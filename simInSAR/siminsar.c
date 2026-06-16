#include "stdio.h"
#include "string.h"

#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"
#include <sys/types.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

/*
  Program to simulate InSAR image including both terrain and motion effects.
*/
static void readMaskFile(char *shelfMaskFile, ShelfMask *shelfMask);

static void readArgs(int argc, char *argv[], sceneStructure *scene,
					 char **demFile, char **displacementFile, char **sceneFile, char **outputFile);
static void usage();

//float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
//char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
//void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
//void *lBuf1, *lBuf2, *lBuf3, *lBuf4;
//char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2, *SEBuf; /* Buffers for offset and azimuth parameter interpolation, and special culling cases. */
//void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4, *offSEBuffSpace;
//void *lBuf1, *lBuf2, *lBuf3, *lBuf4, *lSEBuf;
int32_t llConserveMem = 999; /* Kluge to maintain backwards compat 9/13/06 */

/*
   Global variables definitions
*/
int32_t BufferSize = BUFFERSIZE;			/* Size of nonoverlap region of the buffer */
int32_t BufferLines = 512;					/* # of lines of nonoverlap in buffer */
int32_t HemiSphere = NORTH;
double Rotation = 45.;
double SLat = -91.0;
//char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;

int main(int argc, char *argv[])
{
	sceneStructure scene;
	demStructure llDem;
	xyDEM xyDem;
	xyDEM verticalCorrection;
	xyVEL xyVel;
	void *dem;
	displacementStructure displacements;
	double x, y, lat, lon;
	int32_t i, j, flatFlag;
	ShelfMask *imageMask;
	char *demFile, *displacementFile, *sceneFile, *outputFile;
	GDALAllRegister();
	if (getenv("OMP_NUM_THREADS") == NULL)
		omp_set_num_threads(4);
	/*
	   Read command line args
	*/
	readArgs(argc, argv, &scene, &demFile, &displacementFile, &sceneFile, &outputFile);
	
	/* Added August 2021 to force projections to match dem */
	fprintf(stderr, "outFile %s\n", outputFile);
	
	readXYDEMGeoInfo(demFile, &xyDem, TRUE);
	imageMask = malloc(sizeof(ShelfMask));
	/*
	  Input scene parameters from sceneFile
	*/
	scene.I.stateFlag = TRUE;
	fprintf(stderr, "Parsing scene file (%s)...\n", sceneFile);
	parseSceneFile(sceneFile, &scene);
	/*
	  Compute scene bounding box so we only read the relevant portion of the
	  DEM and velocity map (xyDem.rot/stdLat set by readXYDEMGeoInfo above).
	*/
	double xMin, xMax, yMin, yMax;
	simInSARDEMBounds(&scene, xyDem.rot, xyDem.stdLat, &xMin, &xMax, &yMin, &yMax);
	/*
	  Init and input DEM
	*/
	if (scene.offsetFlag == TRUE || scene.useVelocity == TRUE)
	{
		if (scene.offsetFlag == TRUE)
			fprintf(stderr, "Using offsets\n");
		readXYCropVel(&xyVel, displacementFile, xMin, xMax, yMin, yMax);
	}
	else
	{
		xyVel.xSize = 0;
		xyVel.ySize = 0;
	}

	if (scene.velThresh > 0)
		scene.maskFlag = TRUE;

	fprintf(stderr, "Loading DEM...\n");
	readXYDEMcrop(demFile, &xyDem, xMin, xMax, yMin, yMax);
	dem = (void *)&xyDem;

	if (scene.verticalCorrectionFile != NULL)
	{
		fprintf(stderr, "Loading vertical correction (%s)...\n", scene.verticalCorrectionFile);
		readXYDEMcrop(scene.verticalCorrectionFile, &verticalCorrection, xMin, xMax, yMin, yMax);
		scene.verticalCorrection = &verticalCorrection;
	}
	else
	{
		scene.verticalCorrection = NULL;
	}
	if (scene.maskFlag == TRUE && scene.velThresh == 0)
	{
		readMaskFile(displacementFile, imageMask);
		scene.imageMask = imageMask;
		fprintf(stderr, "++++ %f  %f\n", imageMask->x0, imageMask->y0);
	}
	/*
	  Simulate InSAR image
	*/
	if (scene.geodat2File != NULL)
		simInSARBaselineFromSV(scene.geodat2File, &scene);
	fprintf(stderr, "Running simulation....\n");
	simInSARimage(&scene, dem, &xyVel);
	/*
	  Output image
	*/
	if (scene.heightFlag == TRUE)
		fprintf(stderr, "Outputing DEM\n");
	fprintf(stderr, "Writing results....\n");
	
	outputSimulatedImage(scene, outputFile, demFile, displacementFile);
	fprintf(stderr, "Done\n");
}

int parseBnBpParamsFile(const char *bpParamsFile,
                        double *bn, double *bp,
                        double *dbn, double *dbp)
{
    FILE *fp;
	double x;
    fp = fopen(bpParamsFile, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: could not open %s\n", bpParamsFile);
        return -1;
    }

    if (fscanf(fp, "%lf %lf %lf %lf %lf", bn, bp, dbn, &x, dbp) != 5) {
        fprintf(stderr, "Error: failed to read bn bp dbn dbp from %s\n", bpParamsFile);
        fclose(fp);
        return -1;
    }

    fclose(fp);
    return 0;
}

static void readArgs(int argc, char *argv[], sceneStructure *scene, char **demFile, char **displacementFile,
					 char **sceneFile, char **outputFile)
{
	int32_t filenameArg;
	char *argString;
	char bpParamsFile[2048];
	double bn, bp;
	double bnStart, bnEnd;
	double bpStart, bpEnd;
	double bnMid, dBn, bpMid, dBp;
	double dT;
	int32_t velocityFlag;
	float velThresh;

	int32_t bnFlag = FALSE, bnStartFlag = FALSE, toLLFlag = FALSE;
	int32_t bpFlag = FALSE, bpStartFlag = FALSE;
	int32_t i, n, flatFlag, heightFlag, maskFlag, offsetFlag, velOnlyFlag, tiffFlag;

	if (argc < 5 || argc > 30)
		usage(); /* Check number of args */
	n = argc - 5;
	bn = DEFAULT_SIM_BN;
	bp = DEFAULT_SIM_BP;
	dBn = 0.0;
	dBp = 0.0;
	flatFlag = FALSE;
	heightFlag = FALSE;
	maskFlag = FALSE;
	offsetFlag = FALSE;
	tiffFlag = FALSE;
	bnStart = bn;
	bnEnd = bn;
	bpStart = bp;
	bpEnd = bp;
	velocityFlag = FALSE;
	velOnlyFlag = FALSE;
	velThresh = 0.0;
	dT = 12.0;
	scene->llInput = NULL;
	scene->toLLFlag = FALSE;   /* For offsets */
	scene->saveLLFlag = FALSE; /* for phase/geodat */
	scene->byteOrder = MSB;
	scene->geodat2File = NULL;
	scene->bnArray = NULL;
	scene->bpArray = NULL;
	scene->verticalCorrectionFile = NULL;
	for (i = 1; i <= n; i += 2)
	{
		argString = strchr(argv[i], '-');
		if (strstr(argString, "bnStart") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &bnStart);
			bnStartFlag = TRUE;
			bn = bnStart;
			if (bnFlag == TRUE)
				error("readargs: bnStart/bnEnd incompatible bn\n");
		}
		else if (strstr(argString, "bnEnd") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &bnEnd);
			bnStartFlag = TRUE;
			if (bnFlag == TRUE)
				error("readargs: bnStart/bnEnd incompatible bn\n");
		}
		else if (strstr(argString, "bn") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &bn);
			bnFlag = TRUE;
			if (bnStartFlag == TRUE)
				error("readargs: bn incompatible bnStart/bnEnd\n");
		}
		else if (strstr(argString, "dBn") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &dBn);
			bnFlag = TRUE;
			if (bnStartFlag == TRUE)
				error("readargs: dBn incompatible bnStart/bnEnd\n");
		}
		else if (strstr(argString, "bpStart") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &bpStart);
			bpStartFlag = TRUE;
			bp = bpStart;
			if (bpFlag == TRUE)
				error("readargs: bpStart/bpEnd incompatible bp\n");
		} else if (strstr(argString, "bParamsFile") != NULL)
		{
			sscanf(argv[i + 1], "%s", bpParamsFile);
			bnFlag = TRUE;
			bpFlag = TRUE;
			if (parseBnBpParamsFile(bpParamsFile, &bn, &bp, &dBn, &dBp) != 0)
				error("readArgs: failed to parse bParamsFile %s\n", bpParamsFile);
			if (bnStartFlag == TRUE)
				error("readargs: dBn incompatible bnStart/bnEnd\n");
		}
		else if (strstr(argString, "toLL") != NULL)
		{
			scene->toLLFlag = TRUE;
			scene->llInput = argv[i + 1];
		}
		else if (strstr(argString, "saveLL") != NULL)
		{
			scene->saveLLFlag = TRUE;
		}
		else if (strstr(argString, "bpEnd") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &bpEnd);
			bpStartFlag = TRUE;
			if (bpFlag == TRUE)
				error("readargs: bpStart/bpEnd incompatible bp\n");
		}
		else if (strstr(argString, "dT") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &dT);
		}
		else if (strstr(argString, "bp") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &bp);
			bpStart = bp;
			bpEnd = bp;
			bpFlag = TRUE;
			if (bpStartFlag == TRUE)
				error("readargs: bp incompatible bpStart/bpEnd\n");
		}
		else if (strstr(argString, "dBp") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &dBp);
			bpFlag = TRUE;
			if (bpStartFlag == TRUE)
				error("readargs: bp incompatible bpStart/bpEnd\n");
		}
		else if (strstr(argString, "rPix") != NULL)
			sscanf(argv[i + 1], "%lf", &(scene->I.rangePixelSize));
		else if (strstr(argString, "aPix") != NULL)
			sscanf(argv[i + 1], "%lf", &(scene->I.azimuthPixelSize));
		else if (strstr(argString, "center") != NULL)
		{
			i--;
			fprintf(stderr, "ignoring center flag, obsolete\n");
		}
		else if (strstr(argString, "flat") != NULL)
		{
			i--;
			flatFlag = TRUE;
		}
		else if (strstr(argString, "height") != NULL)
		{
			i--;
			heightFlag = TRUE;
		}
		else if (strstr(argString, "mask") != NULL)
		{
			i--;
			maskFlag = TRUE;
		}
		else if (strstr(argString, "offset") != NULL)
		{
			i--;
			offsetFlag = TRUE;
		}
		else if (strstr(argString, "velThresh") != NULL)
		{
			sscanf(argv[i + 1], "%f", &velThresh);
			velocityFlag = TRUE;
		}
		else if (strstr(argString, "velOnly") != NULL)
		{
			i--;
			velOnlyFlag = TRUE;
			velocityFlag = TRUE;
		}
		else if (strstr(argString, "velocity") != NULL)
		{
			i--;
			velocityFlag = TRUE;
		}
		else if (strstr(argString, "slantRangeDEM") != NULL)
		{
			error("obsolete slantRangeDEM flag used");
		}
		else if (strstr(argString, "xyDEM") != NULL)
		{
			i--;
			fprintf(stderr, "xyDEM flag obsolete: all dems xy");
		}
		else if (strstr(argString, "LSB") != NULL)
		{
			i--;
			scene->byteOrder = LSB;
		}
		else if (strstr(argString, "tiff") != NULL)
		{
			i--;
			tiffFlag = TRUE;
		}
		else if (strstr(argString, "geodat2") != NULL)
		{
			scene->geodat2File = argv[i + 1];
		}
		else if (strstr(argString, "verticalCorrection") != NULL)
		{
			scene->verticalCorrectionFile = argv[i + 1];
		}
		else if (strstr(argString, "ompThreads") != NULL)
		{
			int32_t nThreads = 0;
			sscanf(argv[i + 1], "%d", &nThreads);
			if (nThreads > 0) {
				omp_set_num_threads(nThreads);
				fprintf(stderr, "\033[1;3;34mompThreads set to %d\033[0m\n", nThreads);
			} else
				fprintf(stderr, "\033[1;3;34mompThreads using default (%d)\033[0m\n", omp_get_max_threads());
		}
		else
			usage();
	}
	if (velOnlyFlag == TRUE && (bnFlag == TRUE || bpFlag == TRUE || bnStartFlag == TRUE || bpStartFlag == TRUE ||
	                            scene->geodat2File != NULL))
		error("-velOnly cannot be combined with explicit baseline flags or -geodat2\n");
	if (velOnlyFlag == TRUE)
	{
		bn = 0.0; bp = 0.0;
		bnStart = 0.0; bnEnd = 0.0;
		bpStart = 0.0; bpEnd = 0.0;
		dBn = 0.0; dBp = 0.0;
		fprintf(stderr, "velOnly: baseline set to zero (velocity-only interferogram)\n");
	}
	if (scene->geodat2File != NULL && (bnFlag == TRUE || bpFlag == TRUE || bnStartFlag == TRUE || bpStartFlag == TRUE))
		fprintf(stderr, "WARNING: -geodat2 specified with explicit baseline flags; "
		                "explicit bn/bp values will be ignored\n");
	fprintf(stderr, "offsetFlag %i\n", offsetFlag);
	if (flatFlag == TRUE)
		fprintf(stderr, "flat\n");
	fprintf(stderr, "bn,bp  %f  %f\n", bn, bp);
	*demFile = argv[argc - 4];
	*displacementFile = argv[argc - 3];
	*sceneFile = argv[argc - 2];
	*outputFile = argv[argc - 1];
	scene->bn = bn;
	scene->bp = bp;
	scene->dT = dT;
	fprintf(stderr, "dT = %f", dT);
	scene->flatFlag = flatFlag;
	scene->heightFlag = heightFlag;
	scene->maskFlag = maskFlag;
	scene->offsetFlag = offsetFlag;
	scene->tiffFlag = tiffFlag;
	scene->useVelocity = velocityFlag;
	scene->velThresh = velThresh;

	if (bnStartFlag == TRUE)
	{
		scene->bnStart = bnStart;
		scene->bnEnd = bnEnd;
	}
	else
	{
		scene->bnStart = bn - 0.5 * dBn;
		scene->bnEnd = bn + 0.5 * dBn;
	}
	if (bpStartFlag == TRUE)
	{
		scene->bpStart = bpStart;
		scene->bpEnd = bpEnd;
	}
	else
	{
		scene->bpStart = bp - 0.5 * dBp;
		scene->bpEnd = bp + 0.5 * dBp;
	}
	if (strstr(*outputFile, "stdio"))
		*outputFile = NULL;
	fprintf(stderr, "bn, %f %f\n", scene->bnStart, scene->bnEnd);
	fprintf(stderr, "bp, %f %f\n", scene->bpStart, scene->bpEnd);
	return;
}

static void usage()
{
	error("\n\n%s\n\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n%s\n%s\n%s\n%s\n\n%s\n%s\n%s\n%s\n",
			"Simulate interferogram using a DEM",
			"Usage:",
			"siminsar -LSB -bn bn -dBn dBn -bp bp -dBp dBp ",
			"         -bnStart bnStart -bnEnd bnEnd -bpStart bpStart -bpEnd bpEnd ",
			"         -bParamsFile bParamsFile",
			"         -flat -height -rPix rPix -aPix deltA -velocity -velOnly",
			"         -slantRangeDEM -xyDEM -mask -saveLL -toLL file.dat",
			"         -geodat2 geodat2File -verticalCorrection vcFile -ompThreads N",
			"          demFile displacementFile sceneFile outPutImage",
			"where",
			"   LSB             Output results as LSB [MSB]",
			"   mask            = output a mask using displacement file as mask",
			"   toLL file.dat   = save LL on offset defined grid write outputimage.lat,.lon (can output mask or height with this option)",
			"   saveLL          = save LL on the geodat defined grid and write outputimage.lat and .lon",
			"   bn              = normal component of baseline",
			"   dBn             = change in bn over scene (use instead of bnStart/bnEnd)",
			"   bnStart/bnEnd   = bn at start and end of scene",
			"   bp              = parallel component of baseline",
			"   dBp             = change in bp over scene (use instead of bpStart/bpEnd)",
			"   bpStart/bpEnd   = bp at start and end of scene",
			"   bParamsFile     = file with bn, bp, dBn, dBp for simulation",
			"   dT              = time interval for simulation displacement",
			"   flat            = flattened image",
			"   height          = output height values instead of phase",
			"   rPix            = range single look pixel size",
			"   aPix            = azimuth single look pixel size",
			"   velocity        = use velocity",
			"   verticalCorrection vcFile = xyDEM grid (m/yr) of submergence/emergence rate added to simulated phase",
			"   ompThreads N    = number of OpenMP threads [default: 4]",
			"   slantRangeDEM   = use dem of image size in slant range coords",
			"   xyDEM           = xyDEM file with xyDEM.geodat file",
			"   demFile          = dem file in lat/lon, xy, or slant range format",
			"   displacementFile = velocity file (not .vx,.vy for binary) and .vv. or .*. for .tif or .vrt",
			"   sceneFile        = file with location info",
			"   outPutImage      = simulated interferogram",
			"if outputImage == stdio output is to stdout");
	//exit(0);
}

/*
  Process input file for mosaicDEMs
*/
static void readMaskFile(char *shelfMaskFile, ShelfMask *shelfMask)
{
	fprintf(stderr, "ShelfMaskFile %s\n", shelfMaskFile);

	if (has_extension(shelfMaskFile, ".tif") || has_extension(shelfMaskFile, ".vrt"))
	{
		/* GDAL path for .tif / .vrt mask files */
		GDALDatasetH hDS = GDALOpen(shelfMaskFile, GA_ReadOnly);
		if (hDS == NULL)
			error("readMaskFile: cannot open %s with GDAL\n", shelfMaskFile);

		int32_t nCols = GDALGetRasterXSize(hDS);
		int32_t nRows = GDALGetRasterYSize(hDS);
		double gt[6];
		GDALGetGeoTransform(hDS, gt);

		shelfMask->xSize  = nCols;
		shelfMask->ySize  = nRows;
		shelfMask->deltaX = gt[1] * MTOKM;         /* m → km, positive east */
		shelfMask->deltaY = -gt[5] * MTOKM;        /* gt[5]<0 north-up; store positive */
		shelfMask->x0     = gt[0] * MTOKM;         /* left edge in km */
		shelfMask->y0     = (gt[3] + nRows * gt[5]) * MTOKM; /* bottom edge in km */
		shelfMask->rot    = 0;

		/* Detect hemisphere and standard latitude from EPSG */
		const char *projStr = GDALGetProjectionRef(hDS);
		char *projCopy = strdup(projStr);
		char *projCopyOrig = projCopy;  /* OSRImportFromWkt advances projCopy; free the original */
		OGRSpatialReferenceH hSRS = OSRNewSpatialReference(NULL);
		int32_t epsg = 0;
		if (projCopy != NULL && strlen(projCopy) > 0 &&
			OSRImportFromWkt(hSRS, &projCopy) == OGRERR_NONE)
		{
			const char *epsgCode = OSRGetAuthorityCode(hSRS, NULL);
			if (epsgCode != NULL)
				epsg = atoi(epsgCode);
		}
		OSRDestroySpatialReference(hSRS);
		free(projCopyOrig);

		if (epsg == 3413)
		{
			shelfMask->hemisphere = NORTH;
			shelfMask->stdLat     = 70.0;
			shelfMask->rot        = 45.0;
		}
		else if (epsg == 3031)
		{
			shelfMask->hemisphere = SOUTH;
			shelfMask->stdLat     = 71.0;
			shelfMask->rot        = 0.0;
		}
		else
		{
			/* Default to Greenland if unknown */
			fprintf(stderr, "readMaskFile: unknown EPSG %d, defaulting to Greenland (3413)\n", epsg);
			shelfMask->hemisphere = NORTH;
			shelfMask->stdLat     = 70.0;
			shelfMask->rot        = 45.0;
		}

		fprintf(stderr, "** %i %i \n %f %f \n %f %f \n %f %i %f \n",
				shelfMask->xSize, shelfMask->ySize, shelfMask->deltaX, shelfMask->deltaY,
				shelfMask->x0, shelfMask->y0, shelfMask->rot, shelfMask->hemisphere, shelfMask->stdLat);

		/* Read pixel data: GDAL row 0 = top; ShelfMask row 0 = bottom → flip Y */
		unsigned char *tmp = (unsigned char *)malloc(nCols * nRows);
		GDALRasterBandH hBand = GDALGetRasterBand(hDS, 1);
		GDALRasterIO(hBand, GF_Read, 0, 0, nCols, nRows,
					 tmp, nCols, nRows, GDT_Byte, 0, 0);
		GDALClose(hDS);

		unsigned char *storage = (unsigned char *)malloc(nCols * nRows);
		shelfMask->mask = (unsigned char **)malloc(nRows * sizeof(unsigned char *));
		int32_t i;
		/* mask[i] = bottom-up row i → must hold GDAL row (nRows-1-i) */
		for (i = 0; i < nRows; i++)
		{
			shelfMask->mask[i] = &storage[i * nCols];
			memcpy(shelfMask->mask[i], &tmp[(nRows - 1 - i) * nCols], nCols);
		}
		free(tmp);
	}
	else
	{
		/* GrIMP binary path: read geometry from .geodat sidecar, data from binary file */
		FILE *fp;
		char *geodatFile;
		char line[256];
		int32_t lineCount = 0, eod;
		unsigned char *tmp, *tmp1;
		float dum1, dum2;
		int32_t i;

		geodatFile = (char *)malloc(strlen(shelfMaskFile) + 8);
		geodatFile[0] = '\0';
		geodatFile = strcpy(geodatFile, shelfMaskFile);
		geodatFile = strcat(geodatFile, ".geodat");
		fprintf(stderr, "Shelfmask geodat file %s\n", geodatFile);

		fp = openInputFile(geodatFile);
		if (fp == NULL)
			error("*** readShelf: Error opening %s ***\n", geodatFile);

		lineCount = getDataString(fp, lineCount, line, &eod); /* Skip # 2 line */
		lineCount = getDataString(fp, lineCount, line, &eod);
		sscanf(line, "%f %f\n", &dum1, &dum2);
		shelfMask->xSize = (int)dum1;
		shelfMask->ySize = (int)dum2;
		lineCount = getDataString(fp, lineCount, line, &eod);
		sscanf(line, "%lf %lf\n", &(shelfMask->deltaX), &(shelfMask->deltaY));
		shelfMask->deltaX *= MTOKM;
		shelfMask->deltaY *= MTOKM;
		lineCount = getDataString(fp, lineCount, line, &eod);
		sscanf(line, "%lf %lf\n", &(shelfMask->x0), &(shelfMask->y0));
		fprintf(stderr, "%s\n", line);
		fclose(fp);
		shelfMask->rot = 0;
		shelfMask->hemisphere = SOUTH;
		shelfMask->stdLat = 71.0;

		fprintf(stderr, "** %i %i \n %f %f \n %f %f \n %f %i %f \n",
				shelfMask->xSize, shelfMask->ySize, shelfMask->deltaX, shelfMask->deltaY,
				shelfMask->x0, shelfMask->y0, shelfMask->rot, shelfMask->hemisphere, shelfMask->stdLat);

		shelfMask->mask = (unsigned char **)
			malloc(shelfMask->ySize * sizeof(unsigned char *));
		tmp = (unsigned char *)malloc(shelfMask->xSize * shelfMask->ySize * sizeof(unsigned char));

		fp = openInputFile(shelfMaskFile);
		if (fp == NULL)
			error("*** readShelf: Error opening %s ***\n", shelfMaskFile);
		for (i = 0; i < shelfMask->ySize; i++)
		{
			tmp1 = &(tmp[i * shelfMask->xSize]);
			freadBS(tmp1, sizeof(unsigned char), shelfMask->xSize, fp, BYTEFLAG);
			shelfMask->mask[i] = tmp1;
		}
		fprintf(stderr, "%f %f \n", shelfMask->x0, shelfMask->y0);
	}
}
