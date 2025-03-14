#include "stdio.h"
#include "string.h"

#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"
#include <sys/types.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>

/*
  Program to simulate InSAR image including both terrain and motion effects.
*/
static void readMaskFile(char *shelfMaskFile, ShelfMask *shelfMask);

static void readArgs(int argc, char *argv[], sceneStructure *scene,
					 char **demFile, char **displacementFile, char **sceneFile, char **outputFile);
static void usage();

float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
void *lBuf1, *lBuf2, *lBuf3, *lBuf4;
int32_t llConserveMem = 999; /* Kluge to maintain backwards compat 9/13/06 */

/*
   Global variables definitions
*/
int32_t BufferSize = BUFFERSIZE;			/* Size of nonoverlap region of the buffer */
int32_t BufferLines = 512;					/* # of lines of nonoverlap in buffer */
int32_t HemiSphere = NORTH;
double Rotation = 45.;
double SLat = -91.0;
char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;

int main(int argc, char *argv[])
{
	sceneStructure scene;
	demStructure llDem;
	xyDEM xyDem;
	xyVEL xyVel;
	void *dem;
	displacementStructure displacements;
	double x, y, lat, lon;
	int32_t i, j, flatFlag;
	ShelfMask *imageMask;
	char *demFile, *displacementFile, *sceneFile, *outputFile;
	GDALAllRegister();
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
	  Init and input DEM
	*/
	if (scene.offsetFlag == TRUE || scene.useVelocity == TRUE)
	{
		if (scene.offsetFlag == TRUE)
			fprintf(stderr, "Using offsets\n");
		readXYVel(&xyVel, displacementFile);
	}
	else
	{
		xyVel.xSize = 0;
		xyVel.ySize = 0;
	}
	
	fprintf(stderr, "Loading DEM...\n");
	readXYDEM(demFile, &xyDem);
	dem = (void *)&xyDem;
	if (scene.maskFlag == TRUE)
	{
		readMaskFile(displacementFile, imageMask);
		scene.imageMask = imageMask;
		fprintf(stderr, "++++ %f  %f\n", imageMask->x0, imageMask->y0);
	}
	/*
	  Simulate InSAR image
	*/
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

static void readArgs(int argc, char *argv[], sceneStructure *scene, char **demFile, char **displacementFile,
					 char **sceneFile, char **outputFile)
{
	int32_t filenameArg;
	char *argString;
	double bn, bp;
	double bnStart, bnEnd;
	double bpStart, bpEnd;
	double bnMid, dBn, bpMid, dBp;
	double dT;
	int32_t velocityFlag;

	int32_t bnFlag = FALSE, bnStartFlag = FALSE, toLLFlag = FALSE;
	int32_t bpFlag = FALSE, bpStartFlag = FALSE;
	int32_t i, n, flatFlag, heightFlag, maskFlag, offsetFlag;

	if (argc < 5 || argc > 22)
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
	bnStart = bn;
	bnEnd = bn;
	bpStart = bp;
	bpEnd = bp;
	velocityFlag = FALSE;
	dT = 12.0;
	scene->llInput = NULL;
	scene->toLLFlag = FALSE;   /* For offsets */
	scene->saveLLFlag = FALSE; /* for phase/geodat */
	scene->byteOrder = MSB;
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
		else
			usage();
	}
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
	scene->useVelocity = velocityFlag;

	if (bpStartFlag == TRUE)
	{
		scene->bnStart = bnStart;
		scene->bnEnd = bnEnd;
	}
	else
	{
		scene->bnStart = bn - 0.5 * dBn;
		scene->bnEnd = bn + 0.5 * dBn;
	}
	if (bnStartFlag == TRUE)
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
	error(
		"\n\n%s\n\n%s\n\n%s\n%s\n%s\n%s\n%s\n\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n%s\n%s\n%s\n%s\n\n%s\n%s\n%s\n",
		"Simulate interferogram using a DEM",
		"Usage:",
		"siminsar -LSB -bn bn -dBn dBn -bp bp -dBp dBp ",
		"         -bnStart bnStart -bnEnd -bpStart bpStart -bpEnd ",
		"         -flat -height -rPix rPix -aPix deltA -velocity",
		"         -slantRangeDEM -xyDEM -mask -saveLL -toLL file.dat",
		"          demFile displacementFile sceneFile outPutImage",
		"where",
		"   LSB             Output results as LSB [MSB]",
		"   mask          = output a mask using displaceMent file as mask",
		"   toLL file.dat = save LL on offset defined grid write outputimage.lat,.lon (can output mask or height with this option)",
		"   saveLL  	  =  save LL on the geodat defined grid and write outputimage.lat and .lon",
		"   bn            = normal component of baseline",
		"   dBn           = change in bn over scene (use instead of bnStart/end",
		"   bnStart/bnEnd = bn at start and end of scene",
		"   bp            = parallel component of baseline",
		"   dBn           = change in bp over scene (use instead of bpStart/end",
		"   bpStart/bpEnd = bn at start and end of scene",
		"   dT			  = time interval for simulation displacement",
		"   flat          = flattened image",
		"   height        = output height values instead of phase",
		"   rPix          = range single look pixel size",
		"   aPix          = azimuth single look pixel size",
		"   velocity      = use velocity",
		"   slantRangeDEM = use dem of image size in slant range coords",
		"   xyDEM         = xyDEM file with xyDEM.geodat file",
		"   demFile          = dem file in lat/lon, xy, or slant range format",
		"   displacementFile = velocity file (not .vx,.vy for biary) and .vv. or .*. for .tif or .vrt",
		"   sceneFile        = file with location info",
		"   outPutImage      = simulated interferogram",
		"if outputImage == stdio output is to stdout");
}

/*
  Process input file for mosaicDEMs
*/
static void readMaskFile(char *shelfMaskFile, ShelfMask *shelfMask)
{
	FILE *fp;
	char *geodatFile;
	char line[256];
	int32_t lineCount = 0, eod;
	unsigned char *tmp, *tmp1;
	float dum1, dum2;
	int32_t i;
	/*
	  geodat file name
	*/
	fprintf(stderr, "ShelfMaskFile %s\n", shelfMaskFile);
	geodatFile = (char *)malloc(strlen(shelfMaskFile) + 8);
	geodatFile[0] = '\0';
	geodatFile = strcpy(geodatFile, shelfMaskFile);
	geodatFile = strcat(geodatFile, ".geodat");
	fprintf(stderr, "Shelfmask geodat file %s\n", geodatFile);
	/*
	   Open geodat file
	*/
	fp = openInputFile(geodatFile);
	if (fp == NULL)
		error("*** readShelf: Error opening %s ***\n", geodatFile);
	/*
	  Read parameters
	*/
	lineCount = getDataString(fp, lineCount, line, &eod); /* Skip # 2 line */
	lineCount = getDataString(fp, lineCount, line, &eod);
	sscanf(line, "%f %f\n", &dum1, &dum2); /* read as float in case fp value */
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
	/*
	  Malloc array
	*/
	shelfMask->mask = (unsigned char **)
		malloc(shelfMask->ySize * sizeof(unsigned char *));
	tmp = (unsigned char *)malloc(shelfMask->xSize * shelfMask->ySize * sizeof(unsigned char));
	/*
	   Open shelfFile file
	*/
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
