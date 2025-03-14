#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "geomosaic.h"
#include <sys/types.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include "gdalIO/gdalIO/grimpgdal.h"
#include "gdal.h"
#include "ogr_srs_api.h"

#define PSISAVE 2
#define GAMMACORSAVE 4
#define GAMMASAVE 8
/*
  Mosaic several insar dems with altimetry dem.
*/
static void parseAntPat(char *antPatFile, inputImageStructure *inputImage);
static void readArgs(int argc, char *argv[], char **inputFile, char **demFile, char **outFile, float *fl, int *removePad,
					 int32_t *nearestDate, int32_t *noPower, int32_t *hybridZ, int32_t *rsatFineCal, int32_t *S1Cal, char **date1, char **date2,
					 int32_t *smoothL, int32_t *smoothOut, int32_t *orbitPriority, float *noData, char **driver, int32_t *byteScale);
static void usage();
static void processMosaicDateGeo(outputImageStructure *outputImage, char *date1, char *date2);
static void parseBetaNought(inputImageStructure *inputImage);
static void parseImages(inputImageStructure *inputImages, outputImageStructure *outputImage, char **imageFiles, char **geodatFiles,
						float *weights, char **antPatFiles, int32_t nFiles, int32_t S1Cal, int32_t *maxR, int32_t *maxA, float noData);
static void outputBounds(inputImageStructure *inputImage, outputImageStructure *outputImage, int32_t nFiles);
static void memAllocGeomosaic(inputImageStructure *inputImage, outputImageStructure *outputImage,
							  int32_t maxR, int32_t maxA, int32_t nFiles, int32_t removePad);
static void outputS1Cal(outputImageStructure outputImage, char *outputFile, float **psi, float **gamma, int32_t s1Cal, char *driver);
static void byteScaleImage(outputImageStructure *outputImage, double lowerBound, double upperBound, double scale, double exponent);

typedef struct {
    double lowerBound;
    double upperBound;
    double scale;
    double exponent;
} byteScaleParameters;

byteScaleParameters byteScaleParams = {
	.lowerBound = 0.53,
	.upperBound = 2.4,
	.scale = 0.00015 / 0.13,
	.exponent = 0.2
};

/*
   Global variables definitions
*/
void *lBuf1, *lBuf2, *lBuf3, *lBuf4;
void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
int32_t RangeSize = RANGESIZE;				/* Range size of complex image */
int32_t AzimuthSize = AZIMUTHSIZE;			/* Azimuth size of complex image */
int32_t BufferSize = BUFFERSIZE;			/* Size of nonoverlap region of the buffer */
int32_t BufferLines = 512;					/* # of lines of nonoverlap in buffer */
int32_t HemiSphere = NORTH;
double Rotation = 45.;
double SLat = -91.0;
int32_t nearestDate = -1;
int32_t hybridZ = -1;
int32_t noPower = -1;
int32_t rsatFineCal = FALSE;
int32_t S1Cal = FALSE;
char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
int32_t llConserveMem = 1234;		/* Kluge to maintain backwards compat 9/13/06 */
float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 not use only for mosaic3d compatability */
float *smoothBuf;

int main(int argc, char *argv[])
{
	extern int32_t nearestDate;
	extern int32_t noPower;
	extern char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
	extern float *smoothBuf;
	float **psi, **gamma, *fpImage;
	void *dem;
	xyDEM xyDem;
	outputImageStructure outputImage;
	inputImageStructure *inputImage; /* Linked list of input images */
	float *weights;
	double x1, y1, x2, y2;
	float fl;
	int32_t GTiff, COG, byteScale, dataType;
	int32_t nFiles, maxR, maxA;
	int32_t smoothL, smoothOut, orbitPriority, removePad;
	int32_t i, j; /* LCV */
	float noData;
	char *date1, *date2; /* Date range */
	char *demFile, *inputFile, *outFile;
	char **imageFiles, **geodatFiles, **antPatFiles;
	char tmp[2048];
	char *driver;
	/*
	   Read command line args and compute filenames
	*/
	GDALAllRegister();
	smoothBuf = NULL;
	readArgs(argc, argv, &inputFile, &demFile, &outFile, &fl, &removePad, &nearestDate, &noPower,
			 &hybridZ, &rsatFineCal, &S1Cal, &date1, &date2, &smoothL, &smoothOut, &orbitPriority, &noData, &driver, &byteScale);
	/* This step just reads in the dem projection info, which is then used for the outputs */
	readXYDEMGeoInfo(demFile, &xyDem, TRUE);
	outputImage.slat = xyDem.stdLat;
	processMosaicDateGeo(&outputImage, date1, date2);
	/*
	  read inputfile (uses routine from mosaicDEMS).
	*/
	processInputFileGeo(inputFile, &imageFiles, &geodatFiles, &outputImage, &nFiles, &weights, &antPatFiles);
	/*
	  Malloc input images
	*/
	inputImage = (inputImageStructure *)malloc((size_t)(sizeof(inputImageStructure) * nFiles));
	/*
	  Parse input file and input data for each image.
	*/
	parseImages(inputImage, &outputImage, imageFiles, geodatFiles, weights, antPatFiles, nFiles,
				S1Cal & TRUE, &maxR, &maxA, noData);
	/* Modified Aug 2021 for alternate ps projections */

	fprintf(stderr, "Rotation, slat*  %f %f\n", Rotation, outputImage.slat);
	fprintf(stderr, "maxR, maxA %i %i\n", maxR, maxA);
	/*
	  Find output bounds
	*/
	outputBounds(inputImage, &outputImage, nFiles);
	/*
	   Memory allocation
	*/
	memAllocGeomosaic(inputImage, &outputImage, maxR, maxA, nFiles, removePad);
	/*
	  Init and input DEM
	*/
	x1 = (outputImage.originX - 10.e3) / 1000.;
	y1 = (outputImage.originY - 10.e3) / 1000.;
	x2 = (outputImage.originX + outputImage.deltaX * outputImage.xSize + 10e3) / 1000.;
	y2 = (outputImage.originY + outputImage.deltaY * outputImage.ySize + 10e3) / 1000.;
	fprintf(stderr, "x1, x2, y1, y2, %f %f %f %f\n", x1, x2, y1, y2);
	readXYDEMcrop(demFile, &xyDem, x1, x2, y1, y2);
	dem = (void *)&xyDem;
	/*
	  Do the mosaicking
	*/
	makeGeoMosaic(inputImage, outputImage, dem, nFiles, maxR, maxA, imageFiles, fl, smoothL,
				  smoothOut, orbitPriority, &psi, &gamma);
	fprintf(stderr, "Outputting result....");
	/*
	  Output result
	*/
	if (S1Cal > 0)
		outputS1Cal(outputImage, outFile, psi, gamma, S1Cal, driver);
	else 
	{
		if(driver == NULL)
		{	
			outputGeocodedImage(outputImage, outFile);
		}
		else
		{
			dictNode *summaryMetaData = NULL;
			const char *epsg = getEPSGFromProjectionParams(Rotation, SLat, HemiSphere);
			if(byteScale == FALSE) {
				dataType = GDT_Float32;
			} 
			else
			{
				dataType = GDT_Byte;
				// Constants need to made user definable.
				byteScaleImage(&outputImage, 0.53, 2.4, 0.00015 / 0.13, 0.2);
			}
			outputGeocodedImageTiff(outputImage, outFile, driver, epsg, summaryMetaData, 0., dataType);
		}
		
	}
	/* For now only save incidence angle buffer if S1Cal */
}

static void byteScaleImage(outputImageStructure *outputImage, double lowerBound, double upperBound, double scale, double exponent)
{
	float **originalImage, value;
	unsigned char **scaledImage;
	unsigned char *byteBuff;
	extern byteScaleParameters byteScaleParams;

	originalImage = (float **)outputImage->image;
	byteBuff = (unsigned char *)malloc((size_t)(outputImage->xSize * outputImage->ySize * sizeof(unsigned char)));
	scaledImage  = (unsigned char **)malloc((size_t)(outputImage->ySize * sizeof(unsigned char *)));
	for (int32_t i = 0; i < outputImage->ySize; i++)
	{
		scaledImage[i] = (unsigned char *)&(byteBuff[i * outputImage->xSize]);
		for(int32_t j = 0; j < outputImage->xSize; j++ )
		{
			if(originalImage[i][j] == 0.0)
			{
				scaledImage[i][j] = 0;
			} 
			else
			{
			    value = pow((byteScaleParams.scale * originalImage[i][j]), byteScaleParams.exponent);
				value = (value - byteScaleParams.lowerBound) * (1. /(byteScaleParams.upperBound - byteScaleParams.lowerBound)) * 255.;
				value = min(255, max(1, value));
				scaledImage[i][j] = (unsigned char) value;
			}
		}
	}
	outputImage->image = (void **)scaledImage;
}


static void outputS1Cal(outputImageStructure outputImage, char *outFile, float **psi, float **gamma, int32_t S1Cal, char *driver)
{
	char tmp[2048], *psiFile, *gFile, *sigFile;
	float **tmpFloat;
	int32_t i, j;
	const char *epsg = getEPSGFromProjectionParams(Rotation, SLat, HemiSphere);	
	int32_t dataType = GDT_Float32;
	/*
	  output image
	 */
	sigFile = appendSuffix(outFile, ".sigma0", tmp);
	if(driver == NULL)
	{
		outputGeocodedImage(outputImage, sigFile);
	}
	else
	{
		dictNode *summaryMetaData = NULL;
		outputGeocodedImageTiff(outputImage, sigFile, driver, epsg, summaryMetaData, -30., dataType);
	}
	fprintf(stderr, "%i %i %i\n", S1Cal, S1Cal & PSISAVE, S1Cal & GAMMACORSAVE);
	/*
	  Output gamma
	*/
	if ((S1Cal & GAMMASAVE) > 0)
	{
		gFile = appendSuffix(outFile, ".gamma0", tmp);
		fprintf(stderr, "gammaFile %s\n", gFile);
		tmpFloat = (float **)outputImage.image;
		for (i = 0; i < outputImage.ySize; i++)
		{
			for (j = 0; j < outputImage.xSize; j++)
			{
				/* catch low values */
				if (tmpFloat[i][j] > MINS1DB)
					tmpFloat[i][j] += gamma[i][j];
				if (gamma[i][j] <= MINS1DB)
					tmpFloat[i][j] = MINS1DB;
				tmpFloat[i][j] = max(roundf(tmpFloat[i][j] * 100) / 100., MINS1DB);
			}
		}
		if(driver == NULL)
		{
			outputGeocodedImage(outputImage, gFile);
		}
		else
		{
			dictNode *summaryMetaData = NULL;
			outputGeocodedImageTiff(outputImage, gFile, driver, epsg, summaryMetaData, -30., dataType);
		}
	}
	/* Output psi 	*/
	if ((S1Cal & PSISAVE) > 0)
	{
		outputImage.image = (void **)psi;
		psiFile = appendSuffix(outFile, ".inc", tmp);
		fprintf(stderr, "psiFile %s\n", psiFile);
		outputGeocodedImage(outputImage, psiFile);
	}
	/*  Output save gamma	*/
	if ((S1Cal & GAMMACORSAVE) > 0)
	{
		outputImage.image = (void **)gamma;
		gFile = appendSuffix(outFile, ".gamcor", tmp);
		fprintf(stderr, "gammaCorrFile %s\n", gFile);
		outputGeocodedImage(outputImage, gFile);
	}
}

static void memAllocGeomosaic(inputImageStructure *inputImage, outputImageStructure *outputImage,
							  int32_t maxR, int32_t maxA, int32_t nFiles, int32_t removePad)
{

	extern char *Dbuf1, *Dbuf2;
	extern float *smoothBuf;
	float *buf1, *buf1s;
	int32_t i, j;
	/*
	  Malloc space for input images
	*/
	/*Dbuf1=malloc((size_t)MAXADBUF);  add these 9/13/06 for initlltoimage */
	Dbuf1 = NULL;
	Dbuf2 = malloc((size_t)MAXADBUF);
	buf1 = (float *)malloc((size_t)(sizeof(float) * maxR * maxA));
	fprintf(stderr, "mallocing image buf1 %f\n", (sizeof(float) * maxR * maxA) / 1e6);
	smoothBuf = malloc((size_t)(sizeof(float) * max(maxR, maxA)));
	if (buf1 == NULL)
		error("Unable to malloc input image buffer\n");
	for (i = 0; i < nFiles; i++)
	{
		inputImage[i].image = (void **)malloc((size_t)(inputImage[i].azimuthSize * sizeof(float *)));
		inputImage[i].removePad = removePad;
		inputImage[i].memChan = MEM1; /* added 8/28/2017 to avoid using Abuf*/
		for (j = 0; j < inputImage[i].azimuthSize; j++)
		{
			inputImage[i].image[j] = &(buf1[j * inputImage[i].rangeSize]);
		}
	}
	/*
	  Malloc space for output images
	 */
	outputImage->imageType = POWER;
	outputImage->image = (void **)malloc((size_t)(outputImage->ySize * sizeof(float *)));
	outputImage->scale = (float **)malloc((size_t)(outputImage->ySize * sizeof(float *)));
	buf1 = (float *)malloc((size_t)(outputImage->xSize * outputImage->ySize * sizeof(float)));
	if (buf1 != NULL)
		fprintf(stderr, "Malloced %lu image buffer \n", outputImage->xSize * outputImage->ySize * sizeof(float));
	else
		error("Malloc failed for image buffer of size %lu\n", outputImage->xSize * outputImage->ySize * sizeof(float));

	buf1s = (float *)malloc((size_t)(outputImage->xSize * outputImage->ySize * sizeof(float)));
	if (buf1 != NULL)
		fprintf(stderr, "Malloced %lu scale buffer \n", outputImage->xSize * outputImage->ySize * sizeof(float));
	else
		error("Malloc failed for scale buffer of size %i\n", outputImage->xSize * outputImage->ySize * sizeof(float));
	fprintf(stderr, "mallocing buf2 %f\n", outputImage->xSize * outputImage->ySize * sizeof(float) / 1e6);
	for (i = 0; i < outputImage->ySize; i++)
	{
		outputImage->image[i] = (void *)&(buf1[i * outputImage->xSize]);
		outputImage->scale[i] = (float *)&(buf1s[i * outputImage->xSize]);
	}
}

static void outputBounds(inputImageStructure *inputImage, outputImageStructure *outputImage, int32_t nFiles)
{
	double minX, maxX, minY, maxY, x, y;
	int32_t i, j;
	/* Loop to find min/max for input image */
	minX = 1.e30;
	minY = 1.0e30;
	maxX = -1.e30;
	maxY = -1.0e30;
	for (i = 0; i < nFiles; i++)
	{
		for (j = 1; j < 5; j++)
		{
			lltoxy1(inputImage[i].latControlPoints[j], inputImage[i].lonControlPoints[j], &x, &y, Rotation, outputImage->slat);
			minX = min(x, minX);
			minY = min(y, minY);
			maxX = max(x, maxX);
			maxY = max(y, maxY);
		}
	}
	/* Pad and set as output range */
	minX = (double)((int32_t)minX - 3);
	maxX = (double)((int32_t)maxX + 3);
	minY = (double)((int32_t)minY - 3);
	maxY = (double)((int32_t)maxY + 3);
	fprintf(stderr, "minX, maxX, minY, maxY %f %f %f %f\n", minX, maxX, minY, maxY);
	/* if xSize or ySize == 0, then autosize, otherwise use preset values */
	if (outputImage->xSize == 0 || outputImage->ySize == 0)
	{
		fprintf(stderr, "\n\n*** AUTOSIZING REGION ***\n\n");
		outputImage->xSize = (int)((maxX - minX) / (outputImage->deltaX * MTOKM) + 0.5);
		outputImage->ySize = (int)((maxY - minY) / (outputImage->deltaY * MTOKM) + 0.5);
		outputImage->originX = minX * KMTOM;
		outputImage->originY = minY * KMTOM;
	}
	fprintf(stderr, "%i %i %f %f %f %f\n", outputImage->xSize, outputImage->ySize,
			outputImage->originX, outputImage->originY, outputImage->deltaX, outputImage->deltaY);
}

static void parseImages(inputImageStructure *inputImages, outputImageStructure *outputImage, char **imageFiles, char **geodatFiles,
						float *weights, char **antPatFiles, int32_t nFiles, int32_t S1Cal, int32_t *maxR, int32_t *maxA, float noData)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	inputImageStructure *inputImage;
	double jd;
	float t1, t2, t3, t4;
	int32_t i;
	*maxA = 0;
	*maxR = 0;
	/* outputImage->slat = 70.0;*/
	for (i = 0; i < nFiles; i++)
	{
		inputImage = &(inputImages[i]);
		inputImage->passType = DESCENDING; /* default */
		inputImage->stateFlag = TRUE;
		inputImage->weight = weights[i];
		inputImage->file = imageFiles[i];
		inputImage->noData = noData;
		fprintf(stderr, "Input file %i:: %s %s %f\n", i + 1, geodatFiles[i], imageFiles[i], weights[i]);
		inputImage->betaNought = -1.0;
		if (S1Cal == TRUE)
			parseBetaNought(inputImage);
		parseInputFile(geodatFiles[i], inputImage);
		inputImage->pAnt = NULL;
		inputImage->rAnt = NULL;
		if (antPatFiles[i] != NULL)
		{
			fprintf(stderr, "parse %s\n", antPatFiles[i]);
			parseAntPat(antPatFiles[i], inputImage);
		}
		jd = juldayDouble((int)inputImage->month, (int)inputImage->day, (int)inputImage->year);
		if (jd < outputImage->jd1 || jd > outputImage->jd2)
		{
			inputImage->weight = 0; /* Use weight to indicate outside of range */
		}
		else
		{
			*maxR = max(*maxR, inputImage->rangeSize);
			*maxA = max(*maxA, inputImage->azimuthSize);
			fprintf(stderr, "mR, mA %i %i %i %i\n", *maxR, *maxA, inputImage->rangeSize, inputImage->azimuthSize);
		}
		t1 = secondForSAR(&(inputImage->par)); /* + inputImage.par.azoff/ inputImage.par.prf; why was this still here*/
		t2 = t1 + (inputImage->azimuthSize * inputImage->nAzimuthLooks) / inputImage->par.prf;
		t3 = inputImage->sv.times[0];
		t4 = inputImage->sv.times[inputImage->sv.nState];
		if (t3 > (t1 + 45) || (t4 + 45) < t2)
			fprintf(stdout, "%s  %f %f %lf %lf %lf %lf\n", geodatFiles[i], t1, t2, t3, t4, t1 - t3, t4 - t2);
	}
	if (HemiSphere == NORTH)
		fprintf(stderr, "\n ***** Northern Hemisphere  ******\n");
	else
		fprintf(stderr, "\n ***** Southern Hemisphere  ******\n");
}

static void parseAntPat(char *antPatFile, inputImageStructure *inputImage)
{
	FILE *fp;
	int32_t lineCount, eod;
	char line[256];
	/*
	  Open file for input
	*/
	/* use poly */
	if (strstr(antPatFile, "poly") != NULL)
	{
		inputImage->patSize = RSATFINE;
		return;
	}
	if (strstr(antPatFile, "alos") != NULL)
	{
		inputImage->patSize = ALOS;
		return;
	}

	fp = openInputFile(antPatFile);
	fprintf(stderr, "Ant Pat File %s\n", antPatFile);
	eod = 0;
	lineCount = 0;
	while (!eod)
		lineCount = getDataString(fp, lineCount, line, &eod);
	fprintf(stderr, "linecount = %i\n", lineCount);
	inputImage->rAnt = (float *)malloc((size_t)(sizeof(float) * (lineCount - 1)));
	inputImage->pAnt = (float *)malloc((size_t)(sizeof(float) * (lineCount - 1)));
	inputImage->patSize = lineCount - 1;
	rewind(fp);
	lineCount = 0;
	while (lineCount < inputImage->patSize)
	{
		lineCount = getDataString(fp, lineCount, line, &eod);
		sscanf(line, "%f %f\n", &(inputImage->rAnt[lineCount - 1]), &(inputImage->pAnt[lineCount - 1]));
		/*         fprintf(stderr, "%f %f\n", inputImage->rAnt[lineCount], inputImage->pAnt[lineCount]);*/
	}
}

static void readArgs(int argc, char *argv[], char **inputFile, char **demFile, char **outFile, float *fl, int32_t *removePad,
					 int32_t *nearestDate, int32_t *noPower, int32_t *hybridZ, int32_t *rsatFineCal, int32_t *S1Cal, char **date1, char **date2,
					 int32_t *smoothL, int32_t *smoothOut, int32_t *orbitPriority, float *noData, char **driver, int32_t *byteScale)
{
	int32_t filenameArg;
	char *argString;
	char stringbuf[128];
	char *s;
	int32_t month, day, year;
	int32_t doy[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 333};
	int32_t i, n;

	if (argc < 4 || argc > 28)
	{
		fprintf(stderr, "To many/few args %i", argc);
		usage();
	} /* Check number of args */
	n = argc - 4;
	*fl = 0.0;
	*removePad = 0;
	*nearestDate = -1;
	*noPower = -1;
	*rsatFineCal = FALSE;
	*S1Cal = FALSE;
	*date1 = NULL;
	*date2 = NULL;
	*noData = 0.0;
	*smoothL = 0;
	*smoothOut = 1;
	stringbuf[0] = '\0';
	*orbitPriority = -1;
	*driver = NULL;
	*byteScale = FALSE;
	for (i = 1; i <= n; i++)
	{
		argString = strchr(argv[i], '-');
		if (strstr(argString, "xyDEM") != NULL)
			fprintf(stderr, "xyDEM obsolete flag");
		else if (strstr(argString, "fl") != NULL)
		{
			sscanf(argv[i + 1], "%f", fl);
			i++;
		}
		else if (strstr(argString, "smoothL") != NULL)
		{
			sscanf(argv[i + 1], "%i", smoothL);
			i++;
		}
		else if (strstr(argString, "smoothOut") != NULL)
		{
			sscanf(argv[i + 1], "%i", smoothOut);
			i++;
			if (*smoothOut < 1 || *smoothOut > 1)
				error("smoothOut must be in range 1 (can recompile for up to 3)");
		}
		else if (strstr(argString, "noData") != NULL)
		{
			sscanf(argv[i + 1], "%f", noData);
			i++;
		}
		else if (strstr(argString, "nearestDate") != NULL)
		{
			sscanf(argv[i + 1], "%s", stringbuf);
			i++;
		}
		else if (strstr(argString, "hybridZ") != NULL)
		{
			sscanf(argv[i + 1], "%i", hybridZ);
			i++;
		}
		else if (strstr(argString, "rsatFineCal") != NULL)
		{
			*rsatFineCal = 1;
		} 
		else if (strstr(argString, "COG") != NULL)
		{
			*driver = "COG";
		}	
		else if (strstr(argString, "GTiff") != NULL)
		{
			*driver = "GTiff";
		}	
		else if (strstr(argString, "byteScale") != NULL)
		{
			*byteScale = TRUE;
		}		
		else if (strstr(argString, "S1Cal") != NULL)
		{
			*S1Cal = TRUE;
			*S1Cal |= GAMMASAVE;
		}
		else if (strstr(argString, "S1GammaCorr") != NULL)
		{
			*S1Cal |= GAMMACORSAVE;
		}
		else if (strstr(argString, "S1Psi") != NULL)
		{
			*S1Cal |= PSISAVE;
		}
		else if (strstr(argString, "noPower") != NULL)
		{
			*noPower = TRUE;
		}
		else if (strstr(argString, "ascending") != NULL)
		{
			if (*orbitPriority > -1)
			{
				fprintf(stderr, "\n\nCan't specify both ascending and descending priority \n\n");
				usage();
			}
			*orbitPriority = ASCENDING;
		}
		else if (strstr(argString, "descending") != NULL)
		{
			if (*orbitPriority > -1)
			{
				fprintf(stderr, "\n\nCan't specify both ascending and descending priority \n\n");
				usage();
			}
			*orbitPriority = DESCENDING;
		}
		else if (strstr(argString, "removePad") != NULL)
		{
			sscanf(argv[i + 1], "%i", removePad);
			i++;
		}
		else if (strstr(argString, "date1") != NULL)
		{
			*date1 = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "date2") != NULL)
		{
			*date2 = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "center"))
		{
			fprintf(stderr, "obsolute center flag - can remove");
		}
		else if (strstr(argString, "BSlowerBound") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &byteScaleParams.lowerBound);
			i++;
		}
		else if (strstr(argString, "BSupperBound") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &byteScaleParams.upperBound);
			i++;
		}
		else if (strstr(argString, "BSexponent") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &byteScaleParams.exponent);
			i++;
		}
		else if (strstr(argString, "BSscale") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &byteScaleParams.scale);
			i++;
		}
		else
		{
			fprintf(stderr, "\n\narg %s not parsed \n\n", argString);
			usage();
		} 	
		
	}
	if (strlen(stringbuf) > 8)
	{
		s = strtok(stringbuf, "-");
		if (s != NULL)
			sscanf(s, "%d", &month);
		else
		{
			fprintf(stderr, "error parsing date\n\n");
			usage();
		}
		s = strtok(NULL, "-");
		if (s != NULL)
			sscanf(s, "%d", &day);
		else
		{
			fprintf(stderr, "error parsing date\n\n");
			usage();
		}
		s = strtok(NULL, "-");
		if (s != NULL)
			sscanf(s, "%d", &year);
		else
		{
			fprintf(stderr, "error parsing date\n\n");
			usage();
		}
		if (month > 12 || day > 31 || year < 1990 || year > 2050)
		{
			fprintf(stderr, "invalid day month or year");
			usage();
		}
		fprintf(stderr, "%i %i %i\n", year, month, day);
		*nearestDate = year * 365. + doy[month - 1] + day;
		fprintf(stderr, "nearest data %i %i %i\n", month, day, year);
	}

	if(*driver == NULL && *byteScale == TRUE)
	{
		fprintf(stderr,"\nbyte scale only works with tiff output\n");
		usage();
	}
	if (hybridZ > 0 && nearestDate < 0)
		error("hybrid Z requires a nearest date ");
	*inputFile = argv[argc - 3];
	*demFile = argv[argc - 2];
	*outFile = argv[argc - 1];
	if (*smoothL > 1)
		*smoothL = ((int)(*smoothL / 2)) + 1;
	if (*smoothL > 0 && *smoothOut > 1)
		error("Smoothing of only input or output allowed (not both)");
	fprintf(stderr, "Smooth Length %i\n", *smoothL);
	fprintf(stderr, "Inputfile    = %s\n", *inputFile);
	fprintf(stderr, "High Res DEM = %s\n", *demFile);
	fprintf(stderr, "outFile  = %s\n", *outFile);
	fprintf(stderr, "fl           = %f\n", *fl);
	fprintf(stderr, "removePad     = %i\n", *removePad);
	fprintf(stderr, "nearestDate= %i\n", *nearestDate);
	fprintf(stderr, "noPower= %i\n", *noPower);
	fprintf(stderr, "orbitPriority (<0 none; 0 desc; 1 asc) = %i\n", *orbitPriority);
	if(*byteScale == TRUE)
	{
		fprintf(stderr, "\nByte Scale Parameters\n");
		fprintf(stderr, "\tlowerBound %lf\n", byteScaleParams.lowerBound);
		fprintf(stderr, "\tupperBound %lf\n", byteScaleParams.upperBound);
		fprintf(stderr, "\tscale %lf\n", byteScaleParams.scale);
		fprintf(stderr, "\texponent %lf\n\n", byteScaleParams.exponent);
	}

	if (*orbitPriority >= 0 || (*S1Cal & TRUE) == TRUE)
		if (*fl > 0)
			error("Can use fl with orbitPriority or S1Cal");

	if (*rsatFineCal == TRUE)
		fprintf(stderr, "rsatFineCal = TRUE\n");
	else
		fprintf(stderr, "rsatFineCal = FALSE\n");
	if (*S1Cal == TRUE)
		fprintf(stderr, "S1Cal = TRUE\n");
	else
		fprintf(stderr, "S1Cal = FALSE\n");
	if ((*S1Cal & TRUE) == 0)
		*S1Cal &= FALSE; /* Avoid output of sigma related vars if not cal */
	return;
}

static void parseBetaNought(inputImageStructure *inputImage)
{
	/*
	   Parse beta nought conversion factor from file.
	*/
	FILE *fp;
	float betaNought;
	int32_t lineCount, eod;
	char line[2048];
	char betaNoughtFile[2048], *tmp1;
	int32_t i, lastSlash;
	tmp1 = "betaNought";
	for (i = 0; i < strlen(inputImage->file); i++)
	{
		betaNoughtFile[i] = inputImage->file[i];
		if (betaNoughtFile[i] == '/')
			lastSlash = i;
	}
	for (i = 0; i < strlen(tmp1); i++)
		betaNoughtFile[lastSlash + i + 1] = tmp1[i];
	betaNoughtFile[lastSlash + i + 1] = '\0';
	fprintf(stderr, "parseBetaNought\n");
	fprintf(stderr, "%s\n", betaNoughtFile);
	/* Open and read file */
	fp = openInputFile(betaNoughtFile);
	eod = 0;
	lineCount = 0;
	while (lineCount < 1)
	{
		lineCount = getDataString(fp, lineCount, line, &eod);
	}
	fclose(fp);
	sscanf(line, "%f\n", &betaNought);
	fprintf(stderr, "%s %f\n", line, betaNought);
	if (betaNought < 100 || betaNought > 1000)
		error("invalid BetaNought");
	inputImage->betaNought = betaNought;
}

static void usage()
{
	error("\n\n%s\n\n%s\n\n%s\n%s\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
		  "mosaic images*",
		  "Usage:", " \033[1mgeomosaic -rsatFineCal -S1Cal -noPower -noData -descending -ascending -nearestDate "
					"YYYY:MM:DD -hybridZ zthresh \\ \n\t-date1 MM-DD-YYYY -date2 MM-DD-YYYY -smoothL smoothL -smoothOut"
					"smoothOut -fl fl -removePad pad -xyDEM \\",
		  "\tinputFile Demfile outPutImage\033[0m",
		  "where\n",
		   "\tGTiff =\t\t\t Save to geotiff files will add .tif extension if not present",
		  "\tCOG =\t\t\t Save to Cloud Optimized Geotiff will add .tif extension if not present",
		  "\tbyteScale			write data a scaled byte ouput to look nice using empirically determined params",
		  "\tBSlowerBound		lowerBound for byteScale [0.53]",
		  "\tBSupperBound		upperBound for byteScale [2.4]",
		  "\tBSscale			scale factor for byteScale [0.0011538463]",
		  "\tBSexponent			exponent for byteScale [0.2]",
		  "\tfl               	= feather length",
		  "\tdate1, date2 	  	= range of dates to include in mosaic",
		  "\tdescending      	= put descending on top",
		  "\tascending        	= put ascending on top",
		  "\tsmoothL   		= multilook source images smoothL by smoothL (should be odd or will be made odd)",
		  "\tsmoothOut 	   	= oversample image and then smooth for output (1, 2, 3) [1]  - TEST ONLY smoothOut!=1 will exit",
		  "\tnearestDate            = mosaic only points nearest YYYY:MM:DD",
		  "\thybridZ           	 = if used with nearest Date, date will only be applied for elevations below zthresh",
		  "\trsatFineCal       	 = set to output calibrate rsat fine beam data",
		  "\tnoPower    	   	 = non power data",
		  "\tnoData	    	   	 = no data value",
		  "\tS1Cal  	          	 = set to output calibrated S1 Data",
		  "\tS1Psi  	          	 = if S1Cal set, also output inc angle (outputImage.inc)",
		  "\tS1GammaCorr      = if S1Cal set, also output gamma correction (outputImage.gammacorr)",
		  "\tremovePad         	 = remove first pad lines from first and last col",
		  "\txyDEM              	 = smoothDem is XY type with xyDEM.geodat file",
		  "\tinputFile          	 = file with input params, dem and geodat filenames",
		  "\tdemFile                  = file with nonInsar dem");
}

static void processMosaicDateGeo(outputImageStructure *outputImage, char *date1, char *date2)
{
	int32_t m1, m2, d1, d2, y1, y2;
	/* Process dates and convert to julian dates */
	if (date1 != NULL)
	{
		if (sscanf(date1, "%d-%d-%d", &m1, &d1, &y1) != 3)
			error("invalid date %s\n", date1);
		outputImage->jd1 = juldayDouble(m1, d1, y1);
	}
	else
		outputImage->jd1 = 0.0;
	if (date2 != NULL)
	{
		if (sscanf(date2, "%d-%d-%d", &m2, &d2, &y2) != 3)
			error("invalid date %s\n", date2);
		outputImage->jd2 = (double)juldayDouble(m2, d2, y2) + 0.9999999; /* add 0.999 to ensure end of day */
	}
	else
		outputImage->jd2 = HIGHJD + 1000.; /* way past present */
	if ((date1 != NULL || date2 != NULL) && (outputImage->jd2 < outputImage->jd1))
		error("date 1 follows date2\n");
	if (date1 != NULL)
		fprintf(stderr, "date1= %s %f\n", date1, outputImage->jd1);
	if (date2 != NULL)
		fprintf(stderr, "date2= %s %f\n", date2, outputImage->jd2);
}
