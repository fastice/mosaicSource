#include <stdio.h>
#include "string.h"
#include "mosaicSource/common/common.h"
#include "mosaic3d.h"
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
#include "mosaicSource/landsatMosaic/landSatMosaic.h"
#include "gdalIO/gdalIO/grimpgdal.h"
/*
  Mosaic several insar dems with altimetry dem.


  5/30/07 - Major revision for the 3D part. Modified so that decision to
  derive velocity from two images is base on the track heading rather than asc/desc
  pairs. This required significant changes to the mem conservation scheme. Specifically,
  instead of saving memory seperately for ASC and DESC images, the memory now needs to
  be designated by first or second image in the pair (which will shift). So a memChan
  variable was added to inputImage. If set, this will override the asc/desc allocation.
  So this gets explicitly set before each image is started in the two loops. Likewise
  the image buffers are still pre-allocated, but two global pointers (AImageBuffer and DImageBuffer)
  are used to carry around the locations of the blocks of memory. Prior to reading each image, the
  buffer pointers are set as needed with explicit routine (setBuffer).
*/
static void computeDateRange(char **date1, char **date2, outputImageStructure *outputImage, inputImageStructure *images, vhParams *params);
static void removeOutOfBounds(outputImageStructure *outputImage, inputImageStructure **ascImages, vhParams **ascParams, int32_t *nImages);
static void findOutBounds(outputImageStructure *outputImage, inputImageStructure *ascImages, inputImageStructure *descImages,
						  landSatImage *LSImages, int32_t *autoSize, int32_t writeBlank);
static void mallocOutputImage(outputImageStructure *outputImage);
static void readArgs(int32_t argc, char *argv[], char **inputFile, char **demFile, char **outFileBase, float *fl, char **irregFile,
					 char **shelfMaskFile, double *tieThresh, char **extraTieFile, char **tideFile, int32_t *north, char **landSatFile,
					 int32_t *threeDOffFlag, float *timeThresh, float *timeThreshPhase, char **date1, char **date2,
					 referenceVelocity *refVel, int32_t *statsFlag, outputImageStructure *outputImage, int32_t *writeBlank,
					 char **verticalCorrectionFile, int32_t *GTiff, int32_t *COG);
static void write3Doutput(outputImageStructure outputImage, char *outFileBase);
static void write3DTiffOutput(outputImageStructure outputImage, char *outFileBase, char *driverType, const char *epsg, char *date1,  char *date2);
static int32_t writeMetaFile(inputImageStructure *image, outputImageStructure *outputImage, vhParams *params, char *outFileBase,
							 char *demFile, int32_t writeBlank);
void caldat(int32_t julian, int32_t *mm, int32_t *id, int32_t *iyyy);
static void readReferenceVelMosaic(referenceVelocity *refVel, outputImageStructure *outputImage);
static void processMosaicDate(outputImageStructure *outputImage, char *date1, char *date2);
static void usage();
static void logInputs3d(outputImageStructure *outputImage, char *outFileBase, char *inputFile, char *demFile, char *irregFile,
						char *shelfMaskFile, char *extraTieFile,
						char *tideFile, char *verticalCorrectionFile, float fl, int32_t statsFlag, int32_t threeDOffFlag,
						double tieThresh, referenceVelocity *refVel);
static void logInputFiles3d(outputImageStructure *outputImage, char **geodatFiles, char **phaseFiles, char **baselineFiles,
							char **rOffsetFiles, char **rParamsFiles, char **offsetFiles, char **azParamsFiles, float *nDays, float *weights,
							int32_t *crossFlags, int32_t nFiles, int32_t offsetFlag);
static void get3DProj(inputImageStructure *ascImages, inputImageStructure *descImages, int32_t nAsc, int32_t nDesc,
					  int32_t northFlag, outputImageStructure *outputImage);
static void malloc3DBuffers();
static void init3DImages(outputImageStructure *outputImage, referenceVelocity *refVel);
static void consolodateLists(inputImageStructure **images, vhParams **params, inputImageStructure *ascImages,
							 inputImageStructure *descImages, vhParams *ascParams, vhParams *descParams, int32_t nAsc, int32_t nDesc);
/* moved to mosaic3d.h
   #define MAXOFFBUF 30000000
   #define MAXOFFLENGTH 30000
*/
/*
   Global variables definitions
*/
int32_t RangeSize = RANGESIZE;				/* Range size of complex image */
int32_t AzimuthSize = AZIMUTHSIZE;			/* Azimuth size of complex image */
int32_t BufferSize = BUFFERSIZE;			/* Size of nonoverlap region of the buffer */
int32_t BufferLines = 512;					/* # of lines of nonoverlap in buffer */
int32_t sepAscDesc = TRUE;
int32_t HemiSphere = NORTH;
double Rotation = 45.;
double SLat = -91.;

float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
void *lBuf1, *lBuf2, *lBuf3, *lBuf4;

int32_t llConserveMem = 1234; /* Kluge to maintain backwards compat 9/13/06 */

int main(int argc, char *argv[])
{
	extern char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
	extern void *lBuf1, *lBuf2, *lBuf3, *lBuf4;
	extern void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
	FILE *fp;
	xyDEM dem, verticalCorrection;					  /* Main dem, flow dir dem */
	referenceVelocity refVel;						  /* Reference velocity used for clipping and initializing */
	landSatImage *LSImages;							  /* List of landsat images to include */
	outputImageStructure outputImage;				  /* Output image */
	vhParams *ascParams, *descParams, *params, *tmpP; /* Param lists */
	irregularData *irregDat, *irregTmp;
	inputImageStructure *ascImages, *descImages, *images, *tmp; /* Lists of asc/desc and all  images */
	double tieThresh;											/* Limiting value used for tiepoints */
	char *demFile, *inputFile, *outFileBase;
	char **phaseFiles, **geodatFiles, **baselineFiles;
	float *weights, *nDays;
	float timeThresh, timeThreshPhase;	 /* For mosaicking 3Offsets only use crossing pairs with this many days seperation */
	float fl;							 /* Feathering length */
	char *irregFile;					 /* File with irregularly spaced data for interpolation */
	char **azParamsFiles, **offsetFiles; /* Az paramter and offest files */
	char **rOffsetFiles, **rParamsFiles; /* Range offset file */
	char *extraTieFile, *tideFile;		 /* tide file for tiepoints */
	char *date1, *date2;				 /* Date range */
	char *shelfMaskFile, *landSatFile;	 /* Shelf and landsat file names */
	char *verticalCorrectionFile;
	int32_t i, j;								  /* LCV */
	int32_t northFlag, threeDOffFlag, haveData; /* Flags */
	int32_t nAsc, nDesc, nFiles;			  /* Number of ascending/descending and all files */
	int32_t offsetFlag = TRUE;				  /* Flag to indicate do both offset solution where needed, always true old option removed */
	int32_t statsFlag, *crossFlags;
	int32_t autoSize, writeBlank, count;
	int32_t GTiff, COG;
	const char *epsg=NULL;

	GDALAllRegister();
	/*
	   Read command line args and compute filenames
	*/
	readArgs(argc, argv, &inputFile, &demFile, &outFileBase, &fl, &irregFile, &shelfMaskFile, &tieThresh,
			 &extraTieFile, &tideFile, &northFlag, &landSatFile, &threeDOffFlag,
			 &timeThresh, &timeThreshPhase, &date1, &date2, &refVel, &statsFlag, &outputImage, &writeBlank, &verticalCorrectionFile,
			 &GTiff, &COG);
	/* Added August 2021 to set projection parameters from DEM */
	readXYDEMGeoInfo(demFile, &dem, TRUE);
	/* Removed no offset flag version */
	processMosaicDate(&outputImage, date1, date2);
	/* Write to log */
	logInputs3d(&outputImage, outFileBase, inputFile, demFile, irregFile, shelfMaskFile, extraTieFile, tideFile, verticalCorrectionFile, fl, statsFlag,
				threeDOffFlag, tieThresh, &refVel);
	/*
	  read inputfile
	*/
	getMVhInputFile(inputFile, &phaseFiles, &geodatFiles, &baselineFiles, &offsetFiles, &azParamsFiles, &rOffsetFiles, &rParamsFiles,
					&outputImage, &nDays, &weights, &crossFlags, &nFiles, offsetFlag, outputImage.rOffsetFlag, threeDOffFlag);
	if (outputImage.rOffsetFlag == TRUE)
		fprintf(stderr, "*** USING RANGE/AZIMUTH OFFSETS in solution *** \n");
	fprintf(outputImage.fpLog, ";\n; **** SOURCE DATA **** \n");
	/*
	  Parse inputfiles and set everything up.
	*/
	setup3D(nFiles, phaseFiles, geodatFiles, baselineFiles, offsetFiles, azParamsFiles, rOffsetFiles, rParamsFiles, nDays, weights, crossFlags,
			&ascImages, &descImages, &ascParams, &descParams, &nAsc, &nDesc, offsetFlag, outputImage.rOffsetFlag, threeDOffFlag,
			outputImage.fpLog, &outputImage);
	logInputFiles3d(&outputImage, geodatFiles, phaseFiles, baselineFiles, rOffsetFiles, rParamsFiles, offsetFiles, azParamsFiles,
					nDays, weights, crossFlags, nFiles, offsetFlag);
	
	fprintf(stderr, "%s %s\n", date1, date2); 
	/*
	  Determine hemisphere
	*/
	get3DProj(ascImages, descImages, nAsc, nDesc, northFlag, &outputImage);
	outputImage.slat = dem.stdLat;
	/* Process landsat images */
	LSImages = NULL;
	if (landSatFile != NULL)
	{
		LSImages = parseLSInputs(landSatFile, LSImages, outputImage.jd1, outputImage.jd2, outputImage.timeOverlapFlag);
	}
	/* Find bounding box */
	fprintf(stderr, "nAsc/nDesc %i %i\n", nAsc, nDesc);
	findOutBounds(&outputImage, ascImages, descImages, LSImages, &autoSize, writeBlank);
	/* Remove images that are outside output area */
	removeOutOfBounds(&outputImage, &ascImages, &ascParams, &nAsc);
	removeOutOfBounds(&outputImage, &descImages, &descParams, &nDesc);
	/*
	  Read shelf mask
	*/
	if (shelfMaskFile != NULL)
	{
		fprintf(outputImage.fpLog, "; Starting reading Shelf mask file\n");
		fflush(outputImage.fpLog);
		readShelf(&outputImage, shelfMaskFile);
		fprintf(outputImage.fpLog, "; Finished reading Shelf mask file\n");
		fflush(outputImage.fpLog);
	}
	else
		outputImage.shelfMask = NULL;

	if (verticalCorrectionFile != NULL)
	{
		fprintf(outputImage.fpLog, "; Starting reading vertical correction file\n");
		fflush(outputImage.fpLog);
		readXYDEM(verticalCorrectionFile, &verticalCorrection);
		outputImage.verticalCorrection = &verticalCorrection;
		fprintf(outputImage.fpLog, "; Finished reading vertical correction file\n");
		fflush(outputImage.fpLog);
	}
	else
		outputImage.verticalCorrection = NULL;
	/*
	   read ref vel file for error clipping
	*/
	if (refVel.velFile != NULL && (nAsc + nDesc) > 0)
	{
		readReferenceVelMosaic(&refVel, &outputImage);
	}
	/*
	  Init output image memory
	*/
	malloc3DBuffers();
	mallocOutputImage(&outputImage);
	/* Consolodate lists */
	consolodateLists(&images, &params, ascImages, descImages, ascParams, descParams, nAsc, nDesc);
 	computeDateRange(&date1, &date2, &outputImage, images, params);
	tmpP = params;
	for (tmp = images; tmp != NULL; tmp = tmp->next, tmpP = tmpP->next)
	{
		fprintf(stderr, "%s %s\n", tmpP->offsets.file, tmp->file);
		if (strstr(tmp->file, "nophase") != NULL)
		{
			fprintf(stderr, "Skipping %s :no phase given\n", tmp->file);
			continue;
		}
		fprintf(stderr, "--\n");
	}
	/*
	  Init and input DEM
	*/
	fprintf(outputImage.fpLog, ";\n; About to Read XYDEM %s\n", demFile);
	readXYDEM(demFile, &dem);
	for (tmpP = params; tmpP != NULL; tmpP = tmpP->next)
		tmpP->xydem = dem;
	fprintf(outputImage.fpLog, ";\n; Returned from readXYDEM\n");
	fflush(outputImage.fpLog);
	/*
	  Init values
	*/
	init3DImages(&outputImage, &refVel);

	/*
	  Step 0: Switched to first map since to accomdate discard of large dt.
	  *******************************START Landsat mosaics******************************
	  */
	if (landSatFile != NULL && LSImages != NULL)
	{
		if (statsFlag == TRUE)
			error("Landsat  incompatible with stats flag, which is for speckle tracked offsets only\n");
		makeLandSatMosaic(LSImages, &outputImage, fl);
	}
	/*
	  Step 0: Mosaic using ascending and descending data where possible.
	*/
	if ((nAsc + nDesc) > 0)
	{
		make3DMosaic(images, descImages, params, descParams, &dem, &outputImage, fl, outputImage.no3d, timeThreshPhase);
	}
	/*
	  Step 1:
	*/
	if ((nAsc + nDesc) > 0 && threeDOffFlag == TRUE)
	{
		make3DOffsets(images, params, &dem, &outputImage, fl, timeThresh);
		fprintf(stderr, "End of 3d offsets %i\n", threeDOffFlag);
	}
	/*	  Setup dem stuff	*/

	/*
	   Step 2: Make phase/az offset velocity
	*/
	if (outputImage.noVhFlag == FALSE && (nAsc + nDesc) > 0)
	{
		makeVhMosaic(images, params, &outputImage, fl);
	}
	else
		fprintf(outputImage.fpLog, ";\n; NoVh flag set, not Enterng makeVhMosaic\n;\n");
	/*
	  Step 3: Include fully speckle-tracked data
	*/
	if (outputImage.rOffsetFlag == TRUE && (nAsc + nDesc) > 0)
	{
		speckleTrackMosaic(images, params, &outputImage, fl, &refVel, statsFlag);
	}
	else
		fprintf(outputImage.fpLog, ";\n; rOffset flag False, not Entering speckleTrackMosaic\n;\n");
	/********************************END Landsat mosaics******************************	*/
	/*
	  Step 6: Include irregularly interpolated data, if specified.
	*/
	if (irregFile != NULL)
	{
		fprintf(stderr, "*** incorporating irregularly gridded data from file %s :  ***\n\n", irregFile);
		irregDat = NULL;
		parseIrregFile(irregFile, &irregDat);
		for (irregTmp = irregDat; irregTmp != NULL; irregTmp = irregTmp->next)
		{
			fprintf(stderr, "|%s|\n", irregTmp->file);
			irregTmp->maxLength = 15;
			irregTmp->maxArea = 75.;
		}
		getIrregData(irregDat);
		addIrregData(irregDat, &outputImage, fl);
	}
	/*
	   write meta file
	*/
	haveData = writeMetaFile(images, &outputImage, params, outFileBase, demFile, writeBlank);
	if (landSatFile != NULL)
		haveData = TRUE;
	/*
	  Output result
	*/
	if (outputImage.makeTies == FALSE)
	{
		if (haveData == TRUE || refVel.initMapFlag == TRUE)
		{
			if(COG == FALSE && GTiff == FALSE) {
				write3Doutput(outputImage, outFileBase);
			} 
			else 
			{
				char *driverType;
				if(COG == TRUE) driverType = "COG"; else driverType = "GTiff";
				write3DTiffOutput(outputImage, outFileBase, driverType, epsg, date1, date2);
			}
			
			remove("NoOutput_noDataInRange");
		}
		else
		{
			fp = fopen("NoOutput_noDataInRange", "w");
			error("No output because no data in range, check dates \n");
		}
	}
	else
	{
		writeTieFile(&outputImage, &dem, &verticalCorrection, outFileBase, tieThresh, extraTieFile, tideFile, autoSize);
	}
}

static void computeDateRange(char **date1, char **date2, outputImageStructure *outputImage, inputImageStructure *images, vhParams *params)
{
	double jd1 = outputImage->jd2, jd2 = outputImage->jd1;
	int year, month, day, hour, minute, second;
	if(*date1 == NULL || *date2 == NULL) {
		vhParams *currentParams;
		inputImageStructure *currentImage;
		for (currentImage = images, currentParams = params;
			 currentImage != NULL;
			 currentImage = currentImage->next, currentParams = currentParams->next)
		{
			fprintf(stderr, "jd %lf", currentImage->julDay);
			jd1 = min(jd1, currentImage->julDay);
			jd2 = max(jd2, currentImage->julDay + currentParams->nDays);
		}
	}
	jd_to_date_and_time(jd1, &year, &month, &day, &hour, &minute, &second);
	if(*date1 == NULL) { *date1 = malloc(11); sprintf(*date1, "%4d-%02d-%02d", year, month, day);}
	jd_to_date_and_time(jd2, &year, &month, &day, &hour, &minute, &second);
	if(*date2 == NULL) { *date2 = malloc(11); sprintf(*date2, "%4d-%02d-%02d", year, month, day);}
	fprintf(stderr, "%s %s", *date1, *date2);
}

static void removeOutOfBounds(outputImageStructure *outputImage, inputImageStructure **ascImages, vhParams **ascParams, int32_t *nImages)
{
	/* Go through lists of images and remove ones that are outside the output area */
	int32_t iMin, iMax, jMin, jMax;
	int32_t count, countAdd = 0, countRemove = 0;
	/* Skip empty list */
	if (*ascImages == NULL || *nImages == 0)
		return;
	inputImageStructure *newList, *newListHead, *tmp;
	vhParams *ptmp, *newParams, *newParamsHead;
	/* Init list */
	newList = NULL;
	newParams = NULL;
	newListHead = NULL;
	newParamsHead = NULL;
	count = 0;
	for (tmp = *ascImages, ptmp = *ascParams; tmp != NULL; tmp = tmp->next, ptmp = ptmp->next)
	{
		getRegion(tmp, &iMin, &iMax, &jMin, &jMax, outputImage);
		/* Check in or out of bounds */
		if ((iMin <= iMax && jMin <= jMax))
		{
			if (newList == NULL)
			{
				newList = tmp;
				newListHead = tmp;
				newParams = ptmp;
				newParamsHead = ptmp;
			}
			else
			{
				newList->next = tmp;
				newList = tmp;
				newParams->next = ptmp;
				newParams = ptmp;
			}
			countAdd++;
		}
		else
		{
			*nImages -= 1;
			countRemove++;
		}
	}
	/* Terminate lists */
	if (newList != NULL)
	{
		newList->next = NULL;
		newParams->next = NULL;
	}
	*ascImages = newListHead;
	*ascParams = newParamsHead;
}

static void consolodateLists(inputImageStructure **images, vhParams **params, inputImageStructure *ascImages, inputImageStructure *descImages,
							 vhParams *ascParams, vhParams *descParams, int32_t nAsc, int32_t nDesc)
{
	inputImageStructure *tmp; /* Tmp list for looping */
	vhParams *tmpP;
	fprintf(stderr, "**** nAsc %i, nDesc %i\n", nAsc, nDesc);
	if (nAsc == 0 && nDesc > 0)
	{
		*images = descImages;
		*params = descParams;
	}
	else if (nAsc > 0 && nDesc == 0)
	{
		*images = ascImages;
		*params = ascParams;
	}
	else if (nAsc > 0 && nDesc > 0)
	{
		*images = ascImages;
		/* Skip to end of list */
		for (tmp = ascImages; tmp->next != NULL; tmp = tmp->next)
			;
		/* append */
		tmp->next = descImages;
		*params = ascParams;
		for (tmpP = ascParams; tmpP->next != NULL; tmpP = tmpP->next)
			;
		tmpP->next = descParams;
	}
	else
	{
		*images = NULL;
		*params = NULL;
	}
}

static void init3DImages(outputImageStructure *outputImage, referenceVelocity *refVel)
{
	float **vXimage, **vYimage, **vZimage;
	float **scaleX, **scaleY, **scaleZ;
	float **errorX, **errorY;
	float vxPt, vyPt, exPt, eyPt;
	double x, y;
	int32_t i, j;

	vXimage = (float **)outputImage->image;
	vYimage = (float **)outputImage->image2;
	vZimage = (float **)outputImage->image3;
	errorX = (float **)outputImage->errorX;
	errorY = (float **)outputImage->errorY;
	scaleX = (float **)outputImage->scale;
	scaleY = (float **)outputImage->scale2;
	scaleZ = (float **)outputImage->scale3;
	for (i = 0; i < outputImage->ySize; i++)
	{
		y = (outputImage->originY + i * outputImage->deltaY) * MTOKM;
		for (j = 0; j < outputImage->xSize; j++)
		{
			x = (outputImage->originX + j * outputImage->deltaX) * MTOKM;
			if (refVel->initMapFlag == TRUE)
			{
				/*	fprintf(stderr,"%i %i ",i,j); */
				if (outputImage->timeOverlapFlag == TRUE)
					error("can't use reference map and timeoverlap flag at the same time");
				refVelInterp(x, y, refVel, &vxPt, &vyPt, &exPt, &eyPt);
				if (vxPt > -2e6 && vxPt < 2e6)
				{
					scaleX[i][j] = 1.0 / (exPt * exPt);
					vXimage[i][j] = vxPt;
					scaleY[i][j] = 1.0 / (eyPt * eyPt);
					vYimage[i][j] = vyPt;
					scaleZ[i][j] = 0.0;
					vZimage[i][j] = 0.0;
					errorX[i][j] = exPt * exPt;
					errorY[i][j] = eyPt * eyPt;
				}
			}
			else
			{
				errorX[i][j] = 0.0;
				errorY[i][j] = 0.0;
				scaleX[i][j] = 0.0;
				vXimage[i][j] = -2.e9;
				scaleY[i][j] = 0.0;
				vYimage[i][j] = -2.e9;
				scaleZ[i][j] = 0.0;
				vZimage[i][j] = -2.e9;
			}
		}
	}
}

static void malloc3DBuffers()
{
	extern char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
	extern void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
	extern void *lBuf1, *lBuf2, *lBuf3, *lBuf4;

	Abuf1 = malloc(MAXADBUF);
	Abuf2 = malloc(MAXADBUF2);
	Dbuf1 = malloc(MAXADBUF);
	Dbuf2 = malloc(MAXADBUF2);
	/*
	  Init buffer memory
	*/
	offBufSpace1 = (void *)malloc(MAXOFFBUF);
	offBufSpace2 = (void *)malloc(MAXOFFBUF);
	offBufSpace3 = (void *)malloc(MAXOFFBUF);
	offBufSpace4 = (void *)malloc(MAXOFFBUF);
	/* Used for pointers rows to above buffer space */
	lBuf1 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	lBuf2 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	lBuf3 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	lBuf4 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
}

static void get3DProj(inputImageStructure *ascImages, inputImageStructure *descImages, int32_t nAsc, int32_t nDesc, int32_t northFlag, outputImageStructure *outputImage)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	float lat;

	fprintf(outputImage->fpLog, "; Determining Hemisphere\n;\n");
	if (nAsc != 0)
		lat = ascImages->latControlPoints[1];
	else if (nDesc != 0)
		lat = descImages->latControlPoints[1];
	else
	{
		if (northFlag == FALSE)
		{
			fprintf(stderr, "*** No Asc or Desc images specified - assuming south ***");
			lat = -80.;
		}
		else
		{
			fprintf(outputImage->fpLog, "; *** No Asc or Desc images specified - using north flag ***");
			lat = 80.;
		}
	}
	/* modifed august 2021 to use dem to set projection */
	if (lat < 0)
	{
		HemiSphere = SOUTH; /* Rotation=0.0;  outputImage->slat=71.0;*/
		fprintf(stderr, "\n ***** Southern Hemisphere  ******\n");
		fprintf(outputImage->fpLog, "; ***** Southern Hemisphere  ******\n;\n");
	}
	else
	{
		/* outputImage->slat=70.0;*/
		HemiSphere = NORTH; /* Rotation=0.0;  outputImage->slat=71.0;*/
		fprintf(stderr, "\n ***** Northern Hemisphere  ******\n");
		fprintf(outputImage->fpLog, "; ***** Northern Hemisphere  ******\n;\n");
	}
	fflush(outputImage->fpLog);
}

static void logInputFiles3d(outputImageStructure *outputImage, char **geodatFiles, char **phaseFiles, char **baselineFiles, char **rOffsetFiles,
							char **rParamsFiles, char **offsetFiles, char **azParamsFiles, float *nDays, float *weights, int32_t *crossFlags, int32_t nFiles, int32_t offsetFlag)
{
	int32_t count = 0, i;
	/* loop through and log each file */
	for (i = 0; i < nFiles; i++)
	{
		if (geodatFiles[i] != NULL)
		{
			count++;
			fprintf(outputImage->fpLog, ";\n; Data Take %i\n", count);
			fprintf(outputImage->fpLog, "; phaseFile     : %s\n; Geodat File   : %s\n; Baseline File : %s\n", phaseFiles[i], geodatFiles[i], baselineFiles[i]);
			fprintf(outputImage->fpLog, "; nDays         : %f", nDays[i]);
			fprintf(outputImage->fpLog, "; Weight        : %f\n", weights[i]);
			fprintf(outputImage->fpLog, "; CrossFlag        : %i\n", crossFlags[i]);
			if (offsetFlag == TRUE && offsetFiles[i] != NULL && azParamsFiles[i] != NULL)
				fprintf(outputImage->fpLog, "; Az offsetFile : %s\n; Az Params     : %s\n", offsetFiles[i], azParamsFiles[i]);
			if (outputImage->rOffsetFlag == TRUE && rOffsetFiles[i] != NULL)
				fprintf(outputImage->fpLog, "; Range Offsets : %s\n; Rbaseline     : %s\n", rOffsetFiles[i], rParamsFiles[i]);
			else
				fprintf(outputImage->fpLog, "; Range Offsets : %s\n; Rbaseline     : %s\n", "none", "none");
		}
	}
	fflush(outputImage->fpLog);
}

static void logInputs3d(outputImageStructure *outputImage, char *outFileBase, char *inputFile, char *demFile, char *irregFile, char *shelfMaskFile,
						char *extraTieFile, char *tideFile, char *verticalCorrectionFile, float fl, int32_t statsFlag,
						int32_t threeDOffFlag, double tieThresh, referenceVelocity *refVel)
{
	char *logFile;
	logFile = (char *)malloc(1024);
	logFile[0] = '\0';
	logFile = strcat(logFile, outFileBase);
	logFile = strcat(logFile, ".log");
	outputImage->fpLog = fopen(logFile, "w");
	if (outputImage->fpLog == NULL)
		error("Could not open log file %s", logFile);

	fprintf(outputImage->fpLog, "; Mosaic3d log file\n");
	fprintf(outputImage->fpLog, ";\n; ***Command Line Parameters***\n;\n");
	fprintf(outputImage->fpLog, "; inputFile        : %s\n", inputFile);
	fprintf(outputImage->fpLog, "; demFile          : %s\n", demFile);
	if (irregFile != NULL)
		fprintf(outputImage->fpLog, "; irregFile        : %s\n", irregFile);
	else
		fprintf(outputImage->fpLog, "; irregFile        : %s\n", "None");
	if (shelfMaskFile != NULL)
		fprintf(outputImage->fpLog, "; shelfMaskFile    : %s\n", shelfMaskFile);
	else
		fprintf(outputImage->fpLog, "; shelfMaskFile        : %s\n", "None");
	if (outputImage->makeTies == TRUE && extraTieFile != NULL)
		fprintf(outputImage->fpLog, "; extraTieFile     : %s\n", extraTieFile);
	else
		fprintf(outputImage->fpLog, "; extraTieFile     : %s\n", "None");
	if (outputImage->makeTies == TRUE && tideFile != NULL)
		fprintf(outputImage->fpLog, "; tideFile         : %s\n", tideFile);
	if (verticalCorrectionFile != NULL)
		fprintf(outputImage->fpLog, "; verticalCorrection         : %s\n", verticalCorrectionFile);
	else
		fprintf(outputImage->fpLog, "; verticalCorrection         : %s\n", "None");
	fprintf(outputImage->fpLog, "; OutputFile base  : %s\n", outFileBase);
	fprintf(outputImage->fpLog, "; Feather length   : %f\n", fl);
	fprintf(outputImage->fpLog, "; NoVh Flag        : %i\n", outputImage->noVhFlag);
	fprintf(outputImage->fpLog, "; rOffset Flag     : %i\n", outputImage->rOffsetFlag);
	fprintf(outputImage->fpLog, "; No3D Flag        : %i\n", outputImage->no3d);
	fprintf(outputImage->fpLog, "; Stats Flag        : %i\n", statsFlag);
	fprintf(outputImage->fpLog, "; Maketies Flag    : %i\n", outputImage->makeTies);
	fprintf(outputImage->fpLog, "; No Tide  Flag    : %i\n", outputImage->noTide);
	fprintf(outputImage->fpLog, "; ThreeD Off Flag    : %i\n", threeDOffFlag);
	fprintf(outputImage->fpLog, "; TimeOverlapFlag    : %i\n", outputImage->timeOverlapFlag);
	fprintf(outputImage->fpLog, "; DeltaB    : %i\n", outputImage->deltaB);
	if (refVel->velFile != NULL)
	{
		fprintf(outputImage->fpLog, "; Reference Velocity    : %s\n", refVel->velFile);
		fprintf(outputImage->fpLog, "; ClipFlag    : %f\n", refVel->clipThresh);
		fprintf(outputImage->fpLog, "Include Ref Vel in mosaic %i\n", refVel->initMapFlag);
		fprintf(stderr, "Include Ref Vel in mosaic %i\n", refVel->initMapFlag);
		fprintf(outputImage->fpLog, "Compute errors mode using ref Vel %i\n", refVel->initMapFlag);
	}
	if (outputImage->makeTies == TRUE)
		fprintf(outputImage->fpLog, "; tieThresh        : %lf\n", tieThresh);
	fflush(outputImage->fpLog);
}

static void toSigma(outputImageStructure outputImage)
{
	// convert error variance to sigma
	for (int i = 0; i < outputImage.ySize; i++)
	{
		for (int j = 0; j < outputImage.xSize; j++)
		{
			if (outputImage.errorX[i][j] > 0.0)
			{
				outputImage.errorX[i][j] = (float)sqrt((double)outputImage.errorX[i][j]);
				outputImage.errorY[i][j] = (float)sqrt((double)outputImage.errorY[i][j]);
			}
		}
	}
}

static void filterDT(outputImageStructure outputImage)
{
	/* Added 0.49 to ensure last day of month gets included. */
		double maxdT = (outputImage.jd2 - outputImage.jd1 + 1.0) * 0.5 + 0.49; /* Remove .49 ? */
		fprintf(stderr, "maxdT %f %i\n", maxdT, outputImage.timeOverlapFlag);
		if (maxdT < 0 || maxdT > 20000)
			fprintf(stderr, "invalid dT (<0 or > 20000) in writing ouput");
		/* recast pointers */
		float **vx = (float **)outputImage.image;
		float **vy = (float **)outputImage.image2;
		float **dT = (float **)outputImage.image3;
		for (int i = 0; i < outputImage.ySize; i++)
		{
			for (int j = 0; j < outputImage.xSize; j++)
			{
				/* Ensure only data in date range written */
				if (dT[i][j] > maxdT || dT[i][j] < (-maxdT))
				{
					vx[i][j] = -2.e9;
					vy[i][j] = -2.e9;
					dT[i][j] = -2.e9;
					outputImage.errorX[i][j] = -2.e9;
					outputImage.errorY[i][j] = -2.e9;
				}
			}
		}
}

static void write3Doutput(outputImageStructure outputImage, char *outFileBase)
{
	char *outFileVx, *outFileVy, *outFileVz; /* Output files */
	char *outFileEx, *outFileEy;
	/* Output file names velocity */
	outFileVx = appendSuffix(outFileBase, ".vx", (char *)malloc(strlen(outFileBase) + 4));
	outFileVy = appendSuffix(outFileBase, ".vy", (char *)malloc(strlen(outFileBase) + 4));
	if (outputImage.timeOverlapFlag == FALSE)
	{
		outFileVz = appendSuffix(outFileBase, ".vz", (char *)malloc(strlen(outFileBase) + 4));
	}
	else
	{
		outFileVz = appendSuffix(outFileBase, ".dT", (char *)malloc(strlen(outFileBase) + 4));
	}
	outFileEx = appendSuffix(outFileBase, ".ex", (char *)malloc(strlen(outFileBase) + 4));
	outFileEy = appendSuffix(outFileBase, ".ey", (char *)malloc(strlen(outFileBase) + 4));
	// convert error variance to sigma
	toSigma(outputImage);
	//  Fileter by dT if timeOverlap
	if (outputImage.timeOverlapFlag == TRUE)
	{
		filterDT(outputImage);
	}
	outputGeocodedImage(outputImage, outFileVx);
	free(outputImage.image[0]);
	free(outputImage.image);
	outputImage.image = outputImage.image2;
	outputGeocodedImage(outputImage, outFileVy);
	free(outputImage.image[0]);
	free(outputImage.image);
	outputImage.image = outputImage.image3;
	outputGeocodedImage(outputImage, outFileVz);
	free(outputImage.image[0]);
	free(outputImage.image);
	outputImage.image = (void **)outputImage.errorX;
	outputGeocodedImage(outputImage, outFileEx);
	free(outputImage.image[0]);
	free(outputImage.image);
	outputImage.image = (void **)outputImage.errorY;
	outputGeocodedImage(outputImage, outFileEy);
	free(outputImage.image[0]);
	free(outputImage.image);
}



static void write3DTiffOutput(outputImageStructure outputImage, char *outFileBase, char *driverType, const char *epsg, char *date1,  char *date2)
{
	char *outFileVx, *outFileVy, *outFileVz; /* Output files */
	char *outFileEx, *outFileEy;
	double geoTransform[6];
	float noDataValues[5] = {-2.0e9, -2.0e9, -2.0e9, -2.0e9, -2.0e9};
	/* Output file names velocity */
	outFileVx = appendSuffix(outFileBase, ".vx.tif", (char *)malloc(strlen(outFileBase) + 8));
	outFileVy = appendSuffix(outFileBase, ".vy.tif", (char *)malloc(strlen(outFileBase) + 8));
	if (outputImage.timeOverlapFlag == FALSE)
	{
		outFileVz = appendSuffix(outFileBase, ".vz.tif", (char *)malloc(strlen(outFileBase) + 8));
	}
	else
	{
		outFileVz = appendSuffix(outFileBase, ".dT.tif", (char *)malloc(strlen(outFileBase) + 8));
	}
	outFileEx = appendSuffix(outFileBase, ".ex.tif", (char *)malloc(strlen(outFileBase) + 8));
	outFileEy = appendSuffix(outFileBase, ".ey.tif", (char *)malloc(strlen(outFileBase) + 8));
	// convert error variance to sigma
	toSigma(outputImage);
	//  Fileter by dT if timeOverlap
	if (outputImage.timeOverlapFlag == TRUE)
	{
		filterDT(outputImage);
	}
    // Compute Geotransform
	computeGeoTransform(geoTransform, outputImage.originX, outputImage.originY, outputImage.xSize,
					    outputImage.ySize, outputImage.deltaX, outputImage.deltaY);
	// Set up meta data
	dictNode *summaryMetaData = NULL;
	char *timeStamp = timeStampMeta();
	insert_node(&summaryMetaData, "FirstDate", date1);
	insert_node(&summaryMetaData, "LastDate", date2);
	insert_node(&summaryMetaData, "CreationTime", timeStamp);
	// Get epsg code
	if(epsg == NULL)
	{
		epsg = getEPSGFromProjectionParams(Rotation, SLat, HemiSphere); 
	}	
	// Save files as tiffs
	// Vx
	saveAsGeotiff(outFileVx, (float *)outputImage.image[0], outputImage.xSize,
				outputImage.ySize, geoTransform, epsg, summaryMetaData, driverType, GDT_Float32, noDataValues[0]);
	free(outputImage.image[0]);
	free(outputImage.image);
	// Vy
	saveAsGeotiff(outFileVy, (float *)outputImage.image2[0], outputImage.xSize,
				 outputImage.ySize, geoTransform, epsg, summaryMetaData, driverType, GDT_Float32, noDataValues[1]);
	free(outputImage.image2[0]);
	free(outputImage.image2);
	// Vz
	saveAsGeotiff(outFileVz, (float *)outputImage.image3[0], outputImage.xSize,
				 outputImage.ySize, geoTransform, epsg, summaryMetaData, driverType, GDT_Float32, noDataValues[2]);
	free(outputImage.image3[0]);
	free(outputImage.image3);
	// Ex
	saveAsGeotiff(outFileEx, (float *)outputImage.errorX[0], outputImage.xSize,
				 outputImage.ySize, geoTransform, epsg, summaryMetaData, driverType, GDT_Float32, noDataValues[3]);
	free(outputImage.errorX[0]);
	free(outputImage.errorX);
	// Ey
	saveAsGeotiff(outFileEy, (float *)outputImage.errorY[0], outputImage.xSize,
				 outputImage.ySize, geoTransform, epsg, summaryMetaData, driverType, GDT_Float32, noDataValues[4]);
	free(outputImage.errorY[0]);
	free(outputImage.errorY);	
	//
	// Create VRT
	char *vrtFile = appendSuffix(outFileBase, ".vrt", (char *)malloc(strlen(outFileBase) + 5));
	const char *bands[] = {outFileVx, outFileVy, outFileVz, outFileEx, outFileEy};
	makeTiffVRT(vrtFile, bands, 5, noDataValues, summaryMetaData);
}

static void processMosaicDate(outputImageStructure *outputImage, char *date1, char *date2)
{
	int32_t m1, m2, d1, d2, y1, y2;

	fprintf(stderr, "date1,date2 %s %s\n", date1, date2);
	if (date1 != NULL)
	{
		if (sscanf(date1, "%2d-%2d-%4d", &m1, &d1, &y1) != 3)
			error("invalid date %s\n", date1);
		outputImage->jd1 = juldayDouble(m1, d1, y1);
		fprintf(stderr, "%i %i %i %lf\n", m1, d1, y1, outputImage->jd1);
	}
	else
		outputImage->jd1 = 0.0;

	if (date2 != NULL)
	{
		if (sscanf(date2, "%2d-%2d-%4d", &m2, &d2, &y2) != 3)
			error("invalid date %s\n", date2);
		outputImage->jd2 = (double)juldayDouble(m2, d2, y2) + 0.9999999; /* add 0.999 to ensure end of day */
		fprintf(stderr, "%i %i %i %lf\n", m2, d2, y2, outputImage->jd2);
	}
	else
		outputImage->jd2 = HIGHJD + 1000.; /* way past present */

	if ((date1 != NULL || date2 != NULL) && (outputImage->jd2 < outputImage->jd1))
		error("date 1 follows date2\n");
	fprintf(stderr, "date1,date2 %f %f\n", outputImage->jd1, outputImage->jd2);
}

static int32_t writeMetaFile(inputImageStructure *image, outputImageStructure *outputImage, vhParams *params, char *outFileBase, char *demFile, int32_t writeBlank)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	FILE *fpMeta;
	char *metaFile;
	char sensor[128];
	char *tmp;
	int32_t mm, day, year;
	double dday;
	time_t tproc;
	struct tm *tS, tTmp;
	double x, y, lat, lon;
	size_t sMax = 128;
	char *timeString;
	int32_t hour, minute, second;
	inputImageStructure *imageTmp;
	int32_t haveData;
	/*
	  If single pair result, make a meta file with the date
	*/
	metaFile = (char *)malloc(LINEMAX);
	metaFile[0] = '\0';
	metaFile = strcat(metaFile, outFileBase);
	metaFile = strcat(metaFile, ".meta");
	fpMeta = fopen(metaFile, "w");
	timeString = (char *)malloc(sizeof(char) * 128);
	haveData = writeBlank; /* If writeBlank this will force have data to true */
	for (imageTmp = image; imageTmp != NULL; imageTmp = imageTmp->next)
	{
		if (imageTmp->used == TRUE)
		{ /* Only output images that were used */
			haveData = TRUE;
			/*
				Julian date for central estimate
			*/
			fprintf(fpMeta, "Central Julian Date (CE) for Pair = %10.3lf\n", imageTmp->julDay + params->nDays / 2);
			/*
			  Date for first image
			*/
			caldat((int32_t)(imageTmp->julDay), &mm, &day, &year);
			fprintf(stderr, "mm,day,year %i %i %i %f\n", mm, day, year, imageTmp->julDay);
			/* note uses month as 0 to 11 */
			tTmp.tm_year = year - 1900;
			tTmp.tm_mon = mm - 1;
			tTmp.tm_mday = day;
			strftime(timeString, sMax, "%b:%d:%Y", &tTmp);
			fprintf(fpMeta, "First Image Date (MM:DD:YYYY) = %s\n", timeString);
			/*
			  Date for first second image
			*/
			caldat((int32_t)(imageTmp->julDay + params->nDays), &mm, &day, &year);
			tTmp.tm_year = year - 1900;
			tTmp.tm_mon = mm - 1;
			tTmp.tm_mday = day;
			strftime(timeString, sMax, "%b:%d:%Y", &tTmp);
			fprintf(fpMeta, "Second Image Date (MM:DD:YYYY) = %s\n", timeString);
			/*
			  Nominal time
			*/
			dday = (imageTmp->julDay - (int32_t)imageTmp->julDay);
			hour = (int)(dday * 24);
			minute = (int)((dday * 24 - hour) * 60);
			second = (int)(dday * 86400 - hour * 3600 - minute * 60 + .5);
			tTmp.tm_sec = second;
			tTmp.tm_hour = hour;
			tTmp.tm_min = minute;
			strftime(timeString, sMax, "%H:%M:%S", &tTmp);
			fprintf(fpMeta, "Nominal Time for Pair (HH:MM:SS) = %s\n", timeString);
			/*
			  Sensor
			*/
			sensor[0] = '\0';
			tmp = imageTmp->par.label;
			/*
			  Extra sensor from geodat title. Only implemented now for TSX/TDX
			*/
			if (strstr(tmp, "TDX") != NULL)
				strcat(sensor, "TSX/TDX");
			else if (strstr(imageTmp->par.label, "TSX") != NULL)
				strcat(sensor, "TSX/TDX");
			else
				strcat(sensor, "not_specfied");
			fprintf(fpMeta, "Sensor = %s\n", sensor);
		}

	} /* End for image loop */
	x = outputImage->originX + 0.5 * ((outputImage->xSize - 1) * outputImage->deltaX);
	x *= 0.001;
	y = outputImage->originY + 0.5 * ((outputImage->ySize - 1) * outputImage->deltaY);
	y *= 0.001;
	xytoll1(x, y, HemiSphere, &lat, &lon, Rotation, outputImage->slat);
	if (lon > 180)
		lon -= 360;
	fprintf(fpMeta, "Product Center Latitude  = %10.5lf\n", lat);
	fprintf(fpMeta, "Product Center Longitude = %10.5lf\n", lon);
	/*
		DEM
	*/
	if (strstr(demFile, "gimp1") != NULL)
	{
		fprintf(fpMeta, "DEM version = GIMP DEM V1\n");
	}
	else
	{
		if (strstr(demFile, "gimp2") != NULL)
			fprintf(fpMeta, "DEM version = GIMP DEM V2\n");
		else
			fprintf(fpMeta, "DEM version = Custom\n");
	}
	/*
	  time
	*/
	tproc = time(&tproc);
	tS = localtime(&tproc);
	strftime(timeString, sMax, "%b-%d-%Y-%H:%M:%S\n", tS);
	fprintf(fpMeta, "Production Date/Time = %s\n", timeString);
	if (haveData == FALSE)
		remove(metaFile);
	return (haveData);
}

static void readArgs(int32_t argc, char *argv[], char **inputFile, char **demFile, char **outFileBase, float *fl, char **irregFile, char **shelfMaskFile,
					 double *tieThresh, char **extraTieFile, char **tideFile, int32_t *north, char **landSatFile, int32_t *threeDOffFlag, float *timeThresh, float *timeThreshPhase,
					 char **date1, char **date2, referenceVelocity *refVel, int32_t *statsFlag, outputImageStructure *outputImage, int32_t *writeBlank, char **verticalCorrectionFile,
					int32_t *GTiff, int32_t *COG)
{
	extern int32_t sepAscDesc;
	char *argString;
	float tmp;
	int32_t i, n;
	int32_t noVhFlag, no3d, rOffsetFlag, vzFlag, noTide, timeOverlapFlag;
	int32_t deltaB;

	if (argc < 4 || argc > 30)
	{
		fprintf(stderr, "Arg count > 27: %i\n", argc);
		usage(); /* Check number of args */
	}
	n = argc - 4;
	/*
	  Defaults
	*/
	refVel->clipFlag = FALSE;
	refVel->clipThresh = 100000;
	refVel->initMapFlag = FALSE;
	refVel->velFile = NULL;
	*date1 = NULL;
	*date2 = NULL;
	noTide = FALSE;
	outputImage->makeTies = FALSE;
	*irregFile = NULL;
	rOffsetFlag = FALSE;
	*shelfMaskFile = NULL;
	*extraTieFile = NULL;
	*landSatFile = NULL;
	noVhFlag = FALSE;
	timeOverlapFlag = FALSE;
	no3d = FALSE;
	*threeDOffFlag = FALSE;
	*tieThresh = 100.0; /* Limit on velocity for tie points */
	*fl = 0.0;
	*tideFile = NULL;
	*north = FALSE;
	*timeThresh = 12;
	*timeThreshPhase = 548; /* Allow pairing with in 1.5 years */
	*writeBlank = FALSE;
	vzFlag = VZDEFAULT;
	*statsFlag = FALSE;
	*verticalCorrectionFile = NULL;
	*COG = FALSE;
	*GTiff = FALSE;
	/* Added this flag to sort ignore crossing orbits or like asc/desc types May 6 2014 */
	sepAscDesc = TRUE;
	deltaB = DELTABNONE;
	/*
	  Parse command line args
	*/
	for (i = 1; i <= n; i++)
	{
		argString = strchr(argv[i], '-');
		if (strstr(argString, "center") != NULL)
			fprintf(stderr, "Ignoring obsolete center flag\n");
		else if (strstr(argString, "xyDEM") != NULL)
			fprintf(stderr, "xyDEM flag obsolete - xydem is the default");
		else if (strstr(argString, "writeBlank") != NULL)
			*writeBlank = TRUE;
		else if (strstr(argString, "rOffsets") != NULL)
			rOffsetFlag = TRUE;
		else if (strstr(argString, "offsets") != NULL)
			fprintf(stderr, "ignoring obsolet offsets flag - always enabled\n");
		else if (strstr(argString, "-initMap") != NULL)
			refVel->initMapFlag = TRUE;
		else if (strstr(argString, "makeTies") != NULL)
			outputImage->makeTies = TRUE;
		else if (strstr(argString, "timeOverlap") != NULL)
			timeOverlapFlag = TRUE;
		else if (strstr(argString, "stats") != NULL)
			*statsFlag = TRUE;
		else if (strstr(argString, "COG") != NULL)
			*COG = TRUE;
		else if (strstr(argString, "GTiff") != NULL)
			*GTiff = TRUE;
		else if (strstr(argString, "vzFlag") != NULL)
		{
			if (sscanf(argv[i + 1], "%i\n", &vzFlag) != 1)
				usage();
			i++;
		}
		else if (strstr(argString, "tieThresh") != NULL)
		{
			if (sscanf(argv[i + 1], "%lf\n", tieThresh) != 1)
				usage();
			i++;
		}
		else if (strstr(argString, "extraTies") != NULL)
		{
			*extraTieFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "tideFile") != NULL)
		{
			*tideFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "verticalCorrection") != NULL)
		{
			*verticalCorrectionFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "fl") != NULL)
		{
			sscanf(argv[i + 1], "%f", fl);
			i++;
		}
		else if (strstr(argString, "timeThresh") != NULL)
		{
			sscanf(argv[i + 1], "%f", timeThresh);
			i++;
		}
		else if (strstr(argString, "timePhaseThresh") != NULL)
		{
			sscanf(argv[i + 1], "%f", timeThreshPhase);
			i++;
		}
		else if (strstr(argString, "irreg") != NULL)
		{
			*irregFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "noVh") != NULL)
		{
			noVhFlag = TRUE;
		}
		else if (strstr(argString, "no3d") != NULL)
		{
			no3d = TRUE;
		}
		else if (strstr(argString, "3dOff") != NULL)
		{
			*threeDOffFlag = TRUE;
		}
		else if (strstr(argString, "noSepAscDesc") != NULL)
		{
			sepAscDesc = FALSE;
		}
		else if (strstr(argString, "noTide") != NULL)
		{
			noTide = TRUE;
		}
		else if (strstr(argString, "SVAlongTrack") != NULL)
		{
			deltaB = DELTABQUAD;
		}
		else if (strstr(argString, "SVConst") != NULL)
		{
			deltaB = DELTABCONST;
		}
		else if (strstr(argString, "north") != NULL)
		{
			*north = TRUE;
		}
		else if (strstr(argString, "landSat") != NULL)
		{
			*landSatFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "refVel") != NULL)
		{
			refVel->velFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "clipThresh") != NULL)
		{
			sscanf(argv[i + 1], "%f", &tmp);
			i++;
			refVel->clipFlag = TRUE;
			refVel->clipThresh = tmp;
		}
		else if (strstr(argString, "shelfMask") != NULL)
		{
			*shelfMaskFile = argv[i + 1];
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
		else
			usage();
	}
	if (refVel->initMapFlag == TRUE && refVel->velFile == NULL)
		error("initMap set but no reference velocity provided\n");
	if (*statsFlag == TRUE)
	{
		fprintf(stderr, "setting flags to force speckletracking only since in computeError or stats mode");
		no3d = TRUE;
		*threeDOffFlag = FALSE;
		noVhFlag = TRUE;
		rOffsetFlag = TRUE;
	}
	if (*statsFlag == TRUE)
	{
		printf("stats and timeOverlap flags incompatible, setting timeOverlap flag to False");
		timeOverlapFlag = FALSE;
	}
	if(*COG == TRUE && *GTiff == TRUE)
	{
		error("Select COG or GTiff but not both");
	}
	/* Must have az offsets to do range offsets */
	*inputFile = argv[argc - 3];
	*demFile = argv[argc - 2];
	*outFileBase = argv[argc - 1];
	outputImage->noVhFlag = noVhFlag;
	outputImage->no3d = no3d;
	outputImage->rOffsetFlag = rOffsetFlag;
	outputImage->noTide = noTide;
	outputImage->vzFlag = vzFlag;
	outputImage->deltaB = deltaB;
	outputImage->timeOverlapFlag = timeOverlapFlag;
	return;
}

static void usage()
{
	error("\033[1m\n\n%s\n\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n\n",
		  "mosaic3d: mosaic phase and speckle data to create a velocity mosaic",
		  "Usage:",
		  " mosaic3d -north -GTiff -COG -writeBlank -makeTies -tieThresh -extraTies extraTieFile -date1 MM-DD-YYYY -date2 MM-DD-YYYY -timeOverlap -tideFile tideFile \\",
		  " \t-verticalCorrection vcFile -no3d -3dOff -deltaBQ -deltaBC -noVh -landSat landSatList  -refVel refVelFile -initMap -clipThresh clipThresh -lsClip clipVal  \\",
		  " \t-shelfMask shelfMask -fl fl  -rOffsets -timeThresh timeThresh -irregFile irregFile -vzFlag flag -stats -noSepAscDesc -noTide \\",
		  " \tinputFile demFile outPutImage\n",
		  " where :",
		  "\tnorth =\t\t\t Force northern hemisphere",
		  "\tGTiff =\t\t\t Save to geotiff files",
		  "\tCOG =\t\t\t Save to Cloud Optimized Geotiff",
		  "\twriteBlank =\t\t\t Force writing of results with no data",
		  "\tmakeTies =\t\t Flag to produce tiepoint file instead of usual output ",
		  "\ttieThresh =\t\t Only use values less than tieThresh for tiepoints",
		  "\textraTies =\t\t File with additional ties to include with data file",
		  "\tdate1,date2 =\t\t Time range for data (include only interval in range (date1,date2). format -> MM-DD-YYYY",
		  "\ttimeOverlap =\t\t Default is False so only images full in (date1,date2) are included. \n\t\t\t\t If set,  images that overlap (date1,date2) are included and weighted accordingly",
		  "\ttideFile =\t\t Tidal correction used for creating tiepoints",
		  "\tverticalCorrection =\t\t Vertical velocity correction (sub/emerg vel) in m/yr",
		  "\tno3d =\t\t\t No crossing orbit solution with phase",
		  "\t3dOff =\t\t\t Crossing orbit solution with offsets",
		  "\tSVAlongTrack =\t\tIf available, use along track correction to  state vector solution",
		  "\tSVConst =\t\t If available, use const only correction on top of state vector solution",
		  "\tnoVh =\t\t\t Use only the crossing orbit solution",
		  "\tlandSat =\t\t File containing list of landsat offsets to include in mosaic",
		  "\trefVel =\t\t Velocity mosaic used to clip errors ",
		  "\tinitMap =\t\t Interpolate refVel as startin point for mosaic ",
		  "\tclipThresh =\t\t Use with refVel to clip differences > clipThresh for slow moving regions (< 100 m/yr) ",
		  "\tlsClip (not implemented)=\t\t Clip landsat mosaic using refVelFile ",
		  "\tshelfMask =\t\t Shelf mask file (includes nodata)",
		  "\tfl =\t\t\t Feather length",
		  "\toffsets =\t\t Use offset azimuth offset data for second component ",
		  "\trOffsets =\t\t Use offset data for both components where needed ",
		  "\ttimeThresh =\t\t For 3D crossing offset solution, only use pairs within timeThresh days (def=12)",
		  "\ttimePhaseThresh =\t\t For 3D crossing phase solution, only use pairs within timePhaseThresh days (def=548 days: 1.5 years)",
		  "\tirregFile =\t\t File with list of irregular input data",
		  "\tvzFlag =\t\t Used for changing output in vertical channel 0 for default (vz correction), 1 horizontal 1/sin(psi), 2 for 1/vertical cos(psi), 3 flag for LOS scaled to m/yr, 4 inc angle",
		  "\tstats =\t\t compute unweight mean vx, and vy and standard dev of the inputs (ex,ey) and number of points (vz) - works only for speckleTrack",
		  "\tnoSepAscDesc =\t\t For 3d ignore asc/desc and use heading, default use asc/desc",
		  "\tnoTide =\t\t Don't compute value if on shelf",
		  "\tinputFile =\t\t File with input params, dem and geodat filenames",
		  "\tdemFile =\t\t File with nonInsar dem",
		  "\toutputImage =\t\t Root of output image (e.g., mosaicOffsets)\033[0m");
}

/*
  Estimate output area for image based on location of inputs
*/
static void findOutBounds(outputImageStructure *outputImage, inputImageStructure *ascImages, inputImageStructure *descImages, landSatImage *LSImages, int32_t *autoSize, int32_t writeBlank)
{
	double minX, maxX, minY, maxY, x, y;
	double minXLS, maxXLS, minYLS, maxYLS;
	inputImageStructure *tmp; /* Tmp list for looping */
	landSatImage *LStmp;
	int32_t j;
	/*
	  No data case for tiepoints
	*/
	if (ascImages == NULL && descImages == NULL && LSImages == NULL && writeBlank == FALSE)
	{
		outputImage->xSize = (int)0;
		outputImage->ySize = (int)0;
		outputImage->originX = 0;
		outputImage->originY = 0;
		outputImage->deltaX = 1.;
		outputImage->deltaY = 1.;
	}
	/* If ouput size zero, auto size */
	if (outputImage->xSize == 0 || outputImage->ySize == 0)
	{
		minXLS = (double)LARGEINT;
		maxXLS = -(double)LARGEINT;
		minYLS = (double)LARGEINT;
		maxYLS = -(double)LARGEINT;
		if (LSImages != NULL)
		{
			for (LStmp = LSImages; LStmp != NULL; LStmp = LStmp->next)
			{
				minXLS = min(minXLS, LStmp->matches.x0 * MTOKM);
				maxXLS = max(maxXLS, (LStmp->matches.x0 + LStmp->matches.stepX * LStmp->matches.dx * (LStmp->matches.nx - 1)) * MTOKM);
				minYLS = min(minYLS, LStmp->matches.y0 * MTOKM);
				maxYLS = max(maxYLS, (LStmp->matches.y0 + LStmp->matches.stepY * LStmp->matches.dy * (LStmp->matches.ny - 1)) * MTOKM);
			}
			fprintf(stderr, "LS Bounds  %lf %lf %lf %lf\n", minXLS, maxXLS, minYLS, maxYLS);
		}
		/*
		  Find output bounds
		*/
		fprintf(outputImage->fpLog, ";\n; Entering findOutBounds (mosaic3d.c)\n;\n");
		minX = 1.e30;
		minY = 1.0e30;
		maxX = -1.e30;
		maxY = -1.0e30;
		for (tmp = ascImages; tmp != NULL; tmp = tmp->next)
		{
			fprintf(stderr, "A %f %f %f %f\n", tmp->minX, tmp->minY, tmp->maxX, tmp->maxY);
			minX = min(tmp->minX, minX);
			minY = min(tmp->minY, minY);
			maxX = max(tmp->maxX, maxX);
			maxY = max(tmp->maxY, maxY);
			/*
			for(j=1; j < 5; j++) {
				lltoxy1(tmp->latControlPoints[j],tmp->lonControlPoints[j],&x,&y,
					Rotation,outputImage->slat);
				minX=min(x,minX);minY=min(y,minY);maxX=max(x,maxX); maxY=max(y,maxY);
			}*/
		}
		for (tmp = descImages; tmp != NULL; tmp = tmp->next)
		{
			fprintf(stderr, "D %f %f %f %f\n", tmp->minX, tmp->minY, tmp->maxX, tmp->maxY);
			minX = min(tmp->minX, minX);
			minY = min(tmp->minY, minY);
			maxX = max(tmp->maxX, maxX);
			maxY = max(tmp->maxY, maxY);
			/*
			for(j=1; j < 5; j++) {
				lltoxy1(tmp->latControlPoints[j],tmp->lonControlPoints[j],&x,&y,
					Rotation,outputImage->slat);
				minX=min(x,minX);minY=min(y,minY);maxX=max(x,maxX); maxY=max(y,maxY);
				fprintf(stderr,"%f %f\n", x,y);
				fprintf(stderr,"%f %f %f %f ..\n", tmp->minX, tmp->maxX, tmp->minY, tmp->maxY);
			} */
		}
		if (LSImages != NULL)
		{
			minX = min(minX, minXLS);
			maxX = max(maxX, maxXLS);
			minY = min(minY, minYLS);
			maxY = max(maxY, maxYLS);
		}
		minX = (double)((int32_t)minX - 3);
		maxX = (double)((int32_t)maxX + 3);
		minY = (double)((int32_t)minY - 3);
		maxY = (double)((int32_t)maxY + 3);
		fprintf(stderr, "minX,maxX,minY,maxY %f %f %f %f %i %i", minX, maxX, minY, maxY, outputImage->xSize, outputImage->ySize);
		*autoSize = TRUE;
		fprintf(stderr, "\n\n*** AUTOSIZING REGION *** *\n\n");
		fprintf(outputImage->fpLog, ";   *** AUTOSIZING REGION ***\n;\n");
		fprintf(stderr, "minX,maxX,minY,maxY %f %f %f %f %i %i %f\n", minX, maxX, minY, maxY, outputImage->xSize, outputImage->ySize, outputImage->deltaX);
		outputImage->xSize = (int)((maxX - minX) / (outputImage->deltaX * MTOKM) + 0.5);
		outputImage->ySize = (int)((maxY - minY) / (outputImage->deltaY * MTOKM) + 0.5);
		fprintf(stderr, "minX,maxX,minY,maxY %f %f %f %f %i %i\n", minX, maxX, minY, maxY, outputImage->xSize, outputImage->ySize);
		outputImage->originX = minX * KMTOM;
		outputImage->originY = minY * KMTOM;
	}
	else
	{
		fprintf(stderr, "\n;   *** USER SPECIFIED SIZE ***\n;\n");
		*autoSize = FALSE;
	}

	fprintf(stderr, "dimensions %i %i %f %f %f %f\n\n", outputImage->xSize, outputImage->ySize,
			outputImage->originX, outputImage->originY, outputImage->deltaX, outputImage->deltaY);
	fprintf(outputImage->fpLog, "; Size       : %i %i\n; Origin     : %f %f \n; Spacing    : %f %f\n;\n",
			outputImage->xSize, outputImage->ySize, outputImage->originX, outputImage->originY, outputImage->deltaX, outputImage->deltaY);
	fprintf(outputImage->fpLog, ";\n; Leaving findOutBounds (mosaic3d.c)\n;\n");
	fflush(outputImage->fpLog);
}

static void mallocOutputImage(outputImageStructure *outputImage)
{
	float *buf1, *buf2, *buf3, *buf1s, *buf2s, *buf3s;
	float *bufex, *bufey;
	float **vXimage, **vYimage, **vZimage;
	float **scaleX, **scaleY, **scaleZ;
	int32_t bufSize, i, j, k;
	float *bufx, *bufy, *bufz, *bufs, *bufsx, *bufsy;
	fprintf(outputImage->fpLog, ";\n; Entering  mallocOutputImage (from mosaic3d)\n");
	outputImage->imageType = POWER;
	bufSize = outputImage->ySize * sizeof(float *);
	fprintf(stderr, "%i\n", bufSize);
	outputImage->image = (void **)malloc(bufSize);
	outputImage->image2 = (void **)malloc(bufSize);
	outputImage->image3 = (void **)malloc(bufSize);
	outputImage->scale = (float **)malloc(bufSize);
	outputImage->scale2 = (float **)malloc(bufSize);
	outputImage->scale3 = (float **)malloc(bufSize);
	outputImage->vxTmp = (float **)malloc(bufSize);
	outputImage->sxTmp = (float **)malloc(bufSize);
	outputImage->vyTmp = (float **)malloc(bufSize);
	outputImage->syTmp = (float **)malloc(bufSize);
	outputImage->vzTmp = (float **)malloc(bufSize);
	outputImage->fScale = (float **)malloc(bufSize);
	outputImage->errorX = (float **)malloc(bufSize);
	outputImage->errorY = (float **)malloc(bufSize);

	bufSize = outputImage->xSize * outputImage->ySize * sizeof(float);
	buf1 = (float *)malloc(bufSize);
	buf2 = (float *)malloc(bufSize);
	buf3 = (float *)malloc(bufSize);
	buf1s = (float *)malloc(bufSize);
	buf2s = (float *)malloc(bufSize);
	buf3s = (float *)malloc(bufSize);
	bufx = (float *)malloc(bufSize);
	bufy = (float *)malloc(bufSize);
	bufz = (float *)malloc(bufSize);
	bufs = (float *)malloc(bufSize);
	bufsx = (float *)malloc(bufSize);
	bufsy = (float *)malloc(bufSize);
	bufex = (float *)malloc(bufSize);
	bufey = (float *)malloc(bufSize);

	for (i = 0; i < outputImage->ySize; i++)
	{
		outputImage->image[i] = (void *)&(buf1[i * outputImage->xSize]);
		outputImage->image2[i] = (void *)&(buf2[i * outputImage->xSize]);
		outputImage->image3[i] = (void *)&(buf3[i * outputImage->xSize]);
		outputImage->scale[i] = (float *)&(buf1s[i * outputImage->xSize]);
		outputImage->scale2[i] = (float *)&(buf2s[i * outputImage->xSize]);
		outputImage->scale3[i] = (float *)&(buf3s[i * outputImage->xSize]);
		outputImage->vxTmp[i] = (float *)&(bufx[i * outputImage->xSize]);
		outputImage->vyTmp[i] = (float *)&(bufy[i * outputImage->xSize]);
		outputImage->vzTmp[i] = (float *)&(bufz[i * outputImage->xSize]);
		outputImage->sxTmp[i] = (float *)&(bufsx[i * outputImage->xSize]);
		outputImage->syTmp[i] = (float *)&(bufsy[i * outputImage->xSize]);
		outputImage->fScale[i] = (float *)&(bufs[i * outputImage->xSize]);
		outputImage->errorX[i] = (float *)&(bufex[i * outputImage->xSize]);
		outputImage->errorY[i] = (float *)&(bufey[i * outputImage->xSize]);
	}

	vXimage = (float **)outputImage->image;
	vYimage = (float **)outputImage->image2;
	vZimage = (float **)outputImage->image3;
	scaleX = (float **)outputImage->scale;
	scaleY = (float **)outputImage->scale2;
	scaleZ = (float **)outputImage->scale3;

	for (j = 0; j < outputImage->ySize; j++)
	{
		for (k = 0; k < outputImage->xSize; k++)
		{
			scaleX[j][k] = 0.0;
			vXimage[j][k] = -LARGEINT;
			scaleY[j][k] = 0.0;
			vYimage[j][k] = -LARGEINT;
			scaleZ[j][k] = 0.0;
			vZimage[j][k] = -LARGEINT;
		}
	}
	fprintf(outputImage->fpLog, "; Returned from  mallocOutputImage\n");
	fflush(outputImage->fpLog);
}

static void readReferenceVelMosaic(referenceVelocity *refVel, outputImageStructure *outputImage)
{
	FILE *fp;
	char *geodatFile, *vxFile, *vyFile, *exFile, *eyFile;
	double maxX, maxY;
	char line[1500];
	uint32_t lineCount = 0;
	int32_t eod;
	float *tmp, *tmp1;
	float dum1, dum2;
	uint32_t nx, ny;
	double x0, y0;
	uint32_t xoff, yoff, tail; /* offset into velocity file in samples */
	uint32_t i;
	/*
	  geodat file name
	*/
	fprintf(stderr, "**** Start reading reference velocity file ****\n");
	geodatFile = (char *)malloc(strlen(refVel->velFile) + 11);
	geodatFile[0] = '\0';
	geodatFile = strcpy(geodatFile, refVel->velFile);
	geodatFile = strcat(geodatFile, ".vx.geodat");
	fprintf(stderr, "vel geodat file %s\n", geodatFile);
	/*
	   velocity file names
	*/
	vxFile = (char *)malloc(strlen(refVel->velFile) + 4);
	vxFile[0] = '\0';
	vxFile = strcpy(vxFile, refVel->velFile);
	vxFile = strcat(vxFile, ".vx");
	vyFile = (char *)malloc(strlen(refVel->velFile) + 4);
	vyFile[0] = '\0';
	vyFile = strcpy(vyFile, refVel->velFile);
	vyFile = strcat(vyFile, ".vy");
	fprintf(stderr, "vx,vy file %s %s\n", vxFile, vyFile);
	/*
	  Open geodat file
	*/
	fp = openInputFile(geodatFile);
	if (fp == NULL)
		error("*** readREfVelFile: Error opening %s ***\n", geodatFile);
	/*
	  Read parameters
	*/
	lineCount = getDataString(fp, lineCount, line, &eod); /* Skip # 2 line */
	lineCount = getDataString(fp, lineCount, line, &eod);
	sscanf(line, "%f %f\n", &dum1, &dum2); /* read as float in case fp value */
	nx = (int)dum1;
	ny = (int)dum2;

	lineCount = getDataString(fp, lineCount, line, &eod);
	sscanf(line, "%lf %lf\n", &(refVel->dx), &(refVel->dy));

	lineCount = getDataString(fp, lineCount, line, &eod);
	sscanf(line, "%lf %lf\n", &(x0), &(y0));
	x0 *= KMTOM;
	y0 *= KMTOM;
	fclose(fp);
	fprintf(stderr, "%i %i \n %f %f \n %f %f \n", nx, ny, refVel->dx, refVel->dy, x0, y0);
	refVel->x0 = max(x0, outputImage->originX);
	refVel->y0 = max(y0, outputImage->originY);
	fprintf(stderr, "output grid %lf %lf \n", outputImage->originX, outputImage->originY);
	maxX = min(x0 + refVel->dx * (nx - 1), outputImage->originX + outputImage->deltaX * (outputImage->xSize - 1));
	maxY = min(y0 + refVel->dy * (ny - 1), outputImage->originY + outputImage->deltaY * (outputImage->ySize - 1));
	refVel->nx = (maxX - refVel->x0) / (refVel->dx) + 1;
	refVel->ny = (maxY - refVel->y0) / (refVel->dy) + 1;
	fprintf(stderr, "%f %f  %i %i \n", refVel->x0, refVel->y0, refVel->nx, refVel->ny);
	/* compute file offsets */
	xoff = (uint32_t)((refVel->x0 - x0) / refVel->dx + 0.5);
	yoff = (uint32_t)((refVel->y0 - y0) / refVel->dy + 0.5);
	fprintf(stderr, "xoff,yoff %i %i\n", xoff, yoff);
	/*
	  Malloc array
	*/
	refVel->vx = (float **)malloc(refVel->ny * sizeof(float *));
	tmp = (float *)malloc(refVel->nx * refVel->ny * sizeof(float));
	/*
	  Open vx file
	*/
	if (vxFile == NULL)
		error("*** readVelMosaic: Error opening %s ***\n", vxFile);
	fp = openInputFile(vxFile);
	fseek(fp, (nx * yoff * sizeof(float)), SEEK_SET);
	tail = max(nx - (xoff + refVel->nx), 0);

	if (fp == NULL)
		error("*** readVelMosaic: Error opening %s ***\n", vxFile);
	for (i = 0; i < refVel->ny; i++)
	{
		tmp1 = &(tmp[i * refVel->nx]);
		fseek(fp, (xoff * sizeof(float)), SEEK_CUR);
		freadBS(tmp1, sizeof(float), refVel->nx, fp, FLOAT32FLAG);
		fseek(fp, (tail * sizeof(float)), SEEK_CUR);
		refVel->vx[i] = tmp1;
	}
	fclose(fp);
	/*
	  Malloc array
	*/
	refVel->vy = (float **)malloc(refVel->ny * sizeof(float *));
	tmp = (float *)malloc(refVel->nx * refVel->ny * sizeof(float));
	/*
	  Open vy file
	*/
	if (vyFile == NULL)
		error("*** readVel: Error opening %s ***\n", vyFile);
	fp = openInputFile(vyFile);
	fseek(fp, (nx * yoff * sizeof(float)), SEEK_SET);
	if (fp == NULL)
		error("*** readVel: Error opening %s ***\n", vyFile);
	for (i = 0; i < refVel->ny; i++)
	{
		tmp1 = &(tmp[i * refVel->nx]);
		fseek(fp, (xoff * sizeof(float)), SEEK_CUR);
		freadBS(tmp1, sizeof(float), refVel->nx, fp, FLOAT32FLAG);
		fseek(fp, (tail * sizeof(float)), SEEK_CUR);
		refVel->vy[i] = tmp1;
	}
	fclose(fp);

	/* Error files  if needed */
	if (refVel->initMapFlag == TRUE)
	{
		/*
		   velocity file names
		*/
		exFile = (char *)malloc(strlen(refVel->velFile) + 4);
		exFile[0] = '\0';
		exFile = strcpy(exFile, refVel->velFile);
		exFile = strcat(exFile, ".ex");
		eyFile = (char *)malloc(strlen(refVel->velFile) + 4);
		eyFile[0] = '\0';
		eyFile = strcpy(eyFile, refVel->velFile);
		eyFile = strcat(eyFile, ".ey");
		/*
		  Malloc array
		*/
		refVel->ex = (float **)malloc(refVel->ny * sizeof(float *));
		tmp = (float *)malloc(refVel->nx * refVel->ny * sizeof(float));
		/*
		  Open ex file
		*/
		if (exFile == NULL)
			error("*** readVelMosaic: Error opening %s ***\n", exFile);
		fp = openInputFile(exFile);
		fseek(fp, (nx * yoff * sizeof(float)), SEEK_SET);
		tail = max(nx - (xoff + refVel->nx), 0);

		if (fp == NULL)
			error("*** readVelMosaic: Error opening %s ***\n", exFile);
		for (i = 0; i < refVel->ny; i++)
		{
			tmp1 = &(tmp[i * refVel->nx]);
			fseek(fp, (xoff * sizeof(float)), SEEK_CUR);
			freadBS(tmp1, sizeof(float), refVel->nx, fp, FLOAT32FLAG);
			fseek(fp, (tail * sizeof(float)), SEEK_CUR);
			refVel->ex[i] = tmp1;
		}
		fclose(fp);
		/*
		  Malloc array
		*/
		refVel->ey = (float **)malloc(refVel->ny * sizeof(float *));
		tmp = (float *)malloc(refVel->nx * refVel->ny * sizeof(float));
		/*
		  Open ey file
		*/
		if (eyFile == NULL)
			error("*** readVelMosaic: Error opening %s ***\n", eyFile);
		fp = openInputFile(eyFile);
		fseek(fp, (nx * yoff * sizeof(float)), SEEK_SET);
		tail = max(nx - (xoff + refVel->nx), 0);
		if (fp == NULL)
			error("*** readVelMosaic: Error opening %s ***\n", eyFile);
		for (i = 0; i < refVel->ny; i++)
		{
			tmp1 = &(tmp[i * refVel->nx]);
			fseek(fp, (xoff * sizeof(float)), SEEK_CUR);
			freadBS(tmp1, sizeof(float), refVel->nx, fp, FLOAT32FLAG);
			fseek(fp, (tail * sizeof(float)), SEEK_CUR);
			refVel->ey[i] = tmp1;
		}
		fclose(fp);
		fprintf(stderr, "**** reading error map ****\n");
	}
	fprintf(stderr, "**** End reading reference velocity file ****\n");
}

#define IGREG 2299161

void caldat(int32_t julian, int32_t *mm, int32_t *id, int32_t *iyyy)
{
	int32_t ja, jalpha, jb, jc, jd, je;

	if (julian >= IGREG)
	{
		jalpha = (int32_t)(((float)(julian - 1867216) - 0.25) / 36524.25);
		ja = julian + 1 + jalpha - (int32_t)(0.25 * jalpha);
	}
	else
		ja = julian;
	jb = ja + 1524;
	jc = (int32_t)(6680.0 + ((float)(jb - 2439870) - 122.1) / 365.25);
	jd = (int32_t)(365 * jc + (0.25 * jc));
	je = (int32_t)((jb - jd) / 30.6001);
	*id = jb - jd - (int32_t)(30.6001 * je);
	*mm = je - 1;
	if (*mm > 12)
		*mm -= 12;
	*iyyy = jc - 4715;
	if (*mm > 2)
		--(*iyyy);
	if (*iyyy <= 0)
		--(*iyyy);
}
#undef IGREG
