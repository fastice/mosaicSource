#include "stdio.h"
#include "string.h"
#include "ctype.h"
#include "math.h"
#include <stdlib.h>
#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"
#include "gdalIO/gdalIO/grimpgdal.h"

static GDALRasterBandH getVRTOffsetMeta(char *datFile, int32_t *rO, int32_t *aO, int32_t *nr, int32_t *na, float *deltaR, float *deltaA)
{
	dictNode *metaData = NULL;
	// GDALRasterBandH hBand;
	GDALDatasetH hDS;
	float tmp;
	// Get meta data
	//hBand = GDALGetRasterBand(hDS, 0);
	hDS = GDALOpen(datFile, GDAL_OF_READONLY);
	*nr = GDALGetRasterXSize(hDS);
	*na = GDALGetRasterYSize(hDS);
	readDataSetMetaData(hDS, &metaData);
	*rO = atoi(get_value(metaData, "r0"));
	*aO = atoi(get_value(metaData, "a0"));
	*deltaR = atof(get_value(metaData, "deltaR"));
	*deltaA = atof(get_value(metaData, "deltaA"));
	GDALClose(hDS);
}

/*
 Determine if offset .dat or .vrt given, and read the results
*/
static void parseOffsetParamFile(char *datFile, int32_t *rO, int32_t *aO, int32_t *nr, int32_t *na, float *deltaR, float *deltaA)
{
	int32_t nChar, isVRT = TRUE, isDAT = TRUE;
	char *datTemplate = "dat", *vrtTemplate = "vrt";
	int32_t lineCount, eod, nRead;
	char line[1024];
	GDALDatasetH hDS;
	FILE *fp;
	/*
	Check for dat or vrt suffix.
	*/
	nChar = strlen(datFile) - 3;
	for (int i = 0; i < 3; i++)
	{
		if (datFile[nChar + i] != vrtTemplate[i])
			isVRT = FALSE;
		if (datFile[nChar + i] != datTemplate[i])
			isDAT = FALSE;
	}
	if (isVRT == TRUE)
	{
		fprintf(stderr, "Found VRT\n");
		getVRTOffsetMeta(datFile, rO, aO, nr, na, deltaR, deltaA);
	}
	else if (isDAT == TRUE)
	{
		fprintf(stderr, "Found DAT\n");
		fp = openInputFile(datFile);
		lineCount = getDataString(fp, lineCount, line, &eod);
		nRead = sscanf(line, "%i%i%i%i%f%f", rO, aO, nr, na, deltaR, deltaA);
		fclose(fp);
	}
	else
		error("%s has neither a .dat or .vrt suffix", datFile);
}
/*
  parse scene input file for siminsar.
*/
void parseSceneFile(char *sceneFile, sceneStructure *scene)
{
	FILE *fp;
	int32_t i;
	char *datFile, buf[1024], buf1[1024], *errorFile;
	int32_t rO, aO, nr, na;
	float deltaR, deltaA;
	int32_t lineCount, eod, nRead;
	char line[1024];
	parseInputFile(sceneFile, &(scene->I));
	/*
	  Compute baseline increment size
	*/
	scene->bnStep = (scene->bnEnd - scene->bnStart) / (double)(scene->I.azimuthSize - 1.0);
	scene->bpStep = (scene->bpEnd - scene->bpStart) / (double)(scene->I.azimuthSize - 1.0);

	if (scene->toLLFlag == TRUE)
	{
		/*
		  Read inputfile with description of offset geometry
		*/
		datFile = scene->llInput;
		parseOffsetParamFile(datFile, &rO, &aO, &nr, &na, &deltaR, &deltaA);
		fprintf(stderr, "offset params %i %i %i %i %f %f\n", rO, aO, nr, na, deltaR, deltaA);
		
		scene->aSize = na;
		scene->rSize = nr;
		scene->aO = aO / scene->I.nAzimuthLooks;
		scene->rO = rO / scene->I.nRangeLooks;
		scene->dR = deltaR / scene->I.nRangeLooks;
		scene->dA = deltaA / scene->I.nAzimuthLooks;
		scene->latImage = (double **)malloc(scene->aSize * sizeof(double *));
		scene->lonImage = (double **)malloc(scene->aSize * sizeof(double *));

		for (i = 0; i < scene->aSize; i++)
		{
			scene->lonImage[i] = (double *)malloc(scene->rSize * sizeof(double));
			scene->latImage[i] = (double *)malloc(scene->rSize * sizeof(double));
		}
		if (scene->maskFlag == TRUE || scene->heightFlag == TRUE)
		{
			scene->image = (float **)malloc(scene->aSize * sizeof(float *));
			for (i = 0; i < scene->aSize; i++)
			{
				scene->image[i] = (float *)malloc(scene->rSize * sizeof(float));
			}
		}
		fprintf(stderr, "reading %s\n", scene->llInput);
	}
	else
	{
		scene->aSize = scene->I.azimuthSize;
		scene->rSize = scene->I.rangeSize;
		scene->aO = 0;
		scene->rO = 0;
		scene->dR = 1;
		scene->dA = 1;
		/* Save lat/lon as */
		if (scene->saveLLFlag == TRUE)
		{
			scene->latImage = (double **)malloc(scene->aSize * sizeof(double *));
			scene->lonImage = (double **)malloc(scene->aSize * sizeof(double *));
			for (i = 0; i < scene->aSize; i++)
			{
				scene->lonImage[i] = (double *)malloc(scene->rSize * sizeof(double));
				scene->latImage[i] = (double *)malloc(scene->rSize * sizeof(double));
			}
		}
		/*
		  Init space for image
		*/
		scene->image = (float **)malloc(scene->aSize * sizeof(float *));
		for (i = 0; i < scene->aSize; i++)
		{
			scene->image[i] = (float *)malloc(scene->rSize * sizeof(float));
		}
	}
	return;
}
