#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"
/*
#include "mosaicSource/GeoCodeDEM_p/geocodedem.h"
#include "mosaicSource/computeVH_p/computevh.h"*/
#include "common.h"

static char *dupName(char *phase)
{
	int32_t j;
	char *tmp;
	if (strstr(phase, "none") || strstr(phase, "None") || phase[0] == '\0')
	{
		return (NULL);
	}
	tmp = (char *)malloc(strlen(phase) + 1);
	for (j = 0; j < strlen(phase); j++)
		tmp[j] = phase[j];
	tmp[j] = '\0';
	return tmp;
}
/*
  Process input file for velocity mosaicking
*/
void getMVhInputFile(char *inputFile, char ***phaseFiles, char ***geodatFiles, char ***baselineFiles, char ***offsetFiles,
					 char ***azParamsFiles, char ***rOffsetFiles, char ***rParamsFiles, outputImageStructure *outputImage, float **nDays,
					 float **weights, int32_t **crossFlags, int32_t *nFiles, int32_t offsetFlag, int32_t rOffsetFlag, int32_t threeDOffFlag)
{
	FILE *fp;
	double xo, yo, xs, ys, deltaX, deltaY;
	float weight, nDay;
	int32_t lineCount, eod;
	int32_t i, j;
	char phase[1024], geodat[1024], baseline[1024];
	char offsets[1024], rOffsets[1024], azParams[1024], rParams[1024];
	char line[1024];
	int32_t crossFlag;
	/*
	  Open file for input 
	*/
	fp = openInputFile(inputFile);
	/*
	  Input nr,na,nlr,nla
	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%lf%lf%lf%lf%lf%lf", &xo, &yo, &xs, &ys, &deltaX, &deltaY) != 6)
		error("%s  %i of %s", "getMVhInputFile -- Missing image parameters at line:", lineCount, inputFile);

	outputImage->xSize = (int)(xs / deltaX + 0.5);
	outputImage->ySize = (int)(ys / deltaY + 0.5);
	outputImage->deltaX = deltaX * KMTOM;
	outputImage->deltaY = deltaY * KMTOM;
	outputImage->originX = xo * KMTOM;
	outputImage->originY = yo * KMTOM;
	/*
	  Input ReMajor,ReMinor, Rc, phic, H
	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%i", nFiles) != 1)
		error("%s  %i  of %s", "getMVhInputFile -- Missing geometric parameters at line:", lineCount, inputFile);
	/*
	  Malloc space for arrays of filenames
	*/
	*phaseFiles = (char **)malloc(sizeof(char *) * (*nFiles));
	*geodatFiles = (char **)malloc(sizeof(char *) * (*nFiles));
	*baselineFiles = (char **)malloc(sizeof(char *) * (*nFiles));
	if (offsetFlag == TRUE || threeDOffFlag == TRUE)
	{
		*offsetFiles = (char **)malloc(sizeof(char *) * (*nFiles));
		*azParamsFiles = (char **)malloc(sizeof(char *) * (*nFiles));
	}
	if (rOffsetFlag == TRUE || threeDOffFlag == TRUE)
	{
		*rOffsetFiles = (char **)malloc(sizeof(char *) * (*nFiles));
		*rParamsFiles = (char **)malloc(sizeof(char *) * (*nFiles));
		if (offsetFlag == FALSE)
		{ /* Case where offsets only for 2d */
			*azParamsFiles = (char **)malloc(sizeof(char *) * (*nFiles));
			*offsetFiles = (char **)malloc(sizeof(char *) * (*nFiles));
		}
	}
	*nDays = (float *)malloc(sizeof(float) * (*nFiles));
	*weights = (float *)malloc(sizeof(float) * (*nFiles));
	*crossFlags = (int32_t *)malloc(sizeof(int) * (*nFiles));
	/*
	  Input files
	*/
	for (i = 0; i < *nFiles; i++)
	{
		lineCount = getDataString(fp, lineCount, line, &eod);
		crossFlag = TRUE;
		if ((offsetFlag == FALSE && rOffsetFlag == FALSE) && threeDOffFlag == FALSE)
		{
			if (sscanf(line, "%s %s %s %f %f", phase, geodat, baseline, &nDay, &weight) != 5)
			{
				weight = 1.0;
				if (sscanf(line, "%s %s %s %f", phase, geodat, baseline, &nDay) != 4)
					error("%s  %i  of %s", "getMVhInputFile -- Missing filename at line:", lineCount, inputFile);
			}
		}
		else if ((offsetFlag == TRUE && rOffsetFlag == FALSE) && threeDOffFlag == FALSE)
		{
			/* az but no range offsets */
			if (sscanf(line, "%s %s %s %f %f %s %s", phase, geodat, baseline, &nDay, &weight, offsets, azParams) != 7)
			{
				weight = 1.0;
				if (sscanf(line, "%s %s %s %f %s %s", phase, geodat, baseline, &nDay, offsets, azParams) != 6)
				{
					rOffsets[0] = '\0';
					rParams[0] = '\0';
					offsets[0] = '\0';
					azParams[0] = '\0';
					if (sscanf(line, "%s %s %s %f %f", phase, geodat, baseline, &nDay, &weight) != 5)
						error("%s  %i  of %s", "getMVhInputFile -- Missing filename at line:", lineCount, inputFile);
				}
			}
		}
		else if (rOffsetFlag == TRUE || threeDOffFlag == TRUE)
		{
			if (sscanf(line, "%s %s %s %f %f %s %s %s %s %i", phase, geodat, baseline, &nDay, &weight, offsets, azParams, rOffsets, rParams, &crossFlag) != 10)
			{
				crossFlag = TRUE;
				if (sscanf(line, "%s %s %s %f %f %s %s %s %s", phase, geodat, baseline, &nDay, &weight, offsets, azParams, rOffsets, rParams) != 9)
				{
					weight = 1.0;
					if (sscanf(line, "%s %s %s %f %s %s %s %s", phase, geodat, baseline, &nDay, offsets, azParams,
							   rOffsets, rParams) != 8)
					{
						if (sscanf(line, "%s %s %s %f %f %s %s", phase, geodat, baseline, &nDay, &weight, offsets, azParams) != 7)
						{
							if (sscanf(line, "%s %s %s %f %s %s", phase, geodat, baseline, &nDay, offsets, azParams) != 6)
							{
								offsets[0] = '\0';
								azParams[0] = '\0';
								if (sscanf(line, "%s %s %s %f %f", phase, geodat, baseline, &nDay, &weight) != 5)
									error("%s  %i  of %s", "getMVhInputFile -- Missing filename at line:", lineCount, inputFile);
							}
						}
						rOffsets[0] = '\0';
						rParams[0] = '\0';
					}
				}
			}
		}
		(*phaseFiles)[i] = dupName(phase);
		(*geodatFiles)[i] = dupName(geodat);
		(*baselineFiles)[i] = dupName(baseline);

		if (offsetFlag == TRUE || rOffsetFlag == TRUE || threeDOffFlag == TRUE)
		{
			(*offsetFiles)[i] = dupName(offsets);
			if ((*offsetFiles)[i] != NULL)
				(*azParamsFiles)[i] = dupName(azParams);
		}
		if (rOffsetFlag == TRUE || threeDOffFlag == TRUE)
		{
			(*rOffsetFiles)[i] = dupName(rOffsets);
			if ((*rOffsetFiles)[i] != NULL)
				(*rParamsFiles)[i] = dupName(rParams);
		}
		(*crossFlags)[i] = crossFlag;
		(*weights)[i] = weight;
		(*nDays)[i] = nDay;
	}
	fclose(fp);
}
