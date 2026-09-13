#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"
/*
#include "mosaicSource/GeoCodeDEM_p/geocodedem.h"
#include "mosaicSource/computeVH_p/computevh.h"*/
#include "common.h"

/*
  Append the vertical-correction suffix to a baseline filename. YAML baseline
  files (extension .yaml, see getBaseline.c) must keep ".yaml" as the final
  extension for getBaseline's format-detection check, so the suffix is
  inserted before it (e.g. baseline.26x16.480.yaml -> baseline.26x16.480.surf.yaml).
  Non-yaml baseline files get the suffix appended as before.
*/
extern int32_t useAzIonosphere; /* defined in common/getRegion.c */

static char *appendBaselineSuffix(char *baseline, char *verticalCorrectionSuffix, char *buf)
{
	size_t blen = strlen(baseline);
	if (blen > 5 && strcmp(baseline + blen - 5, ".yaml") == 0)
	{
		strncpy(buf, baseline, blen - 5);
		buf[blen - 5] = '\0';
		strcat(buf, verticalCorrectionSuffix);
		strcat(buf, ".yaml");
		return buf;
	}
	return appendSuffix(baseline, verticalCorrectionSuffix, buf);
}

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
	int32_t i, j, len, ii;
	char phase[2048], phaseTmp[2048], geodat[1024], baseline[1024];
	char offsets[1024], rOffsets[1024], azParams[1024], rParams[1024];
	char line[1024];
	char *verticalCorrectionSuffix=NULL;
	int32_t crossFlag;
	if(outputImage->verticalCorrectionSuffix != NULL)
	{
		len = strlen(outputImage->verticalCorrectionSuffix);
		verticalCorrectionSuffix = (char *)malloc(sizeof(char) * (len+2));
		verticalCorrectionSuffix[0] = '.';
		for(ii=1;  ii <=min(len,2048); ii++) verticalCorrectionSuffix[ii] = outputImage->verticalCorrectionSuffix[ii-1];
		verticalCorrectionSuffix[len+1] = '\0';
	}
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
		// Suffix for vertical correction if its different than default (none)
		// Only the baseline filename gets the suffix -- the phase file is shared
		// across vertical-correction scenarios; only its baseline fit differs.
		(*phaseFiles)[i] = dupName(phase);
		if(verticalCorrectionSuffix != NULL)
		{
			(*baselineFiles)[i] = appendBaselineSuffix(baseline, verticalCorrectionSuffix,
				(char *)malloc(strlen(baseline)+ strlen(verticalCorrectionSuffix) + 1));
		}
		else
		{
			(*baselineFiles)[i] = dupName(baseline);
		}
		/* -iceOnly: select the ice-only PHASE baseline written by tiepoints -iceOnly.
		   Deliberately applied only here, not to rParams (readOffsets.c), so an
		   offsets run is unaffected and existing baselines are never overwritten. */
		if (outputImage->iceOnly == TRUE && (*baselineFiles)[i] != NULL)
		{
			char *bTmp = (*baselineFiles)[i];
			(*baselineFiles)[i] = appendBaselineSuffix(bTmp, ".iceOnly",
				(char *)malloc(strlen(bTmp) + strlen(".iceOnly") + 1));
			free(bTmp);
		}
		/* -flipSquint: same pattern, so tiepoints -flipSquint output is picked up */
		if (outputImage->flipSquint == TRUE && (*baselineFiles)[i] != NULL)
		{
			char *bTmp = (*baselineFiles)[i];
			(*baselineFiles)[i] = appendBaselineSuffix(bTmp, ".flipSquint",
				(char *)malloc(strlen(bTmp) + strlen(".flipSquint") + 1));
			free(bTmp);
		}
		(*geodatFiles)[i] = dupName(geodat);
		if (offsetFlag == TRUE || rOffsetFlag == TRUE || threeDOffFlag == TRUE)
		{
			(*offsetFiles)[i] = dupName(offsets);
			if ((*offsetFiles)[i] != NULL)
			{
				(*azParamsFiles)[i] = dupName(azParams);
				/* -useAzIonosphere: select the azimuth fit that azparams wrote
				   while weighing the ionosphere correction (tieScript writes it
				   to az.est.azIon*.yaml).  Same pattern as -iceOnly above: the
				   production az.est*.yaml is never touched, so turning the flag
				   off reverts instantly and a reference run is always available. */
				if (useAzIonosphere == TRUE && (*azParamsFiles)[i] != NULL)
				{
					char *aTmp = (*azParamsFiles)[i];
					(*azParamsFiles)[i] = appendBaselineSuffix(aTmp, ".azIon",
						(char *)malloc(strlen(aTmp) + strlen(".azIon") + 1));
					free(aTmp);
				}
			}
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
