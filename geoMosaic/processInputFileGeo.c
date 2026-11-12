#include "stdio.h"
#include "string.h"
#include "math.h"
#include "mosaicSource/common/common.h"
/*#include "ers1/GeoCodeDEM_p/geocodedem.h"*/
#include "geomosaic.h"
#include <stdlib.h>

/*
  Process input file for mosaicDEMs
*/
void processInputFileGeo(char *inputFile, char ***insarDEMFiles, char ***demInputFiles,
						 outputImageStructure *outputImage, int *nDEMs, float **weights, char ***antPatFiles,
						 int32_t **tiers)
{
	FILE *fp;
	double xo, yo, xs, ys, deltaX, deltaY;
	int lineCount = 0, eod;
	int i, j;
	float w;
	int tier;
	char *tierTok;
	char insarDEM[1024], demInput[1024], antPat[1024];
	char line[512];
	/*
	  Open file for input
	*/
	fp = openInputFile(inputFile);
	/*
	  Input nr,na,nlr,nla
	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%lf%lf%lf%lf%lf%lf", &xo, &yo, &xs, &ys, &deltaX, &deltaY) != 6)
		error("%s  %i of %s", "processInputFile -- Missing image parameters at line:", lineCount, inputFile);

	outputImage->xSize = (int)(xs / deltaX + 0.5);
	outputImage->ySize = (int)(ys / deltaY + 0.5);
	/* The header is in the projection's OWN working units -- kilometres for a projected grid,
	   DEGREES for a geographic one -- and is stored here in the units the pixel loops use:
	   metres, or degrees unchanged. gridScale is the one place that distinction lives
	   (MTOKM for projected, 1.0 for GP_LATLON), so dividing by it is the inverse of the
	   MTOKM the conversion sites apply. Hard-coding KMTOM here turned a -122.3 degree origin
	   into -122300. */
	{
		double toStored = 1.0 / outputImage->proj.gridScale;
		outputImage->deltaX = deltaX * toStored;
		outputImage->deltaY = deltaY * toStored;
		outputImage->originX = xo * toStored;
		outputImage->originY = yo * toStored;
	}
	/*
	  Input ReMajor,ReMinor, Rc, phic, H
	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%i", nDEMs) != 1)
		error("%s  %i  of %s",
			  "process InputFile -- Missing geometric parameters at line:",
			  lineCount, inputFile);
	/*
	  Malloc space for arrays of filenames
	*/
	*insarDEMFiles = (char **)malloc((size_t)(sizeof(char *) * (*nDEMs)));
	*demInputFiles = (char **)malloc((size_t)(sizeof(char *) * (*nDEMs)));
	*antPatFiles = (char **)malloc((size_t)(sizeof(char *) * (*nDEMs)));
	*weights = (float *)malloc((size_t)(sizeof(float) * (*nDEMs)));
	*tiers = (int32_t *)malloc((size_t)(sizeof(int32_t) * (*nDEMs)));
	/*
	  Input files
	*/
	for (i = 0; i < *nDEMs; i++)
	{
		lineCount = getDataString(fp, lineCount, line, &eod);

		if (sscanf(line, "%s %s %f %s", insarDEM, demInput, &w, antPat) != 4)
		{
			antPat[0] = '\0';
			if (sscanf(line, "%s %s %f", insarDEM, demInput, &w) != 3)
			{
				w = 1.0;
				if (sscanf(line, "%s %s", insarDEM, demInput) != 2)
					error("%s  %i  of %s", "process InputFile -- Missing filename at line:", lineCount, inputFile);
			}
		}
		/*
		  Optional priority tier, a "tier=N" token after the filenames (0, the default, is the
		  highest priority; see makeGeoMosaic). It can land in the antenna-pattern slot, so it
		  is cleared from there. Tiers must not decrease down the list.
		*/
		(*tiers)[i] = 0;
		if (strncmp(antPat, "tier=", 5) == 0)
		{
			antPat[0] = '\0';
		}
		tierTok = strstr(line, "tier=");
		if (tierTok != NULL && tierTok > line && (tierTok[-1] == ' ' || tierTok[-1] == '\t'))
		{
			if (sscanf(tierTok + 5, "%d", &tier) != 1 || tier < 0 || tier > MAXTIER)
			{
				error("process InputFile -- bad tier (0-%i) at line %i of %s", MAXTIER, lineCount, inputFile);
			}
			(*tiers)[i] = tier;
		}
		if (i > 0 && (*tiers)[i] < (*tiers)[i - 1])
		{
			error("process InputFile -- tiers must not decrease (tier %i after %i) at line %i of %s",
				  (*tiers)[i], (*tiers)[i - 1], lineCount, inputFile);
		}
		(*insarDEMFiles)[i] = (char *)malloc((size_t)(strlen(insarDEM) + 1));
		for (j = 0; j < strlen(insarDEM); j++)
			(*insarDEMFiles)[i][j] = insarDEM[j];
		(*insarDEMFiles)[i][j] = '\0';

		(*demInputFiles)[i] = (char *)malloc((size_t)(strlen(demInput) + 1));
		for (j = 0; j < strlen(demInput); j++)
			(*demInputFiles)[i][j] = demInput[j];
		(*demInputFiles)[i][j] = '\0';
		if (strlen(antPat) > 0)
		{
			(*antPatFiles)[i] = (char *)malloc((size_t)(strlen(antPat) + 1));
			for (j = 0; j < strlen(antPat); j++)
				(*antPatFiles)[i][j] = antPat[j];
			(*antPatFiles)[i][j] = '\0';
		}
		else
			(*antPatFiles)[i] = NULL;

		(*weights)[i] = w;
	}
	return;
}
