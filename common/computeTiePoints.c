#include <math.h>
#include <stdlib.h>
#include "string.h"
#include <omp.h>
#include "mosaicSource/common/common.h"
/*#include "tiePoints.h"*/
#include "time.h"
extern double Rotation;
/* Points per warm-start chunk; also the unit of OpenMP work (see the loop below) */
#define TIEPOINTCHUNK 100

/*
  Compute location of  tie points for baselines estimation.
*/

void computeTiePoints(inputImageStructure *inputImage, tiePointsStructure *tiePoints,
					  demStructure dem, int32_t noDEM, char *geodatFile, ShelfMask *shelfMask, int32_t supressOutput)
{
	FILE *fp;
	int32_t i, j;
	conversionDataStructure *cP;
	double x, y, zSp, zWGS;
	double xx, yy;
	double lat, lon;
	double azOrigin;
	double num, den, rho, r;
	double RElip;
	double Re;
	double range, azimuth;
	double deltaZ;
	double vz1;
	int32_t sMask;
	int32_t shelfMaskFlag;
	char tideFile[2048], *tmpS;
	float **fimage;
	clock_t startTime, lastTime, initTime;
	int32_t k = 0;
	startTime = clock();
	initllToImageNew(inputImage);
	if (!supressOutput) fprintf(stderr, "Using new \n");
	initTime = clock();
	cP = &(inputImage->cpAll);
	Re = cP->Re;
	/*shelfMaskFlag = (int32_t)shelfMask;*/

	if (supressOutput == FALSE)
	{
		fprintf(stdout, ";;\n;; H, InSAR parameters  Re, RNear, Rc, nlr,nla\n;;\n");
		fprintf(stdout, " %10.3f %10.3f  %8.3f %8.3f  %i %i %f %f\n",
				inputImage->par.H * MTOKM, cP->Re * MTOKM, cP->RNear * MTOKM, inputImage->par.rc * MTOKM,
				inputImage->nRangeLooks, inputImage->nAzimuthLooks, inputImage->rangePixelSize / inputImage->nRangeLooks, inputImage->par.lambda);
	}
	/*
	  Read spatially varying tide correction
	*/
	tideFile[0] = '\0';
	strcpy(tideFile, geodatFile);
	tmpS = strrchr(tideFile, '/');
	if (tmpS != NULL)
		tmpS[0] = '\0';
	strcat(tideFile, "/tide.difference");
	fp = fopen(tideFile, "r");
	if (fp != NULL)
	{
		inputImage->tideDiffFlag = TRUE;
		readXYDEM(tideFile, &(inputImage->tideDiff));
		fprintf(stderr, ";; Tide difference : %s\n", tideFile);
		if (supressOutput == FALSE)
			fprintf(stdout, ";; Tide difference : %s\n", tideFile);
	}
	else
	{
		inputImage->tideDiffFlag = FALSE;
		if (supressOutput == FALSE)
			fprintf(stdout, ";; Tide difference : %s\n", "none");
	}
	/*
	  Loop over output points — parallelized with OpenMP over fixed chunks of
	  TIEPOINTCHUNK points. llToImageNew warm-starts its Newton iteration from
	  lastTime, and Newton stops anywhere inside its tolerance, so the result
	  depends on the start time. A shared lastTime made that whichever point some
	  thread had just finished, and the output changed from run to run (azimuth by
	  up to ~4e-4 pixels). Instead each thread works on its own copy of inputImage,
	  and lastTime is reset at the start of every chunk (out of range, so
	  llToImageNew starts from centre time); within a chunk the points run in
	  order. The start for every point is then fixed, and the output is identical
	  for any thread count or schedule.
	*/
	{
		int32_t nThreads = omp_get_max_threads();
		int32_t nChunks = (tiePoints->npts + TIEPOINTCHUNK - 1) / TIEPOINTCHUNK;
		int32_t iChunk;
		inputImageStructure *localImgs = (inputImageStructure *)malloc((size_t)nThreads * sizeof(inputImageStructure));
		if (localImgs == NULL)
			error("computeTiePoints: malloc failed for per-thread image copies\n");
		for (j = 0; j < nThreads; j++)
			localImgs[j] = *inputImage;
#pragma omp parallel \
	private(i, iChunk, lat, lon, xx, yy, sMask, deltaZ, vz1, zWGS, zSp, RElip, range, azimuth, azOrigin)
	{
	inputImageStructure *myImg = &localImgs[omp_get_thread_num()];
#pragma omp for schedule(dynamic, 1)
	for (iChunk = 0; iChunk < nChunks; iChunk++)
	{
	myImg->lastTime = -LARGEINT; /* outside sTime..eTime: first point starts from centre time */
	for (i = iChunk * TIEPOINTCHUNK; i < min((iChunk + 1) * TIEPOINTCHUNK, tiePoints->npts); i++)
	{
		lat = tiePoints->lat[i];
		lon = tiePoints->lon[i];
		/* This test should only be needed for lltora cases, since ll values should be valid for tiepoints */
		if ((lat > -90.001 && lat < 90.001) || 1 == 1)
		{
			/* Deal with shelf mask if necessary && shelfMaskFlag != -1 */
			if (shelfMask != NULL)
			{
				lltoxy1(lat, lon, &xx, &yy, Rotation, tiePoints->stdLat);
				sMask = getShelfMask(shelfMask, xx, yy);
				if (sMask == SHELF)
				{
					deltaZ = interpTideDiff(xx, yy, inputImage->tideDiff);
					
					if (deltaZ > MINELEVATION)
					{
						vz1 = deltaZ / ((double)tiePoints->nDays) * 365.25;
						tiePoints->vz[i] += vz1;
					}
				}
			}
			/*
			  Get height for given lat/lon, with spherical correction
			*/
			zWGS = tiePoints->z[i];

			if (tiePoints->imageCoords == FALSE)
			{
				llToImageNew(lat, lon, zWGS, &range, &azimuth, myImg);
				RElip = earthRadius(lat * DTOR, EMINOR, EMAJOR) * KMTOM + tiePoints->z[i];
				/*  Elevation referenced to a sphere */
				zSp = RElip - Re;
				//if((i +11) % 818 == 0) fprintf(stderr, "lat,lon, z  %lf %lf %lf %lf %lf %lf\n", tiePoints->lat[i], tiePoints->lon[i], tiePoints->z[i], range, azimuth, zWGS);
				//if(azimuth > 690) error("stop");
			}
			else
			{
				/* Obsolete */
				range = lat;		   /* Lat variable = range */
				zSp = tiePoints->z[i]; /* Assume r/a/z already spherical */
			}
			tiePoints->z[i] = zSp;
			tiePoints->r[i] = range;
			if (tiePoints->imageCoords == FALSE)
			{
				tiePoints->a[i] = azimuth;
			}
			else
			{
				/* ??? */
				error("check computeTiePoints - for this obsolete line");
				azimuth = inputImage->azimuthSize - lon;
				tiePoints->a[i] = lon;
			}
			azOrigin = -0.5 * inputImage->azimuthPixelSize * inputImage->azimuthSize; // inputImage->nAzimuthLooks;
			tiePoints->x[i] = azOrigin + inputImage->azimuthPixelSize * azimuth;
		}
		else
		{
			tiePoints->a[i] = -LARGEINT;
			tiePoints->r[i] = -LARGEINT;
		}
	}
	}
	}
		free(localImgs);
	}
	lastTime = clock();
	fprintf(stderr, "initTime %f convertTime %f totalTime %f\n", 
		(double)(initTime - startTime) / CLOCKS_PER_SEC, 
		(double)(lastTime - initTime) / CLOCKS_PER_SEC, 
		(double)(lastTime - startTime) / CLOCKS_PER_SEC);
	return;
}
