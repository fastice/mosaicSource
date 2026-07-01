#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"

#define SMR_NODATA_THRESH -1.0e6

/*
   Separable box-filter average with half-widths (wa, wr), skipping no-data pixels
   (value <= SMR_NODATA_THRESH). O(1) per-pixel via a running sum, independent of window size.
*/
static void boxAverage(float **in, float **out, int32_t aSize, int32_t rSize,
                       int32_t wa, int32_t wr)
{
	int32_t i, j, addIdx, remIdx;
	float **tmp;
	double sum;
	int32_t count;

	tmp = (float **)malloc(aSize * sizeof(float *));
	for (i = 0; i < aSize; i++)
		tmp[i] = (float *)malloc(rSize * sizeof(float));

	/* Pass 1: running sum along range, within each azimuth row */
	for (i = 0; i < aSize; i++)
	{
		sum = 0.;
		count = 0;
		for (j = 0; j <= wr && j < rSize; j++)
			if (in[i][j] > SMR_NODATA_THRESH) { sum += in[i][j]; count++; }
		tmp[i][0] = (count > 0) ? (float)(sum / count) : in[i][0];
		for (j = 1; j < rSize; j++)
		{
			addIdx = j + wr;
			remIdx = j - wr - 1;
			if (addIdx < rSize && in[i][addIdx] > SMR_NODATA_THRESH) { sum += in[i][addIdx]; count++; }
			if (remIdx >= 0 && in[i][remIdx] > SMR_NODATA_THRESH) { sum -= in[i][remIdx]; count--; }
			tmp[i][j] = (count > 0) ? (float)(sum / count) : in[i][j];
		}
	}
	/* Pass 2: running sum along azimuth, within each range column */
	for (j = 0; j < rSize; j++)
	{
		sum = 0.;
		count = 0;
		for (i = 0; i <= wa && i < aSize; i++)
			if (tmp[i][j] > SMR_NODATA_THRESH) { sum += tmp[i][j]; count++; }
		out[0][j] = (count > 0) ? (float)(sum / count) : in[0][j];
		for (i = 1; i < aSize; i++)
		{
			addIdx = i + wa;
			remIdx = i - wa - 1;
			if (addIdx < aSize && tmp[addIdx][j] > SMR_NODATA_THRESH) { sum += tmp[addIdx][j]; count++; }
			if (remIdx >= 0 && tmp[remIdx][j] > SMR_NODATA_THRESH) { sum -= tmp[remIdx][j]; count--; }
			out[i][j] = (count > 0) ? (float)(sum / count) : in[i][j];
		}
	}
	for (i = 0; i < aSize; i++)
		free(tmp[i]);
	free(tmp);
}

/*
  Compute scene->radiusImage: per-pixel single-look azimuth-pixel half-width, the largest
  contiguous radius r (1..scene->maxSmoothRadius) for which repeatedly (scene->smoothNIter
  times, box->triangular->Gaussian-ish) box-filtering scene->image with half-width
  (r, round(r*azimuthPixelSize/rangePixelSize)) changes the pixel by no more than
  scene->toleranceImage at that pixel. First-violation cutoff: once a pixel fails at some r,
  its answer is locked at r-1 (the largest *contiguous* safe radius) and never revisited, even
  if a larger r would coincidentally land back inside tolerance.
*/
void computeSmoothRadiusMap(sceneStructure *scene)
{
	int32_t aSize = scene->aSize, rSize = scene->rSize;
	int32_t i, j, r, k, wr, nIter;
	float **smoothed, **current;
	char **locked;
	int64_t nLocked, nTotal;
	double pixRatio;

	scene->radiusImage = (unsigned char **)malloc(aSize * sizeof(unsigned char *));
	smoothed = (float **)malloc(aSize * sizeof(float *));
	current = (float **)malloc(aSize * sizeof(float *));
	locked = (char **)malloc(aSize * sizeof(char *));
	nLocked = 0;
	nTotal = (int64_t)aSize * (int64_t)rSize;
	for (i = 0; i < aSize; i++)
	{
		scene->radiusImage[i] = (unsigned char *)calloc(rSize, sizeof(unsigned char));
		smoothed[i] = (float *)malloc(rSize * sizeof(float));
		current[i] = (float *)malloc(rSize * sizeof(float));
		locked[i] = (char *)calloc(rSize, sizeof(char));
		for (j = 0; j < rSize; j++)
		{
			if (scene->image[i][j] <= SMR_NODATA_THRESH)
			{
				locked[i][j] = TRUE; /* no-data pixels: leave radius at 0 */
				nLocked++;
			}
		}
	}

	pixRatio = scene->I.azimuthPixelSize / scene->I.rangePixelSize;
	nIter = (scene->smoothNIter > 0) ? scene->smoothNIter : 1;

	for (r = 1; r <= scene->maxSmoothRadius && nLocked < nTotal; r++)
	{
		wr = (int32_t)(r * pixRatio + 0.5);
		if (wr < 0) wr = 0;
		/* nIter repeated box-filter passes starting fresh from the original image each step */
		for (i = 0; i < aSize; i++)
			for (j = 0; j < rSize; j++)
				current[i][j] = scene->image[i][j];
		for (k = 0; k < nIter; k++)
		{
			boxAverage(current, smoothed, aSize, rSize, r, wr);
			for (i = 0; i < aSize; i++)
				for (j = 0; j < rSize; j++)
					current[i][j] = smoothed[i][j];
		}
		for (i = 0; i < aSize; i++)
		{
			for (j = 0; j < rSize; j++)
			{
				if (locked[i][j])
					continue;
				if (fabsf(smoothed[i][j] - scene->image[i][j]) <= scene->toleranceImage[i][j])
				{
					scene->radiusImage[i][j] = (unsigned char)r;
				}
				else
				{
					locked[i][j] = TRUE;
					nLocked++;
				}
			}
		}
	}

	for (i = 0; i < aSize; i++)
	{
		free(smoothed[i]);
		free(current[i]);
		free(locked[i]);
	}
	free(smoothed);
	free(current);
	free(locked);
}
