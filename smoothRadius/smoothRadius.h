/*
  smoothradius -- per-pixel smoothing-radius map from merged range/azimuth offsets.
*/
#ifndef SMOOTHRADIUS_H
#define SMOOTHRADIUS_H
#include <stdint.h>

/* No-data test for the offset rasters, matching simoffsets' own convention (-2.e9 fill). */
#define SMR_NODATA_THRESH -1.0e6f

typedef struct smrParamsType
{
	int32_t nr, na;          /* raster dimensions, nr fastest */
	int32_t maxRadiusR;      /* range half-width cap, single-look pixels */
	int32_t maxRadiusA;      /* azimuth half-width cap */
	int32_t nIter;           /* repeated box passes per sweep step */
	int32_t ompThreads;
	float *dr, *da;          /* merged offsets */
	float *tolDr, *tolDa;    /* per-pixel tolerances, same units as dr/da */
	unsigned char *radius;   /* output, 0 = no data */
} smrParams;

void computeSmoothRadiusOffsets(smrParams *p);

#endif
