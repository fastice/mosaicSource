#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "gdal.h"
#include "cpl_conv.h"
#include "mosaicSource/common/common.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include "smoothRadius.h"

/*
  smoothradius -- per-pixel smoothing-radius map from merged range/azimuth offsets.

  Extracted from simoffsets.computeSmoothRadiusMap(), where it dominated the ROFF stage:
  measured at ~873 s per frame on a 3813 x 3502 grid, roughly the whole of the per-frame
  ROFF budget.  The work is 50 sweeps of the full image and cannot be shortened (some pixels
  legitimately reach the cap), so the gain is from compiled code and running the independent
  sweeps in parallel.

  Usage:
    smoothradius [-maxRadiusR n] [-maxRadiusA n] [-nIter n] [-ompThreads n]
                 drFile daFile tolDrFile tolDaFile outputBase

  dr/da/tolDr/tolDa are single-band float rasters (GeoTIFF or VRT).  Writes
  <outputBase>.smr.tif and <outputBase>.smr.vrt, GDT_Byte, 0 = no data.
*/

/*  Required by the common/ objects this links against (llToImageNew, readXYDEM,
    parseInputFile).  smoothradius does no geocoding -- it works purely in the offsets'
    own pixel grid -- so these keep the same defaults the other tools start from. */
int32_t HemiSphere = NORTH;
double Rotation = 45.;
double SLat = -91.0;
int32_t llConserveMem = 999;

static void usage()
{
	fprintf(stderr,
			"\nsmoothradius -- per-pixel smoothing-radius map from merged offsets\n\n"
			"Usage: smoothradius [-maxRadiusR n] [-maxRadiusA n] [-nIter n] [-ompThreads n]\n"
			"                    dr da tolDr tolDa outputBase\n\n"
			"  dr, da           merged range/azimuth offsets (float, tif or vrt)\n"
			"  tolDr, tolDa     per-pixel tolerances, same units as dr/da\n"
			"  outputBase       writes <outputBase>.smr.tif and .smr.vrt (byte, 0 = no data)\n\n"
			"  -maxRadiusR n    range half-width cap in pixels [50]\n"
			"  -maxRadiusA n    azimuth half-width cap in pixels [50]\n"
			"                   equal caps give an isotropic sweep, which is what a\n"
			"                   square-pixel offsets grid wants\n"
			"  -nIter n         box passes per sweep step [3]\n"
			"  -ompThreads n    parallel sweeps [2]; memory is ~0.9 GB + 0.3 GB per thread\n"
			"                   at a 13 Mpixel grid, so raise it with that in mind\n\n");
	exit(1);
}

/* Read a single-band float raster (tif or vrt) into a flat nr x na buffer. */
static float *readFloatRaster(char *file, int32_t *nr, int32_t *na)
{
	float *flat;
	int32_t xSize, ySize, dataType;
	char vrtBuf[2048];
	char *name = checkForVrt(file, vrtBuf);
	/*  checkForVrt returns NULL when there is no .vrt companion; a plain tif is fine, so
	    fall back to the name as given rather than handing GDALOpen a NULL. */
	if (name == NULL)
	{
		name = file;
	}
	/*  readRasterVRT returns a FLAT row-major buffer (allocData -> CPLMalloc), not the
	    float ** row-pointer layout mallocImage() produces.  It also clamps yMax down to
	    ySize-1, so a full read asks for more rows than exist rather than passing a
	    sentinel; -1 would read nothing. */
	flat = (float *)readRasterVRT(name, 1, &xSize, &ySize, &dataType, NULL, NULL, 0, INT32_MAX);
	if (flat == NULL)
	{
		error("smoothradius: could not read %s\n", file);
	}
	if (dataType != GDT_Float32)
	{
		error("smoothradius: %s is not Float32 (type %i)\n", file, dataType);
	}
	if (*nr < 0)
	{
		*nr = xSize;
		*na = ySize;
	}
	else if (xSize != *nr || ySize != *na)
	{
		error("smoothradius: %s is %i x %i, expected %i x %i\n", file, xSize, ySize, *nr, *na);
	}
	return flat;
}

int main(int argc, char **argv)
{
	smrParams p;
	char *drFile, *daFile, *tolDrFile, *tolDaFile, *outBase;
	char tifBuf[2048], vrtBuf[2048];
	char *tifFile, *vrtFile;
	const char *tifFiles[1];
	float noDataArr[1] = {0.0f};
	int32_t i;

	GDALAllRegister();
	p.nr = -1; p.na = -1;
	p.maxRadiusR = 50;
	p.maxRadiusA = 50;
	p.nIter = 3;
	p.ompThreads = 2;

	for (i = 1; i < argc - 5; i++)
	{
		char *arg = argv[i];
		if (strstr(arg, "maxRadiusR") != NULL) { sscanf(argv[++i], "%i", &p.maxRadiusR); }
		else if (strstr(arg, "maxRadiusA") != NULL) { sscanf(argv[++i], "%i", &p.maxRadiusA); }
		else if (strstr(arg, "nIter") != NULL) { sscanf(argv[++i], "%i", &p.nIter); }
		else if (strstr(arg, "ompThreads") != NULL) { sscanf(argv[++i], "%i", &p.ompThreads); }
		else { usage(); }
	}
	if (argc - i != 5)
	{
		usage();
	}
	drFile = argv[i]; daFile = argv[i + 1];
	tolDrFile = argv[i + 2]; tolDaFile = argv[i + 3];
	outBase = argv[i + 4];

	if (p.maxRadiusR > 255 || p.maxRadiusA > 255)
	{
		fprintf(stderr, "smoothradius: WARNING radius cap exceeds the byte output range, "
						"clamping to 255\n");
		if (p.maxRadiusR > 255) p.maxRadiusR = 255;
		if (p.maxRadiusA > 255) p.maxRadiusA = 255;
	}
	if (p.maxRadiusR < 1 || p.maxRadiusA < 1 || p.nIter < 1 || p.ompThreads < 1)
	{
		error("smoothradius: radii, nIter and ompThreads must all be >= 1\n");
	}

	/* The first read fixes the dimensions; the rest must match. */
	p.dr = readFloatRaster(drFile, &p.nr, &p.na);
	p.da = readFloatRaster(daFile, &p.nr, &p.na);
	p.tolDr = readFloatRaster(tolDrFile, &p.nr, &p.na);
	p.tolDa = readFloatRaster(tolDaFile, &p.nr, &p.na);
	p.radius = (unsigned char *)malloc((size_t)p.nr * (size_t)p.na);
	if (p.radius == NULL)
	{
		error("smoothradius: out of memory for the output map\n");
	}
	fprintf(stderr, "smoothradius: %i x %i, maxRadius %i(range)/%i(azimuth), nIter %i, "
					"ompThreads %i\n", p.nr, p.na, p.maxRadiusR, p.maxRadiusA, p.nIter,
			p.ompThreads);

	computeSmoothRadiusOffsets(&p);

	tifFile = appendSuffix(outBase, ".smr.tif", tifBuf);
	writeFlatTiff(tifFile, p.radius, p.nr, p.na, GDT_Byte, 0.0f, NULL);
	tifFiles[0] = tifFile;
	vrtFile = appendSuffix(outBase, ".smr.vrt", vrtBuf);
	makeTiffVRT(vrtFile, tifFiles, 1, noDataArr, NULL);
	fprintf(stderr, "smoothradius: wrote %s\n", tifFile);
	return 0;
}
