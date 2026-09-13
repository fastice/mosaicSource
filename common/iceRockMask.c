#include "stdio.h"
#include "string.h"
#include "math.h"
#include <stdlib.h>
#include "gdal.h"
#include "cpl_conv.h"
#include "ogr_srs_api.h"
#include "mosaicSource/common/common.h"

/*
  Sample a GIMP-style ice/rock/water mask (0=water, 1=rock, 2=ice) at a list of
  lat/lon points.

  Used by tiepoints -iceOnly to drop tie points that do not sit on ice.  Motivation:
  unwrapping across the rock/ice margin can introduce a 2*pi ambiguity; if tie
  points on rock then anchor the baseline fit, that ambiguity is pushed onto the
  whole ice sheet.  Sentinel-1 processing avoids this by using phase on ice only.

  The transform is built from the MASK FILE'S OWN projection via PROJ rather than
  from lltoxy1().  lltoxy1 hard-codes a standard parallel of 70 and takes a
  separate rotation argument; it is close to EPSG:3413 but not guaranteed
  identical, and any small mismatch would misclassify points exactly at the
  margin -- which is the only place the classification matters here.

  Returns a malloc'd array of length n; IRM_OUTSIDE where a point falls off the
  mask.  Returns NULL (with a message) if the mask cannot be read, so the caller
  can decide whether that is fatal.
*/

unsigned char *sampleIceRockMask(char *maskFile, double *lat, double *lon, int32_t n)
{
	GDALDatasetH hDS;
	GDALRasterBandH hBand;
	OGRSpatialReferenceH srcSRS = NULL, dstSRS = NULL;
	OGRCoordinateTransformationH ct = NULL;
	double gt[6], invGt[6];
	unsigned char *out = NULL;
	const char *proj;
	int32_t i, xSize, ySize, nIce = 0, nRock = 0, nWater = 0, nOut = 0;

	if (maskFile == NULL || n <= 0)
		return (NULL);
	GDALAllRegister();
	hDS = GDALOpen(maskFile, GA_ReadOnly);
	if (hDS == NULL)
	{
		fprintf(stderr, "sampleIceRockMask: cannot open %s\n", maskFile);
		return (NULL);
	}
	if (GDALGetGeoTransform(hDS, gt) != CE_None || !GDALInvGeoTransform(gt, invGt))
	{
		fprintf(stderr, "sampleIceRockMask: %s has no usable geotransform\n", maskFile);
		GDALClose(hDS);
		return (NULL);
	}
	xSize = GDALGetRasterXSize(hDS);
	ySize = GDALGetRasterYSize(hDS);
	hBand = GDALGetRasterBand(hDS, 1);
	/* Build lat/lon -> mask CRS transform from the file's own projection */
	proj = GDALGetProjectionRef(hDS);
	if (proj == NULL || strlen(proj) == 0)
	{
		fprintf(stderr, "sampleIceRockMask: %s has no projection\n", maskFile);
		GDALClose(hDS);
		return (NULL);
	}
	srcSRS = OSRNewSpatialReference(NULL);
	dstSRS = OSRNewSpatialReference(NULL);
	OSRSetFromUserInput(srcSRS, "EPSG:4326");
	OSRImportFromWkt(dstSRS, (char **)&proj);
	/* Force lon,lat ordering so the call below is unambiguous across PROJ versions */
	OSRSetAxisMappingStrategy(srcSRS, OAMS_TRADITIONAL_GIS_ORDER);
	OSRSetAxisMappingStrategy(dstSRS, OAMS_TRADITIONAL_GIS_ORDER);
	ct = OCTNewCoordinateTransformation(srcSRS, dstSRS);
	if (ct == NULL)
	{
		fprintf(stderr, "sampleIceRockMask: cannot build transform to the mask CRS\n");
		OSRDestroySpatialReference(srcSRS);
		OSRDestroySpatialReference(dstSRS);
		GDALClose(hDS);
		return (NULL);
	}
	/* Per-point 1x1 reads; GDAL's block cache makes this cheap because tie
	   points are spatially clustered.  Reading a bounding window instead would
	   be up to the full 16620x30000 mask for a long virtual frame. */
	GDALSetCacheMax64((GIntBig)256 * 1024 * 1024);
	out = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
	for (i = 0; i < n; i++)
	{
		double px = lon[i], py = lat[i], pz = 0.0;
		double fx, fy;
		int32_t ix, iy;
		unsigned char v;
		if (!OCTTransform(ct, 1, &px, &py, &pz))
		{
			out[i] = IRM_OUTSIDE;
			nOut++;
			continue;
		}
		fx = invGt[0] + invGt[1] * px + invGt[2] * py;
		fy = invGt[3] + invGt[4] * px + invGt[5] * py;
		ix = (int32_t)floor(fx);
		iy = (int32_t)floor(fy);
		if (ix < 0 || iy < 0 || ix >= xSize || iy >= ySize)
		{
			out[i] = IRM_OUTSIDE;
			nOut++;
			continue;
		}
		if (GDALRasterIO(hBand, GF_Read, ix, iy, 1, 1, &v, 1, 1, GDT_Byte, 0, 0) != CE_None)
		{
			out[i] = IRM_OUTSIDE;
			nOut++;
			continue;
		}
		out[i] = v;
		if (v == IRM_ICE)
			nIce++;
		else if (v == IRM_ROCK)
			nRock++;
		else
			nWater++;
	}
	fprintf(stderr, "sampleIceRockMask: %s -- ice %i, rock %i, water %i, outside %i (of %i)\n",
			maskFile, nIce, nRock, nWater, nOut, n);
	OCTDestroyCoordinateTransformation(ct);
	OSRDestroySpatialReference(srcSRS);
	OSRDestroySpatialReference(dstSRS);
	GDALClose(hDS);
	return (out);
}

/*
  Drop every tie point that is not on ice.  Compacts all per-point arrays in
  place and updates npts.  Returns the number discarded.
*/
int32_t keepIceTiePoints(tiePointsStructure *tiePoints, char *maskFile)
{
	unsigned char *m;
	int32_t i, k = 0, n = tiePoints->npts;
	int32_t hasIon = (tiePoints->hasIonosphere == TRUE && tiePoints->ionPhase != NULL);
	int32_t hasSq = (tiePoints->hasSquintSolution == TRUE && tiePoints->phaseSquint != NULL);

	m = sampleIceRockMask(maskFile, tiePoints->lat, tiePoints->lon, n);
	if (m == NULL)
		error("keepIceTiePoints: could not sample ice/rock mask %s", maskFile);
	for (i = 0; i < n; i++)
	{
		if (m[i] != IRM_ICE)
			continue;
		tiePoints->lat[k] = tiePoints->lat[i];
		tiePoints->lon[k] = tiePoints->lon[i];
		tiePoints->x[k] = tiePoints->x[i];
		tiePoints->y[k] = tiePoints->y[i];
		tiePoints->z[k] = tiePoints->z[i];
		tiePoints->r[k] = tiePoints->r[i];
		tiePoints->a[k] = tiePoints->a[i];
		tiePoints->phase[k] = tiePoints->phase[i];
		tiePoints->delta[k] = tiePoints->delta[i];
		tiePoints->vx[k] = tiePoints->vx[i];
		tiePoints->vy[k] = tiePoints->vy[i];
		tiePoints->vz[k] = tiePoints->vz[i];
		tiePoints->vyra[k] = tiePoints->vyra[i];
		tiePoints->weight[k] = tiePoints->weight[i];
		tiePoints->bsq[k] = tiePoints->bsq[i];
		if (hasIon)
			tiePoints->ionPhase[k] = tiePoints->ionPhase[i];
		if (hasSq)
			tiePoints->phaseSquint[k] = tiePoints->phaseSquint[i];
		k++;
	}
	free(m);
	fprintf(stderr, "\033[1;33mkeepIceTiePoints: kept %i of %i tie points on ice (dropped %i)\033[0m\n",
			k, n, n - k);
	tiePoints->npts = k;
	return (n - k);
}
