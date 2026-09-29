#include "stdio.h"
#include "string.h"
#include <math.h>
#include "mosaicSource/common/common.h"
#include "geomosaic.h"
#include "cRecipes/nrutil.h"
#include <stdlib.h>
#include <omp.h>

/* from lsmosaic, but only need this */
double xyscale(double latctr, int32_t proj);
/*
************************ Geocode image InSAR DEM. **************************
*/
int32_t sepAscDesc;
static void readDEMInput(inputImageStructure *inputImage, char *insarDEMFile);
static float interpolatePowerInputImage(inputImageStructure inputImage, double range, double azimuth);
static float polypat(float theta);
static float polyALOS(float theta);
static void smoothImage(inputImageStructure *inputImage, int32_t smoothL);
static void rsatFinalCal(float **image, int32_t xSize, int32_t ySize);
static void logSigma(float **image, int32_t xSize, int32_t ySize);
static void mallocTmpBuffers(float ***imageTmp, float ***scaleTmp, outputImageStructure *outputImage,
							 float ***psiBuf, float ***psiBufTmp, float ***gBuf, float ***gBufTmp);
static void getGeoMosaicImage(inputImageStructure *inputImage, int32_t *imageDate, int32_t smoothL, int32_t yMin, int32_t yMax);
float applyCorrections(float *value, inputImageStructure *inputImage, double range, double azimuth, double h);
static void geoMosaicScaling(inputImageStructure *inputImage, float **image, float **imageTmp, float **psiBuf,
							 float **psiBufTmp, float **gBuf, float **gBufTmp,
							 float **scale, float **scaleTmp, void *dem, outputImageStructure *outputImage, int32_t orbitPriority,
							 int32_t imageDate, int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax);
static void finalReScale(outputImageStructure *outputImage, float **image, float **scale, int32_t orbitPriority,
						 float **imageAll, float **scaleAll, float **psiBuf, float **gBuf,
						 float **psiBufAll, float **gBufAll);
static float **mallocPlane(outputImageStructure *outputImage, float initValue);
static void applySelection(float **imageTmp, float **scaleTmp, unsigned char **selTmp,
						   inputImageStructure *inputImage, int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax);
static void applyRamp(float **scaleTmp, unsigned char **selTmp,
					  int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax);
static void initImageBuffers(outputImageStructure *outputImage, float **scaleTmp, int32_t orbitPriority, float **psiBuf,
							 float **psiBufTmp, float **gBuf, float **gBufTmp);
static void computeScaleFast(float **inImage, float **scale, int32_t ySize, int32_t xSize, float fl, float weight,
							 double minVal, int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax);

/*
  Compute the lenth of an edge in 3-D
*/
static double edgeLength(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
}

/*
   Given a triangle specified by P1(x1, y1, z1),P2(x2, y2, z2),P3(x3, y3, z3),compute area
*/
static double threeDTriangleArea(double x1, double y1, double z1,
								 double x2, double y2, double z2, double x3, double y3, double z3)
{
	double p12, p23, p31, h1, A;
	p12 = edgeLength(x1, y1, z1, x2, y2, z2);
	p23 = edgeLength(x2, y2, z2, x3, y3, z3);
	p31 = edgeLength(x3, y3, z3, x1, y1, z1);
	h1 = 0.5 * (p12 + p23 + p31);
	A = sqrt(h1 * (h1 - p12) * (h1 - p23) * (h1 - p31));
	return A;
}

/*
   Normalize a vector a V/|V|
*/
static void normVector(double x0, double y0, double z0, double *nx, double *ny, double *nz)
{
	double den;
	den = sqrt(pow(x0, 2.) + pow(y0, 2.) + pow(z0, 2.));
	*nx = x0 / den;
	*ny = y0 / den;
	*nz = z0 / den;
}

/*
  find the normal to three points.
*/
static void threePtNormal(double x1, double y1, double z1,
						  double x2, double y2, double z2, double x3, double y3, double z3,
						  double *nx, double *ny, double *nz)
{
	cross(x2 - x1, y2 - y1, z2 - z1, x3 - x1, y3 - y1, z3 - z1, nx, ny, nz);
	normVector(*nx, *ny, *nz, nx, ny, nz);
}

/*
   Get normal to satellite (look direction) at a point (x0, y0, x0).
   return
   nlx, nly, nlz satellite look direction.
   nsx, nsy, nsz - normalized velocity vector
*/
static void satNorm(double x0, double y0, double z0,
					double azimuth, inputImageStructure *inputImage,
					double *nlx, double *nly, double *nlz, double *nsx, double *nsy, double *nsz)
{
	conversionDataStructure *cp;
	double vsx, vsy, vsz;
	double xs, ys, zs;
	double myTime;
	double den;
	cp = &(inputImage->cpAll);
	myTime = (azimuth * inputImage->nAzimuthLooks) / cp->prf + cp->sTime;
	/* Get the state vector for that time */
	getState(myTime, inputImage, &xs, &ys, &zs, &vsx, &vsy, &vsz);
	/* Norm vector from the point to the satellite */
	normVector(x0 - xs, y0 - ys, z0 - zs, nlx, nly, nlz);
	/* Norm velocity vector */
	normVector(vsx, vsy, vsz, nsx, nsy, nsz);
}

/*
   Return the vector normal to the Earth (nx, ny, nz) at a point (x0, y0, z0)
*/
static void earthNorm(double x0, double y0, double z0, double *nx, double *ny, double *nz)
{
	double den;
	normVector(x0, y0, z0, nx, ny, nz);
}

/*
   get the vector normal to the range as n_earth X sat_trajectory
*/
static void rangePlane(double nex, double ney, double nez, double nsx, double nsy, double nsz, double *nrx, double *nry, double *nrz)
{
	double den;
	cross(nex, ney, nez, nsx, nsy, nsz, nrx, nry, nrz);
	normVector(*nrx, *nry, *nrz, nrx, nry, nrz);
}

/*
  For four points, compute the 3D area as the sum of two triangles.
*/
static double threeDArea(double x1[4], double y1[4], double z1[4])
{
	double A1, A2;
	A1 = threeDTriangleArea(x1[0], y1[0], z1[0], x1[1], y1[1], z1[1], x1[3], y1[3], z1[3]);
	A2 = threeDTriangleArea(x1[1], y1[1], z1[1], x1[2], y1[2], z1[2], x1[3], y1[3], z1[3]);
	return (A1 + A2);
}

static double nrx, nry, nrz, nsx, nsy, nsz, nlx, nly, nlz;
#pragma omp threadprivate(nrx, nry, nrz, nsx, nsy, nsz, nlx, nly, nlz)

/*
  For PS (x,y) point, corresponding to azimuth, compute projection of the pixel around it onto the beta, and gamma
  planes. Return Ab, Ag
*/
void AbAg(double x, double y, double azimuth, inputImageStructure *inputImage,
				 outputImageStructure *outputImage, void *dem, double *Ab, double *Ag, 
				 double *shadow, int32_t recycle)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	double xc[4], yc[4], z[4];
	double latc[4], lonc[4];
	double xi[4] = {-0.5, 0.5, 0.5, -0.5}; /* Use +/- for a smoother result */
	double yi[4] = {-0.5, -0.5, 0.5, 0.5};
	double x1[4], y1[4], z1[4];
	double x1g[4], y1g[4], z1g[4];
	double x1b[4], y1b[4], z1b[4];
	double xp[4], yp[4], zp4[4];
	double nex, ney, nez, nex1, ney1, nez1, nex2, ney2, nez2;
	double vdotn, rdotn;
	double x0, y0, z0;
	int32_t i;
	/* Map pixel corners to ECEF */
	x0 = 0.0;
	y0 = 0.0;
	z0 = 0.0;
	for (i = 0; i < 4; i++)
	{
		xc[i] = x + xi[i] * outputImage->deltaX * MTOKM;
		yc[i] = y + yi[i] * outputImage->deltaY * MTOKM;
		xytoll1(xc[i], yc[i], HemiSphere, &(latc[i]), &(lonc[i]), Rotation, outputImage->slat);
		z[i] = getXYHeight(latc[i], lonc[i], dem, inputImage->cpAll.Re, ELLIPSOIDAL);
		llToECEF(latc[i], lonc[i], z[i], &(x1[i]), &(y1[i]), &(z1[i]));
		/* Avg position */
		x0 += 0.25 * x1[i];
		y0 += 0.25 * y1[i];
		z0 += 0.25 * z1[i];
	}
	/* Normal to look direction  - toward surface  - gamma plane */
	if (recycle == FALSE)
	{
		satNorm(x0, y0, z0, azimuth, inputImage, &nlx, &nly, &nlz, &nsx, &nsy, &nsz);
		/* Normal for plane that contains sat and look vector - beta plane */
		rangePlane(nsx, nsy, nsz, nlx, nly, nlz, &nrx, &nry, &nrz);
	}
	/* Surface normal for each of two triangles forming a pixel */
	if (*shadow > 0)
	{
		threePtNormal(x1[0], y1[0], z1[0], x1[1], y1[1], z1[1], x1[2], y1[2], z1[2], &nex1, &ney1, &nez1);
		threePtNormal(x1[0], y1[0], z1[0], x1[2], y1[2], z1[2], x1[3], y1[3], z1[3], &nex2, &ney2, &nez2);
		/* Shadow check */
		*shadow = min(dot(nex1, ney1, nez1, -nlx, -nly, -nlz), dot(nex2, ney2, nez2, -nlx, -nly, -nlz));
	}
	/* Project points on to gamma and beta planes */
	for (i = 0; i < 4; i++)
	{
		/*  corner points dotted with line of sight vector */
		vdotn = dot(x1[i] - x0, y1[i] - y0, z1[i] - z0, nlx, nly, nlz);
		/* Points projected to line of sight plane */
		x1g[i] = x1[i] - vdotn * nlx;
		y1g[i] = y1[i] - vdotn * nly;
		z1g[i] = z1[i] - vdotn * nlz;
		/* project points on to plane with line of sight and velocity vector */
		rdotn = dot(x1[i] - x0, y1[i] - y0, z1[i] - z0, nrx, nry, nrz);
		x1b[i] = x1[i] - rdotn * nrx;
		y1b[i] = y1[i] - rdotn * nry;
		z1b[i] = z1[i] - rdotn * nrz;
	}
	/* normalization - removed since these are ratioed*/
	*Ag = threeDArea(x1g, y1g, z1g);
	*Ab = threeDArea(x1b, y1b, z1b);
}

/*
   Compute area of overlap for two boxes
*/
static double boxOverlap(double range0, double azimuth0, double range1, double azimuth1, double dx)
{
	double r01, a01, r02, a02;
	double r11, a11, r12, a12;
	double dr, da; /*, dr1,dr2, da1, da2;*/
	/* Compute box corners */
	r01 = range0 - dx;
	r02 = range0 + dx;
	a01 = azimuth0 - dx;
	a02 = azimuth0 + dx;
	r11 = range1 - dx;
	r12 = range1 + dx;
	a11 = azimuth1 - dx;
	a12 = azimuth1 + dx;
	/* No overlap cases, so return 0 */
	if (r12 < r01 || r11 > r02)
		return 0.0;
	if (a12 < a01 || a11 > a02)
		return 0.0;
	dr = min(r12 - r01, r02 - r11);
	da = min(a12 - a01, a02 - a11);
	return max(dr * da / (4 * dx * dx), 0); /* never should compute a negative, but just in case */
}

/*
   This routine is use for RTC  to compute area around a point.
   Basically loop around a point agregate areas that maps into a particular range
   pixel. Uses a fairly lazy way to weight aveverage, which relies on over laps of pixels
   assuming they map squarely into R/D space, which they often will not. Should still
   provide a quick-and-dirty relative weighting that is fine given all of the other uncertainties,
   smoothing etc.
*/
static double areaAboutXY(double range, double azimuth, double x, double y, inputImageStructure *inputImage,
						  outputImageStructure *outputImage, void *dem, double value, double *test, int32_t recycle, double *Ab, double *Ag)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	double shadow, Ab1, Ag1, rat;
	double x1, y1, lat1, lon1, hWGS, range1, azimuth1;
	double fracOverlap;
	int32_t nx = 3, ny = 2;
	double dx;
	int32_t i, j;
	shadow = 1;
	/* Compute areas about about a point */
	AbAg(x, y, azimuth, inputImage, outputImage, dem, Ab, Ag, &shadow, recycle);
	if (shadow < -0.001)
		return (-1.);
	rat = *Ab / *Ag;
	/* If dark region, or high ratio assume rat is representative */
	if (10. * log10(value) < -2. || rat > 0.65)
		return rat;
	/* This part is quite an inefficient way to bin ground samples contributing to the same range bin */
	dx = 0.5;
	/* Each neighbour recomputes look/plane normals from its own azimuth position.
	   To revert to recycling the center-pixel normals (faster but less accurate over
	   steep terrain), replace FALSE with recycle in the AbAg call below and restore:
	       recycle = TRUE;
	   here before the loop. */
	Ab1 = 0.0;
	for (i = -nx; i <= nx; i++)
	{
		x1 = x + i * outputImage->deltaX * MTOKM;
		for (j = -ny; j <= ny; j++)
		{
			if (!((i == 0) && (j == 0)))
			{ /* Don't redo (0,0) case */
				y1 = y + j * outputImage->deltaY * MTOKM;
				xyToLLProj(x1, y1, &lat1, &lon1, &(outputImage->proj));
				hWGS = getXYHeight(lat1, lon1, dem, inputImage->cpAll.Re, ELLIPSOIDAL);
				llToImageNew(lat1, lon1, hWGS, &range1, &azimuth1, inputImage);
				/*
				   Determine how much a range azimuth pixel centered at (range1, azimuth1) would overlap
				   with one at (range, azimuth)
				*/
				fracOverlap = boxOverlap(range, azimuth, range1, azimuth1, dx);
				/* If overlap, then calculate the area projected onto the gamma and beta planes */
				if (fracOverlap > 0)
				{
					shadow = -1;
					AbAg(x1, y1, azimuth1, inputImage, outputImage, dem, &Ab1, &Ag1, &shadow, FALSE);
					*Ab += Ab1 * fracOverlap;
					*Ag += Ag1 * fracOverlap;
				}
			}
		}
	}
	return (*Ab) / (*Ag);
}

/* Compute indices for producing and oversampled result */
static void oversampledxdy(int32_t os, double dx[3], double dy[3])
{
	if (os == 1)
	{
		dx[0] = 0;
		dy[0] = 0;
	}
	else if (os == 2)
	{
		dx[0] = -0.5;
		dy[0] = -0.5;
		dx[1] = 0.5;
		dy[1] = 0.5;
	}
	else if (os == 3)
	{
		dx[0] = -1;
		dy[0] = -1;
		dx[1] = 0;
		dy[1] = 0;
		dx[2] = 1;
		dy[2] = 1;
	}
	else
		error("Invalid oversample - should be between 1 and 3");
}

/*
  Make mosaic of insar and other dems.
*/
/*
  ---- Near/far range incidence selection -------------------------------------------------------

  A cheap geometry-only first pass builds a COARSE buffer (stride angleStride output pixels) of
  the per-pixel minimum (-nearRange) or maximum (-farRange) ellipsoidal incidence angle over all
  inputs. Pass 2 then keeps only inputs within angleTolerance of it. Incidence varies ~0.05 deg/km
  and is monotone across a frame, so a stride of 10 costs ~0.03 deg of quantisation at 100 m.

  Cells are addressed by ABSOLUTE output-lattice index, so a tiled run and an untiled run put the
  cell boundaries in the same places. Pass 1 evaluates each cell at its centre even when that
  falls outside the current tile, so a cell's contents do not depend on the tiling either.

  Pass 1 cannot know where an input has valid DATA (only where it has geometry), so the filter can
  empty a pixel that the plain mosaic would have filled. That is handled downstream by the
  unfiltered fallback accumulator in makeGeoMosaic, not here.
*/

/* Absolute output-lattice pixel index of the grid's first row/column. */
int32_t incCoarseIndex(double origin, double delta, int32_t stride)
{
	return (int32_t)floor(origin / delta + 0.5);
}

/* Floor division; C division truncates toward zero, which would fold cells together at 0. */
static int32_t floorDiv(int32_t a, int32_t b)
{
	return (a >= 0) ? (a / b) : -(((-a) + b - 1) / b);
}

/* Local cell index of output pixel k along one axis (abs0 = absolute index of pixel 0). */
static int32_t cellOf(int32_t abs0, int32_t k, int32_t stride)
{
	return floorDiv(abs0 + k, stride) - floorDiv(abs0, stride);
}

/* Output pixel at the centre of local cell c; may fall outside the grid, which is intended. */
static int32_t cellCentrePixel(int32_t abs0, int32_t c, int32_t stride)
{
	return (floorDiv(abs0, stride) + c) * stride + stride / 2 - abs0;
}

void incCellSpan(incBuffer *incBuf, int32_t isX, int32_t kMin, int32_t kMax, int32_t *cMin, int32_t *cMax)
{
	int32_t abs0 = isX ? incBuf->j0 : incBuf->i0;
	*cMin = cellOf(abs0, kMin, incBuf->stride);
	*cMax = cellOf(abs0, kMax, incBuf->stride);
}

int32_t incCellCentre(incBuffer *incBuf, int32_t isX, int32_t c)
{
	int32_t abs0 = isX ? incBuf->j0 : incBuf->i0;
	return cellCentrePixel(abs0, c, incBuf->stride);
}

static incBuffer *allocIncBuffer(outputImageStructure *outputImage, int32_t stride)
{
	incBuffer *incBuf;
	float *buf;
	int32_t i, j;
	incBuf = (incBuffer *)malloc(sizeof(incBuffer));
	if (incBuf == NULL)
	{
		error("allocIncBuffer: malloc failed");
	}
	incBuf->stride = stride;
	incBuf->j0 = incCoarseIndex(outputImage->originX, outputImage->deltaX, stride);
	incBuf->i0 = incCoarseIndex(outputImage->originY, outputImage->deltaY, stride);
	incBuf->nX = cellOf(incBuf->j0, outputImage->xSize - 1, stride) + 1;
	incBuf->nY = cellOf(incBuf->i0, outputImage->ySize - 1, stride) + 1;
	buf = (float *)malloc((size_t)incBuf->nX * incBuf->nY * sizeof(float));
	incBuf->inc = (float **)malloc((size_t)incBuf->nY * sizeof(float *));
	if (buf == NULL || incBuf->inc == NULL)
	{
		error("allocIncBuffer: malloc failed for %i x %i cells", incBuf->nX, incBuf->nY);
	}
	for (i = 0; i < incBuf->nY; i++)
	{
		incBuf->inc[i] = &(buf[(size_t)i * incBuf->nX]);
		for (j = 0; j < incBuf->nX; j++)
		{
			incBuf->inc[i][j] = INCUNSET;
		}
	}
	fprintf(stderr, "rangeSelect: %i x %i incidence cells, stride %i\n", incBuf->nX, incBuf->nY, stride);
	return incBuf;
}

/*
  Fold psi at output pixel (i1, j1) into the running extremum. Callers must keep one thread per
  cell ROW (the row index depends only on i1), which is what makes this lock free.
*/
void incUpdate(incBuffer *incBuf, int32_t i1, int32_t j1, double psi)
{
	extern int32_t rangeSelect;
	int32_t i, j;
	float cur;
	if (incBuf == NULL)
	{
		return;
	}
	i = cellOf(incBuf->i0, i1, incBuf->stride);
	j = cellOf(incBuf->j0, j1, incBuf->stride);
	if (i < 0 || i >= incBuf->nY || j < 0 || j >= incBuf->nX)
	{
		return;
	}
	cur = incBuf->inc[i][j];
	if (cur <= INCUNSET ||
		(rangeSelect == RANGESELECT_NEAR && psi < cur) ||
		(rangeSelect == RANGESELECT_FAR && psi > cur))
	{
		incBuf->inc[i][j] = (float)psi;
	}
}

float incLookup(incBuffer *incBuf, int32_t i1, int32_t j1)
{
	int32_t i, j;
	if (incBuf == NULL)
	{
		return INCUNSET;
	}
	i = cellOf(incBuf->i0, i1, incBuf->stride);
	j = cellOf(incBuf->j0, j1, incBuf->stride);
	if (i < 0 || i >= incBuf->nY || j < 0 || j >= incBuf->nX)
	{
		return INCUNSET;
	}
	return incBuf->inc[i][j];
}

/*
  Selection weight at (i1, j1), as a byte: 0 rejects the pixel, 255 keeps it at full weight.
  Without -angleRamp this is a hard cut (0 or 255). With it, the weight falls linearly from 255 at
  the extremum to 0 at the tolerance, so an input fades out instead of switching off - which
  matters where two tracks sit near the tolerance boundary and the switch would draw a hairline
  seam. A cell no input reached filters nothing.
*/
int32_t incWeight(incBuffer *incBuf, int32_t i1, int32_t j1, double psi)
{
	extern double angleTolerance;
	extern int32_t angleRamp;
	double d;
	float ext;
	ext = incLookup(incBuf, i1, j1);
	if (ext <= INCUNSET)
	{
		return 255;
	}
	d = fabs(psi - (double)ext);
	if (d > angleTolerance)
	{
		return 0;
	}
	if (angleRamp == FALSE)
	{
		return 255;
	}
	return (int32_t)(255.0 * (1.0 - d / angleTolerance) + 0.5);
}

/*
  Pass 1 for one range/Doppler input: geometry only, no raster read. Evaluates the ellipsoidal
  incidence angle at each coarse cell centre covering the image's region and folds it in.
*/
static void incidenceCoarseRD(inputImageStructure *inputImage, outputImageStructure *outputImage,
							  void *dem, incBuffer *incBuf, int32_t iMin, int32_t iMax,
							  int32_t jMin, int32_t jMax)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	extern double minIncidence;
	extern double maxIncidence;
	double savedLastTime;
	int32_t cMinX, cMaxX, cMinY, cMaxY, ci;
	/* cells covering the region; evaluated at their centres even if those fall outside the tile */
	incCellSpan(incBuf, TRUE, jMin, jMax - 1, &cMinX, &cMaxX);
	incCellSpan(incBuf, FALSE, iMin, iMax - 1, &cMinY, &cMaxY);
	/* llToImageNew warm start: restore it so pass 2 solves exactly as it would have alone */
	savedLastTime = inputImage->lastTime;
#pragma omp parallel for schedule(static) private(ci)
	for (ci = cMinY; ci <= cMaxY; ci++)
	{
		inputImageStructure myImg = *inputImage;
		int32_t cj, i1, j1;
		double x, y, lat, lon, h, hWGS, range, azimuth, aRange, ReH, psiSel;
		i1 = incCellCentre(incBuf, FALSE, ci);
		y = (outputImage->originY + i1 * outputImage->deltaY) * MTOKM;
		for (cj = cMinX; cj <= cMaxX; cj++)
		{
			j1 = incCellCentre(incBuf, TRUE, cj);
			x = (outputImage->originX + j1 * outputImage->deltaX) * MTOKM;
			xyToLLProj(x, y, &lat, &lon, &(outputImage->proj));
			h = getXYHeight(lat, lon, dem, myImg.cpAll.Re, SPHERICAL);
			hWGS = sphericalToWGSElev(h, lat, myImg.cpAll.Re);
			llToImageNew(lat, lon, hWGS, &range, &azimuth, &myImg);
			/*
			  Explicit bounds: llToImageNew only rejects past range +-2000 pixels and checkLL pads
			  the footprint by 50 km, so "not -9999" would let a frame claim tens of km of ground
			  it never imaged and win the extremum there.
			*/
			if (range < 0.0 || azimuth < 0.0 ||
				range > (double)(myImg.rangeSize - 1) || azimuth > (double)(myImg.azimuthSize - 1))
			{
				continue;
			}
			aRange = range * myImg.rangePixelSize + myImg.cpAll.RNear;
			ReH = getReH(&(myImg.cpAll), &myImg, azimuth);
			psiSel = psiRReZReH(aRange, myImg.cpAll.Re + h, ReH) * RTOD;
			/* a gated angle must not set the extremum: pass 2 throws that data away, so the
			   nearest ALLOWED look would otherwise be measured against a discarded one */
			if (psiSel < minIncidence || psiSel > maxIncidence)
			{
				continue;
			}
			incUpdate(incBuf, i1, j1, psiSel);
		}
	}
	inputImage->lastTime = savedLastTime;
}

void makeGeoMosaic(inputImageStructure *inputImage, outputImageStructure outputImage,
				   void *dem, int32_t nFiles, int32_t maxR, int32_t maxA, char **imageFiles, float fl, int32_t smoothL,
				   int32_t smoothOut, int32_t orbitPriority, float ***psiData, float ***gamma, gcovInputs *gcov)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	extern int32_t nearestDate;
	extern int32_t hybridZ;
	extern int32_t rsatFineCal;
	extern int32_t noPower;
	extern int32_t geoLinear;
	extern int32_t S1Cal;
	extern int32_t useSubPixelRTC;
	extern int32_t rangeSelect;
	extern int32_t angleStride;
	extern int32_t angleRamp;
	extern int32_t useIncidence;
	extern double minIncidence;
	extern double maxIncidence;
	double psiE;
	/* near/far range selection: coarse extremum buffer, per-pixel keep flag, and the unfiltered
	   fallback accumulator that makes the filter unable to leave a hole */
	incBuffer *incBuf = NULL;
	unsigned char **selTmp = NULL, *selBuf = NULL;
	float **imageAll = NULL, **scaleAll = NULL, **psiBufAll = NULL, **gBufAll = NULL;
	float **scale, **psiBuf, **psiBufTmp, **gBufTmp, **gBuf;
	FILE *fp;
	double x, y, lat, lon, h, hWGS;
	double range, azimuth;
	float value;
	float **imageTmp, **scaleTmp, **image;
	int32_t imageDate;
	int32_t i, j, i1, j1;
	int32_t iMin, iMax, jMin, jMax;
	float azimuthMin, azimuthMax;
	double Ab, Ag, AbCum, AgCum;
	double test;
	double dx[3], dy[3], invNAvg;
	int32_t recycle, shadow;
	clock_t startTime, lastTime, initTime;
	/* Get indices for oversampling output */
	oversampledxdy(smoothOut, dx, dy);
	/* Compute norm coefficient for average */
	invNAvg = 1.0 / (smoothOut * smoothOut);
	/*
	  Malloc space for tmp buffers
	*/
	mallocTmpBuffers(&imageTmp, &scaleTmp, &outputImage, &psiBuf, &psiBufTmp, &gBuf, &gBufTmp);
	*psiData = psiBuf;
	*gamma = gBuf;
	/* Init image, scale */
	initImageBuffers(&outputImage, scaleTmp, orbitPriority, psiBuf, psiBufTmp, gBuf, gBufTmp);
	image = (float **)outputImage.image;
	scale = (float **)outputImage.scale;
	fprintf(stderr, "\norbitPriority %i\n", orbitPriority);
	/*
	  Pass 1 of the range selection: build the coarse incidence extremum over every input before
	  any accumulation. Geometry only - no raster is read here.
	*/
	if (rangeSelect != RANGESELECT_NONE)
	{
		size_t nPix = (size_t)outputImage.xSize * outputImage.ySize;
		incBuf = allocIncBuffer(&outputImage, angleStride);
		imageAll = mallocPlane(&outputImage, 0.0);
		scaleAll = mallocPlane(&outputImage, 0.0);
		/* the unfiltered accumulation needs its own incidence and gamma-correction planes: it
		   runs after the previous file's filtered pass, so sharing psiBuf/gBuf would let an
		   unselected image supply the gamma correction for selected sigma0 */
		psiBufAll = mallocPlane(&outputImage, 0.0);
		gBufAll = mallocPlane(&outputImage, MINS1DB);
		selBuf = (unsigned char *)malloc(nPix);
		selTmp = (unsigned char **)malloc((size_t)outputImage.ySize * sizeof(unsigned char *));
		if (selBuf == NULL || selTmp == NULL)
		{
			error("makeGeoMosaic: malloc failed for the range-selection buffers");
		}
		for (i = 0; i < outputImage.ySize; i++)
		{
			selTmp[i] = &(selBuf[(size_t)i * outputImage.xSize]);
		}
		for (i = 0; i < nFiles; i++)
		{
			inputImage[i].file = imageFiles[i];
			if (!getRegion(&(inputImage[i]), &iMin, &iMax, &jMin, &jMax, &outputImage) ||
				inputImage[i].weight < 1e-20)
			{
				continue;
			}
			initllToImageNew(&(inputImage[i]));
			incidenceCoarseRD(&(inputImage[i]), &outputImage, dem, incBuf, iMin, iMax, jMin, jMax);
		}
		if (gcov != NULL)
		{
			gcovIncidenceCoarse(gcov, &outputImage, dem, incBuf);
		}
		fprintf(stderr, "rangeSelect pass 1 done\n");
	}
	/*
	  Main mosaicing loop
	*/
	fprintf(stderr, "STARTING MOSAICING\n");
	startTime = clock();
	for (i = 0; i < nFiles; i++)
	{
		/*
		  Process input image
		*/
		inputImage[i].file = imageFiles[i];
		/*  Get bounding box of image; skip if no overlap, file missing, or zero weight */
		if (!getRegion(&(inputImage[i]), &iMin, &iMax, &jMin, &jMax, &outputImage)
		    || inputImage[i].weight < 1e-20)
		{
			fprintf(stderr, "Skip (out of bounds, file missing, or 0 weight)\n");
			continue;
		}
		fprintf(stderr, "*** iMin, iMax, jMin, jMax ++ w %i %i %i %i ++ %f\n", iMin, iMax, jMin, jMax, inputImage[i].weight);
		/* Read image  */
		initllToImageNew(&(inputImage[i])); /* Setup conversions */
		getAzimuthBoundsForXYBox(iMin, iMax, jMin, jMax, &(inputImage[i]), &outputImage, &azimuthMin, &azimuthMax);
		fprintf(stderr, "%f %f\n", azimuthMin, azimuthMax);
		getGeoMosaicImage(&(inputImage[i]), &imageDate, smoothL, (int32_t)azimuthMin, (int32_t)azimuthMax);
		/*
		  Loop over output grid
		*/
		fprintf(stderr, "%s\n", imageFiles[i]);
		/* Each thread gets its own copy of inputImage[i] so that llToImageNew's
		   warm-start write to lastTime stays private and avoids cache-line
		   invalidation across all cores. */
		{
			int nthreads = omp_get_max_threads();
			inputImageStructure *localImgs = (inputImageStructure *)malloc(
				(size_t)nthreads * sizeof(inputImageStructure));
			if (localImgs == NULL)
				error("makeGeoMosaic: malloc failed for per-thread image copies\n");
			{
				int t;
				for (t = 0; t < nthreads; t++)
					localImgs[t] = inputImage[i];
			}
#pragma omp parallel \
			    private(j1, x, y, lat, lon, h, hWGS, range, azimuth, \
			            value, psiE, shadow, recycle, Ab, Ag, AbCum, AgCum, test)
			{
				inputImageStructure *myImg = &localImgs[omp_get_thread_num()];
#pragma omp for schedule(dynamic, 8)
				for (i1 = iMin; i1 < iMax; i1++)
				{
					if ((i1 % 100) == 0) fprintf(stderr, "-- %i\n", i1);
					for (j1 = jMin; j1 < jMax; j1++)
					{
						recycle = FALSE;
						shadow = FALSE;
						/*
						   This loop allows or overampling the result by computing multiple values about x and y.
						   In practice it doesn't help much, is slow, and should therefore be avoided.
						*/
						y = (outputImage.originY + i1 * outputImage.deltaY) * MTOKM;
						x = (outputImage.originX + j1 * outputImage.deltaX) * MTOKM;
						/*
						  Convert x/y stereographic coords to lat/lon
						*/
						extern int32_t linearSubPixelRTC;
						if (!linearSubPixelRTC)
						{
							xytoll1(x, y, HemiSphere, &lat, &lon, Rotation, outputImage.slat);
							/*
							  Get height for given lat/lon
							*/
							h = getXYHeight(lat, lon, dem, myImg->cpAll.Re, SPHERICAL);
							hWGS = sphericalToWGSElev(h, lat, myImg->cpAll.Re);
							/*
							  Convert lat/lon to image coordinates
							*/
							llToImageNew(lat, lon, hWGS, &range, &azimuth, myImg);
						}
						/*
						  Interpolate image
						*/
						if (useSubPixelRTC && (S1Cal & TRUE) == TRUE)
						{
							float  powerVal;
							double AbAcc, AgAcc;
							float  psiAcc;
							extern int32_t jacobianSubPixelRTC;
							extern int32_t linearSubPixelRTC;
							if ((jacobianSubPixelRTC
									? subPixelGammaRTCJacobian(x, y, myImg, &outputImage, dem,
															   &powerVal, &AbAcc, &AgAcc, &psiAcc)
								 : linearSubPixelRTC
									? subPixelGammaRTCLinear(x, y, myImg, &outputImage, dem,
															 &powerVal, &AbAcc, &AgAcc, &psiAcc)
									: subPixelGammaRTC(x, y, myImg, &outputImage, dem,
													   &powerVal, &AbAcc, &AgAcc, &psiAcc)))
								shadow = TRUE;
							value = powerVal;
							psiE  = psiAcc;
							AbCum = AbAcc;
							AgCum = AgAcc;
						}
						else
						{
							value = interpolatePowerInputImage(*myImg, range, azimuth);
							psiE = applyCorrections(&value, myImg, range, azimuth, h) * invNAvg;
							/*
							   RTC corrections for Sentinel.
							*/
							if ((S1Cal & TRUE) == TRUE)
							{
								if (areaAboutXY(range, azimuth, x, y, myImg, &outputImage, dem,
												value, &test, recycle, &Ab, &Ag) < -0.001)
									shadow = TRUE;
								AbCum = Ab;
								AgCum = Ag;
								recycle = TRUE;
							}
						}
						if (value > 0)
						{
							psiBufTmp[i1][j1] = psiE;
							/* Note the sin(psiE) undoes the psiE for sigma nought */
							if (shadow == FALSE && (S1Cal & TRUE) == TRUE)
							{
								double sinPsi = sin(psiE * DTOR);
								gBufTmp[i1][j1] = (sinPsi > 1e-6) ? 10.0 * log10((AbCum / AgCum) / sinPsi) : MINS1DB;
								gBufTmp[i1][j1] = round(gBufTmp[i1][j1] * 100.) / 100.;
								gBufTmp[i1][j1] = min(max(gBufTmp[i1][j1], -29.9), 35.0);
							}
							else
								gBufTmp[i1][j1] = -30.0; /* Negative value indicates shadow */
						}
						else
							value = -LARGEINT;
						/*
						  Range selection: the ellipsoidal incidence angle, computed here rather
						  than taken from psiE, which applyCorrections leaves at 0 on the antenna
						  pattern path. Rejection is deferred to applySelection so the unfiltered
						  fallback can be accumulated first.
						*/
						if (useIncidence == TRUE)
						{
							double aRangeSel = range * myImg->rangePixelSize + myImg->cpAll.RNear;
							double psiSel = psiRReZReH(aRangeSel, myImg->cpAll.Re + h,
													   getReH(&(myImg->cpAll), myImg, azimuth)) * RTOD;
							if (psiSel < minIncidence || psiSel > maxIncidence)
							{
								/* hard gate: drop it from the fallback accumulator too */
								value = -LARGEINT;
								if (rangeSelect != RANGESELECT_NONE)
								{
									selTmp[i1][j1] = 0;
								}
							}
							else if (rangeSelect != RANGESELECT_NONE)
							{
								selTmp[i1][j1] = (unsigned char)incWeight(incBuf, i1, j1, psiSel);
							}
						}
						/* Scaling */
						imageTmp[i1][j1] = value;
						if (value > 0 || (noPower > 0 && value > myImg->noData))
						{
							scaleTmp[i1][j1] = 1;
						}
					} /* End j1 */
				}	  /* End i1 */
			} /* End omp parallel */
			free(localImgs);
		}
		/*
		  With range selection, accumulate the UNFILTERED mosaic first - with its own feathering,
		  since the filtered and unfiltered data edges differ - then apply the filter in place and
		  fall through to the normal accumulation below.
		*/
		if (rangeSelect != RANGESELECT_NONE)
		{
			if (fl > 0 && (iMax > 0 && jMax > 0) && (iMax > iMin && jMax > jMin))
				computeScaleFast(imageTmp, scaleTmp, outputImage.ySize, outputImage.xSize, fl, inputImage[i].weight, (float)0.0,
								 iMin, iMax, jMin, jMax);
			geoMosaicScaling(&(inputImage[i]), imageAll, imageTmp, psiBufAll, psiBufTmp, gBufAll, gBufTmp, scaleAll,
							 scaleTmp, dem, &outputImage, orbitPriority, imageDate, iMin, iMax, jMin, jMax);
			applySelection(imageTmp, scaleTmp, selTmp, &(inputImage[i]), iMin, iMax, jMin, jMax);
		}
		/*
		  Compute scale array for feathering.
		*/
		if (fl > 0 && (iMax > 0 && jMax > 0) && (iMax > iMin && jMax > jMin))
			computeScaleFast(imageTmp, scaleTmp, outputImage.ySize, outputImage.xSize, fl, inputImage[i].weight, (float)0.0,
							 iMin, iMax, jMin, jMax);
		if (rangeSelect != RANGESELECT_NONE && angleRamp == TRUE)
		{
			applyRamp(scaleTmp, selTmp, iMin, iMax, jMin, jMax);
		}
		/*
		  Now sum current result. Falls through if no intersection (iMax&jMax==0)
		*/
		geoMosaicScaling(&(inputImage[i]), image, imageTmp, psiBuf, psiBufTmp, gBuf, gBufTmp, scale,
						 scaleTmp, dem, &outputImage, orbitPriority, imageDate, iMin, iMax, jMin, jMax);
	} /* End for i=0; i < nFiles */
	/*
	  Add already geocoded NISAR GCOV products, accumulated the same way as the range/Doppler products
	*/
	if (gcov != NULL)
	{
		for (i = 0; i < gcov->nFiles; i++)
		{
			inputImageStructure gcovImage;
			if (!gcovToOutputGrid(gcov, i, &outputImage, dem, imageTmp, scaleTmp, psiBufTmp, gBufTmp,
								  &gcovImage, &imageDate, &iMin, &iMax, &jMin, &jMax, incBuf, selTmp, 2))
			{
				continue;
			}
			if (rangeSelect != RANGESELECT_NONE)
			{
				if (fl > 0)
				{
					computeScaleFast(imageTmp, scaleTmp, outputImage.ySize, outputImage.xSize, fl, gcovImage.weight,
									 (geoLinear == TRUE) ? (float)GEOLINEARMIN : (float)0.0,
									 iMin, iMax, jMin, jMax);
				}
				geoMosaicScaling(&gcovImage, imageAll, imageTmp, psiBufAll, psiBufTmp, gBufAll, gBufTmp, scaleAll,
								 scaleTmp, dem, &outputImage, orbitPriority, imageDate, iMin, iMax, jMin, jMax);
				applySelection(imageTmp, scaleTmp, selTmp, &gcovImage, iMin, iMax, jMin, jMax);
			}
			if (fl > 0)
			{
				computeScaleFast(imageTmp, scaleTmp, outputImage.ySize, outputImage.xSize, fl, gcovImage.weight,
								 (geoLinear == TRUE) ? (float)GEOLINEARMIN : (float)0.0,
								 iMin, iMax, jMin, jMax);
			}
			if (rangeSelect != RANGESELECT_NONE && angleRamp == TRUE)
			{
				applyRamp(scaleTmp, selTmp, iMin, iMax, jMin, jMax);
			}
			geoMosaicScaling(&gcovImage, image, imageTmp, psiBuf, psiBufTmp, gBuf, gBufTmp, scale,
							 scaleTmp, dem, &outputImage, orbitPriority, imageDate, iMin, iMax, jMin, jMax);
		}
	}
	/*
	  Final rescaling if not nearestDate
	*/
	finalReScale(&outputImage, image, scale, orbitPriority, imageAll, scaleAll,
				 psiBuf, gBuf, psiBufAll, gBufAll);
	/* Convert values to log db if calibrated  added 11/18/2013 */
	if (rsatFineCal == TRUE)
	{
		rsatFinalCal(image, outputImage.xSize, outputImage.ySize);
	}
	if (rsatFineCal == TRUE || (S1Cal & TRUE) == TRUE)
	{
		logSigma(image, outputImage.xSize, outputImage.ySize);
	}
	lastTime = clock();
	fprintf(stderr, "totalTime %f\n", (double)(lastTime - startTime) / CLOCKS_PER_SEC);
	fprintf(stderr, "RETURN \n");
	return;
}

/*
  1-D squared distance transform (Felzenszwalb & Huttenlocher): d[q] = min_p (q - p)^2 + f[p].
  v (n ints) and z (n + 1 doubles) are workspace. All inputs are integers < 2^53, so the
  results are exact.
*/
static void edt1d(const double *f, int32_t n, double *d, int32_t *v, double *z)
{
	int32_t k = 0, q;
	double s;
	v[0] = 0;
	z[0] = -1e300;
	z[1] = 1e300;
	for (q = 1; q < n; q++)
	{
		s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k])) / (2.0 * q - 2.0 * v[k]);
		while (s <= z[k])
		{
			k--;
			s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k])) / (2.0 * q - 2.0 * v[k]);
		}
		k++;
		v[k] = q;
		z[k] = s;
		z[k + 1] = 1e300;
	}
	k = 0;
	for (q = 0; q < n; q++)
	{
		while (z[k + 1] < q)
		{
			k++;
		}
		d[q] = (double)(q - v[k]) * (q - v[k]) + f[v[k]];
	}
}

/*
  Same feathering weights as common/computeScale.c, in O(pixels) instead of O(edge pixels * fl^2).

  computeScale sets every pixel to weight, then for each data-edge pixel (valid, with an invalid
  pixel in its 3x3 neighbourhood) lowers the pixels within +-fl to the radial kernel
  w1 * min(max(dist, 0.5) / fl, 1), where w1 is the weight of the FIRST call (the kernel is
  cached in rDistSave). With weight == w1 that is weight * min(max(d, 0.5) / fl, 1) for d the
  distance to the nearest edge pixel, which an exact squared-distance transform gives directly;
  the kernel expression is evaluated identically, so results are bit-identical. Only the output
  region [iMin,iMax) x [jMin,jMax) is used by geoMosaicScaling, so the transform covers that
  region plus fl (edges further away cannot reach it); elsewhere scale is weight, as initialised.
  If weight != w1 the original computeScale is called, keeping that case identical by construction.
*/
static void computeScaleFast(float **inImage, float **scale, int32_t ySize, int32_t xSize, float fl, float weight,
							 double minVal, int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax)
{
	extern double **rDistSave;
	static float firstWeight = -1.0;
	int32_t i0, i1, j0, j1, ny, nx, n, i, j, is, il;
	int32_t flInt = (int32_t)fl;
	double cap, dMax, rA;
	float minV, v;
	if (fl == 0)
	{
		computeScale(inImage, scale, ySize, xSize, fl, weight, minVal);
		return;
	}
	if (firstWeight < 0)
	{
		firstWeight = weight;
	}
	if (weight != firstWeight)
	{
		/* computeScale's cached kernel must hold the first weight, exactly as if it had made it */
		if (rDistSave == NULL)
		{
			rDistSave = dmatrix(-fl, fl, -fl, fl);
			fillRadialKernel(rDistSave, fl, firstWeight);
		}
		computeScale(inImage, scale, ySize, xSize, fl, weight, minVal);
		return;
	}
	initFloatMatrix(scale, ySize, xSize, weight);
	if (iMax <= iMin || jMax <= jMin)
	{
		return;
	}
	i0 = max(0, iMin - flInt);
	i1 = min(ySize, iMax + flInt);
	j0 = max(0, jMin - flInt);
	j1 = min(xSize, jMax + flInt);
	ny = i1 - i0;
	nx = j1 - j0;
	n = max(nx, ny);
	/* anything beyond fl is weight; cap keeps the arithmetic small and exact */
	dMax = (double)flInt * flInt;
	cap = 4.0 * (flInt + 2.0) * (flInt + 2.0);
	double *g = (double *)malloc(sizeof(double) * (size_t)nx * ny);
	if (g == NULL)
	{
		error("computeScaleFast: malloc failed");
	}
	/* edge pixels (0) vs others (cap), exactly computeScale's test */
	for (i = i0; i < i1; i++)
	{
		for (j = j0; j < j1; j++)
		{
			g[(size_t)(i - i0) * nx + (j - j0)] = cap;
			if (inImage[i][j] > minVal)
			{
				minV = 1.0e30;
				for (is = max(0, i - 1); is <= min(ySize - 1, i + 1); is++)
				{
					for (il = max(0, j - 1); il <= min(xSize - 1, j + 1); il++)
					{
						minV = min(minV, inImage[is][il]);
					}
				}
				if (minV <= minVal)
				{
					g[(size_t)(i - i0) * nx + (j - j0)] = 0.0;
				}
			}
		}
	}
	/* columns, then rows */
#pragma omp parallel private(i, j)
	{
		double *f = (double *)malloc(sizeof(double) * n);
		double *d = (double *)malloc(sizeof(double) * n);
		double *z = (double *)malloc(sizeof(double) * (n + 1));
		int32_t *vv = (int32_t *)malloc(sizeof(int32_t) * n);
#pragma omp for schedule(static)
		for (j = 0; j < nx; j++)
		{
			for (i = 0; i < ny; i++)
			{
				f[i] = g[(size_t)i * nx + j];
			}
			edt1d(f, ny, d, vv, z);
			for (i = 0; i < ny; i++)
			{
				g[(size_t)i * nx + j] = min(d[i], cap);
			}
		}
#pragma omp for schedule(static)
		for (i = 0; i < ny; i++)
		{
			for (j = 0; j < nx; j++)
			{
				f[j] = g[(size_t)i * nx + j];
			}
			edt1d(f, nx, d, vv, z);
			for (j = 0; j < nx; j++)
			{
				g[(size_t)i * nx + j] = d[j];
			}
		}
		free(f);
		free(d);
		free(z);
		free(vv);
	}
	/* kernel value, evaluated as fillRadialKernel does, then FMIN with the initial weight */
	for (i = i0; i < i1; i++)
	{
		for (j = j0; j < j1; j++)
		{
			double d2 = g[(size_t)(i - i0) * nx + (j - j0)];
			if (d2 <= dMax)
			{
				rA = firstWeight * min(max(sqrt(d2), 0.5) / fl, 1);
				v = (float)rA;
				scale[i][j] = (v < scale[i][j]) ? v : scale[i][j];
			}
		}
	}
	free(g);
}

static void finalReScale(outputImageStructure *outputImage, float **image, float **scale, int32_t orbitPriority,
						 float **imageAll, float **scaleAll, float **psiBuf, float **gBuf,
						 float **psiBufAll, float **gBufAll)
{
	int32_t i1, j1;
	size_t nRescued = 0, nValid = 0;
	if (orbitPriority >= 0)
		return;

	for (i1 = 0; i1 < outputImage->ySize; i1++)
	{
		for (j1 = 0; j1 < outputImage->xSize; j1++)
		{
			if (scale[i1][j1] > 0)
			{
				/* this assumes that dates are larger than 1000, and we won't sum more than 1000*/
				if (scale[i1][j1] < 10000)
					image[i1][j1] /= scale[i1][j1];
				nValid++;
			}
			else if (imageAll != NULL && scaleAll[i1][j1] > 0)
			{
				/*
				  Range selection emptied this pixel: pass 1 sees only geometry, so it can pick an
				  input that turns out to have no data here. Fall back to the unfiltered average
				  so the filter can never remove coverage the plain mosaic would have had.
				*/
				image[i1][j1] = (scaleAll[i1][j1] < 10000) ? imageAll[i1][j1] / scaleAll[i1][j1]
														   : imageAll[i1][j1];
				/* carry the matching incidence and gamma correction across too */
				psiBuf[i1][j1] = psiBufAll[i1][j1];
				gBuf[i1][j1] = gBufAll[i1][j1];
				nRescued++;
				nValid++;
			}
		}
	}
	if (imageAll != NULL)
	{
		fprintf(stderr, "rangeSelect: %lu of %lu valid pixels (%.3f%%) fell back to the "
				"unfiltered average\n", nRescued, nValid,
				(nValid > 0) ? 100.0 * (double)nRescued / (double)nValid : 0.0);
	}
}

/* One full-grid float plane, contiguous with a row-pointer array, filled with initValue. */
static float **mallocPlane(outputImageStructure *outputImage, float initValue)
{
	float **plane, *buf;
	size_t n = (size_t)outputImage->xSize * outputImage->ySize, k;
	int32_t i;
	buf = (float *)malloc(n * sizeof(float));
	plane = (float **)malloc((size_t)outputImage->ySize * sizeof(float *));
	if (buf == NULL || plane == NULL)
	{
		error("mallocPlane: malloc failed for %i x %i", outputImage->xSize, outputImage->ySize);
	}
	for (k = 0; k < n; k++)
	{
		buf[k] = initValue;
	}
	for (i = 0; i < outputImage->ySize; i++)
	{
		plane[i] = &(buf[(size_t)i * outputImage->xSize]);
	}
	return plane;
}

/*
  Scale the feather weight by the -angleRamp selection weight. Must run AFTER computeScaleFast,
  which rebuilds scaleTmp from imageTmp and would otherwise discard it.
*/
static void applyRamp(float **scaleTmp, unsigned char **selTmp,
					  int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax)
{
	int32_t i1, j1;
	for (i1 = iMin; i1 < iMax; i1++)
	{
		for (j1 = jMin; j1 < jMax; j1++)
		{
			if (selTmp[i1][j1] > 0)
			{
				scaleTmp[i1][j1] *= (float)selTmp[i1][j1] / 255.0f;
			}
		}
	}
}

/*
  Drop the pixels the range selection rejected. Called after the unfiltered accumulation, so the
  fallback keeps its own values. scaleTmp is restored for kept pixels because geoMosaicScaling
  zeroes it as it goes and, with fl == 0, nothing else would set it again.
*/
static void applySelection(float **imageTmp, float **scaleTmp, unsigned char **selTmp,
						   inputImageStructure *inputImage, int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax)
{
	extern int32_t noPower;
	extern int32_t geoLinear;
	int32_t i1, j1;
	for (i1 = iMin; i1 < iMax; i1++)
	{
		for (j1 = jMin; j1 < jMax; j1++)
		{
			if (selTmp[i1][j1] == 0)
			{
				imageTmp[i1][j1] = -LARGEINT;
			}
			else if (imageTmp[i1][j1] > 0 || (noPower > 0 && imageTmp[i1][j1] > inputImage->noData) ||
					 (geoLinear == TRUE && imageTmp[i1][j1] > GEOLINEARMIN))
			{
				scaleTmp[i1][j1] = 1;
			}
		}
	}
}

static void geoMosaicScaling(inputImageStructure *inputImage, float **image, float **imageTmp, float **psiBuf, float **psiBufTmp,
							 float **gBuf, float **gBufTmp, float **scale, float **scaleTmp, void *dem, outputImageStructure *outputImage, int32_t orbitPriority,
							 int32_t imageDate, int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	extern int32_t hybridZ;
	extern int32_t nearestDate;
	extern int32_t noPower;
	extern int32_t geoLinear;
	extern int32_t geoMosaicMode;
	double x, y, hWGS;
	double lat, lon;
	int32_t i1, j1;
	for (i1 = iMin; i1 < iMax; i1++)
	{
		y = (outputImage->originY + i1 * outputImage->deltaY) * MTOKM;
		for (j1 = jMin; j1 < jMax; j1++)
		{
			x = (outputImage->originX + j1 * outputImage->deltaX) * MTOKM;
			/* Get elevation */
			if (hybridZ > 0)
			{
				xyToLLProj(x, y, &lat, &lon, &(outputImage->proj));
				hWGS = getXYHeight(lat, lon, dem, inputImage->cpAll.Re, ELLIPSOIDAL);
			}
			else
				hWGS = hybridZ - 1; /* This will force skip */

			if (imageTmp[i1][j1] > 0 || (noPower > 0 && imageTmp[i1][j1] > inputImage->noData) ||
				(geoLinear == TRUE && imageTmp[i1][j1] > GEOLINEARMIN))
			{ /* Points with valid data */
				if (geoMosaicMode == GEOMOSAIC_MIN || geoMosaicMode == GEOMOSAIC_MAX)
				{ /* pixel-wise min or max across all inputs */
					float candidate = imageTmp[i1][j1] * scaleTmp[i1][j1];
					if (scale[i1][j1] <= 0)
					{ /* first valid pixel at this location */
						image[i1][j1] = candidate;
						scale[i1][j1] = 1;
						psiBuf[i1][j1] = psiBufTmp[i1][j1];
						gBuf[i1][j1] = gBufTmp[i1][j1];
					}
					else if ((geoMosaicMode == GEOMOSAIC_MIN && candidate < image[i1][j1]) ||
					         (geoMosaicMode == GEOMOSAIC_MAX && candidate > image[i1][j1]))
					{
						image[i1][j1] = candidate;
						psiBuf[i1][j1] = psiBufTmp[i1][j1];
						gBuf[i1][j1] = gBufTmp[i1][j1];
					}
				}
				/* case for no nearestDate, no orbitPriority, or hWGS override */
				else if ((nearestDate < 0 || (nearestDate > 0 && (int)hWGS > hybridZ)) && orbitPriority < 0)
				{ /* Summing data or non-nearest date or hybridZ*/
					image[i1][j1] += imageTmp[i1][j1] * scaleTmp[i1][j1] * inputImage->weight;
					scale[i1][j1] += scaleTmp[i1][j1];
					psiBuf[i1][j1] = psiBufTmp[i1][j1];
					gBuf[i1][j1] = gBufTmp[i1][j1];
				}
				else if (orbitPriority > -1)
				{ /* orbit priority case */
					if (orbitPriority == ASCENDING)
					{
						/* either put in data if none already, or if not ascending replace */
						if (scale[i1][j1] < -0.1 || scale[i1][j1] == DESCENDING)
						{
							image[i1][j1] = imageTmp[i1][j1] * inputImage->weight;
							scale[i1][j1] = inputImage->passType;
							psiBuf[i1][j1] = psiBufTmp[i1][j1];
							gBuf[i1][j1] = gBufTmp[i1][j1];
						}
					}
					else if (orbitPriority == DESCENDING)
					{
						/* either put in data if none already, or if not descending replace */
						if (scale[i1][j1] < -0.1 || scale[i1][j1] == ASCENDING)
						{
							image[i1][j1] = imageTmp[i1][j1] * inputImage->weight;
							scale[i1][j1] = inputImage->passType;
							psiBuf[i1][j1] = psiBufTmp[i1][j1];
							gBuf[i1][j1] = gBufTmp[i1][j1];
						}
					}
					else
						error("bad orbit priority flag");
				}
				else
				{ /* Nearest date, if orbitPriority not st */
					if ((abs(nearestDate - imageDate) < fabs(nearestDate - scale[i1][j1])) && scaleTmp[i1][j1] > 0)
					{
						image[i1][j1] = imageTmp[i1][j1] * inputImage->weight;
						scale[i1][j1] = imageDate;
						psiBuf[i1][j1] = psiBufTmp[i1][j1];
						gBuf[i1][j1] = gBufTmp[i1][j1];
					}
				} /* end else */
			}	  /* end if imageTmp...*/
			scaleTmp[i1][j1] = 0.0;
		} /* end j1 */
	}	  /* end i1=iMin sum current ... */
}

float applyCorrections(float *value, inputImageStructure *inputImage, double range, double azimuth, double h)
{
	extern int32_t rsatFineCal;
	extern int32_t S1Cal;
	extern int32_t noPower;
	int32_t index;
	float theta, psi;
	double ReH, aRange, drA;
	psi = 0.0;
	if (*value > 0)
	{
		if (inputImage->rAnt != NULL && noPower > 0)
		{
			/* use antenna pattern */
			aRange = range * inputImage->rangePixelSize + inputImage->cpAll.RNear;
			drA = (inputImage->rAnt[1] - inputImage->rAnt[0]);
			index = (aRange - inputImage->rAnt[0]) / drA;
			index = min(max(0, index), inputImage->patSize);
			*value = *value * inputImage->pAnt[index];
		}
		else
		{
			aRange = range * inputImage->rangePixelSize + inputImage->cpAll.RNear;
			ReH = getReH(&(inputImage->cpAll), inputImage, azimuth);
			theta = thetaRReZReH(aRange, (inputImage->cpAll.Re + h), ReH) * RTOD;
			psi = psiRReZReH(aRange, (inputImage->cpAll.Re + h), ReH);
			if ((inputImage->patSize < 0 && noPower < 0))
			{
				/*   use polynomial function - radarsat or alos specific
					 At present this only does fine beam RS, but function could easily be modified for other antenna ppatters
					 calibrated  added 11/18/2013
				*/
				if (rsatFineCal == TRUE)
				{
					*value *= 1.0 / polypat(theta);
					/* ADDED 3/4/15 to reference sigma nought  to ellipsoid
					   removed because it made more stripy  *value*= sin(psi)/sin(inputImage->geoData.psic * DTOR);
					   moved to final part           	   *value *=1.0/27.3;
					   *value-=0.0058;*/
					/* Set noise floor at -30db*/
					/*	      if(*value < 0.001) *value= 0.001;*/
				}
				else
				{
					if (inputImage->patSize == RSATFINE)
						*value *= 1.0 / polypat(theta);
					if (inputImage->patSize == ALOS)
					{
						*value *= 1.0 / polyALOS(theta);
					}
					/* ADDED 3/4/15 to reference sigma nought  to ellipsoid
					 *value*= sin(psi)/sin(inputImage->geoData.psic * DTOR);
					 Removed because it caused additional striping.
					*/
				}
			}
			else if ((S1Cal & TRUE) == TRUE)
			{
				if (inputImage->betaNought > 0)
				{
					*value *= sin(psi) / pow(inputImage->betaNought, 2.0);
				}
				else
				{
					fprintf(stderr, "input image %s", inputImage->file);
					error("S1 Cal, but beta Nought < 0 ");
				}
			}
			else
			{
				/*
				  This equation is for uncalibrated results - primarily sentinel. The first sin(psi) scales
				  for sigma nought. The sin(psi)**2 appears to cosmetically correct for the incidence angle
				*/
				if (noPower < 0)
					*value *= sin(psi) * sin(psi) * sin(psi);
				/*	attempt at lambertian backscatter *value *= cos(37.*DTOR)*cos(37.*DTOR)*sin(psi)/(cos(psi)*cos(psi)); */
			}
		}
	}
	return (float)(psi * RTOD);
}

static void getGeoMosaicImage(inputImageStructure *inputImage, int32_t *imageDate, int32_t smoothL, int32_t yMin, int32_t yMax)
{
	extern int32_t nearestDate;
	float **tmpImage;
	int32_t extraPad, tmpi;
	int32_t doy[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 333};
	int32_t k1, k2, i;
	
	*imageDate = inputImage->year * 365. + doy[inputImage->month - 1] + inputImage->day;
	if (nearestDate > 0)
		fprintf(stderr, "Image date %i, Nearest Date %i %i %i %i\n",
				*imageDate, nearestDate, inputImage->year, inputImage->month, inputImage->day);
	readComplexAsPower(inputImage, yMin, yMax);
	/* Multilook the image */
	if (smoothL > 0)
		smoothImage(inputImage, smoothL);
	if (inputImage->removePad > 0)
	{
		extraPad = 0;
		tmpi = inputImage->rangeSize * inputImage->nRangeLooks;
		/* 1/23/06 This removes a little extra on wider images (should mostly apply to fine
		   beam. This should handle any range shifts */
		if (tmpi > 8500)
			extraPad = 120 / inputImage->nRangeLooks;
		/* Pad around edges */
		for (k2 = 0; k2 < inputImage->azimuthSize; k2++)
			for (k1 = 0; k1 < (inputImage->removePad + extraPad); k1++)
			{
				tmpImage = (float **)inputImage->image;
				tmpImage[k2][k1] = 0.;
				tmpImage[k2][inputImage->rangeSize - k1 - 1] = 0.;
			}
	}
}

static void initImageBuffers(outputImageStructure *outputImage, float **scaleTmp, int32_t orbitPriority, float **psiBuf, float **psiBufTmp, float **gBuf, float **gBufTmp)
{
	float **image, **scale;
	int32_t i1, j1;
	image = (float **)outputImage->image;
	scale = (float **)outputImage->scale;

	for (i1 = 0; i1 < outputImage->ySize; i1++)
	{
		for (j1 = 0; j1 < outputImage->xSize; j1++)
		{
			scale[i1][j1] = 0.0;
			image[i1][j1] = 0.0;
			scaleTmp[i1][j1] = 0.0;
			psiBuf[i1][j1] = 0.0;
			psiBufTmp[i1][j1] = 0.0;
			gBuf[i1][j1] = MINS1DB;
			gBufTmp[i1][j1] = MINS1DB;
			if (orbitPriority > -1)
				scale[i1][j1] = -1;
		}
	}
}

static void mallocTmpBuffers(float ***imageTmp, float ***scaleTmp, outputImageStructure *outputImage, float ***psiBuf, float ***psiBufTmp, float ***gBuf, float ***gBufTmp)
{
	size_t bufSize;
	float *bufI, *bufS, *bufP, *bufPTmp, *bufG, *bufGTmp;
	int32_t i;
	int32_t myBuff;
	myBuff = 0;
	bufSize = outputImage->xSize * outputImage->ySize * sizeof(float);
	bufI = (float *)malloc((size_t)bufSize);
	myBuff += bufSize;
	if (bufI == NULL)
		error("Unable to malloc buffer 1 in mallocTmpBuffers already allocated %f new buf %f\n\n", myBuff / 1e6, bufSize / 1e6);
	bufS = (float *)malloc((size_t)bufSize);
	myBuff += bufSize;
	if (bufS == NULL)
		error("Unable to malloc buffer 2 in mallocTmpBuffers already allocated %f new buf %f\n\n", myBuff / 1e6, bufSize / 1e6);
	bufP = (float *)malloc((size_t)bufSize);
	myBuff += bufSize;
	if (bufP == NULL)
		error("Unable to malloc buffer 3 in mallocTmpBuffers already allocated %f new buf %f\n\n", myBuff / 1e6, bufSize / 1e6);
	bufPTmp = (float *)malloc((size_t)bufSize);
	myBuff += bufSize;
	if (bufPTmp == NULL)
		error("Unable to malloc buffer 4 in mallocTmpBuffers  already allocated %f new buf %f\n\n", myBuff / 1e6, bufSize / 1e6);
	bufG = (float *)malloc((size_t)bufSize);
	myBuff += bufSize;
	if (bufG == NULL)
		error("Unable to malloc buffer 5 in mallocTmpBuffers already allocated %f new buf %f\n\n", myBuff / 1e6, bufSize / 1e6);
	bufGTmp = (float *)malloc((size_t)bufSize);
	myBuff += bufSize;
	if (bufGTmp == NULL)
		error("Unable to malloc buffer 6 in mallocTmpBuffers already allocated %f new buf %f\n\n", myBuff / 1e6, bufSize / 1e6);
	bufSize = sizeof(float *) * outputImage->ySize;
	*imageTmp = (float **)malloc((size_t)bufSize);
	*scaleTmp = (float **)malloc((size_t)bufSize);
	*psiBuf = (float **)malloc((size_t)bufSize);
	*psiBufTmp = (float **)malloc((size_t)bufSize);
	*gBuf = (float **)malloc((size_t)bufSize);
	*gBufTmp = (float **)malloc((size_t)bufSize);

	for (i = 0; i < outputImage->ySize; i++)
	{
		(*imageTmp)[i] = (float *)&(bufI[i * outputImage->xSize]);
		(*scaleTmp)[i] = (float *)&(bufS[i * outputImage->xSize]);
		(*psiBuf)[i] = (float *)&(bufP[i * outputImage->xSize]);
		(*psiBufTmp)[i] = (float *)&(bufPTmp[i * outputImage->xSize]);
		(*gBuf)[i] = (float *)&(bufG[i * outputImage->xSize]);
		(*gBufTmp)[i] = (float *)&(bufGTmp[i * outputImage->xSize]);
	}
}

/* final cal for rsat */
static void rsatFinalCal(float **image, int32_t xSize, int32_t ySize)
{
	int32_t i1, j1;
	for (i1 = 0; i1 < ySize; i1++)
	{
		for (j1 = 0; j1 < xSize; j1++)
		{
			image[i1][j1] *= 1.0 / 27.3;
			image[i1][j1] -= 0.0058;
		}
	}
}

/* compute 10log(sig) for calibrated results. Clip data to -29.9 */
static void logSigma(float **image, int32_t xSize, int32_t ySize)
{
	int32_t i1, j1;
	for (i1 = 0; i1 < ySize; i1++)
	{
		for (j1 = 0; j1 < xSize; j1++)
		{
			/* Limit small values to 0.001001 = -29.92 db , for no data for -30.0 */
			if (image[i1][j1] < 0.001001)
			{
				if (image[i1][j1] > 1.e-9)
					image[i1][j1] = 0.00102;
				else
					image[i1][j1] = MINS1SIG;
			}
			/* Round to 1/100 of dB - in conversion to tif, it may get rounded further */
			image[i1][j1] = roundf((10.0 * log10(image[i1][j1])) * 100) / 100.;
		}
	}
}

/*
  On 3/21/2013 this just does F1 RADARSAT with a 9th degree polynomial that provides
  a tight match to the antenna pattern provided by ASF.
  The pattern is in dB, but the gain is return as a look-angle dependent scale factor.
*/
static float polypat(float theta)
{
	double thetaCentered;
	double logGain;
	float gain;
	thetaCentered = (theta - 33.399360) / 1.880050;
	logGain = 0.0;
	logGain += 0.044834 * pow(thetaCentered, 9);
	logGain += -0.048875 * pow(thetaCentered, 8);
	logGain += -0.338039 * pow(thetaCentered, 7);
	logGain += 0.414746 * pow(thetaCentered, 6);
	logGain += 0.767037 * pow(thetaCentered, 5);
	logGain += -2.244836 * pow(thetaCentered, 4);
	logGain += -1.038396 * pow(thetaCentered, 3);
	logGain += 1.074398 * pow(thetaCentered, 2);
	logGain += 0.935996 * pow(thetaCentered, 1);
	logGain += 1.707515;
	gain = (float)pow(10.0, logGain / 10.0);
	return gain;
}

static float polyALOS(float theta)
{
	double thetaCentered;
	float gain;
	thetaCentered = (theta - 34.0);
	gain = 0.0;
	gain += 0.00006896 * pow(thetaCentered, 4);
	gain += 0.00011602 * pow(thetaCentered, 3);
	gain += -0.05117270 * pow(thetaCentered, 2);
	gain += 0.00516289 * pow(thetaCentered, 1);
	gain += 1.00243524;
	gain = gain * gain; /* tests indicate this is 2way */
	return gain;
}

static float interpolatePowerInputImage(inputImageStructure inputImage, double range, double azimuth)
{
	float **fimage;
	int32_t i, j;
	/*
	  fimage=(float **)inputImage.image;
	  j = (int)round(range);
	  i = (int)round(azimuth);
	  if(range < 0.0 || azimuth < 0.0  ||  j >= inputImage.rangeSize || i >= inputImage.azimuthSize ) return inputImage.noData;
	  return fimage[i][j];
	*/
	fimage = (float **)inputImage.image;
	return bilinearInterp(fimage, range, azimuth, inputImage.rangeSize, inputImage.azimuthSize, inputImage.noData, inputImage.noData);
}

static void smoothImage(inputImageStructure *inputImage, int32_t smoothL)
{
	extern float *smoothBuf;
	float **image;
	float *filt;
	float sum = 0.0, sumF;
	int32_t i, j, k;
	int32_t hw;
	if (smoothL < 1)
		return;
	if (smoothL == 1)
		hw = 1;
	else
		hw = (int)(smoothL / 2);

	float *filtBase = (float *)malloc((size_t)((hw * 2 + 1) * sizeof(float)));
	filt = &(filtBase[hw]);
	if (smoothL > 1)
		for (i = -hw; i <= hw; i++)
			filt[i] = 1;
	else
	{
		filt[-1] = 0.5;
		filt[0] = 1.0;
		filt[1] = 0.5;
	}

	sumF = 0;
	for (i = -hw; i <= hw; i++)
	{
		fprintf(stderr, "filt %f\n", filt[i]);
		sumF += filt[i];
	}
	fprintf(stderr, "hw = %i %f\n", hw, sum);
	if (smoothBuf == NULL)
		error("memory not allocated for smoothImage\n");
	image = (float **)inputImage->image;
	/*
	  smooth in range direction - note borders can be irregular so be careful not to smooth across data/zero transition
	*/
	for (i = 0; i < inputImage->azimuthSize; i++)
	{
		/* Smooth line and save in buf */
		for (j = hw; j < (inputImage->rangeSize - hw); j++)
		{
			smoothBuf[j] = 0;
			sum = 0;
			for (k = -hw; k <= hw; k++)
			{
				if (image[i][j + k] > 1.e-6)
				{
					smoothBuf[j] += image[i][j + k] * filt[k];
					sum += filt[k];
				}
			}
			if (sum < sumF)
				smoothBuf[j] = 0.0;
			else
				smoothBuf[j] /= sum;
		}
		/* Write result back to main buf */
		for (j = hw; j < inputImage->rangeSize - hw; j++)
			image[i][j] = smoothBuf[j];
	}
	fprintf(stderr, "-- %i %i \n", i, j);
	/*
	  smooth in azimuth direction
	*/
	for (j = 0; j < inputImage->rangeSize; j++)
	{
		/* Smooth line and save in buf */
		for (i = hw; i < (inputImage->azimuthSize - hw); i++)
		{
			smoothBuf[i] = 0;
			sum = 0.0;
			for (k = -hw; k <= hw; k++)
			{
				if (image[i + k][j] > 1e-6)
				{
					sum += filt[k];
					smoothBuf[i] += image[i + k][j] * filt[k];
				}
			}
			if (sum < sumF)
				smoothBuf[i] = 0.0;
			else
				smoothBuf[i] /= sum;
		}
		/* Write result back to main buf */
		for (i = hw; i < inputImage->azimuthSize - hw; i++)
			image[i][j] = smoothBuf[i];
	}
	fprintf(stderr, "++ %i %i \n", i, j);
	free(filtBase);
}
