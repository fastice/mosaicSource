#include "stdio.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
/* #include "geotiff/xtiffio.h"   for TIFF */
/* #include "geotiff/geotiffio.h" for GeoTIFF */
#include "cRecipes/nrutil.h"
#include "mosaicSource/common/common.h"
/*#include "mosaicSource/mosaic3D_p/mosaic3d.h"*/
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
#include "mosaicSource/landsatMosaic/landSatMosaic.h"

void computeScaleLS(float **inImage, float **scale, int32_t azimuthSize, int32_t rangeSize, float fl, float weight, double minVal, int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax);
static void readLSOffsetsForMosaic(landSatImage *currentImage);
static void computeFitError(double x, double y, landSatImage *currentImage, double *sigx2, double *sigy2);
static float **LSreadFloatImageforMosaic(char *file, int32_t nx, int32_t ny, float *fBuffer, float **image);
static uint8_t **LSreadByteImageforMosaic(char *file, int32_t nx, int32_t ny, uint8_t *maskB, uint8_t **m);
static void getLSRegion(landSatImage *image, int32_t *iMin, int32_t *iMax, int32_t *jMin, int32_t *jMax, outputImageStructure *outputImage);
static int32_t interpLSdata(double x, double y, matchResult *matches, double *dx, double *dy, double *sx, double *sy);

#define MAXSIG 0.2 /* Maximum ls error */
/*
************************ Mosaic Landsat data. **************************
*/
#define MAXOFFBUFLS 30000000
#define MAXOFFBUFLINESLS 30000
#define MAXMASKBUFLINES 10000
#define MAXLINEBUF 2000
float *fBuf1, *fBuf2, *fBuf3, *fBuf4;
float **fBuf1L, **fBuf2L, **fBuf3L, **fBuf4L;
uint8_t **mask, *maskBuf;
char *lineBuf;

void makeLandSatMosaic(landSatImage *LSImages, outputImageStructure *outputImage, float fl)
{
	extern int HemiSphere;
	extern double Rotation;
	extern float *fBuf1, *fBuf2, *fBuf3, *fBuf4;
	extern float **fBuf1L, **fBuf2L, **fBuf3L, **fBuf4L;
	extern char *lineBuf;
	extern void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
	extern void *lBuf1, *lBuf2, *lBuf3, *lBuf4;
	double sigx2, sigy2; /* sigma^2 from plane fit */
	double x0poly, y0poly;
	float **vXimage, **vYimage, **vZimage;
	float **scaleX, **scaleY, **scaleZ;
	float **vxTmp, **vyTmp, **vzTmp, **fScale, **sxTmp, **syTmp;
	float **errorX, **errorY;
	double x, y;
	double sigmaQx, sigmaQy;
	double scX, scY;
	double vx, vy, ex, ey; /* velocity and error */
	double tmp;
	double lat, lon;
	double pX[MAXFITPARAM], pY[MAXFITPARAM];
	unsigned char sMask;
	ShelfMask *shelfMask; /* Mask with shelf and grounding zone */
	double offx, offy, sx, sy, sx2, sy2;
	double vscalex, vscaley, latscale;
	double sigmaXTieRes, sigmaXImageAvg, sigmaYTieRes, sigmaYImageAvg;
	double tCenter, tOffCenter, deltaOffCenter, tHalfWidth;
	double deltaT;
	double varXGlobal, varYGlobal;
	landSatImage *currentImage;
	int32_t i, j, k, iMin, iMax, jMin, jMax; /* LCV */
	int32_t count;
	int keep;
	size_t lsize;
	/*
	  Malloc tmp space used later
	*/
	lineBuf = (char *)malloc(sizeof(char) * MAXLINEBUF);
	fprintf(stderr, "**** LANDSAT TRACKING SOLUTION ****\n");
	/* This taps into the common memory pool allocated up front in mosaic 3d */
	fBuf1 = (float *)offBufSpace1;
	fBuf2 = (float *)offBufSpace2;
	fBuf3 = (float *)offBufSpace3;
	fBuf4 = (float *)offBufSpace4;
	lsize = MAXOFFBUFLINESLS * sizeof(float *);
	fBuf1L = (float **)lBuf1;
	fBuf2L = (float **)lBuf2;
	fBuf3L = (float **)lBuf3;
	fBuf4L = (float **)lBuf4;
	/* Added for masks Summer 2017 */
	mask = (uint8_t **)malloc(MAXMASKBUFLINES * sizeof(uint8_t *));
	maskBuf = (uint8_t *)malloc(MAXMASKBUFLINES * MAXMASKBUFLINES * sizeof(uint8_t));
	/*
	  Added shelf mask August 2017 - not to be confused with individual masks above
	*/
	shelfMask = outputImage->shelfMask;
	/*
	  Pointers to output images
	*/
	vXimage = (float **)outputImage->image;
	vYimage = (float **)outputImage->image2;
	vZimage = (float **)outputImage->image3;
	scaleX = (float **)outputImage->scale;
	scaleY = (float **)outputImage->scale2;
	scaleZ = (float **)outputImage->scale3;
	vxTmp = (float **)outputImage->vxTmp;
	vyTmp = (float **)outputImage->vyTmp;
	vzTmp = (float **)outputImage->vzTmp;
	sxTmp = (float **)outputImage->sxTmp;
	syTmp = (float **)outputImage->syTmp;
	fScale = (float **)outputImage->fScale;
	errorX = (float **)outputImage->errorX;
	errorY = (float **)outputImage->errorY;
	/*
	  Compute feather scale for existing
	*/
	computeScale((float **)vXimage, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0, (double)(-LARGEINT));
	undoNormalization(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, fScale, FALSE);
	/***********************************Loop over images******************************** */
	count = 1;
	fprintf(stderr, "Outside image loop %i\n", HemiSphere);
	tCenter = (outputImage->jd1 + outputImage->jd2 + 1.) * 0.5; /* Added 1 on Dec 1 to avoid .5 day bias */
	tHalfWidth = (outputImage->jd2 - outputImage->jd1) * 0.5;
	for (currentImage = LSImages; currentImage != NULL; currentImage = currentImage->next)
	{
		deltaT = currentImage->matches.jdLate - currentImage->matches.jdEarly; /* time seperation */
		fprintf(stderr, "current ls image weight %f %f %f %f\n", currentImage->weight, currentImage->matches.jdEarly, currentImage->matches.jdLate, deltaT);
		/*  Get bounding box of image*/
		getLSRegion(currentImage, &iMin, &iMax, &jMin, &jMax, outputImage);
		/* Keeps dT from moving outside of interval - see notes - added 10/9/2019 */
		if (currentImage->weight < 0.5)
		{
			iMax = -1;
			jMax = -1;
		}
		/*  **************** Now loop over output grid*****************  */
		if (iMax > 0 && jMax > 0)
		{
			tOffCenter = (currentImage->matches.jdEarly + currentImage->matches.jdLate) * 0.5;
			deltaOffCenter = tOffCenter - tCenter;
			/*
			  Input offsets
			*/
			readLSOffsetsForMosaic(currentImage);
			x0poly = currentImage->matches.x0 + 0.5 * currentImage->matches.nx * currentImage->matches.stepX * currentImage->matches.dx;
			y0poly = currentImage->matches.y0 + 0.5 * currentImage->matches.ny * currentImage->matches.stepY * currentImage->matches.dy;
			vscalex = currentImage->matches.dx; /* space */
			vscalex *= 365.25 / deltaT;			/* time */
			vscaley = currentImage->matches.dy; /* space */
			vscaley *= 365.25 / deltaT;			/* time */
			/*
			  Sigma tie residual scale to pixels
			*/
			if (currentImage->fitResult.sigmaXRes > 0)
				sigmaXTieRes = (currentImage->fitResult.sigmaXRes / currentImage->matches.dx);
			else
				sigmaXTieRes = 0.0;
			if (currentImage->fitResult.sigmaYRes > 0)
				sigmaYTieRes = (currentImage->fitResult.sigmaYRes / currentImage->matches.dy);
			else
				sigmaYTieRes = 0.0;
			if (currentImage->matches.meanSigmaX > 0)
				sigmaXImageAvg = currentImage->matches.meanSigmaX;
			else
				sigmaXImageAvg = 0;
			if (currentImage->matches.meanSigmaY > 0)
				sigmaYImageAvg = currentImage->matches.meanSigmaY;
			else
				sigmaYImageAvg = 0;
			fprintf(stderr, "sigmaX/Y Image Avg %lf %lf SigmaX/Y fit %lf %lf - dT %lf %lf\n",
					sigmaXImageAvg, sigmaYImageAvg, sigmaXTieRes, sigmaYTieRes, deltaOffCenter, currentImage->weight);
			fprintf(stderr, "sig %f %f  svg %f %f\n", sigmaXImageAvg * vscalex, sigmaYImageAvg * vscaley, sigmaXTieRes * vscalex, sigmaYTieRes * vscaley);
			/*
			  This calcuates the tie residual (represent overall error) - the error due to average image stats (will add back on local image stats below)
			  Avoid going negative and assume at least 100th of a pixel of other error. Modified 3/23/20 to change from 0.01 to 0.0001.
			*/
			sigmaQx = 0.05; /* added 3/24/2020 - empirical fix perhaps for quantization/estimation errror */
			sigmaQy = 0.05;
			varXGlobal = max(sigmaXTieRes * sigmaXTieRes - sigmaXImageAvg * sigmaXImageAvg, 0.000) + sigmaQx * sigmaQx;
			varYGlobal = max(sigmaYTieRes * sigmaYTieRes - sigmaYImageAvg * sigmaYImageAvg, 0.000) + sigmaQy * sigmaQy;
			fprintf(stderr, "var %f %f \n", sqrt(varXGlobal) * vscalex, sqrt(varYGlobal) * vscaley);
			for (k = 0; k < 3; k++)
			{
				pX[k] = currentImage->fitResult.pX[k];
				pY[k] = currentImage->fitResult.pY[k];
			}
			for (i = iMin; i < iMax; i++)
			{
				if ((i % 100) == 0)
					fprintf(stderr, "-- %i \n", i);
				y = (outputImage->originY + i * outputImage->deltaY) * MTOKM;
				for (j = jMin; j < jMax; j++)
				{
					x = (outputImage->originX + j * outputImage->deltaX) * MTOKM;
					vx = NODATA;
					vy = NODATA;
					/* Adding skip via continue if no shelf mask indicates no data */
					if (shelfMask != NULL)
						sMask = getShelfMask(shelfMask, x, y);
					if (sMask == NOSOLUTION)
						continue; /* skip rest of loop */
					/* Ok, shelf mask indicates continue */
					if (interpLSdata(x * KMTOM, y * KMTOM, &(currentImage->matches), &offx, &offy, &sx, &sy) == TRUE)
					{
						/*
						   add global and local noise (see above)Updated to use squared value 3/23/20 since value is squared again below
						*/
						computeFitError((x * KMTOM - x0poly), (y * KMTOM - y0poly), currentImage, &sigx2, &sigy2);
						sx2 = varXGlobal + sx * sx + sigx2;
						sy2 = varYGlobal + sy * sy + sigy2;
						sx2 = min(sx2, MAXSIG * MAXSIG);
						sy2 = min(sy2, MAXSIG * MAXSIG);
						/*
						  Add polynomial corrections, in pixels
						*/
						offx -= pX[0] + pX[1] * (x * KMTOM - x0poly) + pX[2] * (y * KMTOM - y0poly);
						offy -= pY[0] + pY[1] * (x * KMTOM - x0poly) + pY[2] * (y * KMTOM - y0poly);
						/*
						  Compute lat lon
						*/
						xytoll(x, y, HemiSphere, &lat, &lon, Rotation);
						/*
						  computeLatScale
						*/
						latscale = xyscale(lat, currentImage->fitResult.proj);
						vx = offx * vscalex * latscale;
						vy = offy * vscaley * latscale;
						/* Variances */
						ex = (vscalex * vscalex) * sx2;
						ey = (vscaley * vscaley) * sy2;
						scX = 1.0 / (ex);
						scY = 1.0 / (ey);
						/* Increment buffers */
						vxTmp[i][j] = (float)vx * scX;
						vyTmp[i][j] = (float)vy * scY;
						sxTmp[i][j] = scX;
						syTmp[i][j] = scY;
						fScale[i][j] = 1.0; /* Value for zero feathering */
						if (outputImage->timeOverlapFlag == TRUE)
						{
							vzTmp[i][j] = (float)(deltaOffCenter * sqrt(scX * scY));
						}
						else
						{
							vzTmp[i][j] = (float)1;
						}
					}
					else
					{
						vxTmp[i][j] = (float)-LARGEINT;
						fScale[i][j] = 0.0;
					}
				} /* j loop */
			}	  /* i loop */
			/*
			  Compute scale array for feathering.
			*/
			if (fl > 0 && (iMax > 0 && jMax > 0))
			{
				computeScaleLS((float **)vxTmp, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0, (double)(-LARGEINT), iMin, iMax, jMin, jMax);
			}
			/*
			  Now sum current result. Falls through if no intersection (iMax&jMax==0)
			*/
			redoNormalization(currentImage->weight, outputImage, iMin, iMax, jMin, jMax, vXimage, vYimage, vZimage, errorX, errorY,
							  scaleX, scaleY, scaleZ, fScale, vxTmp, vyTmp, vzTmp, sxTmp, syTmp, FALSE);
		}
		count++;
	} /* End image loop */
	/**************************END OF MAIN LOOP ******************************/
	/*
	  Adjust scale
	*/
	endScale(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, FALSE);
	free(lineBuf);
	fprintf(outputImage->fpLog, ";\n; Returning from landsat  Mosaic(.c)\n;\n");
}

static void computeFitError(double x, double y, landSatImage *currentImage, double *sigx2, double *sigy2)
{
	double Cx[MAXFITPARAM][MAXFITPARAM], Cy[MAXFITPARAM][MAXFITPARAM];
	int32_t i, l, k;
	double v[3], tmpx[3], tmpy[3];

	v[0] = 1;
	v[1] = x;
	v[2] = y;
	for (k = 0; k < 3; k++)
	{
		for (l = 0; l < 3; l++)
		{
			Cx[k][l] = currentImage->fitResult.Cx[k][l];
			Cy[k][l] = currentImage->fitResult.Cy[k][l];
		}
	}
	for (i = 0; i < 3; i++)
	{
		tmpx[i] = Cx[i][0] * v[0] + Cx[i][1] * v[1] + Cx[i][2] * v[2];
		tmpy[i] = Cy[i][0] * v[0] + Cy[i][1] * v[1] + Cy[i][2] * v[2];
	}
	/*  Error in pixels	*/
	*sigx2 = tmpx[0] * v[0] + tmpx[1] * v[1] + tmpx[2] * v[2];
	*sigy2 = tmpy[0] * v[0] + tmpy[1] * v[1] + tmpy[2] * v[2];
}

static int32_t lsBinlinear(float **X, int32_t im, int32_t jm, double t, double u, double *dx)
{
	double p1, p2, p3, p4;
	p1 = X[im][jm];
	p2 = X[im][jm + 1];
	p3 = X[im + 1][jm + 1];
	p4 = X[im + 1][jm];
	/* Don't use if all 4pts aren't good - should ensure better quality data */
	if (p1 <= (NODATA + 1) || p2 <= (NODATA + 1) || p3 <= (NODATA + 1) || p4 <= (NODATA + 1))
		return (FALSE);
	*dx = (double)((1.0 - t) * (1.0 - u) * p1 + t * (1.0 - u) * p2 + t * u * p3 + (1.0 - t) * u * p4);
	return TRUE;
}

static int32_t interpLSdata(double x, double y, matchResult *matches, double *dx, double *dy, double *sx, double *sy)
{
	double xi, yi, t, u;
	int32_t jm, im;
	*dx = NODATA;
	*dy = NODATA;
	*sx = NODATA;
	*sy = NODATA;
	/*
	  Image coordinates
	*/
	xi = (x - matches->x0) / (matches->dx * (double)matches->stepX);
	yi = (y - matches->y0) / (matches->dy * (double)matches->stepY);
	/* Integer values, rounded down */
	jm = (int32_t)xi;
	im = (int32_t)yi;
	/* Return if out of bounds, use -2 since  im+1,jm+1 index array,means extreme border points neglected (good thing) */
	if (jm < 0 || im < 0 || jm > (matches->nx - 2) || im > (matches->ny - 2))
		return (FALSE);
	/* Fraction of pixel for interpolation */
	t = (float)(xi - (double)jm);
	u = (float)(yi - (double)im);
	/*
	  interpoints
	*/
	if (lsBinlinear(matches->X, im, jm, t, u, dx) == FALSE)
		return (FALSE);
	if (lsBinlinear(matches->Y, im, jm, t, u, dy) == FALSE)
		return (FALSE);
	if (lsBinlinear(matches->sigmaX, im, jm, t, u, sx) == FALSE)
		return (FALSE);
	if (lsBinlinear(matches->sigmaY, im, jm, t, u, sy) == FALSE)
		return (FALSE);
	return (TRUE); /* Return true since interpolation appears to have worked */
}

/***************************LSreadFloatImage*************************************
input a floating point image file, with size nx by ny
***********************************************************************************/
static void readLSOffsetsForMosaic(landSatImage *currentImage)
{
	extern float *fBuf1, *fBuf2, *fBuf3, *fBuf4;
	extern float **fBuf1L, **fBuf2L, **fBuf3L, **fBuf4L;
	extern uint8_t **mask, *maskBuf;
	uint8_t **lsMask;
	char *offXFile;
	lsFit *fitDat;
	int32_t i, j;
	matchResult *matches;
	fitDat = &(currentImage->fitResult);
	matches = &(currentImage->matches);
	/*
	  Read data X data
	*/
	offXFile = lineBuf;
	/* X offsets */
	offXFile[0] = '\0';
	offXFile = strcpy(offXFile, fitDat->matchFile);
	offXFile = strcat(offXFile, ".dx");
	fprintf(stderr, "Read -  X offsets %s %i %i \n", offXFile, matches->nx, matches->ny);
	matches->X = LSreadFloatImageforMosaic(offXFile, matches->nx, matches->ny, fBuf1, fBuf1L);
	/* Y offsets */
	offXFile[0] = '\0';
	offXFile = strcpy(offXFile, fitDat->matchFile);
	offXFile = strcat(offXFile, ".dy");
	matches->Y = LSreadFloatImageforMosaic(offXFile, matches->nx, matches->ny, fBuf2, fBuf2L);
	/* sigma X offsets */
	offXFile[0] = '\0';
	offXFile = strcpy(offXFile, fitDat->matchFile);
	offXFile = strcat(offXFile, ".sx");
	matches->sigmaX = LSreadFloatImageforMosaic(offXFile, matches->nx, matches->ny, fBuf3, fBuf3L);
	/* sigma Y offsets */
	/*fprintf(stderr,"Read -  sigma Y offsets %s %i %i  %u\n",offXFile,matches->nx,matches->ny,(uint32)fBuf4);*/
	offXFile[0] = '\0';
	offXFile = strcpy(offXFile, fitDat->matchFile);
	offXFile = strcat(offXFile, ".sy");
	matches->sigmaY = LSreadFloatImageforMosaic(offXFile, matches->nx, matches->ny, fBuf4, fBuf4L);
	if (currentImage->maskFile != NULL)
	{
		lsMask = LSreadByteImageforMosaic(currentImage->maskFile, matches->nx, matches->ny, maskBuf, mask);
		for (i = 0; i < matches->ny; i++)
			for (j = 0; j < matches->nx; j++)
				if (lsMask[i][j] == 0)
				{
					matches->X[i][j] = NODATA;
					matches->Y[i][j] = NODATA;
					matches->sigmaX[i][j] = NODATA;
					matches->sigmaY[i][j] = NODATA;
				}
	}
}
/*
  Read LS image
 */
static float **LSreadFloatImageforMosaic(char *file, int32_t nx, int32_t ny, float *fBuffer, float **image)
{
	FILE *fp;
	int32_t i;
	/*    malloc space	*/
	for (i = 0; i < ny; i++)
		image[i] = &(fBuffer[i * nx]);
	/*   open file	*/
	fp = openInputFile(file);
	/*    read file	*/
	freadBS(image[0], sizeof(float), (size_t)(nx * ny), fp, FLOAT32FLAG);
	fclose(fp);
	return (image);
}

static uint8_t **LSreadByteImageforMosaic(char *file, int32_t nx, int32_t ny, uint8_t *maskB, uint8_t **m)
{
	FILE *fp;
	int32_t i;
	/* 	   malloc space  */
	for (i = 0; i < ny; i++)
		m[i] = &(maskB[i * nx]);
	/*    open file	*/
	fp = openInputFile(file);
	/*   read file	*/
	freadBS(m[0], sizeof(int8_t), (size_t)(nx * ny), fp, BYTEFLAG);
	fclose(fp);
	return (m);
}

/*
  Find region where image exists
*/
static void getLSRegion(landSatImage *image, int32_t *iMin, int32_t *iMax, int32_t *jMin, int32_t *jMax, outputImageStructure *outputImage)
{
	extern double Rotation;
	double minX, maxX, minY, maxY;
	double pad;
	/*
	  Loop through points to find max and min locations
	*/
	minX = image->matches.x0;
	maxX = image->matches.x0 + (image->matches.nx - 1) * image->matches.dx * (double)(image->matches.stepX);
	minY = image->matches.y0;
	maxY = image->matches.y0 + (image->matches.ny - 1) * image->matches.dy * (double)(image->matches.stepY);
	minX *= MTOKM;
	minY *= MTOKM;
	maxX *= MTOKM;
	maxY *= MTOKM;
	/*
	  Compute i,j min,max with pad
	*/
	pad = 0;
	*iMin = (int)((minY * KMTOM - outputImage->originY - pad) / outputImage->deltaY);
	*jMin = (int)((minX * KMTOM - outputImage->originX - pad) / outputImage->deltaX);
	*iMax = (int)((maxY * KMTOM - outputImage->originY + pad) / outputImage->deltaY);
	*jMax = (int)((maxX * KMTOM - outputImage->originX + pad) / outputImage->deltaX);
	fprintf(stderr, "%i %i %i %i -- ", *iMin, *iMax, *jMin, *jMax);
	*iMin = max(*iMin, 0);
	*jMin = max(*jMin, 0);
	*iMax = min(outputImage->ySize, *iMax);
	*jMax = min(outputImage->xSize, *jMax);
	fprintf(stderr, "%i %i  %i %i\n ", *iMin, *iMax, *jMin, *jMax);
}

double **rDistSaveLS = NULL;
/*
  Compute scale for feathering.
*/
void computeScaleLS(float **inImage, float **scale, int32_t azimuthSize, int32_t rangeSize, float fl, float weight, double minVal, int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax)
{
	extern double **rDistSaveLS;
	double **rDist, rA;
	int32_t j, k, is, il, s1, s2, l1, l2;
	float minV;
	/*
	  Set up featherin - radial distance from center of kernel
	*/
	if (rDistSaveLS == NULL)
	{
		rDistSaveLS = dmatrix(-fl, fl, -fl, fl);
		fillRadialKernel(rDistSaveLS, fl, weight);
	}
	rDist = rDistSaveLS;
	/*
	  Set initial value for scale array
	*/
	for (j = iMin; j < azimuthSize; j++)
		for (k = 0; k < rangeSize; k++)
			scale[j][k] = weight;
	/*
	  Now loop through weights
	*/
	for (j = iMin; j < iMax; j++)
	{
		for (k = jMin; k < jMax; k++)
		{
			/*
			  Now do adding the feathering at valid points.
			*/
			if (inImage[j][k] > minVal)
			{
				minV = 1000.;
				/* Find edge pixels - check all neighbors and if one is non-valid, its an edge */
				s1 = max(0, j - 1);
				s2 = min(azimuthSize - 1, j + 1);
				l1 = max(0, k - 1);
				l2 = min(rangeSize - 1, k + 1);
				for (is = s1; is <= s2; is++)
					for (il = l1; il <= l2; il++)
					{
						minV = FMIN(minV, (float)(inImage[is][il]));
					}
				/*   Adjust scale for border points
					removed the tile border 6/22, since it was causing small glitches on seam boundaries
					***                if( (minV <= minVal) || (j == 0) || (k == 0) ||
					***                    (j==(azimuthSize-1)) || (k ==(rangeSize-1))  ) {   */
				if ((minV <= minVal))
				{
					s1 = max(0, j - fl);
					s2 = min(azimuthSize - 1, j + fl);
					l1 = max(0, k - fl);
					l2 = min(rangeSize - 1, k + fl);
					for (is = s1; is <= s2; is++)
					{
						for (il = l1; il <= l2; il++)
						{
							scale[is][il] =
								FMIN(rDist[is - j][il - k], scale[is][il]);
						} /* il */
					}	  /* is */
				}		  /* End if(minV... */
			}			  /* End if(inp... */
		}				  /* k */
	}					  /* j */
}
