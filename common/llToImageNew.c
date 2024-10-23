#include "stdio.h"
#include "string.h"
#include <stdlib.h>
#include <math.h>
/*#include "mosaicSource/common/common.h"*/
#include "common.h"
#include "cRecipes/nrutil.h"
#include <unistd.h>
/* added 8/13/16 to handle large numbers of state vectors. This is so polintt uses a 5 point interpolation - make sure not change with updating where its used */

static void computeSatHeightNew(conversionDataStructure *cp, inputImageStructure *inputImage, int32_t memMode);

void initllToImageNew(inputImageStructure *inputImage)
{
	extern int32_t llConserveMem;
	SARData *par;
	conversionDataStructure *cp;
	double latc, lonc, xc, yc, x, y, lon;
	int32_t Hem;
	int32_t i;
	int32_t memMode;

	/* Need to remove this eventually */
	par = &(inputImage->par);
	/* Cnversion Parameters  - some of these are historical and may not get used */
	cp = &(inputImage->cpAll);
	cp->RNear = par->rn;
	cp->RFar = cp->RNear + (inputImage->rangeSize - 1) * inputImage->rangePixelSize;
	cp->RCenter = (cp->RNear + cp->RFar) * 0.5;
	cp->azOff = 0;
	cp->azSize = inputImage->azimuthSize;
	cp->rSize = inputImage->rangeSize;
	/* Use center point to use compute Re */
	latc = inputImage->latControlPoints[0];
	lonc = inputImage->lonControlPoints[0];
	inputImage->minLat = 1000;
	inputImage->maxLat = -1000;
	inputImage->minLon = 1000;
	inputImage->maxLon = -1000;
	for (i = 1; i < 5; i++)
	{
		lon = inputImage->lonControlPoints[i];
		if (lon > 180.)
			lon -= 360.0;
		inputImage->minLat = min(inputImage->minLat, inputImage->latControlPoints[i]);
		inputImage->maxLat = max(inputImage->maxLat, inputImage->latControlPoints[i]);
		inputImage->minLon = min(inputImage->minLon, lon);
		inputImage->maxLon = max(inputImage->maxLon, lon);
	}
	cp->Re = earthRadius(latc * DTOR, EMINOR, EMAJOR) * KMTOM;
	memMode = llConserveMem;
	if (inputImage->isInit == -1111)
		memMode = 0;
	else
		cp->ReH = NULL;
	cp->prf = par->prf;
	cp->sTime = par->hr * 3600. + par->min * 60.0 + par->sec;
	cp->eTime = cp->sTime + (cp->azSize * inputImage->nAzimuthLooks) / cp->prf;
	cp->pixelToAzimuthPixel = 1. / inputImage->azimuthPixelSize;
	cp->toRangePixel = 1.0 / inputImage->rangePixelSize;
	computeSatHeightNew(cp, inputImage, memMode);
	inputImage->isInit = TRUE;
	inputImage->tolerance = 1e-6;
}

static int32_t checkLL(double lat, double lon, inputImageStructure *inputImage)
{
	if (lon > 180.)
		lon -= 360.;
	
	if (lat > (inputImage->maxLat + 0.5) || lat < (inputImage->minLat - 0.5) ||
		lon > (inputImage->maxLon + 0.5) || lon < (inputImage->minLon - 0.5)) {
		//fprintf(stderr, "%f %f %f %f %f %f\n", lat, lon,inputImage->minLat, inputImage->maxLat,inputImage->minLon, inputImage->maxLon );
		return FALSE;
	}
	return TRUE;
}

/* static double lastTime=0.0;*/
/*
  Geolocation algorithm for converting lat,lon,h to range, azimuth coordinates in multi look coordinates.
  Base on technique used in JPL/Caltech ISCE
 */
void llToImageNew(double lat, double lon, double h, double *range, double *azimuth, inputImageStructure *inputImage)
{
	/* extern double lastTime; */
	SARData *par;
	conversionDataStructure *cp; /* Conversion params from input image */
	stateV *sv;					 /* statevectors from input image */
	double sTime, myTime;
	double xt, yt, zt;
	double xs, ys, zs, vsx, vsy, vsz;
	double xs1, ys1, zs1, vsx1, vsy1, vsz1;
	double drx, dry, drz;
	double RNear, rgPixSize;
	double C1, C2, df;
	int32_t tol = 2000;
	int32_t i, n;
	cp = &(inputImage->cpAll);
	sv = &(inputImage->sv);
	/* Avoid extreme values that could give opposite side solution */
	if (checkLL(lat, lon, inputImage) == FALSE)
	{
		*range = -9999.0;
		*azimuth = -9999.0;
		return;
	}
	/* Refine later */
	if (inputImage->lastTime >= (cp->sTime - 10) && inputImage->lastTime <= (cp->eTime + 10))
	{
		myTime = inputImage->lastTime;
	}
	else
		myTime = sv->times[sv->nState / 2];
	llToECEF(lat, lon, h, &xt, &yt, &zt);
	C2 = 0.0; /* Not used for zero dop */
	n = (int32_t)((myTime - sv->times[1]) / (sv->deltaT) + .5);
	n = min(max(0, n - NUSESTATE/2), sv->nState - NUSESTATE);
	for (i = 0; i < 35; i++)
	{
		/* Interpolate postion and velocity */
		polintVec(&(sv->times[n]), &(sv->x[n]), &(sv->y[n]), &(sv->z[n]), &(sv->vx[n]), &(sv->vy[n]), &(sv->vz[n]),
				  myTime, &xs, &ys, &zs, &vsx, &vsy, &vsz);
		/* dr */
		drx = xt - xs;
		dry = yt - ys;
		drz = zt - zs;
		/* setup correction */
		df = dot(drx, dry, drz, vsx, vsy, vsz);
		C1 = -dot(vsx, vsy, vsz, vsx, vsy, vsz);
		myTime -= df / (C1 + C2);
		/* Check for convergence */
		if (fabs(df / C1) < inputImage->tolerance)
			break;
	}
	// Final call 
	polintVec(&(sv->times[n]), &(sv->x[n]), &(sv->y[n]), &(sv->z[n]), &(sv->vx[n]), &(sv->vy[n]), &(sv->vz[n]),
			  myTime, &xs, &ys, &zs, &vsx, &vsy, &vsz);
	*range = (sqrt(dot(drx, dry, drz, drx, dry, drz)) - cp->RNear) * cp->toRangePixel;
	*azimuth = ((myTime - cp->sTime) * cp->prf) / inputImage->nAzimuthLooks;
	/* allow tol ml pixel buffer in case indexing into slc - allow some tol for calcs like heading - other checks will avoid bad coords */
	if (*range < -tol || *range > (inputImage->rangeSize + tol) || *azimuth < -tol || *azimuth > (inputImage->azimuthSize + tol))
	{
		*range = -9999.0;
		*azimuth = -9999.0;
	} else {
	//fprintf(stderr,"%d %f %f %f %f %f %f %e\n", i, myTime, inputImage->lastTime, *range, *azimuth, lat, lon, fabs(df/C1));
}
	inputImage->lastTime = myTime;
}

void llToECEF(double lat, double lon, double h, double *x, double *y, double *z)
{
	double f, latCos, latSin;
	double F2, C, S;
	double dtor;
	dtor = (2.0 * (double)3.141592653589793) / 360.0;
	latCos = cos(lat * dtor);
	latSin = sin(lat * dtor);
	f = -(EMINOR / EMAJOR - 1.0);
	F2 = (1.0 - f) * (1.0 - f);
	C = (double)1.0 / sqrt(latCos * latCos + F2 * latSin * latSin);
	S = C * F2;
	*x = (EMAJOR * 1000.0 * C + h) * latCos * cos(lon * dtor);
	*y = (EMAJOR * 1000.0 * C + h) * latCos * sin(lon * dtor);
	*z = (EMAJOR * 1000.0 * S + h) * latSin;
}

/*
  Compute Re + H for Sat along track
*/
static void computeSatHeightNew(conversionDataStructure *cp, inputImageStructure *inputImage, int32_t memMode)
{
	int32_t satHpoint;
	extern char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
	SARData *par;
	double sTime, pTime;
	double xs[MXST + 1], ys[MXST + 1], zs[MXST + 1]; /* Position state vectors */
	double times[MXST + 1];
	char *buf;
	double x, y, z, e;
	double Re, RNear, RFar, rho, dz;
	double dLatN, myLatN, dLatF, myLatF;
	int32_t az;
	int32_t memType;
	int32_t endofday;
	int32_t i, n;

	satHpoint = 0; /* Unlike earlier program, there are not multiple conversion params, so start at 0 */
	Re = cp->Re;
	RNear = cp->RNear;
	RFar = cp->RFar;
	par = &(inputImage->par);
	if (memMode == 1234)
	{
		/*
		  Kluge added to seperate ascending descending images
		*/
		if (inputImage->passType == DESCENDING)
			memType = 0; /* Default choose by asc/desc */
		else
			memType = 1;
		if (inputImage->memChan == MEM1)
			memType = 0; /* Overide if mem type gets set */
		else if (inputImage->memChan == MEM2)
			memType = 1;
		if (memType == 0)
			buf = Dbuf2;
		else
			buf = Abuf2;
		cp->rNear = (double *)&(buf[satHpoint]);
		satHpoint += sizeof(double) * cp->azSize;
		cp->rFar = (double *)&(buf[satHpoint]);
		satHpoint += sizeof(double) * cp->azSize;
		cp->ReH = (double *)&(buf[satHpoint]);
		satHpoint += sizeof(double) * cp->azSize;
		if (satHpoint > MAXADBUF)
			error("satHpoint exceeds buffer size\n");
	}
	else if (memMode == 999)
	{
		fprintf(stderr, "not conserving memory 125t46\n"); /* not converving, and not explicitly set not to converse with 999*/
		cp->rNear = (double *)malloc((size_t)(sizeof(double) * cp->azSize));
		cp->rFar = (double *)malloc((size_t)(sizeof(double) * cp->azSize));
		cp->ReH = (double *)malloc((size_t)(sizeof(double) * cp->azSize));
	}
	else
	{
		fprintf(stderr, "Using existing buffers\n"); /* not converving, and not explicitly set not to converse with 999*/
	}

	for (i = 1; i <= inputImage->sv.nState; i++)
	{
		xs[i] = inputImage->sv.x[i];
		ys[i] = inputImage->sv.y[i];
		zs[i] = inputImage->sv.z[i];
		times[i] = inputImage->sv.times[i];
		endofday = FALSE;
		if (times[i] > 86400)
			endofday = TRUE;
	}
	sTime = par->hr * 3600. + par->min * 60.0 + par->sec;
	for (i = 0; i < cp->azSize; i++)
	{
		az = (cp->azOff + i);
		pTime = sTime + az * inputImage->nAzimuthLooks / par->prf;
		if (endofday == TRUE && sTime < 7000)
			pTime += 86400;
		if (inputImage->sv.nState > NUSESTATE)
		{
			n = (int32_t)((pTime - times[1]) / (times[2] - times[1]) + .5);
			n = min(max(0, n - 2), inputImage->sv.nState - NUSESTATE);
		}
		else
			n = 0;
		polint(&(times[n]), &(xs[n]), NUSESTATE, pTime, &x, &e);
		polint(&(times[n]), &(ys[n]), NUSESTATE, pTime, &y, &e);
		polint(&(times[n]), &(zs[n]), NUSESTATE, pTime, &z, &e);
		cp->ReH[i] = sqrt(x * x + y * y + z * z);
		/* note works for asc/desc */
	}
}
