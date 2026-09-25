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
static void computeFootprintPolygon(inputImageStructure *inputImage);

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


	stateV *sv = &(inputImage->sv);
	// If sv start is much greater than data start, the data has probably
	// crossed a day boundary and the statevectors have not. 
	// So add 86400 to time to compensate
	if( (sv->times[1] - cp->sTime) > 50000 )
	{
		cp->sTime += 86400;
		fprintf(stderr, 
			"Adjusting data sec of day from %f to %f to match SV start time of %f\n",
			cp->sTime-86400, cp->sTime, sv->times[1]);
	}
	cp->eTime = cp->sTime + (cp->azSize * inputImage->nAzimuthLooks) / cp->prf;
	cp->pixelToAzimuthPixel = 1. / inputImage->azimuthPixelSize;
	cp->toRangePixel = 1.0 / inputImage->rangePixelSize;
	computeSatHeightNew(cp, inputImage, memMode);
	inputImage->isInit = TRUE;
	inputImage->tolerance = 1e-6;
	computeFootprintPolygon(inputImage);
}

/*
  Compute the frame's 4-corner footprint polygon (x/y km) once per inputImage and
  cache it in inputImage->footprintX/Y, instead of recomputing it on every checkLL()
  call (i.e. once per candidate tie/geocoding point -- up to 500k times for tiepoints,
  and far more for mosaic3d's per-pixel geocoding). The corners are re-ordered by
  angle around their centroid before storing: inputImage->latControlPoints[1..4]/
  lonControlPoints[1..4] (parseInputFile.c parseControlPointsGeoJson()'s
  index={0,0,3,1,2} remap) was chosen to make the old axis-aligned min/max computation
  order-independent, not to trace the corners in perimeter order, which the
  point-in-polygon test in checkLL() requires.
*/
static void computeFootprintPolygon(inputImageStructure *inputImage)
{
	/* The footprint polygon and checkLL's test point must be in the SAME projection;
	   which one does not matter for a containment test, but under -epsg the legacy
	   globals no longer describe the output grid, so use the resolved projection for
	   both.  For a polar run grimpDefaultProj() reproduces the old rot/stdLat exactly,
	   including the 70/71 substitution for the -91 sentinel. */
	const grimpProj *proj = grimpDefaultProj();
	double cx[4], cy[4], ang[4], centX, centY;
	int32_t i, j, idx[4] = {0, 1, 2, 3}, tmp;

	for (i = 0; i < 4; i++)
		llToXYProj(inputImage->latControlPoints[i + 1], inputImage->lonControlPoints[i + 1],
				&cx[i], &cy[i], proj);
	centX = (cx[0] + cx[1] + cx[2] + cx[3]) / 4.0;
	centY = (cy[0] + cy[1] + cy[2] + cy[3]) / 4.0;
	for (i = 0; i < 4; i++)
		ang[i] = atan2(cy[i] - centY, cx[i] - centX);
	for (i = 0; i < 3; i++)
		for (j = i + 1; j < 4; j++)
			if (ang[idx[j]] < ang[idx[i]])
			{
				tmp = idx[i];
				idx[i] = idx[j];
				idx[j] = tmp;
			}
	for (i = 0; i < 4; i++)
	{
		inputImage->footprintX[i] = cx[idx[i]];
		inputImage->footprintY[i] = cy[idx[i]];
	}
}

/* Standard ray-casting point-in-polygon test; polyX/polyY must be a simple
   (non-self-intersecting) ring, either winding order. */
static int32_t pointInPolygon(double px, double py, double *polyX, double *polyY, int32_t n)
{
	int32_t i, j, inside = FALSE;
	for (i = 0, j = n - 1; i < n; j = i++)
	{
		if (((polyY[i] > py) != (polyY[j] > py)) &&
			(px < (polyX[j] - polyX[i]) * (py - polyY[i]) / (polyY[j] - polyY[i]) + polyX[i]))
			inside = !inside;
	}
	return inside;
}

static double distToSegment(double px, double py, double x1, double y1, double x2, double y2)
{
	double dx = x2 - x1, dy = y2 - y1;
	double len2 = dx * dx + dy * dy;
	double t, cx, cy;
	if (len2 < 1.0e-9)
		return sqrt((px - x1) * (px - x1) + (py - y1) * (py - y1));
	t = ((px - x1) * dx + (py - y1) * dy) / len2;
	t = max(0.0, min(1.0, t));
	cx = x1 + t * dx;
	cy = y1 + t * dy;
	return sqrt((px - cx) * (px - cx) + (py - cy) * (py - cy));
}

/*
  Reject candidate points that fall outside the frame's actual footprint. Was an
  axis-aligned bounding box (min/max of the 4 corners, +-50km pad) -- for a long,
  curving, near-polar frame that box clips along a constant-x (or constant-y) line
  wherever the true (curved/rotated) footprint diverges from the box edge, cutting off
  real swath data at one end while still letting through off-swath points elsewhere.
  Replaced with a real point-in-polygon test against the 4 corner points (angle-sorted
  and cached once per inputImage in initllToImageNew()/computeFootprintPolygon() above
  -- this function runs per candidate point, up to 500k times for tiepoints and far
  more for mosaic3d's per-pixel geocoding, so the corner conversion/sort must not be
  redone here), with the same 50km pad applied as a distance-to-boundary tolerance
  instead of a box expansion.
*/
static int32_t checkLL(double lat, double lon, inputImageStructure *inputImage)
{
	/* Same projection as computeFootprintPolygon above -- see the note there. */
	const grimpProj *proj = grimpDefaultProj();
	double x, y, dmin;
	int32_t i;

	llToXYProj(lat, lon, &x, &y, proj);

	if (pointInPolygon(x, y, inputImage->footprintX, inputImage->footprintY, 4))
		return TRUE;

	dmin = 1.0e30;
	for (i = 0; i < 4; i++)
	{
		double d = distToSegment(x, y, inputImage->footprintX[i], inputImage->footprintY[i],
								  inputImage->footprintX[(i + 1) % 4], inputImage->footprintY[(i + 1) % 4]);
		if (d < dmin)
			dmin = d;
	}
	return (dmin <= 50.0); /* pad, km */
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
	double C1, C2, df, dT;
	int32_t tol = 2000;
	int32_t i, n;
	cp = &(inputImage->cpAll);
	sv = &(inputImage->sv);
	/* Avoid extreme values that could give opposite side solution */
	if (checkLL(lat, lon, inputImage) == FALSE)
	{
		*range = -9999.0;
		*azimuth = -9999.0;
		//fprintf(stderr, "RETURNING BAD LL lat %f lon %f\n", lat, lon);
		return;
	}
	/* Refine later */
	if (inputImage->lastTime >= (cp->sTime - 10) && inputImage->lastTime <= (cp->eTime + 10) )
	{
		myTime = inputImage->lastTime;
	}
	else
		myTime = (cp->sTime + cp->eTime) * 0.5;
		//myTime = sv->times[sv->nState / 2];
	llToECEF(lat, lon, h, &xt, &yt, &zt);
	C2 = 0.0; /* Not used for zero dop */
	C1 = 0;
	dT = 10000;  // Force intial state vector selection
	for (i = 0; i < 35; i++)
	{
		// if time jumps more than state vector interval, update (added 3/10/25)
		if(fabs(dT > sv->deltaT))
		{
			n = (int32_t)((myTime - sv->times[1]) / (sv->deltaT) + .5);
			n = min(max(0, n - NUSESTATE/2), sv->nState - NUSESTATE);
		}
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
		dT = df / (C1 + C2);
		myTime -= dT;
		/* Check for convergence */
		if (fabs(dT) < inputImage->tolerance)
			break;
	}
	// Final call
	polintVec(&(sv->times[n]), &(sv->x[n]), &(sv->y[n]), &(sv->z[n]), &(sv->vx[n]), &(sv->vy[n]), &(sv->vz[n]),
			  myTime, &xs, &ys, &zs, &vsx, &vsy, &vsz);
	*range = (sqrt(dot(drx, dry, drz, drx, dry, drz)) - cp->RNear) * cp->toRangePixel;
	*azimuth = ((myTime - cp->sTime) * cp->prf) / inputImage->nAzimuthLooks;
	//fprintf(stderr, "\033[33mllToImageNew: lat %f lon %f h %f range %f azimuth %f time %f\n\033[0m", lat, lon, h, *range, *azimuth, myTime);
	/* The zero-Doppler/range solve above (dot(dr,v)=0, range=|dr|) is symmetric under
	   reflection of the target through the orbital plane, so it cannot by itself
	   distinguish a real target from its mirror image on the wrong side of the
	   ground track -- checkLL()'s bounding box is too loose to catch this for long,
	   high-latitude frames. Disambiguate with the same cross-track sign convention
	   already used to steer the forward geolocation in smlocateZD.c
	   (elook = lookDir*acos(...), so sign(sin(elook)) == sign(lookDir)) and to flip
	   the cross-track baseline component in svBase.c's svBnBp(): ph = V x U (U the
	   satellite's own outward unit position vector) is the cross-track direction, and
	   for a correctly-illuminated target dot(target-satellite, ph) always carries the
	   same sign as lookDir. */
	{
		double normS = sqrt(dot(xs, ys, zs, xs, ys, zs));
		double phx, phy, phz;
		cross(vsx, vsy, vsz, xs / normS, ys / normS, zs / normS, &phx, &phy, &phz);
		double side = dot(drx, dry, drz, phx, phy, phz);
		if ((side < 0 && inputImage->lookDir == RIGHT) || (side > 0 && inputImage->lookDir == LEFT))
		{
			*range = -9999.0;
			*azimuth = -9999.0;
			inputImage->lastTime = myTime;
			return;
		}
	}
	/* allow tol ml pixel buffer in case indexing into slc - allow some tol for calcs like heading - other checks will avoid bad coords */
	if (*range < -tol || *range > (inputImage->rangeSize + tol) || *azimuth < -tol || *azimuth > (inputImage->azimuthSize + tol))
	{
		*range = -9999.0;
		*azimuth = -9999.0;
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
