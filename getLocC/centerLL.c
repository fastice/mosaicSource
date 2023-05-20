#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "math.h"
#include "mosaicSource/common/common.h"
#include "cRecipes/nrutil.h"
/*
  Find coordinates for image center in ground range
*/
/* added 8/13/16 to handle large numbers of state vectors. This is so polint32_t uses a 5 point32_t interpolation - make sure not change with updating where its used */
#define NUSESTATE 5

void centerLL(SARData *sarD, stateV *sv, int32_t nla, double *lat, double *lon, double deltaT)
{
	double sTime; /* Starting time in second of day of image */
	double pTime; /* Starting time ''               of pixel */
	int32_t i;	  /* LCV */
	int32_t n;
	double x, y, z, vx, vy, vz; /* Interpolated pos/vel */
	double da;
	double rSlant;
	double ReSat, Re, sch;
	int32_t endofday;
	double rhoNear, rhoFar, rhoMid, rMid;
	double e;
	int32_t nMLA;
	/*
	  Init state vectors
	*/
	nMLA = sarD->nSlpA / nla;
	endofday = FALSE;
	for (i = 1; i <= sv->nState; i++)
	{
		//fprintf(stderr, "sv %f\n", sv->times[i]);
		if (sv->times[i] > 86400)
			endofday = TRUE;
	}
	/*
	  Starting time for image
	*/
	sTime = sarD->hr * 3600. + sarD->min * 60.0 + sarD->sec + deltaT;
	da = ((double)(nMLA)-1.0) * 0.5;
	pTime = sTime + (double)(da * nla) / sarD->prf; /* Center Time */
	fprintf(stderr, "s %f m %f --- %f %i  %f %f\n", sTime, pTime, da, nla, sarD->prf, (double)(da * nla) / sarD->prf);
	/* Wrap time into next  day if near end of sv->times and time straddles day boundary */
	if (endofday == TRUE && sTime < 70000)
		pTime += 86400;
	/*
	   Interpolate state vectors at image center - polint is 1..N, so n=0, will use n=1 as first
	*/
	if (sv->nState > NUSESTATE)
	{
		n = (int32_t)((pTime - sv->times[1]) / (sv->times[2] - sv->times[1]) + .5);
		n = min(max(0, n - 2), sv->nState - NUSESTATE);
	}
	else
		n = 0;
	polintVec(&(sv->times[n]), &(sv->x[n]), &(sv->y[n]), &(sv->z[n]), &(sv->vx[n]), &(sv->vy[n]), &(sv->vz[n]), pTime, &x, &y, &z, &vx, &vy, &vz);
	/*
	  Get slant range for center ground range
	*/
	ReSat = sqrt((x * x + y * y + z * z) / ((x * x + y * y) / pow(EMAJOR * KMTOM, 2.0) + z * z / pow(EMINOR * KMTOM, 2.0)));
	sch = sqrt(x * x + y * y + z * z) - ReSat;
	/* Re target - iterarte 3 times just to get Earth radius correctly*/
	for (i = 0; i < 3; i++)
	{
		Re = earthRadius(*lat * DTOR, EMINOR, EMAJOR) * KMTOM;
		rhoNear = acos((pow(sarD->rn, 2.0) - pow(ReSat + sch, 2.0) - pow(Re, 2.0)) / (-2.0 * Re * (ReSat + sch)));
		rhoFar = acos((pow(sarD->rf, 2.0) - pow(ReSat + sch, 2.0) - pow(Re, 2.0)) / (-2.0 * Re * (ReSat + sch)));
		rSlant = (sarD->rn + sarD->rf) * 0.5;
		rhoMid = rhoRReZReH(rSlant, Re, (ReSat + sch));
		/*  Compute lat/lon values*/
		smlocateZD(x * MTOKM, y * MTOKM, z * MTOKM, vx * MTOKM, vy * MTOKM, vz * MTOKM, rSlant * MTOKM, lat, lon, sarD->lookDir, 0.);
	}
}
