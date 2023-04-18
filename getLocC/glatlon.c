#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "math.h"
#include "mosaicSource/common/common.h"
#include "cRecipes/nrutil.h"
/*
  Find coordinates for slant range array.
*/

static double **naSpace = NULL, **nrSpace = NULL;
static double **latSpace = NULL, **lonSpace = NULL;
/* added 8/13/16 to handle large numbers of state vectors. This is so polint uses a 5 point interpolation - make sure not change with updating where its used */
#define NUSESTATE 5
void glatlon(SARData *sarD, stateV *sv, int32_t nlr, int32_t nla, double ***lat1, double ***lon1, int32_t ma, int32_t mr, double deltaT)
{
	extern double **naSpace, **nrSpace, **latSpace, **lonSpace;
	double sTime; /* Starting time in second of day of image */
	double pTime; /* Starting time ''               of pixel */
	int32_t i, j; /* LCV */
	int32_t n;
	double x, y, z, vx, vy, vz; /* Interpolated pos/vel */
	double **lat, **lon, **na, **nr;
	double da, dr;
	double rSlant;
	int32_t endofday;
	double e;
	int32_t nMLA, nMLR;
	/*
	  Check size
	*/
	nMLA = sarD->nSlpA / nla;
	nMLR = sarD->nSlpR / nlr;
	if (ma < 2 || mr < 2)
		error("\n*** glatlon : invalid size ***\n");
	/*
	  Init state vectors
	*/
	endofday = FALSE;
	for (i = 1; i <= sv->nState; i++)
	{
		if (sv->times[i] > 86400)
			endofday = TRUE;
	}
	/*
	  Set up space
	*/
	if (latSpace == NULL)
		latSpace = dmatrix(0, ma - 1, 0, mr - 1);
	if (lonSpace == NULL)
		lonSpace = dmatrix(0, ma - 1, 0, mr - 1);
	if (naSpace == NULL)
		naSpace = dmatrix(0, ma - 1, 0, mr - 1);
	if (nrSpace == NULL)
		nrSpace = dmatrix(0, ma - 1, 0, mr - 1);
	na = naSpace;
	nr = nrSpace;
	lat = latSpace;
	lon = lonSpace;
	*lat1 = lat;
	*lon1 = lon;
	/*
	  Set up array where points desired
	*/
	da = (double)(nMLA - 1) / (double)(ma - 1);
	dr = (double)(nMLR - 1) / (double)(mr - 1);
	for (i = 0; i < ma; i++)
	{
		for (j = 0; j < mr; j++)
		{
			/* multilook pixels */
			na[i][j] = (double)(i * da);
			nr[i][j] = (double)(j * dr);
		}
	}
	/*
	  Starting time for image
	*/
	sTime = sarD->hr * 3600. + sarD->min * 60.0 + sarD->sec + deltaT;
	/*
	  Compute lat/lon values
	*/
	for (i = 0; i < ma; i++)
	{
		/* Time */
		for (j = 0; j < mr; j++)
		{
			pTime = sTime + (double)(na[i][j] * nla) / sarD->prf + deltaT;
			if (sv->nState > NUSESTATE)
			{
				n = (int32_t)((pTime - sv->times[1]) / (sv->times[2] - sv->times[1]) + .5);
				n = min(max(0, n - 2), sv->nState - NUSESTATE);
			}
			else
				n = 0;
			/* Wrap time into next  day if near end of times and time straddles day boundary */
			if (endofday == TRUE && sTime < 7000)
				pTime += 86400;
			rSlant = sarD->rn + nr[i][j] * sarD->slpR * nlr; /* Slant range */
			polint(&(sv->times[n]), &(sv->x[n]), NUSESTATE, pTime, &x, &e);
			polint(&(sv->times[n]), &(sv->y[n]), NUSESTATE, pTime, &y, &e);
			polint(&(sv->times[n]), &(sv->z[n]), NUSESTATE, pTime, &z, &e);
			polint(&(sv->times[n]), &(sv->vx[n]), NUSESTATE, pTime, &vx, &e);
			polint(&(sv->times[n]), &(sv->vy[n]), NUSESTATE, pTime, &vy, &e);
			polint(&(sv->times[n]), &(sv->vz[n]), NUSESTATE, pTime, &vz, &e);
			smlocateZD(x * MTOKM, y * MTOKM, z * MTOKM, vx * MTOKM, vy * MTOKM, vz * MTOKM, rSlant * MTOKM, &(lat[i][j]), &(lon[i][j]), sarD->lookDir, 0.0);
		}
	} /* End for(i=... */
}
