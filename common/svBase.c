#include "stdio.h"
#include "string.h"
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "cRecipes/nrutil.h"
#include <unistd.h>
/* Buffers for recyclint matrix in az setup */
double **uRows = NULL;
double *uBuffer = NULL;
#define NR_END 1

static void unitTCN(double R[3], double V[3], double T[3], double C[3], double N[3])
{
	double xs, ys, zs, vsx, vsy, vsz;
	double NxV[3];
	double normPos, normNxV;
	/* N -r/|r| */
	normPos = norm(R[0], R[1], R[2]);
	N[0] = -R[0] / normPos;
	N[1] = -R[1] / normPos;
	N[2] = -R[2] / normPos;
	/* C = NxV/|NxV| */
	cross(N[0], N[1], N[2], V[0], V[1], V[2], &(NxV[0]), &(NxV[1]), &(NxV[2]));
	normNxV = norm(NxV[0], NxV[1], NxV[2]);
	C[0] = NxV[0] / normNxV;
	C[1] = NxV[1] / normNxV;
	C[2] = NxV[2] / normNxV;
	/* T = CxN*/
	cross(C[0], C[1], C[2], N[0], N[1], N[2], &(T[0]), &(T[1]), &(T[2]));
}

void svInterpBnBp(inputImageStructure *inputImage, Offsets *offsets, double azimuth, double *bnS, double *bpS)
{
	/* Interpolate Bn, Bp from precomputed values as a function of azimuth offsets coordinate */
	double azTime;
	double mlpAzCenter;
	double floatAzimuth;
	int32_t iAzimuth;
	/* Center of first multi look pixel in sl pixels */
	mlpAzCenter = (inputImage->nAzimuthLooks - 1) * 0.5;
	/* If not alread initialized, the init */
	if (offsets->bnS == NULL || offsets->bpS == NULL)
		svInitBnBp(inputImage, offsets);
	/* Convert multi-look az to single look, then to offset */
	floatAzimuth = (mlpAzCenter + azimuth * inputImage->nAzimuthLooks - offsets->aO) / offsets->deltaA;
	iAzimuth = (int)(floatAzimuth + 0.5);
	/*	iAzimuth=(int)(azimuth+0.5);*/
	iAzimuth = min(max(0, iAzimuth), offsets->na - 1);
	/* Nearest neighbor interpolation */
	*bnS = offsets->bnS[iAzimuth];
	*bpS = offsets->bpS[iAzimuth];
}

/* convert range and time on the lat/lon on the ellipsoid */
static void RTtoLatLon(inputImageStructure *inputImage, double r, double myTime, double *lat, double *lon)
{
	int32_t n;
	stateV *sv;
	double xs, ys, zs, vsx, vsy, vsz;
	sv = &(inputImage->sv);
	n = (int32_t)((myTime - sv->times[1]) / (sv->deltaT) + .5);
	n = min(max(0, n - 2), sv->nState - NUSESTATE);
	polintVec(&(sv->times[n]), &(sv->x[n]), &(sv->y[n]), &(sv->z[n]), &(sv->vx[n]), &(sv->vy[n]), &(sv->vz[n]), myTime, &xs, &ys, &zs, &vsx, &vsy, &vsz);
	smlocateZD(xs * MTOKM, ys * MTOKM, zs * MTOKM, vsx * MTOKM, vsy * MTOKM, vsz * MTOKM, r * MTOKM, lat, lon, (double)(inputImage->lookDir), 0.0);
}

typedef struct
{
	double r;
	double a;
} svdData;
/* Number of points in each dimension to fit azimuth offset function */
#define NAPTS 50
/* set to  4 for first order poly a+b*az +c rg +d rg*az or 5 to add quad in range ... + e *rg^2 */
#define AZORD 5

/* order 1 or order 1+rg^2 - see AZORD  - func for svdfit */
static void order1Fit(void *x, int32_t i, double *afunc, int32_t ma)
{
	svdData *xx;
	double x1, y1;
	xx = (svdData *)x;
	afunc[1] = 1;
	afunc[2] = xx[i].r;
	afunc[3] = xx[i].a;
	afunc[4] = xx[i].r * xx[i].a;
	if (AZORD > 4)
		afunc[5] = xx[i].r * xx[i].r;
	return;
}
/* use fit result from svInitAzParams to  compute az off and ml range/azimuth location */
double svAzOffset(inputImageStructure *inputImage, Offsets *offsets, double range, double azimuth)
{
	double r, a, rPoly, aPoly;
	double dA, a0;
	/*
	  Initialize if needed
	 */
	if (offsets->azInit != TRUE)
		svInitAzParams(inputImage, offsets);
	/*
	  rsescale to rg,az  km referenced to image  enter
	*/
	r = inputImage->par.rn + range * inputImage->nRangeLooks * inputImage->par.slpR;
	a = azimuth * inputImage->nAzimuthLooks * inputImage->par.slpA;
	rPoly = (r - inputImage->par.rc) * MTOKM;
	a0 = (inputImage->azimuthSize * inputImage->nAzimuthLooks * inputImage->par.slpA) * 0.5;
	aPoly = (a - a0) * MTOKM;
	dA = offsets->azFit[0] + offsets->azFit[1] * rPoly + offsets->azFit[2] * aPoly + offsets->azFit[3] * rPoly * aPoly;
	if (AZORD == 5)
		dA += offsets->azFit[4] * rPoly * rPoly;
	return dA;
}

double **dmatrixRecycle(int32_t nrl, int32_t nrh, int32_t ncl, int32_t nch, double **mR, double *mBuf)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch], same as numerical recipies */
{
	int32_t i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	double **m;

	/* allocate pointers to rows */
	/*m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));*/
	m = mR;
	if (!m)
		nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	/* m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));*/
	m[nrl] = (double *)mBuf;
	if (!m[nrl])
		nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for (i = nrl + 1; i <= nrh; i++)
		m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

/*
   Fit a first order 2-d poly to points computed from state vector. Result is used to get sv derived az offsets
   Can do first order or first order + rg^2 - see AZORD above
 */
void svInitAzParams(inputImageStructure *inputImage, Offsets *offsets)
{
	inputImageStructure inputImage2;
	svdData ra[NAPTS * NAPTS + 1];
	double da[NAPTS * NAPTS + 1], sigX[NAPTS * NAPTS + 1];
	double myFit[AZORD + 1];
	double **u, **v, chisq, vSpace[(AZORD + 1) * (AZORD + 1)], *vRows[(AZORD + 1)], w[(AZORD + 1)];
	double ReHBuffer[MAXAZIMUTHLINES], ReH;
	double r0, t0, r, t, dr, dt;
	double lat, lon, sig, mean, da1;
	double range1, range2, azimuth1, azimuth2;
	int32_t ma = AZORD;
	int32_t i, j, n, nPts;
	/* Setup memory recycling */
	if (uRows == NULL)
	{
		fprintf(stderr, "Initializing matrix buffer for svInitAzParams");
		uRows = (double **)malloc((size_t)((NAPTS * NAPTS + NR_END) * sizeof(double *)));
		uBuffer = (double *)malloc((size_t)((NAPTS * NAPTS * AZORD + NR_END) * sizeof(double)));
	}
	/* Increments in rg ad time */
	dr = (inputImage->par.rf - inputImage->par.rn) / (NAPTS - 1);
	dt = (inputImage->azimuthSize * inputImage->nAzimuthLooks / inputImage->par.prf) / (NAPTS - 1);
	/* time at image center */
	t0 = inputImage->cpAll.sTime + (inputImage->azimuthSize * inputImage->nAzimuthLooks / inputImage->par.prf) * 0.5;
	/* read the second image file to get the second set of sv  */
	parseInputFile(offsets->geo2, &inputImage2);
	/* Use temp buff */
	inputImage2.cpAll.ReH = &(ReHBuffer[0]); /* temp storage so init works ok */
	inputImage2.isInit = -1111;				 /* Signal not to allocate memory */
	initllToImageNew(&inputImage2);
	/* Copy the sv to offsets */
	memcpy(&(offsets->sv2), &(inputImage2.sv), sizeof(inputImage2.sv));
	offsets->dt1t2 = inputImage->cpAll.sTime - inputImage2.cpAll.sTime;
	/* Loop over time, and range to compute azimuth shift to produce point used for polynomial fit.*/
	n = 1;
	r0 = inputImage->par.rc;
	for (i = 0; i < NAPTS; i++)
	{
		t = inputImage->cpAll.sTime + i * dt;
		for (j = 0; j < NAPTS; j++)
		{
			r = inputImage->par.rn + j * dr;							 /* absolute range */
			RTtoLatLon(inputImage, r, t, &lat, &lon);					 /* lat/lon in first image */
			llToImageNew(lat, lon, 0, &range1, &azimuth1, inputImage);	 /* convert back to range1,azimuth1 */
			llToImageNew(lat, lon, 0, &range2, &azimuth2, &inputImage2); /* same for second */
			/* only do if point overlaps */
			if (range2 > 0 && range2<inputImage2.rangeSize & azimuth2> 0 && azimuth2 < inputImage2.azimuthSize)
			{
				da[n] = (azimuth2 - azimuth1) * inputImage->nAzimuthLooks * inputImage->par.slpA;
				ra[n].r = (r - r0) * MTOKM;
				sigX[n] = 1.;
				ra[n].a = (t - t0) * inputImage->par.prf * inputImage->par.slpA * MTOKM;
				n++;
			}
		}
	}
	nPts = n - 1;
	if (nPts < 2 * NAPTS)
		error("svInitAzParams : insufficient points to fit");
	/* Updated to recycle memory, 9/15/22 */
	u = dmatrixRecycle(1, nPts, 1, ma, uRows, uBuffer);
	v = dmatrixRecycle(1, ma, 1, ma, vRows, vSpace);
	svdfit((void *)ra, da, sigX, nPts, myFit, ma, u, v, w, &chisq, &order1Fit);
	for (i = 0; i < AZORD; i++)
	{
		offsets->azFit[i] = myFit[i + 1];
	}
	fprintf(stderr, "\nSV Fit for Offsets sqrt(x2/nPts) %f %i\n", sqrt(chisq / nPts), nPts);
	offsets->azInit = TRUE;
}

void svInitBnBp(inputImageStructure *inputImage, Offsets *offsets)
{
	/* Compute bn/bp as a function of along track coordinate for the offsets. Arrays are indexed with offset coordinates */
	inputImageStructure inputImage2;
	double offA, offR, azTime;
	double bTCN[3], bnS, bpS;
	double Re, H, thetaC;
	double ReHBuffer[MAXAZIMUTHLINES], ReH;
	int32_t i;
	Re = inputImage->cpAll.Re;
	H = inputImage->par.H;
	if (inputImage->cpAll.ReH != NULL)
		ReH = inputImage->cpAll.ReH[(int)(inputImage->azimuthSize) / 2];
	else
		error("missing reH in svInitBnBp");
	thetaC = thetaRReZReH(inputImage->cpAll.RCenter, (Re + 0), ReH);
	/* Parse the second input image */
	parseInputFile(offsets->geo2, &inputImage2);
	inputImage2.cpAll.ReH = &(ReHBuffer[0]); /* temp storage so init works ok */
	inputImage2.isInit = -1111;				 /* Signal not to allocate memory */
	initllToImageNew(&inputImage2);
	/* Copy offsets from second image to the second set of state vectors for the offsets */
	memcpy(&(offsets->sv2), &(inputImage2.sv), sizeof(inputImage2.sv));
	offsets->dt1t2 = inputImage->cpAll.sTime - inputImage2.cpAll.sTime;
	/* Compute the constant offset between images */
	svOffsets(inputImage, &inputImage2, offsets, &offR, &offA);
	/* Modified Sept 15, 22 to use offset coord indexing rather than multi-look coords to save memory */
	offsets->bnS = (float *)malloc(sizeof(float) * offsets->na);
	offsets->bpS = (float *)malloc(sizeof(float) * offsets->na);
	offsets->rConst = offR;
	/* Compute baseline along track */
	for (i = 0; i < offsets->na; i++)
	{
		/* Convert to offset coords to azTime. aO is the center of the first offset estimate */
		azTime = inputImage->cpAll.sTime + (offsets->aO + i * offsets->deltaA) / inputImage->par.prf;
		/* Compute baseline */
		svBaseTCN(azTime, offsets->dt1t2, &(inputImage->sv), &(offsets->sv2), bTCN);
		svBnBp(azTime, thetaC, offsets->dt1t2, &(inputImage->sv), &(offsets->sv2), &bnS, &bpS);
		/* Save values */
		offsets->bnS[i] = bnS;
		offsets->bpS[i] = bpS;
	}
}

void svBnBp(double myTime, double theta, double dt1t2, stateV *sv1, stateV *sv2, double *bn, double *bp)
{
	/* Compute the normal and perpendicular components of baseline using TCN solution */
	double bTCN[3];
	svBaseTCN(myTime, dt1t2, sv1, sv2, bTCN);
	*bp = bTCN[2] * cos(theta) + bTCN[1] * sin(theta);
	*bn = bTCN[1] * cos(theta) - bTCN[2] * sin(theta);
}

/* Compute the bulk offsets between images directly from the state vectors. */
void svOffsets(inputImageStructure *image1, inputImageStructure *image2, Offsets *offsets, double *cnstR, double *cnstA)
{
	int32_t i;
	double range1, range2, azimuth1, azimuth2;
	double r1, r2, a1, a2;
	double bn, bp;
	double Re, H, ReH;
	double bTCN[3];
	double azTime, thetaC, theta, thetaD;
	double deltaRo;
	double deltaRb, deltaA;
	int32_t nPts;
	Re = image1->cpAll.Re;
	H = image1->par.H;
	if (image1->cpAll.ReH != NULL)
		ReH = image1->cpAll.ReH[(int)(image1->azimuthSize) / 2];
	else
		ReH = Re + H;
	thetaC = thetaRReZReH(image1->cpAll.RCenter, (Re + 0), (Re + H));
	*cnstR = 0.0;
	*cnstA = 0.0;
	nPts = 0;
	for (i = 0; i < 5; i++)
	{
		llToImageNew(image1->latControlPoints[i], image1->lonControlPoints[i], 0, &range1, &azimuth1, image1);
		llToImageNew(image1->latControlPoints[i], image1->lonControlPoints[i], 0, &range2, &azimuth2, image2);
		Re = earthRadius(image1->latControlPoints[i] * DTOR, EMINOR, EMAJOR) * KMTOM;
		azTime = image1->cpAll.sTime + azimuth1 * image1->nAzimuthLooks / image1->par.prf;
		svBaseTCN(azTime, offsets->dt1t2, &(image1->sv), &(offsets->sv2), bTCN);
		svBnBp(azTime, thetaC, offsets->dt1t2, &(image1->sv), &(offsets->sv2), &bn, &bp);
		r1 = image1->par.rn + range1 * image1->nRangeLooks * image1->par.slpR;
		r2 = image2->par.rn + range2 * image2->nRangeLooks * image2->par.slpR;
		a1 = azimuth1 * image1->nAzimuthLooks * image1->par.slpA;
		a2 = azimuth2 * image2->nAzimuthLooks * image2->par.slpA;
		deltaRb = r2 - r1;
		/* This is the offset for one image to the other */
		deltaRo = (range2 - range1) * image1->nRangeLooks * image1->par.slpR;
		deltaA = a2 - a1;
		/* this ensures that the azimuth difference relative to the early time for the first image */
		if (azimuth1 < 100)
		{
			*cnstR += (deltaRo - deltaRb);
			*cnstA += deltaA;
			nPts++;
		}
	}
	if (nPts > 0)
	{
		*cnstR = *cnstR / (double)nPts;
		*cnstA = *cnstA / (double)nPts;
	}
	else
	{
		*cnstR = 0.0;
		*cnstA = 0.0;
	}
	fprintf(stderr, "\033[1;31m dR,dA %f %f %i\033[0m\n", *cnstR, *cnstA, nPts);
}

void svBaseTCN(double myTime, double dt1t2, stateV *sv1, stateV *sv2, double bTCN[3])
{
	/* Compute baseline in TCN coordinates, given time, timedifference, and state vectors */
	double R1[3], V1[3], R2[3], V2[3];
	double T[3], C[3], N[3];
	double b[3];
	double myTime2, dt, C1;
	int32_t i, j, n;
	/* Interp state vectors for R & V */
	n = (int32_t)((myTime - sv1->times[1]) / (sv1->deltaT) + .5);
	n = min(max(0, n - 2), sv1->nState - NUSESTATE);
	/* Position in master image */
	polintVec(&(sv1->times[n]), &(sv1->x[n]), &(sv1->y[n]), &(sv1->z[n]), &(sv1->vx[n]), &(sv1->vy[n]), &(sv1->vz[n]),
			  myTime, &(R1[0]), &(R1[1]), &(R1[2]), &(V1[0]), &(V1[1]), &(V1[2]));
	/* TCN Vector */
	unitTCN(R1, V1, T, C, N);
	/* Time for second image */
	myTime2 = myTime - dt1t2;
	for (i = 0; i < 5; i++)
	{
		/* interp second state vector */
		n = (int32_t)((myTime2 - sv2->times[1]) / (sv2->deltaT) + .5);
		n = min(max(0, n - 2), sv2->nState - NUSESTATE);
		polintVec(&(sv2->times[n]), &(sv2->x[n]), &(sv2->y[n]), &(sv2->z[n]), &(sv2->vx[n]), &(sv2->vy[n]), &(sv2->vz[n]),
				  myTime2, &(R2[0]), &(R2[1]), &(R2[2]), &(V2[0]), &(V2[1]), &(V2[2]));
		b[0] = R2[0] - R1[0];
		b[1] = R2[1] - R1[1];
		b[2] = R2[2] - R1[2];
		C1 = dot(V2[0], V2[1], V2[2], T[0], T[1], T[2]);
		dt = dot(b[0], b[1], b[2], T[0], T[1], T[2]) / C1;
		myTime2 -= dt;
		/* Iterate until convergence or max 5 tries */
		if (fabs(dt) < 1.0e-11)
			break;
	}
	bTCN[0] = dot(b[0], b[1], b[2], T[0], T[1], T[2]);
	bTCN[1] = dot(b[0], b[1], b[2], C[0], C[1], C[2]);
	bTCN[2] = dot(b[0], b[1], b[2], N[0], N[1], N[2]);
}

/* ********************************** Test Code *****************************/
void svTest(inputImageStructure *inputImage, Offsets *offsets)
{
	inputImageStructure inputImage2;
	double offA, offR, azTime;
	double bTCN1[3], bTCN2[3], bnS, bpS;
	double Re, H, thetaC;
	double ReHBuffer[MAXAZIMUTHLINES], ReH;
	double dt;
	int32_t i;
	Re = inputImage->cpAll.Re;
	dt = 0.01;
	H = inputImage->par.H;
	if (inputImage->cpAll.ReH != NULL)
		ReH = inputImage->cpAll.ReH[(int)(inputImage->azimuthSize) / 2];
	else
		error("missing reH in svInitBnBp");
	thetaC = thetaRReZReH(inputImage->cpAll.RCenter, (Re + 0), ReH);
	//fprintf(stderr, "%f\n", thetaC * RTOD);
	/* Parse the second input image */
	parseInputFile(offsets->geo2, &inputImage2);
	inputImage2.cpAll.ReH = &(ReHBuffer[0]); /* temp storage so init works ok */
	inputImage2.isInit = -1111;				 /* Signal not to allocate memory */
	initllToImageNew(&inputImage2);
	memcpy(&(offsets->sv2), &(inputImage2.sv), sizeof(inputImage2.sv));
	offsets->dt1t2 = inputImage->cpAll.sTime - inputImage2.cpAll.sTime;
	fprintf(stderr, "times %f %f %f\n", inputImage->cpAll.sTime, inputImage2.cpAll.sTime, offsets->dt1t2);
	svOffsets(inputImage, &inputImage2, offsets, &offR, &offA);
	offsets->bnS = (float *)malloc(sizeof(float) * inputImage->azimuthSize);
	offsets->bpS = (float *)malloc(sizeof(float) * inputImage->azimuthSize);
	offsets->rConst = offR;
	for (i = 0; i < inputImage->azimuthSize; i++)
	{
		azTime = inputImage->cpAll.sTime + i * inputImage->nAzimuthLooks / inputImage->par.prf;
		svBaseTCN(azTime - dt * .5, offsets->dt1t2, &(inputImage->sv), &(offsets->sv2), bTCN1);
		svBaseTCN(azTime + dt * .5, offsets->dt1t2, &(inputImage->sv), &(offsets->sv2), bTCN2);
		if (i % 1000 == 0)
		{
			fprintf(stderr, "%lf %lf  %f %f\n", bTCN1[1], bTCN1[2], bTCN2[1], bTCN2[2]);
		}
		fprintf(stdout, "%10.8f %10.8f\n", (bTCN2[1] - bTCN1[1]) / dt, (bTCN2[2] - bTCN1[2]) / dt);
	}
}
