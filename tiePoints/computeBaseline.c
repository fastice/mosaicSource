#include <math.h>
#include "mosaicSource/common/common.h"
#include "tiePoints.h"
#include "cRecipes/nrutil.h"
#include <stdlib.h>
/*
  Estimate baseline parameters.
*/

#define NPARAMSEST 4

void baselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
void noRampBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
void dBpBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
void dBpQBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
void bpBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
void bpbnBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
void bpdBpBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
void bnbpdBpBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma);

double lambda;

void computeBaseline(tiePointsStructure *tiePoints,
					 inputImageStructure inputImage)
{
	double Re, H, RNear, thetaC, dr, rOffset;
	double *a; /* Solution for params */
	double nParams;
	double theta, thetaD, z, r0, Cij;
	extern double lambda;
	double twok;
	double **v, **u, chisq, *sigB, **Cp, **C6;
	double ReH;
	double *y, *sig, *w, *chsq;
	double bn, bp, bSq;
	double Bn, dBn, dBnQ, Bp, dBp, dBpQ;
	double deltaExact, deltaApprox;
	double varP, sigP, meanP, xtmp;
	int32_t nData, ma;
	int32_t azimuth;
	modelValues *x;
	int32_t i, i1, k, j, npts, l1, l2;
	int32_t pIndex[7];
	conversionDataStructure *cP;

	lambda = inputImage.par.lambda;
	fprintf(stderr, "wavelength %f\n", lambda);
	if (tiePoints->quadB == TRUE)
		nParams = 6;
	else if (tiePoints->bpFlag == TRUE)
		nParams = 1;
	else if (tiePoints->bnbpFlag == TRUE)
		nParams = 2;
	else if (tiePoints->bnbpdBpFlag == TRUE)
		nParams = 3;
	else if (tiePoints->bpdBpFlag == TRUE)
		nParams = 2;
	else
		nParams = 4;
	a = dvector(1, nParams);
	if (inputImage.isInit != TRUE)
		initllToImageNew(&inputImage);
	cP = &(inputImage.cpAll);
	ReH = getReH(cP, &inputImage, (inputImage.azimuthSize) / 2);
	twok = 2.0 * 2.0 * PI / lambda;
	Re = tiePoints->Re;
	RNear = tiePoints->RNear;
	thetaC = thetaRReZReH(cP->RCenter, (Re + 0.), ReH);
	fprintf(stderr, "------------------+++ RNear %f %f %f %f\n", RNear, thetaC * RTOD, tiePoints->thetaC * RTOD, ReH);
	/*
	  Range comp params
	*/
	dr = inputImage.rangePixelSize;
	rOffset = 0.0;
	/*
	  Init arrays for least squares fit.
	*/
	npts = 0;
	nData = tiePoints->npts;
	for (i = 0; i < nData; i++)
		if (fabs(tiePoints->phase[i]) < 1.0E6)
			npts++;
	fprintf(stderr, "Npoints %i\n", npts);
	ma = nParams;
	x = (modelValues *)malloc((npts + 1) * sizeof(modelValues));
	y = dvector(1, npts);
	sig = dvector(1, npts);
	u = dmatrix(1, npts, 1, ma);
	v = dmatrix(1, ma, 1, ma);
	w = dvector(1, ma);
	sigB = dvector(1, ma);
	C6 = dmatrix(1, 6, 1, 6);
	Cp = dmatrix(1, ma, 1, ma);
	/*
	  Loop twice, first using flattening value of bsq, and then value
	  from first fit. Should easily converge with just two iterations
	  unless very large error in initial estimate.
	  Changed to use base estimate only. Leaving loop for possible
	  later modification.

	  added second loop on 4/28/14 to iterate on delta and bsq
	*/
	sigP = 10.0; /* Use for first try */
	for (k = 0; k <= 2; k++)
	{
		j = 0;
		varP = 0;
		meanP = 0;
		if (k == 0)
		{
			Bn = tiePoints->BnCorig;
			Bp = tiePoints->BpCorig;
			dBn = tiePoints->dBnorig;
			dBp = tiePoints->dBporig;
			dBnQ = tiePoints->dBnQorig;
			dBpQ = tiePoints->dBpQorig;
		}
		for (i = 0; i < nData; i++)
		{
			if (fabs(tiePoints->phase[i]) < 1.0E6)
			{ /* Use only good points */
				i1 = j + 1;
				z = tiePoints->z[i];
				r0 = RNear + rOffset + tiePoints->r[i] * dr;
				azimuth = (int)tiePoints->a[i];
				ReH = getReH(cP, &inputImage, azimuth);
				theta = acos((r0 * r0 + (ReH) * (ReH)-pow(Re + z, 2.0)) / (2.0 * (ReH)*r0));
				thetaD = theta - thetaC;
				x[i1].x = tiePoints->x[i] / (inputImage.azimuthSize * inputImage.azimuthPixelSize);
				x[i1].thetaD = thetaD;

				bn = Bn + dBn * x[i1].x + dBnQ * x[i1].x * x[i1].x;
				bp = Bp + dBp * x[i1].x + dBpQ * x[i1].x * x[i1].x;
				bSq = bn * bn + bp * bp;
				deltaApprox = -bn * sin(thetaD) - bp * cos(thetaD) + bSq * 0.5 / r0 - (pow(tiePoints->delta[i], 2.0) / (2.0 * r0));
				deltaApprox = -bn * sin(thetaD) - bp * cos(thetaD) + bSq * 0.5 / r0 - (pow(deltaApprox, 2.0) / (2.0 * r0));
				y[i1] = tiePoints->phase[i] - twok * bSq / (2.0 * r0) + pow(deltaApprox, 2.0) * twok / (2.0 * r0);
				xtmp = -twok * sin(thetaD) * (Bn + dBn * x[i1].x + dBnQ * (x[i1].x * x[i1].x)) - twok * cos(thetaD) * (Bp + dBp * x[i1].x + dBpQ * (x[i1].x * x[i1].x));
				varP += (y[i1] - xtmp) * (y[i1] - xtmp);
				meanP += (y[i1] - xtmp);
				/*
				  Subtract off known terms for bpFlag
				*/
				if (tiePoints->bpFlag == TRUE)
				{
					y[i1] -= -twok * sin(thetaD) * bPoly(tiePoints->BnCorig, tiePoints->dBnorig, tiePoints->dBnQorig, x[i1].x);
					y[i1] -= -twok * cos(thetaD) * bPoly(0, tiePoints->dBporig, tiePoints->dBpQorig, x[i1].x);
				}
				else if (tiePoints->bnbpFlag == TRUE)
				{
					y[i1] -= -twok * sin(thetaD) * bPoly(0.0, tiePoints->dBnorig, tiePoints->dBnQorig, x[i1].x);
					y[i1] -= -twok * cos(thetaD) * bPoly(0.0, tiePoints->dBporig, tiePoints->dBpQorig, x[i1].x);
				}
				else if (tiePoints->bpdBpFlag == TRUE)
				{
					y[i1] -= -twok * sin(thetaD) * bPoly(tiePoints->BnCorig, tiePoints->dBnorig, tiePoints->dBnQorig, x[i1].x);
				}
				else if (tiePoints->bnbpdBpFlag == TRUE)
				{
					y[i1] -= -twok * sin(thetaD) * bPoly(0.0, tiePoints->dBnorig, tiePoints->dBnQorig, x[i1].x);
				}
				sig[i1] = sigP;
				j++;
			} /* End if */
		}	  /* End for i */
		if (i1 < 4)
			error("tiepoints: Insufficient Number (%i)  of Valid tie points \n", i1);

		if (tiePoints->dBpFlag == TRUE)
		{
			if (tiePoints->quadB == TRUE)
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &dBpQBaselineCoeffs);
			else if (tiePoints->bpFlag == TRUE)
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &bpBaselineCoeffs);
			else if (tiePoints->bnbpFlag == TRUE)
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &bpbnBaselineCoeffs);
			else if (tiePoints->bpdBpFlag == TRUE)
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &bpdBpBaselineCoeffs);
			else if (tiePoints->bnbpdBpFlag == TRUE)
			{
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &bnbpdBpBaselineCoeffs);
			}
			else
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &dBpBaselineCoeffs);
		}
		else
		{
			if (tiePoints->noRamp == TRUE)
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &noRampBaselineCoeffs);
			else
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &baselineCoeffs);
		}
		svdvar(v, ma, w, Cp);
		if (tiePoints->dBpFlag == FALSE)
			a[3] = twok * a[3] / (inputImage.azimuthSize * inputImage.nAzimuthLooks);
		if (tiePoints->quadB == TRUE)
		{
			Bn = a[1], Bp = a[2];
			dBn = a[4];
			dBp = a[3];
			dBnQ = a[6];
			dBpQ = a[5];
			pIndex[1] = 1;
			pIndex[2] = 2;
			pIndex[3] = 4;
			pIndex[4] = 3;
			pIndex[5] = 6;
			pIndex[6] = 5;
		}
		else if (tiePoints->bpFlag == TRUE)
		{
			Bn = tiePoints->BnCorig;
			Bp = a[1];
			dBn = tiePoints->dBnorig;
			dBp = tiePoints->dBporig;
			dBnQ = tiePoints->dBnQorig;
			dBpQ = tiePoints->dBpQorig;
			pIndex[1] = 0;
			pIndex[2] = 1;
			pIndex[3] = 0;
			pIndex[4] = 0;
			pIndex[5] = 0;
			pIndex[6] = 0;
		}
		else if (tiePoints->bnbpFlag == TRUE)
		{
			Bn = a[1];
			Bp = a[2];
			dBn = tiePoints->dBnorig;
			dBp = tiePoints->dBporig;
			dBnQ = tiePoints->dBnQorig;
			dBpQ = tiePoints->dBpQorig;
			pIndex[1] = 1;
			pIndex[2] = 2;
			pIndex[3] = 0;
			pIndex[4] = 0;
			pIndex[5] = 0;
			pIndex[6] = 0;
		}
		else if (tiePoints->bpdBpFlag == TRUE)
		{
			Bn = tiePoints->BnCorig;
			Bp = a[1];
			dBn = tiePoints->dBnorig;
			dBp = a[2];
			dBnQ = tiePoints->dBnQorig;
			dBpQ = tiePoints->dBpQorig;
			pIndex[1] = 0;
			pIndex[2] = 1;
			pIndex[3] = 0;
			pIndex[4] = 2;
			pIndex[5] = 0;
			pIndex[6] = 0;
		}
		else if (tiePoints->bnbpdBpFlag == TRUE)
		{
			Bn = a[1];
			Bp = a[2];
			dBn = tiePoints->dBnorig;
			dBp = a[3];
			dBnQ = tiePoints->dBnQorig;
			dBpQ = tiePoints->dBpQorig;
			pIndex[1] = 1;
			pIndex[2] = 2;
			pIndex[3] = 0;
			pIndex[4] = 3;
			pIndex[5] = 0;
			pIndex[6] = 0;
		}
		else
		{
			Bn = a[1];
			Bp = a[2];
			dBn = a[4];
			dBp = a[3];
			pIndex[1] = 1;
			pIndex[2] = 2;
			pIndex[3] = 4;
			pIndex[4] = 3;
			pIndex[5] = 0;
			pIndex[6] = 0;
		}
		varP = varP / npts;
		meanP = meanP / npts;
		sigP = sqrt(varP - meanP * meanP);
		fprintf(stderr, "%lf %lf %lf %lf %lf %lf %lf %lf \n", Bn, Bp, dBn, dBp, dBnQ, dBpQ, meanP, sigP);

	} /* end k*/

	fprintf(stdout, ";\n; Ntiepoints/Ngiven used= %i/%i\n;\n", npts, nData);
	/* use sigma=10 as in above */
	fprintf(stdout, "; X2 %f \n", chisq);
	fprintf(stdout, ";*  sigma*sqrt(X2/n)= %f \n", sigP * sqrt(chisq / (double)npts));
	fprintf(stdout, "; Covariance Matrix \n;");
	for (l1 = 1; l1 <= 6; l1++)
	{
		fprintf(stdout, ";* C_%1i ", l1);
		for (l2 = 1; l2 <= 6; l2++)
		{
			Cij = 0;
			if (pIndex[l1] > 0 && pIndex[l1] <= ma && pIndex[l2] > 0 && pIndex[l2] <= ma)
				Cij = Cp[pIndex[l1]][pIndex[l2]];
			fprintf(stdout, " %10.6e ", Cij);
		}
		fprintf(stdout, "\n");
	}
	/*
	  Output results
	*/
	if (tiePoints->dBpFlag == TRUE)
		fprintf(stdout, ";\n; Estimated Baseline\n; Bn,Bp,dBn,dBp\n;\n");
	else
		fprintf(stdout, ";\n; Estimated Baseline\n; Bn,Bp,dBn,omegaA\n;\n");

	if (tiePoints->quadB == TRUE)
	{
		fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f\n&\n", a[1], a[2], a[4], a[3], a[6], a[5]);
	}
	else if (tiePoints->bpFlag == TRUE)
	{
		fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f\n&\n",
				tiePoints->BnCorig, a[1], tiePoints->dBnorig, tiePoints->dBporig, tiePoints->dBnQorig, tiePoints->dBpQorig);
	}
	else if (tiePoints->bnbpFlag == TRUE)
	{
		fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f\n&\n", a[1], a[2], tiePoints->dBnorig, tiePoints->dBporig,
				tiePoints->dBnQorig, tiePoints->dBpQorig);
	}
	else if (tiePoints->bpdBpFlag == TRUE)
	{
		fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f\n&\n",
				tiePoints->BnCorig, a[1], tiePoints->dBnorig, a[2], tiePoints->dBnQorig, tiePoints->dBpQorig);
	}
	else if (tiePoints->bnbpdBpFlag == TRUE)
	{
		fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f\n&\n",
				a[1], a[2], tiePoints->dBnorig, a[3], tiePoints->dBnQorig, tiePoints->dBpQorig);
	}
	else
		fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f\n&\n", a[1], a[2], a[4], a[3]);

	return;
}

void baselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	double twok;
	double thetaD;
	extern double lambda;
	modelValues xx, *xy;
	double xAz;

	twok = 4.0 * PI / lambda;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;
	thetaD = xx.thetaD;

	afunc[1] = -twok * sin(thetaD);
	afunc[2] = -twok * cos(thetaD);
	afunc[3] = twok * xAz;
	afunc[4] = -twok * xAz * sin(thetaD);
	return;
}

void bpBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	extern double lambda;
	double twok;
	double thetaD;

	modelValues xx, *xy;
	double xAz;

	twok = 4.0 * PI / lambda;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;
	thetaD = xx.thetaD;

	afunc[1] = -twok * cos(thetaD);
	return;
}

void noRampBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	/*
	  This version multiplies ramp term by cos(thetad)
	*/
	extern double lambda;
	double twok;
	double thetaD;

	modelValues xx, *xy;
	double xAz;

	twok = 4.0 * PI / lambda;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;

	thetaD = xx.thetaD;
	afunc[1] = -twok * sin(thetaD);
	afunc[2] = -twok * cos(thetaD);
	afunc[3] = twok * xAz * cos(thetaD);
	afunc[4] = -twok * xAz * sin(thetaD);
	return;
}
void dBpBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	/*
	  This version is for computing dBp
	*/
	extern double lambda;
	double twok;
	double thetaD;

	modelValues xx, *xy;
	double xAz;

	twok = 4.0 * PI / lambda;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;

	thetaD = xx.thetaD;
	afunc[1] = -twok * sin(thetaD);
	afunc[2] = -twok * cos(thetaD);

	afunc[3] = -twok * xAz * cos(thetaD);
	afunc[4] = -twok * xAz * sin(thetaD);
	return;
}

void bnbpdBpBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	/*
	  This version is for computing dBp
	*/
	extern double lambda;
	double twok;
	double thetaD;

	modelValues xx, *xy;
	double xAz;

	twok = 4.0 * PI / lambda;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;

	thetaD = xx.thetaD;
	afunc[1] = -twok * sin(thetaD);
	afunc[2] = -twok * cos(thetaD);
	afunc[3] = -twok * xAz * cos(thetaD);
	return;
}

void bpdBpBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	/*
	  This version is for computing dBp
	*/
	extern double lambda;
	double twok;
	double thetaD;

	modelValues xx, *xy;
	double xAz;

	twok = 4.0 * PI / lambda;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;

	thetaD = xx.thetaD;
	afunc[1] = -twok * cos(thetaD);
	afunc[2] = -twok * xAz * cos(thetaD);
	return;
}

void bpbnBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	/*
	  This version is for computing dBp
	*/
	extern double lambda;
	double twok;
	double thetaD;

	modelValues xx, *xy;
	double xAz;

	twok = 4.0 * PI / lambda;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;

	thetaD = xx.thetaD;
	afunc[1] = -twok * sin(thetaD);
	afunc[2] = -twok * cos(thetaD);

	return;
}

void dBpQBaselineCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	/*
	  This version is for computing dBp
	*/
	extern double lambda;
	double twok;
	double thetaD;

	modelValues xx, *xy;
	double xAz;

	twok = 4.0 * PI / lambda;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;

	thetaD = xx.thetaD;
	afunc[1] = -twok * sin(thetaD);
	afunc[2] = -twok * cos(thetaD);

	afunc[3] = -twok * xAz * cos(thetaD);
	afunc[4] = -twok * xAz * sin(thetaD);

	afunc[5] = -twok * xAz * xAz * cos(thetaD);
	afunc[6] = -twok * xAz * xAz * sin(thetaD);

	return;
}
