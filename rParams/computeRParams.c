#include <math.h>
#include "mosaicSource/common/common.h"
#include "rparams.h"
#include "cRecipes/nrutil.h"
#include <stdlib.h>
/*
  Estimate baseline parameters.
*/

/*#define NPARAMSEST 4*/

void rParamsCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
void fewPoints();
void rParamsCoeffsQuad(void *x, int32_t i, double *afunc, int32_t ma);

void rbnbpdBpParamsCoeffs(void *x, int32_t i, double *afunc, int32_t ma);

void rbpdBpParamsCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
void constOnlyParamsCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
void BpCorrectOnlyParamsCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
void rParamsCoeffsQuadCorrect(void *x, int32_t i, double *afunc, int32_t ma);

static double bnQZero(tiePointsStructure *tiePoints, double x)
{
	return bPoly(0.0, tiePoints->dBnorig, tiePoints->dBnQorig, x);
}

static double bnQ(tiePointsStructure *tiePoints, double x)
{
	return bPoly(tiePoints->BnCorig, tiePoints->dBnorig, tiePoints->dBnQorig, x);
}

static double bpQ(tiePointsStructure *tiePoints, double x)
{
	return bPoly(tiePoints->BpCorig, tiePoints->dBporig, tiePoints->dBpQorig, x);
}

static double bpQZero(tiePointsStructure *tiePoints, double x)
{
	return bPoly(0.0, tiePoints->dBporig, tiePoints->dBpQorig, x);
}

void computeRParams(tiePointsStructure *tiePoints, inputImageStructure inputImage, char *baseFile, Offsets *offsets)
{
	double Re, H, RNear, thetaC, dr, rOffset, ReH;
	double *a; /* Solution for params */
	double nParams;
	double theta, thetaD, zSp, r0, Cij;
	double **v, **u, chisq;
	double *y, *sig, *w, *chsq, *sigB, **Cp, **C6;
	int32_t nData, ma;
	double bnO, bpO, deltaO;
	modelValues *x;
	int32_t azimuth;
	int32_t i, i1, k, j, npts, l1, l2, kPrint;
	conversionDataStructure *cP;
	double delta;
	double Bn, dBn, dBnQ, Bp, dBp, dBpQ;
	double bnS, bpS;
	double bn, bp, cnst;
	double varP, sigP, meanP, xtmp;
	double azTime;
	double bTCN[3];
	double bnFix, bpFix;
	int32_t pIndex[7];
	double weightSum, tmp;
	double t;
	nParams = 4.;
	/* Determine number of parameters in the fit */
	if (tiePoints->bnbpdBpFlag == TRUE)
		nParams = 3.;
	else if (tiePoints->quadB == TRUE)
		nParams = 6.;
	else if (tiePoints->bpdBpFlag == TRUE)
		nParams = 2.;
	else if (tiePoints->constOnlyFlag == TRUE)
		nParams = 1.;
	/* Overide any ofther flags for deltab cases */
	if (tiePoints->deltaB == DELTABCONST)
		nParams = 1;
	if (tiePoints->deltaB == DELTABQUAD)
		nParams = 6;
	fprintf(stderr, "*** nParams = %f %i\n", nParams, (int)tiePoints->constOnlyFlag);
	a = dvector(1, nParams);
	/* Added 6/12/07 to adjust ReH along track */
	initllToImageNew(&inputImage);
	cP = &(inputImage.cpAll);
	Re = cP->Re;
	H = inputImage.par.H;
	ReH = getReH(cP, &inputImage, (inputImage.azimuthSize) / 2);
	RNear = cP->RNear;
	/*    Re = tiePoints->Re;     H = tiePoints->H;*/
	RNear = tiePoints->RNear;
	/* Theta C uses H from geodat, which is consistent with that used to define baseline */
	thetaC = thetaRReZReH(cP->RCenter, (Re + 0), (ReH));
	fprintf(stderr, "------++---------RNear %f %f %f %f %f %f\n", RNear, ReH, thetaC * RTOD, Re, H, cP->RCenter);
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
	fprintf(stderr, "nData %i\n", nData);
	for (i = 0; i < nData; i++)
		if (fabs(tiePoints->phase[i]) < 1.0E6)
			npts++;
	ma = nParams;
	fprintf(stderr, "ma %i\n", ma);
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
	*/
	kPrint = 2;
	sigP = 1.0; /* Use for first try */
	cnst = tiePoints->cnstR;
	for (k = 0; k <= kPrint; k++)
	{
		if (k == 0)
		{
			if (tiePoints->deltaB == DELTABNONE)
			{
				Bn = tiePoints->BnCorig;
				Bp = tiePoints->BpCorig;
				dBn = tiePoints->dBnorig;
				dBp = tiePoints->dBporig;
				dBnQ = tiePoints->dBnQorig;
				dBpQ = tiePoints->dBpQorig;
			}
			else
			{
				Bn = 0.0;
				Bp = 0.0;
				dBn = 0.0;
				dBp = 0.0;
				dBnQ = 0.0;
				dBpQ = 0.0;
				bnFix = 0.0;
				bpFix = 0.0;
			}
		}
		varP = 0.0;
		meanP = 0.0;
		fprintf(stderr, "K- %i %f %f %f %f %f %f %f\n", k, Bn, Bp, dBn, dBp, dBnQ, dBpQ, cnst);
		j = 0;
		/* Compute mean weight for renormalization */
		weightSum = 0.0;
		for (i = 0; i < nData; i++)
			if (fabs(tiePoints->phase[i]) < 1.0E6)
				weightSum += tiePoints->weight[i];
		/*
			Loop to setup solution
		*/
		for (i = 0; i < nData; i++)
		{
			if (fabs(tiePoints->phase[i]) < 1.0E6)
			{ /* Use only good points */
				i1 = j + 1;
				zSp = tiePoints->z[i];
				r0 = RNear + rOffset + tiePoints->r[i] * dr;
				azimuth = (int)tiePoints->a[i];
				ReH = getReH(cP, &inputImage, azimuth);
				theta = thetaRReZReH(r0, (Re + zSp), ReH);
				thetaD = theta - thetaC;
				x[i1].x = tiePoints->x[i] / (inputImage.azimuthSize * inputImage.azimuthPixelSize);
				x[i1].thetaD = thetaD;
				/* Baseline line or SV */
				if (tiePoints->deltaB == DELTABNONE)
				{
					bn = bPoly(Bn, dBn, dBnQ, x[i1].x);
					bp = bPoly(Bp, dBp, dBpQ, x[i1].x);
				}
				else
				{
					azTime = cP->sTime + azimuth * inputImage.nAzimuthLooks / inputImage.par.prf;
					svBaseTCN(azTime, offsets->dt1t2, &(inputImage.sv), &(offsets->sv2), bTCN);
					svBnBp(azTime, thetaC, offsets->dt1t2, &(inputImage.sv), &(offsets->sv2), &bnS, &bpS);
					bnFix = bPoly(Bn, dBn, dBnQ, x[i1].x);
					bpFix = bPoly(Bp, dBp, dBpQ, x[i1].x);
					bp = bpS + bpFix;
					bn = bnS + bnFix;
					cnst = tiePoints->cnstR;
				}
				tiePoints->bsq[i] = bp * bp + bn * bn;
				/*
				  This is removing the non-linear term in Equation 7, Jglac 1996, to keep solution linear. It removes
				  the approximate delta^2 based on original solution
				*/
				deltaO = -bn * sin(thetaD) - bp * cos(thetaD) + tiePoints->bsq[i] / (2. * r0);
				y[i1] = tiePoints->phase[i] - tiePoints->bsq[i] / (2 * r0) + deltaO * deltaO / (2.0 * r0);
				/* Topographic contribution in mosaicker, Equation 7 jglac - sol of quad eq*/
				xtmp = (sqrt(pow(r0, 2.0) - 2.0 * r0 * (bn * sin(thetaD) + bp * cos(thetaD)) + bn * bn + bp * bp) - r0 + cnst);
				if (tiePoints->deltaB == DELTABNONE)
				{
					y[i1] -= -cos(thetaD) * tiePoints->BpCorig;
					/* Remove non linear terms + fixed Bp */
					xtmp -= -cos(thetaD) * tiePoints->BpCorig + (tiePoints->bsq[i] / (2 * r0) - deltaO * deltaO / (2.0 * r0));
				}
				else
				{
					xtmp -= (tiePoints->bsq[i] / (2 * r0) - deltaO * deltaO / (2.0 * r0));
				}
				/* Sum mean and variance */
				varP += (y[i1] - xtmp) * (y[i1] - xtmp);
				meanP += (y[i1] - xtmp);
				/*	if(fabs((y[i1]-xtmp)) > 1300) fprintf(stderr,"%i %f %f %f %f %f %f %f\n",i1,y[i1],xtmp,deltaO,bn,bp,r0,cnst   );*/
				/*
				  Subtract known terms for bpFlag
				*/
				if (tiePoints->deltaB == DELTABNONE)
				{
					if (tiePoints->bnbpdBpFlag == TRUE)
						y[i1] -= -sin(thetaD) * bnQZero(tiePoints, x[i1].x);
					else if (tiePoints->bpdBpFlag == TRUE)
						y[i1] -= -sin(thetaD) * bnQ(tiePoints, x[i1].x);
					else if (tiePoints->constOnlyFlag == TRUE)
						y[i1] -= -sin(thetaD) * bnQ(tiePoints, x[i1].x) - cos(thetaD) * bpQZero(tiePoints, x[i1].x);
				}
				else
				{
					y[i1] -= -sin(thetaD) * bnS - cos(thetaD) * bpS + tiePoints->cnstR;
				}
				/* Updated 12/17/21:
				Add capability to weight sigmas to emphasize certain points (e.g., rock w=1) and de-emphasize others (ice that could change w=10).
				The noise may be the same for all points, but the weights help limit points where we think there could be subtle biases. To try and keep
				the mean noise level to that of the data, the normalization will keep the mean sigma consistent with the data. This is not much of a
				problem if most of the points have the same weight. But in cases where say there are 10% low weight points and 90% high weight, this will drive the
				sigmas on the low weight (good) to unreasonably low errors. Either way, it shouldn't affect the solution, but it could skew the covariance matrix
				and hence the baseline error estimates. To avoid this situation, sigmas cannot drop below 0.1.
				*/
				sig[i1] = max(sigP * tiePoints->weight[i] * npts / weightSum, 0.1 * sigP);
				j++;
			} /* End if */
		}	  /* End for i */
		varP = varP / (double)npts;
		meanP = meanP / (double)npts;
		sigP = sqrt(varP - meanP * meanP);
		fprintf(stderr, "k= %i v= %lf  sig= %lf mean = %lf\n", k, varP, sigP, meanP);
		if (i1 < nParams)
		{
			fprintf(stderr, "rParams: Insufficient Number (%i)  of Valid tie points \n", i1);
			fewPoints();
		}
		if (k == kPrint)
		{
			fprintf(stdout, ";\n; Estimated Baseline\n; Bn,Bp,dBn,dBp\n;\n");
			fprintf(stdout, ";\n; Ntiepoints/Ngiven used= %i/%i\n;\n", npts, nData);
		}

		if (tiePoints->deltaB == DELTABNONE)
		{
			if (tiePoints->bnbpdBpFlag == TRUE)
			{ /*    solve for bn,bp,dBp */
				if (npts < 3)
					fewPoints();
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &rbnbpdBpParamsCoeffs);
				Bn = a[1];
				Bp = tiePoints->BpCorig;
				dBn = tiePoints->dBnorig;
				dBp = a[2];
				dBnQ = 0., dBpQ = 0.;
				cnst = a[3];
				pIndex[1] = 1;
				pIndex[2] = 0;
				pIndex[3] = 2;
				pIndex[4] = 0;
				pIndex[5] = 0;
				pIndex[6] = 3;
			}
			else if (tiePoints->bpdBpFlag == TRUE)
			{ /*	  Solve for bp and dBp only 			*/
				if (npts < 2)
					fewPoints();
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &rbpdBpParamsCoeffs);
				Bn = tiePoints->BnCorig;
				Bp = tiePoints->BpCorig;
				dBn = tiePoints->dBnorig;
				dBp = a[1];
				dBnQ = 0., dBpQ = 0.;
				cnst = a[2];
				pIndex[1] = 0;
				pIndex[2] = 0;
				pIndex[3] = 1;
				pIndex[4] = 0;
				pIndex[5] = 0;
				pIndex[6] = 2;
			}
			else if (tiePoints->constOnlyFlag == TRUE)
			{
				if (npts < 2)
					fewPoints();
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &constOnlyParamsCoeffs);
				Bn = tiePoints->BnCorig;
				Bp = tiePoints->BpCorig;
				dBn = tiePoints->dBnorig;
				dBp = tiePoints->dBporig;
				dBnQ = 0., dBpQ = 0.;
				cnst = a[1];
				pIndex[1] = 0;
				pIndex[2] = 0;
				pIndex[3] = 0;
				pIndex[4] = 0;
				pIndex[5] = 0;
				pIndex[6] = 1;
			}
			else if (tiePoints->quadB == TRUE)
			{
				if (npts < 6)
					fewPoints();
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &rParamsCoeffsQuad);
				Bn = a[1];
				Bp = tiePoints->BpCorig;
				dBn = a[3];
				dBp = a[2];
				dBnQ = a[6], dBpQ = a[5];
				cnst = a[4];
				pIndex[1] = 1;
				pIndex[2] = 3;
				pIndex[3] = 2;
				pIndex[4] = 6;
				pIndex[5] = 5;
				pIndex[6] = 4;
			}
			else
			{ /* 	 SOLVE FOR ALL PARAMETERS */
				if (npts < 4)
					fewPoints();
				svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &rParamsCoeffs);
				Bn = a[1];
				Bp = tiePoints->BpCorig;
				dBn = a[3];
				dBp = a[2];
				dBnQ = 0.0, dBpQ = 0.0;
				cnst = a[4];
				pIndex[1] = 1;
				pIndex[2] = 3;
				pIndex[3] = 2;
				pIndex[4] = 0;
				pIndex[5] = 0;
				pIndex[6] = 4;
				fprintf(stderr, " --- %f %f %f %f %f %f %f\n", Bn, Bp, dBn, dBp, dBnQ, dBpQ, cnst);
			}
		}
		else if (tiePoints->deltaB == DELTABCONST)
		{
			if (npts < 2)
				fewPoints();
			svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &BpCorrectOnlyParamsCoeffs);
			Bn = 0;
			Bp = a[1];
			dBn = 0;
			dBp = 0.0;
			dBnQ = 0., dBpQ = 0.;
			cnst = tiePoints->cnstR;
			pIndex[1] = 0;
			pIndex[2] = 0;
			pIndex[3] = 0;
			pIndex[4] = 0;
			pIndex[5] = 0;
			pIndex[6] = 1;
			fprintf(stderr, "fit 1 %f %i\n", a[1], ma);
		}
		else if (tiePoints->deltaB == DELTABQUAD)
		{
			if (npts < 6)
				fewPoints();
			svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, &rParamsCoeffsQuadCorrect);
			cnst = 0.0;
			Bn = a[1];
			Bp = a[4];
			dBn = a[3];
			dBp = a[2];
			dBnQ = a[6], dBpQ = a[5];
			pIndex[1] = 1;
			pIndex[2] = 3;
			pIndex[3] = 2;
			pIndex[4] = 6;
			pIndex[5] = 5;
			pIndex[6] = 4;
		}
		else
			error("computeRParams: invalid deltaB flag ", tiePoints->deltaB);
		svdvar(v, ma, w, Cp);
		fprintf(stderr, " --- %f %f %f %f %f %f %f\n", Bn, Bp, dBn, dBp, dBnQ, dBpQ, cnst);
	}
	/* 12/17/21: Remove chisq since sigP is direct estimate of the variance - note left text X2/n in case other prgrams expect it */
	fprintf(stdout, ";* sigma*sqrt(X2/n)= %lf \n", sigP);
	fprintf(stdout, "; Pseudo erors \n;");
	if (tiePoints->deltaB == DELTABNONE)
		fprintf(stdout, "; Covariance Matrix  Bn, dBn,dBp, dBnQ, dBpQ, const\n;");
	else
		fprintf(stdout, "; Covariance Matrix  Bn, dBn,dBp, dBnQ, dBpQ, Bp\n;");

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
	if (tiePoints->deltaB == DELTABNONE)
	{
		fprintf(stdout, "; Bn Bp dBn dBp const dBnQ dBpQ \n");
		if (tiePoints->bnbpdBpFlag == TRUE)
		{
			fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f %f\n&\n", a[1], tiePoints->BpCorig, tiePoints->dBnorig, a[2], a[3], 0.0, 0.0);
		}
		else if (tiePoints->bpdBpFlag == TRUE)
		{
			fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f %f\n&\n", tiePoints->BnCorig, tiePoints->BpCorig, tiePoints->dBnorig,
					a[1], 0.0, 0.0, a[2]);
		}
		else if (tiePoints->constOnlyFlag == TRUE)
		{
			fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f %f\n&\n", tiePoints->BnCorig, tiePoints->BpCorig, tiePoints->dBnorig,
					tiePoints->dBporig, a[1], 0.0, 0.0);
		}
		else if (tiePoints->quadB == TRUE)
		{
			fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f %f\n&\n", a[1], tiePoints->BpCorig, a[3], a[2], a[4], a[6], a[5]);
		}
		else
		{
			fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f %f\n&\n", a[1], tiePoints->BpCorig, a[3], a[2], a[4], 0.0, 0.0);
		}
	}
	else
	{
		fprintf(stdout, "; Corrections to Bn Bp dBn dBp const dBnQ dBpQ  \n");
		if (tiePoints->deltaB == DELTABCONST)
		{
			fprintf(stdout, "%2.1f  %11.7f  %2.1f  %2.1f %2.1f %2.1f %2.1f\n&\n", 0., a[1], 0., 0., 0., 0., 0.);
		}
		else if (tiePoints->deltaB == DELTABQUAD)
		{
			fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f %f\n&\n", a[1], a[4], a[3], a[2], 0.0, a[6], a[5]);
		}
	}
	return;
}

void fewPoints()
{
	fprintf(stdout, "0. 0. 0. 0. 0. 0.\n&\n");
	exit(-1);
}

void rbpdBpParamsCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	/*
	  This version is for computing dBp
	*/
	double thetaD;

	modelValues xx, *xy;
	double xAz;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;
	thetaD = xx.thetaD;

	afunc[1] = -xAz * cos(thetaD);
	afunc[2] = 1.0;
	return;
}

void constOnlyParamsCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	/*
	  This version is for computing dBp
	*/
	double thetaD;
	modelValues xx, *xy;
	double xAz;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;
	thetaD = xx.thetaD;
	afunc[1] = 1.0;
	return;
}

void BpCorrectOnlyParamsCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	/*
	  This version is for computing dBp
	*/
	double thetaD;
	modelValues xx, *xy;
	double xAz;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;
	thetaD = xx.thetaD;
	afunc[1] = -cos(thetaD);
	return;
}

void rParamsCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	double thetaD;
	modelValues xx, *xy;
	double xAz;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;
	thetaD = xx.thetaD;

	afunc[1] = -sin(thetaD);
	afunc[2] = -xAz * cos(thetaD);
	afunc[3] = -xAz * sin(thetaD);
	afunc[4] = 1.0;
	return;
}

void rParamsCoeffsQuad(void *x, int32_t i, double *afunc, int32_t ma)
{
	double thetaD;
	modelValues xx, *xy;
	double xAz;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;
	thetaD = xx.thetaD;

	afunc[1] = -sin(thetaD);
	afunc[2] = -xAz * cos(thetaD);
	afunc[3] = -xAz * sin(thetaD);
	afunc[4] = 1.0;
	afunc[5] = -xAz * xAz * cos(thetaD);
	afunc[6] = -xAz * xAz * sin(thetaD);
	return;
}

void rParamsCoeffsQuadCorrect(void *x, int32_t i, double *afunc, int32_t ma)
{
	double thetaD;
	modelValues xx, *xy;
	double xAz;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;
	thetaD = xx.thetaD;

	afunc[1] = -sin(thetaD);
	afunc[2] = -xAz * cos(thetaD);
	afunc[3] = -xAz * sin(thetaD);
	afunc[4] = -cos(thetaD);
	afunc[5] = -xAz * xAz * cos(thetaD);
	afunc[6] = -xAz * xAz * sin(thetaD);
	return;
}

void rbnbpdBpParamsCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	double thetaD;
	modelValues xx, *xy;
	double xAz;

	xy = (modelValues *)x;
	xx = xy[i];
	xAz = xx.x;
	thetaD = xx.thetaD;

	afunc[1] = -sin(thetaD);
	afunc[2] = -xAz * cos(thetaD);
	/*     afunc[3] = -xAz * sin(thetaD);*/
	afunc[3] = 1.0;
	return;
}
