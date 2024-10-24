#include <math.h>
#include "mosaicSource/common/common.h"
#include "azparams.h"
#include "cRecipes/nrutil.h"
#include <stdlib.h>
#define AZCONST 800000.
/*
  Estimate baseline parameters.
*/

static void constOnlyCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
static void azCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
static void azCoeffsLinear(void *x, int32_t i, double *afunc, int32_t ma);
static void linearOnlyCoeffs(void *x, int32_t i, double *afunc, int32_t ma);
static void getBaselineRates(double *dbcds, double *dbhds, char *baseFile, double prf, double slPixSize);

static void constantOnlyFit(void *x, double y[], double sig[], int32_t npts, double a[], int32_t ma, double **u, double **v, double w[],
							double *chisq, double dbcds, double dbhds, double result[], int32_t pIndex[], double azconst[]);

static void linearOnlyFit(void *x, double y[], double sig[], int32_t npts, double a[], int32_t ma, double **u, double **v, double w[],
						  double *chisq, double dbcds, double dbhds, double result[], int32_t pIndex[], double azconst[]);

static void azParamsConstFit(void *x, double y[], double sig[], int32_t npts, double a[], int32_t ma, double **u, double **v,
							 double w[], double *chisq, double dbhds, double result[], int32_t pIndex[], double azconst[]);

static void azParamsLinearFit(void *x, double y[], double sig[], int32_t npts, double a[], int32_t ma, double **u, double **v, double w[],
							  double *chisq, double dbhds, double result[], int32_t pIndex[], double azconst[]);

static void computeBaselineRates(inputImageStructure *inputImage, Offsets *offsets, double thetaC, double *dbc, double *dbh)
{
	double bTCNLate[3], bTCNEarly[3];
	conversionDataStructure *cP;
	double azTime1, azTime2;
	double bn1, bn2, bp1, bp2;
	cP = &(inputImage->cpAll);
	// Compute baseline at start and end for image
	azTime1 = cP->sTime + 0 * inputImage->nAzimuthLooks / inputImage->par.prf;
	azTime2 = cP->sTime + inputImage->azimuthSize * inputImage->nAzimuthLooks / inputImage->par.prf;
	svBaseTCN(azTime1, offsets->dt1t2, &(inputImage->sv), &(offsets->sv2), bTCNEarly);
	svBaseTCN(azTime2, offsets->dt1t2, &(inputImage->sv), &(offsets->sv2), bTCNLate);
	// m (prf * slPixSize)
	*dbc = (bTCNLate[1] - bTCNEarly[1]) / (inputImage->azimuthSize * inputImage->nAzimuthLooks * inputImage->par.slpA);
	*dbh = (bTCNLate[2] - bTCNEarly[2]) / (inputImage->azimuthSize * inputImage->nAzimuthLooks * inputImage->par.slpA);
	fprintf(stderr, "--- %e %e \n", *dbc, *dbh);
}

void computeAzParams(tiePointsStructure *tiePoints, inputImageStructure *inputImage, char *baseFile, Offsets *offsets)
{
	double Re, H, RNear, dr;
	double *a; /* Solution for params */
	double theta, theta1, z, r0;
	double thetaC, bn, bp, dbn, dbp;
	double **v, **u, chisq;
	double *y, *sig, *w, *chsq, *sigB, **Cp;
	double azimuth, ReH, Cij;
	aZmodelValues *x;
	double dbhds, dbcds;
	int32_t i, i1, k, j, npts, l1, l2;
	double result[5], azconst[5]; /* fit result, constant used to condition solution */
	double varP, sigP, meanP, xtmp, weightSum;
	conversionDataStructure *cP;
	int32_t nParams;
	int32_t pIndex[5];
	int32_t nData, ma;
	double svSol;

	if (tiePoints->constOnlyFlag == TRUE)
	{
		if (tiePoints->linFlag == FALSE)
			ma = 1;
		else
			ma = 2;
	}
	else
	{
		if (tiePoints->linFlag == FALSE)
			ma = 2;
		else
			ma = 3;
	}
	fprintf(stderr, "ma = %i\n", ma);
	a = dvector(1, ma);
	initllToImageNew(inputImage);
	cP = &(inputImage->cpAll);
	Re = cP->Re;
	H = inputImage->par.H;
	ReH = getReH(cP, inputImage, azimuth);
	RNear = cP->RNear;
	fprintf(stderr, "---------------------RNear %f %f \n", RNear, H);

	/* Fixed 12/1/22 to force zero for sv solutions */
	thetaC = thetaRReZReH(cP->RCenter, (Re + 0), (ReH));
	if (tiePoints->initWithSV == TRUE)
	{
		computeBaselineRates(inputImage, offsets, thetaC, &dbcds, &dbhds);
	}
	else
	{
		getBaselineRates(&dbcds, &dbhds, baseFile, inputImage->par.prf, inputImage->par.slpA);
	}
	fprintf(stderr, "+++ %e %e\n", dbcds, dbhds);
	if (tiePoints->deltaB != DELTABNONE)
	{
		dbcds = 0.0;
		dbhds = 0.0;
	}
	/*
	  Range comp params
	*/
	dr = inputImage->rangePixelSize;
	/*
	  Init arrays for least squares fit.
	*/
	npts = 0;
	nData = tiePoints->npts;
	for (i = 0; i < nData; i++)
		if (fabs(tiePoints->phase[i]) < 1.0E6)
			npts++;
	x = (aZmodelValues *)malloc((npts + 1) * sizeof(aZmodelValues));
	y = dvector(1, npts);
	sig = dvector(1, npts);
	u = dmatrix(1, npts, 1, ma);
	v = dmatrix(1, ma, 1, ma);
	w = dvector(1, ma);
	sigB = dvector(1, ma);
	Cp = dmatrix(1, ma, 1, ma);
	sigP = 10.0;
	/* Compute mean weight for renormalization */
	weightSum = 0.0;
	for (i = 0; i < nData; i++)
		if (fabs(tiePoints->phase[i]) < 1.0E6)
			weightSum += tiePoints->weight[i];
	/*
	  Run 3 times 1) initial estimate with unknown errors, 2) use estimate to determine residual 3) final solution with sigma detemermine by residual
	 */
	result[1] = 0;
	result[2] = dbcds;
	result[3] = dbhds;
	result[4] = 0.0;
	for (k = 0; k < 3; k++)
	{
		j = 0;
		varP = 0;
		meanP = 0;
		/*
		Prep data for solution
		*/
		fprintf(stderr, "result %e %e %e %e\n", result[1], result[2], result[3], result[4]);
		for (i = 0; i < nData; i++)
		{
			if (fabs(tiePoints->phase[i]) < 1.0E6)
			{ /* Use only good points */
				i1 = j + 1;
				z = tiePoints->z[i];
				r0 = RNear + tiePoints->r[i] * dr;
				azimuth = tiePoints->a[i];
				ReH = getReH(cP, inputImage, azimuth);
				theta = thetaRReZReH(r0, (Re + z), ReH);
				x[i1].r = r0;
				x[i1].theta = theta;
				x[i1].x = tiePoints->x[i] / (inputImage->azimuthSize * inputImage->nAzimuthLooks * inputImage->par.slpA);
				/* never solve for dbhds, so remove - will be zero for sv case */
				y[i1] = tiePoints->phase[i] - (-r0 * cos(theta) * dbhds);
				/* If solving as correction to sv model */
				xtmp = result[1] + r0 * sin(theta) * result[2] + result[4] * x[i1].x;
				if (tiePoints->deltaB != DELTABNONE)
				{
					svSol = svAzOffset(inputImage, offsets, tiePoints->r[i], tiePoints->a[i]);
					y[i1] -= svSol;
				}
				/* Const only means don't solve for dbcds - const or linear along track ok */
				if (tiePoints->constOnlyFlag == TRUE)
				{
					y[i1] -= r0 * sin(theta) * dbcds;
					xtmp -= r0 * sin(theta) * dbcds;
				}
				/* Compute mean and variance between model (xtmp) and corrected data (y) */
				varP += (y[i1] - xtmp) * (y[i1] - xtmp);
				meanP += (y[i1] - xtmp);
				/* Updated 12/17/21:
				Add capability to weight sigmas to emphasize certain points (e.g., rock w=1) and de-emphasize others (ice that could change w=10).
				The noise may be the same for all points, but the weights help limit points where we think there could be subtle biases. To try and keep
				the mean noise level to that of the data, the normalization will keep the mean sigma consistent with the data. This is not much of a
				problem if most of the points have the same weight. But in cases where say there are 10% low weight points and 90% high weight, this will drive the
				sigmas on the low weight (good) to unreasonably low errors. Either way, it shouldn't affect the solution, but it could skew the covariance matrix
				and hence the baseline error estimates. To avoid this situation, sigmas cannot drop below 0.1.
				*/
				sig[i1] = max(sigP * tiePoints->weight[i] * npts / weightSum, 0.1 * sigP);
				/* sig[i1] = sigP; */
				j++;
			} /* End if */
		} /* End for i */
		if (i1 < ma)
			error("aparams: Insufficient Number (%i)  of Valid tie points \n", i1);
		/*
		  Solve for parameters
		*/
		varP = varP / npts;
		meanP = meanP / npts;
		sigP = sqrt(varP - meanP * meanP);
		fprintf(stderr, "mean sigma %lf %lf \n", meanP, sigP);
		if (tiePoints->constOnlyFlag == TRUE)
		{
			if (tiePoints->linFlag == FALSE)
			{
				constantOnlyFit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, dbcds, dbhds, result, pIndex, azconst);
			}
			else
			{
				linearOnlyFit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, dbcds, dbhds, result, pIndex, azconst);
			}
		}
		else
		{
			if (tiePoints->linFlag == FALSE)
			{
				azParamsConstFit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, dbhds, result, pIndex, azconst);
			}
			else
			{
				azParamsLinearFit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, dbhds, result, pIndex, azconst);
			}
		}
		/*
		  Compute covariance matrix
		 */
		svdvar(v, ma, w, Cp);
	} /* end k */
	/*
	  Ouput results
	 */
	fprintf(stdout, ";\n; Ntiepoints/Ngiven used= %i/%i\n;\n", npts, nData);
	fprintf(stdout, "; X2 %f \n", chisq);
	/* 12/17/21: Remove chisq since sigP is direct estimate of the variance - note left text X2/n in case other prgrams expect it */
	fprintf(stdout, ";*  sigma*sqrt(X2/n)= %f \n", sigP); /**sqrt(chisq/(double)npts));*/
	fprintf(stdout, "; Covariance Matrix \n;");

	for (l1 = 1; l1 <= 4; l1++)
	{
		fprintf(stdout, ";* C_%1i ", l1);
		for (l2 = 1; l2 <= 4; l2++)
		{
			Cij = 0;
			if (pIndex[l1] > 0 && pIndex[l1] <= ma && pIndex[l2] > 0 && pIndex[l2] <= ma)
				Cij = Cp[pIndex[l1]][pIndex[l2]] * azconst[l1] * azconst[l2];
			fprintf(stdout, " %10.6e ", Cij);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, ";\n");
	fprintf(stdout, "%lf %le %le %lf\n", result[1], result[2], result[3], result[4]);
	return;
}

static void constantOnlyFit(void *x, double y[], double sig[], int32_t npts, double a[], int32_t ma,
							double **u, double **v, double w[], double *chisq, double dbcds, double dbhds, double result[],
							int32_t pIndex[], double azconst[])
{
	svdfit((void *)x, y, sig, npts, a, ma, u, v, w, chisq, constOnlyCoeffs); /* use sigma=10 as in above */
	result[1] = a[1];
	result[2] = dbcds;
	result[3] = dbhds;
	result[4] = 0.0;
	pIndex[1] = 1;
	pIndex[2] = 0;
	pIndex[3] = 0;
	pIndex[4] = 0;
	azconst[1] = 1.0;
	azconst[2] = 1.0;
	azconst[3] = 1.0;
	azconst[4] = 1.0;
}

static void linearOnlyFit(void *x, double y[], double sig[], int32_t npts, double a[], int32_t ma, double **u, double **v,
						  double w[], double *chisq, double dbcds, double dbhds, double result[], int32_t pIndex[], double azconst[])
{
	svdfit((void *)x, y, sig, npts, a, ma, u, v, w, chisq, &linearOnlyCoeffs);
	result[1] = a[1];
	result[2] = dbcds;
	result[3] = dbhds;
	result[4] = a[2];
	pIndex[1] = 1;
	pIndex[2] = 0;
	pIndex[3] = 0;
	pIndex[4] = 2;
	azconst[1] = 1.0;
	azconst[2] = 1.0;
	azconst[3] = 1.0;
	azconst[4] = 1.0;
}

static void azParamsConstFit(void *x, double y[], double sig[], int32_t npts, double a[], int32_t ma, double **u,
							 double **v, double w[], double *chisq, double dbhds, double result[], int32_t pIndex[], double azconst[])
{
	svdfit((void *)x, y, sig, npts, a, ma, u, v, w, chisq, &azCoeffs);
	result[1] = a[1] * AZCONST;
	result[2] = a[2];
	result[3] = dbhds;
	result[4] = 0.0;
	pIndex[1] = 1;
	pIndex[2] = 2;
	pIndex[3] = 0;
	pIndex[4] = 0;
	azconst[1] = AZCONST;
	azconst[2] = 1.0;
	azconst[3] = 1.0;
	azconst[4] = 1.0;
}

static void azParamsLinearFit(void *x, double y[], double sig[], int32_t npts, double a[], int32_t ma, double **u,
							  double **v, double w[], double *chisq,
							  double dbhds, double result[], int32_t pIndex[], double azconst[])
{
	svdfit((void *)x, y, sig, npts, a, ma, u, v, w, chisq, &azCoeffsLinear);
	result[1] = a[1] * AZCONST;
	result[2] = a[2];
	result[3] = dbhds;
	result[4] = a[3] * AZCONST;
	pIndex[1] = 1;
	pIndex[2] = 2;
	pIndex[3] = 0;
	pIndex[4] = 3;
	azconst[1] = AZCONST;
	azconst[2] = 1.0;
	azconst[3] = 1.0;
	azconst[4] = AZCONST;
}

void constOnlyCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	afunc[1] = 1;
	return;
}

void linearOnlyCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	aZmodelValues xx, *xy;

	xy = (aZmodelValues *)x;
	xx = xy[i];
	afunc[1] = 1;
	afunc[2] = xx.x;
	return;
}

void azCoeffs(void *x, int32_t i, double *afunc, int32_t ma)
{
	double theta, r;
	aZmodelValues xx, *xy;

	xy = (aZmodelValues *)x;
	xx = xy[i];

	theta = xx.theta;
	r = xx.r;

	afunc[1] = AZCONST;
	afunc[2] = r * sin(theta);
	return;
}

/* Not used */
void azCoeffsLinear(void *x, int32_t i, double *afunc, int32_t ma)
{
	double theta, r;
	aZmodelValues xx, *xy;

	xy = (aZmodelValues *)x;
	xx = xy[i];

	theta = xx.theta;
	r = xx.r;
	afunc[1] = AZCONST;
	afunc[2] = r * sin(theta);
	afunc[3] = xx.x * AZCONST;
	return;
}

static void getBaselineRates(double *dbcds, double *dbhds, char *baseFile, double prf, double slPixSize)
{
	char line[1024];
	int32_t lineCount, eod;
	double x1, x2, x3;
	FILE *fp;

	fprintf(stderr, "prf,nSingleLook %f %f %s\n", prf, slPixSize, baseFile);
	fp = openInputFile(baseFile);
	lineCount = getDataString(fp, lineCount, line, &eod);
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%lf%lf%lf", &x1, &x2, &x3) != 3)
		error("%s  %i of %s", "readOffsets -- Missing image parameters at line:", lineCount, baseFile);
	*dbcds = x2 / (prf * slPixSize);
	*dbhds = x3 / (prf * slPixSize);
	fprintf(stderr, "dbcds,dbhds, %f %f %f %f\n", *dbcds, *dbhds, x2, x3);
	fclose(fp);
}
