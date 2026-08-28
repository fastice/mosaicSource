#include <math.h>
#include "mosaicSource/common/common.h"
#include "mosaicSource/common/writeTieResidualsGpkg.h"
#include "tiePoints.h"
#include "cRecipes/nrutil.h"
#include <stdlib.h>
#include <string.h>
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

/* Result of fitting one baseline solution (unsquinted or squinted phase array) */
typedef struct
{
	int32_t npts, nData;
	int32_t insufficientPoints; /* TRUE if fewer than 4 valid tie points -- no solution */
	double chisq, sigma;
	double Bn, Bp, dBn, dBp, dBnQ, dBpQ;
	double C[7][7];
	tieResidualsType residuals; /* only populated when fitBaseline() is called with debugFlag */
} BaselineFit;

/*
  Fit one baseline solution against phaseArr (parallel to tiePoints->npts -- either
  tiePoints->phase [unsquinted] or tiePoints->phaseSquint [squinted]). Pulled out of
  computeBaseline() so it can be run twice -- see computeBaseline() below.
*/
static BaselineFit fitBaseline(tiePointsStructure *tiePoints, double *phaseArr,
							   inputImageStructure inputImage, int32_t verbose, int32_t debugFlag)
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
	BaselineFit fit;
	int32_t *origIndex = NULL;
	double *residual = NULL;
	void (*coeffsFn)(void *, int32_t, double *, int32_t) = NULL;

	memset(&fit, 0, sizeof(fit));
	fit.insufficientPoints = FALSE;
	lambda = inputImage.par.lambda;
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
	if (verbose) fprintf(stderr, "------------------+++ RNear %f %f %f %f\n", RNear, thetaC * RTOD, tiePoints->thetaC * RTOD, ReH);
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
		if (fabs(phaseArr[i]) < 1.0E6)
			npts++;
	fprintf(stderr, "%i points wavelength %f\n", npts, lambda);
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
	if (debugFlag)
		origIndex = (int32_t *)malloc(npts * sizeof(int32_t));
	/*
	  Loop twice, first using flattening value of bsq, and then value
	  from first fit. Should easily converge with just two iterations
	  unless very large error in initial estimate.
	  Changed to use base estimate only. Leaving loop for possible
	  later modification.

	  added second loop on 4/28/14 to iterate on delta and bsq
	*/
	sigP = 10.0; /* Use for first try */
	fprintf(stderr, "iter: Bn Bp dBn dBp dBnQ dBpQ meanPhase sigPhase\n");
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
			if (fabs(phaseArr[i]) < 1.0E6)
			{ /* Use only good points */
				i1 = j + 1;
				if (debugFlag)
					origIndex[i1 - 1] = i;
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
				y[i1] = phaseArr[i] - twok * bSq / (2.0 * r0) + pow(deltaApprox, 2.0) * twok / (2.0 * r0);
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
		if (j < 4) {
			/* j (reliably 0-initialized above, incremented per valid point) not i1 --
			   i1 is only ever assigned inside the valid-point branch, so it's
			   uninitialized stack garbage whenever zero points pass, which previously
			   made both this message and the original "if (i1 < 4)" check itself read
			   garbage in that exact case (could in principle even skip this check if
			   the garbage happened to be >= 4). */
			fprintf(stderr, "error:  tiepoints: Insufficient Number (%i)  of Valid tie points \n", j);
			/* No solution -- let computeBaseline() emit the sigma<0 sentinel for both
			   blocks and exit(1) once, rather than doing it per-fit here. */
			fit.insufficientPoints = TRUE;
			fit.npts = j;
			fit.nData = nData;
			return fit;
		}

		if (tiePoints->dBpFlag == TRUE)
		{
			if (tiePoints->quadB == TRUE)
				coeffsFn = &dBpQBaselineCoeffs;
			else if (tiePoints->bpFlag == TRUE)
				coeffsFn = &bpBaselineCoeffs;
			else if (tiePoints->bnbpFlag == TRUE)
				coeffsFn = &bpbnBaselineCoeffs;
			else if (tiePoints->bpdBpFlag == TRUE)
				coeffsFn = &bpdBpBaselineCoeffs;
			else if (tiePoints->bnbpdBpFlag == TRUE)
				coeffsFn = &bnbpdBpBaselineCoeffs;
			else
				coeffsFn = &dBpBaselineCoeffs;
		}
		else
		{
			if (tiePoints->bpFlag == TRUE)
				coeffsFn = &bpBaselineCoeffs;
			else if (tiePoints->bnbpFlag == TRUE)
				coeffsFn = &bpbnBaselineCoeffs;
			else if (tiePoints->noRamp == TRUE)
				coeffsFn = &noRampBaselineCoeffs;
			else
				coeffsFn = &baselineCoeffs;
		}
		svdfit((void *)x, y, sig, npts, a, ma, u, v, w, &chisq, coeffsFn);
		if (debugFlag && k == 2)
		{
			double afunc[7];
			residual = (double *)malloc(npts * sizeof(double));
			for (i1 = 1; i1 <= npts; i1++)
			{
				(*coeffsFn)((void *)x, i1, afunc, ma);
				double model = 0.0;
				for (int32_t m = 1; m <= ma; m++)
					model += afunc[m] * a[m];
				residual[i1 - 1] = y[i1] - model;
			}
		}
		svdvar(v, ma, w, Cp);
		if (tiePoints->dBpFlag == FALSE && !tiePoints->bpFlag && !tiePoints->bnbpFlag)
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

	/* Extract fitted baseline values for output */
	double BnOut, BpOut, dBnOut, dBpOut, dBnQOut, dBpQOut;
	if (tiePoints->quadB == TRUE) {
		BnOut=a[1]; BpOut=a[2]; dBnOut=a[4]; dBpOut=a[3]; dBnQOut=a[6]; dBpQOut=a[5];
	} else if (tiePoints->bpFlag == TRUE) {
		BnOut=tiePoints->BnCorig; BpOut=a[1]; dBnOut=tiePoints->dBnorig; dBpOut=tiePoints->dBporig;
		dBnQOut=tiePoints->dBnQorig; dBpQOut=tiePoints->dBpQorig;
	} else if (tiePoints->bnbpFlag == TRUE) {
		BnOut=a[1]; BpOut=a[2]; dBnOut=tiePoints->dBnorig; dBpOut=tiePoints->dBporig;
		dBnQOut=tiePoints->dBnQorig; dBpQOut=tiePoints->dBpQorig;
	} else if (tiePoints->bpdBpFlag == TRUE) {
		BnOut=tiePoints->BnCorig; BpOut=a[1]; dBnOut=tiePoints->dBnorig; dBpOut=a[2];
		dBnQOut=tiePoints->dBnQorig; dBpQOut=tiePoints->dBpQorig;
	} else if (tiePoints->bnbpdBpFlag == TRUE) {
		BnOut=a[1]; BpOut=a[2]; dBnOut=tiePoints->dBnorig; dBpOut=a[3];
		dBnQOut=tiePoints->dBnQorig; dBpQOut=tiePoints->dBpQorig;
	} else {
		BnOut=a[1]; BpOut=a[2]; dBnOut=a[4]; dBpOut=a[3]; dBnQOut=0.0; dBpQOut=0.0;
	}

	fit.npts = npts;
	fit.nData = nData;
	fit.chisq = chisq;
	fit.sigma = sigP * sqrt(chisq / (double)npts);
	fit.Bn = BnOut;
	fit.Bp = BpOut;
	fit.dBn = dBnOut;
	fit.dBp = dBpOut;
	fit.dBnQ = dBnQOut;
	fit.dBpQ = dBpQOut;
	for (l1 = 1; l1 <= 6; l1++)
		for (l2 = 1; l2 <= 6; l2++) {
			Cij = 0;
			if (pIndex[l1] > 0 && pIndex[l1] <= ma && pIndex[l2] > 0 && pIndex[l2] <= ma)
				Cij = Cp[pIndex[l1]][pIndex[l2]];
			fit.C[l1][l2] = Cij;
		}

	if (debugFlag)
	{
		fit.residuals.npts = npts;
		fit.residuals.origIndex = origIndex;
		fit.residuals.residual = residual;
	}

	return fit;
}

/*
  Write the legacy (non-YAML) text baseline solution -- always the unsquinted fit, matching
  prior behavior (the legacy format has no notion of squint).
*/
static void writeLegacyTextSolution(tiePointsStructure *tiePoints, BaselineFit *fit)
{
	int32_t l1, l2;

	fprintf(stdout, ";\n; Ntiepoints/Ngiven used= %i/%i\n;\n", fit->npts, fit->nData);
	fprintf(stdout, "; X2 %f \n", fit->chisq);
	fprintf(stdout, ";*  sigma*sqrt(X2/n)= %f \n", fit->sigma);
	fprintf(stdout, "; Covariance Matrix \n;");
	for (l1 = 1; l1 <= 6; l1++)
	{
		fprintf(stdout, ";* C_%1i ", l1);
		for (l2 = 1; l2 <= 6; l2++)
			fprintf(stdout, " %10.6e ", fit->C[l1][l2]);
		fprintf(stdout, "\n");
	}
	if (tiePoints->dBpFlag == TRUE)
		fprintf(stdout, ";\n; Estimated Baseline\n; Bn,Bp,dBn,dBp\n;\n");
	else
		fprintf(stdout, ";\n; Estimated Baseline\n; Bn,Bp,dBn,omegaA\n;\n");
	if (tiePoints->quadB == TRUE)
		fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f\n&\n", fit->Bn, fit->Bp, fit->dBn, fit->dBp, fit->dBnQ, fit->dBpQ);
	else if (tiePoints->dBpFlag == TRUE || tiePoints->bpFlag == TRUE || tiePoints->bnbpFlag == TRUE ||
			 tiePoints->bpdBpFlag == TRUE || tiePoints->bnbpdBpFlag == TRUE)
		fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f %f %f\n&\n", fit->Bn, fit->Bp, fit->dBn, fit->dBp, fit->dBnQ, fit->dBpQ);
	else
		fprintf(stdout, "%11.5f  %11.5f  %11.5f  %f\n&\n", fit->Bn, fit->Bp, fit->dBn, fit->dBp);
}

/*
  Write one labeled solution block (noSquint: or squint:) to stdout in YAML.
*/
static void writeYamlSolutionBlock(BaselineFit *fit)
{
	int32_t l1, l2;
	fprintf(stdout, "  nTiepoints: %d\n",               fit->npts);
	fprintf(stdout, "  nTiepointsGiven: %d\n",          fit->nData);
	fprintf(stdout, "  X2: %f\n",                       fit->chisq);
	fprintf(stdout, "  sigma: %.6f  # radians\n",       fit->sigma);
	fprintf(stdout, "  Bn: %.6f  # meters\n",           fit->Bn);
	fprintf(stdout, "  Bp: %.6f  # meters\n",           fit->Bp);
	fprintf(stdout, "  dBn: %.6f  # meters/pixel\n",    fit->dBn);
	fprintf(stdout, "  dBp: %.6f  # meters/pixel\n",    fit->dBp);
	fprintf(stdout, "  dBnQ: %.6f  # meters/pixel^2\n", fit->dBnQ);
	fprintf(stdout, "  dBpQ: %.6f  # meters/pixel^2\n", fit->dBpQ);
	fprintf(stdout, "  C:\n");
	for (l1 = 1; l1 <= 6; l1++) {
		fprintf(stdout, "    - [");
		for (l2 = 1; l2 <= 6; l2++)
			fprintf(stdout, "%s%10.6e", l2 > 1 ? ", " : "", fit->C[l1][l2]);
		fprintf(stdout, "]\n");
	}
}

/*
  Write the top-level ionosphere keys. Always emitted so the output unambiguously records
  the decision (matching rparams, which always writes offsetCorrectionFile: nil when it
  did not use one). The two sigmas only appear when both fits actually ran.
*/
static void writeYamlIonosphereKeys(int32_t haveIonFit, int32_t useIon, char *ionosphereFile,
									double sigmaIon, double sigmaNoIon)
{
	if (useIon == TRUE)
	{
		fprintf(stdout, "ionosphereCorrectionFile: %s\n", ionosphereFile);
	}
	else
	{
		fprintf(stdout, "ionosphereCorrectionFile: nil\n");
	}
	fprintf(stdout, "usingIon: %s\n", useIon ? "True" : "False");
	if (haveIonFit == TRUE)
	{
		fprintf(stdout, "sigmaWithIonCorrection: %f\n", sigmaIon);
		fprintf(stdout, "sigmaWithoutIonCorrection: %f\n", sigmaNoIon);
	}
}

/*
  Write a sigma<0 "no solution" sentinel block -- see getBaseline.c, which treats
  sigma<0 as an unambiguous, machine-readable "no solution" signal for callers.
*/
static void writeYamlSentinelBlock(void)
{
	fprintf(stdout, "  nTiepoints: 0\n");
	fprintf(stdout, "  nTiepointsGiven: 0\n");
	fprintf(stdout, "  X2: -1\n");
	fprintf(stdout, "  sigma: -1  # no solution -- insufficient tie points\n");
	fprintf(stdout, "  Bn: 0.000000  # meters\n");
	fprintf(stdout, "  Bp: 0.000000  # meters\n");
	fprintf(stdout, "  dBn: 0.000000  # meters/pixel\n");
	fprintf(stdout, "  dBp: 0.000000  # meters/pixel\n");
	fprintf(stdout, "  dBnQ: 0.000000  # meters/pixel^2\n");
	fprintf(stdout, "  dBpQ: 0.000000  # meters/pixel^2\n");
	fprintf(stdout, "  C:\n");
	fprintf(stdout, "    - [0, 0, 0, 0, 0, 0]\n");
	fprintf(stdout, "    - [0, 0, 0, 0, 0, 0]\n");
	fprintf(stdout, "    - [0, 0, 0, 0, 0, 0]\n");
	fprintf(stdout, "    - [0, 0, 0, 0, 0, 0]\n");
	fprintf(stdout, "    - [0, 0, 0, 0, 0, 0]\n");
	fprintf(stdout, "    - [0, 0, 0, 0, 0, 0]\n");
}

/*
  Fit and emit the baseline solution(s). When tiePoints->hasSquintSolution (squint
  coefficients were available in the input geodat and -motion/-vr was used), fits both
  the unsquinted (tiePoints->phase) and squinted (tiePoints->phaseSquint) phase arrays
  and writes both as separate labeled YAML blocks; mosaic3d -useSquint then picks which
  block getBaseline() reads. Legacy (non-YAML) text output only ever reflects the
  unsquinted solution, matching prior behavior.
*/
/*
  Write the debug residuals gpkg for one or both (noSquint/squint) fits. Called only
  after both fits have succeeded (or the legacy-only noSquint fit, in non-yaml mode) --
  never for a fit whose insufficientPoints is TRUE, matching the sentinel/exit(1)
  guards already in place above for the real YAML/legacy output.
*/
static void writeDebugResiduals(char *debugFile, tiePointsStructure *tiePoints,
								BaselineFit *noSquintFit, BaselineFit *squintFit,
								int32_t hasSquintSolution)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	OGRDataSourceH ds = openTieResidualsGpkg(debugFile, FALSE);
	writeTieResidualsLayer(ds, hasSquintSolution ? "residuals_noSquint" : "residuals",
							tiePoints, &noSquintFit->residuals,
							"phase_residual_rad", "Baseline fit residual (radians)",
							Rotation, tiePoints->stdLat, HemiSphere);
	free(noSquintFit->residuals.origIndex);
	free(noSquintFit->residuals.residual);
	if (hasSquintSolution) {
		writeTieResidualsLayer(ds, "residuals_squint", tiePoints, &squintFit->residuals,
								"phase_residual_rad",
								"Baseline fit residual (radians, squint-corrected)",
								Rotation, tiePoints->stdLat, HemiSphere);
		free(squintFit->residuals.origIndex);
		free(squintFit->residuals.residual);
	}
	closeTieResidualsGpkg(ds);
}

void computeBaseline(tiePointsStructure *tiePoints,
					 inputImageStructure inputImage, int32_t yamlOutput, int32_t verbose,
					 char *debugFile, int32_t ionosphereMode, double ionSigmaMargin,
					 char *ionosphereFile)
{
	BaselineFit noSquintFit, squintFit, ionFit;
	int32_t debugFlag = (debugFile != NULL);
	int32_t haveIonFit = FALSE, useIon = FALSE;
	double sigmaIon = -1.0, sigmaNoIon = -1.0;
	double *phaseIon = NULL;
	int32_t i;

	noSquintFit = fitBaseline(tiePoints, tiePoints->phase, inputImage, verbose, debugFlag);
	/*
	  Ionosphere: fit a second time against the ionosphere-corrected phase and keep whichever
	  fit is better. Only the noSquint solution participates -- S1 (the only source of these
	  corrections today) carries no squint polynomial, so hasSquintSolution is always FALSE
	  here. If that ever changes, the ion choice would need to be made per squint block.
	*/
	if (ionosphereMode != ION_NONE && tiePoints->hasIonosphere == TRUE)
	{
		phaseIon = (double *)malloc(tiePoints->npts * sizeof(double));
		for (i = 0; i < tiePoints->npts; i++)
		{
			/* ISCE convention: phaseCorrected = phase - ion (topsApp applies its own screen
			   as a*exp(-1.0*J*b), b = topophase.ion). NOTE this is a SUBTRACT, unlike the
			   range-offset ionosphere correction in rParams/getROffsets.c and
			   Mosaic3d/make3DOffsets.c, which is a pre-negated correction that gets added. */
			if (fabs(tiePoints->phase[i]) >= 1.0E6 || tiePoints->ionPhase[i] < -0.98 * LARGEINT)
			{
				phaseIon[i] = (double)(-LARGEINT); /* invalid -- fitBaseline drops it */
			}
			else
			{
				phaseIon[i] = tiePoints->phase[i] - tiePoints->ionPhase[i];
			}
		}
		ionFit = fitBaseline(tiePoints, phaseIon, inputImage, verbose, debugFlag);
		haveIonFit = TRUE;
		sigmaIon = ionFit.insufficientPoints ? -1.0 : ionFit.sigma;
		sigmaNoIon = noSquintFit.insufficientPoints ? -1.0 : noSquintFit.sigma;
		if (ionFit.insufficientPoints == FALSE)
		{
			if (ionosphereMode == ION_FORCE || noSquintFit.insufficientPoints)
			{
				useIon = TRUE;
			}
			else
			{
				useIon = (ionFit.sigma < noSquintFit.sigma * (1.0 - ionSigmaMargin));
			}
		}
		fprintf(stderr, "sigma with ion correction: %f  without: %f  (margin %.3f) -- using %s\n",
				sigmaIon, sigmaNoIon, ionSigmaMargin, useIon ? "with ion" : "without ion");
		/* Keep the winner in noSquintFit so everything below is unchanged. */
		if (useIon == TRUE)
		{
			if (debugFlag && noSquintFit.insufficientPoints == FALSE)
			{
				free(noSquintFit.residuals.origIndex);
				free(noSquintFit.residuals.residual);
			}
			noSquintFit = ionFit;
		}
		else if (debugFlag && ionFit.insufficientPoints == FALSE)
		{
			free(ionFit.residuals.origIndex);
			free(ionFit.residuals.residual);
		}
		free(phaseIon);
	}

	if (!yamlOutput) {
		/* Legacy text format has no notion of squint -- unsquinted solution only,
		   matching prior behavior. */
		if (noSquintFit.insufficientPoints)
			exit(1);
		writeLegacyTextSolution(tiePoints, &noSquintFit);
		if (useIon == TRUE)
		{
			fprintf(stdout, ";* ionosphereCorrectionFile %s\n", ionosphereFile);
		}
		if (debugFlag)
			writeDebugResiduals(debugFile, tiePoints, &noSquintFit, NULL, FALSE);
		return;
	}

	if (tiePoints->hasSquintSolution)
		squintFit = fitBaseline(tiePoints, tiePoints->phaseSquint, inputImage, verbose, debugFlag);

	if (noSquintFit.insufficientPoints || (tiePoints->hasSquintSolution && squintFit.insufficientPoints)) {
		fprintf(stdout, "applyFlatEarth: true\n");
		fprintf(stdout, "hasSquintSolution: %s\n", tiePoints->hasSquintSolution ? "true" : "false");
		fprintf(stdout, "nDays: %.6f  # days\n", tiePoints->nDays);
		writeYamlIonosphereKeys(haveIonFit, useIon, ionosphereFile, sigmaIon, sigmaNoIon);
		fprintf(stdout, "noSquint:\n");
		writeYamlSentinelBlock();
		if (tiePoints->hasSquintSolution) {
			fprintf(stdout, "squint:\n");
			writeYamlSentinelBlock();
		}
		exit(1);
	}

	fprintf(stdout, "applyFlatEarth: true\n");
	fprintf(stdout, "hasSquintSolution: %s\n", tiePoints->hasSquintSolution ? "true" : "false");
	fprintf(stdout, "nDays: %.6f  # days\n", tiePoints->nDays);
	writeYamlIonosphereKeys(haveIonFit, useIon, ionosphereFile, sigmaIon, sigmaNoIon);
	fprintf(stdout, "noSquint:\n");
	writeYamlSolutionBlock(&noSquintFit);
	if (tiePoints->hasSquintSolution) {
		fprintf(stdout, "squint:\n");
		writeYamlSolutionBlock(&squintFit);
	}
	if (debugFlag)
		writeDebugResiduals(debugFile, tiePoints, &noSquintFit, &squintFit, tiePoints->hasSquintSolution);
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
