#include "stdio.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "sys/time.h"
#include "cRecipes/nrutil.h"
#include "mosaicSource/common/common.h"
#include "mosaic3d.h"

/*
************ Crossing-orbit velocity from RANGE OFFSETS -- JOINT (normal equations) ***********

  Drop-in alternative to make3DOffsets() (make3DOffsets.c), identical signature.  THE DEFAULT
  since 2026-08-29; "mosaic3d -legacyPairRange" selects the original.  Offsets analogue of
  make3DMosaicJoint.c.

  Read make3DMosaicJoint.c first: the derivation (per-image sensitivity row from computeA's
  N matrix, the exact (I-AB)^-1 A = (N-B)^-1 collapse, the monotone conditioning gate, and why
  a trace-normalized geometric test is wrong) is identical and is not repeated here.  Only the
  offsets-specific differences are documented below.

  Gated by -jointMaxSigmaRange (its own threshold, independent of the phase round).

  WHAT DIFFERS FROM THE PHASE VERSION
  -----------------------------------
  * The measurement is a range offset already in METRES (interpRangeOffsetInMeters), so the
    velocity scaling is 365.25/(nDays*sin(psi)) with no 4pi/lambda factor.

  * The error budget has four terms, not two:
        sigmaR^2 = interpRangeSigma^2 + demError^2 + sig2Base + rangeAccuracyVar
    interpRangeSigma is Cullst's LOCAL neighbourhood scatter (blind to long-wavelength error by
    construction); rangeAccuracyVar carries the rparams tie-point fit residual, which is the term
    that actually sees ionosphere and orbit ramps.  See mosaicSource/CLAUDE.md "Range-offset
    error budget".

  * IONOSPHERE SIGN IS OPPOSITE.  The offsets correction is a pre-negated correction that is
    ADDED (delta += ionCorrection * rSLPixSize), unlike the phase screen which is SUBTRACTED.
    Do not "fix" this to match the phase path -- see the root CLAUDE.md sign-convention note.

  * NO SQUINT.  make3DOffsets calls computeA(..., FALSE) unconditionally: the zero-Doppler
    condition forces true LOS perpendicular to true velocity at the assigned time regardless of
    squint, so range offsets are self-consistent by construction.  This routine likewise never
    applies the squint heading correction, flag or no flag.  Not an oversight -- see
    mosaicSource/CLAUDE.md "Squint".

  * SERIAL SVD PRE-INIT IS MANDATORY.  svInterpBnBp() lazily initialises global SVD workspace
    (uRows/uBuffer in common/svBase.c); calling it first from inside a parallel region races and
    segfaults.  It is primed serially per image below, exactly as make3DOffsets does.  Do not
    remove this even though it looks redundant.

  * Per-image gates make3DOffsets applies and this keeps: crossFlag, weight < 0.05, a NULL
    rFile, and the rparams sigmaRresidual < 0 "no solution" sentinel (which also latches
    rFile = NULL so later passes skip it).

  * timeThresh becomes a PER-IMAGE window against the mosaic centre, not a pair separation --
    the same reinterpretation the phase version makes for timeThreshPhase.

  KNOWN RISK, STATED UP FRONT
  ---------------------------
  Range offsets carry a heavier tail than phase: the worst 1% of crossing-offset vy pixels hold
  ~48% of the variance.  Reweighting (IRLS) was implemented and MEASURED TO DO NOTHING here --
  reduced chi-square is ~0.85, i.e. the measurements at a pixel already agree with each other to
  within their own sigmas, so there is no inconsistent minority to down-weight.  The tail is
  common-mode error at the pixel, not bad measurements.  Do not re-add reweighting; see
  Documents/crossingOrbitRedundancy.md §13.
*/

static float **mallocZeroImageOff(int32_t nr, int32_t nc);
static void freeImageOff(float **image, int32_t nr);

void make3DOffsetsJoint(inputImageStructure *allImages, vhParams *aParams, xyDEM *dem,
						outputImageStructure *outputImage, float fl, float timeThresh)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	extern int32_t indentRegionOutput;
	inputImageStructure *offImage;
	vhParams *params;
	conversionDataStructure *cp;
	ShelfMask *shelfMask;
	xyDEM *vCorrect;
	double Re, ReH, thetaC, ddum1, ddum2;
	double tCenter, tOffCenter;
	float **vXimage, **vYimage, **vZimage, **errorX, **errorY;
	float **vxTmp, **vyTmp, **vzTmp, **fScale, **sxTmp, **syTmp;
	float **scaleX, **scaleY, **scaleZ;
	/* Normal-equation accumulators; float is safe because cancellation in Nxx*Nyy - Nxy^2 only
	   bites as the solve approaches singular, and those pixels are rejected by the gate. */
	float **Nxx, **Nxy, **Nyy, **bxAcc, **byAcc;
	float **nObs;
	float **sumW, **sumWT;
	float **SddAcc; /* sum w d^2, for chi-square at solve time */
	float **chi2Plane;
	float rSLPixSize, dum;
	float azimuthMin, azimuthMax;
	int32_t iMin, iMax, jMin, jMax;
	int32_t iMinAll, iMaxAll, jMinAll, jMaxAll;
	int32_t ii, nTotal, nUsed;
	int32_t i;
	int64_t nSolved, nRejCond, nRejNoData;
	double nObsSum, nObsMax;
	struct timeval funcStart, funcEnd;
	gettimeofday(&funcStart, NULL);

	fprintf(outputImage->fpLog, ";\n; Entering make3DOffsetsJoint(.c)\n");
	shelfMask = outputImage->shelfMask;
	vCorrect = outputImage->verticalCorrection;
	if (dem->stdLat < 50 || dem->stdLat > 80)
	{
		error("make3DOffsetsJoint invalid slat for dem");
	}
	setupBuffers(outputImage, &vXimage, &vYimage, &vZimage, &scaleX, &scaleY, &scaleZ,
				 &vxTmp, &vyTmp, &vzTmp, &sxTmp, &syTmp, &fScale, &errorX, &errorY);
	computeScale((float **)vXimage, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0, (double)(-LARGEINT));
	undoNormalization(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, fScale, FALSE);

	Nxx = mallocZeroImageOff(outputImage->ySize, outputImage->xSize);
	Nxy = mallocZeroImageOff(outputImage->ySize, outputImage->xSize);
	Nyy = mallocZeroImageOff(outputImage->ySize, outputImage->xSize);
	bxAcc = mallocZeroImageOff(outputImage->ySize, outputImage->xSize);
	byAcc = mallocZeroImageOff(outputImage->ySize, outputImage->xSize);
	nObs = mallocZeroImageOff(outputImage->ySize, outputImage->xSize);
	SddAcc = mallocZeroImageOff(outputImage->ySize, outputImage->xSize);
	sumW = NULL;
	sumWT = NULL;
	chi2Plane = mallocZeroImageOff(outputImage->ySize, outputImage->xSize);
	if (outputImage->timeOverlapFlag == TRUE)
	{
		sumW = mallocZeroImageOff(outputImage->ySize, outputImage->xSize);
		sumWT = mallocZeroImageOff(outputImage->ySize, outputImage->xSize);
	}
	nTotal = 0;
	for (offImage = allImages; offImage != NULL; offImage = offImage->next)
	{
		nTotal++;
	}
	fprintf(stderr, "make3DOffsetsJoint: nTotal Images %i\n", nTotal);
	tCenter = (outputImage->jd1 + outputImage->jd2 + 1.) * 0.5;
	iMinAll = outputImage->ySize;
	jMinAll = outputImage->xSize;
	iMaxAll = 0;
	jMaxAll = 0;
	nUsed = 0;
	int nthreads = omp_get_max_threads();
	inputImageStructure *localImgs = (inputImageStructure *)malloc((size_t)nthreads * sizeof(inputImageStructure));
	if (localImgs == NULL)
	{
		error("make3DOffsetsJoint: malloc failed for per-thread image copies\n");
	}
	/*
	  ****************** PASS 1: accumulate normal equations, one image at a time **************
	*/
	ii = 0;
	params = aParams;
	for (offImage = allImages; offImage != NULL; offImage = offImage->next, params = params->next)
	{
		ii++;
		tOffCenter = offImage->julDay + params->nDays * 0.5;
		if (offImage->crossFlag == FALSE || offImage->weight < 0.05 || params->offsets.rFile == NULL)
		{
			continue;
		}
		/* Per-image temporal window.  NOT the pair gate make3DOffsets applies. */
		if (fabs(offImage->julDay - tCenter) > timeThresh)
		{
			fprintf(stderr, "\t%s skipped: |julDay - tCenter| = %.1f > %.1f\n",
					params->offsets.rFile, fabs(offImage->julDay - tCenter), timeThresh);
			continue;
		}
		indentRegionOutput = FALSE;
		if (!getRegion(offImage, &iMin, &iMax, &jMin, &jMax, outputImage))
		{
			continue;
		}
		if (iMin > iMax || jMin > jMax)
		{
			continue;
		}
		rSLPixSize = offImage->rangePixelSize / offImage->nRangeLooks;
		cp = setupGeoConversions(offImage, &dum, &rSLPixSize, &Re, &ReH, &thetaC, &ddum1, &ddum2);
		getAzimuthBoundsForXYBox(iMin, iMax, jMin, jMax, offImage, outputImage, &azimuthMin, &azimuthMax);
		if (azimuthMin == 0.0f && azimuthMax == 0.0f)
		{
			continue;
		}
		getRParams(&(params->offsets));
		if (params->offsets.sigmaRresidual < 0)
		{
			/* rparams found no solution -- latch rFile NULL so later passes skip it too,
			   the same convention make3DOffsets uses. */
			fprintf(stderr, "make3DOffsetsJoint: %s has no baseline solution (sigma<0) -- skipping\n",
					params->offsets.rFile);
			params->offsets.rFile = NULL;
			continue;
		}
		fprintf(stderr, "\033[1;34moffImage %s  %3.0f days  %4i / %4i\033[0m\n",
				params->offsets.rFile, params->nDays, ii, nTotal);
		/* Only one image resident at a time, so the ASCENDING buffer is used throughout. */
		readRangeOrRangeOffsets(&(params->offsets), ASCENDING, azimuthMin, azimuthMax);
		/* MANDATORY serial pre-init: svInterpBnBp lazily mallocs global SVD workspace and
		   races if first called from threads.  See the header note. */
		if (params->offsets.deltaB != DELTABNONE)
		{
			double bnS, bpS;
			svInterpBnBp(offImage, &(params->offsets), 0.0, &bnS, &bpS);
		}
		nUsed++;
		iMinAll = min(iMinAll, iMin);
		jMinAll = min(jMinAll, jMin);
		iMaxAll = max(iMaxAll, iMax);
		jMaxAll = max(jMaxAll, jMax);
		{
			int t;
			for (t = 0; t < nthreads; t++)
			{
				localImgs[t] = *offImage;
			}
		}
#pragma omp parallel
		{
			int myThread = omp_get_thread_num();
			inputImageStructure *myImg = &localImgs[myThread];
			double lat, lon, x, y, zWGS84, zSp;
			double range, azimuth, Range, myReH, theta, thetaD, psi;
			double delta, demError, sig2Base, sigmaR;
			double hAngle, xyAngle, gamma, savedLastTime;
			double B[2][2], dzdx, dzdy;
			double ax, ay, scale, P, sigma, w, dzdtSubmergence;
			float ionCorrection;
			unsigned char sMask;
			int32_t jj, iRow;
#pragma omp for schedule(dynamic, 8)
			for (iRow = iMin; iRow < iMax; iRow++)
			{
				if ((iRow % 100) == 0)
				{
					fprintf(stderr, "\t--+ %i\n", iRow);
				}
				y = (outputImage->originY + iRow * outputImage->deltaY) * MTOKM;
				for (jj = jMin; jj < jMax; jj++)
				{
					x = (outputImage->originX + jj * outputImage->deltaX) * MTOKM;
					xyToLLProj(x, y, &lat, &lon, &(outputImage->proj));
					zWGS84 = getXYHeight(lat, lon, dem, 0.0, ELLIPSOIDAL);
					if (!(zWGS84 > MINELEVATION))
					{
						continue;
					}
					zSp = sphericalElev(zWGS84, lat, Re);
					llToImageNew(lat, lon, zWGS84, &range, &azimuth, myImg);
					geometryInfo(cp, myImg, azimuth, range, zSp, thetaC, &myReH, &Range, &theta, &thetaD, &psi, zSp);
					delta = interpRangeOffsetInMeters(range, azimuth, &(params->offsets), myImg,
													  Range, thetaD, rSLPixSize, theta, &demError);
					/*  Ionosphere: the offsets correction is PRE-NEGATED and therefore ADDED.
					    Opposite to the phase path.  See the header note. */
					{
						double rangeSLC, azimuthSLC;
						computeSLCFromMLCoords(myImg, range, azimuth, &rangeSLC, &azimuthSLC);
						if (params->offsets.rOffCorrection.rangeOffsetCorrection != NULL)
						{
							ionCorrection = interpolateOffsetIonCorrectionInPixels(
								&(params->offsets.rOffCorrection), rangeSLC, azimuthSLC, -LARGEINT, 0.0);
						}
						else
						{
							ionCorrection = 0.0;
						}
					}
					if (ionCorrection > -0.98 * LARGEINT && delta > -LARGEINT)
					{
						delta += ionCorrection * rSLPixSize;
					}
					sMask = GROUNDED;
					if (shelfMask != NULL)
					{
						sMask = getShelfMask(shelfMask, x, y);
					}
					if (sMask == NOSOLUTION)
					{
						continue;
					}
					if (!(delta > (-LARGEINT + 1)))
					{
						continue;
					}
					if (sMask == GROUNDINGZONE)
					{
						continue;
					}
					if (sMask == SHELF && outputImage->noTide == TRUE)
					{
						continue;
					}
					/*  Error budget.  All four terms are metres of slant range at this point --
					    interpRangeSigma has already been multiplied by rSLPixSize -- so they add
					    in quadrature directly.  The sigma<0 sentinel was filtered above. */
					sig2Base = computeSig2Base(sin(thetaD), cos(thetaD), azimuth, myImg, &(params->offsets));
					sigmaR = interpRangeSigma(range, azimuth, &(params->offsets), myImg, Range, thetaD, rSLPixSize);
					sigmaR = sqrt(sigmaR * sigmaR + demError * demError + sig2Base +
								  rangeAccuracyVar(&(params->offsets)));
					if (sMask == SHELF)
					{
						/* twok = 1.0 for offsets (metres, not phase) */
						interpTideError(&sigmaR, myImg, params, x, y, psi, 1.0);
						delta -= -myImg->tideCorrection * cos(psi) * (double)params->nDays / 365.25;
					}
					if (vCorrect != NULL)
					{
						dzdtSubmergence = interpVCorrect(x, y, vCorrect);
						delta -= -dzdtSubmergence * cos(psi) * (double)params->nDays / 365.25;
					}
					/*  Per-image sensitivity row.  NO squint: offsets are self-consistent under
					    the zero-Doppler condition, which is why make3DOffsets passes FALSE to
					    computeA.  Save/restore lastTime around computeHeading exactly as
					    computeA does -- myImg is the per-thread copy, so this is thread-safe. */
					savedLastTime = myImg->lastTime;
					hAngle = computeHeading(lat, lon, 0.0, myImg, &(myImg->cpAll));
					myImg->lastTime = savedLastTime;
					/* xyAngle = PI/2 + meridian convergence.  For polar stereographic this returns
					   the original atan2(-y,-x), plus PI in the south, bit for bit; for UTM it uses
					   the true convergence.  Algebraically they are the same formula. */
					xyAngle = grimpXYAngle(lat, lon, x, y, &(outputImage->proj));
					gamma = xyAngle - hAngle;
					computeB(x, y, zWGS84, B, &dzdx, &dzdy, psi, psi, (xyDEM *)dem, &(outputImage->proj));
					if (sMask == SHELF)
					{
						/* Zero slope coupling on shelves; dzdx/dzdy kept for vz in pass 2. */
						B[0][0] = 0.0;
						B[0][1] = 0.0;
					}
					ax = cos(gamma) - B[0][0];
					ay = sin(gamma) - B[0][1];
					scale = 365.25 / (params->nDays * sin(psi));
					P = delta * scale;
					sigma = sigmaR * scale;
					if (!(sigma > 0.0))
					{
						continue;
					}
					w = (sumW != NULL) ? offImage->weight / (sigma * sigma) : 1.0 / (sigma * sigma);
					if (!(w > 0.0))
					{
						continue;
					}
					Nxx[iRow][jj] += (float)(w * ax * ax);
					Nxy[iRow][jj] += (float)(w * ax * ay);
					Nyy[iRow][jj] += (float)(w * ay * ay);
					bxAcc[iRow][jj] += (float)(w * ax * P);
					byAcc[iRow][jj] += (float)(w * ay * P);
					SddAcc[iRow][jj] += (float)(w * P * P);
					nObs[iRow][jj] += 1.0f;
					if (sumW != NULL)
					{
						sumW[iRow][jj] += (float)w;
						sumWT[iRow][jj] += (float)(w * (tOffCenter - tCenter));
					}
#pragma omp atomic write
					offImage->used = TRUE;
				} /* j loop */
			}	  /* i loop */
		}		  /* End omp parallel */
	}			  /* End image loop */
	fprintf(stderr, "make3DOffsetsJoint: %i of %i images contributed\n", nUsed, nTotal);
	/*
	  ****************** PASS 2: solve once per pixel *****************************************
	*/
	nSolved = 0;
	nRejCond = 0;
	nRejNoData = 0;
	nObsSum = 0.0;
	nObsMax = 0.0;
	if (nUsed > 0 && iMaxAll > iMinAll && jMaxAll > jMinAll)
	{
#pragma omp parallel for schedule(dynamic, 8) \
	reduction(+ : nSolved, nRejCond, nRejNoData, nObsSum) reduction(max : nObsMax)
		for (i = iMinAll; i < iMaxAll; i++)
		{
			double lat, lon, x, y, zWGS84;
			double nxx, nxy, nyy, bx, by, det, trHalf, disc, lambdaMin;
			double vx, vy, vz, scX, scY, deltaOffCenter;
			double B[2][2], dzdx, dzdy;
			int32_t jj;
			y = (outputImage->originY + i * outputImage->deltaY) * MTOKM;
			for (jj = jMinAll; jj < jMaxAll; jj++)
			{
				nxx = (double)Nxx[i][jj];
				nxy = (double)Nxy[i][jj];
				nyy = (double)Nyy[i][jj];
				bx = (double)bxAcc[i][jj];
				by = (double)byAcc[i][jj];
				det = nxx * nyy - nxy * nxy;
				if (nObs[i][jj] < 1.5)
				{
					nRejNoData++;
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					continue;
				}
				trHalf = 0.5 * (nxx + nyy);
				disc = sqrt((0.5 * (nxx - nyy)) * (0.5 * (nxx - nyy)) + nxy * nxy);
				lambdaMin = trHalf - disc;
				/*  Rejection.  The n-NORMALISED test is the quality criterion --
				    sigmaWorst*sqrt(n) is the effective per-measurement sigma, so it asks
				    whether the DATA is noisy rather than whether coverage is thin.  The
				    absolute test is retained but is really a precision requirement: it
				    falls almost entirely on low-n pixels (see mosaic3d.h). */
				if (!(det > 0.0) || !(lambdaMin > 0.0) ||
					(jointMaxSigmaRange > 0.0 &&
					 !((1.0 / sqrt(lambdaMin)) * sqrt((double)nObs[i][jj]) <= jointMaxSigmaRange)))
				{
					nRejCond++;
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					continue;
				}
				vx = (nyy * bx - nxy * by) / det;
				vy = (nxx * by - nxy * bx) / det;
				scX = det / nyy;
				scY = det / nxx;
				if (jointErrScale != 1.0)
				{
					/* Caller-supplied 1-sigma calibration; a SIGMA scale, so it enters the
					   variance squared.  See mosaic3d.h. */
					double fInf = jointErrScale * jointErrScale;
					scX /= fInf;
					scY /= fInf;
				}
				if (!(scX > -1000. && scX < 1000.) || !(scY > -1000. && scY < 1000.))
				{
					nRejCond++;
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					continue;
				}
				x = (outputImage->originX + jj * outputImage->deltaX) * MTOKM;
				xyToLLProj(x, y, &lat, &lon, &(outputImage->proj));
				zWGS84 = getXYHeight(lat, lon, dem, 0.0, ELLIPSOIDAL);
				computeB(x, y, zWGS84, B, &dzdx, &dzdy, 1.0, 1.0, (xyDEM *)dem, &(outputImage->proj));
				vz = vx * dzdx + vy * dzdy;
				/*  Reduced chi-square: chi2 = Sdd - v.b (cross terms cancel at the solution).
				    ~1 => the measurements agree with each other to within their own sigmas. */
				if (chi2Plane != NULL)
				{
					double chi2 = (double)SddAcc[i][jj] - (vx * bx + vy * by);
					double dof = (double)nObs[i][jj] - 2.0;
					chi2Plane[i][jj] = (dof > 0.0 && chi2 > 0.0) ? (float)(chi2 / dof) : 0.0f;
				}
				vxTmp[i][jj] = (float)(vx * scX);
				vyTmp[i][jj] = (float)(vy * scY);
				if (outputImage->makeTies == TRUE)
				{
					vzTmp[i][jj] = (float)vz;
				}
				else if (outputImage->timeOverlapFlag == TRUE)
				{
					deltaOffCenter = (sumW[i][jj] > 0.0f) ? (double)sumWT[i][jj] / (double)sumW[i][jj] : 0.0;
					vzTmp[i][jj] = (float)(deltaOffCenter * sqrt(scX * scY));
				}
				else
				{
					vzTmp[i][jj] = (float)vz;
				}
				sxTmp[i][jj] = (float)scX;
				syTmp[i][jj] = (float)scY;
				fScale[i][jj] = 1.0;
				nSolved++;
				nObsSum += (double)nObs[i][jj];
				if ((double)nObs[i][jj] > nObsMax)
				{
					nObsMax = (double)nObs[i][jj];
				}
			} /* j loop */
		}	  /* i loop */
		{
			if (fl > 0)
			{
				computeScaleLS((float **)vxTmp, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0,
							   (double)(-LARGEINT), iMinAll, iMaxAll, jMinAll, jMaxAll);
			}
			redoNormalization(1.0, outputImage, iMinAll, iMaxAll, jMinAll, jMaxAll, vXimage, vYimage, vZimage,
							  errorX, errorY, scaleX, scaleY, scaleZ, fScale, vxTmp, vyTmp, vzTmp, sxTmp, syTmp, FALSE);
		}
	}
	fprintf(stderr, "make3DOffsetsJoint: solved %ld  rejCond %ld  rejNoData %ld  meanNobs %.2f  maxNobs %.0f\n",
			(long)nSolved, (long)nRejCond, (long)nRejNoData,
			(nSolved > 0) ? nObsSum / (double)nSolved : 0.0, nObsMax);
	fprintf(outputImage->fpLog, "; make3DOffsetsJoint nImagesUsed : %i of %i\n", nUsed, nTotal);
	fprintf(outputImage->fpLog, "; make3DOffsetsJoint nSolved     : %ld\n", (long)nSolved);
	fprintf(outputImage->fpLog, "; make3DOffsetsJoint rejCond     : %ld\n", (long)nRejCond);
	fprintf(outputImage->fpLog, "; make3DOffsetsJoint meanNobs    : %.3f\n",
			(nSolved > 0) ? nObsSum / (double)nSolved : 0.0);
	fprintf(outputImage->fpLog, "; make3DOffsetsJoint maxSigma    : %f\n", jointMaxSigmaRange);
	fprintf(outputImage->fpLog, "; make3DOffsetsJoint errScale    : %f\n", jointErrScale);
	fprintf(outputImage->fpLog, "; make3DOffsetsJoint note        : joint solve, no pairs -- pairOverCount/rhoOffsets not applicable\n");

	free(localImgs);
	freeImageOff(Nxx, outputImage->ySize);
	freeImageOff(Nxy, outputImage->ySize);
	freeImageOff(Nyy, outputImage->ySize);
	freeImageOff(bxAcc, outputImage->ySize);
	freeImageOff(byAcc, outputImage->ySize);
	if (sumW != NULL)
	{
		freeImageOff(sumW, outputImage->ySize);
		freeImageOff(sumWT, outputImage->ySize);
	}
	/* nObs ownership passes to outputImage for the ".nobs" diagnostic band; mosaic3d frees it.
	   If the phase round already left one, keep that one and drop this -- the bands would
	   otherwise collide and the phase count is the one §9's calibration refers to. */
	if (outputImage->jointNObs == NULL)
	{
		outputImage->jointNObs = nObs;
	outputImage->jointChi2 = chi2Plane;
	}
	else
	{
		freeImageOff(nObs, outputImage->ySize);
		freeImageOff(chi2Plane, outputImage->ySize);
	}
	{
		extern double totalOffsetsIOTime;
		gettimeofday(&funcEnd, NULL);
		fprintf(stderr, "Total offsets I/O time (whole run): %.3f s\n", totalOffsetsIOTime);
		fprintf(stderr, "Total offsets processing time (whole run): %.3f s\n",
				(funcEnd.tv_sec - funcStart.tv_sec) + (funcEnd.tv_usec - funcStart.tv_usec) * 1e-6);
	}
	endScale(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, FALSE);
	fprintf(outputImage->fpLog, ";\n; Returning from make3DOffsetsJoint(.c)\n");
	fflush(outputImage->fpLog);
}

/* mallocImage() allocates per row and does NOT zero. */
static float **mallocZeroImageOff(int32_t nr, int32_t nc)
{
	float **tmp;
	int32_t i, j;
	tmp = mallocImage(nr, nc);
	for (i = 0; i < nr; i++)
	{
		for (j = 0; j < nc; j++)
		{
			tmp[i][j] = 0.0f;
		}
	}
	return tmp;
}

static void freeImageOff(float **image, int32_t nr)
{
	int32_t i;
	if (image == NULL)
	{
		return;
	}
	for (i = 0; i < nr; i++)
	{
		free(image[i]);
	}
	free(image);
}

