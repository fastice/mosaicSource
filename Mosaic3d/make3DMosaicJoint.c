#include "stdio.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "cRecipes/nrutil.h"
#include "mosaicSource/common/common.h"
#include "mosaic3d.h"
#include "sys/time.h"

/*
****************** Estimate 3D velocity from phase -- JOINT (normal-equations) ***************

  Drop-in alternative to make3DMosaic() (make3DMosaic.c) with the identical signature.
  THE DEFAULT since 2026-08-29; with "mosaic3d -legacyPairPhase" the pair path in
  make3DMosaic.c runs exactly as before.  Companion flags: -jointMaxSigma (rejection
  threshold, see CONDITIONING) and -jointErrScale (caller-supplied 1-sigma calibration).

  WHY
  ---
  make3DMosaic() forms every admissible ascending x descending pair at a pixel and accumulates
  each as an independent observation.  A pixel covered by n_A ascending and n_D descending
  products yields n_A*n_D pairs built from only n_A + n_D independent measurements.  The shared
  measurement's error is common-mode across every pair containing it, so averaging drives the
  variance to a floor rather than to zero, while the accumulated weight grows like the pair
  count -- so the crossing-orbit method is over-weighted against the other four methods in the
  mosaic blend.  See mosaicSource/CLAUDE.md "Crossing-pair error over-counting" (which ships
  -rhoPhase/-rhoOffsets as a partial fix and notes it "corrects the reported error, NOT the
  weighting"), and Documents/crossingOrbitRedundancy.md for the full derivation.

  This routine forms NO pairs.  Each measurement enters the solution exactly once, so
  double-counting is impossible by construction.

  THE MATH
  --------
  computeA() (common/initRoutines.c) is exactly N^-1 for

      N = [ cos(beta)        sin(beta)       ]     beta       = xyAngle - H_A
          [ cos(alpha+beta)  sin(alpha+beta) ]     alpha+beta = xyAngle - H_D

  (verify: A[0][0] = (cos b - cos a cos(a+b))/sin^2 a = sin(a+b)/sin a, and likewise for the
  other three entries).  Each ROW depends on exactly ONE image's heading, so a per-image
  sensitivity row exists and no pairing is needed to obtain it.

  computeVxy() computes v = (I - AB)^-1 A P.  With A = N^-1 that collapses EXACTLY (not to
  first order) to

      v = (I - N^-1 B)^-1 N^-1 P = [N (I - N^-1 B)]^-1 P = (N - B)^-1 P

  so row i of the system is P_i = a_i . v with

      a_i     = ( cos(gamma_i) - dzdx/tan(psi_i),  sin(gamma_i) - dzdy/tan(psi_i) )
      gamma_i = xyAngle - H_i
      P_i     = phase_i * 365.25 / (twok_i * nDays_i * sin(psi_i))
      sigma_i = phaseError_i * (that same scale factor)

  Accumulate over all images, each entering once, with w_i = weight_i / sigma_i^2:

      Nxx += w ax ax    Nxy += w ax ay    Nyy += w ay ay
      bx  += w ax P     by  += w ay P

  and solve once per pixel:

      det = Nxx*Nyy - Nxy^2
      vx  = (Nyy*bx - Nxy*by) / det        var(vx) = Nyy/det   ->  scX = det/Nyy
      vy  = (Nxx*by - Nxy*bx) / det        var(vy) = Nxx/det   ->  scY = det/Nxx

  The covariance is N^-1 exactly, which is what the mosaic accumulator already expects from
  computeVxy() (it wants 1/variance in scaleX/scaleY).

  CONDITIONING -- and why it is NOT a generalization of the |alpha| gate
  ----------------------------------------------------------------------
  The obvious move is a scale-free version of the existing gates (computeA()'s |alpha| > 0.8
  and computeVxy()'s detC > 0.5), e.g.

      cond = det(N) / (trace(N)/2)^2        <-- WRONG, do not reinstate

  which does reduce to sin^2(alpha) at n = 2 with equal weights.  It is still wrong, because
  it is NOT MONOTONE under adding data.  Measured, at a perfect 90-degree crossing:

      nA=1,  nD=1   cond 1.000   accept        worst-direction sigma 1.000
      nA=10, nD=1   cond 0.331   REJECT        worst-direction sigma 1.000  <-- same solution!
      nA=20, nD=1   cond 0.181   REJECT        worst-direction sigma 1.000

  Adding ascending images leaves the achieved precision unchanged (or better) while driving
  cond down, because trace grows with every row but det(N) only grows through CROSSING pairs.
  A pixel can therefore flip accept -> reject purely by acquiring more data.  On the first
  Greenland test run this gate discarded 22.9% of otherwise-solvable pixels.

  Any normalized, purely geometric measure has this defect: normalizing away "how much data"
  necessarily turns the test into AVERAGE geometry quality, and redundant same-direction data
  always drags an average down while never hurting the estimate.  Monotonicity and scale-free
  geometry are incompatible here.

  So the gate is on ACHIEVED PRECISION instead, which is monotone by construction because
  N grows by a positive-semidefinite rank-1 term with every measurement:

      lambdaMin = smallest eigenvalue of N
      sigmaWorst = 1/sqrt(lambdaMin)        worst-direction formal sigma, m/yr
      reject if sigmaWorst * sqrt(n) > -jointMaxSigma (default 35 m/yr)

  max over unit u of u^T N^-1 u = 1/lambdaMin, so this is exactly "the formal error in the
  worst-constrained direction", in physical units, and it cannot be made to fail by adding
  data.  It deliberately does NOT reproduce the |alpha| > 0.8 gate -- that gate is geometry-only
  and so cannot be monotone.  Note the SOLUTION still reduces exactly to the pair solution at
  n = 2 (see THE MATH above); only the rejection rule differs.

  The threshold is intentionally permissive.  This is an inverse-variance mosaic: a poorly
  constrained pixel gets a correspondingly huge reported sigma and is down-weighted
  automatically wherever it is averaged, which is more informative than a hard reject.

  The Lagrange identity behind the original idea still holds and is still the reason the joint
  solve sees geometry no pair does:

      det(N) = sum_{i<j} w_i w_j |a_i|^2 |a_j|^2 sin^2(psi_ij)

  i.e. three tracks at modest separation jointly constrain a pixel that no single pair could.

  DIFFERENCES FROM make3DMosaic() -- read before comparing outputs
  ---------------------------------------------------------------
  * timeThreshPhase is a PER-IMAGE window against the mosaic centre, not a pair separation.
    An image contributes if |julDay - tCenter| <= timeThreshPhase.  A pair that passed the old
    |jdA - jdD| <= thresh test can now fail if BOTH images sit far from tCenter on the same side.
  * sepAscDesc is meaningless here (no pairs to restrict); any image whose geometry helps the
    conditioning contributes, including two ascending images.
  * pairOverCount / rhoPhase / inflatePairOverCount() are not used and not needed.
  * getIntersect() and computeSceneAlpha() are pair-overlap constructs; per-image getRegion()
    bounds replace them.
  * Feathering: computeScaleLS() runs once over the union region rather than once per pair.
  * timeOverlapFlag folds image->weight into w_i in place of the pair's sqrt(wA*wD), and takes
    the weighted-mean time offset via the sumW/sumWT accumulators.  Outside that mode
    image->weight is ignored, exactly as make3DMosaic ignores it (combWeight = 1.0).

  TWO THINGS TO KNOW BEFORE TRUSTING THE OUTPUT
  ---------------------------------------------
  * Cov = N^-1 assumes the measurements are INDEPENDENT.  They are not: two date-pairs sharing
    an acquisition inherit that frame's baseline residual and ionospheric screen, so Cov(d) is
    not diagonal.  Removing the combinatorial over-count does not remove the shared-frame
    correlation underneath it.

    MEASURED, Greenland 2026-08-28 (Documents/crossingOrbitRedundancy.md §8): against the
    pair baseline at rho = 0.6, std(d) vs Sentinel-1 moved only 3.014 -> 2.960 (1.8% better)
    while RMS(e) fell 2.995 -> 0.994.  So k = std(d)/RMS(e) went 1.01 -> 2.98: THE REPORTED
    ERROR IS ~3x TOO SMALL.  The velocity was already right; only the error bar changed.
    (crossingOrbitRedundancyReview.md §2 predicted k ~ 2.3 by extrapolation.)

    => DO NOT SHIP THIS AS A PRODUCT ERROR BUDGET until the shared-frame term is modelled:
    Cov(d) = D + sum_f sigma_f^2 u_f u_f^T, diagonal-plus-low-rank, inverted by Woodbury and
    used as the GLS weight.  sigma_f already exists as vhParam->sigma (phase) /
    offsets.sigmaRresidual (offsets).  Not implemented here.
  * NO OUTLIER REJECTION -- and measurement says none is needed.  Huber IRLS was implemented and
    had NO effect on 4 of 4 real cases.  The ".chi2" band explains why: reduced chi-square ~0.85,
    so the measurements at a pixel already agree with each other to within their own sigmas.
    There is no inconsistent minority to down-weight; the tail is common-mode error at the pixel.
    Do not re-add reweighting -- see Documents/crossingOrbitRedundancy.md §13.

  Interaction with singleImageFastPath (mosaic3d.c): that path requires no3d == TRUE, and the
  no3d early return below happens BEFORE setupBuffers()/computeScale(), exactly as in
  make3DMosaic -- so the buffers it leaves NULL are never touched here either.

  computePhiZM3d(), computePhiFlatEarthM3d() and setBuffer() below are verbatim copies of the
  statics in make3DMosaic.c -- they cannot be shared without editing that file.  KEEP IN SYNC
  if the originals change.
*/

/* The gate defaults (JOINTMAXSIGMAPHASEDEF / JOINTMAXSIGMARANGEDEF) live in mosaic3d.h,
   where mosaic3d.c's definitions of the globals can see them too. */

static void setBuffer(inputImageStructure *inputImage, float *buf);
static double computePhiZM3d(double *thetaD, double z, double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
							 double Range, double Re, double ReH, double ReHfixed, double thetaC, double thetaCfixedReH, double *phaseError);
static double computePhiFlatEarthM3d(double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
									 double Range, double Re, double ReHfixed, double thetaCfixedReH, double *phaseError);
static float **mallocZeroImage(int32_t nr, int32_t nc);
static void freeImage(float **image, int32_t nr);

void make3DMosaicJoint(inputImageStructure *ascImages, inputImageStructure *descImages,
					   vhParams *ascParams, vhParams *descParams, xyDEM *dem,
					   outputImageStructure *outputImage, float fl, int32_t no3d, float timeThreshPhase)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	extern float *AImageBuffer;
	extern float *AIonBuffer;
	extern int32_t indentRegionOutput;
	extern int32_t useSquint;
	conversionDataStructure *cp;			 /* coordinate conversion info for the current image */
	inputImageStructure *allImages, *phaseImage; /* consolidated asc+desc list, current image */
	vhParams *params;
	ShelfMask *shelfMask;
	xyDEM *vCorrect;
	double Re, ReH, thetaC, ReHfixed, thetaCfixedReH; /* geometry params for the current image */
	double twok;									  /* 4pi/lambda for the current image */
	double tCenter, tOffCenter;						  /* mosaic centre time, image centre time */
	float geoTolerance = 1e-3;						  /* as make3DMosaic: a few metres is plenty */
	float **vXimage, **vYimage, **vZimage, **errorX, **errorY;
	float **vxTmp, **vyTmp, **vzTmp, **fScale, **sxTmp, **syTmp;
	float **scaleX, **scaleY, **scaleZ;
	/* Normal-equation accumulators, one plane each, over the whole output grid.
	   float (not double) is safe: cancellation in Nxx*Nyy - Nxy^2 only bites as cond -> 0, and
	   those pixels are rejected by the conditioning gate anyway.  det is formed in double. */
	float **Nxx, **Nxy, **Nyy, **bxAcc, **byAcc;
	float **nObs;				 /* contributing measurements per pixel -- diagnostic */
	float **SddAcc; /* sum w d^2, for chi-square at solve time */
	float **chi2Plane;
	float **sumW, **sumWT;		 /* timeOverlapFlag only; NULL otherwise */
	float dum1, dum2;
	float azimuthMin, azimuthMax;
	int32_t iMin, iMax, jMin, jMax;					 /* current image's region */
	int32_t iMinAll, iMaxAll, jMinAll, jMaxAll;		 /* union over all contributing images */
	int32_t ii, nTotal, nUsed;
	int32_t i;
	int64_t nSolved, nRejCond, nRejNoData;
	double nObsSum, nObsMax;
	struct timeval funcStart, funcEnd;
	gettimeofday(&funcStart, NULL);

	fprintf(outputImage->fpLog, ";\n; Entering make3DMosaicJoint(.c)\n");
	if (no3d == TRUE)
	{
		fprintf(outputImage->fpLog, ";\n; no3d flag set, returning;\n; Returning from make3DMosaicJoint(.c)\n");
		return;
	}
	shelfMask = outputImage->shelfMask;
	vCorrect = outputImage->verticalCorrection;
	if (dem->stdLat < 50 || dem->stdLat > 80)
	{
		error("mosaic3d invalid slat for dem");
	}
	/*
	  Pointers to output images
	*/
	setupBuffers(outputImage, &vXimage, &vYimage, &vZimage, &scaleX, &scaleY, &scaleZ,
				 &vxTmp, &vyTmp, &vzTmp, &sxTmp, &syTmp, &fScale, &errorX, &errorY);
	/*
	  Compute feather scale for existing, and undo normalization.  Identical to make3DMosaic --
	  this routine is still one accumulation round among several (Landsat -> here -> offsets ->
	  vh -> speckle), so it must compose the same way.
	*/
	computeScale((float **)vXimage, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0, (double)(-LARGEINT));
	undoNormalization(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, fScale, FALSE);
	/*
	  Allocate and zero the accumulators
	*/
	Nxx = mallocZeroImage(outputImage->ySize, outputImage->xSize);
	Nxy = mallocZeroImage(outputImage->ySize, outputImage->xSize);
	Nyy = mallocZeroImage(outputImage->ySize, outputImage->xSize);
	bxAcc = mallocZeroImage(outputImage->ySize, outputImage->xSize);
	byAcc = mallocZeroImage(outputImage->ySize, outputImage->xSize);
	nObs = mallocZeroImage(outputImage->ySize, outputImage->xSize);
	SddAcc = mallocZeroImage(outputImage->ySize, outputImage->xSize);
	sumW = NULL;
	sumWT = NULL;
	chi2Plane = mallocZeroImage(outputImage->ySize, outputImage->xSize);
	if (outputImage->timeOverlapFlag == TRUE)
	{
		sumW = mallocZeroImage(outputImage->ySize, outputImage->xSize);
		sumWT = mallocZeroImage(outputImage->ySize, outputImage->xSize);
	}
	allImages = ascImages; /* consolodateLists() has already appended desc to asc */
	params = ascParams;
	nTotal = 0;
	for (phaseImage = allImages; phaseImage != NULL; phaseImage = phaseImage->next)
	{
		nTotal++;
	}
	fprintf(stderr, "make3DMosaicJoint: nTotal Images %i\n", nTotal);
	tCenter = (outputImage->jd1 + outputImage->jd2 + 1.) * 0.5; /* as make3DMosaic */
	iMinAll = outputImage->ySize;
	jMinAll = outputImage->xSize;
	iMaxAll = 0;
	jMaxAll = 0;
	nUsed = 0;
	int nthreads = omp_get_max_threads();
	inputImageStructure *localImgs = (inputImageStructure *)malloc((size_t)nthreads * sizeof(inputImageStructure));
	if (localImgs == NULL)
	{
		error("make3DMosaicJoint: malloc failed for per-thread image copies\n");
	}
	/*
	  ****************** PASS 1: accumulate normal equations, one image at a time **************
	*/
	ii = 0;
	params = ascParams;
	for (phaseImage = allImages; phaseImage != NULL; phaseImage = phaseImage->next, params = params->next)
	{
		ii++;
		tOffCenter = phaseImage->julDay + params->nDays * 0.5;
		if (strstr(phaseImage->file, "nophase") != NULL)
		{
			continue;
		}
		/* Per-image temporal window.  NOT the pair gate -- see the header note. */
		if (fabs(phaseImage->julDay - tCenter) > timeThreshPhase)
		{
			fprintf(stderr, "\t%s skipped: |julDay - tCenter| = %.1f > %.1f\n",
					phaseImage->file, fabs(phaseImage->julDay - tCenter), timeThreshPhase);
			continue;
		}
		indentRegionOutput = FALSE;
		if (!getRegion(phaseImage, &iMin, &iMax, &jMin, &jMax, outputImage))
		{
			continue;
		}
		if (iMin > iMax || jMin > jMax)
		{
			continue;
		}
		{
			double sq = (useSquint && phaseImage->hasSquintPolynomial)
							? evaluateSquint(phaseImage, phaseImage->rangeSize * 0.5, phaseImage->azimuthSize * 0.5)
							: 0.0;
			fprintf(stderr, "\033[1;34mphaseImage %s %3.0f -- %5.3f -- %4i / %4i  squint %s %5.2f\033[0m\n",
					phaseImage->file, params->nDays, phaseImage->par.lambda, ii, nTotal,
					useSquint ? "on" : "off", sq);
		}
		/*  Set buffer, memory channel, and read image.  Only one image is resident at a time, so
		    the D-side buffers (DImageBuffer/DIonBuffer, MEM2) are never needed here. */
		setBuffer(phaseImage, AImageBuffer);
		cp = setupGeoConversions(phaseImage, &dum1, &dum2, &Re, &ReH, &thetaC, &ReHfixed, &thetaCfixedReH);
		phaseImage->tolerance = geoTolerance;
		getAzimuthBoundsForXYBox(iMin, iMax, jMin, jMax, phaseImage, outputImage, &azimuthMin, &azimuthMax);
		if (azimuthMin == 0.0f && azimuthMax == 0.0f)
		{
			continue;
		}
		getMosaicInputImage(phaseImage, azimuthMin, azimuthMax);
		getIonospherePhaseImage(phaseImage, AIonBuffer, azimuthMin, azimuthMax);
		phaseImage->memChan = MEM1;
		twok = (4.0 * PI) / phaseImage->par.lambda;
		nUsed++;
		iMinAll = min(iMinAll, iMin);
		jMinAll = min(jMinAll, jMin);
		iMaxAll = max(iMaxAll, iMax);
		jMaxAll = max(jMaxAll, jMax);
		{
			int t;
			for (t = 0; t < nthreads; t++)
			{
				localImgs[t] = *phaseImage;
			}
		}
#pragma omp parallel
		{
			int myThread = omp_get_thread_num();
			inputImageStructure *myImg = &localImgs[myThread];
			double lat, lon, x, y, zWGS84, zSp;
			double range, azimuth, Range, myReH, theta, thetaD, psi;
			double phase, ionPhase, phiZ, phaseError;
			double hAngle, xyAngle, gamma, savedLastTime;
			double B[2][2], dzdx, dzdy;
			double ax, ay, scale, P, sigma, w, dzdtSubmergence;
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
					interpPhaseImage(myImg, range, azimuth, &phase);
					sMask = GROUNDED;
					if (shelfMask != NULL)
					{
						sMask = getShelfMask(shelfMask, x, y);
					}
					if (sMask == NOSOLUTION)
					{
						continue;
					}
					if (!(phase > -LARGEINT))
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
					/*
					  Phase due to topography (or the flat-earth orbit ramp for ISCE/NISAR).
					  Identical to make3DMosaic, single-sided.
					*/
					if (params->applyFlatEarth)
					{
						phiZ = computePhiFlatEarthM3d(azimuth, params, myImg, Range, Re, ReHfixed, thetaCfixedReH, &phaseError);
					}
					else
					{
						phiZ = computePhiZM3d(&thetaD, zSp, azimuth, params, myImg, Range, Re, myReH, ReHfixed, thetaC, thetaCfixedReH, &phaseError);
					}
					phase = phase - phiZ;
					/*  Ionosphere correction, still in radians, before the velocity scaling.  The
					    file holds the ionospheric phase in the same sign convention as the phase
					    image, so it is SUBTRACTED -- unlike the range-offset correction in
					    make3DOffsets.c, which is pre-negated and added. */
					if (myImg->ionospherePhase != NULL)
					{
						interpIonPhaseImage(myImg, range, azimuth, &ionPhase);
						if (ionPhase > -0.98 * LARGEINT)
						{
							phase -= ionPhase;
						}
					}
					/*  Tide correction */
					if (sMask == SHELF)
					{
						interpTideError(&phaseError, myImg, params, x, y, psi, twok);
						phase -= -myImg->tideCorrection * cos(psi) * twok * (double)params->nDays / 365.25;
					}
					/*  Submergence / SMB correction.  A per-measurement additive term:
					    d_i -> d_i - u_iz * dzdt, with cos(psi) == u_iz. */
					if (vCorrect != NULL)
					{
						dzdtSubmergence = interpVCorrect(x, y, vCorrect);
						phase -= -dzdtSubmergence * cos(psi) * twok * (double)params->nDays / 365.25;
					}
					/*
					  Per-image sensitivity row.  computeA() builds N from the two images'
					  headings; each of its rows depends on one heading only, so the row for this
					  image alone is (cos(gamma), sin(gamma)) with gamma = xyAngle - hAngle.
					  computeHeading() calls llToImageNew() at lat+-dlat, clobbering the
					  Newton warm-start cache -- save/restore exactly as computeA() does.
					  myImg is the per-thread copy, so this is thread-safe.
					*/
					savedLastTime = myImg->lastTime;
					hAngle = computeHeading(lat, lon, 0.0, myImg, &(myImg->cpAll));
					if (useSquint && myImg->hasSquintPolynomial)
					{
						double sqRange, sqAzimuth;
						llToImageNew(lat, lon, 0.0, &sqRange, &sqAzimuth, myImg);
						hAngle += evaluateSquint(myImg, sqRange, sqAzimuth) * DTOR;
					}
					myImg->lastTime = savedLastTime;
					/* xyAngle = PI/2 + meridian convergence.  For polar stereographic this returns
					   the original atan2(-y,-x), plus PI in the south, bit for bit; for UTM it uses
					   the true convergence.  Algebraically they are the same formula. */
					xyAngle = grimpXYAngle(lat, lon, x, y, &(outputImage->proj));
					gamma = xyAngle - hAngle;
					/*
					  Slope coupling.  computeB() is called with this image's psi in both slots so
					  row 0 is this image's own row -- reuses limitSlope()/badZ()/interpXYDEM()
					  verbatim rather than reimplementing the slope.
					*/
					computeB(x, y, zWGS84, B, &dzdx, &dzdy, psi, psi, (xyDEM *)dem, &(outputImage->proj));
					if (sMask == SHELF)
					{
						/* Zero slope coupling on ice shelves -- shelf-interior slope should be
						   near zero, so a large DEM reading is almost always a stale rift or
						   calving front.  Matches make3DOffsets.c / speckleTrackMosaic.c.
						   dzdx/dzdy themselves are left alone for the vz calculation in pass 2. */
						B[0][0] = 0.0;
						B[0][1] = 0.0;
					}
					/*  v = (N - B)^-1 P, so the row is N_i - B_i */
					ax = cos(gamma) - B[0][0];
					ay = sin(gamma) - B[0][1];
					/*  Scale phase to velocity (m/yr) */
					scale = 365.25 / (twok * params->nDays * sin(psi));
					P = phase * scale;
					sigma = phaseError * scale;
					if (!(sigma > 0.0))
					{
						continue;
					}
					/*  make3DMosaic applies image->weight ONLY via combWeight = sqrt(wA*wD) in
					    timeOverlap mode; the normal path uses combWeight = 1.0 and ignores it
					    (speckleTrackMosaic.c:141 likewise warns if weight != 1 outside that
					    mode).  Match that, so the two routines agree when all weights are 1. */
					w = (sumW != NULL) ? phaseImage->weight / (sigma * sigma) : 1.0 / (sigma * sigma);
					if (!(w > 0.0))
					{
						continue;
					}
					/*  Accumulate.  #pragma omp for partitions by row, so each pixel is touched
					    by exactly one thread within this image, and images are processed
					    sequentially -- no atomics needed. */
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
					phaseImage->used = TRUE;
				} /* j loop */
			}	  /* i loop */
		}		  /* End omp parallel */
	}			  /* End image loop */
	fprintf(stderr, "make3DMosaicJoint: %i of %i images contributed\n", nUsed, nTotal);
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
					/* 0 or 1 measurements -- two unknowns need two independent looks */
					nRejNoData++;
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					continue;
				}
				/*  Worst-direction formal sigma: max over unit u of u^T N^-1 u = 1/lambdaMin.
				    Monotone under adding measurements (N gains a PSD rank-1 term each time), so
				    a pixel can never be rejected for having MORE data -- see CONDITIONING at the
				    top for why a trace-normalized geometric test cannot have that property. */
				trHalf = 0.5 * (nxx + nyy);
				disc = sqrt((0.5 * (nxx - nyy)) * (0.5 * (nxx - nyy)) + nxy * nxy);
				lambdaMin = trHalf - disc;
				/*  Rejection.  The n-NORMALISED test is the quality criterion --
				    sigmaWorst*sqrt(n) is the effective per-measurement sigma, so it asks
				    whether the DATA is noisy rather than whether coverage is thin.  The
				    absolute test is retained but is really a precision requirement: it
				    falls almost entirely on low-n pixels (see mosaic3d.h). */
				if (!(det > 0.0) || !(lambdaMin > 0.0) ||
					(jointMaxSigma > 0.0 &&
					 !((1.0 / sqrt(lambdaMin)) * sqrt((double)nObs[i][jj]) <= jointMaxSigma)))
				{
					nRejCond++;
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					continue;
				}
				vx = (nyy * bx - nxy * by) / det;
				vy = (nxx * by - nxy * bx) / det;
				/*  Covariance is N^-1, so var(vx) = Nyy/det and var(vy) = Nxx/det. The mosaic
				    accumulator wants 1/variance, matching computeVxy()'s scaleX/scaleY. */
				scX = det / nyy;
				scY = det / nxx;
				/*  Divide the WEIGHT, not just the reported error: v = sum(v*sc)/sum(sc), so
				    scaling both sc and v*sc leaves this round's velocity untouched while
				    inflating the reported sigma AND correctly reducing this round's weight
				    against the other methods in the blend.  That second effect is the one
				    inflatePairOverCount() explicitly cannot achieve. */
				if (jointErrScale != 1.0)
				{
					/* jointErrScale is a SIGMA scale, so it enters the variance squared.
					   Caller-supplied 1-sigma calibration; see mosaic3d.h. */
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
				/*  Slopes for vz.  Recomputed here rather than stored -- two more full-grid
				    planes would cost more than one interpXYDEM stencil per solved pixel.  psi is
				    irrelevant since only dzdx/dzdy are used, so pass 1.0 twice. */
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
					/* Weighted mean time offset over the contributing images, in place of the
					   pair average 0.5*(tOffCenterA + tOffCenterD - 2*tCenter). */
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
		/*
		  Feather and fold this round into the mosaic accumulators.  One call, over the union
		  region, rather than make3DMosaic's one call per pair.
		*/
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
	free(localImgs);
	fprintf(stderr, "make3DMosaicJoint: solved %ld  rejCond %ld  rejNoData %ld  meanNobs %.2f  maxNobs %.0f\n",
			(long)nSolved, (long)nRejCond, (long)nRejNoData,
			(nSolved > 0) ? nObsSum / (double)nSolved : 0.0, nObsMax);
	fprintf(outputImage->fpLog, "; make3DMosaicJoint nImagesUsed : %i of %i\n", nUsed, nTotal);
	fprintf(outputImage->fpLog, "; make3DMosaicJoint nSolved     : %ld\n", (long)nSolved);
	fprintf(outputImage->fpLog, "; make3DMosaicJoint rejCond     : %ld\n", (long)nRejCond);
	fprintf(outputImage->fpLog, "; make3DMosaicJoint meanNobs    : %.3f\n",
			(nSolved > 0) ? nObsSum / (double)nSolved : 0.0);
	fprintf(outputImage->fpLog, "; make3DMosaicJoint maxSigma    : %f\n", jointMaxSigma);
	fprintf(outputImage->fpLog, "; make3DMosaicJoint note        : joint solve, no pairs -- pairOverCount/rho not applicable\n");
	/*
	  Free accumulators
	*/
	freeImage(Nxx, outputImage->ySize);
	freeImage(Nxy, outputImage->ySize);
	freeImage(Nyy, outputImage->ySize);
	freeImage(bxAcc, outputImage->ySize);
	freeImage(byAcc, outputImage->ySize);
	/* nObs is NOT freed: ownership passes to outputImage so mosaic3d can write the
	   ".nobs" diagnostic band.  mosaic3d frees it after writing. */
	outputImage->jointNObs = nObs;
	outputImage->jointChi2 = chi2Plane;
	if (sumW != NULL)
	{
		freeImage(sumW, outputImage->ySize);
		freeImage(sumWT, outputImage->ySize);
	}
	{
		extern double totalPhaseIOTime;
		gettimeofday(&funcEnd, NULL);
		fprintf(stderr, "Total phase I/O time (whole run): %.3f s\n", totalPhaseIOTime);
		fprintf(stderr, "Total phase processing time (whole run): %.3f s\n",
				(funcEnd.tv_sec - funcStart.tv_sec) + (funcEnd.tv_usec - funcStart.tv_usec) * 1e-6);
	}
	endScale(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, FALSE);
	fprintf(outputImage->fpLog, ";\n; Returning from make3DMosaicJoint(.c)\n");
	fflush(outputImage->fpLog);
}

/*
  mallocImage() mallocs row by row and does NOT zero -- the accumulators must start at 0.
*/
static float **mallocZeroImage(int32_t nr, int32_t nc)
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

static void freeImage(float **image, int32_t nr)
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

/*
  ***** Verbatim copies of the statics in make3DMosaic.c -- KEEP IN SYNC. *****
*/

/*
  Compute flat-earth baseline phase for ISCE/NISAR products (topo already removed).
  Mirrors the changeflat formula so the orbit-error ramp is corrected inline.
*/
static double computePhiFlatEarthM3d(double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
									 double Range, double Re, double ReHfixed,
									 double thetaCfixedReH, double *phaseError)
{
	double normAzimuth, imageLength;
	double bn, bp, bSq, delta;
	double theta, thetaDFlat;
	double sinThetaD, cosThetaD;
	double v[7], tmpV[7];
	double sig2Base;
	double twok;
	double xsq;
	int32_t i, j;

	twok = 4.0 * PI / phaseImage->par.lambda;
	imageLength = (double)phaseImage->azimuthSize;
	normAzimuth = (azimuth - 0.5 * imageLength) / imageLength;
	xsq = normAzimuth * normAzimuth;
	bn = vhParam->Bn + normAzimuth * vhParam->dBn + xsq * vhParam->dBnQ;
	bp = vhParam->Bp + normAzimuth * vhParam->dBp + xsq * vhParam->dBpQ;
	bSq = bn * bn + bp * bp;

	/* Flat-earth look angle (same formula as computePhiZM3d line for thetaDFlat) */
	theta = acos((Range * Range + ReHfixed * ReHfixed - Re * Re) / (2.0 * ReHfixed * Range));
	thetaDFlat = theta - thetaCfixedReH;

	sinThetaD = sin(thetaDFlat);
	cosThetaD = cos(thetaDFlat);
	v[1] = -twok * sinThetaD;
	v[2] = -twok * cosThetaD;
	v[3] = -twok * sinThetaD * normAzimuth;
	v[4] = -twok * cosThetaD * normAzimuth;
	v[5] = -twok * sinThetaD * xsq;
	v[6] = -twok * cosThetaD * xsq;
	for (i = 1; i <= 6; i++)
	{
		tmpV[i] = 0;
		for (j = 1; j <= 6; j++)
			tmpV[i] += vhParam->C[i][j] * v[j];
	}
	sig2Base = 0.0;
	for (j = 1; j <= 6; j++)
		sig2Base += tmpV[j] * v[j];

	delta = -bn * sinThetaD - bp * cosThetaD + bSq * 0.5 / Range;
	*phaseError = sqrt(sig2Base + min(6 * PI, vhParam->sigma) * min(6 * PI, vhParam->sigma));
	return delta * twok;
}

/*
  Compute phase due to topography
*/
static double computePhiZM3d(double *thetaD, double z, double azimuth, vhParams *vhParam, inputImageStructure *phaseImage, double Range, double Re,
							 double ReH, double ReHfixed, double thetaC, double thetaCfixedReH, double *phaseError)
{
	double normAzimuth, imageLength;
	double bn, bp, bSq, delta;
	double theta, thetaDFlat;
	double phiZ;
	double xsq;
	double sig2Base;
	double twok;
	int32_t i, j;
	double v[7], tmpV[7];
	double sinThetaD, cosThetaD;
	/*
	  Compute baseline
	*/
	twok = 4.0 * PI / phaseImage->par.lambda;
	imageLength = (double)phaseImage->azimuthSize;
	normAzimuth = (azimuth - 0.5 * imageLength) / imageLength;
	xsq = normAzimuth * normAzimuth;
	bn = vhParam->Bn + normAzimuth * vhParam->dBn + xsq * vhParam->dBnQ;
	bp = vhParam->Bp + normAzimuth * vhParam->dBp + xsq * vhParam->dBpQ;
	bSq = bn * bn + bp * bp;
	/*
	   Compute look angles for nonflat surface
	*/
	theta = acos((Range * Range + ReH * ReH - pow((Re + z), 2.0)) / (2.0 * ReH * Range));
	*thetaD = theta - thetaC;
	sinThetaD = sin(*thetaD);
	cosThetaD = cos(*thetaD);
	/*
	  Vector used to comptue baseline error
	*/
	v[1] = -twok * sinThetaD;
	v[2] = -twok * cosThetaD;
	v[3] = -twok * sinThetaD * normAzimuth;
	v[4] = -twok * cosThetaD * normAzimuth;
	v[5] = -twok * sinThetaD * xsq;
	v[6] = -twok * cosThetaD * xsq;
	/*
	  C*v
	*/
	for (i = 1; i <= 6; i++)
	{
		tmpV[i] = 0;
		for (j = 1; j <= 6; j++)
			tmpV[i] += vhParam->C[i][j] * v[j];
	}
	/*
	  sigma^2= vT * C*v
	*/
	sig2Base = 0.0;
	for (j = 1; j <= 6; j++)
		sig2Base += tmpV[j] * v[j];
	/*
	  Assume error more than 1/2 of a fringe, reflects tie point error rather than phase noise
	*/
	*phaseError = sqrt(sig2Base + min(6 * PI, vhParam->sigma) * min(6 * PI, vhParam->sigma)); /* note this is returning phase error as sigma */
	/*
	   This delta uses a varying ReH because it is the computationally correct version
	*/
	delta = sqrt(pow(Range, 2.0) - 2.0 * Range * (bn * sinThetaD + bp * cosThetaD) + bSq) - Range;
	/*
	  Compute thetaD for a flat surface
	  This delta uses a fixed ReH because ultimately it is adding back what was subtracted
	  interferogram phase = phi - phi_flat=(phit + phiv) - phiflat
	  phiz=phithat - phiflat
	  phase - phiz= (phit+phiv)-phiflat - (phithat - phiflat) = phiv + (phit-phithat)
	*/
	theta = acos((Range * Range + ReHfixed * ReHfixed - Re * Re) / (2.0 * ReHfixed * Range));
	thetaDFlat = theta - thetaCfixedReH;
	/*
	   Substract flat earth phase
	*/
	delta -= -bn * sin(thetaDFlat) - bp * cos(thetaDFlat) + bSq * 0.5 / Range;
	phiZ = delta * twok;
	return phiZ;
}

static void setBuffer(inputImageStructure *inputImage, float *buf)
{
	int32_t i;
	for (i = 0; i < inputImage->azimuthSize; i++)
	{
		inputImage->image[i] = &(buf[i * inputImage->rangeSize]);
	}
}
