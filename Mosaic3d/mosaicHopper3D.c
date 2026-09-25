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
********** THE 3D HOPPER: EVERY OBSERVABLE, THREE COMPONENTS, ONE ACCUMULATOR ****************

  "mosaic3d -hopper3D".  DEFAULT OFF.  The 3-component relative of mosaicHopper.c: PHASE, RANGE
  OFFSETS and AZIMUTH OFFSETS in a single per-pixel normal-equation system, each row carrying its
  own frame's sigma, but solving for (vx, vy, vz) instead of (vx, vy) -- and deciding PER PIXEL
  which of the two answers to keep.

      PHASE          u = ( sin(psi) cos(gPh),  sin(psi) sin(gPh), -cos(psi) )
                     P = phase * 365.25/(twok nDays)
      RANGE OFFSET   u = ( sin(psi) cos(g),    sin(psi) sin(g),   -cos(psi) )
                     P = dr    * 365.25/nDays
      AZIMUTH OFFSET u = (-sin(g),             cos(g),             0        )
                     P = da    * 365.25/nDays

  NOTE what is NOT here that IS in mosaicHopper: no slope in the rows, and no 1/sin(psi) in the
  data or the sigma.  Both are deliberate and both are load-bearing -- see PROJECTION below.

  THE 2D/3D TOGGLE
  ----------------
  The surface-parallel (2D) and unconstrained (3D) solutions are NOT two different accumulations.
  With C = [[1,0],[0,1],[sx,sy]] -- the surface-parallel constraint vz = s.vh written as a map
  from the 2D subspace into 3-space -- the identity

      N2 = C^T N3 C        b2 = C^T b3

  is EXACT (derivation: Documents/twoDthreeDProjection.md; verified against -speckleTrackJoint at
  median 1.1e-5 m/yr on real data).  So one pass over the data yields BOTH systems, and the choice
  is deferred to solve time -- which is the only place it can be made, since solvability depends
  on lambdaMin of the ASSEMBLED matrix and that is not known until every track has contributed.

  The two sin(psi) factors -- one taken out of the data into the row, one appearing in the weight
  as w2D = sin^2(psi) w3D -- are precisely what make the identity hold.  This is why the 3D form
  must NOT carry 1/sin(psi) the way the 2D form does.

  Per pixel:
      lam = lambdaMin(N3)
      3D  if lam > 0 && nObs >= 3 && nRange >= 2
             && (1/sqrt(lam)) * sqrt(nRange) <= hopper3DMaxSigma
      2D  otherwise  (projected; NO pixel is lost to the toggle)

  nRange counts PHASE + RANGE rows, i.e. rows with vertical sensitivity.  Azimuth rows have
  u_z == 0 and cannot inform the vertical, so normalising by total nObs would credit the gate with
  measurements that say nothing about what it is deciding.

  Unlike mosaicTrue3D there is NO hard requirement for an azimuth row.  NISAR-left plus S1-right
  give antiparallel look directions that make a range-only 3D solve well conditioned (measured:
  33,812/33,812 pixels solved with zero azimuth rows).  lambdaMin already tests the geometry; a
  row-type requirement second-guesses it.

  hopper3DMaxSigma defaults to 35 and that default is ARBITRARY.  N3 is not on the same scale as
  N2 -- the 3D rows are unit vectors while w3D = w2D/sin^2(psi) -- and lambdaMin(N3) is dominated
  by the vertical direction, which is exactly what the toggle is deciding.  The ".mode" band
  exists so the split can be swept and measured rather than assumed.

  GEOMETRY CONVENTION
  -------------------
  xyAngle = atan2(-y,-x) and computeB() for the slope, matching mosaicHopper.c and
  mosaicTrue3DPhase.c.  NOT mosaicTrue3D.c's computeXYangle() + xyGetZandSlope(), which
  re-projects through the DEM's own rot/stdLat and clamps slope at 0.12 with a psScale correction
  instead of 0.25 with none.  That choice is what lets the reduction test against mosaicHopper be
  meaningful; it also means the test against mosaicTrue3D can only be a distributional one.

  Row suppression is per observable, never per image -- a frame with no phase still contributes
  its offsets.  (mosaicHopper got this wrong until 2026-08-30 and silently built S1 mosaics from
  57% of the available offsets.)

  Derived from mosaicHopper.c for the geometry, threading, squint, tide, SMB and buffer handling;
  the 3-component rows and the lambdaMin3/invSym3 helpers come from mosaicTrue3D.c.
  See Documents/hopper3DPlan.md.
*/
/*
  ACCUMULATOR PRECISION -- selectable, one line.

      #define HOPPER3DACC double   (default)
      #define HOPPER3DACC float

  double by default: a 3x3 determinant cancels harder than a 2x2, the vertical column carries the
  least directional diversity of the three, and SddAcc in particular is a difference of large
  near-equal sums -- in float it collapsed chi2 to zero throughout mosaicTrue3DPhase's output.

  float halves the accumulator footprint (80 -> 40 B/px, i.e. 176 -> 136 B/px for the module) and
  makes the forced-2D mode bit-comparable with mosaicHopper, which accumulates in float.  Reduction
  test A measures the difference directly: at double it is median 1.2e-05 m/yr against mosaicHopper;
  at float it should fall to ~0.

  Deliberately a COMPILE-TIME switch rather than a CLI flag, matching TRUE3DACC in mosaicTrue3D.c:
  a runtime choice would need either two accumulation code paths (which then drift apart) or a
  branch in the innermost loop.  Flip it here and rebuild.
*/

#define HOPPER3DACC double

static void setBufferHop3(inputImageStructure *inputImage, float *buf);
static double computePhiZM3dHop3(double *thetaD, double z, double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
							 double Range, double Re, double ReH, double ReHfixed, double thetaC, double thetaCfixedReH, double *phaseError);
static double computePhiFlatEarthM3dHop3(double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
									 double Range, double Re, double ReHfixed, double thetaCfixedReH, double *phaseError);
static float **mallocZeroImageHop3(int32_t nr, int32_t nc);
static HOPPER3DACC **mallocZeroDoubleHop3(int32_t nr, int32_t nc);
static void freeImageHop3(float **image, int32_t nr);
static void freeDoubleHop3(HOPPER3DACC **image, int32_t nr);
/*  Shared with mosaicTrue3D.c -- closed-form symmetric-3x3 smallest eigenvalue and inverse.
    Declared non-static there precisely so the 3D solvers do not each carry a copy. */
double lambdaMin3(double a11, double a12, double a13, double a22, double a23, double a33);
int32_t invSym3(double a11, double a12, double a13, double a22, double a23, double a33,
				double Cinv[3][3], double *det);

void mosaicHopper3D(inputImageStructure *ascImages, inputImageStructure *descImages,
					   vhParams *ascParams, vhParams *descParams, xyDEM *dem,
					   outputImageStructure *outputImage, float fl, int32_t no3d, float timeThreshPhase,
					   referenceVelocity *refVel)
{
	extern double hopper3DMaxSigma;
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
	int32_t hasPhase, framePhaseOK, useRangeRow, useAzRow;
	float azSLPixSize, rSLPixSize;
	double dumRe, dumReH, dumThetaC;
	float geoTolerance = 1e-3;						  /* as make3DMosaic: a few metres is plenty */
	float **vXimage, **vYimage, **vZimage, **errorX, **errorY;
	float **vxTmp, **vyTmp, **vzTmp, **fScale, **sxTmp, **syTmp;
	float **scaleX, **scaleY, **scaleZ;
	/*  Normal-equation accumulators, one plane each over the whole output grid.  ALL DOUBLE,
	    unlike the 2D hopper's float planes.  A 3x3 determinant cancels harder than a 2x2 and the
	    vertical column carries the least directional diversity, so the margin float leaves is not
	    enough.  SddAcc in particular MUST be double: chi2 = Sdd - v.b is a difference of large
	    near-equal sums, and in float it collapsed to zero throughout the core of
	    mosaicTrue3DPhase's output. */
	HOPPER3DACC **Nxx, **Nxy, **Nxz, **Nyy, **Nyz, **Nzz;
	HOPPER3DACC **bxAcc, **byAcc, **bzAcc;
	HOPPER3DACC **SddAcc;										/* sum w d^2, for chi-square at solve time */
	float **nObs, **nObsPh, **nObsRg, **nObsAz;				 /* contributing measurements per pixel -- diagnostic */
	/*  -gateNEff: sum(w) and sum(w^2) per pixel, for the weighted effective row count
	    (sum w)^2/sum(w^2) used in place of the raw count in the rejection gate. */
	float **gW, **gW2;
	float **chi2Plane;
	float **vZ3D, **errZ, **modeBand;			 /* solved vertical, its VARIANCE, and 3-vs-2 per pixel */
	float **sumW, **sumWT;		 /* timeOverlapFlag only; NULL otherwise */
	float dum1, dum2;
	float azimuthMin, azimuthMax;
	int32_t iMin, iMax, jMin, jMax;					 /* current image's region */
	int32_t iMinAll, iMaxAll, jMinAll, jMaxAll;		 /* union over all contributing images */
	int32_t ii, nTotal, nUsed;
	int32_t i;
	int64_t nSolved, nRejCond, nRejNoData, nSolved3D, nSolved2D, nRejClip;
	double nObsSum, nObsMax;
	struct timeval funcStart, funcEnd;
	gettimeofday(&funcStart, NULL);

	fprintf(outputImage->fpLog, ";\n; Entering mosaicHopper3D(.c)\n");
	if (no3d == TRUE)
	{
		fprintf(outputImage->fpLog, ";\n; no3d flag set, returning;\n; Returning from mosaicHopper3D(.c)\n");
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
	Nxx = mallocZeroDoubleHop3(outputImage->ySize, outputImage->xSize);
	Nxy = mallocZeroDoubleHop3(outputImage->ySize, outputImage->xSize);
	Nxz = mallocZeroDoubleHop3(outputImage->ySize, outputImage->xSize);
	Nyy = mallocZeroDoubleHop3(outputImage->ySize, outputImage->xSize);
	Nyz = mallocZeroDoubleHop3(outputImage->ySize, outputImage->xSize);
	Nzz = mallocZeroDoubleHop3(outputImage->ySize, outputImage->xSize);
	bxAcc = mallocZeroDoubleHop3(outputImage->ySize, outputImage->xSize);
	byAcc = mallocZeroDoubleHop3(outputImage->ySize, outputImage->xSize);
	bzAcc = mallocZeroDoubleHop3(outputImage->ySize, outputImage->xSize);
	SddAcc = mallocZeroDoubleHop3(outputImage->ySize, outputImage->xSize);
	nObs = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
	nObsPh = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
	nObsRg = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
	nObsAz = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
	gW = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
	gW2 = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
	chi2Plane = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
	vZ3D = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
	errZ = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
	modeBand = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
	sumW = NULL;
	sumWT = NULL;
	if (outputImage->timeOverlapFlag == TRUE)
	{
		sumW = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
		sumWT = mallocZeroImageHop3(outputImage->ySize, outputImage->xSize);
	}
	allImages = ascImages; /* consolodateLists() has already appended desc to asc */
	params = ascParams;
	nTotal = 0;
	for (phaseImage = allImages; phaseImage != NULL; phaseImage = phaseImage->next)
	{
		nTotal++;
	}
	fprintf(stderr, "mosaicHopper: nTotal Images %i\n", nTotal);
	tCenter = (outputImage->jd1 + outputImage->jd2 + 1.) * 0.5; /* as make3DMosaic */
	iMinAll = outputImage->ySize;
	jMinAll = outputImage->xSize;
	iMaxAll = 0;
	jMaxAll = 0;
	int32_t timeThreshPhaseWarned = FALSE;   /* -timePhaseThresh is a pair gate; warn once, do not filter */
	nUsed = 0;
	int nthreads = omp_get_max_threads();
	inputImageStructure *localImgs = (inputImageStructure *)malloc((size_t)nthreads * sizeof(inputImageStructure));
	if (localImgs == NULL)
	{
		error("mosaicHopper: malloc failed for per-thread image copies\n");
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
		/*  Partial-temporal-overlap weights are carried in the input file (setupquarters writes
		    a fraction < 1 for frames only partly inside the mosaic window).  They are honoured
		    ONLY when -timeOverlap is set; without it they would be silently discarded and the
		    mosaic would weight a 4%-overlap frame the same as a full one.  speckleTrackMosaic
		    and speckleTrackMosaicJoint have always treated that as fatal -- match them.
		    Added 2026-08-30 after an S1 run (weights 0.042..1.0) silently ignored them here
		    while the 2D path correctly refused. */
		if (fabs(phaseImage->weight - 1.0) > 0.01 && outputImage->timeOverlapFlag == FALSE)
		{
			error("non unity weight, but overlap flag not set\n");
		}
		/*  These two gates are PHASE-ROW gates, not image gates.  They used to `continue`, which
		    threw away the frame's range and azimuth offsets as well -- contradicting the header
		    above and, on a Sentinel-1 archive where 43% of frames are "nophase", silently
		    building the mosaic from 57% of the available offsets.  That is why hop_s1 landed
		    within 0.002 m/yr of the phase-only product.  Suppress the phase row and carry on;
		    the `hasPhase == FALSE && no offset rows` test further down is what drops an image
		    that genuinely has nothing to contribute. */
		framePhaseOK = (strstr(phaseImage->file, "nophase") == NULL);
		/*  -timePhaseThresh DOES NOT APPLY TO THE HOPPER.  See the identical note in
		    mosaicHopper.c: it is a PAIR gate limiting the separation between the two images of
		    a crossing pair, which keeps the pairwise solvers from forming a combinatorial number
		    of pairs that are then averaged for no gain.  The hopper accumulates each image once
		    and forms no pairs, so the temporal extent is set by the DATE RANGE (-date1/-date2).
		    Until 2026-09-07 this reused the parameter as a per-image window against the mosaic
		    centre date -- a different quantity, harmless only because the templates carry 10000.
		    Now: warn once that it does not apply, and do not filter. */
		if (timeThreshPhaseWarned == FALSE && timeThreshPhase > 0.0)
		{
			timeThreshPhaseWarned = TRUE;
			fprintf(stderr, "; note: -timePhaseThresh (%.0f) does not apply to the hopper -- it is a "
					"pair gate, and the hopper forms no pairs.  Temporal extent is set by "
					"-date1/-date2.\n", timeThreshPhase);
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
		setBufferHop3(phaseImage, AImageBuffer);
		/*  ONE call only.  dum1/dum2 ARE azSLPixSize/rSLPixSize -- calling setupGeoConversions
	    a second time for the offsets corrupts the conversion state the phase path is
	    using, which silently zeroed every row (symptom: 223 images contribute, nObs 0). */
		cp = setupGeoConversions(phaseImage, &azSLPixSize, &rSLPixSize, &Re, &ReH, &thetaC, &ReHfixed, &thetaCfixedReH);
		/*  DO NOT relax the geocoding tolerance here.  make3DMosaicJoint sets 1e-3 as a SPEED
		    optimisation for the phase path ("a few metres is plenty"); the default from
		    initLLtoImage is 1e-6.  Relaxing it makes the OFFSETS rows less precisely geolocated
		    than make3DOffsetsJoint's -- measured at 0.045 m/yr median, which is exactly what the
		    reduction test against that solver caught.  Making it conditional on hasPhase is worse
		    still: the tolerance would then differ between frames within one mosaic.  So the
		    hopper keeps the tight default for every row, and pays the extra Newton iterations. */
		getAzimuthBoundsForXYBox(iMin, iMax, jMin, jMax, phaseImage, outputImage, &azimuthMin, &azimuthMax);
		if (azimuthMin == 0.0f && azimuthMax == 0.0f)
		{
			continue;
		}
		hasPhase = (phaseImage->file != NULL) && (framePhaseOK == TRUE) && (noPhaseRows == FALSE);
		/*  GUARDED READ.  A corrupt or truncated input for THIS product (a half-written
		    rBaseline, a bad .dat) used to exit() the whole run from deep inside error().  With
		    thousands of products that is the wrong trade: while errorRecoveryJmp is armed,
		    error() longjmps back here, the product is recorded and skipped, and the list is
		    printed at the end of the run (printFailedProducts).  The guard is disarmed as soon
		    as the reads are done -- it MUST be NULL before the parallel region below. */
		jmp_buf readJmp;
		if (setjmp(readJmp) != 0)
		{
			errorRecoveryJmp = NULL;
			recordFailedProduct(params->offsets.rFile != NULL ? params->offsets.rFile : phaseImage->file,
								errorRecoveryMsg);
			continue;
		}
		errorRecoveryJmp = &readJmp;
		if (hasPhase == TRUE)
		{
			getMosaicInputImage(phaseImage, azimuthMin, azimuthMax);
			getIonospherePhaseImage(phaseImage, AIonBuffer, azimuthMin, azimuthMax);
		}
		phaseImage->memChan = MEM1;
		/*  OFFSETS for the SAME frame.  Read into their own buffers alongside the phase, so one
		    pass over the image yields all three observables.  Row suppression below is per
		    observable -- a frame missing any one still contributes the others. */
		useRangeRow = FALSE; useAzRow = FALSE;
		if (params->offsets.rFile != NULL)
		{
			readOffsetDataAndParams(&(params->offsets), azimuthMin, azimuthMax);
			if (params->offsets.deltaB != DELTABNONE)
			{
				double bnS, bpS;
				svAzOffset(phaseImage, &(params->offsets), 0.0, 0.0);
				svInterpBnBp(phaseImage, &(params->offsets), 0.0, &bnS, &bpS);
			}
			useRangeRow = (params->offsets.sigmaRresidual >= 0) && (noRangeRows == FALSE);
			useAzRow = (params->offsets.file != NULL && params->offsets.da != NULL &&
						params->offsets.sigmaAresidual >= 0 &&
						params->offsets.sigmaAresidual <= outputImage->sigmaAThresh &&
						/*  -sigmaAThreshVel: the azparams residual is METRES over the pair
						    interval; scale it to m/yr so one threshold means the same thing for
						    12- and 24-day pairs.  Inert when < 0. */
						(sigmaAThreshVel < 0.0 ||
						 params->offsets.sigmaAresidual * 365.25 / (double)params->nDays
							 <= sigmaAThreshVel) &&
						noAzimuthRows == FALSE);
		}
		errorRecoveryJmp = NULL;
		if (hasPhase == FALSE && useRangeRow == FALSE && useAzRow == FALSE)
		{
			continue;
		}
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
			int32_t pixHasPhase;
			double lat, lon, x, y, zWGS84, zSp;
			double range, azimuth, Range, myReH, theta, thetaD, psi;
			double phase, ionPhase, phiZ, phaseError;
			double hAngle, xyAngle, gamma, savedLastTime;
			double sinPsi, cosPsi;
			double urx, ury, urz, uphx, uphy;	/* 3-component LOS rows: offsets and phase */
			double gammaPh, hAngleSq, scale, P, sigma, w, dzdtSubmergence;
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
					/*  ONLY read the phase image when this frame actually has one.  For a
					    "nophase" frame getMosaicInputImage() was never called, so the buffer
					    holds whatever the previous frame left behind -- or nothing at all on the
					    first such frame.  Reading it segfaults. */
					phase = (double)(-LARGEINT);
					if (hasPhase == TRUE)
					{
						interpPhaseImage(myImg, range, azimuth, &phase);
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
					/*  A pixel with no valid PHASE can still have valid OFFSETS, so this must
					    not 'continue' -- it only disables the phase row.  (Getting this wrong
					    is what made the first build accumulate nothing.) */
					pixHasPhase = (hasPhase == TRUE) && (phase > -LARGEINT);
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
					/*  Skipped entirely for a pixel with no phase: these read the BASELINE
					    solution, which a "nophase" frame need not have, and the result would be
					    unused anyway.

					    computePhiZM3dHop3 WRITES its thetaD argument, so it gets a private copy
					    -- every offsets solver (make3DOffsets, make3DOffsetsJoint,
					    speckleTrackMosaic) feeds interpRangeOffsetInMeters and computeSig2Base
					    the thetaD that geometryInfo() produced, and the offsets rows below must
					    match them.  No effect on the flat-earth (NISAR) path. */
					phiZ = 0.0;
					phaseError = 0.0;
					if (pixHasPhase == TRUE)
					{
						if (params->applyFlatEarth)
						{
							phiZ = computePhiFlatEarthM3dHop3(azimuth, params, myImg, Range, Re, ReHfixed, thetaCfixedReH, &phaseError);
						}
						else
						{
							double thetaDPhase = thetaD;
							phiZ = computePhiZM3dHop3(&thetaDPhase, zSp, azimuth, params, myImg, Range, Re, myReH, ReHfixed, thetaC, thetaCfixedReH, &phaseError);
						}
					}
					/*  ---- SHARED GEOMETRY (was deleted by an over-wide edit; restored) -----
					    hAngle/gamma and the slope matrix are common to all three observables, so
					    they are computed ONCE per pixel here.  computeHeading() calls
					    llToImageNew() at lat+-dlat and clobbers the Newton warm-start cache, so
					    save/restore lastTime exactly as computeA() does.  myImg is the per-thread
					    copy, so this is thread-safe. */
					savedLastTime = myImg->lastTime;
					hAngle = computeHeading(lat, lon, 0.0, myImg, &(myImg->cpAll));
					hAngleSq = hAngle;
					if (useSquint && myImg->hasSquintPolynomial)
					{
						double sqRange, sqAzimuth;
						llToImageNew(lat, lon, 0.0, &sqRange, &sqAzimuth, myImg);
						hAngleSq += evaluateSquint(myImg, sqRange, sqAzimuth) * DTOR;
					}
					myImg->lastTime = savedLastTime;
					/* xyAngle = PI/2 + meridian convergence.  For polar stereographic this returns
					   the original atan2(-y,-x), plus PI in the south, bit for bit; for UTM it uses
					   the true convergence.  Algebraically they are the same formula. */
					xyAngle = grimpXYAngle(lat, lon, x, y, &(outputImage->proj));
					/*  TWO headings, deliberately.  Squint applies to the PHASE row only.
					    make3DOffsets/make3DOffsetsJoint pass FALSE to computeA unconditionally
					    because the zero-Doppler condition makes range/azimuth OFFSETS
					    self-consistent regardless of squint (root CLAUDE.md, "Squint").  Using
					    the squinted heading for the offsets rows rotates them by ~1.6 deg, worth
					    ~1.8 m/yr at 60 m/yr -- which is exactly what the reduction test against
					    make3DOffsetsJoint caught. */
					gammaPh = xyAngle - hAngleSq;
					gamma = xyAngle - hAngle;
					/*  ---- THE THREE-COMPONENT LOS ROWS ----------------------------------
					    NO SLOPE, and NO 1/sin(psi) folded into the data or the sigma below.
					    The slope enters only at solve time, through C in pass 2; taking the
					    sin(psi) out of the data and putting it in the row is what makes
					    N2 = C^T N3 C hold exactly (Documents/twoDthreeDProjection.md sec 3).

					    The vertical component is -cos(psi): a target moving UP shortens the
					    slant range.  Had it been +cos(psi) the projected row would read
					    cos(gamma) + cot(psi) sx and the reduction test would fail by ~3% of
					    velocity rather than by 1e-5.

					    Two headings, deliberately: squint applies to the PHASE row only.
					    make3DOffsets passes FALSE to computeA unconditionally because the
					    zero-Doppler condition makes range/azimuth OFFSETS self-consistent
					    regardless of squint (root CLAUDE.md, "Squint").  Using the squinted
					    heading for the offsets rows rotates them ~1.6 deg -- worth ~1.8 m/yr
					    at 60 m/yr, which is what the 2D hopper's reduction test caught. */
					sinPsi = sin(psi);
					cosPsi = cos(psi);
					urx = sinPsi * cos(gamma);
					ury = sinPsi * sin(gamma);
					urz = -cosPsi;
					uphx = sinPsi * cos(gammaPh);
					uphy = sinPsi * sin(gammaPh);
					/* uphz == urz: squint rotates the row about the vertical, so the vertical
					   component is common to both LOS rows. */
					if (pixHasPhase == TRUE)
					{
						phase = phase - phiZ;
						if (myImg->ionospherePhase != NULL)
						{
							interpIonPhaseImage(myImg, range, azimuth, &ionPhase);
							/*  Phase screen shares the phase image's sign convention -> SUBTRACT.
							    Opposite to the offsets correction below, which is pre-negated and
							    ADDED.  Do not "fix" for consistency. */
							if (ionPhase > -0.98 * LARGEINT) { phase -= ionPhase; }
						}
						if (sMask == SHELF)
						{
							interpTideError(&phaseError, myImg, params, x, y, psi, twok);
							phase -= -myImg->tideCorrection * cos(psi) * twok * (double)params->nDays / 365.25;
						}
						if (vCorrect != NULL)
						{
							dzdtSubmergence = interpVCorrect(x, y, vCorrect);
							phase -= -dzdtSubmergence * cos(psi) * twok * (double)params->nDays / 365.25;
						}
						/*  NO sin(psi) here -- it lives in the row.  See the row comment above. */
						scale = 365.25 / (twok * params->nDays);
						P = phase * scale;
						sigma = phaseError * scale;
						if (sigma > 0.0 && fabs(phase) < 0.9 * LARGEINT)
						{
							w = (sumW != NULL) ? phaseImage->weight / (sigma * sigma) : 1.0 / (sigma * sigma);
							if (w > 0.0)
							{
								Nxx[iRow][jj] += w * uphx * uphx;
								Nxy[iRow][jj] += w * uphx * uphy;
								Nxz[iRow][jj] += w * uphx * urz;
								Nyy[iRow][jj] += w * uphy * uphy;
								Nyz[iRow][jj] += w * uphy * urz;
								Nzz[iRow][jj] += w * urz * urz;
								bxAcc[iRow][jj] += w * uphx * P;
								byAcc[iRow][jj] += w * uphy * P;
								bzAcc[iRow][jj] += w * urz * P;
								SddAcc[iRow][jj] += w * P * P;
								nObs[iRow][jj] += 1.0f;
								nObsPh[iRow][jj] += 1.0f;
								gW[iRow][jj] += (float)w;
								gW2[iRow][jj] += (float)(w * w);
								obsDumpRecord(iRow, jj, "phase", phaseImage->file,
											  phaseImage->julDay, (double)params->nDays,
											  uphx, uphy, urz, P, sigma, w);
								/*  dT accumulated INSIDE the weight block, not once per pixel at
								    the bottom.  mosaicHopper.c:506-510 does the latter, which
								    reads a stale `w` from the previous pixel whenever this
								    pixel's rows were all rejected.  Follow mosaicTrue3D.c. */
								if (sumW != NULL)
								{
									sumW[iRow][jj] += (float)w;
									sumWT[iRow][jj] += (float)(w * (tOffCenter - tCenter));
								}
							}
						}
					}
					/*  ---- ROW 2: RANGE OFFSET (same LOS direction, own sigma) ---------- */
					if (useRangeRow == TRUE)
					{
						double demErrR, sig2BaseR, sigmaR, dr, scaleR;
						float ionCorr;
						dr = interpRangeOffsetInMeters(range, azimuth, &(params->offsets), myImg,
													   Range, thetaD, rSLPixSize, theta, &demErrR);
						if (params->offsets.rOffCorrection.rangeOffsetCorrection != NULL)
						{
							double rSLC, aSLC;
							computeSLCFromMLCoords(myImg, range, azimuth, &rSLC, &aSLC);
							ionCorr = interpolateOffsetIonCorrectionInPixels(&(params->offsets.rOffCorrection),
																			 rSLC, aSLC, -LARGEINT, 0.0);
						}
						else { ionCorr = 0.0; }
						if (dr > (-LARGEINT + 1) && fabs(dr) < 13.0E4)
						{
							if (ionCorr > -0.98 * LARGEINT) { dr += ionCorr * rSLPixSize; }
							sigmaR = interpRangeSigma(range, azimuth, &(params->offsets), myImg, Range, thetaD, rSLPixSize);
							sig2BaseR = computeSig2Base(sin(thetaD), cos(thetaD), azimuth, myImg, &(params->offsets));
							sigmaR = sqrt(sigmaR * sigmaR + demErrR * demErrR + sig2BaseR +
										  rangeAccuracyVar(&(params->offsets)));
							dr -= shelfMaskCorrection(myImg, params, sMask, x, y, psi, &sigmaR);
							if (vCorrect != NULL)
							{
								dzdtSubmergence = interpVCorrect(x, y, vCorrect);
								dr -= -dzdtSubmergence * cos(psi) * (double)params->nDays / 365.25;
							}
							/*  Again NO sin(psi): it is in the row. */
							scaleR = 365.25 / (double)params->nDays;
							P = dr * scaleR;
							sigma = sigmaR * scaleR;
							if (sigma > 0.0)
							{
								w = (sumW != NULL) ? phaseImage->weight / (sigma * sigma) : 1.0 / (sigma * sigma);
								if (w > 0.0)
								{
									Nxx[iRow][jj] += w * urx * urx;
									Nxy[iRow][jj] += w * urx * ury;
									Nxz[iRow][jj] += w * urx * urz;
									Nyy[iRow][jj] += w * ury * ury;
									Nyz[iRow][jj] += w * ury * urz;
									Nzz[iRow][jj] += w * urz * urz;
									bxAcc[iRow][jj] += w * urx * P;
									byAcc[iRow][jj] += w * ury * P;
									bzAcc[iRow][jj] += w * urz * P;
									SddAcc[iRow][jj] += w * P * P;
									nObs[iRow][jj] += 1.0f;
									nObsRg[iRow][jj] += 1.0f;
									gW[iRow][jj] += (float)w;
									gW2[iRow][jj] += (float)(w * w);
									obsDumpRecord(iRow, jj, "range", phaseImage->file,
												  phaseImage->julDay, (double)params->nDays,
												  urx, ury, urz, P, sigma, w);
									if (sumW != NULL)
									{
										sumW[iRow][jj] += (float)w;
										sumWT[iRow][jj] += (float)(w * (tOffCenter - tCenter));
									}
								}
							}
						}
					}
					/*  ---- ROW 3: AZIMUTH OFFSET (along-track, NO slope term) ----------- */
					if (useAzRow == TRUE)
					{
						double sigmaA, sig2Off, da, kAz, aax, aay;
						da = interpAzOffset(range, azimuth, &(params->offsets), myImg, Range, theta, azSLPixSize);
/* azimuth ionosphere: no-op unless the az fit recorded one and -useAzIonosphere is set */
if (da > -0.98 * LARGEINT)
	da += azIonCorrectionMeters(&(params->offsets), myImg, range, azimuth, azSLPixSize);
						if (fabs(da) < 10.0e4)
						{
							sigmaA = interpAzSigma(range, azimuth, &(params->offsets), myImg, Range, theta, azSLPixSize);
							sig2Off = computeSig2AzParam(sin(theta), cos(theta), azimuth, Range, myImg, &(params->offsets));
							sigmaA = sqrt(sigmaA * sigmaA + sig2Off + azimuthAccuracyVar(&(params->offsets)));
							kAz = 365.25 / (double)params->nDays;
							/*  Along-track row: NO vertical component.  Azimuth offsets sense
							    horizontal motion only, so u_z == 0 and the row contributes
							    nothing to Nxz/Nyz/Nzz/bz.  That is also why C^T u_a is the same
							    in 2D and 3D -- the projection leaves this row untouched
							    (Documents/twoDthreeDProjection.md sec 3.5) -- and why the
							    gate below normalises by nRange, not nObs. */
							aax = -sin(gamma);
							aay = cos(gamma);
							P = (double)da * kAz;
							sigma = sigmaA * kAz;
							if (sigma > 0.0)
							{
								w = (sumW != NULL) ? phaseImage->weight / (sigma * sigma) : 1.0 / (sigma * sigma);
								if (w > 0.0)
								{
									Nxx[iRow][jj] += w * aax * aax;
									Nxy[iRow][jj] += w * aax * aay;
									Nyy[iRow][jj] += w * aay * aay;
									bxAcc[iRow][jj] += w * aax * P;
									byAcc[iRow][jj] += w * aay * P;
									SddAcc[iRow][jj] += w * P * P;
									nObs[iRow][jj] += 1.0f;
									nObsAz[iRow][jj] += 1.0f;
									gW[iRow][jj] += (float)w;
									gW2[iRow][jj] += (float)(w * w);
									/*  az component is identically 0 for azimuth rows -- see the
									    row comment above; it is dumped as 0 so the Python side can
									    form the same 3x3 without special-casing. */
									obsDumpRecord(iRow, jj, "azimuth", phaseImage->file,
												  phaseImage->julDay, (double)params->nDays,
												  aax, aay, 0.0, P, sigma, w);
									if (sumW != NULL)
									{
										sumW[iRow][jj] += (float)w;
										sumWT[iRow][jj] += (float)(w * (tOffCenter - tCenter));
									}
								}
							}
						}
					}
#pragma omp atomic write
					phaseImage->used = TRUE;
				} /* j loop */
			}	  /* i loop */
		}		  /* End omp parallel */
	}			  /* End image loop */
	fprintf(stderr, "mosaicHopper: %i of %i images contributed\n", nUsed, nTotal);
	/*
	  ****************** PASS 2: solve once per pixel *****************************************
	*/
	nSolved = 0;
	nRejCond = 0;
	nRejNoData = 0;
	nSolved3D = 0;
	nSolved2D = 0;
	nRejClip = 0;
	nObsSum = 0.0;
	nObsMax = 0.0;
	if (nUsed > 0 && iMaxAll > iMinAll && jMaxAll > jMinAll)
	{
#pragma omp parallel for schedule(dynamic, 8) \
	reduction(+ : nSolved, nRejCond, nRejNoData, nSolved3D, nSolved2D, nRejClip, nObsSum) reduction(max : nObsMax)
		for (i = iMinAll; i < iMaxAll; i++)
		{
			double lat, lon, x, y, zWGS84;
			double n11, n12, n13, n22, n23, n33, b1, b2, b3;
			double nxx, nxy, nyy, bx, by, det, trHalf, disc, lambdaMin;
			double vx, vy, vz, scX, scY, deltaOffCenter, gateCap, gateSpd;
			double sx, sy, lam3, det3, Cinv[3][3], nRange, nGate, gateN;
			double B[2][2], dzdx, dzdy;
			int32_t use3D;
			unsigned char sMask2;
			int32_t jj;
			y = (outputImage->originY + i * outputImage->deltaY) * MTOKM;
			for (jj = jMinAll; jj < jMaxAll; jj++)
			{
				if (nObs[i][jj] < 1.5)
				{
					/* 0 or 1 measurements -- two unknowns need two independent looks */
					nRejNoData++;
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					continue;
				}
				n11 = Nxx[i][jj]; n12 = Nxy[i][jj]; n13 = Nxz[i][jj];
				n22 = Nyy[i][jj]; n23 = Nyz[i][jj]; n33 = Nzz[i][jj];
				b1 = bxAcc[i][jj]; b2 = byAcc[i][jj]; b3 = bzAcc[i][jj];
				/*  ---- THE SLOPE, and therefore C --------------------------------------
				    Recomputed here rather than carried in two full-grid planes -- the same
				    trade mosaicHopper.c records, and it is now a correctness requirement too,
				    since the accumulator is deliberately slope-free.  psi is irrelevant
				    because only dzdx/dzdy are used, so pass 1.0 twice.

				    SHELF carve-out: mosaicHopper zeroes the slope coupling on shelf pixels
				    (a large DEM slope over shelf interior is almost always a stale rift or
				    calving front), so C must use sx=sy=0 there or every shelf pixel disagrees
				    with the 2D reference.  Note the reported surface-parallel vz below uses
				    the UNzeroed dzdx/dzdy -- reproducing mosaicHopper.c:607-608 exactly,
				    including its internal inconsistency. */
				x = (outputImage->originX + jj * outputImage->deltaX) * MTOKM;
				xyToLLProj(x, y, &lat, &lon, &(outputImage->proj));
				zWGS84 = getXYHeight(lat, lon, dem, 0.0, ELLIPSOIDAL);
				computeB(x, y, zWGS84, B, &dzdx, &dzdy, 1.0, 1.0, (xyDEM *)dem);
				sMask2 = GROUNDED;
				if (shelfMask != NULL) { sMask2 = getShelfMask(shelfMask, x, y); }
				sx = dzdx;
				sy = dzdy;
				if (sMask2 == SHELF) { sx = 0.0; sy = 0.0; }
				/*  Dumped AFTER the shelf carve-out, so it is the sx/sy the projection
				    actually uses -- not the raw DEM slope. */
				obsDumpPixel(i, jj, sx, sy);
				/*  ---- DECIDE: 3D or its surface-parallel projection --------------------
				    nRange counts rows with vertical sensitivity (phase + range offsets).
				    Azimuth rows have u_z == 0, so normalising by total nObs would credit the
				    gate with measurements that say nothing about the vertical -- the very
				    thing being decided. */
				nRange = (double)nObsPh[i][jj] + (double)nObsRg[i][jj];
				lam3 = lambdaMin3(n11, n12, n13, n22, n23, n33);
				/*  hopper3DMaxSigma sign convention:
				        < 0   force the 2D projection everywhere (reduction test A)
				        == 0  gate off -- 3D wherever N3 is invertible
				        > 0   the n-normalised gate */
				use3D = (hopper3DMaxSigma >= 0.0) &&
						(lam3 > 0.0) && (nObs[i][jj] >= 2.5) && (nRange >= 1.5) &&
						(hopper3DMaxSigma == 0.0 ||
						 (1.0 / sqrt(lam3)) * ((gateAbsolute == TRUE) ? 1.0 : sqrt(nRange))
							 <= hopper3DMaxSigma);
				if (use3D == TRUE)
				{
					if (invSym3(n11, n12, n13, n22, n23, n33, Cinv, &det3) == FALSE)
					{
						use3D = FALSE;
					}
				}
				if (use3D == TRUE)
				{
					/*  Unconstrained 3-component solve.  Cov = N3^-1, so var(vx) = Cinv[0][0]
					    and the mosaic accumulator's 1/variance is its reciprocal. */
					vx = Cinv[0][0] * b1 + Cinv[0][1] * b2 + Cinv[0][2] * b3;
					vy = Cinv[1][0] * b1 + Cinv[1][1] * b2 + Cinv[1][2] * b3;
					vz = Cinv[2][0] * b1 + Cinv[2][1] * b2 + Cinv[2][2] * b3;
					if (!(Cinv[0][0] > 0.0) || !(Cinv[1][1] > 0.0))
					{
						nRejCond++;
						vxTmp[i][jj] = (float)-LARGEINT;
						fScale[i][jj] = 0.0;
						continue;
					}
					scX = 1.0 / Cinv[0][0];
					scY = 1.0 / Cinv[1][1];
					/*  errZ carries the VARIANCE, not the sigma.  mosaic3d.c:942-945 takes the
					    sqrt on output; storing a sigma here would ship its square root. */
					vZ3D[i][jj] = (float)vz;
					errZ[i][jj] = (float)Cinv[2][2];
					modeBand[i][jj] = 3.0f;
					/*  chi2 for the 3D solve: 3 parameters, so dof = n - 3. */
					if (chi2Plane != NULL)
					{
						double chi2 = SddAcc[i][jj] - (vx * b1 + vy * b2 + vz * b3);
						double dof = (double)nObs[i][jj] - 3.0;
						chi2Plane[i][jj] = (dof > 0.0 && chi2 > 0.0) ? (float)(chi2 / dof) : 0.0f;
					}
					nSolved3D++;
				}
				else
				{
					/*  Surface-parallel fallback, by PROJECTION of the same accumulator:
					        N2 = C^T N3 C      b2 = C^T b3
					    with C = [[1,0],[0,1],[sx,sy]].  Written out rather than looped, so the
					    correspondence with the derivation is checkable by eye. */
					nxx = n11 + 2.0 * sx * n13 + sx * sx * n33;
					nxy = n12 + sx * n23 + sy * n13 + sx * sy * n33;
					nyy = n22 + 2.0 * sy * n23 + sy * sy * n33;
					bx = b1 + sx * b3;
					by = b2 + sy * b3;
					det = nxx * nyy - nxy * nxy;
					trHalf = 0.5 * (nxx + nyy);
					disc = sqrt((0.5 * (nxx - nyy)) * (0.5 * (nxx - nyy)) + nxy * nxy);
					lambdaMin = trHalf - disc;
					/*  mosaicHopper's gate, verbatim -- jointMaxSigma and nObs, NOT
					    hopper3DMaxSigma and nRange.  This branch has to reproduce the 2D
					    solver, and a different gate would silently change which pixels
					    survive. */
					/*  -gateNEff replaces the raw row count with the WEIGHTED effective count
					    (sum w)^2/sum(w^2).  The gate reads sigmaWorst*sqrt(n) as the
					    per-measurement sigma, which holds only when the rows carry comparable
					    weight.  They do not: azimuth rows sit ~4 orders of magnitude below phase
					    (measured 7e-6 vs 4e-2), contribute nothing to lambdaMin, yet pad n and
					    inflate the gate ~1.7x -- enough to reject fast-ice pixels whose solutions
					    agree with GPS to 0.2-4%.  Identical to the raw count when weights are
					    homogeneous.  Documents/biasInvestigation/currentState.md. */
					nGate = (double)nObs[i][jj];
					if (gateNEff == TRUE && gW2[i][jj] > 0.0)
					{
						nGate = ((double)gW[i][jj] * (double)gW[i][jj]) / (double)gW2[i][jj];
					}
					/*  -gateAbsolute drops the sqrt(n) factor: the test becomes a plain cap on the
					    worst-direction formal sigma.  See mosaic3d.c -- sigma_worst carries the
					    geometric dilution as well as the noise, so the n-normalised form rejects
					    geometry-limited pixels no matter how good their data is. */
					gateN = (gateAbsolute == TRUE) ? 1.0 : sqrt(nGate);
					if (!(det > 0.0) || !(lambdaMin > 0.0))
					{
						nRejCond++;
						vxTmp[i][jj] = (float)-LARGEINT;
						fScale[i][jj] = 0.0;
						continue;
					}
					vx = (nyy * bx - nxy * by) / det;
					vy = (nxx * by - nxy * bx) / det;
					/*  Sigma cap, applied AFTER the solve so -gateSpeedFrac can scale it with the
					    speed: effective cap = max(jointMaxSigma, gateSpeedFrac*|v|).  With
					    gateSpeedFrac 0 this is identical to the pre-solve test it
					    replaces -- the solve has no side effects, and det/lambdaMin are already
					    validated above. */
					gateCap = jointMaxSigma;
					if (gateSpeedFrac > 0.0)
					{
						gateSpd = sqrt(vx * vx + vy * vy);
						if (gateSpeedFrac * gateSpd > gateCap)
						{
							gateCap = gateSpeedFrac * gateSpd;
						}
					}
					if (gateCap > 0.0 && !((1.0 / sqrt(lambdaMin)) * gateN <= gateCap))
					{
						nRejCond++;
						vxTmp[i][jj] = (float)-LARGEINT;
						fScale[i][jj] = 0.0;
						continue;
					}
					scX = det / nyy;
					scY = det / nxx;
					/*  Reported vz is the surface-parallel PREDICTION from the UNzeroed slope,
					    matching mosaicHopper.c:607-608.  No formal vertical error exists for a
					    constrained solve, so errZ stays 0 (written as no-data). */
					vz = vx * dzdx + vy * dzdy;
					vZ3D[i][jj] = (float)vz;
					errZ[i][jj] = 0.0f;
					modeBand[i][jj] = 2.0f;
					if (chi2Plane != NULL)
					{
						double chi2 = SddAcc[i][jj] - (vx * bx + vy * by);
						double dof = (double)nObs[i][jj] - 2.0;
						chi2Plane[i][jj] = (dof > 0.0 && chi2 > 0.0) ? (float)(chi2 / dof) : 0.0f;
					}
					nSolved2D++;
				}
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
				/*  vz, the slope and chi2 are all set inside the two branches above -- the 3D
				    branch has a SOLVED vertical, the 2D branch a slope-derived one, and their
				    degrees of freedom differ (n-3 vs n-2).  Nothing further to compute here. */
				/*  -maxChi2: blunder screen, applied to whichever branch ran.  chi2 asks whether
				    the rows agree with EACH OTHER -- the failure the formal sigma cannot see.
				    Placed after both branches so the .chi2 band still records the value. */
				if (maxChi2 > 0.0 && chi2Plane != NULL && (double)chi2Plane[i][jj] > maxChi2)
				{
					nRejCond++;
					/*  The branches above already counted this pixel as solved; undo that so
					    the "solved N (3D a / 2D b)" line still adds up. */
					if (use3D == TRUE) { nSolved3D--; } else { nSolved2D--; }
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					continue;
				}
				/*  Reference-velocity clip (-clipVel/-clipThresh), applied to whichever branch
				    ran, for the same reason as -maxChi2 above: it needs the final (vx, vy).
				    Inert unless refVel->clipFlag is set; a pixel with no reference value is
				    kept.  Decrements the branch counter so the "solved N (3D a / 2D b)" line
				    still adds up, exactly as the chi2 screen does. */
				if (clipVelChecked(x, y, vx, vy, refVel) == TRUE)
				{
					nRejClip++;
					if (use3D == TRUE) { nSolved3D--; } else { nSolved2D--; }
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					continue;
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
	fprintf(stderr, "mosaicHopper3D: solved %ld (3D %ld / 2D %ld)  rejCond %ld  rejClip %ld  rejNoData %ld  meanNobs %.2f  maxNobs %.0f\n",
			(long)nSolved, (long)nSolved3D, (long)nSolved2D, (long)nRejCond, (long)nRejClip, (long)nRejNoData,
			(nSolved > 0) ? nObsSum / (double)nSolved : 0.0, nObsMax);
	fprintf(outputImage->fpLog, "; mosaicHopper3D nImagesUsed : %i of %i\n", nUsed, nTotal);
	fprintf(outputImage->fpLog, "; mosaicHopper3D nSolved     : %ld\n", (long)nSolved);
	fprintf(outputImage->fpLog, "; mosaicHopper3D nSolved3D   : %ld\n", (long)nSolved3D);
	fprintf(outputImage->fpLog, "; mosaicHopper3D nSolved2D   : %ld\n", (long)nSolved2D);
	fprintf(outputImage->fpLog, "; mosaicHopper3D frac3D      : %.4f\n",
			(nSolved > 0) ? (double)nSolved3D / (double)nSolved : 0.0);
	fprintf(outputImage->fpLog, "; mosaicHopper3D rejCond     : %ld\n", (long)nRejCond);
	fprintf(outputImage->fpLog, "; mosaicHopper3D rejClip     : %ld\n", (long)nRejClip);
	fprintf(outputImage->fpLog, "; mosaicHopper3D meanNobs    : %.3f\n",
			(nSolved > 0) ? nObsSum / (double)nSolved : 0.0);
	fprintf(outputImage->fpLog, "; mosaicHopper3D maxSigma3D  : %f\n", hopper3DMaxSigma);
	fprintf(outputImage->fpLog, "; mosaicHopper3D jointMaxSigma : %f\n", jointMaxSigma);
	fprintf(outputImage->fpLog, "; mosaicHopper3D note        : joint solve, no pairs -- pairOverCount/rho not applicable\n");
	/*
	  Free accumulators
	*/
	freeDoubleHop3(Nxx, outputImage->ySize);
	freeDoubleHop3(Nxy, outputImage->ySize);
	freeDoubleHop3(Nxz, outputImage->ySize);
	freeDoubleHop3(Nyy, outputImage->ySize);
	freeDoubleHop3(Nyz, outputImage->ySize);
	freeDoubleHop3(Nzz, outputImage->ySize);
	freeDoubleHop3(bxAcc, outputImage->ySize);
	freeDoubleHop3(byAcc, outputImage->ySize);
	freeDoubleHop3(bzAcc, outputImage->ySize);
	freeDoubleHop3(SddAcc, outputImage->ySize);
	freeImageHop3(nObsPh, outputImage->ySize);
	freeImageHop3(nObsRg, outputImage->ySize);
	freeImageHop3(nObsAz, outputImage->ySize);
	freeImageHop3(gW, outputImage->ySize);
	freeImageHop3(gW2, outputImage->ySize);
	/*  nObs/chi2Plane/vZ3D/errZ/modeBand are NOT freed here: ownership passes to outputImage so
	    mosaic3d can write the .nobs/.chi2/.vz3d/.ez/.mode bands, and mosaic3d frees them after.
	    Guarded rather than assigned outright (mosaicTrue3D.c:708-717 does the same) -- an
	    unconditional assignment leaks whichever round handed a plane over first. */
	if (outputImage->jointNObs == NULL) { outputImage->jointNObs = nObs; }
	else { freeImageHop3(nObs, outputImage->ySize); }
	if (outputImage->jointChi2 == NULL) { outputImage->jointChi2 = chi2Plane; }
	else { freeImageHop3(chi2Plane, outputImage->ySize); }
	if (outputImage->vZ3D == NULL) { outputImage->vZ3D = vZ3D; }
	else { freeImageHop3(vZ3D, outputImage->ySize); }
	if (outputImage->errorZ == NULL) { outputImage->errorZ = errZ; }
	else { freeImageHop3(errZ, outputImage->ySize); }
	if (outputImage->hopperMode == NULL) { outputImage->hopperMode = modeBand; }
	else { freeImageHop3(modeBand, outputImage->ySize); }
	if (sumW != NULL)
	{
		freeImageHop3(sumW, outputImage->ySize);
		freeImageHop3(sumWT, outputImage->ySize);
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
  DOUBLE accumulator planes.  mallocImage() only supplies float, so these are allocated
  row-by-row in the same shape (a row-pointer array over independently malloc'd rows) and freed
  the same way by freeDoubleHop3().  Zeroed on allocation -- the accumulators must start at 0.
*/
static HOPPER3DACC **mallocZeroDoubleHop3(int32_t nr, int32_t nc)
{
	HOPPER3DACC **tmp;
	int32_t i, j;
	tmp = (HOPPER3DACC **)malloc((size_t)nr * sizeof(HOPPER3DACC *));
	if (tmp == NULL)
	{
		error("mosaicHopper3D: malloc failed for double accumulator plane\n");
	}
	for (i = 0; i < nr; i++)
	{
		tmp[i] = (HOPPER3DACC *)malloc((size_t)nc * sizeof(HOPPER3DACC));
		if (tmp[i] == NULL)
		{
			error("mosaicHopper3D: malloc failed for double accumulator row\n");
		}
		for (j = 0; j < nc; j++)
		{
			tmp[i][j] = 0.0;
		}
	}
	return tmp;
}

static void freeDoubleHop3(HOPPER3DACC **image, int32_t nr)
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
  mallocImage() mallocs row by row and does NOT zero -- the accumulators must start at 0.
*/
static float **mallocZeroImageHop3(int32_t nr, int32_t nc)
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

static void freeImageHop3(float **image, int32_t nr)
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
static double computePhiFlatEarthM3dHop3(double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
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

	/* Flat-earth look angle (same formula as computePhiZM3dHop3 line for thetaDFlat) */
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
static double computePhiZM3dHop3(double *thetaD, double z, double azimuth, vhParams *vhParam, inputImageStructure *phaseImage, double Range, double Re,
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

static void setBufferHop3(inputImageStructure *inputImage, float *buf)
{
	int32_t i;
	for (i = 0; i < inputImage->azimuthSize; i++)
	{
		inputImage->image[i] = &(buf[i * inputImage->rangeSize]);
	}
}
