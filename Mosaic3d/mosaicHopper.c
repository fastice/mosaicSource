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
********** THE HOPPER: EVERY OBSERVABLE IN ONE PER-FRAME-WEIGHTED SOLVE *********************

  "mosaic3d -hopper".  DEFAULT OFF.  Puts PHASE, RANGE OFFSETS and AZIMUTH OFFSETS into a single
  per-pixel normal-equation system, each row carrying its own frame's sigma:

      PHASE          u = (cos g, sin g) - cot(psi) s,   P = phase * 365.25/(twok nDays sin psi)
      RANGE OFFSET   u = same direction,                P = dr    * 365.25/(nDays sin psi)
      AZIMUTH OFFSET u = (-sin g, cos g),               P = da    * 365.25/nDays

  (with -true3D each becomes the 3-component form of mosaicTrue3D.c: the first two gain a
  -cos(psi) vertical component and drop the slope term, the azimuth row gains a zero.)

  WHAT QUESTION THIS ANSWERS
  --------------------------
  Today a combined product is built by running SEPARATE ROUNDS -- phase, then crossing offsets,
  then speckle -- each producing its own weighted mosaic, which are then renormalised and averaged
  by redoNormalization()/endScale().  Averaging mosaics is NOT the same as weighting every frame
  individually: the round-level average uses each round's ACCUMULATED weight, so a frame that is
  excellent in a round with poor overall coverage is diluted, and a poor frame in a strong round
  rides along.  The hopper removes the intermediate step -- there is one weighting, at the frame.

  The bar it has to clear is high.  Measured on the 1600 m grid against the multi-year S1
  reference, std(d) in m/yr (lower is better):

      NISAR   opt_t10000_both (two rounds averaged)  3.930   <- best
              joint_phase                            4.728
              stj_ra (range+azimuth speckle)         6.815
              jointoff (range crossing)              8.503
      S1      s1_bo_joint (two rounds averaged)      2.242   <- best
              s1stj_ra                               2.677
              s1_rg_joint                            3.567
              s1_ph_joint                            4.803

  Note phase and range offsets SWAP RANKS between the two sensors, so no fixed blend is right for
  both -- which is the strongest argument for weighting at the frame.

  WHY THE THREE OBSERVABLES CAN SHARE AN ACCUMULATOR
  --------------------------------------------------
  They are independent in every way that matters.  Speckle-matching error and interferometric
  phase noise share no mechanism, and the IONOSPHERIC term is ANTI-CORRELATED between them --
  phase advance against group delay, equal and opposite in range units -- so combining them should
  partially CANCEL it rather than double-count.  The only genuinely shared term is the baseline
  solution, which is far smaller.  (An earlier version of this comment claimed the opposite; it
  was wrong.)

  Each image contributes up to three rows.  A frame with no phase ("nophase" in the input file)
  still contributes its offsets; a frame whose azparams failed still contributes range and phase;
  a frame outside timeThreshPhase still contributes its offsets.  Rows are suppressed
  individually, never the image -- the same rule as speckleTrackMosaicJoint.c.

  This was NOT true before 2026-08-30: the nophase and timeThreshPhase tests were image-level
  `continue`s, so on an archive where many frames carry offsets but no phase the mosaic was built
  from a fraction of the offsets.  Measured on Sentinel-1: 2099 of 4886 frames in one sector were
  "nophase" and contributed nothing at all.

  Derived from make3DMosaicJoint.c, which already consolidates the asc/desc lists and carries the
  flat-earth/terrain baselines, squint, tide, SMB and buffer handling; the offsets rows are lifted
  from mosaicTrue3D.c.  See Documents/hopper.md.
*/
static void setBufferHop(inputImageStructure *inputImage, float *buf);
static double computePhiZM3dHop(double *thetaD, double z, double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
							 double Range, double Re, double ReH, double ReHfixed, double thetaC, double thetaCfixedReH, double *phaseError);
static double computePhiFlatEarthM3dHop(double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
									 double Range, double Re, double ReHfixed, double thetaCfixedReH, double *phaseError);
static float **mallocZeroImageHop(int32_t nr, int32_t nc);
static void freeImageHop(float **image, int32_t nr);

void mosaicHopper(inputImageStructure *ascImages, inputImageStructure *descImages,
					   vhParams *ascParams, vhParams *descParams, xyDEM *dem,
					   outputImageStructure *outputImage, float fl, int32_t no3d, float timeThreshPhase,
					   referenceVelocity *refVel)
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
	int32_t hasPhase, framePhaseOK, useRangeRow, useAzRow;
	float azSLPixSize, rSLPixSize;
	double dumRe, dumReH, dumThetaC;
	float geoTolerance = 1e-3;						  /* as make3DMosaic: a few metres is plenty */
	float **vXimage, **vYimage, **vZimage, **errorX, **errorY;
	float **vxTmp, **vyTmp, **vzTmp, **fScale, **sxTmp, **syTmp;
	float **scaleX, **scaleY, **scaleZ;
	/* Normal-equation accumulators, one plane each, over the whole output grid.
	   float (not double) is safe: cancellation in Nxx*Nyy - Nxy^2 only bites as cond -> 0, and
	   those pixels are rejected by the conditioning gate anyway.  det is formed in double. */
	float **Nxx, **Nxy, **Nyy, **bxAcc, **byAcc;
	float **nObs, **nObsPh, **nObsRg, **nObsAz;				 /* contributing measurements per pixel -- diagnostic */
	/*  Weighted effective row count for the rejection gate: (sum w)^2 / sum(w^2).  See gateNEff
	    in mosaic3d.c -- must match mosaicHopper3D.c or the forced-2D equivalence breaks. */
	float **gW, **gW2;
	float **SddAcc; /* sum w d^2, for chi-square at solve time */
	float **chi2Plane;
	float **sumW, **sumWT;		 /* timeOverlapFlag only; NULL otherwise */
	float dum1, dum2;
	float azimuthMin, azimuthMax;
	int32_t iMin, iMax, jMin, jMax;					 /* current image's region */
	int32_t iMinAll, iMaxAll, jMinAll, jMaxAll;		 /* union over all contributing images */
	int32_t ii, nTotal, nUsed;
	int32_t i;
	int64_t nSolved, nRejCond, nRejNoData, nRejClip;
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
	Nxx = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	Nxy = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	Nyy = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	bxAcc = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	byAcc = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	nObs = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	nObsPh = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	nObsRg = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	nObsAz = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	gW = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	gW2 = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	SddAcc = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	sumW = NULL;
	sumWT = NULL;
	chi2Plane = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
	if (outputImage->timeOverlapFlag == TRUE)
	{
		sumW = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
		sumWT = mallocZeroImageHop(outputImage->ySize, outputImage->xSize);
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
		/*  -timePhaseThresh DOES NOT APPLY TO THE HOPPER.  It is a PAIR gate: it limits the
		    separation between the two images of a crossing pair, which is what stops the
		    pairwise solvers forming a combinatorial number of pairs that then get averaged
		    together for no gain.  The hopper forms no pairs -- it accumulates each image once --
		    so there is nothing for the threshold to limit, and the temporal extent of the
		    solution is set by the DATE RANGE (-date1/-date2) alone.
		    Until 2026-09-07 this code reused the parameter as a per-image window against the
		    mosaic centre date, which is a different quantity entirely: with a 1 Jun 2025 -
		    30 Nov 2026 mosaic the centre is ~1 Mar 2026, so passing the offsets value of 37
		    would have discarded every phase row outside late Jan - early Apr 2026.  It only
		    went unnoticed because the templates carry 10000, which never triggers it.
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
		setBufferHop(phaseImage, AImageBuffer);
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
			double B[2][2], dzdx, dzdy;
			double ax, ay, axPh, ayPh, gammaPh, hAngleSq, scale, P, sigma, w, dzdtSubmergence;
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
					    first such frame.  Reading it segfaults.  (This surfaced the moment the
					    image-level `continue` above became a phase-row suppression: before that,
					    a nophase frame never reached this line.) */
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
					  Identical to make3DMosaic, single-sided.  Skipped entirely for a frame with
					  no phase: these read the BASELINE solution, which a "nophase" frame need not
					  have, and the result would be unused anyway.
					*/
					phiZ = 0.0;
					phaseError = 0.0;
					if (pixHasPhase == TRUE)
					{
						if (params->applyFlatEarth)
						{
							phiZ = computePhiFlatEarthM3dHop(azimuth, params, myImg, Range, Re, ReHfixed, thetaCfixedReH, &phaseError);
						}
						else
						{
							/*  computePhiZM3dHop WRITES its thetaD argument (*thetaD = theta -
							    thetaC).  Give it a private copy: every offsets solver in the
							    codebase -- make3DOffsets, make3DOffsetsJoint,
							    speckleTrackMosaic -- feeds interpRangeOffsetInMeters and
							    computeSig2Base the thetaD that geometryInfo() produced, and the
							    offsets rows below must match them.  Letting the phase call
							    overwrite it made the hopper's offsets rows diverge from the
							    solver they are supposed to reduce to.  No effect on the
							    flat-earth (NISAR) path, which does not touch thetaD -- which is
							    why the NISAR reduction test passed at 0.000 regardless. */
							double thetaDPhase = thetaD;
							phiZ = computePhiZM3dHop(&thetaDPhase, zSp, azimuth, params, myImg, Range, Re, myReH, ReHfixed, thetaC, thetaCfixedReH, &phaseError);
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
					computeB(x, y, zWGS84, B, &dzdx, &dzdy, psi, psi, (xyDEM *)dem);
					if (sMask == SHELF) { B[0][0] = 0.0; B[0][1] = 0.0; }
					/*  ---- ROW 1: PHASE ------------------------------------------------
					    Slope coupling is shared by the two LOS observables, so B[][] computed
					    above serves both.  Each row is accumulated INDEPENDENTLY with its own
					    frame sigma -- that is the whole point of the hopper. */
					ax = cos(gamma) - B[0][0];
					ay = sin(gamma) - B[0][1];
					axPh = cos(gammaPh) - B[0][0];
					ayPh = sin(gammaPh) - B[0][1];
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
						scale = 365.25 / (twok * params->nDays * sin(psi));
						P = phase * scale;
						sigma = phaseError * scale;
						if (sigma > 0.0 && fabs(phase) < 0.9 * LARGEINT)
						{
							w = (sumW != NULL) ? phaseImage->weight / (sigma * sigma) : 1.0 / (sigma * sigma);
							if (w > 0.0)
							{
								Nxx[iRow][jj] += (float)(w * axPh * axPh);
								Nxy[iRow][jj] += (float)(w * axPh * ayPh);
								Nyy[iRow][jj] += (float)(w * ayPh * ayPh);
								bxAcc[iRow][jj] += (float)(w * axPh * P);
								byAcc[iRow][jj] += (float)(w * ayPh * P);
								SddAcc[iRow][jj] += (float)(w * P * P);
								nObs[iRow][jj] += 1.0f;
								nObsPh[iRow][jj] += 1.0f;
								gW[iRow][jj] += (float)w;
								gW2[iRow][jj] += (float)(w * w);
								/*  az = 0: this solver accumulates a 2x2 directly, so there is no
								    vertical component and no surface-parallel projection.  The
								    PIXEL record emits sx=sy=0 to match, making C the identity. */
								obsDumpRecord(iRow, jj, "phase", phaseImage->file,
											  phaseImage->julDay, (double)params->nDays,
											  axPh, ayPh, 0.0, P, sigma, w);
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
							scaleR = 365.25 / ((double)params->nDays * sin(psi));
							P = dr * scaleR;
							sigma = sigmaR * scaleR;
							if (sigma > 0.0)
							{
								w = (sumW != NULL) ? phaseImage->weight / (sigma * sigma) : 1.0 / (sigma * sigma);
								if (w > 0.0)
								{
									Nxx[iRow][jj] += (float)(w * ax * ax);
									Nxy[iRow][jj] += (float)(w * ax * ay);
									Nyy[iRow][jj] += (float)(w * ay * ay);
									bxAcc[iRow][jj] += (float)(w * ax * P);
									byAcc[iRow][jj] += (float)(w * ay * P);
									SddAcc[iRow][jj] += (float)(w * P * P);
									nObs[iRow][jj] += 1.0f;
									nObsRg[iRow][jj] += 1.0f;
									gW[iRow][jj] += (float)w;
									gW2[iRow][jj] += (float)(w * w);
									obsDumpRecord(iRow, jj, "range", phaseImage->file,
												  phaseImage->julDay, (double)params->nDays,
												  ax, ay, 0.0, P, sigma, w);
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
							aax = -sin(gamma);
							aay = cos(gamma);
							P = (double)da * kAz;
							sigma = sigmaA * kAz;
							if (sigma > 0.0)
							{
								w = (sumW != NULL) ? phaseImage->weight / (sigma * sigma) : 1.0 / (sigma * sigma);
								if (w > 0.0)
								{
									Nxx[iRow][jj] += (float)(w * aax * aax);
									Nxy[iRow][jj] += (float)(w * aax * aay);
									Nyy[iRow][jj] += (float)(w * aay * aay);
									bxAcc[iRow][jj] += (float)(w * aax * P);
									byAcc[iRow][jj] += (float)(w * aay * P);
									SddAcc[iRow][jj] += (float)(w * P * P);
									nObs[iRow][jj] += 1.0f;
									nObsAz[iRow][jj] += 1.0f;
									gW[iRow][jj] += (float)w;
									gW2[iRow][jj] += (float)(w * w);
									obsDumpRecord(iRow, jj, "azimuth", phaseImage->file,
												  phaseImage->julDay, (double)params->nDays,
												  aax, aay, 0.0, P, sigma, w);
								}
							}
						}
					}
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
	fprintf(stderr, "mosaicHopper: %i of %i images contributed\n", nUsed, nTotal);
	/*
	  ****************** PASS 2: solve once per pixel *****************************************
	*/
	nSolved = 0;
	nRejClip = 0;
	nRejCond = 0;
	nRejNoData = 0;
	nObsSum = 0.0;
	nObsMax = 0.0;
	if (nUsed > 0 && iMaxAll > iMinAll && jMaxAll > jMinAll)
	{
#pragma omp parallel for schedule(dynamic, 8) \
	reduction(+ : nSolved, nRejCond, nRejNoData, nRejClip, nObsSum) reduction(max : nObsMax)
		for (i = iMinAll; i < iMaxAll; i++)
		{
			double lat, lon, x, y, zWGS84;
			double nxx, nxy, nyy, bx, by, det, trHalf, disc, lambdaMin;
			double vx, vy, vz, scX, scY, deltaOffCenter, nGate, gateN, gateCap, gateSpd;
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
				/*  sx=sy=0: no surface-parallel projection in this solver -- see the row comment
				    above.  Emitted before the gates so a dump point still gets a PIXEL record
				    even where the pixel is ultimately rejected. */
				obsDumpPixel(i, jj, 0.0, 0.0);
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
				computeB(x, y, zWGS84, B, &dzdx, &dzdy, 1.0, 1.0, (xyDEM *)dem);
				vz = vx * dzdx + vy * dzdy;
				/*  Reduced chi-square: chi2 = Sdd - v.b (cross terms cancel at the solution).
				    ~1 => the measurements agree with each other to within their own sigmas. */
				if (chi2Plane != NULL)
				{
					double chi2 = (double)SddAcc[i][jj] - (vx * bx + vy * by);
					double dof = (double)nObs[i][jj] - 2.0;
					chi2Plane[i][jj] = (dof > 0.0 && chi2 > 0.0) ? (float)(chi2 / dof) : 0.0f;
					/*  -maxChi2: blunder screen.  Rows that disagree with EACH OTHER, which is
					    the one failure the formal sigma cannot see -- a single frame with an
					    unwrapping error gives a small sigma and a wildly wrong velocity.  Applied
					    after chi2Plane is stored so the diagnostic band still records why. */
					if (maxChi2 > 0.0 && (double)chi2Plane[i][jj] > maxChi2)
					{
						nRejCond++;
						vxTmp[i][jj] = (float)-LARGEINT;
						fScale[i][jj] = 0.0;
						continue;
					}
				}
				/*  Reference-velocity clip (-clipVel/-clipThresh).  Same rule the legacy speckle
				    round has always applied, now available under the hopper: reject a solved
				    pixel that disagrees with the reference map by more than clipThresh, and
				    only where one of the two is below 100 m/yr.  Inert unless refVel->clipFlag
				    is set, and a pixel with no reference value is kept.  Applied here, after
				    the solve and the chi2 screen, because it needs the final (vx, vy). */
				if (clipVelChecked(x, y, vx, vy, refVel) == TRUE)
				{
					nRejClip++;
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
	fprintf(stderr, "mosaicHopper: solved %ld  rejCond %ld  rejClip %ld  rejNoData %ld  meanNobs %.2f  maxNobs %.0f\n",
			(long)nSolved, (long)nRejCond, (long)nRejClip, (long)nRejNoData,
			(nSolved > 0) ? nObsSum / (double)nSolved : 0.0, nObsMax);
	fprintf(outputImage->fpLog, "; mosaicHopper nImagesUsed : %i of %i\n", nUsed, nTotal);
	fprintf(outputImage->fpLog, "; mosaicHopper nSolved     : %ld\n", (long)nSolved);
	fprintf(outputImage->fpLog, "; mosaicHopper rejCond     : %ld\n", (long)nRejCond);
	fprintf(outputImage->fpLog, "; mosaicHopper rejClip     : %ld\n", (long)nRejClip);
	fprintf(outputImage->fpLog, "; mosaicHopper meanNobs    : %.3f\n",
			(nSolved > 0) ? nObsSum / (double)nSolved : 0.0);
	fprintf(outputImage->fpLog, "; mosaicHopper jointMaxSigma : %f\n", jointMaxSigma);
	fprintf(outputImage->fpLog, "; mosaicHopper note        : joint solve, no pairs -- pairOverCount/rho not applicable\n");
	/*
	  Free accumulators
	*/
	freeImageHop(Nxx, outputImage->ySize);
	freeImageHop(Nxy, outputImage->ySize);
	freeImageHop(Nyy, outputImage->ySize);
	freeImageHop(bxAcc, outputImage->ySize);
	freeImageHop(byAcc, outputImage->ySize);
	freeImageHop(gW, outputImage->ySize);
	freeImageHop(gW2, outputImage->ySize);
	/* nObs is NOT freed: ownership passes to outputImage so mosaic3d can write the
	   ".nobs" diagnostic band.  mosaic3d frees it after writing. */
	outputImage->jointNObs = nObs;
	outputImage->jointChi2 = chi2Plane;
	if (sumW != NULL)
	{
		freeImageHop(sumW, outputImage->ySize);
		freeImageHop(sumWT, outputImage->ySize);
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
static float **mallocZeroImageHop(int32_t nr, int32_t nc)
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

static void freeImageHop(float **image, int32_t nr)
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
static double computePhiFlatEarthM3dHop(double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
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

	/* Flat-earth look angle (same formula as computePhiZM3dHop line for thetaDFlat) */
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
static double computePhiZM3dHop(double *thetaD, double z, double azimuth, vhParams *vhParam, inputImageStructure *phaseImage, double Range, double Re,
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

static void setBufferHop(inputImageStructure *inputImage, float *buf)
{
	int32_t i;
	for (i = 0; i < inputImage->azimuthSize; i++)
	{
		inputImage->image[i] = &(buf[i * inputImage->rangeSize]);
	}
}
