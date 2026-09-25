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
********** TRUE 3-COMPONENT VELOCITY FROM INTERFEROMETRIC PHASE ******************************

  "mosaic3d -true3DPhase".  Solves (vx, vy, vz) with NO surface-parallel constraint, from PHASE.
  DEFAULT OFF.

  Derived from make3DMosaicJoint.c by the same sin(psi) multiplication used for the offsets in
  mosaicTrue3D.c.  The 2D phase row is a = (cos g, sin g) - cot(psi)(dz/dx, dz/dy) with
  P2D = phase * 365.25/(twok nDays sin psi).  Multiplying a.vh = P2D through by sin(psi) and
  substituting vz = s.vh gives

      RANGE-LIKE row   u = ( sin(psi) cos g,  sin(psi) sin g,  -cos(psi) )
      P = phase * 365.25/(twok nDays),   sigma = phaseError * 365.25/(twok nDays)

  i.e. the 1/sin(psi) comes OUT of both the data and the sigma.  NO slope term in the row -- the
  slope enters only when the 3x3 is projected back to the surface-parallel subspace.

  WHY A PHASE VERSION EXISTS
  --------------------------
  Two reasons, and the second is the important one.

  1. PRECISION.  Phase sigma is ~1.24 m/yr against ~11.5 for range offsets.  Phase carries no
     azimuth component so the geometry is weaker (simulated cond 128 vs 40), but the noise is
     ~9x smaller, and simulation puts sigma_vz at ~1.2 m/yr -- BETTER than the offsets 3D solve
     despite the worse conditioning.

  2. A DIFFERENT IONOSPHERE PATH.  This is the discriminating test.  mosaicTrue3D returns a NISAR
     vz of -19 m/yr where the truth is ~-0.2, and that survives every explanation tried: it is
     not the SMB correction, not conditioning (measured per-pixel psi spread 1.33 deg, cond 40),
     not independent noise and not per-track noise (both give ZERO bias in simulation with the
     real geometry).  Only a ~15 m/yr COMMON-MODE range term reproduces it.  The offsets
     ionosphere correction IS applied (216 of 223 frames) and is pre-negated and ADDED; the phase
     screen is split-spectrum and SUBTRACTED, through entirely separate code.  So:

        phase 3D gives a sensible vz  -> the fault is in the offsets ionosphere correction
        phase 3D also gives ~-19      -> the fault is upstream of both

  NOTE ON THE VERTICAL AND PHASE.  Phase measures the same line of sight as a range offset, so it
  has the same -cos(psi) vertical sensitivity and is subject to the same common-mode ambiguity.
  A phase result close to truth would NOT mean phase is immune -- only that ITS ionosphere
  handling leaves a smaller residual.

  Everything else -- flat-earth vs terrain baselines, squint, tide, SMB, the buffer handling and
  the asc/desc list walk -- is inherited verbatim from make3DMosaicJoint.c.
*/
static void setBuffer(inputImageStructure *inputImage, float *buf);
static double computePhiZM3d(double *thetaD, double z, double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
							 double Range, double Re, double ReH, double ReHfixed, double thetaC, double thetaCfixedReH, double *phaseError);
static double computePhiFlatEarthM3d(double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
									 double Range, double Re, double ReHfixed, double thetaCfixedReH, double *phaseError);
static double **mallocZeroDoubleP3(int32_t nr, int32_t nc)
{
	double **tmp;
	int32_t i, j;
	tmp = (double **)malloc((size_t)nr * sizeof(double *));
	if (tmp == NULL) { error("mosaicTrue3DPhase: malloc failed\n"); }
	for (i = 0; i < nr; i++)
	{
		tmp[i] = (double *)malloc((size_t)nc * sizeof(double));
		if (tmp[i] == NULL) { error("mosaicTrue3DPhase: malloc failed row %i\n", i); }
		for (j = 0; j < nc; j++) { tmp[i][j] = 0.0; }
	}
	return tmp;
}

static void freeDoubleP3(double **image, int32_t nr)
{
	int32_t i;
	if (image == NULL) { return; }
	for (i = 0; i < nr; i++) { free(image[i]); }
	free(image);
}

static float **mallocZeroImageP3(int32_t nr, int32_t nc);
static double **mallocZeroDoubleP3(int32_t nr, int32_t nc);
static void freeDoubleP3(double **image, int32_t nr);
static void freeImage(float **image, int32_t nr);

void mosaicTrue3DPhase(inputImageStructure *ascImages, inputImageStructure *descImages,
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
	/*  3x3 normal equations, double (a 3x3 determinant cancels harder than a 2x2 and the
	    vertical column has far less directional diversity).  Matches mosaicTrue3D.c. */
	double **Nxx, **Nxy, **Nxz, **Nyy, **Nyz, **Nzz, **bxAcc, **byAcc, **bzAcc;
	float **vZ3D, **errZ;
	float **nObs;				 /* contributing measurements per pixel -- diagnostic */
	/*  MUST be double, not float: chi2 = Sdd - v.b is a difference of two large nearly
	    equal sums, so float here destroys it (symptom: chi2 clamps to 0).  N and b are
	    double for the 3x3, and Sdd has to match them. */
	double **SddAcc;
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
	Nxx = mallocZeroDoubleP3(outputImage->ySize, outputImage->xSize);
	Nxz = mallocZeroDoubleP3(outputImage->ySize, outputImage->xSize);
	Nyz = mallocZeroDoubleP3(outputImage->ySize, outputImage->xSize);
	Nzz = mallocZeroDoubleP3(outputImage->ySize, outputImage->xSize);
	bzAcc = mallocZeroDoubleP3(outputImage->ySize, outputImage->xSize);
	vZ3D = mallocZeroImageP3(outputImage->ySize, outputImage->xSize);
	errZ = mallocZeroImageP3(outputImage->ySize, outputImage->xSize);
	{
		int32_t ii, kk;
		for (ii = 0; ii < outputImage->ySize; ii++)
			for (kk = 0; kk < outputImage->xSize; kk++)
			{ vZ3D[ii][kk] = (float)-LARGEINT; errZ[ii][kk] = (float)-LARGEINT; }
	}
	Nxy = mallocZeroDoubleP3(outputImage->ySize, outputImage->xSize);
	Nyy = mallocZeroDoubleP3(outputImage->ySize, outputImage->xSize);
	bxAcc = mallocZeroDoubleP3(outputImage->ySize, outputImage->xSize);
	byAcc = mallocZeroDoubleP3(outputImage->ySize, outputImage->xSize);
	nObs = mallocZeroImageP3(outputImage->ySize, outputImage->xSize);
	SddAcc = mallocZeroDoubleP3(outputImage->ySize, outputImage->xSize);
	sumW = NULL;
	sumWT = NULL;
	chi2Plane = mallocZeroImageP3(outputImage->ySize, outputImage->xSize);
	if (outputImage->timeOverlapFlag == TRUE)
	{
		sumW = mallocZeroImageP3(outputImage->ySize, outputImage->xSize);
		sumWT = mallocZeroImageP3(outputImage->ySize, outputImage->xSize);
	}
	allImages = ascImages; /* consolodateLists() has already appended desc to asc */
	params = ascParams;
	nTotal = 0;
	for (phaseImage = allImages; phaseImage != NULL; phaseImage = phaseImage->next)
	{
		nTotal++;
	}
	fprintf(stderr, "mosaicTrue3DPhase: nTotal Images %i\n", nTotal);
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
		error("mosaicTrue3DPhase: malloc failed for per-thread image copies\n");
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
			double ax, ay, az, scale, P, sigma, w, dzdtSubmergence;
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
					/*  3-component row.  NO slope term: the slope enters only when the 3x3 is
					    projected back to the surface-parallel subspace.  B[][] is still computed
					    above (unchanged) because pass 2 uses dzdx/dzdy for the vz comparison. */
					ax = sin(psi) * cos(gamma);
					ay = sin(psi) * sin(gamma);
					az = -cos(psi);
					/*  Scale phase to velocity (m/yr) */
					/*  NO 1/sin(psi): the row carries sin(psi) explicitly. */
					scale = 365.25 / (twok * params->nDays);
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
					Nxx[iRow][jj] += w * ax * ax;
					Nxy[iRow][jj] += w * ax * ay;
					Nxz[iRow][jj] += w * ax * az;
					Nyy[iRow][jj] += w * ay * ay;
					Nyz[iRow][jj] += w * ay * az;
					Nzz[iRow][jj] += w * az * az;
					bxAcc[iRow][jj] += w * ax * P;
					byAcc[iRow][jj] += w * ay * P;
					bzAcc[iRow][jj] += w * az * P;
					SddAcc[iRow][jj] += w * P * P;
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
	fprintf(stderr, "mosaicTrue3DPhase: %i of %i images contributed\n", nUsed, nTotal);
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
			double n11, n12, n13, n22, n23, n33, bx, by, bz, det, lam;
			double Cinv[3][3];
			double vx, vy, vz, scX, scY, deltaOffCenter;
			double B[2][2], dzdx, dzdy;
			int32_t jj;
			y = (outputImage->originY + i * outputImage->deltaY) * MTOKM;
			for (jj = jMinAll; jj < jMaxAll; jj++)
			{
				if (nObs[i][jj] < 2.5)
				{
					/* THREE unknowns now, so three independent looks are the minimum. */
					nRejNoData++;
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					continue;
				}
				n11 = Nxx[i][jj]; n12 = Nxy[i][jj]; n13 = Nxz[i][jj];
				n22 = Nyy[i][jj]; n23 = Nyz[i][jj]; n33 = Nzz[i][jj];
				bx = bxAcc[i][jj]; by = byAcc[i][jj]; bz = bzAcc[i][jj];
				/*  Worst-direction formal sigma, 1/sqrt(lambdaMin(N3)), n-normalised.  Monotone
				    under adding data, exactly as in the 2D solvers.  lambdaMin3/invSym3 are
				    shared with mosaicTrue3D.c so the two cannot drift apart. */
				lam = lambdaMin3(n11, n12, n13, n22, n23, n33);
				if (!(lam > 0.0) ||
					(true3DMaxSigma > 0.0 &&
					 !((1.0 / sqrt(lam)) * sqrt((double)nObs[i][jj]) <= true3DMaxSigma)) ||
					invSym3(n11, n12, n13, n22, n23, n33, Cinv, &det) == FALSE ||
					!(Cinv[0][0] > 0.0) || !(Cinv[1][1] > 0.0) || !(Cinv[2][2] > 0.0))
				{
					nRejCond++;
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					continue;
				}
				vx = Cinv[0][0] * bx + Cinv[0][1] * by + Cinv[0][2] * bz;
				vy = Cinv[1][0] * bx + Cinv[1][1] * by + Cinv[1][2] * bz;
				vz = Cinv[2][0] * bx + Cinv[2][1] * by + Cinv[2][2] * bz;
				vZ3D[i][jj] = (float)vz;
				errZ[i][jj] = (float)Cinv[2][2];   /* VARIANCE; sqrt applied on output */
				scX = 1.0 / Cinv[0][0];
				scY = 1.0 / Cinv[1][1];
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
					double chi2 = SddAcc[i][jj] - (vx * bx + vy * by + vz * bz);
					double dof = (double)nObs[i][jj] - 3.0;
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
	fprintf(stderr, "mosaicTrue3DPhase: solved %ld  rejCond %ld  rejNoData %ld  meanNobs %.2f  maxNobs %.0f\n",
			(long)nSolved, (long)nRejCond, (long)nRejNoData,
			(nSolved > 0) ? nObsSum / (double)nSolved : 0.0, nObsMax);
	fprintf(outputImage->fpLog, "; mosaicTrue3DPhase nImagesUsed : %i of %i\n", nUsed, nTotal);
	fprintf(outputImage->fpLog, "; mosaicTrue3DPhase nSolved     : %ld\n", (long)nSolved);
	fprintf(outputImage->fpLog, "; mosaicTrue3DPhase rejCond     : %ld\n", (long)nRejCond);
	fprintf(outputImage->fpLog, "; mosaicTrue3DPhase meanNobs    : %.3f\n",
			(nSolved > 0) ? nObsSum / (double)nSolved : 0.0);
	fprintf(outputImage->fpLog, "; mosaicTrue3DPhase maxSigma    : %f\n", jointMaxSigma);
	fprintf(outputImage->fpLog, "; mosaicTrue3DPhase note        : joint solve, no pairs -- pairOverCount/rho not applicable\n");
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
	outputImage->vZ3D = vZ3D;
	outputImage->errorZ = errZ;
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
static float **mallocZeroImageP3(int32_t nr, int32_t nc)
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
