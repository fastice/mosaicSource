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
********** TRUE 3-COMPONENT VELOCITY FROM RANGE + AZIMUTH OFFSETS ****************************

  "mosaic3d -true3D".  Solves (vx, vy, vz) per pixel with NO surface-parallel constraint.
  DEFAULT OFF.  Offsets only -- no phase rows in this version.

  Read speckleTrackMosaicJoint.c first: the per-image loop, the gates, the ionosphere sign, the
  serial SVD pre-init and the chi-square band are identical and are not re-documented here.

  WHY
  ---
  The surface-parallel constraint vz = s.vh is cheap to carry and expensive to be wrong about.
  Measured (Documents/speckleTrackJoint.md §6), Greenland geometry, 48 range + 31 azimuth rows:

     cost of DROPPING the constraint   ~1.3 m/yr of horizontal noise  (8.1 vs 6.8)
     cost of KEEPING it with a DEM slope error of 0.02 / 0.05 / 0.10, on a 7000 m/yr glacier:
                                        121 / 309 / 618 m/yr of horizontal BIAS

  It also turns out to be MORE robust to a systematic range error, not less: a bias common to
  every frame is geometrically vertical (asc and desc both carry -cos psi), so the 3-component
  solve parks it in vz instead of leaking it into the horizontal.  Only an asc/desc-ANTIsymmetric
  bias hurts, and that hurts the constrained solve identically because it is geometrically
  indistinguishable from real horizontal motion.

  THE ROWS -- DERIVED FROM THE EXISTING CODE, NOT ASSUMED
  -------------------------------------------------------
  make3DOffsetsJoint builds  a = (cos g, sin g) - cot(psi) (dz/dx, dz/dy)  with
  P2D = delta * 365.25/(nDays sin psi).  Multiply a.vh = P2D through by sin(psi) and substitute
  vz = s.vh:

      sin(psi) (cos g vx + sin g vy) - cos(psi) vz  =  delta * 365.25/nDays

  so the three-component rows are

      RANGE    u_r = ( sin(psi) cos g,  sin(psi) sin g,  -cos(psi) )
      AZIMUTH  u_a = (      -sin g,          cos g,           0    )

  with P = delta * 365.25/nDays and sigma = sigmaR|A * 365.25/nDays -- i.e. the 1/sin(psi) comes
  OUT of both the data and the sigma, which is where it belongs.

  ** THE VERTICAL COMPONENT IS NEGATIVE. **  That sign is forced by the algebra above, not
  chosen, and it is the single easiest thing in this file to get wrong.  It is verified by the
  projection identity below, which is the reduction test: with

      C = [[1,0],[0,1],[sx,sy]],   N2 = C^T N3 C,   b2 = C^T b3

  the projected solve must reproduce -speckleTrackJoint EXACTLY.  If the sign of -cos(psi) were
  flipped, the projection would still be self-consistent but would disagree with the 2D solver.
  -true3DProject forces that projection and exists solely for this test.

  AZIMUTH ROWS MATTER MORE HERE THAN ANYWHERE ELSE
  ------------------------------------------------
  Without them the 3x3 is poorly conditioned (cond 124 vs 15 for NISAR geometry) and at a single
  incidence angle it is exactly singular -- there is no third direction.  It is the ACROSS-swath
  spread of psi that makes range-only 3D possible at all, and azimuth that makes it good.
  Consequently azimuthAccuracyVar() is not optional here: an over-trusted azimuth row in a 3x3
  corrupts the horizontal directly, rather than merely inflating a weight.  See
  Documents/speckleTrackJoint.md §4, where the azimuth budget was measured 9x too small.

  vz HAS ITS OWN BUFFER
  ---------------------
  vZ3D/errorZ (geocode.h), NOT the vzTmp/image3 plane.  Nine modules write vzTmp and -timeOverlap
  repurposes it as the ".dT" band (mosaic3d.c:984), so there is no vz plane to "reserve" in the
  mode every calibration run uses.  vx/vy go through redoNormalization()/endScale() weighted
  exactly like any other round; only vz/ez are written directly, since this is their only writer.

  ACCUMULATORS ARE double.  A 3x3 determinant cancels harder than a 2x2 and the vertical column
  has far less directional diversity than the horizontal ones.  To measure what that choice
  actually costs, flip TRUE3DACC to float below and rebuild -- deliberately a compile-time switch
  rather than a CLI flag, so there is exactly one accumulation code path and no chance of the two
  drifting apart.  The result of that comparison belongs in Documents/true3DPlan.md.
*/

#define TRUE3DACC double

static TRUE3DACC **mallocZeroDouble3D(int32_t nr, int32_t nc);
static void freeDouble3D(TRUE3DACC **image, int32_t nr);
static float **mallocZeroFloat3D(int32_t nr, int32_t nc);
static void freeFloat3D(float **image, int32_t nr);
double lambdaMin3(double a11, double a12, double a13, double a22, double a23, double a33);
int32_t invSym3(double a11, double a12, double a13, double a22, double a23, double a33,
					   double C[3][3], double *det);

void mosaicTrue3D(inputImageStructure *images, vhParams *params, outputImageStructure *outputImage,
				  float fl, referenceVelocity *refVel)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	extern int32_t indentRegionOutput;
	inputImageStructure *currentImage;
	vhParams *currentParams;
	conversionDataStructure *cP;
	ShelfMask *shelfMask;
	xyDEM *vCorrect;
	xyDEM *vzDem;
	double Re, ReH, thetaC, ddum1, ddum2;
	double tCenter, tOffCenter;
	float **vXimage, **vYimage, **vZimage, **errorX, **errorY;
	float **vxTmp, **vyTmp, **vzTmp, **fScale, **sxTmp, **syTmp;
	float **scaleX, **scaleY, **scaleZ;
	/* Normal equations, 6 + 3 planes.  double by default -- see the header note. */
	TRUE3DACC **Nxx, **Nxy, **Nxz, **Nyy, **Nyz, **Nzz, **bxAcc, **byAcc, **bzAcc, **SddAcc;
	/*  Time-weighted centre offset for the ".dT" band under -timeOverlap.  REQUIRED:
	    filterDT() (mosaic3d.c) nulls vx/vy/ex/ey wherever |dT| exceeds the half-window,
	    so a solver that leaves dT unset silently erases its own output. */
	TRUE3DACC **sumW, **sumWT;
	/*  DIAGNOSTIC (-true3DDiag): weighted mean RANGE RATE, split by pass direction.
	    This is the discriminating test for the unexplained NISAR vz.  The vertical row
	    component is -cos(psi), the SAME sign for ascending and descending, so:
	      * a common-mode range signal (same sign both passes) is indistinguishable from vz
	        and will drive it;
	      * genuine horizontal motion appears with OPPOSITE sign in the two passes.
	    Comparing mean(P_asc) with mean(P_desc) therefore separates a real data-level
	    common-mode term from a geometry or sign error, which no test so far distinguishes. */
	TRUE3DACC **sumWPasc, **sumWasc, **sumWPdesc, **sumWdesc;
	/*  DIAGNOSTIC 2: per-pixel weighted mean and spread of psi over the RANGE rows.
	    The vertical is separable ONLY through the variation of cos(psi) among the rows at a
	    pixel -- with identical psi the 3x3 conditioning collapses (cond 28 -> 1182 for real
	    Greenland azimuths, verified from first principles).  Nothing measured so far looks at
	    this: per-frame MLIncidenceCenter describes pointing repeatability, not the per-pixel
	    spread, which also depends on where the pixel falls in each frame's swath. */
	TRUE3DACC **sumWpsi, **sumWpsi2;
	float **nObs, **nObsAz, **chi2Plane;
	float **vZ3D, **errZ;
	/*  Per-pixel xy surface slope, captured in pass 1 from xyGetZandSlope()'s (dz/dr, dz/da)
	    rotated back to xy.  It is image-INDEPENDENT (the rotation is orthogonal and exact), so
	    last-writer-wins is well defined.  Needed only by the -true3DProject reduction test, and
	    it must be THIS slope -- not computeB()'s -- or the test compares against the wrong 2D
	    solver: xyGetZandSlope clamps at 0.12 with a psScale correction and a low-elevation
	    zeroing, computeB clamps at 0.25 with neither. */
	float **sxPlane, **syPlane;
	float azSLPixSize, rSLPixSize;
	float azimuthMin, azimuthMax;
	int32_t iMin, iMax, jMin, jMax;
	int32_t iMinAll, iMaxAll, jMinAll, jMaxAll;
	int32_t *jRowMin, *jRowMax;
	int32_t count, total, nUsed, useAzRow;
	int32_t i;
	int64_t nSolved, nRejCond, nRejNoData, nClip;
	double nObsSum, nObsAzSum, nObsMax;
	struct timeval funcStart, funcEnd;
	gettimeofday(&funcStart, NULL);

	fprintf(stderr, "**** TRUE 3D SOLUTION (range + azimuth, no surface-parallel) ****\n");
	fprintf(outputImage->fpLog, ";\n; Entering mosaicTrue3D(.c)\n;\n");

	if (outputImage->outputRAFlag == TRUE)
	{
		error("mosaicTrue3D: -outputRA is not supported with -true3D\n");
	}
	shelfMask = outputImage->shelfMask;
	vCorrect = outputImage->verticalCorrection;

	setupBuffers(outputImage, &vXimage, &vYimage, &vZimage, &scaleX, &scaleY, &scaleZ,
				 &vxTmp, &vyTmp, &vzTmp, &sxTmp, &syTmp, &fScale, &errorX, &errorY);
	computeScale((float **)vXimage, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0, (double)(-LARGEINT));
	undoNormalization(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, fScale, FALSE);

	Nxx = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	Nxy = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	Nxz = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	Nyy = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	Nyz = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	Nzz = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	bxAcc = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	byAcc = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	bzAcc = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	SddAcc = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	nObs = mallocZeroFloat3D(outputImage->ySize, outputImage->xSize);
	nObsAz = mallocZeroFloat3D(outputImage->ySize, outputImage->xSize);
	chi2Plane = mallocZeroFloat3D(outputImage->ySize, outputImage->xSize);
	vZ3D = mallocZeroFloat3D(outputImage->ySize, outputImage->xSize);
	errZ = mallocZeroFloat3D(outputImage->ySize, outputImage->xSize);
	sumW = NULL;
	sumWT = NULL;
	sumWPasc = sumWasc = sumWPdesc = sumWdesc = NULL;
	sumWpsi = sumWpsi2 = NULL;
	if (true3DDiag == TRUE)
	{
		sumWPasc = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
		sumWasc = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
		sumWPdesc = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
		sumWdesc = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
		sumWpsi = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
		sumWpsi2 = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	}
	if (outputImage->timeOverlapFlag == TRUE)
	{
		sumW = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
		sumWT = mallocZeroDouble3D(outputImage->ySize, outputImage->xSize);
	}
	sxPlane = mallocZeroFloat3D(outputImage->ySize, outputImage->xSize);
	syPlane = mallocZeroFloat3D(outputImage->ySize, outputImage->xSize);
	for (i = 0; i < outputImage->ySize; i++)
	{
		int32_t k;
		for (k = 0; k < outputImage->xSize; k++)
		{
			vZ3D[i][k] = (float)-LARGEINT;
			errZ[i][k] = (float)-LARGEINT;
		}
	}

	total = 0;
	for (currentImage = images; currentImage != NULL; currentImage = currentImage->next)
	{
		total++;
	}
	tCenter = (outputImage->jd1 + outputImage->jd2 + 1.) * 0.5;
	iMinAll = outputImage->ySize;
	jMinAll = outputImage->xSize;
	iMaxAll = 0;
	jMaxAll = 0;
	nUsed = 0;
	count = 0;
	vzDem = NULL;
	int nthreads = omp_get_max_threads();
	inputImageStructure *localImgs = (inputImageStructure *)malloc((size_t)nthreads * sizeof(inputImageStructure));
	jRowMin = (int32_t *)malloc((size_t)outputImage->ySize * sizeof(int32_t));
	jRowMax = (int32_t *)malloc((size_t)outputImage->ySize * sizeof(int32_t));
	if (localImgs == NULL || jRowMin == NULL || jRowMax == NULL)
	{
		error("mosaicTrue3D: malloc failed\n");
	}
	/*
	  ****************** PASS 1: accumulate 3x3 normal equations ******************************
	*/
	currentParams = params;
	for (currentImage = images; currentImage != NULL; currentImage = currentImage->next, currentParams = currentParams->next)
	{
		count++;
		tOffCenter = currentImage->julDay + currentParams->nDays * 0.5;
		/*  Partial-temporal-overlap weights are carried in the input file (setupquarters writes
		    a fraction < 1 for frames only partly inside the mosaic window).  They are honoured
		    ONLY when -timeOverlap is set; without it they would be silently discarded and the
		    mosaic would weight a 4%-overlap frame the same as a full one.  speckleTrackMosaic
		    and speckleTrackMosaicJoint have always treated that as fatal -- match them.
		    Added 2026-08-30 after an S1 run (weights 0.042..1.0) silently ignored them here
		    while the 2D path correctly refused. */
		if (fabs(currentImage->weight - 1.0) > 0.01 && outputImage->timeOverlapFlag == FALSE)
		{
			error("non unity weight, but overlap flag not set\n");
		}
		if (currentParams->offsets.rFile == NULL || currentImage->weight < 0.00001)
		{
			continue;
		}
		fprintf(stderr, "\033[1;31m3D %i of %i\033[0m %s\n", count, total, currentParams->offsets.file);
		indentRegionOutput = FALSE;
		cP = setupGeoConversions(currentImage, &azSLPixSize, &rSLPixSize, &Re, &ReH, &thetaC, &ddum1, &ddum2);
		getRegion(currentImage, &iMin, &iMax, &jMin, &jMax, outputImage);
		if (iMin > iMax || jMin > jMax)
		{
			continue;
		}
		getAzimuthBoundsForXYBox(iMin, iMax, jMin, jMax, currentImage, outputImage, &azimuthMin, &azimuthMax);
		if (azimuthMin == 0.0f && azimuthMax == 0.0f)
		{
			continue;
		}
		getRowBounds(currentImage, outputImage, iMin, iMax, jMin, jMax, jRowMin, jRowMax);
		/*  GUARDED READ -- see mosaicHopper.c: a corrupt input for this product is recorded
		    and skipped instead of exiting.  Disarmed before the parallel region. */
		jmp_buf readJmp;
		if (setjmp(readJmp) != 0)
		{
			errorRecoveryJmp = NULL;
			recordFailedProduct(currentParams->offsets.rFile, errorRecoveryMsg);
			continue;
		}
		errorRecoveryJmp = &readJmp;
		readOffsetDataAndParams(&(currentParams->offsets), azimuthMin, azimuthMax);
		/* MANDATORY serial pre-init -- both routines lazily malloc global SVD workspace. */
		if (currentParams->offsets.deltaB != DELTABNONE)
		{
			double bnS, bpS;
			svAzOffset(currentImage, &(currentParams->offsets), 0.0, 0.0);
			svInterpBnBp(currentImage, &(currentParams->offsets), 0.0, &bnS, &bpS);
		}
		errorRecoveryJmp = NULL;
		if (currentParams->offsets.sigmaRresidual < 0)
		{
			fprintf(stderr, "mosaicTrue3D: %s no range baseline solution (sigma<0) -- skipping\n",
					currentParams->offsets.rParamsFile);
			currentParams->offsets.rFile = NULL;
			continue;
		}
		useAzRow = TRUE;
		if (currentParams->offsets.file == NULL || currentParams->offsets.da == NULL)
		{
			useAzRow = FALSE;
		}
		else if (currentParams->offsets.sigmaAresidual < 0 ||
				 currentParams->offsets.sigmaAresidual > outputImage->sigmaAThresh)
		{
			fprintf(stderr, "\tazimuth row suppressed (sigmaAresidual %f)\n", currentParams->offsets.sigmaAresidual);
			useAzRow = FALSE;
		}
		if (noAzimuthRows == TRUE)
		{
			useAzRow = FALSE;
		}
		nUsed++;
		if (vzDem == NULL)
		{
			vzDem = &(currentParams->xydem);
		}
		iMinAll = min(iMinAll, iMin);
		jMinAll = min(jMinAll, jMin);
		iMaxAll = max(iMaxAll, iMax);
		jMaxAll = max(jMaxAll, jMax);
		{
			int t;
			for (t = 0; t < nthreads; t++)
			{
				localImgs[t] = *currentImage;
			}
		}
#pragma omp parallel
		{
			inputImageStructure *myImg = &localImgs[omp_get_thread_num()];
			double lat, lon, x, y, zSp, zWGS84;
			double range, azimuth, Range, myReH, theta, thetaD, psi;
			double dzda, dzdr, demError, sig2Base, sig2Off, sigmaR, sigmaA;
			double hAngle, xyAngle, gamma, savedLastTime;
			double urx, ury, urz, aax, aay, kScale, sdzdx, sdzdy;
			double Pr, Pa, sigR, sigA, wR, wA, dzdtSubmergence;
			float da, dr, ionCorrection;
			unsigned char sMask;
			int32_t jj, iRow, drValid, daValid;
#pragma omp for schedule(dynamic, 8)
			for (iRow = iMin; iRow < iMax; iRow++)
			{
				if ((iRow % 100) == 0)
				{
					fprintf(stderr, "\t--+ %i\n", iRow);
				}
				y = (outputImage->originY + iRow * outputImage->deltaY) * MTOKM;
				for (jj = jRowMin[iRow]; jj < jRowMax[iRow]; jj++)
				{
					x = (outputImage->originX + jj * outputImage->deltaX) * MTOKM;
					xytoll1(x, y, HemiSphere, &lat, &lon, Rotation, outputImage->slat);
					xyGetZandSlope(lat, lon, x, y, &zSp, &zWGS84, &dzda, &dzdr, cP, currentParams, myImg);
					if (!(zWGS84 > (MINELEVATION + 1) && zWGS84 < 10000.0))
					{
						continue;
					}
					llToImageNew(lat, lon, zWGS84, &range, &azimuth, myImg);
					geometryInfo(cP, myImg, azimuth, range, zSp, thetaC, &myReH, &Range, &theta, &thetaD, &psi, zSp);
					da = interpAzOffset(range, azimuth, &(currentParams->offsets), myImg, Range, theta, azSLPixSize);
/* azimuth ionosphere: no-op unless the az fit recorded one and -useAzIonosphere is set */
if (da > -0.98 * LARGEINT)
	da += azIonCorrectionMeters(&(currentParams->offsets), myImg, range, azimuth, azSLPixSize);
					dr = interpRangeOffsetInMeters(range, azimuth, &(currentParams->offsets), myImg,
												   Range, thetaD, rSLPixSize, theta, &demError);
					if (currentParams->offsets.rOffCorrection.rangeOffsetCorrection != NULL)
					{
						double rangeSLC, azimuthSLC;
						computeSLCFromMLCoords(myImg, range, azimuth, &rangeSLC, &azimuthSLC);
						ionCorrection = interpolateOffsetIonCorrectionInPixels(&(currentParams->offsets.rOffCorrection),
																			   rangeSLC, azimuthSLC, -LARGEINT, 0.0);
					}
					else
					{
						ionCorrection = 0.0;
					}
					if (ionCorrection > -0.98 * LARGEINT)
					{
						dr += ionCorrection * rSLPixSize;
					}
					sMask = GROUNDED;
					if (shelfMask != NULL)
					{
						sMask = getShelfMask(shelfMask, x, y);
					}
					if (sMask == NOSOLUTION || sMask == GROUNDINGZONE)
					{
						continue;
					}
					if (sMask == SHELF && outputImage->noTide == TRUE)
					{
						continue;
					}
					drValid = (fabs(dr) < 13.0E4);
					daValid = (useAzRow == TRUE) && (fabs(da) < 10.0e4);
					if (drValid == FALSE && daValid == FALSE)
					{
						continue;
					}
					savedLastTime = myImg->lastTime;
					hAngle = computeHeading(lat, lon, 0.0, myImg, cP);
					myImg->lastTime = savedLastTime;
					computeXYangle(lat, lon, &xyAngle, currentParams->xydem);
					gamma = xyAngle - hAngle;
					/*  Capture the xy slope for the -true3DProject reduction test.  Exactly the
					    quantity speckleTrackMosaicJoint builds its 2D row from, including the
					    SHELF zeroing, so the projection compares like with like. */
					rotateFlowDirectionToXY(dzdr, dzda, &sdzdx, &sdzdy, xyAngle, hAngle);
					if (sMask == SHELF)
					{
						sdzdx = 0.0;
						sdzdy = 0.0;
					}
					sxPlane[iRow][jj] = (float)sdzdx;
					syPlane[iRow][jj] = (float)sdzdy;
					/* NOTE: no 1/sin(psi).  The 3D row carries sin(psi) explicitly. */
					kScale = 365.25 / (double)currentParams->nDays;
					if (drValid == TRUE)
					{
						sigmaR = interpRangeSigma(range, azimuth, &(currentParams->offsets), myImg, Range, thetaD, rSLPixSize);
						sig2Base = computeSig2Base(sin(thetaD), cos(thetaD), azimuth, myImg, &(currentParams->offsets));
						sigmaR = sqrt(sigmaR * sigmaR + demError * demError + sig2Base +
									  rangeAccuracyVar(&(currentParams->offsets)));
						dr -= shelfMaskCorrection(myImg, currentParams, sMask, x, y, psi, &sigmaR);
						if (vCorrect != NULL)
						{
							/*  NOTE: with vz solved rather than assumed, an SMB vertical
							    correction is arguably double-counting.  Applied anyway so the
							    data vector matches every other solver; the correct treatment is
							    to leave it OFF for -true3D and let the solve find the vertical.
							    Flagged in Documents/true3DPlan.md, not decided here. */
							dzdtSubmergence = interpVCorrect(x, y, vCorrect);
							dr -= -dzdtSubmergence * cos(psi) * (double)currentParams->nDays / 365.25;
						}
						urx = sin(psi) * cos(gamma);
						ury = sin(psi) * sin(gamma);
						urz = -cos(psi);
						Pr = (double)dr * kScale;
						sigR = sigmaR * kScale;
						if (sigR > 0.0)
						{
							wR = (outputImage->timeOverlapFlag == TRUE) ? currentImage->weight / (sigR * sigR)
																		: 1.0 / (sigR * sigR);
							if (wR > 0.0)
							{
								Nxx[iRow][jj] += wR * urx * urx;
								Nxy[iRow][jj] += wR * urx * ury;
								Nxz[iRow][jj] += wR * urx * urz;
								Nyy[iRow][jj] += wR * ury * ury;
								Nyz[iRow][jj] += wR * ury * urz;
								Nzz[iRow][jj] += wR * urz * urz;
								bxAcc[iRow][jj] += wR * urx * Pr;
								byAcc[iRow][jj] += wR * ury * Pr;
								bzAcc[iRow][jj] += wR * urz * Pr;
								SddAcc[iRow][jj] += wR * Pr * Pr;
								nObs[iRow][jj] += 1.0f;
								if (sumWpsi != NULL)
								{
									sumWpsi[iRow][jj] += wR * psi;
									sumWpsi2[iRow][jj] += wR * psi * psi;
								}
								if (sumWPasc != NULL)
								{
									if (myImg->passType == ASCENDING)
									{
										sumWPasc[iRow][jj] += wR * Pr;
										sumWasc[iRow][jj] += wR;
									}
									else
									{
										sumWPdesc[iRow][jj] += wR * Pr;
										sumWdesc[iRow][jj] += wR;
									}
								}
								if (sumW != NULL)
								{
									sumW[iRow][jj] += wR;
									sumWT[iRow][jj] += wR * (tOffCenter - tCenter);
								}
							}
						}
					}
					if (daValid == TRUE)
					{
						sigmaA = interpAzSigma(range, azimuth, &(currentParams->offsets), myImg, Range, theta, azSLPixSize);
						sig2Off = computeSig2AzParam(sin(theta), cos(theta), azimuth, Range, myImg, &(currentParams->offsets));
						sigmaA = sqrt(sigmaA * sigmaA + sig2Off + azimuthAccuracyVar(&(currentParams->offsets)));
						sigA = sigmaA * kScale;
						/* Azimuth row: horizontal, no vertical sensitivity, no slope term. */
						aax = -sin(gamma);
						aay = cos(gamma);
						Pa = (double)da * kScale;
						if (sigA > 0.0)
						{
							wA = (outputImage->timeOverlapFlag == TRUE) ? currentImage->weight / (sigA * sigA)
																		: 1.0 / (sigA * sigA);
							if (wA > 0.0)
							{
								Nxx[iRow][jj] += wA * aax * aax;
								Nxy[iRow][jj] += wA * aax * aay;
								Nyy[iRow][jj] += wA * aay * aay;
								bxAcc[iRow][jj] += wA * aax * Pa;
								byAcc[iRow][jj] += wA * aay * Pa;
								SddAcc[iRow][jj] += wA * Pa * Pa;
								nObs[iRow][jj] += 1.0f;
								nObsAz[iRow][jj] += 1.0f;
								if (sumW != NULL)
								{
									sumW[iRow][jj] += wA;
									sumWT[iRow][jj] += wA * (tOffCenter - tCenter);
								}
							}
						}
					}
				} /* j */
			}	  /* i */
		}		  /* omp parallel */
		currentImage->used = TRUE;
	}			  /* image loop */
	fprintf(stderr, "mosaicTrue3D: %i of %i images contributed\n", nUsed, total);
	/*
	  ****************** PASS 2: solve 3x3 (or its 2D projection) per pixel *******************
	*/
	nSolved = 0;
	nRejCond = 0;
	nRejNoData = 0;
	nClip = 0;
	nObsSum = 0.0;
	nObsAzSum = 0.0;
	nObsMax = 0.0;
	if (nUsed > 0 && iMaxAll > iMinAll && jMaxAll > jMinAll)
	{
#pragma omp parallel for schedule(dynamic, 8) \
	reduction(+ : nSolved, nRejCond, nRejNoData, nClip, nObsSum, nObsAzSum) reduction(max : nObsMax)
		for (i = iMinAll; i < iMaxAll; i++)
		{
			double lat, lon, x, y, zWGS84;
			double n11, n12, n13, n22, n23, n33, bx, by, bz;
			double Cinv[3][3], det, lam, vx, vy, vz, scX, scY, nRange;
			double B[2][2], dzdx, dzdy;
			int32_t jj;
			y = (outputImage->originY + i * outputImage->deltaY) * MTOKM;
			for (jj = jMinAll; jj < jMaxAll; jj++)
			{
				/*  THREE unknowns need three rows -- but not any three.  A 3-component solve
				    REQUIRES at least one AZIMUTH row: azimuth is the only observable that
				    constrains the along-track horizontal directly, and without it the third
				    dimension rests solely on the spread of psi across the swath, which is
				    marginal (measured per-pixel psi std 1.33 deg; cond 40 with azimuth, 128
				    without).  A track that happens to lack azimuth still contributes its range
				    row -- the requirement is one azimuth row AT THE PIXEL, not per image. */
				nRange = (double)nObs[i][jj] - (double)nObsAz[i][jj];
				/*  -noAzimuthRows is an explicit request for a range-only solve, so it
				    overrides the azimuth requirement.  Only meaningful when the geometry
				    supplies the third dimension some other way -- e.g. NISAR (left-looking)
				    combined with S1 (right-looking), whose antiparallel look pairs make
				    range-only 3D well conditioned (cond 28 vs 1182 non-antiparallel). */
				if (nObs[i][jj] < 2.5 || nRange < 1.5 ||
					(noAzimuthRows == FALSE && nObsAz[i][jj] < 0.5))
				{
					nRejNoData++;
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					continue;
				}
				n11 = Nxx[i][jj]; n12 = Nxy[i][jj]; n13 = Nxz[i][jj];
				n22 = Nyy[i][jj]; n23 = Nyz[i][jj]; n33 = Nzz[i][jj];
				bx = bxAcc[i][jj]; by = byAcc[i][jj]; bz = bzAcc[i][jj];
				x = (outputImage->originX + jj * outputImage->deltaX) * MTOKM;
				if (true3DProject == TRUE)
				{
					/*  REDUCTION TEST PATH.  Project onto the surface-parallel subspace with
					    C = [[1,0],[0,1],[sx,sy]] and solve the resulting 2x2, which must
					    reproduce -speckleTrackJoint exactly.  Not a production path. */
					double sx, sy, m11, m12, m22, p1, p2, d2;
					sx = (double)sxPlane[i][jj];
					sy = (double)syPlane[i][jj];
					m11 = n11 + 2.0 * sx * n13 + sx * sx * n33;
					m12 = n12 + sx * n23 + sy * n13 + sx * sy * n33;
					m22 = n22 + 2.0 * sy * n23 + sy * sy * n33;
					p1 = bx + sx * bz;
					p2 = by + sy * bz;
					d2 = m11 * m22 - m12 * m12;
					if (!(d2 > 0.0))
					{
						nRejCond++;
						vxTmp[i][jj] = (float)-LARGEINT;
						fScale[i][jj] = 0.0;
						continue;
					}
					vx = (m22 * p1 - m12 * p2) / d2;
					vy = (m11 * p2 - m12 * p1) / d2;
					vz = sx * vx + sy * vy;
					scX = d2 / m22;
					scY = d2 / m11;
					errZ[i][jj] = (float)-LARGEINT;
					vZ3D[i][jj] = (float)vz;
				}
				else
				{
					lam = lambdaMin3(n11, n12, n13, n22, n23, n33);
					/*  n-normalisation uses the RANGE count, not the total.  Azimuth rows have
					    ZERO vertical sensitivity (u_a = (-sin g, cos g, 0)), so they cannot
					    inform the worst-constrained direction, which for this geometry is
					    vertical or near it.  Normalising by the total would credit the vertical
					    with measurements that say nothing about it. */
					if (!(lam > 0.0) ||
						(true3DMaxSigma > 0.0 &&
						 !((1.0 / sqrt(lam)) * sqrt(nRange) <= true3DMaxSigma)) ||
						invSym3(n11, n12, n13, n22, n23, n33, Cinv, &det) == FALSE)
					{
						nRejCond++;
						vxTmp[i][jj] = (float)-LARGEINT;
						fScale[i][jj] = 0.0;
						continue;
					}
					vx = Cinv[0][0] * bx + Cinv[0][1] * by + Cinv[0][2] * bz;
					vy = Cinv[1][0] * bx + Cinv[1][1] * by + Cinv[1][2] * bz;
					vz = Cinv[2][0] * bx + Cinv[2][1] * by + Cinv[2][2] * bz;
					if (!(Cinv[0][0] > 0.0) || !(Cinv[1][1] > 0.0) || !(Cinv[2][2] > 0.0))
					{
						nRejCond++;
						vxTmp[i][jj] = (float)-LARGEINT;
						fScale[i][jj] = 0.0;
						continue;
					}
					scX = 1.0 / Cinv[0][0];
					scY = 1.0 / Cinv[1][1];
					vZ3D[i][jj] = (float)vz;
					errZ[i][jj] = (float)Cinv[2][2]; /* VARIANCE; sqrt applied on output */
				}
				if (jointErrScale != 1.0)
				{
					double fInf = jointErrScale * jointErrScale;
					scX /= fInf;
					scY /= fInf;
					if (errZ[i][jj] > (-LARGEINT + 1))
					{
						errZ[i][jj] *= (float)fInf;
					}
				}
				if (!(scX > -1000. && scX < 1000.) || !(scY > -1000. && scY < 1000.))
				{
					nRejCond++;
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					vZ3D[i][jj] = (float)-LARGEINT;
					errZ[i][jj] = (float)-LARGEINT;
					continue;
				}
				if (refVel->clipFlag == TRUE && clipVel((float)x, (float)y, (float)vx, (float)vy, refVel) == TRUE)
				{
					nClip++;
					vxTmp[i][jj] = (float)-LARGEINT;
					fScale[i][jj] = 0.0;
					vZ3D[i][jj] = (float)-LARGEINT;
					errZ[i][jj] = (float)-LARGEINT;
					continue;
				}
				{
					double chi2 = SddAcc[i][jj] - (vx * bx + vy * by + vz * bz);
					double dof = (double)nObs[i][jj] - (true3DProject == TRUE ? 2.0 : 3.0);
					chi2Plane[i][jj] = (dof > 0.0 && chi2 > 0.0) ? (float)(chi2 / dof) : 0.0f;
				}
				/*  vx/vy through the normal weighted-accumulation path, like any other round. */
				vxTmp[i][jj] = (float)(vx * scX);
				vyTmp[i][jj] = (float)(vy * scY);
				sxTmp[i][jj] = (float)scX;
				syTmp[i][jj] = (float)scY;
				/*  The SHARED vz/dT plane keeps its normal meaning.  Giving true3D its own
				    vZ3D plane does not mean this one stops needing valid contents: under
				    -timeOverlap it is the ".dT" band and filterDT() erases vx/vy/ex/ey wherever
				    it is out of range, so leaving it unset wipes the whole product. */
				if (outputImage->timeOverlapFlag == TRUE)
				{
					double deltaOffCenter = (sumW[i][jj] > 0.0) ? sumWT[i][jj] / sumW[i][jj] : 0.0;
					vzTmp[i][jj] = (float)(deltaOffCenter * sqrt(scX * scY));
				}
				else
				{
					/*  Not the slope-derived vz the other solvers write -- the SOLVED one, so
					    ".vz" means vertical velocity under -true3D.  ".vz3d" duplicates it so the
					    band is unambiguous even in the -timeOverlap case above. */
					vzTmp[i][jj] = (float)vz;
				}
				fScale[i][jj] = 1.0;
				nSolved++;
				nObsSum += (double)nObs[i][jj];
				nObsAzSum += (double)nObsAz[i][jj];
				if ((double)nObs[i][jj] > nObsMax)
				{
					nObsMax = (double)nObs[i][jj];
				}
			} /* j */
		}	  /* i */
		if (fl > 0)
		{
			computeScaleLS((float **)vxTmp, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0,
						   (double)(-LARGEINT), iMinAll, iMaxAll, jMinAll, jMaxAll);
		}
		redoNormalization(1.0, outputImage, iMinAll, iMaxAll, jMinAll, jMaxAll, vXimage, vYimage, vZimage,
						  errorX, errorY, scaleX, scaleY, scaleZ, fScale, vxTmp, vyTmp, vzTmp, sxTmp, syTmp, FALSE);
	}
	fprintf(stderr, "mosaicTrue3D: solved %ld  rejCond %ld  rejNoData %ld  clipped %ld  "
					"meanNobs %.2f  meanNobsAz %.2f  maxNobs %.0f\n",
			(long)nSolved, (long)nRejCond, (long)nRejNoData, (long)nClip,
			(nSolved > 0) ? nObsSum / (double)nSolved : 0.0,
			(nSolved > 0) ? nObsAzSum / (double)nSolved : 0.0, nObsMax);
	fprintf(outputImage->fpLog, "; true3D nImagesUsed : %i of %i\n", nUsed, total);
	fprintf(outputImage->fpLog, "; true3D nSolved     : %ld\n", (long)nSolved);
	fprintf(outputImage->fpLog, "; true3D rejCond     : %ld\n", (long)nRejCond);
	fprintf(outputImage->fpLog, "; true3D meanNobs    : %.3f\n", (nSolved > 0) ? nObsSum / (double)nSolved : 0.0);
	fprintf(outputImage->fpLog, "; true3D meanNobsAz  : %.3f\n", (nSolved > 0) ? nObsAzSum / (double)nSolved : 0.0);
	fprintf(outputImage->fpLog, "; true3D maxSigma    : %f\n", true3DMaxSigma);
	fprintf(outputImage->fpLog, "; true3D projectMode : %i\n", true3DProject);

	free(localImgs);
	free(jRowMin);
	free(jRowMax);
	freeDouble3D(Nxx, outputImage->ySize); freeDouble3D(Nxy, outputImage->ySize);
	freeDouble3D(Nxz, outputImage->ySize); freeDouble3D(Nyy, outputImage->ySize);
	freeDouble3D(Nyz, outputImage->ySize); freeDouble3D(Nzz, outputImage->ySize);
	freeDouble3D(bxAcc, outputImage->ySize); freeDouble3D(byAcc, outputImage->ySize);
	freeDouble3D(bzAcc, outputImage->ySize); freeDouble3D(SddAcc, outputImage->ySize);
	freeFloat3D(nObsAz, outputImage->ySize);
	if (sumW != NULL)
	{
		freeDouble3D(sumW, outputImage->ySize);
		freeDouble3D(sumWT, outputImage->ySize);
	}
	freeFloat3D(sxPlane, outputImage->ySize);
	freeFloat3D(syPlane, outputImage->ySize);
	if (outputImage->jointNObs == NULL)
	{
		outputImage->jointNObs = nObs;
		outputImage->jointChi2 = chi2Plane;
	}
	else
	{
		freeFloat3D(nObs, outputImage->ySize);
		freeFloat3D(chi2Plane, outputImage->ySize);
	}
	/* Ownership passes to outputImage; mosaic3d writes and frees them. */
	outputImage->vZ3D = vZ3D;
	outputImage->errorZ = errZ;
	if (true3DDiag == TRUE)
	{
		/*  Overwrite the (otherwise unused in diag mode) planes with the pass-split means so
		    they reach disk through the existing writer: .vz3d <- mean P_asc, .ez <- mean P_desc. */
		int32_t ii, kk;
		for (ii = 0; ii < outputImage->ySize; ii++)
		{
			for (kk = 0; kk < outputImage->xSize; kk++)
			{
				TRUE3DACC sw = sumWasc[ii][kk] + sumWdesc[ii][kk];
				if (sw > 0.0)
				{
					double mp = sumWpsi[ii][kk] / sw;
					double vp = sumWpsi2[ii][kk] / sw - mp * mp;
					vZ3D[ii][kk] = (float)(mp * RTOD);
					errZ[ii][kk] = (float)((vp > 0.0 ? sqrt(vp) : 0.0) * RTOD);
				}
				else
				{
					vZ3D[ii][kk] = (float)-LARGEINT;
					errZ[ii][kk] = (float)-LARGEINT;
				}
			}
		}
		fprintf(stderr, "mosaicTrue3D: -true3DDiag -- .vz3d holds weighted MEAN psi (deg), .ez holds weighted STD psi (deg)\n");
		freeDouble3D(sumWPasc, outputImage->ySize); freeDouble3D(sumWasc, outputImage->ySize);
		freeDouble3D(sumWPdesc, outputImage->ySize); freeDouble3D(sumWdesc, outputImage->ySize);
		freeDouble3D(sumWpsi, outputImage->ySize); freeDouble3D(sumWpsi2, outputImage->ySize);
	}
	{
		gettimeofday(&funcEnd, NULL);
		fprintf(stderr, "mosaicTrue3D total time: %.3f s\n",
				(funcEnd.tv_sec - funcStart.tv_sec) + (funcEnd.tv_usec - funcStart.tv_usec) * 1e-6);
	}
	endScale(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, FALSE);
	fprintf(outputImage->fpLog, ";\n; Returning from mosaicTrue3D(.c)\n;\n");
	fflush(outputImage->fpLog);
}

/*
  Smallest eigenvalue of a symmetric 3x3 by the closed-form trigonometric method (Smith 1961).
  No iteration and no workspace, so it is safe inside a parallel region -- unlike the cRecipes
  SVD routines, which carry the lazy-global-workspace hazard documented for svInterpBnBp.
*/
double lambdaMin3(double a11, double a12, double a13, double a22, double a23, double a33)
{
	double p1, q, p2, p, r, phi, b11, b12, b13, b22, b23, b33, detB;
	p1 = a12 * a12 + a13 * a13 + a23 * a23;
	if (p1 <= 0.0)
	{
		double m = a11;
		if (a22 < m) { m = a22; }
		if (a33 < m) { m = a33; }
		return m;
	}
	q = (a11 + a22 + a33) / 3.0;
	p2 = (a11 - q) * (a11 - q) + (a22 - q) * (a22 - q) + (a33 - q) * (a33 - q) + 2.0 * p1;
	p = sqrt(p2 / 6.0);
	if (!(p > 0.0))
	{
		return q;
	}
	b11 = (a11 - q) / p; b22 = (a22 - q) / p; b33 = (a33 - q) / p;
	b12 = a12 / p; b13 = a13 / p; b23 = a23 / p;
	detB = b11 * (b22 * b33 - b23 * b23) - b12 * (b12 * b33 - b23 * b13) + b13 * (b12 * b23 - b22 * b13);
	r = detB / 2.0;
	if (r <= -1.0) { phi = M_PI / 3.0; }
	else if (r >= 1.0) { phi = 0.0; }
	else { phi = acos(r) / 3.0; }
	/* eig1 = q + 2p cos(phi) is the largest; the smallest is at phi + 2pi/3. */
	return q + 2.0 * p * cos(phi + (2.0 * M_PI / 3.0));
}

/* Inverse of a symmetric positive-definite 3x3 by cofactors.  FALSE if not invertible. */
int32_t invSym3(double a11, double a12, double a13, double a22, double a23, double a33,
				double C[3][3], double *det)
{
	double c11, c12, c13, c22, c23, c33, d;
	c11 = a22 * a33 - a23 * a23;
	c12 = a13 * a23 - a12 * a33;
	c13 = a12 * a23 - a13 * a22;
	d = a11 * c11 + a12 * c12 + a13 * c13;
	*det = d;
	if (!(fabs(d) > 0.0))
	{
		return FALSE;
	}
	c22 = a11 * a33 - a13 * a13;
	c23 = a13 * a12 - a11 * a23;
	c33 = a11 * a22 - a12 * a12;
	C[0][0] = c11 / d; C[0][1] = c12 / d; C[0][2] = c13 / d;
	C[1][0] = c12 / d; C[1][1] = c22 / d; C[1][2] = c23 / d;
	C[2][0] = c13 / d; C[2][1] = c23 / d; C[2][2] = c33 / d;
	return TRUE;
}

static TRUE3DACC **mallocZeroDouble3D(int32_t nr, int32_t nc)
{
	TRUE3DACC **tmp;
	int32_t i, j;
	tmp = (TRUE3DACC **)malloc((size_t)nr * sizeof(TRUE3DACC *));
	if (tmp == NULL) { error("mosaicTrue3D: malloc failed for accumulator rows\n"); }
	for (i = 0; i < nr; i++)
	{
		tmp[i] = (TRUE3DACC *)malloc((size_t)nc * sizeof(TRUE3DACC));
		if (tmp[i] == NULL) { error("mosaicTrue3D: malloc failed for accumulator row %i\n", i); }
		for (j = 0; j < nc; j++) { tmp[i][j] = 0.0; }
	}
	return tmp;
}

static void freeDouble3D(TRUE3DACC **image, int32_t nr)
{
	int32_t i;
	if (image == NULL) { return; }
	for (i = 0; i < nr; i++) { free(image[i]); }
	free(image);
}

static float **mallocZeroFloat3D(int32_t nr, int32_t nc)
{
	float **tmp;
	int32_t i, j;
	tmp = mallocImage(nr, nc);
	for (i = 0; i < nr; i++)
	{
		for (j = 0; j < nc; j++) { tmp[i][j] = 0.0f; }
	}
	return tmp;
}

static void freeFloat3D(float **image, int32_t nr)
{
	int32_t i;
	if (image == NULL) { return; }
	for (i = 0; i < nr; i++) { free(image[i]); }
	free(image);
}
