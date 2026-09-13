#include "stdio.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "cRecipes/nrutil.h"
#include "mosaicSource/common/common.h"
#include "mosaic3d.h"
#include <time.h>

static double now()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Non-static: speckleTrackMosaicJoint.c applies the same clip in its solve pass.
   Prototyped in mosaic3d.h so the two cannot drift apart. */
int clipVel(float x, float y, float vx, float vy, referenceVelocity *refVel);
/*
  Pure speckle tracking solution.

  Ionosphere correction:
    If rparams recorded a ;* offsetCorrectionFile line in the baseline file,
    readOffsets loads the named GeoTIFF (written by estimateIonosphere.py on
    the native ROFF grid in SLC pixels) into offsets.rOffCorrection.  The
    correction file is linked to range.offsets.vrt via the VRT metadata key
    ionosphereRangeOffsetCorrection, stamped there by SetupNISAR.py.
    At each output pixel the correction is interpolated in SLC pixel coords
    and subtracted from the metre-domain range offset (correction × SLC pixel
    size in metres).
*/
void speckleTrackMosaic(inputImageStructure *images, vhParams *params, outputImageStructure *outputImage, float fl,
						referenceVelocity *refVel, int statsFlag)
{
	extern int HemiSphere;
	extern double Rotation;
	double lat, lon;
	double ddum1, ddum2;
	vhParams *currentParams;
	conversionDataStructure *cP;
	inputImageStructure *currentImage;
	xyDEM *vCorrect;
	double ReH, Re;
	float **vXimage, **vYimage, **vZimage;
	float **scaleX, **scaleY, **scaleZ;
	double sigmaR, sigmaA; /* Sigmas for offsets */
	double dzda, dzdr;	   /* Slopes in azimuth and range */
	double range, azimuth; /* range,azimuth coords */
	double Range;		   /* Slant Range */
	double x, y, zSp, zWGS84;
	double thetaC, thetaD; /* Center look angle and deviation from center */
	double theta;		   /* Look angle */
	double psi, cotanpsi;  /* Incidence angle and its cotan */
	double hAngle;		   /* Heading angle */
	double va, vr;		   /* Components of vel in azimuth and range */
	double xyAngle;		   /* Angle from north */
	double scaleDr;		   /* Scale factor for range displacement */
	double er, ea, ex, ey; /* Relative errors */
	double demError;
	double sig2Base, sig2Off;
	ShelfMask *shelfMask; /* Mask with shelf and grounding zone */
	unsigned char sMask;
	float azSLPixSize, rSLPixSize;
	uint32_t noData;
	double scX, scY;
	float **vxTmp, **vyTmp, **vzTmp, **fScale, **sxTmp, **syTmp;
	float **errorX, **errorY;
	float da, dr, ionCorrection;
	float azimuthMin, azimuthMax;
	double vx, vy, vz, dzdtSubmergence;
	double dzdx, dzdy;
	double tCenter, tOffCenter, deltaOffCenter;
	int32_t i, j;
	int iMin, iMax, jMin, jMax;
	int32_t *jRowMin, *jRowMax; /* per-row tightened column bounds, see getRowBounds() */
	int count, total; /* Current image counter - info only */
	int validData;
	int32_t drValid, daValid, goodPixel; /* independent range/azimuth offset validity,
	                                         RA diagnostic mode only -- see outputRAFlag below */

	fprintf(stderr, "**** SPECKLE TRACKING SOLUTION ****\n");
	fprintf(outputImage->fpLog, ";\n; Entering speckleTrackMosaic(.c)\n;\n");

	shelfMask = outputImage->shelfMask;
	vCorrect = outputImage->verticalCorrection;
	/* Count total */
	total = 0;
	for (currentImage = images; currentImage != NULL; currentImage = currentImage->next)
		total++;
	/*
	  Pointers to output images
	*/
	if (outputImage->singleImageFastPath == TRUE)
	{
		/* mallocOutputImage() left the multi-image accumulation buffers NULL for
		   this path -- only vxTmp/vyTmp/vzTmp/errorX/errorY exist, and the
		   per-pixel loop below writes final values directly into them (no
		   feathering/weighting/redoNormalization/endScale needed with one
		   contributing image). */
		vXimage = NULL; vYimage = NULL; vZimage = NULL;
		scaleX = NULL; scaleY = NULL; scaleZ = NULL;
		sxTmp = NULL; syTmp = NULL; fScale = NULL;
		vxTmp = outputImage->vxTmp;
		vyTmp = outputImage->vyTmp;
		vzTmp = outputImage->vzTmp;
		errorX = outputImage->errorX;
		errorY = outputImage->errorY;
	}
	else
	{
		setupBuffers(outputImage, &vXimage, &vYimage, &vZimage, &scaleX, &scaleY, &scaleZ,
					 &vxTmp, &vyTmp, &vzTmp, &sxTmp, &syTmp, &fScale, &errorX, &errorY);
		/*
		  Compute feather scale for existing ,Then undo the prior normalization so errors are all weighted.
		*/
		computeScale((float **)vXimage, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0, (double)(-LARGEINT));
		undoNormalization(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, fScale, statsFlag);
	}
	/*
	  Loop over images
	*/
	count = 1;
	currentParams = params;
	tCenter = (outputImage->jd1 + outputImage->jd2 + 1.) * 0.5; /* Added 1 on Dec 1 to avoid .5 day bias */
	int nthreads = omp_get_max_threads();
	inputImageStructure *localImgs = (inputImageStructure *)malloc(
		(size_t)nthreads * sizeof(inputImageStructure));
	if (localImgs == NULL)
		error("speckleTrackMosaic: malloc failed for per-thread image copies\n");
	jRowMin = (int32_t *)malloc((size_t)outputImage->ySize * sizeof(int32_t));
	jRowMax = (int32_t *)malloc((size_t)outputImage->ySize * sizeof(int32_t));
	if (jRowMin == NULL || jRowMax == NULL)
		error("speckleTrackMosaic: malloc failed for jRowMin/jRowMax\n");
	for (currentImage = images; currentImage != NULL; currentImage = currentImage->next)
	{   fprintf(stderr, "Adding image %i of %i: %s\n", count, total, currentParams->offsets.file);
		/* Compute central time, and delta from nominal*/
		tOffCenter = currentImage->julDay + currentParams->nDays * 0.5;
		deltaOffCenter = tOffCenter - tCenter;
		/* Error check weight */
	
		if (fabs(currentImage->weight - 1.0) > 0.01 && outputImage->timeOverlapFlag == FALSE)
			error("non unity weight, but overlap flag not set\n");
		if (currentParams->offsets.rFile == NULL || currentImage->weight < 0.00001)
		{
			currentParams = currentParams->next;
			continue; /* Skip if no data*/
		}
		fprintf(stderr, "\033[1;31mRIGHT(+)/LEFT(-): %i of %i\033[0m %s\n", count, total, currentParams->offsets.file);
		count++;
		/*
		  Conversions initialization and get bounding box
		*/
		cP = setupGeoConversions(currentImage, &azSLPixSize, &rSLPixSize, &Re, &ReH, &thetaC, &ddum1, &ddum2);
		getRegion(currentImage, &iMin, &iMax, &jMin, &jMax, outputImage);
		double t0 = now();
		getAzimuthBoundsForXYBox(iMin, iMax, jMin, jMax, currentImage, outputImage, &azimuthMin, &azimuthMax);
		double t1 = now();
		//azimuthMin =0; azimuthMax= LARGEINT;
		fprintf(stderr, "Time to get azimuth bounds: %f seconds\n", t1-t0);
		fprintf(stderr, "\033[35mRegion i %i to %i j %i to %i azimuth range %f to %f\033[0m	\n", iMin, iMax, jMin, jMax, azimuthMin, azimuthMax);
		/*
		  Read Offset
		*/
		if (iMin <= iMax && jMin <= jMax)
		{
			double t0 = now();
			/*  GUARDED READ -- see mosaicHopper.c: a corrupt input for this product is
			    recorded and skipped instead of exiting.  Disarmed before the parallel region. */
			jmp_buf readJmp;
			if (setjmp(readJmp) != 0)
			{
				errorRecoveryJmp = NULL;
				recordFailedProduct(currentParams->offsets.rFile, errorRecoveryMsg);
				currentParams = currentParams->next;
				continue;
			}
			errorRecoveryJmp = &readJmp;
			readOffsetDataAndParams(&(currentParams->offsets), azimuthMin, azimuthMax);
			double t1 = now();
			fprintf(stderr, "Time read offsets: %f seconds\n", t1-t0);
			/* Prime lazy-init routines in the serial section so threads never race on
			   shared workspace (uRows/uBuffer for svInitAzParams; bnS/bpS malloc in
			   svInitBnBp). After these calls azInit==TRUE and bnS/bpS are allocated. */
			if (currentParams->offsets.deltaB != DELTABNONE) {
				double bnS, bpS;
				svAzOffset(currentImage, &(currentParams->offsets), 0.0, 0.0);
				svInterpBnBp(currentImage, &(currentParams->offsets), 0.0, &bnS, &bpS);
			}
			errorRecoveryJmp = NULL;
		}
		else
		{
			iMax = iMin - 1;
			jMax = jMin - 1;
		}
		/* Tighten the per-row column range to the swath's actual (possibly diagonal)
		   footprint instead of getRegion()'s single axis-aligned bbox — see
		   common/getRegion.c:getRowBounds() for the geometry. Falls back to [jMin,jMax)
		   on any row it can't tighten, so this only ever saves work, never drops data. */
		if (iMax > iMin && jMax > jMin)
			getRowBounds(currentImage, outputImage, iMin, iMax, jMin, jMax, jRowMin, jRowMax);
		if(currentParams->offsets.sigmaAresidual > outputImage->sigmaAThresh)
		{
			fprintf(stderr, "Skipping sigmaAresidual > sigmaAThresh: %f > %f\n",
				currentParams->offsets.sigmaAresidual, outputImage->sigmaAThresh);
			currentParams = currentParams->next;
			continue;
		}
		if(currentParams->offsets.sigmaAresidual < 0)
		{
			/* azparams found no solution (sigma<0 sentinel -- see fewPointsAz() in
			   computeAzparams.c); not caught by the sigmaAThresh check above since
			   negative is never > a positive threshold. */
			fprintf(stderr, "Skipping %s: no azimuth baseline solution (sigmaAresidual<0)\n",
				currentParams->offsets.azParamsFile);
			currentParams = currentParams->next;
			continue;
		}
		if(currentParams->offsets.sigmaRresidual < 0)
		{
			/* rparams found no solution (sigma<0 sentinel -- see fewPoints() in
			   computeRParams.c).  Speckle tracking needs both components: the range
			   and azimuth offsets of one image are solved into a single velocity
			   here, so a missing range baseline is as disqualifying as a missing
			   azimuth one and the image is skipped rather than mosaicked with the
			   zeroed baseline the sentinel block carries.  make3DOffsets is the
			   other case -- it forms velocity from the range offsets of two
			   crossing passes and never reads azparams, so there a missing azimuth
			   solution costs nothing and only this same range check applies. */
			fprintf(stderr, "Skipping %s: no range baseline solution (sigmaRresidual<0)\n",
				currentParams->offsets.rParamsFile);
			currentParams = currentParams->next;
			continue;
		}
		/*
		  Now loop over output grid
		*/
		da = 0.0;
		/* Refresh per-thread copies with current image's warm-start state. */
		{
			{
				int t;
				for (t = 0; t < nthreads; t++)
					localImgs[t] = *currentImage;
			}
#pragma omp parallel \
			private(j, x, y, lat, lon, zSp, zWGS84, dzda, dzdr, \
			        range, azimuth, Range, theta, thetaD, psi, cotanpsi, ReH, \
			        hAngle, va, vr, xyAngle, scaleDr, \
			        sigmaA, sigmaR, sig2Base, sig2Off, demError, \
			        vx, vy, vz, dzdtSubmergence, dzdx, dzdy, \
			        ex, ey, er, ea, scX, scY, \
			        da, dr, ionCorrection, sMask, noData, validData, \
			        drValid, daValid, goodPixel)
			{
				inputImageStructure *myImg = &localImgs[omp_get_thread_num()];
#pragma omp for schedule(dynamic, 8)
				for (i = iMin; i < iMax; i++)
				{
					if ((i % 100) == 0) {
						fprintf(stderr, "-- %i  %f \n", i, myImg->weight);
					}
					y = (outputImage->originY + i * outputImage->deltaY) * MTOKM;
					for (j = jRowMin[i]; j < jRowMax[i]; j++)
					{
						/*
						  Convert x/y stereographic coords to lat/lon
						*/
						x = (outputImage->originX + j * outputImage->deltaX) * MTOKM;
						xytoll1(x, y, HemiSphere, &lat, &lon, Rotation, outputImage->slat);
						/*
						  Get slope and elevation
						*/
						xyGetZandSlope(lat, lon, x, y, &zSp, &zWGS84, &dzda, &dzdr, cP, currentParams, myImg);
						validData = FALSE;
						/* Validate on zWGS84 (true elevation), not zSp -- zSp additionally
						   carries a latitude-dependent spherical-earth correction
						   (earthRadius(lat)*KMTOM - cP->Re, xyGetZandSlope.c) that's several
						   km for tracks spanning a wide latitude range relative to their
						   reference radius, and was spuriously failing this elevation sanity
						   check for real, low-elevation ice-sheet pixels. zSp itself is still
						   correct and still used below (geometryInfo() needs cP->Re + zSp) --
						   only the validity gate was wrong. Matches make3DMosaic.c/
						   make3DOffsets.c, which already validate zWGS84 directly. */
						if (zWGS84 > (MINELEVATION + 1) && zWGS84 < 10000.0)
						{ /* If valid z ....*/
							llToImageNew(lat, lon, zWGS84, &range, &azimuth, myImg);
							/* Note use theta c fixed, which is referenced to baseline */
							geometryInfo(cP, myImg, azimuth, range, zSp, thetaC, &ReH, &Range, &theta, &thetaD, &psi, zSp);
							cotanpsi = 1.0 / tan(psi);
							// Get azimuth and range components from the offset field. Note these values come back as meters
							da = interpAzOffset(range, azimuth, &(currentParams->offsets), myImg, Range, theta, azSLPixSize);
/* azimuth ionosphere: no-op unless the az fit recorded one and -useAzIonosphere is set */
if (da > -0.98 * LARGEINT)
	da += azIonCorrectionMeters(&(currentParams->offsets), myImg, range, azimuth, azSLPixSize);
							dr = interpRangeOffsetInMeters(range, azimuth, &(currentParams->offsets), myImg, Range, thetaD, rSLPixSize, theta, &demError);
							if (currentParams->offsets.rOffCorrection.rangeOffsetCorrection != NULL)
							{
								double rangeSLC, azimuthSLC;
								computeSLCFromMLCoords(myImg, range, azimuth, &rangeSLC, &azimuthSLC);
								ionCorrection = interpolateOffsetIonCorrectionInPixels(&(currentParams->offsets.rOffCorrection),
																					   rangeSLC, azimuthSLC, -LARGEINT, 0.0);
							}
							else
								ionCorrection = 0.0;
							if (ionCorrection > -0.98 * LARGEINT)
							{
								dr += ionCorrection * rSLPixSize;
							}
							/*
							  If shelf mask, get mask value
							*/
							sMask = GROUNDED;
							if (shelfMask != NULL)
								sMask = getShelfMask(shelfMask, x, y);
							if (sMask == NOSOLUTION)
							{
								da = -LARGEINT;
								dr = -LARGEINT;
							};
							/*
							  Process only good points. In RA diagnostic mode (outputRAFlag),
							  range and azimuth offsets are validated independently -- a bad
							  azimuth offset should not discard an otherwise-good range offset
							  (and vice versa), since vr and va are reported as separate
							  components, not combined into a single xy vector. xy mode still
							  requires both, since rotating to map coordinates needs both
							  components together.
							*/
							drValid = (fabs(dr) < 13.0E4);
							daValid = (fabs(da) < 10.0e4);
							/* The independent-validity relaxation below is restricted to
							   singleImageFastPath: vxTmp/vyTmp are written directly as final
							   output there (no accumulation). In the general multi-image path,
							   redoNormalization() (common/scalingFunctions.c) gates accumulation
							   of BOTH vx and vy on vxTmp alone, so a vx-only-valid pixel would
							   silently pull a stale/sentinel vyTmp into vYimage -- not safe
							   without reworking that shared function (also used by
							   make3DMosaic.c/makeVhMosaic.c/make3DOffsets.c). */
							if (outputImage->outputRAFlag && outputImage->singleImageFastPath)
								goodPixel = (drValid || daValid) && sMask != GROUNDINGZONE && (!(sMask == SHELF && outputImage->noTide == TRUE));
							else
								goodPixel = drValid && daValid && sMask != GROUNDINGZONE && (!(sMask == SHELF && outputImage->noTide == TRUE));
							if (goodPixel)
							{
								/* Compute sigma for velocity error estimate. Note that the sigmas come back as meters;
								/* Moved inside of if statement 3/1/16 */
								sigmaA = interpAzSigma(range, azimuth, &(currentParams->offsets), myImg, Range, theta, azSLPixSize);
								sig2Off = computeSig2AzParam(sin(theta), cos(theta), azimuth, Range, myImg, &(currentParams->offsets));
								sigmaA = sqrt(sigmaA * sigmaA + sig2Off);
								sigmaR = interpRangeSigma(range, azimuth, &(currentParams->offsets), myImg, Range, thetaD, rSLPixSize);
								sig2Base = computeSig2Base(sin(thetaD), cos(thetaD), azimuth, myImg, &(currentParams->offsets));
								/* Add the rparams tie-point fit residual in quadrature -- see the
								   matching comment in make3DOffsets.c.  The .sr band alone
								   describes local matching noise and cannot represent
								   long-wavelength error. */
								sigmaR = sqrt(sigmaR * sigmaR + demError * demError + sig2Base +
											  rangeAccuracyVar(&(currentParams->offsets)));
								/*
								  SHELF MASK CORRECTION HERE
								*/
								dr -= shelfMaskCorrection(myImg, currentParams, sMask, x, y, psi, &sigmaR);
								if (vCorrect != NULL)
								{
									dzdtSubmergence = interpVCorrect(x, y, vCorrect);
									dr -= -dzdtSubmergence * cos(psi) * (double)currentParams->nDays / 365.25;
								}
								/*
								   Compute velocity
								*/
								hAngle = computeHeading(lat, lon, 0., myImg, cP);
								/*
								   Compute flow direction in xy coords from dem and angle of x from north
								*/
								computeXYangle(lat, lon, &xyAngle, currentParams->xydem);
								/*
								  Note va for left sign flip done in azOffset
								*/
								scaleDr = (365.25 / (double)currentParams->nDays) * (1.0 / sin(psi));
								va = da * (365.25 / (double)currentParams->nDays);
								/*
								   Turn off slope correction for shelf to avoid shelf front or rift artifacts for now this is the default (as of 10/13/17)
								   slopes on shelves, should be small (especially relative to the 3% quoted error.
								*/
								if (sMask == SHELF)
								{
									vr = (dr * scaleDr + va * cotanpsi * 0.0) / (1.0 - cotanpsi * 0.0);
								}
								else if (outputImage->outputRAFlag && !daValid)
								{
									/* da failed its sanity check -- va is untrustworthy, so drop the
									   cross term (RA diagnostic mode only; xy mode never reaches here
									   since goodPixel already required both valid above). dzdr itself
									   comes from the DEM, not from da, so the denominator is unaffected. */
									vr = (dr * scaleDr) / (1.0 - cotanpsi * dzdr);
								}
								else
								{
									vr = (dr * scaleDr + va * cotanpsi * dzda) / (1.0 - cotanpsi * dzdr);
								}
								ea = sigmaA * (365.25 / (double)currentParams->nDays);
								/* If pixel already done and azimuth offsets used,
								   assume azimuth offsets have already been used so multiply sqrt(2)
								   to avoid double averaging.
								   This only applies if phase is being used too.
								*/
								if (outputImage->noVhFlag == FALSE && vXimage[i][j] > (-LARGEINT + 1))
									ea *= 1.41421;
								er = sigmaR * scaleDr;
								/*
								   Rotate velocity to xy coordinates, or keep as range/azimuth
								*/
								if (outputImage->outputRAFlag)
								{
									/* Independent per-component validity -- see goodPixel above. A
									   pixel can reach here with only one of dr/da valid; the other
									   component's slot gets the standard -LARGEINT no-data sentinel. */
									vx = drValid ? vr : (double)-LARGEINT;
									vy = daValid ? va : (double)-LARGEINT;
									dzdx = dzdr;
									dzdy = dzda;
									ex = drValid ? er * er : (double)-LARGEINT;
									ey = daValid ? ea * ea : (double)-LARGEINT;
								}
								else
								{
									rotateFlowDirectionToXY(vr, va, &vx, &vy, xyAngle, hAngle);
									rotateFlowDirectionToXY(dzdr, dzda, &dzdx, &dzdy, xyAngle, hAngle);
									/* NOTE THIS RETURNS VARIANCES */
									errorsToXY(er, ea, &ex, &ey, xyAngle, hAngle);
								}
								/*
								   Clip data
								*/
								noData = FALSE;
								if (refVel->clipFlag == TRUE)
									noData = clipVel(x, y, vx, vy, refVel);
								/*
								  Compute vertical velocity. Needs both vx and vy -- in the
								  independent-validity RA fast path one of them may be the
								  -LARGEINT sentinel, in which case vz isn't meaningful either.
								  (drValid && daValid is unconditionally true on every other path,
								  since goodPixel already required both there.)
								*/
								vz = (drValid && daValid) ? (vx * dzdx + vy * dzdy) : (double)-LARGEINT;

								/*vx=dzdx; vy=dzdy;*/
								if (noData == FALSE)
								{
									currentImage->used = TRUE;
									if (statsFlag == FALSE)
									{
										scX = 1.0 / (ex);
										scY = 1.0 / (ey);
									}
									else
									{
										scX = 1.0;
										scY = 1.0;
									}
									if (outputImage->singleImageFastPath == TRUE)
									{
										/* No weighting/averaging to do with one contributing
										   image -- vx/ex are already the final values (see
										   mallocOutputImage()'s singleImageFastPath comment). */
										vxTmp[i][j] = (float)vx;
										vyTmp[i][j] = (float)vy;
										errorX[i][j] = (float)ex;
										errorY[i][j] = (float)ey;
									}
									else
									{
										vxTmp[i][j] = (float)vx * scX;
										vyTmp[i][j] = (float)vy * scY;
										fScale[i][j] = 1; /* Value for zero feathering */
									}
									validData = TRUE;
									/*
									   If overlap flag = true, then use the vz buff for deltaT
									*/
									if (outputImage->makeTies == TRUE)
									{
										vzTmp[i][j] = vz;
									}
									else if (outputImage->timeOverlapFlag == TRUE)
									{
										vzTmp[i][j] = (float)(deltaOffCenter * sqrt(scX * scY));
									}
									else if (statsFlag == FALSE)
									{
										if (outputImage->vzFlag == VZDEFAULT)
											vzTmp[i][j] = (float)vz; /* vz ; */
										else if (outputImage->vzFlag == VZHORIZONTAL)
										{
											vzTmp[i][j] = (float)(dr * scaleDr);
										} /* scaled by sin(psi) for h */
										else if (outputImage->vzFlag == VZVERTICAL)
										{
											vzTmp[i][j] = (float)(dr * scaleDr * sin(psi) / cos(psi));
										} /* undo h by * sin, then make vert /cos */
										else if (outputImage->vzFlag == VZLOS)
										{
											vzTmp[i][j] = (float)(dr * scaleDr * sin(psi));
										} /* undo h by * sin */
										else if (outputImage->vzFlag == VZINC)
										{
											vzTmp[i][j] = (float)psi * RTOD;
										} /* undo h by * sin */
									}
									else
									{
										vzTmp[i][j] = 1.0;
									}
									if (outputImage->singleImageFastPath == FALSE)
									{
										sxTmp[i][j] = scX;
										syTmp[i][j] = scY;
									}
								}
							} /* end fabs(dr) < 13.0E4 && fabs(da)... */
							/* Write incidence angle for all valid-elevation pixels */
							if (outputImage->vzFlag == VZINC)
								vzTmp[i][j] = (float)psi * RTOD;
						} /* end if valid z */
						/* Mark as no data if not valid data */
						if (validData == FALSE)
						{
							if (outputImage->singleImageFastPath == TRUE)
							{
								/* No later gated accumulation pass exists to keep
								   vyTmp/vzTmp/errorX/errorY quarantined for invalid
								   pixels the way redoNormalization()'s vxTmp-gate does
								   in the normal path -- reset all five explicitly. */
								vxTmp[i][j] = (float)-LARGEINT;
								vyTmp[i][j] = (float)-LARGEINT;
								vzTmp[i][j] = (float)-LARGEINT;
								errorX[i][j] = (float)-LARGEINT;
								errorY[i][j] = (float)-LARGEINT;
							}
							else
							{
								vxTmp[i][j] = (float)-LARGEINT;
								fScale[i][j] = 0.0;
							}
						}
					} /* j loop */
				}	  /* i loop */
			} /* End omp parallel */
		}

		if (outputImage->singleImageFastPath == FALSE)
		{
			/*
			  Compute scale array for feathering.
			*/
			if (fl > 0 && iMax > 0 && jMax > 0 && statsFlag == FALSE)
				computeScale((float **)vxTmp, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0, (double)(-LARGEINT));
			/*
			  Now sum current result. Falls through if no intersection (iMax&jMax==0)
			*/
			redoNormalization(currentImage->weight, outputImage, iMin, iMax, jMin, jMax, vXimage, vYimage, vZimage, errorX, errorY,
							  scaleX, scaleY, scaleZ, fScale, vxTmp, vyTmp, vzTmp, sxTmp, syTmp, statsFlag);
		}
		/* singleImageFastPath: vxTmp/vyTmp/vzTmp/errorX/errorY were already written as
		   final values directly in the per-pixel loop above -- nothing to accumulate. */
		/*  Update image pointer to move to  next image*/
		currentParams = currentParams->next;
	} /* End image loop */
	free(localImgs);
	free(jRowMin);
	free(jRowMax);
	/*   ********************END OF MAIN LOOP *****************************
		Adjust scale
	*/
	if (outputImage->singleImageFastPath == FALSE)
	{
		endScale(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, statsFlag);
	}
	else
	{
		/* vxTmp/vyTmp/vzTmp/errorX/errorY already hold final values (written directly
		   in the per-pixel loop) -- alias them in as image/image2/image3 so the rest
		   of mosaic3d.c (GeoTIFF/binary writer) needs no changes at all. */
		outputImage->image = (void **)outputImage->vxTmp;
		outputImage->image2 = (void **)outputImage->vyTmp;
		outputImage->image3 = (void **)outputImage->vzTmp;
	}
	fprintf(outputImage->fpLog, ";\n; Returning from speckleTrackMosaic(.c)\n;\n");
}

/*********************************************************************************************************
Input baseline info.
**********************************************************************************************************/

/*
   if difference between velocity map and reference map exceeds some value, return noData=TRUE
   Only apply to speeds < 100 m/yr. This option has not been used except for special cases.
*/
/*
  Return-checked reference-velocity clip, for the hopper solvers.

  Same rule as clipVel() below -- reject when the solved velocity differs from the reference by
  more than refVel->clipThresh, and only where one of the two is slower than 100 m/yr, so fast
  ice is never clipped against a reference that may simply be from a different epoch.

  It differs from clipVel() in one respect, deliberately: it honours refVelInterp()'s return
  value.  refVelInterp() returns FALSE without writing its outputs when the pixel is outside the
  reference grid or the reference is no-data, so a caller that ignores the return compares
  against uninitialised stack.  Here that case means "no reference to judge against", and the
  pixel is KEPT.  clipVel() is left as-is so the legacy path is bit-for-bit unchanged.

  Pure: reads only refVel, writes nothing shared.  Safe to call from inside an OpenMP region.
*/
int clipVelChecked(double x, double y, double vx, double vy, referenceVelocity *refVel)
{
	float vxPt, vyPt, exPt, eyPt;
	double dv, refSpeed, newSpeed;

	if (refVel == NULL || refVel->clipFlag != TRUE)
	{
		return (FALSE);
	}
	if (refVelInterp(x, y, refVel, &vxPt, &vyPt, &exPt, &eyPt) != TRUE)
	{
		return (FALSE); /* no reference value here -- nothing to clip against */
	}
	dv = sqrt((vx - (double)vxPt) * (vx - (double)vxPt) + (vy - (double)vyPt) * (vy - (double)vyPt));
	refSpeed = sqrt((double)vxPt * (double)vxPt + (double)vyPt * (double)vyPt);
	newSpeed = sqrt(vx * vx + vy * vy);
	if (dv > refVel->clipThresh && (refSpeed < 100.0 || newSpeed < 100.0))
	{
		return (TRUE);
	}
	return (FALSE);
}

int clipVel(float x, float y, float vx, float vy, referenceVelocity *refVel)
{
	float vxPt, vyPt, exPt, eyPt, exy;
	/*
	   Interp reference vel
	*/
	refVelInterp(x, y, refVel, &vxPt, &vyPt, &exPt, &eyPt);
	/* Compute difference */
	exy = sqrt((vx - vxPt) * (vx - vxPt) + (vy - vyPt) * (vy - vyPt));
	/* Only do if either new vel or ref vel < 100 m/yr */
	if (exy > refVel->clipThresh && (sqrt(vxPt * vxPt + vyPt * vyPt) < 100 || sqrt(vx * vx + vy * vy) < 100))
		return (TRUE);
	return (FALSE);
}

#define NODATA -2000000000
static int32_t refBinlinear(float **X, int32_t im, int32_t jm, double t, double u, float *dx)
{
	double p1, p2, p3, p4;
	p1 = X[im][jm];
	p2 = X[im][jm + 1];
	p3 = X[im + 1][jm + 1];
	p4 = X[im + 1][jm];
	/* Don't use if all 4pts aren't good - should ensure better quality data */
	if (p1 <= (NODATA + 1) || p2 <= (NODATA + 1) || p3 <= (NODATA + 1) || p4 <= (NODATA + 1))
		return (FALSE);
	*dx = (double)((1.0 - t) * (1.0 - u) * p1 + t * (1.0 - u) * p2 + t * u * p3 + (1.0 - t) * u * p4);
	return TRUE;
}

unsigned char refVelInterp(double x, double y, referenceVelocity *refVel, float *vxPt, float *vyPt, float *exPt, float *eyPt)
{
	int32_t im, jm;
	double t, u, xi, yi;

	if (refVel->velFile == NULL)
		return (FALSE); /* no refVel, so return */

	xi = ((x * KMTOM - refVel->x0) / refVel->dx + 0.5);
	yi = ((y * KMTOM - refVel->y0) / refVel->dy + 0.5);

	jm = (int32_t)xi;
	im = (int32_t)yi;

	if (jm < 0 || im < 0 || jm >= refVel->nx || im >= refVel->ny)
		return (FALSE); /* outside of bounds -> no refVel */
	if (refVel->vx[im][jm] < (NODATA + 1) || refVel->vy[im][jm] < (NODATA + 1))
		return (FALSE); /* no refVel value */
	/*
	  nearest neigbbor on border
	*/
	if (jm == (refVel->nx - 1) || im == (refVel->ny - 1))
	{
		if (refVel->vx[im][jm] < (NODATA + 1))
			return (FALSE);
		*vxPt = refVel->vx[im][jm];
		*vyPt = refVel->vy[im][jm];
		if (refVel->initMapFlag == TRUE)
		{
			*exPt = refVel->ex[im][jm];
			*eyPt = refVel->ey[im][jm];
		}
		return (TRUE);
	}

	t = (float)(xi - (double)jm);
	u = (float)(yi - (double)im);
	if (refBinlinear(refVel->vx, im, jm, t, u, vxPt) == FALSE)
		return FALSE;
	if (refBinlinear(refVel->vy, im, jm, t, u, vyPt) == FALSE)
		return FALSE;
	/* not doing initMap so errors not needed */
	if (refVel->initMapFlag == FALSE)
		return (TRUE);
	if (refBinlinear(refVel->ex, im, jm, t, u, exPt) == FALSE)
		return FALSE;
	if (refBinlinear(refVel->ey, im, jm, t, u, eyPt) == FALSE)
		return FALSE;
	return (TRUE);
}
