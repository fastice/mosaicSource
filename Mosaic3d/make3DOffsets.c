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
  Compute velocities from range/range offsets data.

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
void make3DOffsets(inputImageStructure *allImages, vhParams *aParams, xyDEM *dem, outputImageStructure *outputImage, float fl, float timeThresh)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	extern int32_t sepAscDesc;
	extern int32_t indentRegionOutput;
	inputImageStructure *aOffImage, *dOffImage; /* Nominal asc/desc image */
	vhParams *dParams;							/* Velocity params */
	conversionDataStructure *aCp, *dCp;			/* asc/desc coordinate conversion info */
	ShelfMask *shelfMask;
	xyDEM *vCorrect;
	double lat, lon, x, y, zWGS84; /* lat/lon - x,y,z coords */
	double aRange, dRange;		   /* asc/desc absolute range */
	double arange, drange;		   /* asc/desc range coordinates in image coords */
	double aAzimuth, dAzimuth;	   /* azimuth coordinates */
	double aDelta, dDelta;
	double aReH, dReH, aRe, dRe; /* Asc/desc Earth radii and radii + alt */
	double A[2][2], B[2][2];
	double aDemError, dDemError;
	;
	double aSigmaR, dSigmaR;
	double dzdx, dzdy;						   /* Slopes for computing vertical velocity and 3 d solution */
	double aP, dP, aPe, dPe;				   /* scaled offsets and offset  error */
	double aThetaC, dThetaC, aThetaD, dThetaD; /* Asc/desc center look angle, and dev from center */
	double aTheta, dTheta, aPsi, dPsi;		   /* Asc/desc look angle, inc angle */
	double aSig2Base, dSig2Base;
	double tCenter, tOffCenterA, tOffCenterD, deltaOffCenter;
	double combWeight;
	double ddum1, ddum2;
	double scX, scY;
	double aZSp, dZSp;					/* asc/desc elevations corrected to local sphere */
	double vx, vy, vz, dzdtSubmergence; /* velocity solution */
	float **vXimage, **vYimage, **vZimage, **errorX, **errorY;
	;															 /* velocity and error buffers */
	float **vxTmp, **vyTmp, **vzTmp, **fScale, **sxTmp, **syTmp; /* Temp solutions */
	float **scaleX, **scaleY, **scaleZ;							 /*  scale buffers */
	float dRSLPixSize, aRSLPixSize;
	float aIonCorrection, dIonCorrection;
	float azimuthMin, azimuthMax;
	float dum;
	int32_t iMin, iMax, jMin, jMax; /* range in pixels over which to compute solutions */
	int32_t aa, dd, nTotal;			/* Counters for asc/desc images and total numer of images*/
	int32_t validData, Aset;		/* Flags to indicate a valide solution, and A updates */
	int32_t nCrossing;				/* Count of crossing orbit pairs found for aOffImage */
	int32_t i, j;
	unsigned char sMask;
	struct timeval funcStart, funcEnd;
	gettimeofday(&funcStart, NULL);
	shelfMask = outputImage->shelfMask;
	vCorrect = outputImage->verticalCorrection;
	fprintf(stderr, "fl = %f\n", (double)fl);
	fprintf(outputImage->fpLog, ";\n; Entering make3DOffsets(.c)\n");
	A[0][0] = 0;
	A[1][0] = 0;
	A[0][1] = 0;
	A[1][1] = 0;
	/*
	  Pointers to output images
	*/
	setupBuffers(outputImage, &vXimage, &vYimage, &vZimage, &scaleX, &scaleY, &scaleZ,
				 &vxTmp, &vyTmp, &vzTmp, &sxTmp, &syTmp, &fScale, &errorX, &errorY);
	/*
	  Compute feather scale for existing
	*/
	computeScale((float **)vXimage, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0, (double)(-LARGEINT));
	if (dem->stdLat < 50 || dem->stdLat > 80)
		error("mosaic3doff invalid slat for dem");
	/*
	  Init array. This undoes the prior normalization so errors are all weighted.
	*/
	undoNormalization(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, fScale, FALSE);
	nTotal = 0;
	for (aOffImage = allImages; aOffImage != NULL; aOffImage = aOffImage->next)
	{
		nTotal++;
	}
	fprintf(stderr, "nTotal Images %i", nTotal);
	/*
	   MAIN LOOP Loop over ascending images
	*/
	aa = 0;
	tCenter = (outputImage->jd1 + outputImage->jd2 + 1.) * 0.5; /* Added 1 on Dec 1 to avoid .5 day bias */
	int nthreads = omp_get_max_threads();
	inputImageStructure *localAImgs = (inputImageStructure *)malloc(
		(size_t)nthreads * sizeof(inputImageStructure));
	inputImageStructure *localDImgs = (inputImageStructure *)malloc(
		(size_t)nthreads * sizeof(inputImageStructure));
	if (localAImgs == NULL || localDImgs == NULL)
		error("make3DOffsets: malloc failed for per-thread image copies\n");
	for (aOffImage = allImages; aOffImage != NULL; aOffImage = aOffImage->next, aParams = aParams->next)
	{
		aa++;
		/*
			Skip if no cross flag set
		*/
		tOffCenterA = aOffImage->julDay + aParams->nDays * 0.5;
		/* skip if crossFlag False or weight too small (<5%) */
		if (aOffImage->crossFlag == FALSE || aOffImage->weight < 0.05 || aParams->offsets.rFile == NULL)
			continue;
		/* Skip if no overlap or file missing */
		indentRegionOutput = FALSE;
		if (!getRegion(aOffImage, &iMin, &iMax, &jMin, &jMax, outputImage))
			continue;
		fprintf(stderr,"\033[1;34maOffImage->rangeFile %s\033[0m\n", aParams->offsets.rFile);
		/*
		  Setup conversions
		 */
		aRSLPixSize = aOffImage->rangePixelSize / aOffImage->nRangeLooks;
		/*
		  Setup conversion parameters
		*/
		aCp = setupGeoConversions(aOffImage, &dum, &aRSLPixSize, &aRe, &aReH, &aThetaC, &ddum1, &ddum2);
		// This is going to read the full data take for the outer loop image.
		getAzimuthBoundsForXYBox(iMin, iMax, jMin, jMax, aOffImage, outputImage, &azimuthMin, &azimuthMax);
		getRParams(&(aParams->offsets));
		if (aParams->offsets.sigmaRresidual < 0)
		{
			/* rparams found no solution (sigma<0 sentinel) -- set rFile = NULL so this
			   and every future iteration skip via the existing rFile==NULL check above
			   (line ~115), the same convention already used for images with no range
			   offset processing requested at all. */
			fprintf(stderr, "make3DOffsets: %s has no baseline solution (sigma<0) -- skipping\n",
			        aParams->offsets.rFile);
			aParams->offsets.rFile = NULL;
			continue;
		}
		readRangeOrRangeOffsets(&(aParams->offsets), ASCENDING, azimuthMin, azimuthMax);
		
		/*
		 **** SECOND LOOP ***
		 */
		dd = 0;
		/* All images are in a list so once an image is processed, start with rest of images in the loop */
		dParams = aParams->next; /* Start at next element below outer loop */
		indentRegionOutput = TRUE; /* tab-indent getRegion/computeSceneAlpha/etc. output for the
		                              rest of this outer iteration, so it visually groups with the
		                              crossing-pair inner loop's own tab-indented prints below */
		nCrossing = 0;
		for (dOffImage = aOffImage->next; dOffImage != NULL; dOffImage = dOffImage->next, dParams = dParams->next)
		{
			dd++;
			tOffCenterD = dOffImage->julDay + dParams->nDays * 0.5;
			if ((dOffImage->passType == aOffImage->passType && sepAscDesc == TRUE) || dOffImage->crossFlag == FALSE)
				continue;
			//error("Could not process image %i",fabs(aOffImage->julDay - dOffImage->julDay) > timeThresh);
			/*
				Check images close enough in time
			*/
			if (fabs(aOffImage->julDay - dOffImage->julDay) > timeThresh || dOffImage->weight < 0.05 || dParams->offsets.rFile == NULL)
				continue;
			//fprintf(stderr, "time JD %f %f\n", aOffImage->julDay, dOffImage->julDay);
			/* Skip if no overlap or file missing */
			if (!getRegion(dOffImage, &iMin, &iMax, &jMin, &jMax, outputImage))
				continue;
			fprintf(stderr,"\t\033[38;5;208mdOffImage->rangeFile %s\033[0m\n", dParams->offsets.rFile);
			/*
			  Get region of  possible intersection - pass if not interect
			*/
			getIntersect(dOffImage, aOffImage, &iMin, &iMax, &jMin, &jMax, outputImage);
			//fprintf(stderr, "%i %i %i %i\n", iMin, iMax, jMin, jMax); error("STOP");
			// No intersect, then skip rest of this this loop iteration.
			if (iMax == 0 && jMax == 0)
				continue;
			/*
			  Init conversion stuff
			*/
			dCp = setupGeoConversions(dOffImage, &dum, &dRSLPixSize, &dRe, &dReH, &dThetaC, &ddum1, &ddum2);
			/*
			   Compute approximate heading by sampling overlap region
			*/
			computeSceneAlpha(outputImage, aOffImage, dOffImage, aCp, dCp, dem, &iMin, &iMax, &jMin, &jMax);
			if (iMax == 0 && jMax == 0)
				continue;
			/* 2026-06-16: also skip when computeSceneAlpha returns inverted bounds (not caught by 0,0 sentinel) */
			if (iMin > iMax || jMin > jMax)
				continue;
			/*
			  Read in descending image if needed (i.e., nozero intersect).
			*/
			getAzimuthBoundsForXYBox(iMin, iMax, jMin, jMax, dOffImage, outputImage, &azimuthMin, &azimuthMax);
			/* 2026-06-16: skip read when getAzimuthBoundsForXYBox finds no intersection (azimuthMin=azimuthMax=0) */
			if (azimuthMin == 0.0f && azimuthMax == 0.0f)
				continue;
			getRParams(&(dParams->offsets));
			if (dParams->offsets.sigmaRresidual < 0)
			{
				/* see matching comment at the ASCENDING getRParams() call above */
				fprintf(stderr, "\tmake3DOffsets: %s has no baseline solution (sigma<0) -- skipping\n",
				        dParams->offsets.rFile);
				dParams->offsets.rFile = NULL;
				continue;
			}
			readRangeOrRangeOffsets(&(dParams->offsets), DESCENDING, azimuthMin, azimuthMax);

			//error("STOP");
			
			/*
			  Loop over output grid and compute velocities
			*/
			nCrossing++;
			fprintf(stderr, "\t---- \033[1mAsc %i / %i\033[0m \033[1mDes %i\033[0m\n", aa, nTotal, dd);
			/* Prime svInitBnBp in serial before threads race on bnS/bpS malloc */
			if (aParams->offsets.deltaB != DELTABNONE) {
				double bnS, bpS;
				svInterpBnBp(aOffImage, &(aParams->offsets), 0.0, &bnS, &bpS);
				svInterpBnBp(dOffImage, &(dParams->offsets), 0.0, &bnS, &bpS);
			}
			{
				int t;
				for (t = 0; t < nthreads; t++) { localAImgs[t] = *aOffImage; localDImgs[t] = *dOffImage; }
			}
#pragma omp parallel \
			private(j, x, y, lat, lon, zWGS84, \
			        aZSp, dZSp, arange, drange, aAzimuth, dAzimuth, \
			        aDelta, dDelta, aReH, dReH, aRange, dRange, \
			        aTheta, dTheta, aThetaD, dThetaD, aPsi, dPsi, \
			        aSig2Base, dSig2Base, aSigmaR, dSigmaR, \
			        aDemError, dDemError, aIonCorrection, dIonCorrection, \
			        aP, dP, aPe, dPe, vx, vy, vz, scX, scY, \
			        dzdx, dzdy, dzdtSubmergence, deltaOffCenter, sMask, validData)
			{
				int myThread = omp_get_thread_num();
				inputImageStructure *myAImg = &localAImgs[myThread];
				inputImageStructure *myDImg = &localDImgs[myThread];
				double A[2][2], B[2][2];
				int32_t rowAset;
#pragma omp for schedule(dynamic, 8)
				for (i = iMin; i < iMax; i++)
				{
					rowAset = FALSE; /* force a fresh A on the first valid-data pixel of every
					                    row -- j==jMin alone isn't enough since that column's
					                    data may itself be invalid, leaving A never refreshed */
					y = (outputImage->originY + i * outputImage->deltaY) * MTOKM;
					for (j = jMin; j < jMax; j++)
					{
						/*
						  Convert x/y stereographic coords to lat/lon
						*/
						x = (outputImage->originX + j * outputImage->deltaX) * MTOKM;
						xytoll1(x, y, HemiSphere, &lat, &lon, Rotation, dem->stdLat);
						zWGS84 = getXYHeight(lat, lon, dem, 0.0, ELLIPSOIDAL);
						validData = FALSE;
						/*
						   Process points where elevation is known
						*/
						if (zWGS84 > MINELEVATION)
						{
							/*  Convert elevations to spherical reference	*/
							aZSp = sphericalElev(zWGS84, lat, aRe);
							dZSp = sphericalElev(zWGS84, lat, dRe);
							/*
							  Compute range azimuth positions
							*/
							llToImageNew(lat, lon, zWGS84, &arange, &aAzimuth, myAImg);
							geometryInfo(aCp, myAImg, aAzimuth, arange, aZSp, aThetaC, &aReH, &aRange, &aTheta, &aThetaD, &aPsi, aZSp);
							llToImageNew(lat, lon, zWGS84, &drange, &dAzimuth, myDImg);
							geometryInfo(dCp, myDImg, dAzimuth, drange, dZSp, dThetaC, &dReH, &dRange, &dTheta, &dThetaD, &dPsi, dZSp);
							/*  Interpolate range offsets */
							dDelta = interpRangeOffsetInMeters(drange, dAzimuth, &(dParams->offsets), myDImg, dRange, dThetaD, dRSLPixSize, dTheta, &dDemError);
							aDelta = interpRangeOffsetInMeters(arange, aAzimuth, &(aParams->offsets), myAImg, aRange, aThetaD, aRSLPixSize, aTheta, &aDemError);
							{
								double dRangeSLC, dAzimuthSLC;
								computeSLCFromMLCoords(myDImg, drange, dAzimuth, &dRangeSLC, &dAzimuthSLC);
								if (dParams->offsets.rOffCorrection.rangeOffsetCorrection != NULL)
									dIonCorrection = interpolateOffsetIonCorrectionInPixels(&(dParams->offsets.rOffCorrection), dRangeSLC, dAzimuthSLC, -LARGEINT, 0.0);
								else
									dIonCorrection = 0.0;
							}
							if (dIonCorrection > -0.98 * LARGEINT && dDelta > -LARGEINT)
								dDelta += dIonCorrection * dRSLPixSize;
							{
								double aRangeSLC, aAzimuthSLC;
								computeSLCFromMLCoords(myAImg, arange, aAzimuth, &aRangeSLC, &aAzimuthSLC);
								if (aParams->offsets.rOffCorrection.rangeOffsetCorrection != NULL)
									aIonCorrection = interpolateOffsetIonCorrectionInPixels(&(aParams->offsets.rOffCorrection), aRangeSLC, aAzimuthSLC, -LARGEINT, 0.0);
								else
									aIonCorrection = 0.0;
							}
							if (aIonCorrection > -0.98 * LARGEINT && aDelta > -LARGEINT)
								aDelta += aIonCorrection * aRSLPixSize;
						}
						else
						{
							aDelta = -LARGEINT;
							dDelta = -LARGEINT;
						} /* End if (z >... */
						/*
						  If shelf mask, get mask value
						*/
						sMask = GROUNDED;
						if (shelfMask != NULL)
							sMask = getShelfMask(shelfMask, x, y);
						if (sMask == NOSOLUTION)
						{
							aDelta = -LARGEINT;
							dDelta = -LARGEINT;
						};
						/*
						  If there is valid offsets data from both images then compute velocity
						*/
						if (aDelta > (-LARGEINT + 1) && dDelta > (-LARGEINT + 1) && zWGS84 > MINELEVATION && sMask != GROUNDINGZONE && (!(sMask == SHELF && outputImage->noTide == TRUE)))
						{
							/*
							  Compute error due to baseline
							*/
							aSig2Base = computeSig2Base(sin(aThetaD), cos(aThetaD), aAzimuth, myAImg, &(aParams->offsets));
							dSig2Base = computeSig2Base(sin(dThetaD), cos(dThetaD), dAzimuth, myDImg, &(dParams->offsets));
							aSigmaR = interpRangeSigma(arange, aAzimuth, &(aParams->offsets), myAImg, aRange, aThetaD, aRSLPixSize);
							dSigmaR = interpRangeSigma(drange, dAzimuth, &(dParams->offsets), myDImg, dRange, dThetaD, dRSLPixSize);
							aSigmaR = sqrt(aSigmaR * aSigmaR + aDemError * aDemError + aSig2Base);
							dSigmaR = sqrt(dSigmaR * dSigmaR + dDemError * dDemError + dSig2Base);
							/*
							  Tide corrections
							*/
							if (sMask == SHELF)
							{
								/* Interp tide errors, set twok (last param) as 1.0 for offsets */
								interpTideError(&aSigmaR, myAImg, aParams, x, y, aPsi, 1.0);
								interpTideError(&dSigmaR, myDImg, dParams, x, y, dPsi, 1.0);
								aDelta -= -myAImg->tideCorrection * cos(aPsi) * (double)aParams->nDays / 365.25;
								dDelta -= -myDImg->tideCorrection * cos(dPsi) * (double)dParams->nDays / 365.25;
							} /* ENd if(smask... */
							if (vCorrect != NULL)
							{
								dzdtSubmergence = interpVCorrect(x, y, vCorrect);
								aDelta -= -dzdtSubmergence * cos(aPsi) * (double)aParams->nDays / 365.25;
								dDelta -= -dzdtSubmergence * cos(dPsi) * (double)dParams->nDays / 365.25;
							}
							/*  Update A every 3rd pixel; rowAset guarantees a fresh A on the first
							    valid-data pixel of every row (regardless of which chunk a thread was on
							    previously, avoiding stale/uninitialized carry-over) -- using j==jMin alone
							    isn't sufficient since that column's own data may be invalid, in which case
							    A would never get refreshed for the row. */
							if ((j % 3) == 0 || rowAset == FALSE)
							{
								/* Offsets are self-consistent regardless of squint by construction --
								   the zero-Doppler condition forces true LOS perpendicular to true
								   velocity at the assigned time, independent of squint (see
								   mosaicSource/CLAUDE.md "Squint"). So this path never applies the
								   correction, flag or no flag -- not an oversight. */
								computeA(lat, lon, x, y, myAImg, myDImg, A, FALSE);
								rowAset = TRUE;
							}
							/*
							  Only pursue solution if sufficient difference in angles for 3d solution
							*/
							if (A[0][0] != -LARGEINT)
							{
								/*  Compute B (note B is really C in the TGARS paper	*/
								computeB(x, y, zWGS84, B, &dzdx, &dzdy, aPsi, dPsi, (xyDEM *)dem);
								/*  Scale offsets for velocity computation (scale for m/yr)	*/
								aP = 365.25 * aDelta / ((double)(aParams->nDays) * sin(aPsi));
								dP = 365.25 * dDelta / ((double)(dParams->nDays) * sin(dPsi));
								aPe = 365.25 * aSigmaR / (aParams->nDays * sin(aPsi));
								dPe = 365.25 * dSigmaR / (dParams->nDays * sin(dPsi));
								/*  Compute velocity */
								computeVxy(aP, dP, aPe, dPe, A, B, &vx, &vy, &scX, &scY);
								/*  Compute vertical velocity	*/
								vz = vx * dzdx + vy * dzdy;
								/*  Update output arrays */
								vxTmp[i][j] = vx * scX;
								vyTmp[i][j] = vy * scY;
								validData = TRUE;
								if (outputImage->makeTies == TRUE)
								{
									vzTmp[i][j] = vz;
								}
								else if (outputImage->timeOverlapFlag == TRUE)
								{
									/* For lack of better option, use the average of the two data takes */
									deltaOffCenter = 0.5 * (tOffCenterA + tOffCenterD - 2.0 * tCenter);
									vzTmp[i][j] = (float)(deltaOffCenter * sqrt(scX * scY));
								}
								else
								{
									vzTmp[i][j] = vz;
								}
								sxTmp[i][j] = scX; /* This is summing up 1/sigma^2*/
								syTmp[i][j] = scY;
								fScale[i][j] = 1.0; /* Value for zero feathering */
							}
						}
						if (validData == FALSE)
						{
							vxTmp[i][j] = (float)-LARGEINT;
							fScale[i][j] = 0.0;
						}
					} /* j loop */
					if ((i % 100) == 0)
						fprintf(stderr, "\t--+ %i\n", i);
				} /* i loop */
			} /* End omp parallel */
			/*
			  Compute scale array for feathering.
			*/
			if (fl > 0 && (iMax > 0 && jMax > 0))
				computeScaleLS((float **)vxTmp, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0, (double)(-LARGEINT), iMin, iMax, jMin, jMax);
			/*
			  Now sum current result. Falls through if no intersection (iMax&jMax==0)
			*/
			if (outputImage->timeOverlapFlag == TRUE)
			{
				combWeight = sqrt(aOffImage->weight * dOffImage->weight);
				fprintf(stderr, "\t\033[1mComb weight = %lf |Ta-Td| %lf\033[0m\n", combWeight, fabs(aOffImage->julDay - dOffImage->julDay));
			}
			else
				combWeight = 1.0;

			redoNormalization(combWeight, outputImage, iMin, iMax, jMin, jMax, vXimage, vYimage, vZimage, errorX, errorY,
							  scaleX, scaleY, scaleZ, fScale, vxTmp, vyTmp, vzTmp, sxTmp, syTmp, FALSE);
		} /* End desc loop */
		if (nCrossing == 0)
			fprintf(stderr, "No crossing orbit within timeThresh days\n\n");
		else
			fprintf(stderr, "Found %i crossing orbit(s)\n\n", nCrossing);
	}	  /* End asc loop */
	free(localAImgs);
	free(localDImgs);
	/**************************END OF MAIN LOOP ******************************/
	fprintf(stderr, "Out of main loop\n");
	{
		extern double totalOffsetsIOTime;
		gettimeofday(&funcEnd, NULL);
		fprintf(stderr, "Total offsets I/O time (whole run): %.3f s\n", totalOffsetsIOTime);
		fprintf(stderr, "Total offsets processing time (whole run): %.3f s\n",
		        (funcEnd.tv_sec - funcStart.tv_sec) + (funcEnd.tv_usec - funcStart.tv_usec) * 1e-6);
	}
	/*
	  Adjust scale
	*/
	endScale(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, FALSE);
	fprintf(outputImage->fpLog, ";\n; Returning from make3DOffsets(.c)\n");
	fflush(outputImage->fpLog);
}
