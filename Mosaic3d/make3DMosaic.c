#include "stdio.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "cRecipes/nrutil.h"
#include "mosaicSource/common/common.h"
#include "mosaic3d.h"
#include "sys/time.h"

static void setBuffer(inputImageStructure *inputImage, float *buf);
static double computePhiZM3d(double *thetaD, double z, double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
							 double Range, double Re, double ReH, double ReHfixed, double thetaC, double thetaCfixedReH, double *phaseError);
static double computePhiFlatEarthM3d(double azimuth, vhParams *vhParam, inputImageStructure *phaseImage,
									 double Range, double Re, double ReHfixed, double thetaCfixedReH, double *phaseError);
/*
************************ Estimate 3D velocity from phase **************************
*/

/* convert range and time on the lat/lon on the ellipsoid */
static void RTtoLatLon(inputImageStructure *inputImage, double r, double myTime, double *lat, double *lon)
{
	int32_t n;
	stateV *sv;
	double xs, ys, zs, vsx, vsy, vsz;
	sv = &(inputImage->sv);
	n = (long int)((myTime - sv->times[1]) / (sv->deltaT) + .5);
	n = min(max(0, n - 2), sv->nState - NUSESTATE);
	polintVec(&(sv->times[n]), &(sv->x[n]), &(sv->y[n]), &(sv->z[n]), &(sv->vx[n]), &(sv->vy[n]), &(sv->vz[n]), myTime, &xs, &ys, &zs, &vsx, &vsy, &vsz);
	smlocateZD(xs * MTOKM, ys * MTOKM, zs * MTOKM, vsx * MTOKM, vsy * MTOKM, vsz * MTOKM, r * MTOKM, lat, lon, (double)(inputImage->lookDir), 0.0);
}

static void printLatLon(inputImageStructure *inputImage)
{
	int32_t i1, j1;
	double dr, dt, lat, lon;
	dr = (inputImage->par.rf - inputImage->par.rn) / (1);
	dt = (inputImage->azimuthSize * inputImage->nAzimuthLooks / inputImage->par.prf) / (1);
	fprintf(stderr, "%f %f\n", dr, dt);
	for (i1 = 0; i1 < 2; i1++)
		for (j1 = 0; j1 < 2; j1++)
		{
			fprintf(stderr, "%f %f \n", inputImage->par.rc, inputImage->cpAll.sTime);
			RTtoLatLon(inputImage, inputImage->par.rc + i1 * dr, inputImage->cpAll.sTime + j1 * dt, &lat, &lon); /* lat/lon in first image */
			fprintf(stderr, "%f %f \n", lat, lon);
		}
	RTtoLatLon(inputImage, inputImage->par.rc + 0.5 * dr, inputImage->cpAll.sTime + 0.5 * dt, &lat, &lon); /* lat/lon in first image */
	fprintf(stderr, "%f %f \n", lat, lon);
}

void make3DMosaic(inputImageStructure *ascImages, inputImageStructure *descImages,
				  vhParams *ascParams, vhParams *descParams, xyDEM *dem, outputImageStructure *outputImage, float fl, int32_t no3d, float timeThreshPhase)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	extern float *AImageBuffer, *DImageBuffer;
	extern float *AIonBuffer, *DIonBuffer;
	extern int32_t sepAscDesc;
	extern int32_t indentRegionOutput;
	conversionDataStructure *aCp, *dCp;							/* asc/desc coordinate conversion info */
	inputImageStructure *allImages, *aPhaseImage, *dPhaseImage; /* list of all images, and individual asc/desc images */
	vhParams *aParams, *dParams;
	ShelfMask *shelfMask;
	xyDEM *vCorrect;
	double lat, lon, x, y, zWGS84; /* lat/lon - x,y,z coords */
	double aRange, dRange;		   /* asc/desc absolute range */
	double arange, drange;		   /* asc/desc range coordinates in image coords */
	double aAzimuth, dAzimuth;
	double aPhase, dPhase;
	double aIonPhase, dIonPhase;
	double aReH, dReH, aRe, dRe;								   /* Asc/desc Earth radii and radii + alt */
	double aReHfixed, aThetaCfixedReH, dReHfixed, dThetaCfixedReH; /* Fixed geo params */
	double scX, scY;											   /* X,Y scale factors */
	double A[2][2], B[2][2];									   /* 3d solution matrices */
	double phaseErrorA, phaseErrorD;							   /* asc/desc phase errors */
	double dzdx, dzdy;											   /* Slopes for computing vertical velocity and 3 d solution */
	double aPhiZ, dPhiZ;										   /* Phase due to topopgraphy */
	double aP, dP, aPe, dPe;									   /* scaled phases and phase error */
	double aThetaC, dThetaC, aThetaD, dThetaD;					   /* Asc/desc center look angle, and dev from center */
	double aTheta, dTheta, aPsi, dPsi;							   /* Asc/desc look angle, inc angle */
	double tCenter, tOffCenterA, tOffCenterD, deltaOffCenter;	   /* Variables for tracking time offsets */
	double combWeight;											   /* Weight based on time overlap */
	double aZSp, dZSp;											   /* asc/desc elevations corrected to local sphere */
	double twokA, twokD;										   /* 4pi/lambda */
	double scaleA, scaleD;
	double dzdtSubmergence;
	double vx, vy, vz;		   /* velocity solution */
	float geoTolerance = 1e-3; /* Tolerance for geocoding 1e-3 should give a few meters, good enough for velocity */
	float **vXimage, **vYimage, **vZimage, **errorX, **errorY;
	;															 /* velocity and error buffers */
	float **vxTmp, **vyTmp, **vzTmp, **fScale, **sxTmp, **syTmp; /* Temp solutions */
	extern int32_t pairOverCount; /* -pairOverCount; default off, see mosaic3d.c */
	float **nAtmp, **nDtmp;   /* per-pixel outer-image and pair counts for the over-counting correction */
	float **errorX0, **errorY0; /* errorX/errorY accumulators as this round found them */
	unsigned char *aContrib;  /* flat [ySize×xSize]: did current aPhaseImage contribute at this pixel? */
	float **scaleX, **scaleY, **scaleZ;							 /*  scale buffers */
	float dum1, dum2;	
	float azimuthMin, azimuthMax;										 /* Placeholder dummys for function calls */
	int32_t validData;											 /* Flag to indicate a valid solution */
	int32_t iMin, iMax, jMin, jMax;								 /* range in pixels over which to compute solutions */
	int32_t aa, dd, nTotal;										 /* Counters for asc/desc images and total numer of images*/
	int32_t i, j, i1, j1, count;
	unsigned char sMask;
	struct timeval start, stop;
	struct timeval funcStart, funcEnd;
	gettimeofday(&funcStart, NULL);

	fprintf(outputImage->fpLog, ";\n; Entering make3DMosaic(.c)\n");
	if (no3d == TRUE)
	{
		fprintf(outputImage->fpLog, ";\n; no3d flag set, returning;\n; Returning from make3DMosaic(.c)\n");
		return;
	}
	A[0][0] = 0;
	A[1][0] = 0;
	A[0][1] = 0;
	A[1][1] = 0;
	shelfMask = outputImage->shelfMask;
	vCorrect = outputImage->verticalCorrection;
	/*
	  Pointers to output images
	*/
	setupBuffers(outputImage, &vXimage, &vYimage, &vZimage, &scaleX, &scaleY, &scaleZ,
				 &vxTmp, &vyTmp, &vzTmp, &sxTmp, &syTmp, &fScale, &errorX, &errorY);
	nAtmp = NULL;
	nDtmp = NULL;
	aContrib = NULL;
	if (pairOverCount == TRUE)
	{
		int i1, j1;
		nAtmp = mallocImage(outputImage->ySize, outputImage->xSize);
		nDtmp = mallocImage(outputImage->ySize, outputImage->xSize);
		for (i1 = 0; i1 < outputImage->ySize; i1++)
			for (j1 = 0; j1 < outputImage->xSize; j1++)
				nAtmp[i1][j1] = nDtmp[i1][j1] = 0.0f;
		aContrib = (unsigned char *)calloc(
			(size_t)outputImage->ySize * outputImage->xSize, sizeof(unsigned char));
		if (aContrib == NULL)
			error("make3DMosaic: calloc failed for aContrib\n");
	}
	if (dem->stdLat < 50 || dem->stdLat > 80)
		error("mosaic3d invalid slat for dem");
	/*
	  Compute feather scale for existing, and undo normalization
	*/
	computeScale((float **)vXimage, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0, (double)(-LARGEINT));
	undoNormalization(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, fScale, FALSE);
	/* Snapshot the error accumulators as this round inherits them, so the pair
	   over-counting correction below inflates only what this round adds -- see
	   inflatePairOverCount() in common/scalingFunctions.c. */
	if (pairOverCount == TRUE)
	{
		int i1, j1;
		errorX0 = mallocImage(outputImage->ySize, outputImage->xSize);
		errorY0 = mallocImage(outputImage->ySize, outputImage->xSize);
		for (i1 = 0; i1 < outputImage->ySize; i1++)
			for (j1 = 0; j1 < outputImage->xSize; j1++)
			{
				errorX0[i1][j1] = errorX[i1][j1];
				errorY0[i1][j1] = errorY[i1][j1];
			}
	}
	allImages = ascImages; /* Added 5/30 to avoid using asc/desc */
	aParams = ascParams;
	nTotal = 0;
	for (aPhaseImage = allImages; aPhaseImage != NULL; aPhaseImage = aPhaseImage->next)
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
		error("make3DMosaic: malloc failed for per-thread image copies\n");
	for (aPhaseImage = allImages; aPhaseImage != NULL; aPhaseImage = aPhaseImage->next, aParams = aParams->next)
	{
		aa++;
		/* Use for calculating time skew 09/21/17 */
		tOffCenterA = aPhaseImage->julDay + aParams->nDays * 0.5;
		/*
			Skip inner loop if outer image is a nophase
		*/
		if (strstr(aPhaseImage->file, "nophase") != NULL)
			continue;
		/* Skip if no overlap or file missing */
		indentRegionOutput = FALSE;
		if (!getRegion(aPhaseImage, &iMin, &iMax, &jMin, &jMax, outputImage))
			continue;
		{
			extern int32_t useSquint;
			double aSq = (useSquint && aPhaseImage->hasSquintPolynomial)
				? evaluateSquint(aPhaseImage, aPhaseImage->rangeSize * 0.5, aPhaseImage->azimuthSize * 0.5)
				: 0.0;
			if (!useSquint)
				fprintf(stderr, "\033[1;34maPhaseImage %s %3.0f -- %5.3f -- %4i squint off\033[0m\n",
						aPhaseImage->file, aParams->nDays, aPhaseImage->par.lambda, aa);
			else if (aPhaseImage->hasSquintPolynomial)
				fprintf(stderr, "\033[1;34maPhaseImage %s %3.0f -- %5.3f -- %4i squint on sq=%5.2f\033[0m\n",
						aPhaseImage->file, aParams->nDays, aPhaseImage->par.lambda, aa, aSq);
			else
				fprintf(stderr, "\033[1;34maPhaseImage %s %3.0f -- %5.3f -- %4i squint on sq=\033[31m%5.2f\033[1;34m\033[0m\n",
						aPhaseImage->file, aParams->nDays, aPhaseImage->par.lambda, aa, aSq);
		}
		/*  Set buffer, memory channel for sharedmem, and read image		*/
		setBuffer(aPhaseImage, AImageBuffer);
		/*  Setup conversion parameters		*/
		aCp = setupGeoConversions(aPhaseImage, &dum1, &dum2, &aRe, &aReH, &aThetaC, &aReHfixed, &aThetaCfixedReH);
		aPhaseImage->tolerance = geoTolerance;
		getAzimuthBoundsForXYBox(iMin, iMax, jMin, jMax, aPhaseImage, outputImage, &azimuthMin, &azimuthMax);
		getMosaicInputImage(aPhaseImage, azimuthMin, azimuthMax);
		getIonospherePhaseImage(aPhaseImage, AIonBuffer, azimuthMin, azimuthMax);
		aPhaseImage->memChan = MEM1;
		twokA = (4.0 * PI) / aPhaseImage->par.lambda;
		/*
		   SECOND LOOP: Loop over descending images:  Changed 05/30/07 to search all images below one from aPhase loop
		*/
		dd = 0;
		dParams = aParams->next; /* Start at next element below outer loop */
		indentRegionOutput = TRUE; /* tab-indent getRegion/computeSceneAlpha/etc. output for the
		                              rest of this outer iteration, so it visually groups with the
		                              crossing-pair inner loop's own tab-indented prints below */
		for (dPhaseImage = aPhaseImage->next; dPhaseImage != NULL; dPhaseImage = dPhaseImage->next, dParams = dParams->next)
		{
			dd++;
			/* Use for calculating time skew 09/21/17 */
			tOffCenterD = dPhaseImage->julDay + dParams->nDays * 0.5;
			/* Skip if to far in time */
			if (fabs(aPhaseImage->julDay - dPhaseImage->julDay) > timeThreshPhase)
				continue;
			/*   Check inner image is not a nophase and only proceed with this iteration of the inner loop if so. */
			if (strstr(dPhaseImage->file, "nophase") != NULL)
				continue;
			/*
			   Moved up 9/14/2016 to avoid init ll - had to modify getRegion to use lat/lon from geodat
			*/
			if (!getRegion(dPhaseImage, &iMin, &iMax, &jMin, &jMax, outputImage))
			{
				/* No overlap or file missing — flag as nophase for future outer loops */
				dPhaseImage->file = strdup("nophase");
				continue;
			}
			{
				extern int32_t useSquint;
				double dSq = (useSquint && dPhaseImage->hasSquintPolynomial)
					? evaluateSquint(dPhaseImage, dPhaseImage->rangeSize * 0.5, dPhaseImage->azimuthSize * 0.5)
					: 0.0;
				if (!useSquint)
					fprintf(stderr, "\t\033[38;5;208mdPhaseImage %s %3.0f -- %5.3f -- %4i squint off\033[0m\n",
							dPhaseImage->file, dParams->nDays, dPhaseImage->par.lambda, dd);
				else if (dPhaseImage->hasSquintPolynomial)
					fprintf(stderr, "\t\033[38;5;208mdPhaseImage %s %3.0f -- %5.3f -- %4i squint on sq=%5.2f\033[0m\n",
							dPhaseImage->file, dParams->nDays, dPhaseImage->par.lambda, dd, dSq);
				else
					fprintf(stderr, "\t\033[38;5;208mdPhaseImage %s %3.0f -- %5.3f -- %4i squint on sq=\033[31m%5.2f\033[38;5;208m\033[0m\n",
							dPhaseImage->file, dParams->nDays, dPhaseImage->par.lambda, dd, dSq);
			}
			if (dPhaseImage->passType == aPhaseImage->passType && sepAscDesc == TRUE)
				continue;
			/*
			  Get region of  possible intersection;   This provides final iMin/iMax.. to loop over
			*/
			getIntersect(dPhaseImage, aPhaseImage, &iMin, &iMax, &jMin, &jMax, outputImage);
			if (iMax == 0 && jMax == 0)
				continue; /* Skip this image if no data */
			/*
			  Init conversion stuff
			*/
			dPhaseImage->memChan = MEM2;
			dCp = setupGeoConversions(dPhaseImage, &dum1, &dum2, &dRe, &dReH, &dThetaC, &dReHfixed, &dThetaCfixedReH);
			dPhaseImage->tolerance = geoTolerance;
			/* fprintf(stderr,"%i %i %i %i", iMin, iMax, jMin, jMax); */
			twokD = (4.0 * PI) / dPhaseImage->par.lambda;
			/*   Compute approximate heading by sampling overlap region  - set iMax,jMax zero if no good solution */
			computeSceneAlpha(outputImage, aPhaseImage, dPhaseImage, aCp, dCp, dem, &iMin, &iMax, &jMin, &jMax);
			if (iMax == 0 && jMax == 0)
				continue; /* no data in range, so skip */
			/* 2026-06-16: also skip when computeSceneAlpha returns inverted bounds (not caught by 0,0 sentinel) */
			if (iMin > iMax || jMin > jMax)
				continue;
			/*  Read in descending image if needed (i.e., nozero intersect).	*/
			setBuffer(dPhaseImage, DImageBuffer);
			getAzimuthBoundsForXYBox(iMin, iMax, jMin, jMax, dPhaseImage, outputImage, &azimuthMin, &azimuthMax);
			/* 2026-06-16: skip read when getAzimuthBoundsForXYBox finds no intersection (azimuthMin=azimuthMax=0) */
			if (azimuthMin == 0.0f && azimuthMax == 0.0f)
				continue;
			getMosaicInputImage(dPhaseImage, azimuthMin, azimuthMax);
			getIonospherePhaseImage(dPhaseImage, DIonBuffer, azimuthMin, azimuthMax);
			/*
			  Loop over output grid and compute velocities
			*/
			fprintf(stderr, "\t---- \033[1mAsc %i / %i\033[0m \033[1mDes %i\033[0m\n", aa, nTotal, dd);
			gettimeofday(&start, NULL);
			count = 0;
			{
				int t;
				for (t = 0; t < nthreads; t++) { localAImgs[t] = *aPhaseImage; localDImgs[t] = *dPhaseImage; }
			}
			double dbgSumPeA = 0, dbgSumPeA2 = 0, dbgSumPeD = 0, dbgSumPeD2 = 0;
			int64_t dbgN = 0;
#pragma omp parallel \
			private(j, x, y, lat, lon, zWGS84, \
			        aZSp, dZSp, arange, drange, aAzimuth, dAzimuth, \
			        aPhase, dPhase, aIonPhase, dIonPhase, aReH, dReH, aRange, dRange, \
			        aTheta, dTheta, aThetaD, dThetaD, aPsi, dPsi, \
			        aPhiZ, dPhiZ, phaseErrorA, phaseErrorD, \
			        aP, dP, aPe, dPe, scaleA, scaleD, \
			        vx, vy, vz, scX, scY, dzdx, dzdy, \
			        dzdtSubmergence, deltaOffCenter, sMask, validData) \
			reduction(+: dbgSumPeA, dbgSumPeA2, dbgSumPeD, dbgSumPeD2, dbgN)
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
					                    row -- see make3DOffsets.c */
					if ((i % 100) == 0)
						fprintf(stderr, "\t--+ %i\n", i);
					/* y - coordinate */
					y = (outputImage->originY + i * outputImage->deltaY) * MTOKM;
					for (j = jMin; j < jMax; j++)
					{
						/*  x-coordinate, then convert x/y stereographic coords to lat/lon	*/
						x = (outputImage->originX + j * outputImage->deltaX) * MTOKM;
						xyToLLProj(x, y, &lat, &lon, &(outputImage->proj));
						zWGS84 = getXYHeight(lat, lon, dem, 0.0, ELLIPSOIDAL);
						/*
						   Process points where elevation is known
						*/
						validData = FALSE; /* Assume didn't work until successful */
						if (zWGS84 > MINELEVATION)
						{
							/*   Convert elevations to spherical reference	*/
							aZSp = sphericalElev(zWGS84, lat, aRe);
							dZSp = sphericalElev(zWGS84, lat, dRe);
							/*
							  Compute range azimuth position
							*/
							llToImageNew(lat, lon, zWGS84, &arange, &aAzimuth, myAImg);
							geometryInfo(aCp, myAImg, aAzimuth, arange, aZSp, aThetaC, &aReH, &aRange, &aTheta, &aThetaD, &aPsi, aZSp);
							llToImageNew(lat, lon, zWGS84, &drange, &dAzimuth, myDImg);
							geometryInfo(dCp, myDImg, dAzimuth, drange, dZSp, dThetaC, &dReH, &dRange, &dTheta, &dThetaD, &dPsi, dZSp);
							/*  Interpolate Phase*/
							interpPhaseImage(myAImg, arange, aAzimuth, &aPhase);
							interpPhaseImage(myDImg, drange, dAzimuth, &dPhase);
							/*  If shelf mask, get mask value */
							sMask = GROUNDED;
							if (shelfMask != NULL)
								sMask = getShelfMask(shelfMask, x, y);
							if (sMask == NOSOLUTION)
							{
								aPhase = -LARGEINT;
								dPhase = -LARGEINT;
							};
							/*
							  If there is valid phase data from both images then compute velocity
							*/
							if (aPhase > -LARGEINT && dPhase > -LARGEINT && zWGS84 > MINELEVATION && sMask != GROUNDINGZONE && (!(sMask == SHELF && outputImage->noTide == TRUE)))
							{
								/*
								  Compute phase due to topography.  Everything is looped through
								  and only the pairs where there is a good angular seperation are used. In general, one will be asc and one will be desc, but they could be flipped.
								  As a consequence, this next step has to look at the flag to see if it should flip the azimuth coordinate when computing the baseline.
								*/
								if (aParams->applyFlatEarth) {
									aPhiZ = computePhiFlatEarthM3d(aAzimuth, aParams, myAImg, aRange, aRe, aReHfixed, aThetaCfixedReH, &phaseErrorA);
									dPhiZ = computePhiFlatEarthM3d(dAzimuth, dParams, myDImg, dRange, dRe, dReHfixed, dThetaCfixedReH, &phaseErrorD);
								} else {
									aPhiZ = computePhiZM3d(&aThetaD, aZSp, aAzimuth, aParams, myAImg, aRange, aRe, aReH, aReHfixed, aThetaC, aThetaCfixedReH, &phaseErrorA);
									dPhiZ = computePhiZM3d(&dThetaD, dZSp, dAzimuth, dParams, myDImg, dRange, dRe, dReH, dReHfixed, dThetaC, dThetaCfixedReH, &phaseErrorD);
								}
								aPhase = aPhase - aPhiZ;
								dPhase = dPhase - dPhiZ;
								/*  Ionosphere corrections, still in radians and before the
								    velocity scaling below. The file holds the ionospheric phase
								    itself, in the same sign convention as the phase image, so it
								    is SUBTRACTED -- unlike the range-offset ionosphere correction
								    in make3DOffsets.c, which is a pre-negated correction that is
								    added. See tiePoints/computeBaseline.c for the derivation. */
								if (myAImg->ionospherePhase != NULL)
								{
									interpIonPhaseImage(myAImg, arange, aAzimuth, &aIonPhase);
									if (aIonPhase > -0.98 * LARGEINT)
									{
										aPhase -= aIonPhase;
									}
								}
								if (myDImg->ionospherePhase != NULL)
								{
									interpIonPhaseImage(myDImg, drange, dAzimuth, &dIonPhase);
									if (dIonPhase > -0.98 * LARGEINT)
									{
										dPhase -= dIonPhase;
									}
								}
								/*  Tide corrections	*/
								if (sMask == SHELF)
								{
									/* update tide correct, and compute phaseImage->tideCorrection */
									interpTideError(&phaseErrorA, myAImg, aParams, x, y, aPsi, twokA);
									interpTideError(&phaseErrorD, myDImg, dParams, x, y, dPsi, twokD);
									aPhase -= -myAImg->tideCorrection * cos(aPsi) * twokA * (double)aParams->nDays / 365.25;
									dPhase -= -myDImg->tideCorrection * cos(dPsi) * twokD * (double)dParams->nDays / 365.25;
								} /* ENd if(smask... */
								/* Submergence corrections */
								if (vCorrect != NULL)
								{
									dzdtSubmergence = interpVCorrect(x, y, vCorrect);
									aPhase -= -dzdtSubmergence * cos(aPsi) * twokA * (double)aParams->nDays / 365.25;
									dPhase -= -dzdtSubmergence * cos(dPsi) * twokD * (double)dParams->nDays / 365.25;
								}
								/*  Update A every 3rd pixel; rowAset guarantees a fresh A on the first
								    valid-data pixel of every row -- see make3DOffsets.c */
								if ((j % 3) == 0 || rowAset == FALSE)
								{
									extern int32_t useSquint;
									computeA(lat, lon, x, y, myAImg, myDImg, A, useSquint, &(outputImage->proj));
									rowAset = TRUE;
								}
								/*
								  Only pursue solution if sufficient difference  in angles for 3d solution
								*/
								if (A[0][0] != -LARGEINT)
								{
									/*  Compute B (note B is really C in the TGARS paper	*/
									computeB(x, y, zWGS84, B, &dzdx, &dzdy, aPsi, dPsi, (xyDEM *)dem);
									if (sMask == SHELF)
									{
										/* Zero slope coupling in the crossing-pair solve on ice shelves --
										   see make3DOffsets.c for the full explanation (matches
										   speckleTrackMosaic.c's 10/13/17 precedent). */
										B[0][0] = 0.0; B[0][1] = 0.0; B[1][0] = 0.0; B[1][1] = 0.0;
									}
									/*  Scale phases for velocity computation (scale for m/yr)	*/
									scaleA = 365.25 / (twokA * aParams->nDays * sin(aPsi));
									scaleD = 365.25 / (twokD * dParams->nDays * sin(dPsi));
									aP = aPhase * scaleA;
									dP = dPhase * scaleD;
									aPe = phaseErrorA * scaleA;
									dPe = phaseErrorD * scaleD;
									dbgSumPeA  += phaseErrorA;
									dbgSumPeA2 += phaseErrorA * phaseErrorA;
									dbgSumPeD  += phaseErrorD;
									dbgSumPeD2 += phaseErrorD * phaseErrorD;
									dbgN++;
									/*  Compute velocity */
									computeVxy(aP, dP, aPe, dPe, A, B, &vx, &vy, &scX, &scY);
									/*  Compute vertical velocity	*/
									vz = vx * dzdx + vy * dzdy;
									/*
									  Update output arrays
									*/
									if (!(scX > -1000. && scX < 1000.))
										error("invalid velocity %f %f %f %f %f %f\n", vx, vy, phaseErrorA, phaseErrorD, aPe, dPe);
									vxTmp[i][j] = vx * scX;
									vyTmp[i][j] = vy * scY;

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
									if (pairOverCount == TRUE)
									{
										nDtmp[i][j] += 1.0f;
										aContrib[i * outputImage->xSize + j] = 1;
									}
									fScale[i][j] = 1.0; /* Value for zero feathering */
#pragma omp atomic write
									aPhaseImage->used = TRUE;
#pragma omp atomic write
									dPhaseImage->used = TRUE;
									validData = TRUE;
								} /* else fprintf(stderr,"LARGEA\n"); */
							}
						}
						if (validData == FALSE)
						{
							vxTmp[i][j] = (float)-LARGEINT;
							fScale[i][j] = 0.0;
						}
					} /* j loop */
				}	  /* i loop */
			} /* End omp parallel */
			if (dbgN > 0) {
				double mA = dbgSumPeA / dbgN, mD = dbgSumPeD / dbgN;
				double sA = sqrt(dbgSumPeA2 / dbgN - mA * mA);
				double sD = sqrt(dbgSumPeD2 / dbgN - mD * mD);
				double radToCmA = aPhaseImage->par.lambda / (4 * PI) * 100.0;
				double radToCmD = dPhaseImage->par.lambda / (4 * PI) * 100.0;
				fprintf(stderr,
					"\tPHASEERR n=%ld  peA: %.4f+/-%.4f rad (%.3f cm)  peD: %.4f+/-%.4f rad (%.3f cm)\n",
					(long)dbgN, mA, sA, mA * radToCmA, mD, sD, mD * radToCmD);
			}
				  /*
					Compute scale array for feathering.
				  */
			gettimeofday(&stop, NULL);
			/* fprintf(stderr,"T  %lf\n", (double)(stop.tv_usec - start.tv_usec)/1e6 + (double)(stop.tv_sec - start.tv_sec));*/
			if (fl > 0 && (iMax > 0 && jMax > 0))
				computeScaleLS((float **)vxTmp, fScale, outputImage->ySize, outputImage->xSize, fl, (float)1.0, (double)(-LARGEINT), iMin, iMax, jMin, jMax);
			/*
			  Now sum current result. Falls through if no intersection (iMax&jMax==0)
			*/
			if (outputImage->timeOverlapFlag == TRUE)
			{
				combWeight = sqrt(aPhaseImage->weight * dPhaseImage->weight);
				fprintf(stderr, "\t\033[1mComb weight = %lf |Ta-Td| %lf\033[0m\n", combWeight, fabs(aPhaseImage->julDay - dPhaseImage->julDay));
			}
			else
				combWeight = 1.0;
			redoNormalization(combWeight, outputImage, iMin, iMax, jMin, jMax, vXimage, vYimage, vZimage, errorX, errorY,
							  scaleX, scaleY, scaleZ, fScale, vxTmp, vyTmp, vzTmp, sxTmp, syTmp, FALSE);
			/* ******************************
			   Use end of goto used to skip inner loop for nophase */
		} /* End desc loop */
		if (pairOverCount == TRUE)
		{
			int i1, j1;
			for (i1 = 0; i1 < outputImage->ySize; i1++)
				for (j1 = 0; j1 < outputImage->xSize; j1++)
					if (aContrib[i1 * outputImage->xSize + j1]) {
						nAtmp[i1][j1] += 1.0f;
						aContrib[i1 * outputImage->xSize + j1] = 0;
					}
		}
	}	  /* End asc loop */
	free(localAImgs);
	free(localDImgs);
	if (aContrib != NULL)
		free(aContrib);
	/**************************END OF MAIN LOOP ******************************/
	fprintf(stderr, "Out of main loop\n");
	{
		extern double totalPhaseIOTime;
		gettimeofday(&funcEnd, NULL);
		fprintf(stderr, "Total phase I/O time (whole run): %.3f s\n", totalPhaseIOTime);
		fprintf(stderr, "Total phase processing time (whole run): %.3f s\n",
		        (funcEnd.tv_sec - funcStart.tv_sec) + (funcEnd.tv_usec - funcStart.tv_usec) * 1e-6);
	}
	/* Correct for pair over-counting.  Must run before endScale(), on this round's own
	   error contribution only -- see inflatePairOverCount() in common/scalingFunctions.c. */
	if (pairOverCount == TRUE)
	{
		extern double rhoPhase; /* defined in common/getRegion.c */
		int i1;
		inflatePairOverCount(outputImage, errorX, errorY, errorX0, errorY0, nAtmp, nDtmp, rhoPhase);
		for (i1 = 0; i1 < outputImage->ySize; i1++)
		{
			free(nAtmp[i1]);
			free(nDtmp[i1]);
			free(errorX0[i1]);
			free(errorY0[i1]);
		}
		free(nAtmp);
		free(nDtmp);
		free(errorX0);
		free(errorY0);
	}
	/*
	  Adjust scale
	*/
	endScale(outputImage, vXimage, vYimage, vZimage, errorX, errorY, scaleX, scaleY, scaleZ, FALSE);
	fprintf(outputImage->fpLog, ";\n; Returning from make3DOffs(.c)\n");
	fflush(outputImage->fpLog);
}

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
	for (i = 1; i <= 6; i++) {
		tmpV[i] = 0;
		for (j = 1; j <= 6; j++)
			tmpV[i] += vhParam->C[i][j] * v[j];
	}
	sig2Base = 0.0;
	for (j = 1; j <= 6; j++)
		sig2Base += tmpV[j] * v[j];

	delta = -bn * sinThetaD - bp * cosThetaD + bSq * 0.5 / Range;
	*phaseError = sqrt(sig2Base + min(6*PI, vhParam->sigma) * min(6*PI, vhParam->sigma));
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
	*phaseError = sqrt(sig2Base + min(6*PI, vhParam->sigma) * min(6*PI, vhParam->sigma)); /* note this is returning phase error as sigma */
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
