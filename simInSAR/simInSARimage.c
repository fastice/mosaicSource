#include "stdio.h"
#include "string.h"
#include <sys/types.h>
#include <math.h>
#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"
#include <omp.h>
/*
  Function to simulate InSAR image including both terrain and motion effects.

  01/04/00 Modified to handle block updates of Re, ReH. This is kluged a
  little but since blocks are update by range and rg found uniquely for
  that range it should be no problem.
*/

static void withHeight(inputImageStructure *inputImage, double *lat, double *lon, double azimuth, double rsl, double zWGS)
{
	int32_t n;
	stateV *sv;
	conversionDataStructure *cp;
	double xs, ys, zs, vsx, vsy, vsz;
	double myTime;
	cp = &(inputImage->cpAll);
	myTime = (azimuth * inputImage->nAzimuthLooks) / cp->prf + cp->sTime;
	sv = &(inputImage->sv);
	n = (int32_t)((myTime - sv->times[1]) / (sv->deltaT) + .5);
	n = min(max(0, n - 2), sv->nState - NUSESTATE);
	polintVec(&(sv->times[n]), &(sv->x[n]), &(sv->y[n]), &(sv->z[n]), &(sv->vx[n]), &(sv->vy[n]), &(sv->vz[n]), myTime, &xs, &ys, &zs, &vsx, &vsy, &vsz);

	smlocateZD(xs * MTOKM, ys * MTOKM, zs * MTOKM, vsx * MTOKM, vsy * MTOKM, vsz * MTOKM, rsl * MTOKM, lat, lon, (double)(inputImage->lookDir), zWGS * MTOKM);
}

/*
  Compute flow direction and angle relative to north or south.
*/
static void computeXYslope(double x, double y, double *dzdx, double *dzdy, xyDEM xydem)
{
	double z, z1x, z1y;
	double dx, dy;

	dx = xydem.deltaX;
	dy = xydem.deltaY;
	z = interpXYDEM(x, y, xydem);
	z1x = interpXYDEM(x + dx, y, xydem);
	z1y = interpXYDEM(x, y + dy, xydem);
	if (z <= MINELEVATION || z1x <= MINELEVATION || z1y <= MINELEVATION ||
		z > 10000.0 || z1x > 10000.0 || z1y > 10000.0)
	{
		*dzdx = 0;
		*dzdy = 0.0;
		return;
	}
	/*
	  Compute gradient
	*/
	*dzdx = (z1x - z) / (KMTOM * dx);
	*dzdy = (z1y - z) / (KMTOM * dy);
	return;
}

/************************STATIC ROUTINE DECLARATIONS**************************/

unsigned char getShelfMask(ShelfMask *imageMask, double x, double y);
/************************END STATIC ROUTINES DECLARATIONS*********************/

void interpXYVel(double x, double y, xyVEL *xyVel, double *vx, double *vy)
{
	double xi, yi;
	int32_t i, j, iSize, jSize;
	double t, u, p1, p2, p3, p4;
	double z;

	xi = (x - xyVel->x0) / xyVel->deltaX;
	yi = (y - xyVel->y0) / xyVel->deltaY;
	*vx = bilinearInterp(xyVel->vx, xi, yi, xyVel->xSize, xyVel->ySize, -LARGEINT, -LARGEINT);
	*vy = bilinearInterp(xyVel->vy, xi, yi, xyVel->xSize, xyVel->ySize, -LARGEINT, -LARGEINT);
}

void simInSARDEMBounds(sceneStructure *scene, double rot, double stdLat,
                       double *xMin, double *xMax, double *yMin, double *yMax)
{
	/* Compute x/y bounding box of the simulation area by sampling 9 corner/mid
	   points, so callers can crop the DEM and velocity map before reading them. */
	inputImageStructure img = scene->I;   /* local copy; initllToImageNew writes cpAll */
	double aOff, rOff, Re, ReH, iFloat, range, rhoSp, rg, lat, lon, x, y;
	int ii, jj;
	int iSamp[3], jSamp[3];

	initllToImageNew(&img);
	aOff = 0.5 * (img.nAzimuthLooks - 1) / img.nAzimuthLooks;
	rOff = 0.5 * (img.nRangeLooks - 1) * img.par.slpR;
	Re   = img.cpAll.Re;

	*xMin =  1e15;  *xMax = -1e15;
	*yMin =  1e15;  *yMax = -1e15;

	iSamp[0] = 0;  iSamp[1] = scene->aSize / 2;  iSamp[2] = scene->aSize - 1;
	jSamp[0] = 0;  jSamp[1] = scene->rSize / 2;  jSamp[2] = scene->rSize - 1;

	for (ii = 0; ii < 3; ii++) {
		iFloat = iSamp[ii] * scene->dA + scene->aO - aOff - img.cpAll.azOff;
		ReH    = getReH(&img.cpAll, &img, iFloat);
		for (jj = 0; jj < 3; jj++) {
			range = img.cpAll.RNear + (scene->rO + jSamp[jj] * scene->dR)
			        * img.rangePixelSize - rOff;
			rhoSp = rhoRReZReH(range, Re, ReH);
			rg    = Re * rhoSp - 100.;
			groundRangeToLLNew(rg, iFloat, &lat, &lon, &img, FALSE);
			lltoxy1(lat, lon, &x, &y, rot, stdLat);
			if (x < *xMin) *xMin = x;
			if (x > *xMax) *xMax = x;
			if (y < *yMin) *yMin = y;
			if (y > *yMax) *yMax = y;
		}
	}
	/* 50 km margin covers terrain relief and numerical rounding */
	*xMin -= 50.;  *xMax += 50.;
	*yMin -= 50.;  *yMax += 50.;
	fprintf(stderr, "Scene DEM crop bounds: x=[%.1f,%.1f] y=[%.1f,%.1f] km\n",
	        *xMin, *xMax, *yMin, *yMax);
}

void simInSARimage(sceneStructure *scene, void *dem, xyVEL *xyVel)
{
	inputImageStructure *inputImage;
	xyDEM *xyDem;
	double R, theta, thetaD, delta, thetaC, rhoSp, rhoWGS, RCenter, RNear, RWGS, R0, ReWGS, ReH, Re;
	double lat, lon;
	double baselineSquared;
	double toRho, r, rg, drgStep, range;
	double deltaToPhase;
	double iFloat;
	double hSp, hWGS;
	double x1, y1;
	double rOff, aOff;
	double vx, vy, vr, psi, dzdx, dzdy, dzdr, dzdtSubmergence;
	double speed, Xv;
	double hAngle, xyAngle;
	int64_t nIterTotal, nIter;
	int32_t recycle;
	int32_t iLoop, jLoop, i, iIndex;
	fprintf(stderr, "SIM INSAR IMAGE\n");
	/*	  Dem type	*/
	xyDem = (xyDEM *)dem;
	/*
	  Init coord conversions
	*/
	inputImage = &(scene->I);
	initllToImageNew(inputImage);
	deltaToPhase = 4.0 * PI / scene->I.par.lambda;
	fprintf(stderr, "Delta to phase %f lambda %f\n", deltaToPhase, scene->I.par.lambda);
	fprintf(stderr, "xyVel xSize=%d vx=%p vy=%p\n", xyVel->xSize, (void*)xyVel->vx, (void*)xyVel->vy);
	Re = inputImage->cpAll.Re;
	RCenter = inputImage->cpAll.RCenter;
	ReH = inputImage->cpAll.Re + inputImage->par.H;
	thetaC = thetaRReZReH(RCenter, (Re + 0), ReH);
	/* Constants for phase computation */
	fprintf(stderr, "Re,thetaC,RNear,RFar, %f %f   %f %f %f %i \n ",
			Re, thetaC, inputImage->cpAll.RNear, inputImage->cpAll.RFar, inputImage->rangePixelSize, inputImage->nRangeLooks);
	fprintf(stderr, "R/A size %i, %i\n", scene->rSize, scene->aSize);
	/* Corrections from ml pixels to sl pixels */
	rOff = 0.5 * (inputImage->nRangeLooks - 1) * inputImage->par.slpR;
	/* First part in slp so need to divide nlooks to get in mlp pixels */
	aOff = 0.5 * (inputImage->nAzimuthLooks - 1) / inputImage->nAzimuthLooks;
	RNear = inputImage->cpAll.RNear;
	nIterTotal = 0;

	/* Per-thread copies of inputImageStructure for thread-safe coordinate conversion.
	   groundRangeToLLNew and withHeight cache state in cpAll; threads must not share it. */
	int nthreads = omp_get_max_threads();
	inputImageStructure *localImgs =
		(inputImageStructure *)malloc((size_t)nthreads * sizeof(inputImageStructure));
	if (localImgs == NULL) error("simInSARimage: malloc localImgs failed\n");
	{
		int t;
		for (t = 0; t < nthreads; t++) localImgs[t] = *inputImage;
	}

	/*
	  Loop over output grid.
	*/
#pragma omp parallel \
	private(iFloat, i, iIndex, ReH, rhoSp, rg, range, recycle, \
	        nIter, jLoop, R, R0, ReWGS, rhoWGS, hWGS, hSp, lat, lon, \
	        theta, thetaD, delta, baselineSquared, drgStep, \
	        x1, y1, vx, vy, vr, psi, dzdx, dzdy, dzdr, dzdtSubmergence, hAngle, xyAngle) \
	reduction(+: nIterTotal)
	{
		int myThread = omp_get_thread_num();
		inputImageStructure *myImg = &localImgs[myThread];
		conversionDataStructure *myCP = &(myImg->cpAll);
		double localBn, localBp;

#pragma omp for schedule(dynamic, 8)
		for (iLoop = 0; iLoop < scene->aSize; iLoop++)
		{
			iFloat = iLoop * scene->dA + scene->aO - aOff - myCP->azOff;
			i = (int32_t)(iFloat + 0.5);
			iIndex = iLoop;
			ReH = getReH(myCP, myImg, iFloat);
			rhoSp = rhoRReZReH(RNear, (Re + 0), ReH);
			rg = Re * rhoSp - 100;															 /* Initial ground range on spherical earth */
			range = inputImage->cpAll.RNear + scene->rO * inputImage->rangePixelSize - rOff; /* Starting range */
			if (scene->bnArray != NULL)
			{
				localBn = scene->bnArray[iIndex];
				localBp = scene->bpArray[iIndex];
			}
			else
			{
				localBn = scene->bnStart + (double)iIndex * scene->bnStep;
				localBp = scene->bpStart + (double)iIndex * scene->bpStep;
			}
			if (i == 0)
				fprintf(stderr, "%f %f \n", localBn, localBp);
			if (myThread == 0 && i % 100 == 1)
				fprintf(stderr, "%i\n", i);
			/* Pre-compute heading at 3 range positions for quadratic interpolation in jLoop.
			   Runs the full terrain-convergence loop for each anchor pixel so the heading is
			   computed at terrain-corrected lat/lon (same as the per-pixel original path). */
			double hAngle_row[3] = {0., 0., 0.};
			if (xyVel->xSize > 0 && i >= 0 && i < myCP->azSize)
			{
				int jRef[3] = {0, scene->rSize / 2, scene->rSize - 1};
				int k;
				for (k = 0; k < 3; k++)
				{
					int jj = jRef[k];
					double range_k = myCP->RNear + (scene->rO + jj * scene->dR) * inputImage->rangePixelSize - rOff;
					double rhoSp_k = rhoRReZReH(range_k, Re + 0., ReH);
					double rg_k = Re * rhoSp_k - 100;
					double lat_k, lon_k, hWGS_k, ReWGS_k, rhoWGS_k, R_k, R0_k, drgStep_k;
					int32_t recycle_k = FALSE;
					int iter_k = 0;
					while (1) {
						R0_k = groundRangeToLLNew(rg_k, iFloat, &lat_k, &lon_k, myImg, recycle_k);
						recycle_k = TRUE;
						ReWGS_k = earthRadius(lat_k * DTOR, EMINOR, EMAJOR) * KMTOM;
						rhoWGS_k = rhoRReZReH(R0_k, (ReWGS_k + 0), ReH);
						hWGS_k = getXYHeight(lat_k, lon_k, xyDem, Re, ELLIPSOIDAL);
						R_k = slantRange(rhoWGS_k, ReWGS_k + hWGS_k, ReH);
						drgStep_k = (range_k - R_k) * 0.5;
						if (R_k > (range_k - 0.1)) break;
						rg_k += drgStep_k;
						if (iter_k > 8000) break;
						iter_k++;
					}
					withHeight(myImg, &lat_k, &lon_k, iFloat, range_k, hWGS_k);
					hAngle_row[k] = computeHeading(lat_k, lon_k, 0., myImg, myCP);
					if (hAngle_row[k] > 9000.)
						hAngle_row[k] = k > 0 ? hAngle_row[k - 1] : 0.;
				}
			}
			/* Loop over range */
			for (jLoop = 0; jLoop < scene->rSize; jLoop++)
			{
				/* Initialize iteration */
				delta = 0.0;
				nIter = 0;
				recycle = FALSE;
				if (i >= 0 && i < myCP->azSize)
				{
					while (1)
					{
						/* For rg on spherical earth, compute lat/lon for the corresponding ellipsoidal earth
						   and return the corresponding slant range */
						R0 = groundRangeToLLNew(rg, iFloat, &lat, &lon, myImg, recycle);
						recycle = TRUE; /* flags to reuse previously computed values that are still valid */
						/* Compute WGS radius at lat/lon */
						ReWGS = earthRadius(lat * DTOR, EMINOR, EMAJOR) * KMTOM;
						/* Compute rho for the ellipsoidal radius */
						rhoWGS = rhoRReZReH(R0, (ReWGS + 0), ReH);
						/* Get height for lat/lon referenced to the ellipsoid */
						hWGS = getXYHeight(lat, lon, xyDem, Re, ELLIPSOIDAL);
						/* Compute height corrected R */
						R = slantRange(rhoWGS, ReWGS + hWGS, ReH);
						/* Not the most efficient iteration, but good enough */
						drgStep = (range - R) * 0.5;
						if (R > (range - 0.1))
							break;
						/* Speed convergence if way out */
						rg += drgStep; /* Increment ground range */
						/* Avoid infinite loop */
						if (nIter > 8000)
						{
							fprintf(stderr, "%f %f %f\n", Re, rhoSp, ReH);
							break;
						}
						nIter++;
					}
					nIterTotal += nIter;
					withHeight(myImg, &lat, &lon, iFloat, range, hWGS);

					if (scene->maskFlag == TRUE && scene->toLLFlag == FALSE)
					{
						lltoxy1(lat, lon, &x1, &y1, xyDem->rot, xyDem->stdLat);
						if (scene->velThresh > 0) {
							interpXYVel(x1, y1, xyVel, &vx, &vy);
							if (vx > -1.99e9 && vy > -1.99e9)
								scene->image[iIndex][jLoop] = (hypotf((float)vx, (float)vy) > scene->velThresh) ? 0.0f : 1.0f;
							else
								scene->image[iIndex][jLoop] = 0.0f;
						} else {
							scene->image[iIndex][jLoop] = getShelfMask(scene->imageMask, x1, y1);
						}
					}
					else if (scene->heightFlag == TRUE) //&& scene->toLLFlag == FALSE
					{	// Save height only
						scene->image[iIndex][jLoop] = (float)(hWGS);
					}
					else if (scene->toLLFlag == TRUE || scene->saveLLFlag == TRUE)
					{
						// Save lat/lon
						scene->latImage[iIndex][jLoop] = lat;
						scene->lonImage[iIndex][jLoop] = lon;
						// Also save mask
						if (scene->maskFlag == TRUE)
						{
							lltoxy1(lat, lon, &x1, &y1, xyDem->rot, xyDem->stdLat);
							scene->image[iIndex][jLoop] = getShelfMask(scene->imageMask, x1, y1);
						}
					}
					else
					{ 	// Phase case
						/* One lltoxy1 for both the hSp DEM lookup and the velocity lookups below. */
						lltoxy1(lat, lon, &x1, &y1, xyDem->rot, xyDem->stdLat);
						/* Use h on single spherical reference for phase
						  even though locally spherical value was used for location */
						hSp = getXYHeightXY(x1, y1, lat, xyDem, Re, SPHERICAL);
						/*
						  Compute look angle for height h and range R.
						  Note I had to make a decision here whether to use ReH or Re + H
						  to compute phases. I went with the less precise Re+H to maintain
						  internal consistency through out the code. This will should have
						  the effect of introducing and erroneous phase ramp, but since the
						  code is consistent, it should be cancelled via the effective baseline.
						  If I ever choose to update the code to carry Re+h throughout,
						  then this line should replace the theta calculation.
						  theta = acos( ( range*range + ReH*ReH - (Re+h)*(Re+h) )/
						  (2.0*ReH*range) );
						*/
						theta = thetaRReZReH(range, (Re + hSp), ReH);
						thetaD = theta - thetaC;
						/*
						  Add component of delta from topography in meters.
						*/
						baselineSquared = pow(localBn, 2.0) + pow(localBp, 2.0);
						delta += (sqrt(pow(range, 2.0) - 2.0 * range * (localBn * sin(thetaD) + localBp * cos(thetaD)) +
									   baselineSquared) -
								  range);
						/*
						   Flatten phase image if flatFlag set. Include possible baseline estimation error.
						*/
						if (scene->flatFlag == TRUE)
						{
							theta = thetaRReZReH(range, (Re + 0.), ReH);
							thetaD = theta - thetaC;
							delta -= -localBn * sin(thetaD) - localBp * cos(thetaD) + baselineSquared * 0.5 / range;
						}
						/*
						  Get displacement
						*/
						vx = 0.0;
						vy = 0.0;
						vr = 0.0;
						if (xyVel->xSize > 0 || scene->verticalCorrection != NULL)
							psi = psiRReZReH(range, (Re + hSp), ReH);
						if (xyVel->xSize > 0)
						{
							/* x1, y1 already computed above — reuse for velocity lookup */
							//fprintf(stderr, "vel lookup: lat=%f lon=%f x1=%f y1=%f\n", lat, lon, x1, y1);
							interpXYVel(x1, y1, xyVel, &vx, &vy);
							/* Compute rotation angles */
							computeXYangleXY(x1, y1, &xyAngle);
							/* Quadratic interpolation of pre-computed per-row heading */
							{
								double t = scene->rSize > 1 ? (double)jLoop / (scene->rSize - 1) : 0.5;
								hAngle = hAngle_row[0] * 2.0 * (t - 0.5) * (t - 1.0)
								       + hAngle_row[1] * (-4.0 * t * (t - 1.0))
								       + hAngle_row[2] * 2.0 * t * (t - 0.5);
							}
							/* Rotate to ground range */
							double vr_horiz = vx * cos(hAngle - xyAngle) - vy * sin(hAngle - xyAngle);
							vr = vr_horiz;
							/* Get slopes */
							computeXYslope(x1, y1, &dzdx, &dzdy, *xyDem);
							/* Rotate to ground range direction */
							dzdr = dzdx * cos(hAngle - xyAngle) - dzdy * sin(hAngle - xyAngle);

							vr = vr * sin(psi) - vr * dzdr * cos(psi);
							//if (vx > -1.99e9 && vy > -1.99e9)
							//	fprintf(stderr, "DIAG: hAngle=%.4f xyAngle=%.4f rotAngle=%.4f psi=%.4f(deg=%.1f) vr_horiz=%.4f dzdr=%.6f vr_slant=%.4f vx=%.2f vy=%.2f\n",
							//		hAngle, xyAngle, hAngle-xyAngle, psi, psi*180./PI, vr_horiz, dzdr, vr, vx, vy);
							/*
							  Speed-adaptive smoothing tolerance (-minTol/-percentSpeed/-maxTol):
							  reuses sin(psi)/deltaToPhase/dT, the same factors already used above
							  to scale the real velocity into a phase contribution, applied instead
							  to the clamped tolerance value. Approximate (uses the dominant sin(psi)
							  projection, ignores the smaller slope term) by design.
							*/
							if (scene->smoothRadiusFlag == TRUE && vx > -1.99e9 && vy > -1.99e9)
							{
								speed = hypot(vx, vy);
								Xv = fmin(fmax(scene->minTol, scene->percentSpeed / 100. * speed), scene->maxTol);
								scene->toleranceImage[iIndex][jLoop] =
									(float)(Xv * (scene->dT / 365.) * deltaToPhase * sin(psi));
							}
						}
						/*
						  Convert range diff to phase diff.
						*/
						if(vx > -1.99e9 && vy > -1.99e9)
						{
							//fprintf(stderr, "%f\n", scene->image[iIndex][jLoop] );
							scene->image[iIndex][jLoop] = (float)(delta * deltaToPhase);
							//fprintf(stderr, "%f %f %f %f %f %f\n", delta, vr, scene->dT, deltaToPhase, vx, vy);
							scene->image[iIndex][jLoop] += vr * scene->dT / 365. * deltaToPhase;
							/*
							  Vertical (submergence/emergence) contribution to phase. mosaic3d's
							  inverse operation (make3DMosaic.c) ADDS dzdtSubmergence*cos(psi)*twok*nDays/365.25
							  to the observed phase to remove the vertical-motion contribution before
							  solving for horizontal velocity; here we SUBTRACT the same term to add
							  that contribution into the simulated phase.
							*/
							if (scene->verticalCorrection != NULL)
							{
								dzdtSubmergence = interpVCorrect(x1, y1, scene->verticalCorrection);
								scene->image[iIndex][jLoop] -= dzdtSubmergence * cos(psi) * scene->dT / 365. * deltaToPhase;
							}
						} else
						{
							scene->image[iIndex][jLoop] = -LARGEINT;
						}

					} /* End else */
					range += inputImage->rangePixelSize * scene->dR;
				}
				else
				{
					/* No data cases because outside az range */
					if (scene->maskFlag == TRUE && scene->toLLFlag == FALSE)
					{
						scene->image[iIndex][jLoop] = 0;
					}
					else if (scene->heightFlag == TRUE && scene->toLLFlag == FALSE)
					{
						scene->image[iIndex][jLoop] = (float)-9999.0;
					}
					else if (scene->toLLFlag == TRUE)
					{
						scene->latImage[iIndex][jLoop] = -LARGEINT;
						scene->lonImage[iIndex][jLoop] = -LARGEINT;
						if (scene->maskFlag == TRUE)
						{
							scene->image[iIndex][jLoop] = 0.0;
						}
					}
					else
					{
						scene->image[iIndex][jLoop] = -2.0e9;
					}
				}
			} /* end for jLoop ... */
		}	  /* end for iLoop .... */
	} /* end omp parallel */

	free(localImgs);
	fprintf(stderr, "nIt Avg %f\n", (double)nIterTotal / (scene->rSize * scene->aSize));
	fprintf(stderr, "bn bp %f %f \n", scene->bn, scene->bp);
	return;
}
