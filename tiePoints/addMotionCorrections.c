#include <math.h>
#include "mosaicSource/common/common.h"
#include "tiePoints.h"

/*
  Compute motion correction.
*/
void addMotionCorrections(inputImageStructure inputImage,  tiePointsStructure *tiePoints)
{   
	extern int HemiSphere;
	int i,j;
	conversionDataStructure *cP;
	double x,y,z, zWGS;
	double lat,lon;
	double azOrigin;
	double num,den,rho1,rho2 ;
	double range1,range2,azimuth1,azimuth2;
	double vyra;
	double lambda1;
	double da,dgr;
	double rot;
	double ReH, Re, RNear;
	double theta,psi,r1;
	double gRange1,gRange2;
	double hAngle,xyAngle,rotAngle;
	double x1,y1;
	double tSign;
	double deltaT;
	double twok;
	double dlat;
	double hAngle1;
	double phiDisplacement;
	float **fimage;
	/*
	  Setup stuff for Shusun Li's conversion routines
	*/
	if(inputImage.isInit != TRUE)   initllToImageNew( &inputImage);	
	lambda1=inputImage.par.lambda;
	cP = &(inputImage.cpAll);
	Re = cP->Re;
	RNear = cP->RNear;
	fprintf(stderr,"---------------------RNear %f\n",RNear);
	twok = (4.0 * PI / lambda1);
	deltaT = tiePoints->nDays/365.25;
	/*
	  Loop over tie points
	*/
	for(i=0; i < tiePoints->npts; i++) {
		lat = tiePoints->lat[i];
		lon = tiePoints->lon[i];
		z = tiePoints->z[i];
		/*
		  Convert lat/lon to image coordinates
		*/                     
		z = sphericalToWGSElev(z, lat, Re); /* 0.0;*/
		llToImageNew(lat, lon, z ,&range1, &azimuth1, &inputImage);
		hAngle = computeHeading(lat, lon, z, &inputImage, &(inputImage.cpAll));
		/*
		  Compute angle for ps coordinats from north
		*/
		computeXYangleNoDem(lat, lon, &xyAngle, tiePoints->stdLat);
		/*
		  Compute angle for CW rotation of ps to radar
		*/
		rotAngle = hAngle - xyAngle;
		/*
		  Compute range displacement vyra by ccw rot from ps to ra
		*/         
		if(tiePoints->vrFlag==FALSE) {
			vyra = tiePoints->vx[i] * cos(rotAngle) - tiePoints->vy[i] * sin(rotAngle);
		} else vyra=tiePoints->vx[i];

		z = tiePoints->z[i];
		zWGS = sphericalToWGSElev(z, lat, Re);
		llToImageNew(lat, lon, zWGS, &range1, &azimuth1, &inputImage);
		ReH = getReH(cP, &inputImage, azimuth1);
		r1 = RNear + inputImage.nRangeLooks * RangePixelSize * range1;
		theta = thetaRReZReH(r1, (Re+z), ReH);
		psi = psiRReZReH(r1, (Re + z), ReH);
		/*
		  Compute phase correction for motion
		*/
		tSign = 1.0;
		if(tiePoints->timeReverseFlag == TRUE) tSign = -1.0;
		phiDisplacement = tSign * twok * deltaT * (vyra * sin(psi) - tiePoints->vz[i] * cos(psi));
		tiePoints->phase[i] -= phiDisplacement;
		tiePoints->vyra[i] = vyra;
	}

	return;
}
