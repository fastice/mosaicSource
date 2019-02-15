#include <math.h>
#include "mosaicSource/common/common.h"
#include "rparams.h"

/*
  Compute motion correction.
*/
void addVelCorrections(inputImageStructure *inputImage, tiePointsStructure *tiePoints)
{   
	extern int HemiSphere;
	int i,j;
	conversionDataStructure *cP;
	double x,y;
	double zWGS,zSp;    
	double lat,lon;
	double range1,azimuth1;
	double vyra;
	double rot;
	double H, Re, RNear,ReH;
	double r1, theta,psi;
	double hAngle,xyAngle,rotAngle;
	double tSign;
	double deltaT;
	double displacement;
	float **fimage;
	/* Setup conversions */
	cP = &(inputImage->cpAll);
	H = inputImage->par.H;
	Re = cP->Re;
	RNear = inputImage->par.rn;
	fprintf(stderr,"---------------------RNear,H,Re %f %f %f\n",RNear,H,Re);
	deltaT = tiePoints->nDays/365.25;
	/*
	  Loop over tie points
	*/
	for(i=0; i < tiePoints->npts; i++) {
		if(tiePoints->phase[i] > (-LARGEINT/10)) {
			lat = tiePoints->lat[i];
			lon = tiePoints->lon[i];
			 /*
			  Compute angle for CW rotation of ps to radar
			*/
			hAngle=computeHeading( lat, lon,0,inputImage, cP);
			computeXYangleNoDem(lat,lon,&xyAngle,tiePoints->stdLat);
			rotAngle = hAngle - xyAngle;
			/*
			  Compute range displacement vyra by ccw rot from ps to ra
			*/         
			if(tiePoints->vrFlag==FALSE) {
				vyra = tiePoints->vx[i] * cos(rotAngle) - tiePoints->vy[i] * sin(rotAngle);
			} else vyra=tiePoints->vx[i];
			zSp=tiePoints->z[i];
			zWGS = tiePoints->z[i] - Re + earthRadius(tiePoints->lat[i]*DTOR ,EMINOR,EMAJOR) * 1000.0 ;
			/*fprintf(stderr,"%f %f %f %f\n",Re, earthRadius(tiePoints->lat[i]*DTOR ,EMINOR,EMAJOR) * 1000.0 ,zSp,zWGS);*/
			llToImageNew(lat,lon,zWGS,&range1,&azimuth1,inputImage);
			if(cP->ReH == NULL) error("geometryInfo : ReH[] not defined");
			if(azimuth1 < cP->azSize && azimuth1 >= 0) ReH = cP->ReH[(int)azimuth1];
			r1 = RNear + inputImage->nRangeLooks * RangePixelSize*range1;
			psi = psiRReZReH(r1,(Re+zSp),ReH);
			/*
			  Compute phase correction for motion
			*/
			tSign =1.0;
			if(tiePoints->timeReverseFlag == TRUE) tSign= -1.0;
			displacement = tSign*deltaT*(  vyra*sin(psi)  -tiePoints->vz[i]*cos(psi)  );
			tiePoints->phase[i] -=  displacement;
			tiePoints->vyra[i] =  vyra;
			/*
			if(fabs(vyra) > 10)
			fprintf(stderr,"vyra %f %f %f\n",vyra,tiePoints->vx[i],tiePoints->vy[i]);*/
		}
	}
	return;
}
