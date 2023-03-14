#include "math.h"
#include "common.h"
static  void computeXYslope(double lat,double lon, double *dzdx, double *dzdy,  xyDEM xydem,double *slopeMag);

static double xyScaleApprox(double lat, double stdLat)
{
/*
  This is a much simpler formula than the xyscale in the rest of the code, but it agrees to at least 3 or 4 digits, which is close enough for scaling derivatives
  Assumes lat is same hemisphere as slat
*/
	return (1.0 + sin(fabs(lat) * DTOR))/(1.0 + sin(stdLat*DTOR));
}

void xyGetZandSlope(double lat, double  lon,  double x, double y, double *zSp, double *zWGS84,double *da, double *dr,
		    conversionDataStructure *cP,vhParams *vhParam, inputImageStructure *image)
{
	double gRange,azimuth;
	double range,z1,azimuth1,gRange1;
	double dzdx,dzdy,slopeMag;
	double psScale;
	double hAngle, xyAngle;
	double lati,loni,xi,yi;
	double zOrig;
	int i,j;   
	/*
	  default values in case of early return
	*/
	*da = 0.0;	*dr = 0.0;
	/* 
	   Conversion from output grid to perhaps alternate xy projection  added 9/11/96
	*/
	lltoxy1(lat, lon, &xi, &yi, vhParam->xydem.rot, vhParam->xydem.stdLat);
	psScale = xyScaleApprox(lat, vhParam->xydem.stdLat);
	/* 
	   Get height
	*/
	*zWGS84 = interpXYDEM(xi, yi, vhParam->xydem); 
	zOrig = *zWGS84;
	if( *zWGS84 <= (MINELEVATION+1) || *zWGS84 > 10000.0) return;
	/*
	  Ellipsoidal to spherical earth correction
	*/
	*zSp = *zWGS84 + (earthRadius(lat*DTOR,image->par.ReMinor, image->par.ReMajor) * KMTOM  - cP->Re);
	/*
	  Compute flow direction in xy coords from dem and angle of x from north 
	*/
	computeXYslope(lat, lon, &dzdx, &dzdy, vhParam->xydem, &slopeMag);
	/*
		Correct for fact dx, dy are polar stereo coords not real distances
	*/
	dzdx = dzdx / psScale;
	dzdy = dzdy / psScale;
	slopeMag = slopeMag / psScale;
	/*fprintf(stderr,"%f\n",psScale);*/
	/*
	  Avoid really large slopes ???
	*/ 
	dzdx = limitSlope(dzdx, 0.12);
	dzdy = limitSlope(dzdy, 0.12);
	/* if(slopeMag > 0.1) {dzdx= dzdx/slopeMag * 0.1; dzdy =dzdy/slopeMag * 0.1;} */
	/*
	  Take care of slopes in low areas near ice shelf fronts
	*/
	if(slopeMag > 0.05 && zOrig < 75) {dzdx=0; dzdy=0;}
	if(slopeMag > 0.03 && zOrig < 50) {dzdx=0; dzdy=0;}
	/* 
	   Compute angle between range direction and north 
	*/        
	hAngle = computeHeading(lat, lon, 0, image, cP);
	if(hAngle > (2.0001*PI)) { *zSp = 1000000.; *zWGS84=1000000;  return; }
	/* 
	   Rotate dzdx,dzdy to dr,da 
	*/
	computeXYangle(lat, lon, &xyAngle, vhParam->xydem);
	rotateFlowDirectionToRA(dzdx, dzdy, da, dr, xyAngle, hAngle);
	return;
}

/*
  Compute flow direction and angle relative to north or south. 
*/
static  void computeXYslope(double lat,double lon, double *dzdx, double *dzdy,  xyDEM xydem,double *slopeMag)
{   
	extern int HemiSphere;
	double x,y,z0x, z0y, z1x,z1y;
	double dx,dy;

	dx = xydem.deltaX ;
	dy = xydem.deltaY;
	lltoxy1(lat, lon, &x, &y, xydem.rot, xydem.stdLat);
	z0x = interpXYDEM(x - dx/2.,  y, xydem);
	z1x = interpXYDEM(x + dx/2. ,y, xydem);
	z0y = interpXYDEM(x, y - dy/2., xydem);		
	z1y = interpXYDEM(x, y + dy/2, xydem);	
	if( z0x <= MINELEVATION || z1x <= MINELEVATION ||
	    z0y <= MINELEVATION || z1y <= MINELEVATION ||
	    z0x > 10000.0 || z0y > 10000.0 ||
	    z1x > 10000.0 || z1y > 10000.0 ) {
		*dzdx=0; *dzdy = 0.0; *slopeMag = 0.0; return;
	}
	/*
	  Compute gradient
	*/
	*dzdx = (z1x - z0x) / (KMTOM * dx);
	*dzdy = (z1y - z0y) / (KMTOM * dy);
	/*
	  compute slope mag
	*/
	*slopeMag = sqrt(*dzdx * (*dzdx) + *dzdy * (*dzdy));
	return;
}
