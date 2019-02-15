#include "stdio.h"
#include"string.h"
#include <math.h>
#include "mosaicSource/common/common.h"
/*
  Compute height for lat/lon from DEM.
  Modified 2/7/96 to add mode for converting from ellipsoidal to spherical
  heights.

  lat, lon inputs should be in degrees.
  Re in KM
*/    

double getXYHeight(double lat, double lon, xyDEM *xydem, double Re,int heightFlag)
{
	double xi,yi;
	double x,y, rDist;
	int i,j, iSize, jSize;
	double t,u, p1, p2, p3, p4;
	double z;
	/*
	  Convert lat/lon to xy for dem
	*/
	lltoxy1(lat,lon,&x,&y,xydem->rot,xydem->stdLat);
	xi =  ((x-xydem->x0)/xydem->deltaX);
	yi =  ((y-xydem->y0)/xydem->deltaY);
	/*

	 */
	j = (int)xi;
	i = (int)yi;
	t = xi - (double)j;
	u = yi - (double)i;
	iSize = xydem->ySize; jSize= xydem->xSize;
	/*
	  If location of DEM return 0.
	*/
	if(i < 0 || i >= iSize || j < 0 || j >= jSize) {
		if(heightFlag == ELLIPSOIDAL) return 0;
		return earthRadius(lat*DTOR,EMINOR,EMAJOR)*KMTOM - Re;
	}
	/*
	  Border pixels
	*/
	if( i==(iSize-1) || j==(jSize-1) ) return xydem->z[i][j];
	/*
	  Set up interpolation   
	*/
	p1 = (double)xydem->z[i][j];  
	p2 = (double)xydem->z[i][j+1];
	p3 = (double)xydem->z[i+1][j+1];
	p4 = (double)xydem->z[i+1][j];    
	/*
	  Use bilinear interpolation to compute height.
	*/
	z = ((1.0 - t)*(1.0 - u) * p1 + t * (1.0 - u) * p2 +  t * u * p3 +                (1.0 - t) * u* p4);

	if(heightFlag == ELLIPSOIDAL)  return z;
	/*
	  Radial distance from earth center to point
	*/
	rDist = z + earthRadius(lat*DTOR,EMINOR,EMAJOR)*KMTOM;
	/*
	  Height reference to spherical earth
	*/
	z = rDist - Re;

	return z;
}
