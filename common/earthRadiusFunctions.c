#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "float.h"
#include "common.h"

/*
  Compute radial distance from earth center to elipsoid surface
  lat should be in radians
  rp = rminor
  re=rmajor
*/
double earthRadius(double lat, double rp, double re)
{
	double earthRad;
	double n, x, z;
	/*    fprintf(stderr,"rp,re,lat %f %f %f",rp,re,lat); */
	n = (re * re) / sqrt(pow(re * cos(lat), 2.0) + pow(rp * sin(lat), 2.0));
	x = n * cos(lat);
	z = pow((rp / re), 2.0) * n * sin(lat);
	earthRad = sqrt(x * x + z * z);

	return earthRad;
}

/*
  earth radius for WGS84
 */
double earthRadiusWGS84(double lat)
{
	double rp = EMINOR, re = EMAJOR;
	return earthRadius(lat, rp, re);
}

/*
  Curvature for earth at lat
 */
double earthRadiusCurvatureWGS84(double lat)
{
	double earthRad;
	double n, x, z;
	double rp = EMINOR, re = EMAJOR;
	n = (re * re) / sqrt(pow(re * cos(lat), 2.0) + pow(rp * sin(lat), 2.0));
	return n;
}
