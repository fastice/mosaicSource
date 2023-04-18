#include "stdio.h"
#include "string.h"
#include "math.h"
#include <stdlib.h>
#include "geotiff/xtiffio.h"   /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "clib/standard.h"
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"

double xyscale(double latctr, int32 proj)
/*

  Use this scale factor to correct distance

  converted from old idl function

  ;  This function calculates the scaling factor for a polar stereographic
  ;  projection (ie. SSM/I grid) to correct area calculations. The scaling
  ;  factor is defined (from Snyder, 1982, Map Projections used by the U.S.
  ;  Geological Survey) as:
  ;
  ;    k = (mc/m)*(t/tc), where:
  ;
  ;    m = cos(lat)/sqrt(1 - e2*sin(lat)^2)
  ;    t = tan(Pi/4 - lat/2)/((1 - e*sin(lat))/(1 + e*sin(lat)))^(e/2)
  ;    e2 = 0.006693883 is the earth eccentricity (Hughes ellipsoid)
  ;    e = sqrt(e2)
  ;    mc = m at the reference latitude (70 degrees)
  ;    tc = t at the reference latitude (70 degrees)
  ;
  ;  The ratio mc/tc is precalculated and stored in the variable m70_t70.
  ;
*/
{
	double clat, lat, m70_t70, k, m, e2, e, slat, t;
	/*
	  if( not keyword_set(south)) then m70_t70 = 1.9332279d0 $
	  else m70_t70=1.9390299d  ; for 71 deg
	*/
	if (latctr < 0)
		latctr = -latctr; /* Ensure postive value - assumes you wouldn't use PS in wrong hemisphere */
	if (proj == NSIDCNORTH)
	{
		m70_t70 = 1.9332279;
	}
	else if (proj == NSIDCSOUTH)
	{
		m70_t70 = 1.9390299;
	}
	else
		error("xyscale: Invalid proj %i\n", proj);
	/*e2 = 0.006693883d0 HUGHES */
	/* WGS 84 https://en.wikipedia.org/wiki/Geodetic_datum */
	e2 = 0.00669437999014;
	e = sqrt(e2);
	lat = (latctr * DTOR);
	slat = sin(lat);
	clat = cos(lat);
	m = clat / sqrt(1.0 - e2 * slat * slat);
	t = tan(PI / 4.0 - lat / 2.0) / pow((1.0 - e * slat) / (1.0 + e * slat), (e / 2.0));
	k = m70_t70 * t / m;
	return (1.0 / k);
}
