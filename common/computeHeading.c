#include <math.h>
#include <stdlib.h>
#include "common.h"

/*
   Compute the track heading at lat/lon/z
   Note this gives the angle of the cross track (not along track)
   direction from north (CW).
*/
double computeHeading(double lat, double lon, double z, inputImageStructure *inputImage,
					  conversionDataStructure *cP)
{
	extern int32_t HemiSphere;
	double range1, range2, azimuth1, azimuth2;
	double r1, r2;
	double rho1, rho2;
	double zSp1, zSp2;
	double dlat;
	double gRange1, gRange2;
	double da, dgr;
	double hAngle;
	double Re, RNear;
	double ReH;
	int32_t index;

	/* modified 7/31/02 - made dlat symmetrical,+- 0.1 instead of +0.22,  and change rang limits to 500 - this should avoid zero rhos */
	if (inputImage->lookDir == LEFT)
		dlat = -0.0025;
	else
		dlat = 0.0025;
	Re = cP->Re;
	RNear = cP->RNear;
	z = 0.0;
	zSp1 = z - earthRadius((lat - dlat) * DTOR, EMINOR, EMAJOR) * KMTOM + cP->Re;
	zSp2 = z - earthRadius((lat + dlat) * DTOR, EMINOR, EMAJOR) * KMTOM + cP->Re;
	/* Compute range/azimuth coords for point dlat to north and south of lat */
	llToImageNew(lat - dlat, lon, z, &range1, &azimuth1, inputImage);
	llToImageNew(lat + dlat, lon, z, &range2, &azimuth2, inputImage);
	/* Handle case where along track variation in lookup */
	if (cP->ReH != NULL)
	{
		index = min(max(0, (int)azimuth1), cP->azSize - 1);
		ReH = cP->ReH[index];
	}
	else
	{ /*ReH = cP->Re + cP->H; */
		error("computeHeading - obsolte ");
	}
	if (range1 < -2000 || range1 > (inputImage->rangeSize + 2000))
		return (9999.);
	/*
	  Compute azimuth displacement for latitude displacement
	*/
	da = (azimuth2 - azimuth1) * inputImage->azimuthPixelSize;
	/*
	  Compute ground range  displacement for latitude displacement
	  Do computations on a spherical earth
	*/
	r1 = RNear + inputImage->rangePixelSize * range1;
	rho1 = rhoRReZReH(r1, (Re + zSp1), ReH);
	/*	rho1 = acos( min((r1*r1 - ReH*ReH - (Re+zSp1)*(Re+zSp1) ) / (-2.0*(Re+zSp1)*ReH ),1.) );*/
	gRange1 = Re * rho1;
	r2 = RNear + inputImage->rangePixelSize * range2;
	rho2 = rhoRReZReH(r2, (Re + zSp2), ReH);
	/*rho2 = acos( min((r2*r2 - ReH*ReH - (Re+zSp2)*(Re+zSp2) ) / (-2.0*(Re+zSp2)*ReH ),1.) );*/
	gRange2 = Re * rho2;
	/* delta ground range */
	dgr = gRange2 - gRange1;
	/*
	  Compute track heading
	*/
	if (inputImage->lookDir == RIGHT)
	{
		hAngle = atan2(da, dgr);
	}
	else if (inputImage->lookDir == LEFT)
	{
		hAngle = atan2(da, -dgr);
	}
	else
		error("*** computeHeading: invalid lookDir ***");
	return hAngle;
}
