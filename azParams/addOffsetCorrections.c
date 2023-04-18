#include <math.h>
#include "mosaicSource/common/common.h"
#include "azparams.h"

/*
  Compute motion correction.
*/
void addOffsetCorrections(inputImageStructure *inputImage, tiePointsStructure *tiePoints)
{
	extern int32_t HemiSphere;
	int32_t i, j;
	conversionDataStructure *cP;
	double lat, lon;
	double va;
	double hAngle, xyAngle, rotAngle;
	double x1, y1;
	double deltaT;

	initllToImageNew(inputImage);
	cP = &(inputImage->cpAll);
	deltaT = tiePoints->nDays / 365.25;
	/*
	  Loop over tie points
	*/
	for (i = 0; i < tiePoints->npts; i++)
	{
		lat = tiePoints->lat[i];
		lon = tiePoints->lon[i];
		/*
		  Compute angle for CW rotation of ps to radar
		*/
		hAngle = computeHeading(lat, lon, 0, inputImage, cP);
		computeXYangleNoDem(lat, lon, &xyAngle, tiePoints->stdLat);
		rotAngle = hAngle - xyAngle;
		/*
		  Compute vyra by ccw rot from ps to ra
		*/
		va = tiePoints->vx[i] * sin(rotAngle) + tiePoints->vy[i] * cos(rotAngle);
		/*
		  Compute displacement
		*/
		tiePoints->phase[i] -= va * deltaT;
	}

	return;
}
