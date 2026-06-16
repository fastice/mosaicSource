#include <math.h>
#include "common.h"

/*
  Compute flat-earth baseline phase for ISCE/NISAR products where topography
  has already been removed. Applies the same formula as changeflat/reflattenImage.c
  (flat-earth term only), so mosaic3d can correct the orbit-error ramp inline
  without a changeflat intermediate file.
*/
void computePhiFlatEarth(double *phiFlat, double azimuth, vhParams *vhParam,
						 inputImageStructure *phaseImage, double Range,
						 double ReHfixed, double Re, double thetaCfixed, double *phaseError)
{
	double normAzimuth, imageLength;
	double bn, bp, bSq, delta;
	double thetaFlat, thetaDFlat;
	double twok;
	double xsq;

	twok = 4.0 * PI / phaseImage->par.lambda;
	imageLength = (double)phaseImage->azimuthSize;
	normAzimuth = (azimuth - 0.5 * imageLength) / imageLength;
	xsq = normAzimuth * normAzimuth;
	bn = bPoly(vhParam->Bn, vhParam->dBn, vhParam->dBnQ, normAzimuth);
	bp = bPoly(vhParam->Bp, vhParam->dBp, vhParam->dBpQ, normAzimuth);
	bSq = bn * bn + bp * bp;

	/* Flat-earth look angle — identical to computePhiZ.c:29 and changeflat:145 */
	thetaFlat = thetaRReZReH(Range, Re, ReHfixed);
	thetaDFlat = thetaFlat - thetaCfixed;

	/* Flat-earth range delta only (no terrain-height term) */
	delta = -bn * sin(thetaDFlat) - bp * cos(thetaDFlat) + bSq * 0.5 / Range;
	*phiFlat = delta * twok;

	/* PI/4 unwrap error — same as computePhiZ when Bn=0 */
	*phaseError = PI / 4.0;
}
