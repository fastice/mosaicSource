#include <math.h>
#include "common.h"
/*
  Compute phase (phiZ) for a given elevation. z looks like its in thetaD, z not actually used??
*/
void computePhiZ(double *phiZ, double azimuth, vhParams *vhParam, inputImageStructure *phaseImage, double thetaD, double Range,
				 double ReH, double ReHfixed, double Re, double thetaCfixed, double *phaseError)
{
	double normAzimuth, imageLength;
	double bn, bp, bSq, delta;
	double xsq;
	double twok;
	double thetaFlat, thetaDFlat;
	twok = 4.0 * PI / phaseImage->par.lambda;
	/* geometry stuff */
	imageLength = (double)phaseImage->azimuthSize;
	normAzimuth = (azimuth - 0.5 * imageLength) / imageLength;
	xsq = normAzimuth * normAzimuth;
	/* Baseline */
	bn = bPoly(vhParam->Bn, vhParam->dBn, vhParam->dBnQ, normAzimuth);
	bp = bPoly(vhParam->Bp, vhParam->dBp, vhParam->dBpQ, normAzimuth);
	bSq = bn * bn + bp * bp;
	/* Compute delta and theta; complete phase */
	delta = sqrt(pow(Range, 2.0) - 2.0 * Range * (bn * sin(thetaD) + bp * cos(thetaD)) + bSq) - Range;
	/*
	  Phase has already had flat earth removed, so remove that.
	  This step needs to use the same azimuth-independent thetaC used in changeflat.
	*/
	thetaFlat = thetaRReZReH(Range, (Re + 0), ReHfixed);
	thetaDFlat = thetaFlat - thetaCfixed;
	delta -= -bn * sin(thetaDFlat) - bp * cos(thetaDFlat) + bSq * 0.5 / Range;
	/* Changed 06/06/07 to make baseline dependent phase error */
	/* Assumes PI/4 phase error, 30 elevation error */
	*phaseError = sqrt(pow((double)(PI / 4.), 2.0) + pow((double)((fabs(bn) * 30.0 / (Range * sin(thetaFlat))) * twok), 2.0));
	/* Scale delta to phase */
	*phiZ = delta * twok;
	return;
}
