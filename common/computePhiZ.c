#include <math.h>
#include "common.h"
/*
  Compute phase (phiZ) for a given elevation. z looks like its in thetaD, z not actually used??
*/
void computePhiZ(double *phiZ,double z,  double azimuth,vhParams *vhParam,  inputImageStructure *phaseImage, double thetaD,double Range,
		 double ReH,double ReHfixed, double Re, double thetaC,double *phaseError)
{
	double normAzimuth,imageLength;
	double bn, bp,bSq, delta;
	double xsq;
	double theta, thetaDFlat;
	double lambda;
	lambda=phaseImage->par.lambda;
	/* geometry stuff */
	imageLength = (double)phaseImage->azimuthSize;
	normAzimuth = (azimuth-0.5*imageLength)/imageLength;
	xsq = normAzimuth*normAzimuth;
	/* Baseline */
	bn = bPoly( vhParam->Bn,  vhParam->dBn, vhParam->dBnQ, normAzimuth);
	bp = bPoly( vhParam->Bp,  vhParam->dBp, vhParam->dBpQ, normAzimuth);
	bSq = bn*bn + bp*bp;
	/* Compute delta and theta */
	delta = sqrt(  pow(Range,2.0) -2.0*Range*(bn*sin(thetaD) + bp*cos(thetaD)) + bSq  )  - Range;
	theta=thetaRReZReH(Range, (Re+0), ReHfixed);
	/*theta = acos( (Range*Range + ReHfixed*ReHfixed - Re*Re)/( 2.0*(ReHfixed)*Range ) );*/
	/* Changed 06/06/07 to make baseline dependent phase error */
	*phaseError = sqrt( pow((double)(PI/4.),2.0) + pow( (double) ((fabs(bn)*30.0 /(Range*sin(theta))) *4.0*PI/lambda),2.0) );
	/*
	  update flatten delta
	*/
	thetaDFlat = theta - thetaC;
	delta -=  -bn * sin(thetaDFlat) - bp * cos(thetaDFlat) + bSq * 0.5 / Range;
	/* Scale delta to phase */     
	*phiZ = delta * 4.0 * PI / lambda;
	return;
}

    
