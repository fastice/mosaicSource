#include "math.h"
#include "common.h"

/*
       This file contains interpAzOffsets, interpAzSigma, interpRangeOffsets, interpRangeSigma, which are used to interpolate the offsets fields.
*/

static void computeCoordsForInterp(inputImageStructure *inputImage ,Offsets *offsets, double *range,double *azimuth,double *imageLength,double *normAzimuth);

/*
   Interpolate azimuth offsets map for velocity generation and apply baseline/geometry corrections
*/
float interpAzOffset(double range,double azimuth,Offsets *offsets, inputImageStructure *inputImage, double Range,double theta, float azSLPixSize)
{
	float result;
	float zeroOffset;
	double alongTrack;
	double imageLength, normAzimuth;
	double azimuthML, rangeML;
	/*
	  Compute azimuth in image coord stuff
	*/
	azimuthML=azimuth; /* compute coords will convert input azimuth to range offsets, so save */
	rangeML=range; /* compute coords will convert input azimuth to range offsets, so save */		
	computeCoordsForInterp(inputImage,offsets,&range,&azimuth,&imageLength,&normAzimuth);
        alongTrack=normAzimuth*offsets->doffdx;
	/*
		Interpolate 
	*/
	result = bilinearInterp((float **)offsets->da,range,azimuth,offsets->nr,offsets->na,-0.9999*LARGEINT,(float)-LARGEINT);
	if (result < -0.9999*LARGEINT) return -LARGEINT;
	if(offsets->deltaB != DELTABNONE) {
		zeroOffset=svAzOffset(inputImage,offsets,rangeML,azimuthML) ; /* Offset from SV */
	} else zeroOffset=0.0;
	
	zeroOffset += (float)offsets->c1 + Range*sin(theta)*offsets->dbcds -   Range * cos(theta)*offsets->dbhds;
	/* Apply scaling corrections */
	if(inputImage->lookDir==LEFT) result *=-1.0;
	result *= azSLPixSize;
	result -=zeroOffset;
	result -=alongTrack;

	return result;

}


/*
   Interpolate azimuth sigma  map for velocity generation
*/
float interpAzSigma(double range,double azimuth,Offsets *offsets,  inputImageStructure *inputImage, double Range,double theta, float azSLPixSize)
{
	float result;
	double alongTrack;
	double imageLength, normAzimuth;
	/*
	  Compute azimuth in image coord stuff
	*/
	computeCoordsForInterp(inputImage,offsets,&range,&azimuth,&imageLength,&normAzimuth);
	alongTrack=normAzimuth*offsets->doffdx;	
	/*
	 Interpolate
	*/
	result = bilinearInterp((float **)offsets->sa,range,azimuth,offsets->nr,offsets->na,-0.99999*LARGEINT,(float)-LARGEINT);
	/*if (result < -0.99999*LARGEINT) return -LARGEINT;	*/
	/* Added 10/20/2017 to cap errors at 1/10 of a pixel - sigma streaks could make it larger */
	if(result > 0.1 || result < 0) result=0.1;		
	result = sqrt(result*result + offsets->sigmaStreaks *offsets->sigmaStreaks);
	result *= azSLPixSize;
	return result;
}

/*
   Interpolate range offset  map for velocity generation and apply baseline/geometry corrections
*/
float interpRangeOffset(double range,double azimuth,Offsets *offsets,   inputImageStructure *inputImage,
			double Range,double thetaD, float rSLPixSize, double theta, double *demError)
{
	float result;
	float zeroOffset;
	double bn,bp, bSq;
	double bnS,bpS;
	double  xsq;
	double imageLength, normAzimuth, azimuthML;
	/*	  Compute azimuth in image coord stuff	*/
	azimuthML=azimuth; /* compute coords will convert input azimuth to range offsets, so save */
	computeCoordsForInterp(inputImage,offsets,&range,&azimuth,&imageLength,&normAzimuth);
	/* 
	   Baseline or deltaBaseline (SV case)
	*/
	xsq=normAzimuth*normAzimuth;
	bn = offsets->bn + offsets->dBn * normAzimuth + offsets->dBnQ * xsq;
	bp = offsets->bp + offsets->dBp * normAzimuth + offsets->dBpQ * xsq;
	/* 
	 State vector baseline ?
	*/
	if(offsets->deltaB != DELTABNONE) {
		svInterpBnBp(inputImage,offsets,azimuthML,&bnS,&bpS);
		bn = bnS + bn;
		bp = bpS + bp;
	}
	bSq = bn*bn + bp*bp;
	/*
	  Note this can be derived from Eq 7, JGlac 1996, page 566. Solution for quadratic equation
	*/   
	zeroOffset = sqrt( pow(Range,2.0) -2.0*Range*(bn*sin(thetaD) + bp*cos(thetaD)) + bSq)  - Range + offsets->rConst;
	/*
 	if(range > 0 && range < inputImage->rangeSize && azimuth > 0 && azimuth < inputImage->azimuthSize) {
	fprintf(stderr,"NEW zereOffsets,thetaD,bn,bp %f %f %f %f %f %f %f\n",zeroOffset,thetaD*RTOD,bn,bp,range,azimuth,Range);
	error("stop");
	}*/
	/* Changed 3/1/16 from 30 to 15 for the nominal DEM error */
	*demError= fabs(bn) * 15.0 /(Range*sin(theta));
	/*
	  Interpolate 
	*/
	result = bilinearInterp((float **)offsets->dr,range,azimuth,offsets->nr,offsets->na,-0.9999*LARGEINT,(float)-LARGEINT);
	if (result < -0.9999*LARGEINT) return -LARGEINT;
	/*    
	   apply corrections 
	*/
	result *= rSLPixSize; /* This puts offset in meters */
	result -=zeroOffset;  /* Substract the offset in meters */
	return result;
}

/*
   Interpolate range sigma  map for velocity generation
*/
float interpRangeSigma(double range,double azimuth,Offsets *offsets,  inputImageStructure *inputImage,   double Range,double thetaD, float rSLPixSize)
{
	float result;
	double imageLength, normAzimuth;
	/*
	  Compute azimuth in image coord stuff
	*/
	computeCoordsForInterp(inputImage,offsets,&range,&azimuth,&imageLength,&normAzimuth);
	/*
	 Interpolate
	*/
	result = bilinearInterp((float **)offsets->sr,range,azimuth,offsets->nr,offsets->na,-0.9999*LARGEINT,(float)-LARGEINT);
	
	/*if (result < -0.9999*LARGEINT) return -LARGEINT;	*/
	/* 
		If sigma range is set, it represents a minimum error 
	*/
	if(result > 0.1 || result < 0.0) result=0.1;			
	result=max(offsets->sigmaRange,result);
	result *= rSLPixSize; /* put in units of meters */
	return result;
}

static void computeCoordsForInterp(inputImageStructure *inputImage ,Offsets *offsets,
							    double *range,double *azimuth,double *imageLength,double *normAzimuth) {
	/*
	  Compute azimuth in image coord stuff
	*/
	/* added 10/30/2013 to fix azimuth problem -((inputImage->nAzimuthLooks-1)*0.5)  . */
	*imageLength = (double)inputImage->azimuthSize;   
	*normAzimuth = (*azimuth-0.5* (*imageLength))/ (*imageLength);
	*range=  (*range * inputImage->nRangeLooks     - ((inputImage->nRangeLooks  -1)*(-0.5)) -offsets->rO)/offsets->deltaR;
	*azimuth=(*azimuth * inputImage->nAzimuthLooks - ((inputImage->nAzimuthLooks-1)*(-0.5)) -offsets->aO)/offsets->deltaA;
	return;
}
