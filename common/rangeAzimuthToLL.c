#include "mosaicSource/common/common.h"

int rangeAzimuthToLL(double *rg, double range,double iFloat,double rhoSp,double ReH, double Re, double *lat, double *lon,double *hWGS,inputImageStructure *inputImage,xyDEM *xyDem,double tol,double step ) {
	double R0,ReWGS,rhoWGS, drgStep, R;
	int recycle, nIter;
	recycle=FALSE;
	nIter=0;
	while(1) {
		/* For rg on spherical earth, compute lat/lon for the corresponding ellipsoidal earth 
		   and return the corresponding slant range */
		R0=groundRangeToLLNew(*rg,iFloat,lat,lon,inputImage,recycle);
		recycle=TRUE; /* flags to reuse previously computed values that are still valid */
		/* Compute WGS radious at lat/lon */
		ReWGS=earthRadius(*lat*DTOR,EMINOR,EMAJOR)*KMTOM;
		/* Comput rho for the ellipsodal radius */
		rhoWGS=rhoRReZReH( R0,(ReWGS+0), ReH);
		/*  Get height for lat/lon referenced to the ellipsoid */
		*hWGS= getXYHeight(*lat,*lon,xyDem,Re, ELLIPSOIDAL);					
		/*  Compute height corrected R    	*/
		R=slantRange( rhoWGS, ReWGS+ (*hWGS),  ReH);
		/* Not the most efficient iteration, but good enough */
		drgStep=(range-R)*step;
	/*	fprintf(stderr,"R0 %9.1f %9.1f  %9.1f  %8.3f %5.0f\n",R0,R,range,drgStep,*hWGS);		*/
		if(R > (range-tol))  break;				
		/* Speed convergence if way out */
		*rg += drgStep;    /* Increment ground range */
		/* Avoid infinite loop */
		if(nIter > 8000){fprintf(stderr,"%f %f %f\n",Re,rhoSp,ReH); break; }
		nIter++;
	} /* Endwhile */

	
	return nIter;
}

	
