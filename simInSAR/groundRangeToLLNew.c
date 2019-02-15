#include "stdio.h"
#include"string.h"
#include <math.h>
#include <stdlib.h>
#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"


static void getState(double myTime,inputImageStructure *inputImage, double *xs,double *ys, double *zs, double *vsx, double *vsy, double *vsz) {
	int n;
	stateV *sv;  /* statevectors from input image */
	
	sv=&(inputImage->sv);
	n=(long int)((myTime - sv->times[1])/(sv->deltaT)+.5);
	n=min(max(0,n-2), sv->nState-NUSESTATE);	
	polintVec(&(sv->times[n]), &(sv->x[n]), &(sv->y[n]),&(sv->z[n]),  &(sv->vx[n]), &(sv->vy[n]),&(sv->vz[n]), myTime, xs,ys ,zs, vsx,vsy ,vsz);
}
/*
  Compute params for groundrange/azimuth to slantrange.
  Adapted from asf code by Shusun Li.

*/
static double xs,ys,zs,vsx,vsy,vsz;
double  groundRangeToLLNew( double groundRange, float azimuth,double *lat,double *lon,inputImageStructure *inputImage,int recycle)
{ 
	conversionDataStructure *cp;

	double rho, ReH,h,rsl;
	double x,y;
	double myTime;
	int i;
	int azIndex, rIndex,azimuthInt;
	cp = &(inputImage->cpAll);
	azimuthInt=min(  max( 0,(int)(azimuth+0.5) ),   inputImage->azimuthSize-1);
	myTime=(azimuth* inputImage->nAzimuthLooks)/cp->prf + cp->sTime;
	if(recycle == FALSE) 
		getState( myTime,inputImage, &xs,&ys, &zs, &vsx, &vsy, &vsz);
		
/*	fprintf(stderr,"%f %f %f %f %f %f %f \n",myTime, xs,ys, zs, vsx, vsy, vsz);*/
	rho = groundRange/cp->Re;
	ReH=cp->ReH[azimuthInt];
	/* Locate point on with that range on an ellipsoidal earth */
	h=0;
	rsl = sqrt( (ReH*ReH) + (cp->Re+h)*(cp->Re+h) - 2.0*(cp->Re+h)*(ReH)*cos(rho) );
	smlocateZD( xs*MTOKM,  ys*MTOKM,  zs*MTOKM, vsx*MTOKM,vsy*MTOKM, vsz*MTOKM, rsl*MTOKM,lat,lon,(double)(inputImage->lookDir),0.0);
	return rsl;
/*	fprintf(stderr,"%f %f \n",*lat,*lon);*/
}



