#include <math.h>
/*#include "mosaicSource/GeoCodeDEM_p/geocodedem.h"*/
#include "cRecipes/nrutil.h"
#include "common.h"
/* 
	Save kernel memory to recyle on repeated calls 
*/
double **rDistSave=NULL;
/*
  Compute scale for feathering.
*/
void  computeScale(float **inImage,float **scale,  int azimuthSize,int rangeSize,float fl,  float weight,double minVal)
{
	extern double **rDistSave;
	double **rDist,rA; 			/* Radial distance kernel */
	int j,k,is,il,s1,s2,l1,l2;		/* Loop variabiles */
	float minV;
	/* 
	No feathring case - set weight to all 1's
	*/
	if(fl == 0 ) {
		initFloatMatrix(scale,azimuthSize,rangeSize,1.0);
		return;
	 }
	/*
	  Set up featherin - radial distance from center of kernel.
	  Allocate space if first time, or recycle previoius
	*/
	if(rDistSave==NULL)  rDistSave = dmatrix(-fl, fl,-fl,fl);
	rDist = rDistSave;
	/* 
	   Compute radial distance for kernel 
	*/
	for(is = 0; is <= fl; is++) {
		for(il = 0; il <= fl; il++) {
			rA=weight * min(max( sqrt((double)(is * is) + (double)(il * il)),0.5 )/fl,1);
			rDist[il][is]  = rA; rDist[il][-is] = rA;
			rDist[-il][is] = rA; rDist[-il][-is] = rA;           
		} /* end for(il... */
	} /* end for(is... */
	/*
	  Set initial value for scale array with specified weight value
	*/
	initFloatMatrix(scale,(long int) azimuthSize,(long int)rangeSize,weight);
	/*
	  Now loop through weights
	*/
	for(j=0; j <  azimuthSize; j++) {
		for(k=0; k <  rangeSize; k++) {
			/*
			  Adding the feathering at valid points.
			*/ 
			if(inImage[j][k] > minVal ) {
				/* 
					Find edge pixels - check all neighbors and if one is non-valid, its an edge 
				*/
				minV=1.0e30;  /* will get reduced if its an edge */				
				s1=max(0,j-1); s2= min( azimuthSize-1,j+1);
				l1=max(0,k-1); l2= min( rangeSize-1,k+1);
				/* find minV for point and its neigbors */
				for(is=s1; is <= s2; is++) 
					for(il=l1; il <= l2; il++) {
						minV=FMIN(minV,(float)(inImage[is][il]));
					}
				/* 
					if minV <= minVal its a border pixel, so set scale to distance from border. Use minvalue where kernels overlap.
				 */
				if( (minV <= minVal)  ) { 
					s1=max(0,j-fl); s2= min( azimuthSize-1,j+fl);
					l1=max(0,k-fl); l2= min( rangeSize-1,k+fl);
					/* set the scale to the rad dist, or a less value of already set */
					for(is=s1; is <= s2; is++) {
						for(il=l1; il <= l2; il++) {
							scale[is][il] = 	FMIN(rDist[is-j][il-k], scale[is][il]);
						} /* il */
					} /* is */
				} /* End if(minV... */
			} /* End if(inp... */
		} /* k */
	} /* j */
}
  
