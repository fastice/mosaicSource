#include "common.h"

/*
  Find region where image exists
*/
void getRegion(inputImageStructure *image, int *iMin,int *iMax, int *jMin,int *jMax,  outputImageStructure *outputImage)
{
	extern double Rotation;
	double xa1,ya1;
	double minX,maxX,minY,maxY;
	double pad;
	int i;
	/*
	  Loop through points to find max and min locations
	*/    
	for(i=1; i <=4; i++) {
		/* changed 10/6/05 to avoid detailed setup and memory leak
		   lltoxy1(image->cpAll.latControlPoints[i],
		   image->cpAll.lonControlPoints[i],&xa1,&ya1,
		   Rotation,outputImage->slat);
		   fprintf(stderr,"%f %f %f %f\n",image->cpAll.latControlPoints[i],image->cpAll.lonControlPoints[i],xa1,ya1);
		*/
		lltoxy1(image->latControlPoints[i],image->lonControlPoints[i],&xa1,&ya1,Rotation,outputImage->slat);
		/* Commented out print statements 7/31/2015 */
		fprintf(stderr,"%f %f %f %f\n",image->latControlPoints[i],image->lonControlPoints[i],xa1,ya1);
		fprintf(stderr,"%f %f\n",Rotation,outputImage->slat);
		if(i==1) {minX=xa1; maxX=xa1; minY=ya1; maxY=ya1;}
		else {
			minX=min(minX,xa1); maxX=max(maxX,xa1);
			minY=min(minY,ya1); maxY=max(maxY,ya1);
		}
	}

	fprintf(stderr,"%f %f %f %f\n",minX, maxX, minY, maxY);
	/*
	  Compute i,j min,max with pad
	*/ 
	pad=15000.; 
	*iMin=(int)((minY*KMTOM-outputImage->originY - pad)/outputImage->deltaY);
	*jMin=(int)((minX*KMTOM-outputImage->originX - pad)/outputImage->deltaX);
	*iMax=(int)((maxY*KMTOM-outputImage->originY + pad)/outputImage->deltaY);
	*jMax=(int)((maxX*KMTOM-outputImage->originX + pad)/outputImage->deltaX);
	*iMin=max(*iMin,0); *jMin=max(*jMin,0);
	*iMax=min(outputImage->ySize,*iMax); *jMax=min(outputImage->xSize,*jMax);
}
