#include "stdio.h"
#include"string.h"
#include <math.h>
#include "mosaicSource/common/common.h"
/*
  Compute shift for registering images using corner points.
*/
    void coarseRegister(inputImageStructure *inputImage1,inputImageStructure *inputImage2)
{   
   double lat,lon;
   double range1,range2;
   double azimuth1,azimuth2;
   double deltaR,deltaA;
   int i,nPts;
   fprintf(stderr,"coarse reg\n");


   fprintf(stderr,"Setup conversions 1\n");
   initllToImageNew(inputImage1); /* Setup conversions */
   fprintf(stderr,"Setup conversions 2\n");
   initllToImageNew(inputImage2); /* Setup conversions */
   fprintf(stderr,"Compute locations\n");
   deltaR=0.0; deltaA=0.0; nPts=0;
   for(i = 0; i < 5 ; i++ ) {
	   lat=inputImage1->latControlPoints[i];
	   lon=inputImage1->lonControlPoints[i];
	   llToImageNew(lat,lon,0,&range1,&azimuth1,inputImage1);
	   llToImageNew(lat,lon,0,&range2,&azimuth2,inputImage2);
	   /*
	   fprintf(stderr,"lat,lon %f %f\n",lat,lon);
	   fprintf(stderr,"r2-r1 %f a2-a1 %f\n",(range2-range1)*inputImage1->nRangeLooks,(azimuth2-azimuth1)*inputImage1->nAzimuthLooks);*/
	   /*	   fprintf(stderr,"r1,a1 %f %f\n",range1,azimuth1);
		   fprintf(stderr,"r2,a2 %f %f\n",range2,azimuth2);*/
	   if(range2 > 0 && range2 < inputImage2->rangeSize & azimuth2 > 0 && azimuth2 < inputImage2->azimuthSize) {
		   deltaR += (range2-range1)*inputImage1->nRangeLooks;
		   deltaA += (azimuth2-azimuth1)*inputImage1->nAzimuthLooks;
		   nPts++;
	   }
   }
   deltaR = deltaR/(double)nPts;
   deltaA = deltaA/(double)nPts;   
   /*  
      Took out 8/30/00
     if(inputImage1->passType==DESCENDING) deltaA *=-1;
   */
   fprintf(stdout,"%i %i\n",(int)deltaR,(int) deltaA); 
}




