#include "stdio.h"
#include"string.h"
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "mosaicSource/common/common.h"
/*
   Output geocode image. Writes two files one for image, and xxx.geodat 
   with image header info
*/
    void outputGeocodedImage(outputImageStructure outputImage,char *outputFile)
{   
    FILE *imageFP, *imageDatFP;
    int i;
    int pSize;
    char *outputFileDat, *tmp;
/*
    Open image outputfile
*/
    imageFP = fopen(outputFile,"w");
        if(imageFP == NULL) 
            error("*** outputGeocodedImage: Error opening %s ***\n",outputFile);
/*
    Output image data
*/
/*    if(outputImage.imageType == COMPLEX) pSize =sizeof(ers1Complex);
    else */
    pSize = sizeof(float);
    for(i=0; i < outputImage.ySize; i++)  
         if(outputImage.imageType == COMPLEX) fwriteBS(outputImage.image[i],outputImage.xSize*pSize,1,imageFP,INT16FLAG);
         else fwriteBS(outputImage.image[i],outputImage.xSize*pSize,1,imageFP,FLOAT32FLAG);
    fclose(imageFP);
/*
    Form header filename by adding .geodat suffix to outputFile
*/
    tmp = (char *) malloc((size_t)512); 
    tmp[0] = '\0'; 
    tmp = strcat(tmp,outputFile);   /* Compute input filename */
    tmp = strcat(tmp,".geodat");
    outputFileDat = (char *) malloc((size_t)(strlen(tmp)+1));
    outputFileDat = strcpy(outputFileDat,tmp);
    free(tmp);
/*
    Open image outputfile
*/
    imageDatFP = fopen(outputFileDat,"w");
        if(imageDatFP == NULL) error(
            "*** outputGeocodedImage: Error opening %s ***\n",outputFileDat);
/*
    Output image header info.
*/
    fprintf(imageDatFP,"# 2\n;\n;  Image size (pixels) nx ny \n;\n");
    fprintf(imageDatFP,"%i  %i\n",outputImage.xSize,outputImage.ySize);
    fprintf(imageDatFP,";\n;  Pixel size (m) deltaX deltaY \n;\n");
    fprintf(imageDatFP,"%f  %f\n",outputImage.deltaX,outputImage.deltaY);
    fprintf(imageDatFP,";\n;  Origin, lower left corner (km) Xo  Yo \n;\n");
    fprintf(imageDatFP,"%f  %f\n&\n",
            outputImage.originX*MTOKM,outputImage.originY*MTOKM);
    fclose(imageDatFP);
    return;
}
