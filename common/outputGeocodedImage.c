#include "stdio.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "mosaicSource/common/common.h"

size_t fwriteOptionalBS(void *ptr, size_t nitems, size_t size, FILE *fp, int32_t flags, int32_t byteOrder) {
   if(byteOrder == LSB) return fwrite(ptr, nitems, size, fp);
   else return fwriteBS(ptr, nitems, size, fp, flags);
}

/*
   Output geocode image. Writes two files one for image, and xxx.geodat
   with image header info

*/
void outputGeocodedImage(outputImageStructure outputImage, char *outputFile)
{
    FILE *imageFP, *imageDatFP;
    int32_t i;
    int pSize;
    char *outputFileDat, *tmp;
    /*
        Open image outputfile
    */
    imageFP = fopen(outputFile, "w");
    if (imageFP == NULL)
        error("*** outputGeocodedImage: Error opening %s ***\n", outputFile);
    /*
        Output image data
    */
    /*    if(outputImage.imageType == COMPLEX) pSize =sizeof(ers1Complex);
        else */
    pSize = sizeof(float);
    for (i = 0; i < outputImage.ySize; i++)
        if (outputImage.imageType == COMPLEX)
            fwriteBS(outputImage.image[i], outputImage.xSize * pSize, 1, imageFP, INT16FLAG);
        else
            fwriteBS(outputImage.image[i], outputImage.xSize * pSize, 1, imageFP, FLOAT32FLAG);
    fclose(imageFP);
    /*
        Form header filename by adding .geodat suffix to outputFile
    */
    tmp = (char *)malloc((size_t)512);
    tmp[0] = '\0';
    tmp = strcat(tmp, outputFile); /* Compute input filename */
    tmp = strcat(tmp, ".geodat");
    outputFileDat = (char *)malloc((size_t)(strlen(tmp) + 1));
    outputFileDat = strcpy(outputFileDat, tmp);
    free(tmp);
    /*
        Open image outputfile
    */
    imageDatFP = fopen(outputFileDat, "w");
    if (imageDatFP == NULL)
        error(
            "*** outputGeocodedImage: Error opening %s ***\n", outputFileDat);
    /*
        Output image header info.
    */
    fprintf(imageDatFP, "# 2\n;\n;  Image size (pixels) nx ny \n;\n");
    fprintf(imageDatFP, "%i  %i\n", outputImage.xSize, outputImage.ySize);
    fprintf(imageDatFP, ";\n;  Pixel size (m) deltaX deltaY \n;\n");
    fprintf(imageDatFP, "%f  %f\n", outputImage.deltaX, outputImage.deltaY);
    fprintf(imageDatFP, ";\n;  Origin, lower left corner (km) Xo  Yo \n;\n");
    fprintf(imageDatFP, "%f  %f\n&\n",
            outputImage.originX * MTOKM, outputImage.originY * MTOKM);
    fclose(imageDatFP);
    return;
}


void outputGeocodedImageTiff(outputImageStructure outputImage, char *outputFile, char *driverType, const char *epsg,
                            dictNode *summaryMetaData, float noDataValue, int32_t dataType)
{
    double geoTransform[6];
    char *outputFileTiff;
    if(hasSuffix(outputFile, ".tif"))
    {
        outputFileTiff = outputFile;
    }
    else 
    {
        outputFileTiff = appendSuffix(outputFile, ".tif", (char *)malloc(strlen(outputFile) + 5));
    }
    
    if (outputImage.imageType == COMPLEX) 
    {
        error("Tiff output for complex not currently supported");
    }
    // Compute Geotransform
	computeGeoTransform(geoTransform, outputImage.originX, outputImage.originY, outputImage.xSize,
					    outputImage.ySize, outputImage.deltaX, outputImage.deltaY);
	// Set up meta data
	char *timeStamp = timeStampMeta();
	insert_node(&summaryMetaData, "CreationTime", timeStamp);
    // Write to tiff file
    saveAsGeotiff(outputFileTiff, (float *)outputImage.image[0], outputImage.xSize,
				outputImage.ySize, geoTransform, epsg, summaryMetaData, driverType, dataType, noDataValue);
}