#include "stdio.h"
#include "string.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include <libgen.h>
#define STR_BUFFER_SIZE 1024
#define STR_BUFF(fmt, ...) ({                                    \
    char *__buf = (char *)calloc(STR_BUFFER_SIZE, sizeof(char)); \
    snprintf(__buf, STR_BUFFER_SIZE, fmt, ##__VA_ARGS__);        \
    __buf;                                                       \
})


static void popuplateMeta(dictNode **metaData, sceneStructure scene) {
	insert_node(metaData, "r0", STR_BUFF("%i", (int) (scene.rO * scene.I.nRangeLooks)));
	insert_node(metaData, "a0", STR_BUFF("%i", (int) (scene.aO * scene.I.nAzimuthLooks)));
	insert_node(metaData, "deltaR", STR_BUFF("%i", (int)(scene.dR * scene.I.nRangeLooks)));
	insert_node(metaData, "deltaA", STR_BUFF("%i",  (int)(scene.dA * scene.I.nAzimuthLooks)));
	insert_node(metaData, "sigmaRange", STR_BUFF("%f", 0.0));
	insert_node(metaData, "sigmaStreaks", STR_BUFF("%f", 0.0));
}
static void outputLL(sceneStructure scene, char *outputFile) 
{
	dictNode *metaData = NULL;
	FILE *imageFP;
	char *file, buf1[2048], buf2[2048], bufvrt[2048];
	double **images[2] = {scene.latImage, scene.lonImage};
	char *byteSwapOption;
	char *suffixes[2] = {".lat", ".lon"};
	char *bandNames[2] = {"lat", "lon"};
	char *bandFiles[2] = {buf1, buf2};
	int32_t dataTypes[2] = {GDT_Float64, GDT_Float64};
	int32_t i, k;

	for(k=0; k < 2; k++) {
		file = appendSuffix(outputFile, suffixes[k], bandFiles[k]);
		fprintf(stderr, "writing %s\n", buf1);
		imageFP = fopen(file, "w");
		for (i = 0; i < scene.aSize; i++)
		{
			fwriteOptionalBS(images[k][i], scene.rSize, sizeof(double), imageFP, FLOAT64FLAG, scene.byteOrder);
		}
		fclose(imageFP);
	}
	file = appendSuffix(outputFile, ".ll.vrt", bufvrt);
	if(scene.byteOrder == MSB) byteSwapOption = "ByteOrder=MSB"; else byteSwapOption = "ByteOrder=LSB";
	popuplateMeta(&metaData, scene);
	writeSingleVRT(scene.rSize, scene.aSize, metaData, file, bandFiles, bandNames, dataTypes, byteSwapOption, 2);
}

static void outputSimImage(sceneStructure scene, char *outputFile)
{
	dictNode *metaData = NULL;
	char *buf, buf1[2048], buf2[2048];
	char *file, *fileVRT;
	char *bandNames[1];
	char *bandFiles[1];
	char *byteSwapOption;
	int32_t i, j, dataTypes[1];
	FILE *imageFP;
	// Setup file name
	if(scene.maskFlag == TRUE && (scene.saveLLFlag == TRUE || scene.toLLFlag)) {
		file = appendSuffix(outputFile, ".mask", buf1);
		bandNames[0] = "Mask";
		dataTypes[0] = GDT_Byte;
		byteSwapOption = NULL;
	}
	else
	{
		file = outputFile;
		bandNames[0] = "Phase";
		if(scene.heightFlag == TRUE) bandNames[0] = "Height" ;
		dataTypes[0] = GDT_Float32;
		if(scene.byteOrder == MSB) byteSwapOption = "ByteOrder=MSB"; else byteSwapOption = "ByteOrder=LSB";
	}
	bandFiles[0] = file;
	// Open file
	imageFP = fopen(file, "w");
	if (imageFP == NULL)
		error("*** outputSimulatedImage: Error opening %s ***\n", outputFile);
	/*
		Output image data
	*/
	if (scene.maskFlag == TRUE)
		buf = (char *)malloc(scene.rSize * sizeof(char));
	for (i = 0; i < scene.aSize; i++)
	{
		if (scene.maskFlag == FALSE)
		{	// Floating point cases
			fwriteOptionalBS(scene.image[i], scene.rSize, sizeof(float), imageFP, FLOAT32FLAG, scene.byteOrder);
		}
		else // Mask
		{
			// Convert data to byte
			dataTypes[0] = GDT_Byte;	
			for (j = 0; j < scene.rSize; j++)
				buf[j] = (char)scene.image[i][j];
			fwriteBS(buf, scene.rSize, sizeof(char), imageFP, BYTEFLAG);
		}
	}
	/*  free(buf);*/
	if (scene.maskFlag == TRUE) free(buf);
	fclose(imageFP);
	
	fileVRT = appendSuffix(file, ".vrt", buf2);
	fprintf(stderr, "VRT File: %s\n", fileVRT);
	popuplateMeta(&metaData, scene);
	writeSingleVRT(scene.rSize, scene.aSize, metaData, fileVRT, bandFiles, bandNames, dataTypes, byteSwapOption, 1);
}

/*
  Output simulated image. Writes two files one for image, and xxx.simdat
  with image header info
*/
void outputSimulatedImage(sceneStructure scene, char *outputFile, char *demFile, char *displacementFile)
{
	FILE *imageFP, *imageDatFP;
	int32_t i, j;
	int32_t pSize;
	char *outputFileDat, *tmp;
	char *buf, buf1[2048];
	/*
	  Open image outputfile
	*/
	if( (scene.saveLLFlag == FALSE && scene.toLLFlag) || scene.maskFlag == TRUE)
		outputSimImage(scene, outputFile);
	if(scene.saveLLFlag == TRUE || scene.toLLFlag == TRUE) outputLL(scene, outputFile);
	//
	//  Form header filename by adding .simdat suffix to outputFile
	//
	outputFileDat = appendSuffix(outputFile, ".simdat", buf1);
	// Open image outputfile
	imageDatFP = fopen(outputFileDat, "w");
	if (imageDatFP == NULL)
		error("*** outputSimulatedImage: Error opening %s ***\n", outputFileDat);
	// Output image header info.
	fprintf(imageDatFP, "# 2\n;\n;  Image size (pixels) nx ny \n;\n");
	fprintf(imageDatFP, "%i  %i\n", scene.I.rangeSize, scene.I.azimuthSize);
	fprintf(imageDatFP, "; Baseline start/end\n");
	fprintf(imageDatFP, "%f  %f\n&\n", scene.bnStart, scene.bnEnd);
	fprintf(imageDatFP, "; demFile :        %s\n", demFile);
	fprintf(imageDatFP, "; displacemtnFile: %s\n", displacementFile);
	fclose(imageDatFP);
	return;
}
