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
	insert_node(metaData, "lambda", STR_BUFF("%f", scene.I.par.lambda));
	insert_node(metaData, "dT", STR_BUFF("%f", scene.dT));
	if (scene.verticalCorrectionFile != NULL)
		insert_node(metaData, "verticalCorrection", STR_BUFF("%s", scene.verticalCorrectionFile));
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
	GDALDataType dataTypes[2] = {GDT_Float64, GDT_Float64};
	int32_t i, k;
	popuplateMeta(&metaData, scene);
	if (scene.tiffFlag) {
		/* GeoTIFF path: write .lat.tif + .lon.tif, then .ll.vrt via makeTiffVRT */
		char tifBuf1[2048], tifBuf2[2048];
		const char *tifFiles[2];
		char *tifSuffixes[2] = {".lat.tif", ".lon.tif"};
		char *tifBufs[2] = {tifBuf1, tifBuf2};
		float noDataArr[2] = {-2.e9f, -2.e9f};
		size_t nPx = (size_t)scene.aSize * scene.rSize;
		double *flat = (double *)malloc(nPx * sizeof(double));
		if (flat == NULL) error("outputLL: malloc failed for tif buffer\n");
		for (k = 0; k < 2; k++) {
			for (i = 0; i < scene.aSize; i++)
				memcpy(flat + (size_t)i * scene.rSize, images[k][i], scene.rSize * sizeof(double));
			tifFiles[k] = appendSuffix(outputFile, tifSuffixes[k], tifBufs[k]);
			fprintf(stderr, "writing tif %s\n", tifFiles[k]);
			writeFlatTiff(tifFiles[k], flat, scene.rSize, scene.aSize, GDT_Float64, -2.e9f, metaData);
		}
		free(flat);
		file = STR_BUFF("%s.ll.vrt", outputFile);
		makeTiffVRT(file, tifFiles, 2, noDataArr, metaData);
		free(file);
		return;
	}
	for(k=0; k < 2; k++) {
		file = appendSuffix(outputFile, suffixes[k], bandFiles[k]);
		fprintf(stderr, "writing %s\n", file);
		imageFP = fopen(file, "w");
		for (i = 0; i < scene.aSize; i++)
		{
			fwriteOptionalBS(images[k][i], scene.rSize, sizeof(double), imageFP, FLOAT64FLAG, scene.byteOrder);
		}
		fclose(imageFP);
	}
	file = STR_BUFF("%s.ll.vrt", outputFile);
	if(scene.byteOrder == MSB) byteSwapOption = "ByteOrder=MSB"; else byteSwapOption = "ByteOrder=LSB";
	writeSingleVRT(scene.rSize, scene.aSize, metaData, file, bandFiles, bandNames, dataTypes, byteSwapOption, -2.0e9, 2);
}

static void outputSimImage(sceneStructure scene, char *outputFile)
{
	dictNode *metaData = NULL;
	char *buf, buf1[2048], buf2[2048];
	char *file, *fileVRT;
	char *bandNames[1];
	char *bandFiles[1];
	char *byteSwapOption;
	int32_t i, j;
	GDALDataType dataTypes[1];
	FILE *imageFP;
	
	// Setup file name
	if(scene.maskFlag == TRUE) {
		if(scene.saveLLFlag == TRUE || scene.toLLFlag)
			file = appendSuffix(outputFile, ".mask", buf1);
		else
			file = appendSuffix(outputFile, "\0", buf1);
		// Malloc buff for conversion to byte
		buf = (char *)malloc(scene.rSize * sizeof(char));
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
	/* GeoTIFF path for the mask only: write <file>.tif + <file>.vrt, matching the
	   "<source>.tif" naming the Python writers use. The phase/height branch keeps the
	   raw path -- it feeds a different set of consumers that have no tiff support. */
	if (scene.tiffFlag && scene.maskFlag == TRUE)
	{
		char tifBuf[2048];
		const char *tifFile;
		const char *tifFiles[1];
		float noDataArr[1] = {0.0f};
		size_t nPx = (size_t)scene.aSize * scene.rSize;
		unsigned char *flat = (unsigned char *)malloc(nPx);
		if (flat == NULL)
			error("outputSimImage: malloc failed for tif buffer\n");
		for (i = 0; i < scene.aSize; i++)
			for (j = 0; j < scene.rSize; j++)
				flat[(size_t)i * scene.rSize + j] = (unsigned char)scene.image[i][j];
		tifFile = appendSuffix(file, ".tif", tifBuf);
		fprintf(stderr, "writing tif %s\n", tifFile);
		popuplateMeta(&metaData, scene);
		writeFlatTiff(tifFile, flat, scene.rSize, scene.aSize, GDT_Byte, 0.0f, metaData);
		free(flat);
		free(buf);
		tifFiles[0] = tifFile;
		fileVRT = STR_BUFF("%s.vrt", file);
		/* Name the band explicitly: the tif is <root>.mask.tif, so makeTiffVRT's
		   filename-derived name would be "mask", not the "Mask" the raw path
		   (writeSingleVRT below) writes and the Python readers ask for. */
		makeTiffVRTNamed(fileVRT, tifFiles, (const char **)bandNames, 1, noDataArr, metaData);
		free(fileVRT);
		return;
	}
	// Open file
	imageFP = fopen(file, "w");
	if (imageFP == NULL)
		error("*** outputSimulatedImage: Error opening %s ***\n", outputFile);
	/*
		Output image data
	*/		
	for (i = 0; i < scene.aSize; i++)
	{
		if (scene.maskFlag == FALSE)
		{	// Floating point cases
			fwriteOptionalBS(scene.image[i], scene.rSize, sizeof(float), imageFP, FLOAT32FLAG, scene.byteOrder);
		}
		else // Mask
		{
			// Convert data to byte
			for (j = 0; j < scene.rSize; j++)
				buf[j] = (unsigned char)scene.image[i][j];
			fwriteBS(buf, scene.rSize, sizeof(char), imageFP, BYTEFLAG);
		}
	}
	
	fprintf(stderr,"\n+\n");
	if (scene.maskFlag == TRUE) free(buf);
	fclose(imageFP);
	// Now write VRT
	fileVRT = STR_BUFF("%s.vrt", file); 
	// fileVRT = appendSuffix(file, ".vrt", buf2);
	popuplateMeta(&metaData, scene);
	writeSingleVRT(scene.rSize, scene.aSize, metaData, fileVRT, bandFiles, bandNames, dataTypes, byteSwapOption, -2.e9, 1);
}

/*
  Write scene.radiusImage (per-pixel azimuth-pixel smoothing radius, GDT_Byte) computed by
  computeSmoothRadiusMap(). Modeled on outputSimImage()'s mask-writing branch, since both are
  single-byte rasters.
*/
static void outputSmoothRadius(sceneStructure scene, char *outputFile)
{
	dictNode *metaData = NULL;
	char *file, buf1[2048];
	char *fileVRT;
	char *bandNames[1];
	char *bandFiles[1];
	GDALDataType dataTypes[1];
	FILE *imageFP;
	int32_t i;

	popuplateMeta(&metaData, scene);
	if (scene.tiffFlag)
	{
		char tifBuf[2048];
		const char *tifFile;
		const char *tifFiles[1];
		float noDataArr[1] = {0.0f};
		size_t nPx = (size_t)scene.aSize * scene.rSize;
		unsigned char *flat = (unsigned char *)malloc(nPx);
		if (flat == NULL)
			error("outputSmoothRadius: malloc failed for tif buffer\n");
		for (i = 0; i < scene.aSize; i++)
			memcpy(flat + (size_t)i * scene.rSize, scene.radiusImage[i], scene.rSize);
		tifFile = appendSuffix(outputFile, ".smr.tif", tifBuf);
		fprintf(stderr, "writing tif %s\n", tifFile);
		writeFlatTiff(tifFile, flat, scene.rSize, scene.aSize, GDT_Byte, 0.0f, metaData);
		free(flat);
		tifFiles[0] = tifFile;
		fileVRT = STR_BUFF("%s.smr.vrt", outputFile);
		makeTiffVRT(fileVRT, tifFiles, 1, noDataArr, metaData);
		return;
	}

	file = appendSuffix(outputFile, ".smr", buf1);
	fprintf(stderr, "writing %s\n", file);
	imageFP = fopen(file, "w");
	if (imageFP == NULL)
		error("*** outputSmoothRadius: Error opening %s ***\n", file);
	for (i = 0; i < scene.aSize; i++)
		fwriteBS(scene.radiusImage[i], scene.rSize, sizeof(unsigned char), imageFP, BYTEFLAG);
	fclose(imageFP);

	bandFiles[0] = file;
	bandNames[0] = "SmoothRadius";
	dataTypes[0] = GDT_Byte;
	fileVRT = STR_BUFF("%s.vrt", file);
	writeSingleVRT(scene.rSize, scene.aSize, metaData, fileVRT, bandFiles, bandNames, dataTypes, NULL, -2.e9, 1);
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
	//fprintf(stderr, "%i %i %i\n", scene.saveLLFlag, scene.toLLFlag, scene.maskFlag );
	//if( (scene.saveLLFlag == FALSE && scene.toLLFlag) )
	// If either saveLL (geodat grid) or toLLFlag (offset grid)
	if(scene.saveLLFlag == TRUE || scene.toLLFlag == TRUE) 
	{
		outputLL(scene, outputFile);
		// If either mask or height flag set, output these variables
		fprintf(stderr, "%i %i\n", scene.maskFlag, scene.heightFlag);
	   	if(scene.maskFlag == TRUE || scene.heightFlag == TRUE)
		{
			fprintf(stderr, "Output height or mask...\n");
			outputSimImage(scene, outputFile);
	   	}	
	} else {
		fprintf(stderr, "Output...");
		outputSimImage(scene, outputFile);
	}
	if (scene.smoothRadiusFlag == TRUE)
	{
		fprintf(stderr, "Output smoothing-radius map...\n");
		outputSmoothRadius(scene, outputFile);
	}
	//
	//  Form header filename by adding .simdat suffix to outputFile
	//
	buf1[0] = '\0';
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
