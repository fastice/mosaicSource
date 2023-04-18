#include <stdio.h>

typedef struct landSatImageType
{
	matchResult matches;
	lsFit fitResult;
	char *maskFile;
	double weight;
	struct landSatImageType *next;
	int32_t nImages;
} landSatImage;

landSatImage *parseLSInputs(char *inputFile, landSatImage *LSImages, double jd1, double jd2, int32_t timeOverlapFlag);
void makeLandSatMosaic(landSatImage *LSImages, outputImageStructure *outputImage, float fl);
