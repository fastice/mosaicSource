
typedef struct landSatImageType {
	matchResult matches;
	lsFit fitResult;
	char *maskFile;
	double weight;
	struct  landSatImageType *next;
	int32 nImages;
} landSatImage;


landSatImage  *parseLSInputs(char *inputFile, landSatImage *LSImages, double jd1, double jd2,  int timeOverlapFlag);
void makeLandSatMosaic(landSatImage *LSImages,outputImageStructure *outputImage,float fl);
