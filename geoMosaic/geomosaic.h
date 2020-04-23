#define HIGHJD 1.e7
#define LOWJD 0.0
/* -30dB min for S1 */
#define MINS1DB -30.0
#define MINS1SIG 0.001
/*
   Process input file for mosaicDEMs
*/
    void  processInputFile(char *inputFile,char ***insarDEMFiles, char ***demInputFiles,   outputImageStructure *outputImage, int *nDEMs);

    void makeGeoMosaic(inputImageStructure *inputImage,  outputImageStructure outputImage,
                       void *dem, int nFiles, int maxR, int maxA,   char **imageFiles,float fl, int smoothL,int smoothOut, int orbitPriority, float ***psiBuf, float ***gBuf);
     void  processInputFileGeo(char *inputFile,char ***insarDEMFiles, char ***demInputFiles,  outputImageStructure *outputImage, int *nDEMs,
                           float **weights, char ***antPatFiles);

#define RSATFINE -2
#define ALOS -3
