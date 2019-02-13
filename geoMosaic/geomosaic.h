#define HIGHJD 1.e7
#define LOWJD 0.0

/*
   Process input file for mosaicDEMs
*/
    void  processInputFile(char *inputFile,char ***insarDEMFiles, char ***demInputFiles,   outputImageStructure *outputImage, int *nDEMs);

    void makeGeoMosaic(inputImageStructure *inputImage,  outputImageStructure outputImage,
                       void *dem, int nFiles, int maxR, int maxA,   char **imageFiles,float fl, int smoothL, int orbitPriority);
     void  processInputFileGeo(char *inputFile,char ***insarDEMFiles, char ***demInputFiles,  outputImageStructure *outputImage, int *nDEMs,
                           float **weights, char ***antPatFiles);

#define RSATFINE -2
#define ALOS -3
