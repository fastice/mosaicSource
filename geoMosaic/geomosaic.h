#define HIGHJD 1.e7
#define LOWJD 0.0
/* -30dB min for S1 */
#define MINS1DB -30.0
#define MINS1SIG 0.001
/* geoMosaicMode values */
#define GEOMOSAIC_AVERAGE 0
#define GEOMOSAIC_MIN     1
#define GEOMOSAIC_MAX     2
/*
   Process input file for mosaicDEMs
*/
void processInputFile(char *inputFile, char ***insarDEMFiles, char ***demInputFiles, outputImageStructure *outputImage, int *nDEMs);

void makeGeoMosaic(inputImageStructure *inputImage, outputImageStructure outputImage,
                   void *dem, int nFiles, int maxR, int maxA, char **imageFiles, float fl, int smoothL, int smoothOut, int orbitPriority, float ***psiBuf, float ***gBuf);
void processInputFileGeo(char *inputFile, char ***insarDEMFiles, char ***demInputFiles, outputImageStructure *outputImage, int *nDEMs,
                         float **weights, char ***antPatFiles);

void readComplexAsPower(inputImageStructure *inputImage, int32_t yMin, int32_t yMax);
void  AbAg(double x, double y, double azimuth,
           inputImageStructure *inputImage, outputImageStructure *outputImage,
           void *dem, double *Ab, double *Ag, double *shadow, int32_t recycle);
float applyCorrections(float *value, inputImageStructure *inputImage,
                       double range, double azimuth, double h);
int32_t subPixelGammaRTC(double x, double y,
                          inputImageStructure *inputImage,
                          outputImageStructure *outputImage,
                          void *dem,
                          float *power, double *Ab, double *Ag, float *psiE);
int32_t subPixelGammaRTCLinear(double x, double y,
                                inputImageStructure *inputImage,
                                outputImageStructure *outputImage,
                                void *dem,
                                float *power, double *Ab, double *Ag, float *psiE);
int32_t subPixelGammaRTCJacobian(double x, double y,
                                  inputImageStructure *inputImage,
                                  outputImageStructure *outputImage,
                                  void *dem,
                                  float *power, double *Ab, double *Ag, float *psiE);

#define RSATFINE -2
#define ALOS -3
