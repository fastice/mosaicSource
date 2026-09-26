#define HIGHJD 1.e7
#define LOWJD 0.0
/* -30dB min for S1 */
#define MINS1DB -30.0
#define MINS1SIG 0.001
/* geoMosaicMode values */
#define GEOMOSAIC_AVERAGE 0
#define GEOMOSAIC_MIN     1
#define GEOMOSAIC_MAX     2
/* S1Cal output bit flags */
#define PSISAVE 2
#define GAMMACORSAVE 4
#define GAMMASAVE 8
/* -calOutput values: which calibrated quantities to output */
#define CALOUTPUT_BOTH   0
#define CALOUTPUT_SIGMA0 1
#define CALOUTPUT_GAMMA0 2
/* -nearRange / -farRange: keep only inputs whose incidence angle is within angleTolerance
   of the per-pixel minimum (NEAR) or maximum (FAR) over all inputs */
#define RANGESELECT_NONE 0
#define RANGESELECT_NEAR 1
#define RANGESELECT_FAR  2
/* coarse incidence buffer: a cell no input reached filters nothing */
#define INCUNSET -999.0
/* Already geocoded NISAR GCOV inputs, from the -gcov yaml file */
typedef struct
{
	char polarization[16]; /* covariance term, e.g. HHHH */
	char frequency[4];	   /* A or B */
	int32_t useMask;	   /* drop samples flagged invalid/fill in the GCOV mask */
	/* Directory holding rtcGammaToSigmaFactor outside the granule (yaml key factorFrom).
	   Empty (the default) means the factor lives in the granule, as it does in an archive
	   product. Slim products share one factor per track/frame/grid across cycles; the
	   downloader leaves a per-granule symlink here named exactly like the granule, so this
	   code only ever opens <factorDir>/<granule basename> and never has to derive the key. */
	char factorDir[2048];
	int32_t nFiles;
	char **files;
	float *weights;
} gcovInputs;
/*
   Process input file for mosaicDEMs
*/
void processInputFile(char *inputFile, char ***insarDEMFiles, char ***demInputFiles, outputImageStructure *outputImage, int *nDEMs);

void makeGeoMosaic(inputImageStructure *inputImage, outputImageStructure outputImage,
                   void *dem, int nFiles, int maxR, int maxA, char **imageFiles, float fl, int smoothL, int smoothOut, int orbitPriority, float ***psiBuf, float ***gBuf,
                   gcovInputs *gcov);
/* Coarse incidence-angle buffer shared by the range-selection passes. Indices are ABSOLUTE
   output-lattice cell indices (see incCoarseIndex), so a tiled run and an untiled run agree. */
typedef struct
{
	float **inc;      /* nY x nX, INCUNSET where no input reached the cell */
	int32_t nX, nY;   /* cell counts covering the output grid */
	int32_t stride;   /* output pixels per cell */
	int32_t j0, i0;   /* absolute cell index of the grid's first cell */
} incBuffer;

int32_t incCoarseIndex(double origin, double delta, int32_t stride);
void incUpdate(incBuffer *incBuf, int32_t i1, int32_t j1, double psi);
float incLookup(incBuffer *incBuf, int32_t i1, int32_t j1);
int32_t incWeight(incBuffer *incBuf, int32_t i1, int32_t j1, double psi);
void gcovIncidenceCoarse(gcovInputs *gcov, outputImageStructure *outputImage, void *dem, incBuffer *incBuf);
/* cell span covering output pixels [kMin,kMax] and the output pixel at a cell's centre; isX
   selects the axis. The centre may fall outside the grid, which is what keeps a cell's value
   independent of how the mosaic is tiled. */
void incCellSpan(incBuffer *incBuf, int32_t isX, int32_t kMin, int32_t kMax, int32_t *cMin, int32_t *cMax);
int32_t incCellCentre(incBuffer *incBuf, int32_t isX, int32_t c);
void readGCOVYaml(char *yamlFile, gcovInputs *gcov);
void gcovBounds(gcovInputs *gcov, double *minX, double *maxX, double *minY, double *maxY);
int32_t gcovToOutputGrid(gcovInputs *gcov, int32_t iFile, outputImageStructure *outputImage, void *dem,
                         float **imageTmp, float **scaleTmp, float **psiBufTmp, float **gBufTmp,
                         inputImageStructure *gcovImage, int32_t *imageDate,
                         int32_t *iMin, int32_t *iMax, int32_t *jMin, int32_t *jMax,
                         incBuffer *incBuf, unsigned char **selTmp, int32_t pass);
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
