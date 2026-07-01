
/*#include "ers1/getLocC_p/parfile.h"*/

#define DEFAULT_SIM_BN 0.0
#define DEFAULT_SIM_BP 0.0
#define LL 1
#define LR 2
#define UL 3
#define UR 4


/*
  #define REMAJOR 6378.137
  #define REMINOR 6356.7523142
  #define THETAC 20.355
  #define H_ALT 789.00
*/
#define SLANTRANGEDEM 3
#define NORMALDEM 1

typedef struct latLonPairType
{
	double lat;
	double lon;
} latLonPair;

typedef struct xyPairType
{
	double x;
	double y;
} xyPair;

typedef struct sceneStructureType
{
	int32_t flatFlag;
	int32_t heightFlag;
	int32_t maskFlag;
	int32_t offsetFlag;
	int32_t toLLFlag;
	int32_t saveLLFlag;
	char *llInput;
	double deltaBn;
	double deltaBp;
	double bn;
	double bnStart;
	double bnEnd;
	double bnStep;
	double bp;
	double bpStart;
	double bpEnd;
	double bpStep;
	double dT;
	/* size info added 2/3/17 */
	int32_t aSize;
	int32_t rSize;
	float aO;
	float rO;
	float dR;
	float dA;
	inputImageStructure I;
	ShelfMask *imageMask;
	double width;
	double length;
	float **image;
	double **latImage;
	double **lonImage;
	int32_t useVelocity;
	int32_t byteOrder;
	float velThresh;      /* > 0 enables velThresh mode */
	char  *geodat2File;   /* second geodat for SV-derived baseline; NULL = use explicit bn/bp */
	float *bnArray;       /* per-azimuth bn from SV; NULL when not used */
	float *bpArray;       /* per-azimuth bp from SV; NULL when not used */
	int32_t tiffFlag;     /* -tiff: write lat/lon as GeoTIFF instead of binary */
	char  *verticalCorrectionFile; /* -verticalCorrection vcFile; NULL = no correction */
	xyDEM *verticalCorrection;     /* loaded from verticalCorrectionFile, NULL if not used */
	/* Variable smoothing-radius map (-minTol/-percentSpeed/-maxTol) */
	int32_t smoothRadiusFlag;      /* TRUE if -minTol/-percentSpeed/-maxTol given */
	double minTol;                 /* m/yr floor for adaptive tolerance */
	double percentSpeed;           /* percent (e.g. 1.0 = 1%) of local speed for adaptive tolerance */
	double maxTol;                 /* m/yr ceiling for adaptive tolerance */
	int32_t maxSmoothRadius;       /* pixels; sweep cap, clamped to <= 255 */
	int32_t smoothNIter;           /* repeated box-filter passes per sweep step (Gaussian-ish) */
	float **toleranceImage;        /* per-pixel native-units (radian) tolerance, same size as image */
	unsigned char **radiusImage;   /* per-pixel azimuth-pixel smoothing radius, output of computeSmoothRadiusMap */
} sceneStructure;

typedef struct displacementStructureType
{
	int32_t coordType;
	int32_t size1;
	int32_t size2;
	double minC1;
	double minC2;
	double maxC1;
	double maxC2;
	double deltaC1;
	double deltaC2;
	double **dR;
} displacementStructure;

/*
  parse scene input file for siminsar.
*/
void parseSceneFile(char *sceneFile, sceneStructure *scene);
/*
  Compute per-azimuth bn/bp arrays in scene from two sets of state vectors.
  Call after parseSceneFile so scene->aSize is set.
*/
void simInSARBaselineFromSV(char *geodat2File, sceneStructure *scene);

/*
  Initialize conversion matrices and constants for groundRangeToll
*/
void initGroundRangeToLLNew(inputImageStructure *inputImage);
/*
  Function to simulate InSAR image including both terrain and motion effects.
*/
void simInSARDEMBounds(sceneStructure *scene, double rot, double stdLat,
                       double *xMin, double *xMax, double *yMin, double *yMax);
void simInSARimage(sceneStructure *scene, void *dem, xyVEL *xyVel);
/*
  Output simulated image. Writes two files one for image, and xxx.simdat
  with image header info
*/
void outputSimulatedImage(sceneStructure scene, char *outputFile, char *demFile, char *displacementFile);
/*
  Input displacement map
*/
void getDisplacementMap(char *displacementFile, displacementStructure *displacements, demStructure dem, sceneStructure scene);
/*
  Compute displacement for lat/lon from Displacement map.
*/
double getDisplacement(double lat, double lon, displacementStructure *displacements);
/*
  Input slant range dem.
*/
void getSlantRangeDEM(char *demFile, demStructure *dem, sceneStructure scene);
/*
  Compute speed-threshold mask from velocity map and write as GeoTIFF.
*/
void velThreshMask(xyVEL *xyVel, float velThresh, char *outputFile);
/*
  Compute scene->radiusImage from scene->image and scene->toleranceImage (see
  computeSmoothRadius.c for the sweep algorithm). Only meaningful when
  scene->smoothRadiusFlag is set.
*/
void computeSmoothRadiusMap(sceneStructure *scene);
