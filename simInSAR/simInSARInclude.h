
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
  Initialize conversion matrices and constants for groundRangeToll
*/
void initGroundRangeToLLNew(inputImageStructure *inputImage);
/*
  Function to simulate InSAR image including both terrain and motion effects.
*/
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
