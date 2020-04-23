#include <stdio.h>
#include "clib/standard.h"
#include "cRecipes/cRecipes.h"
#include "geocode.h"
#define MAXOFFBUF 72000000
#define MAXOFFLENGTH 30000

#define MINELEVATION -1000
#define MINVELOCITY 5
#define NUSESTATE 5
#define GROUNDED 0
#define SHELF 1
#define GROUNDINGZONE 2
#define NOSOLUTION 5
/* elipsoid */
#define F 1.0d/298.257222563
#define EMINOR 6356.7523142
#define EMAJOR 6378.137
/*
  Constants
*/
#define RANGEONLY 0
#define RANGEANDAZIMUTH 1
#define LARGEINT 2000000000
#define DESCENDING 0
#define ASCENDING 1
#define INIT 2
#define UPDATE 3
#define RETURN 0 
#define CONTINUE 1
#define OK 0
#define MAXLAYERS 10000
#define COMPLEX 0
#define POWER 1
#define AMPLITUDE 2
#define CORRELATION 3
#define PHASE 4
#define FLOAT 5
#define UNWRAPPED 6
#define VZDEFAULT 0
#define VZHORIZONTAL 1
#define VZVERTICAL 2
#define VZLOS 3
#define VZINC 4
#define  AZONLY 0
#define RGANDAZ 1
#define  RGONLY 2
/* default values - not really used */
#define RANGESIZE 2048
#define AZIMUTHSIZE 12800
#define AZIMUTHPIXELSIZE 3.8999999
#define RANGEPIXELSIZE 7.8999
#define BUFFERSIZE 8388608
#define LAMBDAERS1 0.05656
#define MAXAZIMUTHLINES 200000
#define DELTABQUAD 2
#define DELTABCONST 1
#define DELTABNONE 0
extern double RangePixelSize;    /* Size in m of range pixel */
extern double AzimuthPixelSize;    /* Size in m of azimuth pixel */



typedef struct irregType {
	char *file;
	double *x;
	double *y;
	double *vx;
	double *vy;
	double *vz;
	int *link;
	int nData;
	int nTri;
	double maxLength;
	double maxArea;
	struct irregType *next;
} irregularData;


typedef struct OffsetsType {
	int nr;
	int na;
	int rO;
	int aO;
	float deltaA;
	float deltaR;
	/* Azimuth params */
	double c1;
	double dbcds;
	double dbhds;
	double doffdx;
	/* Range params */
	double bn;
	double bp;
	double dBn;
	double dBp;
	double dBnQ;
	double dBpQ;
	double rConst;
	double Cr[7][7];
	double Ca[7][7];
	/* Pixels spacing */
	float **da;
	float **dr;
	float **sa;
	float **sr;
	char *file;
	char *azParamsFile;
	char *rFile;
	char *rParamsFile;
	double sigmaStreaks;
	double sigmaRange;
	double sigmaRresidual;
	double sigmaAresidual;
	char *geo1;
	char *geo2;
	stateV sv2;
	double dt1t2;
	int deltaB;
	float *bnS;
	float *bpS;
	double rOffS;
	double azFit[10];
	int azInit;
} Offsets;

typedef struct vhParamsType {
	double nDays;
	int noSlant;
	int vh3D;
	int acrossAngle;
	double angleLimit;
	double Bn; 
	double Bp;
	double dBn; 
	double dBp; 
	double dBnQ; 
	double dBpQ;
	double C[7][7];
	double sigma;
	int offsetFlag; 
	int rOffsetFlag; 
	char *xyDEMFile;
	char *altXYDEMFile;
	Offsets offsets;
	xyDEM xydem;
	struct vhParamsType *next;
} vhParams;


#define MAXTIEPOINTS 500000

typedef struct tiePointsType {
	int npts;
	int noRamp;
	int dBpFlag;
	int quadB;
	int motionFlag;
	int vrFlag;
	int imageCoords;
	int bnbpFlag;    /* Flag to estimate only the center baseline components*/
	int bpFlag;      /* Flag to estimate only the bp center component */
	int bnbpdBpFlag;
	int bpdBpFlag;
	int timeReverseFlag;
	int quiet; /* Don't echo tiepoints to baseline solution */
	double BnCorig;      /* Original params */
	double BpCorig;
	double dBnorig;
	double dBporig;
	double dBnQorig;
	double dBpQorig;
	double deltaB;
	double stdLat;
	double nDays;
	double H;
	double Re;
	double thetaC;
	double RNear;
	double *bsq; /* Modified 12/09/94 to vector */
	double *lat;
	 double *lon;
	double *x;
	double *y;
	double *z;
	double *r;
	double *a;
	double *phase;
	double *delta;
	double *vx;
	double *vy;
	double *vz;
	double *vyra;
	double cnstA; /* const dead reckoned az offset */
	double cnstR;	/* const dead reckoned rg offset */
	/*
	  Parameters below are used exclusively for azParams.
	*/
	int constOnlyFlag;
	int linFlag;
} tiePointsStructure;


typedef struct unwrapPhaseType {
	int rangeSize;
	int azimuthSize;
	float **phase;
} unwrapPhaseStructure;

typedef struct modelValuesType {
	double x;
	double thetaD;
} modelValues;

typedef struct aZmodelValuesType {
    double r;
    double theta;
    double x;
} aZmodelValues;


/*
  Process input file for mosaicDEMs

  void  processInputFile(char *inputFile,char ***insarDEMFiles, char ***demInputFiles,  outputImageStructure *outputImage, int *nDEMs);
*/

/*
  Make mosaic of insar and other dems.
*/
void makeDEMMosaic(inputImageStructure *inputImage,  outputImageStructure outputImage, void *dem, int nDEMs, int maxR, int maxA,  char **insarDEMFiles);
/*
  Compute scale for feathering.
*/
void  computeScale(float **inImage,float **scale,  int azimuthSize,int rangeSize,float fl,  float weight,double minVal);
void fillRadialKernel(double **rDist, float fl, float weight);
/*
  void  getInputFile(char *inputFile,char ***insarDEMFiles,  char ***demInputFiles,  outputImageStructure *outputImage,  float **weights, int *nDEMs);
*/
float interpolateInSARDEM( inputImageStructure *slantRangeDEM, double **yImage, double gRange, double azimuth, double deltaGroundRange );
double computeHeading(double lat,double lon,double z, inputImageStructure *inputImage, conversionDataStructure *cP);

/* Initialize a floating point matrix */
void initFloatMatrix(float **x,long int nr,long int nc,float initValue);
int getDataStringSpecial(FILE *fp, int lineCount, char *line,int *eod, char special, int *specialFound);
void getBaseline(char *baselineFile, vhParams *params);
void rotateFlowDirectionToXY(double drn,double dan, double *dxn, double *dyn, double xyAngle, double hAngle);
void rotateFlowDirectionToRA(double dxn,double dyn,  double *dan, double *drn,  double xyAngle, double hAngle);
void xyGetZandSlope(double lat, double  lon,  double x, double y, double *zSp,double *zWGS84, double *da, double *dr,
		    conversionDataStructure *cP,vhParams *vhParam, inputImageStructure *slantRangeDEM);

void computeXYvh(inputImageStructure *phaseImage, inputImageStructure *slantRangeDEM, vhParams *vhParam, outputImageStructure outputImage );
void computeXYdata(double lat,double lon, double *dxn, double *dyn,double *xyAngle,  xyDEM xydem,double *slopeMag);
void errorsToXY(double er,double ea, double *ex, double *ey,  double xyAngle, double hAngle);
double interpXYDEM(double x, double y, xyDEM xydem);
void getRegion(inputImageStructure *image, int *iMin,int *iMax, int *jMin,int *jMax,  outputImageStructure *outputImage);
unsigned char getShelfMask(ShelfMask *shelfMask,double x, double y);
void computePhiZ(double *phiZ,double z,  double azimuth,vhParams *vhParam,  inputImageStructure *phaseImage, double thetaD,double Range,	 double ReH,double ReHfixed, double Re, double thetaC,double *phaseError);
double interpTideDiff(double x, double y, xyDEM xydem);
void  getIrregData(irregularData *irregData);
void addIrregData(irregularData *irregDat,  outputImageStructure *outputImage, float fl);
void  parseIrregFile(char *irregFile,irregularData **irregDat);
void  getMVhInputFile(char *inputFile,char ***phaseFiles, char ***geodatFiles, char ***baselineFiles,   char ***offsetFiles, char ***azParamsFiles, char ***rOffsetFiles,char ***rParamsFiles, 
		      outputImageStructure *outputImage,  float **nDays, float **weights, int **crossFlags, int *nFiles,   int offsetFlag, int rOffsetFlag, int threeDOffFlag);
float interpAzOffset(double range,double azimuth,Offsets *offsets,  inputImageStructure *inputImage,  double Range,double theta, float azSLPixSize);
float interpAzSigma(double range,double azimuth,Offsets *offsets,  inputImageStructure *inputImage, double Range,double theta, float azSLPixSize);
float interpRangeOffset(double range,double azimuth,Offsets *offsets,  inputImageStructure *inputImage,double Range,double thetaD, float rSLPixSize, double theta, double *demError);
float interpRangeSigma(double range,double azimuth,Offsets *offsets,  inputImageStructure *inputImage, double Range,double thetaD, float rSLPixSize);
float bilinearInterp(float **fimage,double range,double azimuth,int nr,int na,float minvalue, float noData);

void computeXYangle(double lat,double lon,double *xyAngle, xyDEM xydem);
void computeXYangleNoDem(double lat,double lon,double *xyAngle,double stdLat);
	void endScale(outputImageStructure *outputImage,float **vXimage,float **vYimage,float **vZimage,float **errorX,float **errorY, float **scaleX,float **scaleY,float **scaleZ,int statsFlag);
void readBothOffsets( Offsets *offsets);

void undoNormalization(outputImageStructure *outputImage,float **vXimage,float **vYimage,float **vZimage,float **errorX,float **errorY, float **scaleX,float **scaleY,float **scaleZ,float **fScale,int statsFlag);
void redoNormalization(float myWeight,outputImageStructure *outputImage,int iMin,int iMax,int jMin,int jMax,
		       float **vXimage,float **vYimage,float **vZimage,float **errorX,float **errorY, float **scaleX,float **scaleY,float **scaleZ,float **fScale,
       		       float **vxTmp,float **vyTmp,float **vzTmp, float **sxTmp,float **syTmp,  int statsFlag);

void setupBuffers(outputImageStructure *outputImage, float ***vXimage, float ***vYimage, float ***vzImage,float ***scaleX,float ***scaleY, float ***scaleZ,
		  float ***vxTmp, float ***vyTmp, float ***vzTemp,float ***sxTmp,float ***syTmp, float ***fScale,float ***errorX,float ***errorY);

conversionDataStructure *setupGeoConversions (inputImageStructure *currentImage, float *azSLPixSize, float *rSLPixSize,double *Re, double *ReH,double *thetaC,double *ReHfixed,double *thetaCfixed);
float shelfMaskCorrection(inputImageStructure *currentImage, vhParams *currentParams,unsigned char sMask,double x,double y,double psi, double*sigmaR);

void getAzParams( Offsets *offsets);
void getRParams( Offsets *offsets);
void readOffsetDataAndParams(Offsets *offsets);
void readRangeOrRangeOffsets( Offsets *offsets,int orbitType);
void readOffsets( Offsets *offsets);
void readAzimuthOffsets( Offsets *offsets);
void readRangeOffsets( Offsets *offsets);
void readOldPar(char *parFile, SARData *sarD, stateV *sv);
void geometryInfo(conversionDataStructure *cP,inputImageStructure *currentImage,double azimuth,double range, double z, double thetaC,double *ReH,double *Range,double *theta,double *thetaD, double *phi,double zSp);

double computeSig2AzParam(double sinTheta,double cosTheta, double azimuth, double Range, inputImageStructure *inputImage,Offsets *offsets);
double computeSig2Base(double sinThetaD,double cosThetaD, double azimuth, inputImageStructure *inputImage,Offsets *offsets);

void interpTideError(double *phaseError,  inputImageStructure *phaseImage,	vhParams *params,double x,double y, double dPsi, double twok);

double sphericalElev(double z,double lat,double Re);
double sphericalToWGSElev(double z,double lat,double Re);
void computeA(double lat,double lon,double x,double y,  inputImageStructure *aPhaseImage,  inputImageStructure *dPhaseImage,  double A[2][2]); 
void computeB(double x, double y, double z, double B[2][2], double *dzdx, double *dzdy,double aPsi, double dPsi,xyDEM *xydem);
void computeVxy(double aP, double dP,double aPe, double dPe,double A[2][2], double B[2][2],	 double *vx, double *vy, double *scaleX,double *scaleY);
void getMosaicInputImage(inputImageStructure *inputImage);
void errorsToXY(double er,double ea, double *ex, double *ey, double xyAngle, double hAngle);
void getIntersect(inputImageStructure *dPhaseImage,inputImageStructure *aPhaseImage,  int *iMin,int *iMax,int *jMin,int *jMax, outputImageStructure *outputImage);
void computeSceneAlpha(outputImageStructure *outputImage,	inputImageStructure *aOffImage, inputImageStructure *dOffImage,
		       conversionDataStructure *aCp,  conversionDataStructure *dCp, xyDEM *dem,  int *iMin,int *iMax,int *jMin,int *jMax);

/*ers1Complex **mallocComplexImage(int nr,int nc);*/

float **mallocImage(int nr,int nc);
void centerLLNew( inputImageStructure *inputImage, double *lat, double *lon,double deltaT);
long julday(int mm, int id, int iyyy);
double juldayDouble(int mm, int id, int iyyy);

double earthRadius(double lat, double rp, double re);
double earthRadiusCurvatureWGS84(double lat);
double earthRadiusWGS84(double lat);
void interpPhaseImage(inputImageStructure *inputImage, double range,double azimuth,  double *phase);
void polintVec(double xa[], double y1[],double y2[],double y3[],double y4[],double y5[],double y6[],
	       double x, double *yr1,double *yr2,double *yr3,double *yr4,double *yr5,double *yr6 );
double getHeight(double lat, double lon, demStructure *dem, double Re,int heightFlag);
double rhoRReZReH(double R,double ReZ, double ReH);
double thetaRReZReH(double R,double ReZ, double ReH);
double psiRReZReH(double R,double ReZ, double ReH);
double slantRange(double rho,double ReZ, double ReH);
void readShelf(outputImageStructure *outputImage,char *shelfMaskFile);
void initOutputImage(outputImageStructure *outputImage, inputImageStructure inputImage);
void getDEM(char *demFile, void *dem1);
void smlocateZD(double rxs,  double rys,  double rzs, double rvsx, double rvsy, double rvsz, double rsl,double *lat, double *lon, double lookDir, double trgalt);
void computeTiePoints(inputImageStructure *inputImage, tiePointsStructure *tiePoints,
                      demStructure dem,int noDEM,char *geodatFile,ShelfMask *shelfMask, int supressOutput);
void readTiePoints(FILE *fp, tiePointsStructure *tiePoints, int noDEM);
double dot(double x1,double y1,double z1,double x2, double y2, double z2);
void cross(double a1,double a2,double a3, double b1,double b2,double b3, double *c1,double *c2,double *c3);
double norm(double x1,double y1,double z1);
void svBaseTCN(double myTime, double dt1t2, stateV *sv1,stateV *sv2,double bTCN[3]);
void svBnBp(double myTime, double theta, double dt1t2, stateV *sv1,stateV *sv2,double *bn,double *bp);
void svOffsets(inputImageStructure *image1,inputImageStructure *image2, Offsets *offsets, double *cnstR, double *cnstA);
void svInitBnBp(inputImageStructure *inputImage, Offsets *offsets );
void svInterpBnBp(inputImageStructure *inputImage, Offsets *offsets, double azimuth,double *bnS,double *bpS );
void svTest(inputImageStructure *inputImage, Offsets *offsets );
void svInitAzParams(inputImageStructure *inputImage, Offsets *offsets );
double svAzOffset(inputImageStructure *inputImage, Offsets *offsets,double range, double azimuth );
void setTiePointsMapProjectionForHemisphere(tiePointsStructure *tiePoints);
double getReH( conversionDataStructure *cP, inputImageStructure *inputImage,double azimuth);
double bPoly(double b0,double b1,double b2,double x);
double secondForSAR(SARData *par);
char *appendSuffix(char *file, char *suffix,char *buf );
int rangeAzimuthToLL(double *rg, double range,double iFloat,double rhoSp,double ReH, double Re, double *lat, double *lon,double *hWGS,inputImageStructure *inputImage,xyDEM *xyDem,double tol,double step );
double groundRangeToLLNew(double groundRange, double azimuth, double *lat,double *lon, inputImageStructure *inputImage,int recycle);
void llToECEF(double lat,double lon,double h, double *x,double *y,double *z);
void getState(double myTime,inputImageStructure *inputImage, double *xs,double *ys, double *zs, double *vsx, double *vsy, double *vsz);
