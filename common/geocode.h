#define MEM1 7777
#define MEM2 8888
#define MAXADBUF 150000000
#define MAXADBUF2 5000000
#define MXST 75
#define MNST 5

#define LSB 0
#define MSB 1
/*
  This is the include file for using geocode routines with other programs.
*/

typedef struct shelfMaskType
{
	double x0;
	double y0;
	int32_t xSize;
	int32_t ySize;
	double deltaX;
	double deltaY;
	double stdLat;
	int32_t hemisphere;
	double rot;
	unsigned char **mask;
} ShelfMask;

typedef struct stateVType
{
	int32_t nState;
	double t0;
	double deltaT;
	double x[1000]; /* Note state vectors store as 1:nState */
	double y[1000];
	double z[1000];
	double vx[1000];
	double vy[1000];
	double vz[1000];
	double times[1000];
} stateV;

typedef struct conversionDataType
{
	double cf_x[10];
	double cf_y[10];
	double a2[3][3];
	double pixelSize;
	double pixelToAzimuthPixel;
	double toRangePixel;
	double toAzimuthPixel;
	double RNear;
	double RCenter;
	double RFar;
	double *rNear;
	double *rFar;
	double Re;
	double swathWidth;
	double latControlPoints[6];
	double lonControlPoints[6];
	double *ReH;
	double sTime;
	double eTime;
	double prf;
	int32_t azOff;
	int32_t azSize;
	int32_t rSize;
	int32_t passType;
} conversionDataStructure;

typedef struct xyDEMTYPE
{
	double x0;
	double y0;
	int32_t xSize;
	int32_t ySize;
	double deltaX;
	double deltaY;
	double stdLat;
	int32_t hemisphere;
	double rot;
	float **z;
} xyDEM;

typedef struct xyVELTYPE
{
	double x0;
	double y0;
	int32_t xSize;
	int32_t ySize;
	double deltaX;
	double deltaY;
	double stdLat;
	int32_t hemisphere;
	double rot;
	float **vx;
	float **vy;
} xyVEL;

typedef struct sarDataType
{
	char *label;
	double lambda;
	double prf;
	int32_t hr;
	int32_t min;
	double sec;
	int32_t year;
	int32_t month;
	int32_t day;
	double slpR;
	double slpA;
	int32_t nSlpR; /* n single look range pixels */
	int32_t nSlpA; /* n single look azimuth pixels */
	double fd[4];
	double rawRange[3]; /* Not Used */
	double rn;
	double rc;
	double rf;
	double H;
	double ReMajor;
	double ReMinor;
	double rangeError;
	double lookDir;
	double echoTD; /* not currently used */
} SARData;

typedef struct inputImageType
{
	char *file;
	int32_t stateFlag;
	int32_t llInit;
	int32_t passType;
	int32_t rangeSize;
	int32_t azimuthSize;
	int32_t nRangeLooks;
	int32_t nAzimuthLooks;
	int32_t imageType;
	int32_t nAzGrid;	/* Number of lat/lon regions along track used for geoocoding */
	int32_t nRGrid;		/* Number of lat/lon region across track used for geododing */
	int32_t nAzPixGrid; /* Nominal azimuth size for region on grid */
	int32_t nRPixGrid;	/* Nominal range size for region on grid */
	int32_t lookDir;
	double azimuthPixelSize;
	double rangePixelSize;
	double *latControlPoints; /* Maintain for backwards compatability */
	double *lonControlPoints;
	double minLat;
	double maxLat;
	double minLon;
	double maxLon;
	double minX;
	double maxX;
	double minY;
	double maxY;
	conversionDataStructure **conversionData;
	conversionDataStructure cpAll;
	SARData par;
	/*	PROC_PAR par;         Structure with some of the parameters from the par file (those which are in the geodat */
	float weight;		   /* Weight to apply to overall image */
	float **scale;		   /* Scale array */
	double tideCorrection; /* Tidal correction for velocity estimation vz m/yr */
	int32_t tideDiffFlag;
	xyDEM tideDiff;
	int32_t removePad; /* strip off first and last removePad columns */
	float *rAnt;	   /* range and pattern for antenna pattern correction */
	float *pAnt;
	float betaNought; /* beta nought conversion for S1 radiometric correction */
	int32_t patSize;
	int32_t year;
	int32_t month;
	int32_t day;
	double julDay;
	int32_t orbit;
	int32_t memChan;   /* Kluge added 05/31/07 to remove asc/desc dependence for mem init */
	int32_t used;	   /* set to indicate whether data from this image was used */
	int32_t crossFlag; /* Flag use to indicate don't use for crossing orbit solution */
	int32_t useNew;
	stateV sv;
	void **image;
	int32_t isInit;
	float noData;
	double lastTime; /* Last time used for geocoding */
	double tolerance; /* Tolerange for geocoding. Set to 1e-6 in initLLtoImage */
	struct inputImageType *next;
} inputImageStructure;

typedef struct outputImageType
{
	int32_t xSize;
	int32_t ySize;
	int32_t noMem; /* Flag set to true when only output parameters needed */
	int32_t imageType;
	double deltaX;
	double deltaY;
	double originX;
	double originY;
	double slat; /* this has not been fully implemented yet */
	int32_t makeTies;
	void **image;
	void **image2; /* Second Image for special applications (i.e vx,vy) */
	void **image3; /* Third image for 3-d velocities */
	float **scale; /* Scale array for mosaicing */
	float **scale2;
	float **scale3;
	/*  Arrays for carrying around temp space in velocity mosaicking */
	float **vxTmp;
	float **vyTmp;
	float **vzTmp;
	float **fScale;
	float **sxTmp;
	float **syTmp;
	/* Error arrays */
	float **errorX;
	float **errorY;
	float **exTmp;
	float **eyTmp;
	/*
	  Mask for ice shelf
	*/
	ShelfMask *shelfMask;
	xyDEM *verticalCorrection;
	/*  flags for velocity work */
	int32_t noVhFlag;
	int32_t no3d;
	int32_t vzFlag;
	int32_t rOffsetFlag;
	int32_t noTide;
	int32_t timeOverlapFlag;
	FILE *fpLog;
	char *date1;
	char *date2;
	double jd1;
	double jd2;
	int32_t deltaB; /* used to indicate sv baselines */
} outputImageStructure;

typedef struct demType
{
	int32_t coordType;
	double minLat;
	double minLon;
	double maxLat;
	double maxLon;
	double deltaLat;
	double deltaLon;
	int32_t latSize;
	int32_t lonSize;
	float **demData;
} demStructure;

#define NCONTROLPOINTS 6
#define NORTH 0
#define EAST 1
#define SOUTH 2
#define WEST 3

#define ELLIPSOIDAL 0
#define SPHERICAL 1
#define RIGHT 1.0
#define LEFT -1.0
#define ROT 45.00
#define HEMISPHERE NORTH
#define MTOKM 0.001
#define KMTOM 1000.0
#define NINCIDENCEANGLES 1024
#define DELTAGROUNDRANGE 10.
#define NEXTRAANGLES 100
#define XY 0
#define LATLON 1
/*
  Parse input file for geocode.
*/
void parseInputFile(char *inputFile, inputImageStructure *inputImage);
void parseGeojson(char *inputFile, inputImageStructure *inputImage);
// Debug routine
void dumpParsedInput(inputImageStructure *inputImage, FILE *fp);
/*
  Input phase or power image for geocode
*/
void getGeoCodeInputImage(char *imageFile, inputImageStructure *inputImage,
						  int32_t imageFlag);
/*
  Input dem
*/
void getDEM(char *demFile, void *dem);
/*
  Init output image
*/
void initOutputImage(outputImageStructure *outputImage,
					 inputImageStructure inputImage);
/*
  Geocode image using  DEM for terrain correction.
*/
void geoCodeImage(inputImageStructure inputImage,
				  outputImageStructure outputImage,
				  void *dem);
// Write LSB or MSB data
size_t fwriteOptionalBS(void *ptr, size_t nitems, size_t size, FILE *fp, int32_t flags, int32_t byteOrder);
/*
  Output geocode image. Writes two files one for image, and xxx.geodat
  with image header info
*/
void outputGeocodedImage(outputImageStructure outputImage, char *outputFile);
/*
   Convert lat/lon and h to image coordinates using alogrithm by Shusun Li.
*/

/*
  Compute height for lat/lon from DEM.
*/
double getHeight(double lat, double lon, demStructure *dem, double Re, int32_t heightFlag);
/*
  Interpolate input image interpolateInputImage.c
*/

float interpolateFloatInputImage(inputImageStructure inputImage, double range, double azimuth);

/*
  Input xy (PS) DEM.
*/
void readXYDEMGeoInfo(char *xyFile, xyDEM *xydem, int32_t resetProjection);
void readXYDEMcrop(char *xyFile, xyDEM *xydem, float xmin, float xmax, float ymin, float ymax);
void readXYDEM(char *xyFile, xyDEM *xydem);
void readXYVel(xyVEL *xyvel, char *xyFile);

void xytoll1(double x, double y, int32_t hemi, double *alat, double *alon,
			 double dlam, double slat);
void xytoll(double x, double y, int32_t hemi, double *alat, double *alon,
			double dlam);

void lltoxy1(double alat, double alon, double *x, double *y, double dlam, double slat);
void lltoxy(double alat, double alon, double *x, double *y, double dlam);

/*    double getXYDEM(double lat, double lon, xyDEM *dem,
	  double Re,int32_t heightFlag);
*/

double getXYHeight(double lat, double lon, xyDEM *xydem,
				   double Re, int32_t heightFlag);

/*
  Memory allocation routines
*/
float **mallocImage(int32_t nr, int32_t na);
/*ers1Complex **mallocComplexImage(int32_t nr,int32_t na);*/

/*---***    void findRange(double out_y, double h, int32_t azimuth, double *R,
  conversionDataStructure *cp );*/
void findRange(double out_y, double h, int32_t azimuth, double *R,
			   conversionDataStructure *cp, double ReAll);

void asfConversions(double lat, double lon, double *out_x, double *out_y,
					conversionDataStructure *conversionData);

void initllToImageNew(inputImageStructure *inputImage);
void llToImageNew(double lat, double lon, double h, double *range, double *azimuth, inputImageStructure *inputImage);
