#define MEM1 7777
#define MEM2 8888
#define   MAXADBUF 150000000
#define MAXADBUF2 5000000
#define MXST 75
#define MNST 5
/*
  This is the include file for using geocode routines with other programs.
*/

typedef struct shelfMaskType {
	double x0;
	double y0;
	int xSize;
	int ySize;
	double deltaX;
	double deltaY;
	double stdLat;
	int hemisphere;
	double rot;
	unsigned char **mask;
} ShelfMask;


typedef struct stateVType {
	int nState;
	double t0;
	double deltaT;
	double x[1000]; 	/* Note state vectors store as 1:nState */
	double y[1000];
	double z[1000];
	double vx[1000];
	double vy[1000];
	double vz[1000];
	double times[1000];
} stateV;

typedef struct conversionDataType {
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
	int azOff;
	int azSize;
	int rSize;
	int passType;
} conversionDataStructure;

typedef struct xyDEMTYPE {
	double x0;
	double y0;
	int xSize;
	int ySize;
	double deltaX;
	double deltaY;
	double stdLat;
	int hemisphere;
	double rot;
	float **z;
} xyDEM;

typedef struct xyVELTYPE {
	double x0;
	double y0;
	int xSize;
	int ySize;
	double deltaX;
	double deltaY;
	double stdLat;
	int hemisphere;
	double rot;
	float **vx;
	float **vy;
} xyVEL;

typedef struct sarDataType {
	char *label;
	double lambda;
	double prf;
	int hr;
	int min;
	double sec;
	int year;
	int month;
	int day;
	double slpR;
	double slpA;
	int nSlpR;  /* n single look range pixels */
	int nSlpA;  /* n single look azimuth pixels */
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

typedef struct inputImageType {
	char *file;
	int stateFlag;
	int llInit;
	int passType;
	int rangeSize;
	int azimuthSize;
	int nRangeLooks;
	int nAzimuthLooks;
	int imageType;
	int nAzGrid; /* Number of lat/lon regions along track used for geoocoding */
	int nRGrid;   /* Number of lat/lon region across track used for geododing */
	int nAzPixGrid;  /* Nominal azimuth size for region on grid */
	int nRPixGrid;   /* Nominal range size for region on grid */
	int lookDir;
	double azimuthPixelSize;
	double rangePixelSize;
	double *latControlPoints; /* Maintain for backwards compatability */
	double *lonControlPoints;
	double minLat;
	double maxLat;
	double minLon;
	double maxLon;	
	conversionDataStructure **conversionData;
	conversionDataStructure cpAll;
	SARData par;
/*	PROC_PAR par;         Structure with some of the parameters from the par file (those which are in the geodat */
	float weight;      /* Weight to apply to overall image */
	float **scale;      /* Scale array */
	double tideCorrection; /* Tidal correction for velocity estimation vz m/yr */
	int tideDiffFlag;
	xyDEM tideDiff;
	int removePad;  /* strip off first and last removePad columns */
	float *rAnt; /* range and pattern for antenna pattern correction */
	float *pAnt;
	float betaNought; /* beta nought conversion for S1 radiometric correction */
	int patSize;
	int year;
	int month;
	int day;
	double julDay;
	int orbit;
	int memChan; /* Kluge added 05/31/07 to remove asc/desc dependence for mem init */
	int used;  /* set to indicate whether data from this image was used */
	int crossFlag ;  /* Flag use to indicate don't use for crossing orbit solution */
	int useNew;
	stateV sv;
	void **image;
	int isInit;
	float noData;
	double lastTime;  /* Last time used for geocoding */
	float tolerance; /* Tolerange for geocoding. Set to 1e-6 in initLLtoImage */
	struct inputImageType *next;
} inputImageStructure;

typedef struct outputImageType {
	int xSize;
	int ySize;
	int noMem; /* Flag set to true when only output parameters needed */
	int imageType;
	double deltaX;
	double deltaY;
	double originX;
	double originY;
	double slat;    /* this has not been fully implemented yet */
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
	/*  flags for velocity work */
	int noVhFlag;  
	int no3d;
	int vzFlag;
	int rOffsetFlag;
	int noTide;
	int timeOverlapFlag;
	FILE *fpLog;
	char *date1;
	char *date2;
	double jd1;
	double jd2;
	int deltaB; /* used to indicate sv baselines */
} outputImageStructure;


typedef struct demType {
	int coordType;
	double minLat;
	double minLon;
	double maxLat;
	double maxLon;
	double deltaLat;
	double deltaLon;
	int latSize;
	int lonSize;
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
#define ROT  45.00
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
/*
  Input phase or power image for geocode
*/
void getGeoCodeInputImage(char *imageFile, inputImageStructure *inputImage,
			  int imageFlag);
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
/*
  Output geocode image. Writes two files one for image, and xxx.geodat 
  with image header info
*/
void outputGeocodedImage(outputImageStructure outputImage,char *outputFile);
/* 
   Convert lat/lon and h to image coordinates using alogrithm by Shusun Li.
*/

/*
  Compute height for lat/lon from DEM.
*/
double getHeight(double lat, double lon, demStructure *dem, double Re,int heightFlag);
/*
  Interpolate input image interpolateInputImage.c
*/

float interpolateFloatInputImage(inputImageStructure inputImage, double range,double azimuth);

/*
  Input xy (PS) DEM.
*/
void readXYDEMGeoInfo(char *xyFile, xyDEM *xydem, int resetProjection);
void readXYDEMcrop( char *xyFile,xyDEM *xydem,float xmin,float xmax,float ymin, float ymax);
void readXYDEM( char *xyFile,xyDEM *xydem);
void readXYVel(xyVEL *xyvel, char *xyFile);

void xytoll1(double x, double y,int hemi, double *alat,double *alon,
	     double dlam,double slat);
void xytoll(double x, double y,int hemi, double *alat,double *alon,
	    double dlam);
   
void lltoxy1(double alat,double alon,double *x, double *y, double dlam, double slat);
void lltoxy(double alat,double alon,double *x, double *y, double dlam);

/*    double getXYDEM(double lat, double lon, xyDEM *dem,
      double Re,int heightFlag);
*/

double getXYHeight(double lat, double lon, xyDEM *xydem,
		   double Re,int heightFlag);

/*
  Memory allocation routines
*/
float **mallocImage(int nr,int na);
/*ers1Complex **mallocComplexImage(int nr,int na);*/

/*---***    void findRange(double out_y, double h, int azimuth, double *R, 
  conversionDataStructure *cp );*/
void findRange(double out_y, double h, int azimuth, double *R, 
	       conversionDataStructure *cp,double ReAll );

void asfConversions(double lat, double lon,double *out_x,double *out_y,
		    conversionDataStructure *conversionData );

void initllToImageNew( inputImageStructure *inputImage);
void llToImageNew(double lat,double lon, double h, double *range, double *azimuth, inputImageStructure *inputImage);
