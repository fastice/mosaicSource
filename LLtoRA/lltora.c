#include "stdio.h"
#include"string.h"
#include "source/common/common.h"
#include <sys/types.h>
#include <sys/time.h>
#include <math.h>
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "time.h"
/*
  Estimate r,a for given lat lon z 
  This program is uses some of the routines for geocode, 
  which means there is alot of unused junk to initialize everything correctly.
*/
static void  readArgs(int argc,char *argv[],int *passType, char **demFile,char **inputFile, char **tiePointFile, char **outFile  ) ;
static void usage();
void  computeScaleLS(float **inImage,float **scale,  int32 azimuthSize,int32 rangeSize,float fl,   float weight,double minVal,int32 iMin,int32 iMax,int32 jMin,int32 jMax);
void readLLinput(FILE *fp, tiePointsStructure *tiePoints, char *DEM) ;
/* 
   Global variables definitions, many are not used but are so they are defined in various functions
*/
int RangeSize=RANGESIZE;       /* Range size of complex image */
int AzimuthSize=AZIMUTHSIZE;   /* Azimuth size of complex image */
int BufferSize=BUFFERSIZE;     /* Size of nonoverlap region of the buffer */
int BufferLines = 512;         /* # of lines of nonoverlap in buffer */
double RangePixelSize=RANGEPIXELSIZE; /* Range PixelSize */
double AzimuthPixelSize=AZIMUTHPIXELSIZE; /* Azimuth PixelSize */
int HemiSphere=NORTH;
double Rotation=45.;
int sepAscDesc=TRUE;   

float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
char *Abuf1,*Abuf2,*Dbuf1,*Dbuf2;
int llConserveMem=999; /* Kluge to maintain backwards compat 9/13/06 */
float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
char *Abuf1,*Abuf2,*Dbuf1,*Dbuf2;
void *offBufSpace1,*offBufSpace2,*offBufSpace3,*offBufSpace4;
void *lBuf1,*lBuf2,*lBuf3,*lBuf4;

void main(int argc, char *argv[])
{   
	FILE *tiePointFp, *fpOut;
	demStructure dem;
	tiePointsStructure tiePoints;
	inputImageStructure inputImage;       /* input image */
	char   *inputFile, *tiePointFile  ;
	char *outputFile, *DemFile;
	int imageFlag,passType;
	int bufferSize;
	int i,j;                               /* LCV */
	float dr,da;
	ShelfMask *shelfMask=NULL;
	uint32 size[2];
	float dum[3];
	/* 
	   Read command line args and compute filenames
	*/
	readArgs(argc,argv, &passType,&DemFile,&inputFile,&tiePointFile,&outputFile);
	if(inputFile !=NULL) fprintf(stderr,"inputFile  %s \n",inputFile);
	if(DemFile !=NULL) fprintf(stderr,"DEM %s \n", DemFile);
	if(outputFile !=NULL) fprintf(stderr,"output %s \n",outputFile);
	/*	  Parse input file	*/
	inputImage.passType = passType;
	inputImage.stateFlag = TRUE;
	parseInputFile(inputFile, &inputImage);
	inputImage.imageType = PHASE;
	inputImage.useNew=TRUE;
	fprintf(stderr,"useNew %i\n",inputImage.useNew);
	/*	  Input tiepoints	*/
	tiePointFp = openInputFile(tiePointFile);
	if(outputFile == NULL) {
		readTiePoints(tiePointFp, &tiePoints,TRUE);
	} else		{
		readLLinput(tiePointFp, &tiePoints, DemFile) ;
	}
	if(tiePoints.lat[0] < 0) {
		fprintf(stderr,"**** SOUTHERN HEMISPHERE ****");
		HemiSphere=SOUTH; tiePoints.stdLat=71; Rotation=0.0;
	} else fprintf(stderr,"**** NORTHERN HEMISPHERE ****");
	tiePoints.imageCoords = FALSE;
	/*
	  Call to set up stuff, don't really use the outputimage
	*/ 
	computeTiePoints(&inputImage,&tiePoints,dem,TRUE,inputFile,shelfMask,TRUE);
	/*
	  Output results for checking to sterr
	*/
	dr=(inputImage.nRangeLooks-1)/2.0;
	da=(inputImage.nAzimuthLooks-1)/2.0;
	fprintf(stderr,"%f %f\n",dr,da);
	if(outputFile == NULL) {
		fprintf(stdout,"# 5\n");
		for(i=0; i < tiePoints.npts; i++) 
			fprintf(stdout," %10.3f %10.3f %10.7f %10.7f %7.1f\n",
				(tiePoints.r[i])*inputImage.nRangeLooks+dr,tiePoints.a[i]*inputImage.nAzimuthLooks+da, tiePoints.lat[i], tiePoints.lon[i],tiePoints.z[i]);   
		fprintf(stdout,"&\n");
	} else { /* Output results to binary output file */
		fpOut=fopen(outputFile,"w");
		size[0]=tiePoints.npts;
		size[1]=3;
	        fwriteBS( (void *)size,sizeof(uint32),2, fpOut,INT32FLAG);
		for(i=0; i < tiePoints.npts; i++)  {
			dum[0]=tiePoints.r[i]*inputImage.nRangeLooks;	dum[1]=tiePoints.a[i]*inputImage.nAzimuthLooks;	dum[2]=tiePoints.z[i];
			if(dum[0] < 0 || dum[0] > (inputImage.nRangeLooks*inputImage.rangeSize)-1) { dum[0]=-9999; dum[1]=-9999;}
			if(dum[1] < 0 || dum[1] > (inputImage.nAzimuthLooks*inputImage.azimuthSize-1)) {dum[0]=-9999; dum[1]=-9999;}			
		        fwriteBS( (void *)dum,sizeof(float),3, fpOut,FLOAT32FLAG);
		}
	}
	return; 
}

static void  readArgs(int argc,char *argv[],int *passType, char **demFile,char **inputFile, char **tiePointFile, char **outFile  ) 
{
	int filenameArg;
	char *argString;
	if( argc < 3 || argc > 5 ) usage();        /* Check number of args */ 
	*passType = DESCENDING;  

	if(*passType < 0) *passType = DESCENDING;
	*demFile=NULL;
	if(argc==5) {
		*outFile=argv[argc-1];
		*tiePointFile = argv[argc-2];
		*demFile = argv[argc-3];
		*inputFile = argv[argc-4];
	} else if(argc==4 ) {
		*outFile=argv[argc-1];
		*tiePointFile = argv[argc-2];
		*inputFile = argv[argc-3];
	}else {
		*outFile=NULL;
		*tiePointFile = argv[argc-1];
		*inputFile = argv[argc-2];
	}
	return;
}

 
static void usage()
{ 
	error("\n\n%s\n%s\n\n%s\n\n%s\n\n%s\n%s\n%s\n\n%s\n%s\n",
	      "Compute single look image pixels locations for lat lon z",  
	      "Output is to stdout range,azimuth,lat,lon,z ",
	      "Usage:",
	      " lltora    geoInputFile  llFile",
	      " lltora    geoInputFile  llFile outfile (binary input mode - 3 columns in)",
	      " lltora    geoInputFile  DEM llFile outfile (binary input mode - 2 columns in)",
	      "where",
	      "  geoInputFile = geodat file",
	      "  llFile       = (lat,lon,z)" );
}

/*
  Read binary file for tiepoints.
*/
void readLLinput(FILE *fp, tiePointsStructure *tiePoints, char *DEM)
{
	double lat, lon;
	double dum[3];
	uint32 size[2];
	double xx,yy;
	uint32 i;
	xyDEM  dem;
     
	freadBS( ( void *)size, sizeof(size[0]),2,fp,INT32FLAG);
	fprintf(stderr,"Size %i %i \n",size[0],size[1]);
	if(size[1] != 2 && size[1] != 3) error("Invalid size in first line");
	/* Modified 4/02/07 to fix crash with declared in structure */
	tiePoints->lat = (double *)malloc(sizeof(double)*size[0]);  tiePoints->lon = (double *)malloc(sizeof(double)*size[0]);
	tiePoints->x = (double *)malloc(sizeof(double)*size[0]);   tiePoints->y = (double *)malloc(sizeof(double)*size[0]);
	tiePoints->z = (double *)malloc(sizeof(double)*size[0]);
	tiePoints->r = (double *)malloc(sizeof(double)*size[0]);  tiePoints->a = (double *)malloc(sizeof(double)*size[0]);
	/* End 4/2/7 fix */
	if(size[1] == 2) 	readXYDEM(DEM,&dem);
	tiePoints->npts=0;
	for(i=0; i < size[0];  i++) {
		/* Read line */
		freadBS( ( void *)dum, sizeof(dum[0]),size[1],fp,FLOAT64FLAG);
		/* set hemisphere */
		if(i==0) {
			if(dum[0] < 0) {
				fprintf(stderr,"**** SOUTHERN HEMISPHERE ****");
				HemiSphere=SOUTH;	   tiePoints->stdLat=71;		Rotation=0.0;
			} else {
				fprintf(stderr,"**** NORTHERN HEMISPHERE ****");
				HemiSphere=NORTH;	   tiePoints->stdLat=70;		Rotation=45.;
			} 
		}
		if(size[1] == 2) {
			lltoxy1(dum[0],dum[1],&xx,&yy,Rotation,tiePoints->stdLat);
			dum[2]=interpXYDEM(xx, yy,dem);
			if( dum[2] < 0) {dum[0]=-1; dum[1]=-1;}
		}

		tiePoints->lat[tiePoints->npts] = dum[0];
		tiePoints->lon[tiePoints->npts] = dum[1];
		tiePoints->z[tiePoints->npts] = dum[2];    
		(tiePoints->npts)++;
	}
	if(tiePoints->npts > size[0]) /* Too many tiepts ? */
		error("readTiePoints -- size[0]=%i exceeded %i",   size[0], tiePoints->npts);
}

