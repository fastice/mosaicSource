#include "stdio.h"
#include"string.h"
#include"stdlib.h"
#include "mosaicSource/common/common.h"
#include "azparams.h"
#include <sys/types.h>
#include <sys/time.h>
#include "math.h"
/*
  Estimate baseline using tiepoints.

  This program uses some of the routines for geocode, 
  which means there is alot of unused junk to initialize everything correctly.
*/

static void  readArgs(int argc,char *argv[],char **geodatFile, char **tiePointFile,char **offsetFile, char **baselineFile, tiePointsStructure *tiePoints);
static void usage();  

char *Abuf1,*Abuf2,*Dbuf1,*Dbuf2;
int llConserveMem=999; /* NO mem conserve Kluge to maintain backwards compat 9/13/06 */
#define MAXOFFBUF 72000000
#define MAXOFFLENGTH 30000
void *offBufSpace1,*offBufSpace2,*offBufSpace3,*offBufSpace4;
void *lBuf1,*lBuf2,*lBuf3,*lBuf4;
int sepAscDesc;
/* 
   Global variables definitions
*/
int RangeSize=RANGESIZE;       /* Range size of complex image */
int AzimuthSize=AZIMUTHSIZE;   /* Azimuth size of complex image */
int BufferSize=BUFFERSIZE;     /* Size of nonoverlap region of the buffer */
int BufferLines = 512;         /* # of lines of nonoverlap in buffer */
double RangePixelSize=RANGEPIXELSIZE; /* Range PixelSize */
double AzimuthPixelSize=AZIMUTHPIXELSIZE; /* Azimuth PixelSize */
int HemiSphere=NORTH;
double Rotation=45.;
double SLat=-91.0;
float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 not use only for mosaic3d compatability */


void main(int argc, char *argv[])
{   
	FILE *tiePointFp;
	demStructure dem;
	tiePointsStructure tiePoints;
	double nDays,stdLat;
	inputImageStructure inputImage;       /* input image */
	outputImageStructure outputImage;
	char *demFile, *geodatFile, *tiePointFile, *offsetFile, *baselineFile;
	char *outputFile;
	Offsets offsets;
	int imageFlag,passType,noDEM,noRamp,dBpFlag,timeReverseFlag;
	int bufferSize;
	int imageCoords;
	int constOnlyFlag, linFlag;
	int i,j;                               /* LCV */
	Abuf1=NULL; Abuf2=NULL; Dbuf1=NULL;Dbuf2=NULL;
	/* Used for pointers rows to above buffer space */
	lBuf1=(void *)malloc(sizeof(float *)*MAXOFFLENGTH);
	lBuf2=(void *)malloc(sizeof(float *)*MAXOFFLENGTH);
	offBufSpace1=(void *)malloc(MAXOFFBUF);
	offBufSpace2=(void *)malloc(MAXOFFBUF);	 
	/* 
	   Read command line args and compute filenames
	*/
	readArgs(argc,argv,&geodatFile,&tiePointFile,&offsetFile,&baselineFile, &tiePoints);
	/*
	  Parse input file
	*/
	inputImage.stateFlag = TRUE;
	noDEM=TRUE;
	offsets.file=offsetFile;
	parseInputFile(geodatFile, &inputImage); 
	/*
	  Set inputs
	*/
	if(constOnlyFlag == TRUE) fprintf(stderr,"\n(****Constant only fit*****\n");
	if(linFlag == TRUE) fprintf(stderr,"\n(****Including linear term fit*****\n");
	/*
	  Input tiepoints
	*/
	tiePointFp = openInputFile(tiePointFile);
	tiePoints.motionFlag=TRUE;
	readTiePoints(tiePointFp, &tiePoints,noDEM);
	/*
	  Determine hemisphere
	*/
	setTiePointsMapProjectionForHemisphere(&tiePoints);	
	tiePoints.imageCoords = FALSE;
	/*
	  Call to set up stuff, don't really use the outputimage
	*/
	outputImage.noMem = TRUE;
	outputImage.originX=LARGEINT;     outputImage.originY=LARGEINT;
	initOutputImage(&outputImage,inputImage);
	/*
	  Compute image coords and get z from dem if necessary
	*/
	computeTiePoints(&inputImage,&tiePoints,dem,noDEM,geodatFile,NULL,FALSE);
	/*
	  Extract phases from phase file.
	*/
	getOffsets(offsetFile,&tiePoints,inputImage,&offsets);
	fprintf(stderr,"%s %s %i\n",offsets.geo1,offsets.geo2,(int)tiePoints.deltaB);
	if(offsets.geo1 != NULL && offsets.geo2 != NULL && tiePoints.deltaB != DELTABNONE) {
		svInitAzParams(&inputImage,&offsets );
		fprintf(stderr,"sv fit %f %f %f %f\n",offsets.azFit[0],offsets.azFit[1],offsets.azFit[2],offsets.azFit[3]);
	} else { tiePoints.cnstA=0.0; tiePoints.cnstR=0.0; }
	/*
	  Add baseline corrections for previously removed baselines
	*/
	addOffsetCorrections( &inputImage,&tiePoints);
	/*
	  Output results for checking to sterr
	*/
	if(tiePoints.quiet == FALSE ) {
		for(i=0; i < tiePoints.npts; i++) 
			if( fabs(tiePoints.phase[i]) < 200000)  
				fprintf(stderr,"%8.1f %8.1f %8.1f ---  %7.2f %7.2f --- %f ---- %f\n",
					tiePoints.x[i],	tiePoints.y[i],tiePoints.z[i],tiePoints.r[i],tiePoints.a[i],tiePoints.phase[i],tiePoints.vyra[i]);
		fprintf(stderr,"\n");
	}
	/*
	  Estimate baseline solution and output to stdout
	*/
	computeAzParams(&tiePoints,&inputImage,baselineFile,&offsets); 
	return; 
}





static void  readArgs(int argc,char *argv[],char **geodatFile, char **tiePointFile,char **offsetFile, char **baselineFile, tiePointsStructure *tiePoints)
{
	int filenameArg;
	char *argString;
	double nDays=-1;
	int linFlag=FALSE,constOnlyFlag=FALSE;
	int deltaB=DELTABNONE;
	int i,n;
	if( argc < 5 || argc > 18 ) usage();        /* Check number of args */ 
	n = argc - 5; 
	tiePoints->quiet=FALSE;
	for(i=1; i <= n; i++) {
		argString = strchr(argv[i],'-');  
		if(strstr(argString,"constOnly") != NULL ) {
			constOnlyFlag = TRUE;
		} else if(strstr(argString,"linear") != NULL ) {
			linFlag = TRUE;
		} else if(strstr(argString,"nDays") != NULL) {
			sscanf(argv[i+1],"%lf",&nDays); i++;
		}  else if(strstr(argString,"useSV") != NULL ) {
			deltaB=DELTABCONST;
		}  else if(strstr(argString,"quiet") != NULL ) {
			tiePoints->quiet=TRUE;			
		} else if(i != n) usage();  
	}
	if(nDays < 0) error("nDays not specified in azparams");
	*geodatFile = argv[argc-4];
	*tiePointFile = argv[argc-3];
	*offsetFile = argv[argc-2];
	*baselineFile = argv[argc-1];
	tiePoints->deltaB=deltaB;
	tiePoints->linFlag=linFlag;
	tiePoints->constOnlyFlag=constOnlyFlag;
	tiePoints->nDays = nDays; 	
return;
}

 
static void usage()
{ 
	error(
	      "\n\n%s\n%s\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
	      "Compute parameters to calibrate azimuth offsets",
	      "Usage:",
	      " azparams  -constOnly -linear -useSV -quiet -nDays nDays geodatFile tiepointsFile offsetFile "
	      "baselineFile",
	      "where",
	      "  constOnly    = do not estimate baseline dependent terms, can be combined with linear fit",
	      "  linear       = add linear along track fit to either constOnly or the baseline parameter solution",	      
	      "  useSV     = as any of the four other options except it adds a correction after an SV determined offset removed",
	      "  quiet          = don't echo tiepoints to solution",
	      "  nDays        = temporals basline in days (default=24)",
	      "  geoDatFile   = geo param file",
	      "  tiepointFile = Tiepoint location file in (lat,lon,z,vx,vy,vz)",
	      "  offsetFile  =  azimuth offset file (offsetFile.dat must also exist)",
	      "  baselineFile = cw state vector baseline file ");
}
