#include "stdio.h"
#include "stdlib.h"
#include"string.h"
#include "mosaicSource/common/common.h"
#include "rparams.h"
#include <sys/types.h>
#include <sys/time.h>


/*
  Estimate range params using tiepoints.

  This program uses some of the routines for geocode, 
  which means there is alot of unused junk to initialize everything correctly.
*/

static void  readArgs(int argc,char *argv[],char **geodatFile,    char **tiePointFile,char **offsetFile, char **baselineFile,tiePointsStructure *tiepoints, char **shelfMaskFile);
static void setMapProjectionForHemisphere(tiePointsStructure *tiePoints);

char *Abuf1,*Abuf2,*Dbuf1,*Dbuf2;
int llConserveMem=999; /* NO mem conserve Kluge to maintain backwards compat 9/13/06 */
void *offBufSpace1,*offBufSpace2,*offBufSpace3,*offBufSpace4;
void *lBuf1,*lBuf2,*lBuf3,*lBuf4;
static void usage();
#define MAXOFFBUF 72000000
#define MAXOFFLENGTH 30000
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
int sepAscDesc=TRUE;
float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 not use only for mosaic3d compatability */


void main(int argc, char *argv[])
{   
	extern char *Abuf1,*Abuf2,*Dbuf1,*Dbuf2;
	FILE *tiePointFp;
	demStructure dem;
	tiePointsStructure tiePoints;
	double nDays,stdLat;
	inputImageStructure inputImage, inputImage2;       /* input image */
	outputImageStructure outputImage;
	char *demFile, *geodatFile, *tiePointFile, *offsetFile, *baselineFile;
	double deltaR,deltaA;
	char *outputFile;
	char *shelfMaskFile;
	Offsets offsets;
	int imageFlag,passType,noDEM,noRamp,dBpFlag,timeReverseFlag;
	int bufferSize;
	int imageCoords;
	int linFlag;
	int i,j;                               /* LCV */
	Abuf1=NULL; Abuf2=NULL; Dbuf1=NULL;Dbuf2=NULL;
	/* Used for pointers rows to above buffer space */
	lBuf3=(void *)malloc(sizeof(float *)*MAXOFFLENGTH);
	lBuf4=(void *)malloc(sizeof(float *)*MAXOFFLENGTH);
	offBufSpace3=(void *)malloc(MAXOFFBUF);
	offBufSpace4=(void *)malloc(MAXOFFBUF);	 	
	/* 
	   Read command line args and compute filenames
	*/
	readArgs(argc,argv,&geodatFile,&tiePointFile,&offsetFile,&baselineFile,&tiePoints,&shelfMaskFile);
	/*
	  Parse input file
	*/
	noDEM=TRUE;
	offsets.rFile=offsetFile;
	parseInputFile(geodatFile, &inputImage); 
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
	/*
	  Call to set up stuff, don't really use the outputimage
	*/
	outputImage.noMem = TRUE;
	outputImage.originX=LARGEINT;
	outputImage.originY=LARGEINT;
	initOutputImage(&outputImage,inputImage);
	/*
	  Read shelf mask if needed
	*/
	outputImage.fpLog=stderr;
	if(shelfMaskFile != NULL) readShelf(&outputImage,shelfMaskFile);  else outputImage.shelfMask=NULL;
	/*
	  Compute image coords and get z from dem if necessary
	*/
	computeTiePoints(&inputImage,&tiePoints,dem,noDEM,geodatFile,outputImage.shelfMask,FALSE);
	/*
	  Extract offsets from phase file.
	*/
	fprintf(stderr,"%s\n",offsets.rFile);
	getROffsets(offsetFile,&tiePoints,inputImage,&offsets);
	fprintf(stderr,"%s %s\n",offsets.geo1,offsets.geo2);
	if(offsets.geo1 != NULL && offsets.geo2 != NULL && tiePoints.deltaB != DELTABNONE) {
		parseInputFile(offsets.geo2, &inputImage2);
		initllToImageNew( &inputImage2);
		memcpy(&(offsets.sv2),&(inputImage2.sv),sizeof(inputImage2.sv));
		offsets.dt1t2=inputImage.cpAll.sTime-inputImage2.cpAll.sTime;
		fprintf(stderr,"times %f %f %f\n",inputImage.cpAll.sTime,inputImage2.cpAll.sTime,offsets.dt1t2);
		svOffsets(&inputImage,&inputImage2,  &offsets,&(tiePoints.cnstR),&(tiePoints.cnstA));
	} else if( tiePoints.deltaB == DELTABNONE ) { tiePoints.cnstA=0.0; tiePoints.cnstR=0.0; }
	else error("SV baselines but geodats not specified in .dat file ");
	
	fprintf(stderr,"Rg/Az offsets %10.5f %10.5f\n",tiePoints.cnstR,tiePoints.cnstA);
	/*
	  Get baseline info
	*/
	getBaselineFile(baselineFile,&tiePoints,inputImage);
	/*
	 remove velocity components
	*/
	addVelCorrections(&inputImage,&tiePoints);
	/*
	  Output results for checking to sterr
	*/
	/* for(i=0; i < tiePoints.npts; i++) 
	   if( fabs(tiePoints.phase[i]) < 200000)  
	   fprintf(stderr,"%8.1f %8.1f %8.1f ---  %7.2f %7.2f --- %f ---- %f %f %f\n",
	   tiePoints.x[i],
	   tiePoints.y[i],tiePoints.z[i],tiePoints.r[i],tiePoints.a[i],
	   tiePoints.phase[i],tiePoints.vyra[i],tiePoints.vx[i],tiePoints.vy[i]);
	   fprintf(stderr,"\n");*/

	/*
	  Estimate baseline solution and output to stdout
	*/
	computeRParams(&tiePoints,inputImage,baselineFile,&offsets);

	return; 
}



static void  readArgs(int argc,char *argv[],char **geodatFile,    char **tiePointFile,char **offsetFile, char **baselineFile,tiePointsStructure *tiePoints, char **shelfMaskFile)

{
	int filenameArg;
	char *argString;
        double nDays;
        int bnbpFlag, bpFlag, bnbpdBpFlag, bpdBpFlag,constOnlyFlag,quadB,deltaB;
	int i,n;

	if( argc < 5 || argc > 22 ) usage();        /* Check number of args */ 
	n = argc - 5; 
	nDays = 24;
	quadB=FALSE;
	bpFlag = FALSE;
	bnbpFlag = FALSE;
	bnbpdBpFlag = FALSE;
	bpdBpFlag = FALSE;
	constOnlyFlag=FALSE;
	deltaB=DELTABNONE;
	shelfMaskFile=NULL;
	tiePoints->quiet=FALSE;
	for(i=1; i <= n; i++) {
		argString = strchr(argv[i],'-');  
		if(strstr(argString,"nDays") != NULL) {
			sscanf(argv[i+1],"%lf",&nDays); i++;
		} else if(strstr(argString,"shelfMask") != NULL) {
			*shelfMaskFile =  argv[i+1]; i++;
		} else if(strstr(argString,"bnbpOnly") != NULL ) {
			bnbpFlag = TRUE;
			if(bpFlag==TRUE || bpdBpFlag==TRUE || bnbpdBpFlag == TRUE) 
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -"
				      "incompatable");
		} else if(strstr(argString,"bnbpdBpOnly") != NULL ) {
			bnbpdBpFlag = TRUE;
			if(bpFlag==TRUE || bpdBpFlag==TRUE || bnbpFlag == TRUE) 
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -incompatable");
		}  else if(strstr(argString,"bpdBpOnly") != NULL ) {
			bpdBpFlag = TRUE;
			if(bpFlag==TRUE || bnbpFlag==TRUE || bnbpdBpFlag == TRUE) 
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB - incompatable");
		}  else if(strstr(argString,"quadB") != NULL ) {
			quadB = TRUE;
			if(bpFlag==TRUE || bnbpFlag==TRUE || bnbpdBpFlag == TRUE) 
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB - incompatable");
		}  else if(strstr(argString,"constOnly") != NULL ) {
			constOnlyFlag=TRUE;
			if(bpFlag==TRUE || bnbpFlag==TRUE || bnbpdBpFlag == TRUE) 
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB - incompatable");
		}  else if(strstr(argString,"deltaBQ") != NULL ) {
			deltaB=DELTABQUAD;
		}  else if(strstr(argString,"deltaBC") != NULL ) {
			deltaB=DELTABCONST;			
		}  else if(strstr(argString,"quiet") != NULL ) {
			tiePoints->quiet=TRUE;
		}  else if(i != n) usage();  
	}
	*geodatFile = argv[argc-4];
	*tiePointFile = argv[argc-3];
	*offsetFile = argv[argc-2];
	*baselineFile = argv[argc-1];
	/*
	  Set inputs
	*/
	tiePoints->constOnlyFlag = constOnlyFlag;
	tiePoints->linFlag = TRUE;
	tiePoints->quadB = quadB;
	tiePoints->dBpFlag=TRUE;
	tiePoints->bnbpFlag = bnbpFlag;
	tiePoints->bpFlag = bpFlag;
	tiePoints->bnbpdBpFlag = bnbpdBpFlag;
	tiePoints->bpdBpFlag = bpdBpFlag;
	tiePoints->nDays = nDays; 
	tiePoints->vrFlag = FALSE;
	tiePoints->deltaB = deltaB;
	if(tiePoints->constOnlyFlag == TRUE) fprintf(stderr,"\n(****Constant only fit*****\n");
	if( tiePoints->linFlag == TRUE) fprintf(stderr,"\n(****Including linear term fit*****\n");	
	return;
}

 
static void usage()
{ 
	error(
	      "\n\n%s\n%s\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
	      "Compute parameters to calibrate azimuth offsets",
	      "Usage:",
	      " rparams  -shelfMaskFile shelfMaskFile -deltaBQ -deltaBC -quadB -bpdBpOnly -bnbpdBpOnly -nDays nDays geodatFile tiepointsFile offsetFile "
	      "baselineFile",
	      "where",
	      "  shelfMask   	 = shelfMask file to indicate use tidal corrections",
	      "  constOnly   	 = estimate only the constant term",
	      "  deltaBQ   	 = estimate the quadratic correction to the state vector baseline",
	      "  deltaBC   	 = estimate the correction to Bp component of the state vector baseline",	      	      
	      "  quadB    	 = estimate  the quadratic terms",
	      "  bnbpdBpOnly    = estimate only bn,bp,dBp",
	      "  quiet          = don't echo tiepoints to solution",
	      "  nDays        	= temporals basline in days (default=24)",
	      "  geoDatFile  	= geo param file",
	      "  tiepointFile 	= Tiepoint location file in (lat,lon,z,vx,vy,vz)",
	      "  offsetFile  	= azimuth offset file (offsetFile.dat must also exist)",
	      "  baselineFile 	= cw state vector baseline file ");
}


