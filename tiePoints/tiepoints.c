#include "stdio.h"
#include"string.h"
#include "mosaicSource/common/common.h"
#include "tiePoints.h"
#include <sys/types.h>
#include <sys/time.h>
#include <math.h>
int sepAscDesc;
/*
  Estimate baseline using tiepoints.

  This program uses some of the routines for geocode, 
  which means there is alot of unused junk to initialize everything correctly.
*/

static void  readArgs(int argc,char *argv[],int *imageFlag,int *passType,
		      int *noDEM,int *noRamp,char **demFile,char **inputFile,
		      char **tiePointFile,char **phaseFile,
		      char **baselineFile,int *dBpFlag,
		      int *imageCoords, int *motionFlag,int *timeReverseFlag,
		      double *nDays,int *quadB, double *stdLat, 
		      int *bnbpFlag, int *bpFlag, int *bnbpdBpFlag,
		      int *bpdBpFlag, int *vrFlag, char **shelfMaskFile);

static void usage();

char *Abuf1,*Abuf2,*Dbuf1,*Dbuf2;
int llConserveMem=999;  /* mem conserve Kluge to maintain backwards compat 9/13/06 */
/* 
   Global variables definitions
*/
int RangeSize=RANGESIZE;       /* Range size of complex image */
int AzimuthSize=AZIMUTHSIZE;   /* Azimuth size of complex image */
int BufferSize=BUFFERSIZE;     /* Size of nonoverlap region of the buffer */
int BufferLines = 512;         /* # of lines of nonoverlap in buffer */
double RangePixelSize=RANGEPIXELSIZE; /* Range PixelSize */
double AzimuthPixelSize=AZIMUTHPIXELSIZE; /* Azimuth PixelSize */
char *shelfMaskFile;
int HemiSphere=NORTH;
double Rotation=45.;
double SLat=-91.0;    
float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 not use only for mosaic3d compatability */

void *offBufSpace1,*offBufSpace2,*offBufSpace3,*offBufSpace4;
void *lBuf1,*lBuf2,*lBuf3,*lBuf4;


void main(int argc, char *argv[])
{   
	FILE *tiePointFp;
	demStructure dem;
	tiePointsStructure tiePoints;
	double nDays,stdLat;
	inputImageStructure inputImage;       /* input image */
	outputImageStructure outputImage;
	char *demFile, *inputFile, *tiePointFile, *phaseFile, *baselineFile;
	char *outputFile;
	int imageFlag,passType,noDEM,noRamp,dBpFlag,timeReverseFlag;
	int bufferSize;
	int imageCoords;
	int motionFlag, quadB,bnbpFlag,bpFlag,bnbpdBpFlag,bpdBpFlag,vrFlag;
	int i,j;                               /* LCV */
	/* 
	   Read command line args and compute filenames
	*/
	fprintf(stderr,"1\n");
	readArgs(argc,argv,&imageFlag,&passType,&noDEM,&noRamp,&demFile,
		 &inputFile,&tiePointFile,&phaseFile,&baselineFile,
		 &dBpFlag,&imageCoords,&motionFlag,&timeReverseFlag,&nDays,&quadB,
		 &stdLat,&bnbpFlag, &bpFlag, &bnbpdBpFlag,&bpdBpFlag, &vrFlag, &shelfMaskFile);
	/*
	  Parse input file
	*/

	inputImage.imageType = imageFlag;
	inputImage.passType = passType;
	inputImage.stateFlag = TRUE;
	parseInputFile(inputFile, &inputImage);

	if(passType == ASCENDING) fprintf(stderr,"Ascending pass\n");
	else fprintf(stderr,"Descending pass\n");

	tiePoints.motionFlag = motionFlag;
	tiePoints.timeReverseFlag = timeReverseFlag;
	tiePoints.nDays = nDays;
	tiePoints.quadB = quadB;
	tiePoints.bnbpFlag = bnbpFlag;
	tiePoints.bpFlag = bpFlag;
	tiePoints.bnbpdBpFlag = bnbpdBpFlag;
	tiePoints.bpdBpFlag = bpdBpFlag;
	tiePoints.stdLat = stdLat;
	tiePoints.vrFlag = vrFlag;

	if(quadB == TRUE) fprintf(stderr,"\n(****Quadratic fit*****\n");
	else if(bnbpFlag==TRUE) fprintf(stderr,"\n(****bn,bp fit*****\n");
	else if(bpFlag==TRUE) fprintf(stderr,"\n(****bp only fit*****\n");
	else if(bnbpdBpFlag==TRUE) fprintf(stderr,"\n(***bn,bp,dBp only fit****\n");
	else if(bpdBpFlag==TRUE) fprintf(stderr,"\n(****bp,dBp only fit*****\n");
	else fprintf(stderr,"\n(****Linear fit****\n*");

	if(motionFlag == TRUE) fprintf(stderr, "\n\n***Using %f day motion tiepoints*** \n\n",tiePoints.nDays);
	if(timeReverseFlag == TRUE) fprintf(stderr, "\n\n***Using timeReverse Flag *** \n\n");
	/*
	  Input tiepoints
	*/
	tiePointFp = openInputFile(tiePointFile);
	readTiePoints(tiePointFp, &tiePoints,noDEM);

	if(tiePoints.lat[0] < 0) {
		fprintf(stderr,"**** SOUTHERN HEMISPHERE ****");
		HemiSphere=SOUTH;
		tiePoints.stdLat=71;
		Rotation=0.0;
	} else fprintf(stderr,"**** NORTHERN HEMISPHERE ****");

	tiePoints.noRamp = noRamp;
	tiePoints.dBpFlag = dBpFlag;
	tiePoints.imageCoords = imageCoords;
	if(tiePoints.noRamp == TRUE) fprintf(stderr,"NO RAMP\n");
	else fprintf(stderr,"RAMP\n");
	/*
	  Call to set up stuff, don't really use the outputimage
	*/
	outputImage.noMem = TRUE;
	initOutputImage(&outputImage,inputImage);
	outputImage.fpLog=stderr;
	/*
	  Read shelf mask if needed
	*/
	if(shelfMaskFile != NULL) {
		readShelf(&outputImage,shelfMaskFile);
	} else outputImage.shelfMask=NULL;
	/*
	  Compute image coords and get z from dem if necessary
	*/
	computeTiePoints(&inputImage,&tiePoints,dem,noDEM,inputFile,outputImage.shelfMask,FALSE);
	/*
	  Extract phases from phase file.
	*/
	getPhases(phaseFile,&tiePoints,inputImage);
	/*
	  Add baseline corrections for previously removed baselines
	*/
	addBaselineCorrections( baselineFile, &tiePoints,inputImage);
	/*
	  Motion corrections
	*/
	if(motionFlag == TRUE || vrFlag==TRUE) {
		fprintf(stderr, "Before motion corrections\n");
		addMotionCorrections(inputImage,&tiePoints);
		fprintf(stderr, "After motion corrections\n");		
	 }
	/*
	  Output results for checking to sterr
	*/
	for(i=0; i < tiePoints.npts; i++) 
		if( fabs(tiePoints.phase[i]) < 200000)  fprintf(stderr,"%8.1f %8.1f %8.1f ---  %7.2f %7.2f --- %f ---- %f\n",
								tiePoints.x[i],
								tiePoints.y[i],tiePoints.z[i],tiePoints.r[i],tiePoints.a[i],
								tiePoints.phase[i],tiePoints.vyra[i]);
	fprintf(stderr,"\n");

	/*
	  Estimate baseline solution and output to stdout
	*/
	computeBaseline(&tiePoints,inputImage);
	return; 
}





static void  readArgs(int argc,char *argv[],int *imageFlag,int *passType,
		      int *noDEM,int *noRamp,char **demFile,char **inputFile,
		      char **tiePointFile,char **phaseFile, char **baselineFile,int *dBpFlag,
		      int *imageCoords, int *motionFlag,int *timeReverseFlag,
		      double *nDays,int *quadB, double *stdLat,  int *bnbpFlag, int *bpFlag, int *bnbpdBpFlag,
		      int *bpdBpFlag, int *vrFlag,char **shelfMaskFile) 
{
	int filenameArg;
	char *argString;
	int i,n;

	if( argc < 5 || argc > 20 ) usage();        /* Check number of args */ 

	*imageFlag = DESCENDING;  /* Default */
	*passType = -1;  
	n = argc - 5;
	*noDEM = TRUE;
	*noRamp = FALSE;
	*bpFlag = FALSE;
	*bnbpFlag = FALSE;
	*bpdBpFlag=FALSE;
	*bnbpdBpFlag=FALSE;
	*dBpFlag = FALSE;
	*imageCoords = FALSE;
	*motionFlag = FALSE;
	*timeReverseFlag = FALSE;
	*shelfMaskFile=NULL;
	*nDays = 3.0;
	*stdLat = 70.0;
	*quadB=FALSE;
	*vrFlag=FALSE;
	for(i=1; i <= n; i++) {
		argString = strchr(argv[i],'-');  
		if(strstr(argString,"descending") != NULL && *passType < 0)
			*passType = DESCENDING;
		else if(strstr(argString,"ascending") != NULL && *passType < 0)
			*passType= ASCENDING;
		else if(strstr(argString,"noDEM") != NULL ) 
			*noDEM = TRUE;
		else if(strstr(argString,"noRamp") != NULL ) 
			*noRamp = TRUE;
		else if(strstr(argString,"shelfMask") != NULL) {
			*shelfMaskFile =  argv[i+1]; i++;
		} else if(strstr(argString,"quadB") != NULL ) {
			*quadB = TRUE;
			if(*bpFlag==TRUE || *bnbpFlag==TRUE ||
			   *bpdBpFlag==TRUE || *bnbpdBpFlag == TRUE) 
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -"
				      "incompatable");
		} else if(strstr(argString,"bnbpOnly") != NULL ) {
			*bnbpFlag = TRUE;
			if(*bpFlag==TRUE || *quadB==TRUE ||
			   *bpdBpFlag==TRUE || *bnbpdBpFlag == TRUE) 
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -"
				      "incompatable");
		} else if(strstr(argString,"bpOnly") != NULL ) {
			*bpFlag = TRUE;
			if(*quadB==TRUE || *bnbpFlag==TRUE || 
			   *bpdBpFlag==TRUE || *bnbpdBpFlag == TRUE) 
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -"
				      "incompatable");
		}  else if(strstr(argString,"bnbpdBpOnly") != NULL ) {
			*bnbpdBpFlag = TRUE;
			if(*bpFlag==TRUE || *quadB==TRUE ||
			   *bpdBpFlag==TRUE || *bnbpFlag == TRUE) 
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -"
				      "incompatable");
		} else if(strstr(argString,"bpdBpOnly") != NULL ) {
			*bpdBpFlag = TRUE;
			if(*bpFlag==TRUE || *quadB==TRUE ||
			   *bnbpFlag==TRUE || *bnbpdBpFlag == TRUE) 
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -"
				      "incompatable");
		}  else if(strstr(argString,"imageCoords") != NULL ) 
			*imageCoords = TRUE;
		else if(strstr(argString,"vr") != NULL ) 
			*vrFlag = TRUE;
		else if(strstr(argString,"dBp") != NULL ) 
			*dBpFlag = TRUE;
		else if(strstr(argString,"center") != NULL)
			fprintf(stderr,"ignoring obsolete center flag\n");
		else if(strstr(argString,"motion") != NULL)
			*motionFlag = TRUE;
		else if(strstr(argString,"timeReverse") != NULL)
			*timeReverseFlag = TRUE;
		else if(strstr(argString,"nDays") != NULL) {
			sscanf(argv[i+1],"%lf",nDays); i++; }
		else if(strstr(argString,"stdLat") != NULL) {
			sscanf(argv[i+1],"%lf",stdLat); i++; }
		else if(strstr(argString,"rPix") != NULL) {
			sscanf(argv[i+1],"%lf",&RangePixelSize); i++; }
		else if(strstr(argString,"aPix") != NULL) {
			sscanf(argv[i+1],"%lf",&AzimuthPixelSize); i++; }
		else if(i != n) usage();  
	}
	if(*dBpFlag == TRUE && *noRamp == TRUE) {
		fprintf(stderr,"*** -dBp incompatible with -noRamp ***");
		usage();
	} else if(*imageCoords == TRUE && *motionFlag == TRUE) {
		fprintf(stderr,"*** motionFlag incompatible with imageCoords ***");
		usage();
	}
	if(*passType < 0) *passType = DESCENDING;
	*imageFlag = PHASE;
	*demFile = NULL;
	if(*noDEM == FALSE) *demFile = argv[argc-5];
	*inputFile = argv[argc-4];
	*tiePointFile = argv[argc-3];
	*phaseFile = argv[argc-2];
	*baselineFile = argv[argc-1];
	return;
}

 
static void usage()
{ 
	error(
	      "\n\n%s\n%s\n%s\n\n%s\n\n%s\n%s\n%s\n\n%s\n\n%s\n%s\n%s\n%s\n\n"
	      "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
	      "Compute baseline for given tiepoints.",
	      "At least 4 tiepoints must be given in the tiepoints file",
	      "Output is to stdout",
	      "Usage:",
	      " tiepoints  -vr -imageCoords -motion  -dBp -passType -noRamp -noDEM \\",
	      "            -shelfMask shelfMask -timeReverse -nDays nDays -quadB -rPix rPix -aPix aPix \\",
	      "            dem geoInput tiepointsFile uwPhase baselineFile",
	      "where",
	      "  bpOnly      estimate only bp",
	      "  bnbpOnly    estimate only bn,bp",
	      "  bpdBpOnly   estimate only bp,dBp",
	      "  bnbpdBpOnly estimate only bn,bp,dBp",   
	      "  imageCoords  = tiepoints in imageCoords instead of lat/lon",
	      "  passType     = DESCENDING (default) or ASCENDING",
	      "  noRamp       = if set, treat apply cos(thetad) depedence to omegaA",
	      "  dBp          = use deltaBp in place of omegaA", 
	      "  noDEM        = set flag if no DEM used (obsolete - DEM NEVER USED)",
	      "  shelfMask    = shelfMask file to indicate use tidal corrections",
	      "  timeReverse  = flag to reverse time when order of orbits switched",
	      "  nDays        = temporals basline in days (default=3)",
	      "  quadB        = flag for quadratic fit",
	      "  rPix, aPix   = Range/Azimuth single look pixel size",
	      "  dem          = (OMIT if dem not used) dem file for tiepoint elevations",
	      "  motion       = moving tiepoints, tiepoint file lat,lon,z,vx,vy,vz (m/yr)",
	      "  vr           = moving tiepoints, tiepoint file lat,lon,z,vr,vz (m/yr)",
	      "  geoInputFile = geo params file",
	      "  tiepointFile = Tiepoint location file in (lat,lon) or (lat,lon,z)",
	      "  uwPhaseFile  = unwrapped phase image",
	      "  baselineFile = baseline used to unwrap phases");
}
