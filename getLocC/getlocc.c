#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#undef MNST
#include "math.h"
#include "mosaicSource/common/common.h"

/*
  Create geodat file from .par file from clw processor.
*/

static void readArgs(int argc,char *argv[], int *nlr, int *nla,  int *noffset, double *squint, char **parFile,
		     char **geodatFile, int *passType, double *lookDir, int *squintTime, int *ersFix, int *surveyFlag,float *lambda);    
static void usage();
static void echoSarD(FILE *fp, SARData *sarD,stateV *sv,int nlr,int nla);
/* These routines are largely fixed to this program and not used elsewhere - so not bothering with include file*/

void correctTime(SARData *sarD, double *squint, double noffset, double *tskew,double *toffset, int squintTime);
void centerLL(SARData *sarD, stateV *sv, int nla,double *lat, double *lon,double deltaT);
void glatlon( SARData *sarD, stateV *sv, int nlr, int nla, double ***lat1,double ***lon1,int ma, int mr, double deltaT);     
/* 
   Global variables definitions (NOT USED, NEED FOR LINKING ERS CODE
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
void *offBufSpace1,*offBufSpace2,*offBufSpace3,*offBufSpace4;
void *lBuf1,*lBuf2,*lBuf3,*lBuf4;
char *Abuf1,*Abuf2,*Dbuf1,*Dbuf2;
int llConserveMem=999; /* Kluge to maintain backwards compat 9/13/06 */

void main(int argc, char *argv[])
{   
	FILE *fp;
	int nlr, nla;  /* Number of range and azimuth looks */
	int noffset;
	double theta,phic,rearth;
	char *parFile, *geodatFile; /* Input and output files */
	int passType;  /* 0 -> descending, 1 -> ascending */
	double tskew,toffset, squint, tdum1,tdum2,tdum3;
	double **lat,**lon;
	double latc=45.0,lonc;
	double lookDir;
	double slrpix;
	float lambda;
	int ersFix;
	int squintTime;
	int surveyFlag;
	int i;
	int nMLR,nMLA;
	double deltaAzML,deltaRgML;
	double mldT;
	SARData sarD;
	stateV sv;
    
	tskew=0.;
	readArgs(argc,argv,&nlr,&nla,&noffset,&squint, &parFile, &geodatFile, &passType,&lookDir, &squintTime,&ersFix,&surveyFlag,&lambda);
	sarD.lambda = lambda;
	/* Read old style gamm cw par file - no keywords */
	readOldPar(parFile,&sarD,&sv);
	nMLR=sarD.nSlpR/nlr; /* Mulitlook size */
	nMLA=sarD.nSlpA/nla;
	deltaAzML=nla * sarD.slpA;
	deltaRgML=nlr * sarD.slpR;	
	/*
	  Make sure range convention obeyed.
	*/
	sarD.lookDir=lookDir;
	/* Sort of obsolete, but just in case */
	if(surveyFlag==TRUE) sarD.prf /=8.;
	/*
	  Make sure ranges centered on multilook pixels rather than sl pixels
	*/
	sarD.rn +=(nlr-1)*sarD.slpR*0.5;
	sarD.rf =sarD.rn + (nMLR-1) * deltaRgML;
	sarD.rc = (sarD.rn + sarD.rf)*0.5;
        /* Nearly obsolete, keeping just in case */
	if(ersFix==TRUE) {
		fprintf(stderr,"*** APPLYING ERS FIX *** \n");     
		sarD.rn -= 54.0;  sarD.rc -= 54.0;      sarD.rf -= 54.0;             
	}
	/* Start writing geodat file */
	fp = fopen(geodatFile,"w");
	echoSarD(fp, &sarD,&sv,nlr,nla);
	/*
	  Correct starting time for squint and offset
	*/
	fprintf(stderr,"; Skew offset (s), squint (deg) : %f  %lf\n",tskew,squint);
	correctTime(&sarD,&squint,(double)noffset,&tskew,&toffset,squintTime);
	fprintf(fp,"; Offset of first recordin complex image (s) : %f\n",toffset);
	fprintf(stderr,"; Skew offset (s), squint (deg) : %f  %lf\n",tskew,squint);
	fprintf(fp,"; Skew offset (s), squint (deg) : %f  %lf\n",tskew,squint);
	/* Indicate pass type - asc/desc */
	if(passType == DESCENDING) fprintf(fp,";\n; Descending Pass \n");
	else if(passType == ASCENDING) fprintf(fp,";\n; Ascending Pass \n");
	else error("getlocc: Invalid Pass Type");
	/*
	  Range/azimuth sizes
	*/
	fprintf(fp,";\n; rangesize,azimuthsize,nrangelooks,nazimuthlooks\n;\n");
	fprintf(fp,"%i  %i  %i  %i\n",nMLR,nMLA,nlr,nla);
	/* Center point ll */
	centerLL(&sarD,&sv,nla,&latc,&lonc,0.0);
	/*
	  Note added mldT add, then subtract to be consistent time to glat/lon corrected  for center of multiLook pixel.
	*/
	mldT=0.5*(nla-1);
	tdum1=0; tdum2=0; tdum3=0;
	correctTime(&sarD,&tdum2,mldT,&tdum3,&tdum1,squintTime);
	glatlon(&sarD,&sv,nlr,nla,&lat,&lon,3,3,0.0);
	correctTime(&sarD,&tdum2,-mldT,&tdum3,&tdum1,squintTime);    /* Undo correction so time is SLC not ML */
	/*
	  Compute look and incidence angles
	*/
	rearth = earthRadius(latc*DTOR,EMINOR,EMAJOR) * KMTOM;
	theta = sarD.rc * sarD.rc + 2.0 * sarD.H * rearth + sarD.H * sarD.H;
	theta = theta/(2.0 * sarD.rc * (rearth + sarD.H));
	theta = acos(theta);
	fprintf(stdout,"theta/latc %f %f\n",theta*RTOD,latc);
	phic = ((rearth+sarD.H)/rearth) * sin(theta);
	phic = asin(phic)*RTOD;
	/* Ouput info */
	fprintf(fp,";\n; ReMajor, ReMinor, Rc, phic, h\n;\n");
	fprintf(fp,"%12.6f %12.6f %13.7f %8.4f %13.7f\n",  EMAJOR,EMINOR,sarD.rc*MTOKM,phic,sarD.H*MTOKM);    
	fprintf(fp,";\n; ll,lr,ul,ur\n;\n");
        /* Order corner points for ll,lr.. convention */    
	if(passType == DESCENDING) {
		fprintf(fp,"%13.8f %13.8f\n",lat[2][2],lon[2][2]);
		fprintf(fp,"%13.8f %13.8f\n",lat[2][0],lon[2][0]);
		fprintf(fp,"%13.8f %13.8f\n",lat[0][2],lon[0][2]);
		fprintf(fp,"%13.8f %13.8f\n",lat[0][0],lon[0][0]);
	} else if(passType == ASCENDING) {
		fprintf(fp,"%13.8f %13.8f\n",lat[0][0],lon[0][0]);
		fprintf(fp,"%13.8f %13.8f\n",lat[0][2],lon[0][2]);
		fprintf(fp,"%13.8f %13.8f\n",lat[2][0],lon[2][0]);
		fprintf(fp,"%13.8f %13.8f\n",lat[2][2],lon[2][2]);
	} else error("*** INVALID PASS TYPE ***\n");
	/* Output Center point */
	fprintf(fp,"%13.8f %13.8f\n",latc,lonc);
	/*
	  Other parameters
	*/
	fprintf(fp,";\n; Range/azimuth single look pixel sizes \n;\n");
	fprintf(fp,"%f  %f\n",sarD.slpR,sarD.slpA);
	if(passType == DESCENDING) fprintf(fp,";\ndescending\n");
	else fprintf(fp,";\nascending\n");
	/* Look dir */
	fprintf(fp,";\n; Look direction\n;\n");
	if(lookDir < 0) fprintf(fp,"left\n");
	else fprintf(fp,"right\n");
	/*
	  Write flag to indicate state vectors present - the notion they wouldn't be is 20 years obsolete
	*/
	fprintf(fp,";\n; Flag to indicate state vectors and associated data\n;\n");
	fprintf(fp,"state\n");
	fprintf(fp,"; time after squint and skew corrections\n %i %i %f\n", sarD.hr,sarD.min,sarD.sec);
	fprintf(fp,"; prf \n%f\n",sarD.prf);
	fprintf(fp,"; wavelength\n%f\n",sarD.lambda);
	/*
	  Write state vectors
	*/
	fprintf(fp,"; number of state vectors\n%i\n",sv.nState);
	fprintf(fp,"; time of first vector \n%f\n",sv.t0);
	fprintf(fp,"; state vector interval \n%f\n",sv.deltaT);
	fprintf(fp,"; state vectors \n");
	for(i=1; i <= sv.nState; i++) {
		fprintf(fp,"%.9e  %.9e  %.9e\n", sv.x[i],sv.y[i],sv.z[i]);
		fprintf(fp,"%.9e  %.9e  %.9e\n", sv.vx[i],sv.vy[i],sv.vz[i]);
	}
	return; 
}


static void readArgs(int argc,char *argv[], int *nlr, int *nla, int *noffset, double *squint, char **parFile,
		     char **geodatFile, int *passType,double *lookDir, int *squintTime, int *ersFix, int *surveyFlag, float *lambda)
{
	int i;
	if( argc < 9 || argc > 13 ) usage();        /* Check number of args */ 
	*squintTime=FALSE;
	*ersFix=FALSE;
	*surveyFlag=FALSE;
	*lambda=LAMBDAERS1;
	if(argc>=10){
		for(i=1; i < argc-8; i++) {
			if(strstr(argv[i],"squintTime") != NULL) *squintTime=TRUE;
			else if(strstr(argv[i],"ersFix") != NULL) *ersFix=TRUE;
			else  if(strstr(argv[i],"survey") != NULL) *surveyFlag=TRUE;
			else if(strstr(argv[i],"lambda")!= NULL) {i++; sscanf(argv[i],"%f",lambda);} 
			else usage();
		}
	}
	sscanf(argv[argc-8],"%i",nlr);
	sscanf(argv[argc-7],"%i",nla);
	sscanf(argv[argc-6],"%i",noffset);
	fprintf(stderr,"%i %i %i\n",*nlr,*nla,*noffset);
	sscanf(argv[argc-5],"%lf",squint);
	*parFile = argv[argc-4];
	*geodatFile = argv[argc-3];
	sscanf(argv[argc-2],"%i",passType);
	sscanf(argv[argc-1],"%lf",lookDir);
	if(*lookDir < 0.0 ) *lookDir=-1.0;
	else *lookDir=1.0;
	return;
}

 
static void usage()
{ 
	error("getloc -lambda lambda -ersFix -squintTime nlr nla noffset squintAngle parfile outfile\n\n"
	      " lambda = wavelength\n"
	      " passtype = 0 for desc and 1 for right\n"
	      " squintTime flag to specify squint in time\n"
	      "lookDir = +1.0 for right and -1.0 for left");
}



/* This program writes the largely comment information to the geodat file */
static void echoSarD(FILE *fp, SARData *sarD,stateV *sv,int nlr,int nla)
{
	char *months[12]={"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"};  
	fprintf(fp,"; Image name: %s\n",sarD->label);
	fprintf(fp,"; Image date: %02d %s %4i\n",sarD->day,months[sarD->month-1],sarD->year);
	fprintf(fp,"; Image time: %i %i %f\n",sarD->hr,sarD->min,sarD->sec);
	fprintf(fp,"; Nominal center lat,lon: %f %f\n",0.,0.);
	fprintf(fp,"; track direction: %f\n",0.);
	fprintf(fp,"; S/C altitude: %f\n",sarD->H);
	fprintf(fp,"; Average height above terrain: %f\n",0.0);
	fprintf(fp,"; Vel along track: %f\n",0.0);
	fprintf(fp,"; PRF :   %f\n",sarD->prf);
	fprintf(fp,"; near/cen/far range : %f %f %f\n",sarD->rn,sarD->rc,sarD->rf);
	fprintf(fp,"; Range pixel spacing :   %f\n",sarD->slpR*(double)(nlr));
	fprintf(fp,"; Number of looks (rg,az) :   %i %i\n",nlr,nla);
	fprintf(fp,"; Azimuth pixel spacing :   %f\n",sarD->slpA*(double)nla);
	fprintf(fp,"; Number of pixels (rg,az) :  %i  %i\n",sarD->nSlpR/nlr,sarD->nSlpA/nla);
	fprintf(fp,"; Number of state vectors :   %i\n",sv->nState);
	fprintf(fp,"; Start time of state vectors :   %f\n",sv->t0); 
	fprintf(fp,"; Interval between 2 state vectors :   %f\n",sv->deltaT);  
	fprintf(fp,"; Look direction  :   %f\n",sarD->lookDir);  
}
