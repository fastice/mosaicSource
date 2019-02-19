#define _XOPEN_SOURCE
#include <time.h>
#include "clib/standard.h"
#include "stdio.h"
#include <stdlib.h>
#include <string.h>
#include "mosaicSource/common/common.h"

#define MIN(a,b)  ( ( (a) < (b) ) ? (a) : (b) )

static void getLabel(FILE *fpPar,SARData *sarD,int quiet) {
	char buffer[513];	
	/* Label */
	if(fgets(buffer,512,fpPar) == NULL) error("error reading par title");	/* title of the scene and processing     parameters */
	buffer[strlen(buffer)-1]='\0';/* strip off newline character */    
	sarD->label=malloc(strlen(buffer)+1);
	strcpy(sarD->label,buffer);
	if( quiet == FALSE ) fprintf(stderr,"sarD->label : %s\n",sarD->label);
}

static void getParDate(FILE *fpPar, SARData *sarD,int quiet) {
	struct tm tp;
	char buffer[20];	
	char *cp;
	char *months[12]={"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"};  
	if(fgets(buffer,16,fpPar) == NULL) error("error reading date");	
	buffer[strlen(buffer)-1]='\0';
	cp=(char *)strptime(buffer,"%e %B %Y",&tp);
	sarD->day=tp.tm_mday;
	sarD->year=tp.tm_year+1900;
	sarD->month=tp.tm_mon+1;
	if(quiet == FALSE) fprintf(stderr,"  %i %i %i %s\n",sarD->year,sarD->month,sarD->day,months[tp.tm_mon]);
}

static void getParTime(FILE *fpPar, SARData *sarD,int quiet) {
	char buffer[25];	
	if(fgets(buffer,24,fpPar) == NULL) error("error reading date");	
	buffer[strlen(buffer)-1]='\0';
	fprintf(stderr,"%s\n",buffer);
	if(sscanf(buffer,"%i%i%lf",&(sarD->hr),&(sarD->min),&(sarD->sec)) != 3) error("parsing date");
	if(quiet == FALSE) fprintf(stderr,"  %i %i %f\n",sarD->hr,sarD->min,sarD->sec);
}

static void getSatH(FILE *fpPar, SARData *sarD,int quiet) {
	char buffer[25];	
	if(fgets(buffer,24,fpPar) == NULL) error("error reading date");	
	buffer[strlen(buffer)-1]='\0';
	fprintf(stderr,"%s\n",buffer);
	if( sscanf(buffer,"%lf",&(sarD->H) ) != 1) error("parsing sat height");
	if(quiet == FALSE) fprintf(stderr,"Sat Height  %lf\n",sarD->H);
}

static void skipLine(FILE *fpPar,int quiet) {
	char buffer[512];
	if(fgets(buffer,512,fpPar) == NULL) error("error reading date");
	if(quiet == FALSE) fprintf(stderr,"skip:  %s",buffer);	
}

static void getPRF(FILE *fpPar, SARData *sarD, int quiet) {
	if(fscanf(fpPar,"%lf",&(sarD->prf)) != 1) error("parsing date");
	if(quiet == FALSE) fprintf(stderr,"prf %lf\n",sarD->prf);	
}

static void getDoppler(FILE *fpPar, SARData *sarD, int quiet) {
	if(fscanf(fpPar,"%le%le%le%le",&(sarD->fd[0]),&(sarD->fd[1]),&(sarD->fd[2]),&(sarD->fd[3]))    != 4) error("error parsing Doppler");
	if(quiet == FALSE) fprintf(stderr,"Doppler: %le %le %le %le\n",sarD->fd[0],sarD->fd[1], sarD->fd[2], sarD->fd[3] );	
}

static void getEchoTD(FILE *fpPar, SARData *sarD, int quiet) {
	if(fscanf(fpPar,"%lf",&(sarD->echoTD) )    != 1) error("error parsing time delay");
	if(quiet == FALSE) fprintf(stderr,"TD: %le \n",sarD->echoTD );	
}
static void getRawRange(FILE *fpPar, SARData *sarD, int quiet) {
	if(fscanf(fpPar,"%lf%lf%lf",&(sarD->rawRange[0]),&(sarD->rawRange[1]),&(sarD->rawRange[2]))   != 3) error("error parsing Raw ranges");
	if(quiet == FALSE) fprintf(stderr,"Raw Range (n/c/f): %le %le %le \n",sarD->rawRange[0],sarD->rawRange[1], sarD->rawRange[2]); 
}
static void getSlcRange(FILE *fpPar, SARData *sarD, int quiet) {
	if(fscanf(fpPar,"%lf%lf%lf",&(sarD->rn),&(sarD->rc),&(sarD->rf))   != 3) error("error parsing slc ranges");
	if(quiet == FALSE) fprintf(stderr,"slc Range (n/c/f): %le %le %le \n",(sarD->rn),(sarD->rc), (sarD->rf)); 
}

static void getRangeRes(FILE *fpPar, SARData *sarD, int quiet) {
	double dum;
	if(fscanf(fpPar,"%lf%lf",&(sarD->slpR),&(dum) )  != 2) error("error parsing slc range res");
	if(quiet == FALSE) fprintf(stderr,"Range pixel %lf \n",sarD->slpR); 
}

static void getAzimuthRes(FILE *fpPar, SARData *sarD, int quiet) {
	double dum;
	if(fscanf(fpPar,"%lf%lf",&(sarD->slpA),&(dum) )  != 2) error("error parsing slc az res");
	if(quiet == FALSE) fprintf(stderr,"Azimuth pixel %lf \n",sarD->slpA); 
}

static void getSlcSize(FILE *fpPar, SARData *sarD, int quiet) {
	double dum;
	if(fscanf(fpPar,"%d%d",&(sarD->nSlpR),&(sarD->nSlpA) )  != 2) error("error parsing slc size");
	if(quiet == FALSE) fprintf(stderr,"Range/Azimuth pixels  %d %d \n",sarD->nSlpR,sarD->nSlpA); 
}

static void getSlcNState(FILE *fpPar, stateV *sv, int quiet) {
	if(fscanf(fpPar,"%d",&(sv->nState) )  != 1) error("error parsing n State");
	if(quiet == FALSE) fprintf(stderr,"N state %d \n",sv->nState);
}

static void getSlcT0State(FILE *fpPar, stateV *sv, int quiet) {
	if(fscanf(fpPar,"%lf",&(sv->t0) )  != 1) error("error parsing t0 state");
	if(quiet == FALSE) fprintf(stderr,"t0 state %f \n",sv->t0);
}

static void getSlcDeltaTState(FILE *fpPar, stateV *sv, int quiet) {
	if(fscanf(fpPar,"%lf",&(sv->deltaT) )  != 1) error("error parsing deltaT state");
	if(quiet == FALSE) fprintf(stderr,"t0 state %f \n",sv->deltaT);
}

static void getSlcStateV(FILE *fpPar, stateV *sv, int quiet) {
	int i;
	for(i=1; i <= sv->nState; i++ ) {
		sv->times[i] = sv->t0 + (i-1)*sv->deltaT;		
		if(fscanf(fpPar,"%lf%lf%lf",&(sv->x[i]),&(sv->y[i]),&(sv->z[i])  )  != 3) error("error parsing pos state");
		if(fscanf(fpPar,"%lf%lf%lf",&(sv->vx[i]),&(sv->vy[i]),&(sv->vz[i])  )  != 3) error("error parsing vel state");
		if(quiet == FALSE) fprintf(stderr,"state %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e \n",
					   sv->x[i],sv->y[i],sv->z[i],sv->vx[i],sv->vy[i],sv->vz[i]);
	}
}

	

/*
  read old par file 
*/
void readOldPar(char *parFile, SARData *sarD, stateV *sv)
{
	char buffer[1024];
	int quiet;
	int nState;
	int i;
	char *ret;
	int fret;
	FILE *fpPar;

	quiet=TRUE; /* Set to FALSE for debugging */
	fprintf(stderr,"-- %s\n",parFile);
	fpPar = openInputFile(parFile);
	/* Label */
	getLabel(fpPar,sarD,quiet);
	/* date */
	getParDate(fpPar,sarD,quiet);
	/* time */
	getParTime(fpPar,sarD,quiet);	
	/* skip pol/channel, and lat/lon */
	skipLine(fpPar,quiet); /* not used */
	skipLine(fpPar,quiet); /* lat lon computed elsewhere */
	skipLine(fpPar,quiet); /* track angle computed locally */	
        getSatH(fpPar,sarD,quiet);
	skipLine(fpPar,quiet); /* no need for avg terrain*/
	skipLine(fpPar,quiet); /* no need for center pos */
	skipLine(fpPar,quiet); /* no need for center vel */
	skipLine(fpPar,quiet); /* no need for center acc */
	/* prf */
	getPRF(fpPar,sarD,quiet);
	/* doppler */
	getDoppler(fpPar,sarD,quiet);
	/* echo time delay, not used currently*/
	getEchoTD(fpPar,sarD,quiet);
	skipLine(fpPar,quiet); /* first skip to clear sscanf*/
	skipLine(fpPar,quiet); /* no current need for range offsets*/			
	/* raw data ranges, not used currently*/
	getRawRange(fpPar,sarD,quiet);
	/* slc data ranges */
	getSlcRange(fpPar,sarD,quiet);
	getRangeRes(fpPar,sarD,quiet);
	skipLine(fpPar,quiet); /* flush res line */
	skipLine(fpPar,quiet); /* flush res line */
	skipLine(fpPar,quiet); /* flush res line */
	skipLine(fpPar,quiet); /* flush res line */
	skipLine(fpPar,quiet); /* flush res line */
	skipLine(fpPar,quiet); /* flush res line */
	skipLine(fpPar,quiet); /* flush res line */
	skipLine(fpPar,quiet); /* flush res line */
	skipLine(fpPar,quiet); /* flush res line */
	getAzimuthRes(fpPar,sarD,quiet);
	getSlcSize(fpPar,sarD,quiet);
	/*
	  State vectors
	*/
	getSlcNState(fpPar,sv,quiet);
	getSlcT0State(fpPar,sv,quiet);
	getSlcDeltaTState(fpPar,sv,quiet);
	getSlcStateV(fpPar,sv,quiet);			
	fclose(fpPar);
	return;
}



