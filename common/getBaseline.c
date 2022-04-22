#include "stdio.h"
#include "string.h"
#include "common.h"
#include "math.h"

/*
  Input baseline info.
*/
void getBaseline(char *baselineFile, vhParams *params, int noPhase)
{
	FILE *fp;
	double BnEst,BpEst,dBnEst,dBpEst,dBpEstQ,dBnEstQ;
	int  nBaselines;
	char *tmp;
	double dum;
	double fdum1,fdum2,fdum3,fdum4,fdum5,fdum6;
	int lineCount=0, eod, special;
	int i,j;
	int ci;
	char line[256];
	/*
	  Input parm info
	*/
	for(i=1; i <=6; i++)    for(j=1; j <=6; j++) params->C[i][j]=0;
	params->sigma=PI; /* Default value */
	if( strstr(baselineFile,"nobaseline") != NULL || noPhase == TRUE) {
		/* No baseline used for this solution so return */
		params->Bn  = 0;    params->Bp = 0;
		params->dBn = 0;   params->dBp  = 0;
		params->dBnQ = 0;  params->dBpQ  = 0;
		return;
	}
	fp = openInputFile(baselineFile);
	/*
	  Skip past initial data lines
	*/
	lineCount=getDataString(fp,lineCount,line,&eod);
	lineCount=getDataString(fp,lineCount,line,&eod);
	if( sscanf(line,"%i",&nBaselines) != 1)   error("getBaseline -- Missing baseline params at line: %i of %s",lineCount, baselineFile);
	if(nBaselines < 1 || nBaselines > 3)   error("getBaseline -- invalid number of baselines at line %i of %s\n",  lineCount, baselineFile);
	if(nBaselines > 1) lineCount=getDataString(fp,lineCount,line,&eod);
	if(nBaselines > 2) lineCount=getDataString(fp,lineCount,line,&eod);
	/*
	  Read covariance matrix if there is one. 
	*/
	special=TRUE;
	while(special==TRUE) {
		lineCount=getDataStringSpecial(fp,lineCount,line,&eod,'*',&special);
		ci=0;
		tmp=strstr(line,"C_1");   if(tmp !=NULL) {  tmp+=3; ci=1;}
		if(tmp==NULL)  { tmp=strstr(line,"C_2");   if(tmp !=NULL) { tmp+=3; ci=2;}}
		if(tmp==NULL)  { tmp=strstr(line,"C_3");   if(tmp !=NULL) { tmp+=3; ci=3;}}
		if(tmp==NULL) { tmp=strstr(line,"C_4");   if(tmp !=NULL) { tmp+=3; ci=4;}}
		if(tmp==NULL) {tmp=strstr(line,"C_5");   if(tmp !=NULL) { tmp+=3; ci=5;}}
		if(tmp==NULL) { tmp=strstr(line,"C_6");   if(tmp !=NULL) { tmp+=3; ci=6;}}

		if(ci > 0)  {
			sscanf(tmp,"%lf%lf%lf%lf%lf%lf",&fdum1,&fdum2,&fdum3,&fdum4,&fdum5,&fdum6);
			params->C[ci][1]=fdum1; 		    params->C[ci][2]=fdum2;
			params->C[ci][3]=fdum3; 		    params->C[ci][4]=fdum4;
			params->C[ci][5]=fdum5; 		    params->C[ci][6]=fdum6;
		} else if( tmp=strstr(line,"sigma*sqrt(X2/n)=")) {
			tmp+=strlen("sigma*sqrt(X2/n)=");  sscanf(tmp,"%lf",&fdum1); params->sigma=fdum1;
		}
	}
	fprintf(stderr,"sigma*sqrt(X2/n) = %lf\n",params->sigma);
	/*
	  Input baseline estimated with tiepoints.
	*/
	if( sscanf(line,"%lf%lf%lf%lf%lf%lf",   &dum,&dum,&dum,&dum,&dum,&dum) == 6 ) {
		fprintf(stderr,"\n***BASELINE FIT WAS QUADRATIC***\n");
		sscanf(line,"%lf%lf%lf%lf%lf%lf",&BnEst,&BpEst,&dBnEst,&dBpEst, &dBnEstQ,&dBpEstQ);
	} else {
		if( sscanf(line,"%lf%lf%lf%lf",&BnEst,&BpEst,&dBnEst,&dBpEst) != 4)
			error("%s %i", "getBaseline  -- Missing baseline params at line:",lineCount);
		dBnEstQ = 0.0;
		dBpEstQ = 0.0;
	} 
	/* populate baseline params structure */
	params->Bn  = BnEst; 
	params->Bp = BpEst;
	params->dBn = dBnEst; 
	params->dBp  = dBpEst;
	params->dBnQ = dBnEstQ; 
	params->dBpQ  = dBnEstQ;
	fclose(fp);

}

