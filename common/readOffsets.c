#include "stdio.h"
#include"string.h"
#include <stdlib.h>
#include "common.h"
#include <libgen.h>

static void readCov(FILE *fp, int n,double C[7][7], double *sigmaResidual ,char *line);
static void readOffsetFile(float **data,int nr, int na,char *offsetFile);

/*
   Read the offset data and paramter files
 */
void readOffsetDataAndParams(Offsets *offsets) {
	readBothOffsets(offsets);
	fprintf(stderr,"Read az est\n");
	getAzParams(offsets);
	fprintf(stderr,"Read aBaseline\n");
	getRParams(offsets);
	fprintf(stderr,"Offsets and paramters read\n");	
	if(offsets->deltaB != DELTABNONE && offsets->geo2==NULL)
		error("offsets deltaB set but no second geodat for %s\n",offsets->rFile);
}

static void initOffsetBuffers(Offsets *offsets,int mode) {
	extern void *offBufSpace1,*offBufSpace2,*offBufSpace3,*offBufSpace4;
	extern void *lBuf1,*lBuf2,*lBuf3,*lBuf4;
	float *fBuf1,*fBuf2,*fBuf3,*fBuf4;
	int i,nr,na;
	nr=offsets->nr;
	na=offsets->na;
	if(mode != RGONLY) {
		offsets->da = (float **)lBuf1;
		offsets->sa = (float **)lBuf2;
		fBuf1=(float *)offBufSpace1;
		fBuf2=(float *)offBufSpace2;
	}
	if(mode != AZONLY) {
		offsets->dr = (float **)lBuf3;	
		offsets->sr = (float **)lBuf4;
		fBuf3=(float *)offBufSpace3;
		fBuf4=(float *)offBufSpace4;		
	}
	if(nr * na * 4 > MAXOFFBUF) error("Offsets Buffer Size exceeds MAXOFFBUF of %i\n",MAXOFFBUF);
	for(i=0; i < na; i++) {
		if(mode != RGONLY) {		
			offsets->da[i]=&(fBuf1[i*nr]);
			offsets->sa[i]=&(fBuf2[i*nr]);
		}
		if(mode != AZONLY) {		
			offsets->dr[i]=&(fBuf3[i*nr]);
			offsets->sr[i]=&(fBuf4[i*nr]);
		}
	}
}

static char *mergePath(char *file1,char *path) {
	char fileTmp[1024],*tmp,*merged;
	fileTmp[0]='\0';
	tmp=strcat(fileTmp, path);
	tmp=strcat(fileTmp, "/");
	tmp=strcat(fileTmp,file1);
	merged=(char *)malloc(strlen(fileTmp)+1);
	merged[0]='\0';
        return(strcat(merged,fileTmp));
}
/*
  Read the params
*/
static void readOffsetParams(char *datFile, Offsets *offsets) {
	FILE *fp;
	int rO,aO,nr,na;
	float deltaA,deltaR;
	double sigmaS,sigmaR;
	char line[1024],*tmp;
	char file1[512],file2[512],*path;
	int lineCount, eod;
	int nRead;

	/* Read params file */
	fp = openInputFile(datFile);
	lineCount=getDataString(fp,lineCount,line,&eod);
	nRead=sscanf(line,"%i%i%i%i%f%f%lf%lf",&rO,&aO,&nr,&na,&deltaR,&deltaA,&sigmaS,&sigmaR);

	if(nRead != 6 && nRead != 7 && nRead !=8)
		error("%s  %i of %s","readOffsets -- Missing image parameters at line:", lineCount,datFile); 
	else if(nRead == 6) { 
		fprintf(stderr,"**** WARNING-MISSING SIGMA STREAKS for %s\n",datFile);
		sigmaS=0.0;
		sigmaR=0.0;
	} else if(nRead ==7) sigmaR=0.0;
	/* load param in structure */
        offsets->rO=rO;
	offsets->aO=aO;
	offsets->deltaA = deltaA;
	offsets->deltaR = deltaR;
	offsets->sigmaStreaks=sigmaS;
	offsets->sigmaRange=sigmaR;
	fprintf(stderr,"SigmaStreaks/Range = %f %f %i\n",sigmaS,sigmaR,nRead);
	offsets->nr=nr;
	offsets->na=na;
	fprintf(stderr,"nr na %i %i\n",nr,na);
	/*
	   read geodat files if the exist
	 */
	tmp=fgets(line,1024,fp);
	offsets->geo1=NULL;
	offsets->geo2=NULL;	
	if(tmp != NULL) {
		sscanf(line,"%s %s %s",file1,file2,datFile);
	        path=dirname(datFile);
		offsets->geo1=mergePath(file1,path);
		offsets->geo2=mergePath(file2,path);
		fprintf(stderr,"geo1 %s geo2 %s\n",offsets->geo1,offsets->geo2);
	}
	fclose(fp);			  
}

/* 
	read a single offset file
*/
static void readOffsetFile(float **data,int nr, int na,char *offsetFile) {
	FILE *fp;
	int i;
	fprintf(stderr,"--- Reading offset file %s",offsetFile);
	fp = openInputFile(offsetFile);
	for(i=0; i < na; i++) freadBS(data[i],sizeof(float),nr,fp, FLOAT32FLAG);
	fclose(fp);
	fprintf(stderr,"- done - \n");	
}

/* 
	read offsets and error files 
*/
void readBothOffsets( Offsets *offsets) 
{
	char *datFile, buf[1024],bufa[1024],bufd[1024];
	char *eFileA,*eFileR;
	char *file;
	int i;
	/*
	  Read inputfile
	*/
	file=offsets->file;
	datFile=appendSuffix(file,".dat",buf);
	eFileA=appendSuffix(file,".sa",bufa);
	eFileR=appendSuffix(offsets->rFile,".sr",bufd);		
	/*
	 read offset param
	*/
	readOffsetParams(datFile,offsets);
	/*
	setup buffers
	*/
	initOffsetBuffers(offsets,RGANDAZ);
	/*
	  Read files
	*/
	readOffsetFile(offsets->da ,offsets->nr,offsets->na,offsets->file);
	readOffsetFile(offsets->dr ,offsets->nr,offsets->na,offsets->rFile);
	readOffsetFile(offsets->sa ,offsets->nr,offsets->na,eFileA);
	readOffsetFile(offsets->sr ,offsets->nr,offsets->na,eFileR);
}

static char *RgOffsetsParamName(char *rParamsFile,char *newFile,int deltaB) {
	char *suffix[3] = {"",".deltabp",".quad"};
	if(deltaB > DELTABQUAD || deltaB < DELTABNONE) error("invalide deltaB flag %i",deltaB);
	return(appendSuffix(rParamsFile,suffix[deltaB],newFile));
}


void getRParams( Offsets *offsets)
{
	FILE *fp;
	double bn,bp,dBn,dBp,rConst,dBnQ,dBpQ;
	int  nBaselines;
	double dum;
	int ci;
	int lineCount=0, eod, special;
	int i,j;
	char line[256], *tmp, paramFile[1024];
	/*
	  Input parm info
	*/
	RgOffsetsParamName(offsets->rParamsFile,paramFile,offsets->deltaB);
	fp = fopen(paramFile,"r");
	if(fp == NULL) {
		/* If quad or const doesn't exist revert to original */
		fp=openInputFile(offsets->rParamsFile);
		offsets->deltaB=DELTABNONE;
	}
	fprintf(stderr,"deltaB %i\n",offsets->deltaB);
	for(i=1; i<=6; i++) for(j=1; j<=6; j++) offsets->Cr[i][j]=0.0;
	/*
	  Skip past initial data lines
	*/
	lineCount=getDataString(fp,lineCount,line,&eod);
	lineCount=getDataString(fp,lineCount,line,&eod);
	if( sscanf(line,"%i",&nBaselines) != 1)   error("%s  %i", "getBaseline -- Missing baseline params at line:",lineCount);
	if(nBaselines < 1 || nBaselines > 3)   error("getBaseline -- invalid number of baselines at line %i\n",  lineCount);

	if(nBaselines > 1) lineCount=getDataString(fp,lineCount,line,&eod);
	if(nBaselines > 2) lineCount=getDataString(fp,lineCount,line,&eod);
	/*
	  Read covariance matrix if there is one. 
	*/
	readCov(fp,6,offsets->Cr, &(offsets->sigmaRresidual),line);
	for(i=1; i <=6; i++) fprintf(stderr,"%le %le %le %le %le %le \n",
		(offsets->Cr[i][1]),(offsets->Cr[i][2]),(offsets->Cr[i][3]),(offsets->Cr[i][4]),(offsets->Cr[i][5]),(offsets->Cr[i][6]));	
	fprintf(stderr,"range sigma*sqrt(X2/n) = %lf (m)\n",offsets->sigmaRresidual);
	/*
	  Input baseline estimated with tiepoints.
	*/
	/*lineCount=getDataString(fp,lineCount,line,&eod);*/
	dBpQ=0.0;
	dBnQ=0.0;
	if(sscanf(line,"%lf%lf%lf%lf%lf%lf%lf",&bn,&bp,&dBn,&dBp,&rConst,&dBnQ,&dBpQ) != 7) {
		if( sscanf(line,"%lf%lf%lf%lf%lf",&bn,&bp,&dBn,&dBp,&rConst) != 5) 
			error( "\n\ngetRoffsets:Invalid range offset baseline file\nFile: %s\nLine: %s", offsets->rParamsFile,line);
	}
	offsets->bn=bn; offsets->bp=bp; 
	offsets->dBn=dBn; offsets->dBp=dBp; offsets->rConst=rConst;  
	offsets->dBnQ=dBnQ; offsets->dBpQ=dBpQ;

	fprintf(stderr,"bn %f %f %f bp %f %f %f off %f\n",offsets->bn,offsets->dBn,
		offsets->dBnQ,offsets->bp,offsets->dBp,offsets->dBpQ,offsets->rConst);
	fclose(fp);

}


/*
  This reads the range offsets, but uses the azimuth  offsets buffer for the asc and the range for the descending
*/
void readRangeOrRangeOffsets( Offsets *offsets,int orbitType) 
{
	extern void *offBufSpace1,*offBufSpace2,*offBufSpace3,*offBufSpace4;
	extern void *lBuf1,*lBuf2,*lBuf3,*lBuf4;
	char *datFile,*file ,buf[1024],bufd[1024];
	char *eFileR;
	float *fBufR,*fBufER;
	int i;
	/*
	  Read inputfile
	*/
	file=offsets->file;
	datFile=appendSuffix(file,".dat",buf);
	readOffsetParams(datFile,offsets);
	/* 
	   Set up for offsets type 
	*/
	eFileR=offsets->rFile;   
	eFileR=appendSuffix(eFileR,".sr",bufd);
	fprintf(stderr,"--- Reading offset file %s %i %i\n",eFileR,offsets->nr,offsets->na);

	if(orbitType == ASCENDING) {
		offsets->dr = (float **)lBuf1;
		offsets->sr = (float **)lBuf2;
		fBufR=(float *)offBufSpace1;
		fBufER=(float *)offBufSpace2;
	} else {
		offsets->dr = (float **)lBuf3;
		offsets->sr = (float **)lBuf4;
		fBufR=(float *)offBufSpace3;
		fBufER=(float *)offBufSpace4;
	}
	if(offsets->nr * offsets->na * 4 > MAXOFFBUF) error("Offsets Buffer Size exceeds MAXOFFBUF of %i\n",MAXOFFBUF);
	for(i=0; i < offsets->na; i++) {
		offsets->dr[i]=&(fBufR[i*offsets->nr]);
		offsets->sr[i]=&(fBufER[i*offsets->nr]);
	}
	/*
	  Read files
	*/
	readOffsetFile(offsets->dr ,offsets->nr,offsets->na,offsets->rFile);
	readOffsetFile(offsets->sr ,offsets->nr,offsets->na,eFileR);		      
	fprintf(stderr,"error File %s\n",eFileR);
}
	      
static char *AzOffsetsParamName(char *aParamsFile,char *newFile,int deltaB) {
	char *suffix[3] = {"",".const",".svlinear"};
	if(deltaB > DELTABQUAD || deltaB < DELTABNONE) error("invalide deltaB flag %i",deltaB);
	return(appendSuffix(aParamsFile,suffix[deltaB],newFile));
}


/*
  Input azimuth parameter info. Modifed 3/1/16 to read in covariance matrix
*/
void getAzParams( Offsets *offsets)
{
	FILE *fp;
	int lineCount=0, eod;
	int i, j;
	char line[256], paramFile[1024];
	/*
	  Open az param file
	*/
	AzOffsetsParamName(offsets->azParamsFile,paramFile,offsets->deltaB);
	fp = fopen(paramFile,"r");
	if(fp == NULL) {
		/* If quad or const doesn't exist revert to original */
		fp=openInputFile(offsets->azParamsFile);
		offsets->deltaB=DELTABNONE;
		if(fp == NULL) error("getAzParams: Error opening %s\n",offsets->azParamsFile);
	}
	/*
	  Skip past initial data lines
	*/
	lineCount=getDataString(fp,lineCount,line,&eod);   
	/*
	  Read covariance matrix if there is one. 
	*/
	for(i=1; i<=4; i++)  for(j=1; j<=4; j++) offsets->Ca[i][j]=0.0;
	readCov(fp,4,offsets->Ca, &(offsets->sigmaAresidual),line);
	fprintf(stderr,"azimuth sigma*sqrt(X2/n) = %lf (m)\n",offsets->sigmaAresidual);
	/* read params from fit */
	if( sscanf(line,"%lf%lf%lf%lf",  &(offsets->c1),&(offsets->dbcds),&(offsets->dbhds), &(offsets->doffdx)) != 4) {
		error("getAzParams  -- Missing baseline params at line: %s %i\n %s",paramFile,lineCount,line); 
	}
	fclose(fp);
}


/*
   Read 4x4 or 6x6 cov matrix for the az/rg params file format
*/
static void readCov(FILE *fp, int n,double C[7][7], double *sigmaResidual,char *line) {
	int lineCount=0, eod;
	int count;
	int special, ci;
	char  *tmp;
	double fdum1,fdum2,fdum3,fdum4,fdum5,fdum6;
	count =0;
	special=TRUE;
	while(special == TRUE) {
		lineCount=getDataStringSpecial(fp,lineCount,line,&eod,'*',&special);
		ci=0;
		tmp=strstr(line,"C_1");   if(tmp !=NULL) {  tmp+=3; ci=1;}
		if(tmp==NULL)  { tmp=strstr(line,"C_2");   if(tmp !=NULL) { tmp+=3; ci=2;}}
		if(tmp==NULL)  { tmp=strstr(line,"C_3");   if(tmp !=NULL) { tmp+=3; ci=3;}}
		if(tmp==NULL) { tmp=strstr(line,"C_4");   if(tmp !=NULL) { tmp+=3; ci=4;}}
		if(tmp==NULL) {tmp=strstr(line,"C_5");   if(tmp !=NULL) { tmp+=3; ci=5;}}
		if(tmp==NULL) { tmp=strstr(line,"C_6");   if(tmp !=NULL) { tmp+=3; ci=6;}}
		if(ci > 0)  {
			fdum1=0;fdum2=0;fdum3=0.0;fdum4=0.0;fdum5=0.0; fdum6=0.0;				
			if( n==4 ) sscanf(tmp,"%lf%lf%lf%lf",&fdum1,&fdum2,&fdum3,&fdum4);
			else sscanf(tmp,"%lf%lf%lf%lf%lf%lf",&fdum1,&fdum2,&fdum3,&fdum4,&fdum5,&fdum6);
			C[ci][1]=fdum1; 		    C[ci][2]=fdum2;
			C[ci][3]=fdum3; 		    C[ci][4]=fdum4;
			if( n == 6 ) { C[ci][5]=fdum5;  C[ci][6]=fdum6; }
		} else if( tmp=strstr(line,"sigma*sqrt(X2/n)=")) {
			tmp+=strlen("sigma*sqrt(X2/n)=");  sscanf(tmp,"%lf",sigmaResidual); 
		}
		/* Make sure no infintite loop if file wrong */
		count++;
		if( count > 8) error("Problem reading covariance matrix");
	}
}


/*
   Read azimuth offsets
 */
void readOffsets( Offsets *offsets) 
{
	char *datFile, buf[1024],bufa[1024];
	char *eFileA,*file;
	/*
	  Read inputfile
	*/
	file=offsets->file;
	datFile=appendSuffix(file,".dat",buf);
	eFileA=appendSuffix(file,".sa",bufa);
	/*
	 read offset param
	*/
	readOffsetParams(datFile,offsets);
	/*
	setup buffers
	*/
	initOffsetBuffers(offsets,AZONLY);
	/*
	  Read files
	*/
	readOffsetFile(offsets->da ,offsets->nr,offsets->na,offsets->file);
	readOffsetFile(offsets->sa ,offsets->nr,offsets->na,eFileA);
}

/*
   Read range  offsets (for now no sigma)
 */
void readRangeOffsets( Offsets *offsets) 
{
	char *datFile, buf[1024],bufa[1024];
	char *eFileR,*file;
	/*
	  Read inputfile
	*/
	file=offsets->rFile;
	datFile=appendSuffix(file,".dat",buf);
	/*
	 read offset param
	*/
	readOffsetParams(datFile,offsets);
	/*
	setup buffers
	*/
	initOffsetBuffers(offsets,RGONLY);
	/*
	  Read files
	*/
	readOffsetFile(offsets->dr ,offsets->nr,offsets->na,offsets->rFile);
}	

/*
   Read azimuth offsets
 */
void readAzimuthOffsets( Offsets *offsets) 
{
	char *datFile, buf[1024],bufa[1024];
	char *eFileA,*file;
	/*
	  Read inputfile
	*/
	file=offsets->file;
	datFile=appendSuffix(file,".dat",buf);
	/*
	 read offset param
	*/
	readOffsetParams(datFile,offsets);
	/*
	setup buffers
	*/
	initOffsetBuffers(offsets,AZONLY);
	/*
	  Read files
	*/
	readOffsetFile(offsets->da ,offsets->nr,offsets->na,offsets->file);
}

/*
  Input phase or power image for geocode
*/
void getMosaicInputImage(inputImageStructure *inputImage)
{   
	FILE *fp;
	float **fimage;
	float *imageLine;
	int i,j;
	/*
	  Open image
	*/    
	fprintf(stderr,"Reading %s --- ",inputImage->file);
	fimage=(float **)inputImage->image;
	if(strstr(inputImage->file,"nophase") == NULL) {
		fp = fopen(inputImage->file,"r");
		if(fp == NULL) error("*** getPhaseOrPowerImage: Error opening %s ***\n",
				     inputImage->file);
	} else fp = NULL;
	imageLine = fimage[0];

	if(fp != NULL) {
		freadBS(imageLine,sizeof(float),inputImage->rangeSize*inputImage->azimuthSize,fp,FLOAT32FLAG);
	} else { /* nophase case */
		for(i=0; i < inputImage->azimuthSize; i++) 
			for(j=0; j < inputImage->rangeSize; j++) fimage[i][j]=-2.0e9;
	}

	fprintf(stderr,"completed \n");
	if(fp != NULL) fclose(fp);
	return;
}
