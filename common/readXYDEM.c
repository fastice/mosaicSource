#include "stdio.h"
#include <stdlib.h>
#include "string.h"
#include "mosaicSource/common/common.h"
#include "math.h"
#include "clib/standard.h"  
#include <sys/types.h>

void readXYDEMGeoInfo(char *xyFile, xyDEM *xydem, int resetProjection);

void readXYDEMcrop(char *xyFile, xyDEM *xydem,float xmin,float xmax,float ymin,float ymax)
{  
	extern int HemiSphere;
	extern double Rotation;
	extern double SLat;	
	FILE *fp;
	float **image;
	int lineCount=0, eod;
	int i,j;
	float dum1,dum2;
	float xmx,ymx,ymn,xmn;
	long int x1,x2,y1,y2;
	xyDEM Tmp, *xyTemp;
	char *geodatFile;
	char line[256];
	off_t offset1;
	xyTemp=&Tmp; /* temp dem */
	fprintf(stderr,"READING *** cropped *** XYDEM\n");
	/*
	  Defaults
	*/
	xydem->rot=Rotation;
	xydem->hemisphere=HemiSphere;
	xydem->stdLat = 70.0;
	/*
	  Compute geodat file name
	*/
	geodatFile = (char *) malloc((size_t)256); 
	geodatFile[0] ='\0'; 
	geodatFile = strcat(geodatFile,xyFile);   /* filename */
	geodatFile = strcat(geodatFile,".geodat");
	fp = openInputFile(geodatFile);
	if(fp == NULL) 
		error("*** getXYDEM: Error opening %s ***\n",xyFile);
	/*
	  Read parameters for DEM
	*/
	lineCount=getDataString(fp,lineCount,line,&eod); /* Skip # 2 line */
	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%f %f\n",&dum1,&dum2);
	xyTemp->xSize=(int)dum1;
	xyTemp->ySize=(int)dum2;

	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%lf %lf\n",&(xyTemp->deltaX),&(xyTemp->deltaY));
	xyTemp->deltaX *= MTOKM;  xyTemp->deltaY *= MTOKM;

	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%lf %lf\n",&(xyTemp->x0),&(xyTemp->y0));
	fprintf(stderr,"%s\n",line);
	/*
	  Now compute parameters for cropped dem
	*/
	xmn=xyTemp->x0; /* lower left corner defines min extent of input dem */
	ymn=xyTemp->y0; 

	xmx=xmn + xyTemp->deltaX * (float)(xyTemp->xSize-1); /* upper right corner defines max extent */
	ymx=ymn + xyTemp->deltaY * (float)(xyTemp->ySize-1);

	x1=(long)( (max(xmn,xmin) -xyTemp->x0)/xyTemp->deltaX + 0.5);
	y1=(long)( (max(ymn,ymin) -xyTemp->y0)/xyTemp->deltaY + 0.5);

	x2=x1+(long)((min(xmx,xmax) - max(xmn,xmin))/xyTemp->deltaX + 0.5)+1;
	y2=y1+(long)((min(ymx,ymax) - max(ymn,ymin))/xyTemp->deltaY + 0.5)+1;
	xydem->xSize=x2-x1+1;
	xydem->ySize=y2-y1+1;
	xydem->x0=xyTemp->x0 + x1 * xyTemp->deltaX;
	xydem->y0=xyTemp->y0 + y1 * xyTemp->deltaY;
	xydem->deltaX=xyTemp->deltaX;
	xydem->deltaY=xyTemp->deltaY;

	fprintf(stderr,"input dem limits %f %f %f %f\n",xmn,xmx,ymn,ymx);
	fprintf(stderr,"requested dem limits %f %f %f %f\n",xmin,xmax,ymin,ymax);
	fprintf(stderr,"output dem limits %f %f %f %f\n",max(xmn,xmin),min(xmx,xmax),max(ymn,ymin),min(ymx,ymax) );
	fprintf(stderr,"pixel limits %li %li %li %li\n",x1,x2,y1,y2);

	lineCount=getDataString(fp,lineCount,line,&eod);
	fprintf(stderr,"%s\n",line);
	/* Hemisphere */
	if(eod != TRUE) {
		fprintf(stderr,"Hemispher\n");
		if(strstr(line,"SOUTH")) {
			HemiSphere=SOUTH;  
		} else HemiSphere=NORTH;
		lineCount=getDataString(fp,lineCount,line,&eod);
	}
	/* Rotation */

	if(eod != TRUE) {
		fprintf(stderr,"Rot\n");
		sscanf(line,"%lf",&(xydem->rot)); 
		lineCount=getDataString(fp,lineCount,line,&eod);

	}
	/* Std lat */
	if(eod != TRUE) {
		fprintf(stderr,"Stdlat\n");
		sscanf(line,"%lf",&(xydem->stdLat)); 
	}
	/*
	  Make sure stdlat is -
	*/
	/*    if(HemiSphere==SOUTH) xydem->stdLat=-fabs(xydem->stdLat);*/
	fprintf(stderr,"Dem Size %i %i\n",xydem->xSize,xydem->ySize);
	fprintf(stderr,"Dem Res %f %f\n",xydem->deltaX,xydem->deltaY);
	fprintf(stderr,"Dem Origin %f %f\n",xydem->x0,xydem->y0);
	if(HemiSphere == SOUTH) {
		fprintf(stderr,"SOUTHER HEMISPHERE PROJECTION\n");
		xydem->stdLat=71;
	}
	else  fprintf(stderr,"NORTHERN HEMISPHERE PROJECTION\n");
	fprintf(stderr,"Rotation = %f\n",xydem->rot);
	fprintf(stderr,"Standard Lat = %f\n",xydem->stdLat);
	fclose(fp);
	/*
	  Open input file
	*/
	fp = fopen(xyFile,"r");
	if(fp == NULL) 
		error("*** getXYDEM: Error opening %s ***\n",xyFile);
	/*
	  Read input file
	*/ 
	image = (float **)mallocImage(xydem->ySize,xydem->xSize);
	j=0;
	for(i=y1; i <= y2; i++) {


		offset1 = ((off_t)xyTemp->xSize * (off_t) i + (off_t)x1 )* (off_t)sizeof(float);
		/*
		  fprintf(stderr,"%i %i\n",xydem->xSize,xydem->ySize);
		  fprintf(stderr,"offset - %i %i %lli %i %i \n",j,i,offset1, xydem->xSize,sizeof(off_t));
		  fprintf(stderr,"--- %i\n",sizeof(off_t));
		*/
		/* Changed 11/4/2016 from fseeko64 */
		fseeko(fp,offset1,SEEK_SET);
		freadBS(image[j],sizeof(float),xydem->xSize,fp,FLOAT32FLAG);
		j++;
	}

	xydem->z = image;
	fprintf(stderr,"end readXYDEMcrop\n");
	fclose(fp);
} 

void readXYDEMGeoInfo(char *xyFile, xyDEM *xydem, int resetProjection) {
	extern int HemiSphere;
	extern double Rotation;
	extern double SLat;
	FILE *fp;
	float **image;
	int lineCount=0, eod;
	int i;
	float dum1,dum2;
	char *geodatFile;
	char line[256];
	/*
	  Defaults
	*/
	xydem->rot=Rotation;
	xydem->hemisphere=HemiSphere;
	xydem->stdLat = 70.0;
	/*
	  Compute geodat file name
	*/
	geodatFile = (char *) malloc((size_t)256); 
	geodatFile[0] ='\0'; 
	geodatFile = strcat(geodatFile, xyFile);   /* filename */
	geodatFile = strcat(geodatFile, ".geodat");
	fp = openInputFile(geodatFile);
	if(fp == NULL) 
		error("*** getXYDEM: Error opening %s ***\n",xyFile);
	/*
	  Read parameters
	*/
	lineCount=getDataString(fp,lineCount,line,&eod); /* Skip # 2 line */
	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%f %f\n",&dum1,&dum2);
	xydem->xSize = (int)dum1;
	xydem->ySize = (int)dum2;

	lineCount = getDataString(fp, lineCount, line, &eod);
	sscanf(line,"%lf %lf\n",&(xydem->deltaX), &(xydem->deltaY));
	xydem->deltaX *= MTOKM;  xydem->deltaY *= MTOKM;
	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%lf %lf\n",&(xydem->x0),&(xydem->y0));
	fprintf(stderr,"%s\n",line);
	lineCount=getDataString(fp,lineCount,line,&eod);
	fprintf(stderr,"%s\n",line);
	/* Hemisphere */
	if(eod != TRUE) {
		fprintf(stderr,"Hemispher\n");
		if(strstr(line,"SOUTH")) {
			HemiSphere=SOUTH;  
		} else HemiSphere=NORTH;
		lineCount=getDataString(fp,lineCount,line,&eod);
	}
	/* Rotation */
	if(eod != TRUE) {
		fprintf(stderr,"Rot\n");
		sscanf(line,"%lf",&(xydem->rot)); 
		lineCount=getDataString(fp,lineCount,line,&eod);
	}
	/* Std lat */
	if(eod != TRUE) {
		fprintf(stderr,"Stdlat\n");
		sscanf(line,"%lf",&(xydem->stdLat)); 
	}
	/*
	  Make sure stdlat is -
	*/
	fprintf(stderr,"Dem Size %i %i\n",xydem->xSize,xydem->ySize);
	fprintf(stderr,"Dem Res %f %f\n",xydem->deltaX,xydem->deltaY);
	fprintf(stderr,"Dem Origin %f %f\n",xydem->x0,xydem->y0);
	if(HemiSphere == SOUTH) {
		fprintf(stderr,"SOUTHER HEMISPHERE PROJECTION\n");
	} else  fprintf(stderr,"NORTHERN HEMISPHERE PROJECTION\n");
	fprintf(stderr,"Rotation = %f\n",xydem->rot);
	fprintf(stderr,"Standard Lat = %f\n",xydem->stdLat);
	fclose(fp);
	/* Force projection to match dem */
	if(resetProjection == TRUE) {
		Rotation = xydem->rot;
		SLat = xydem->stdLat;
	}
}


void readXYDEM(char *xyFile, xyDEM *xydem)
{  
	extern int HemiSphere;
	extern double Rotation;
 
	FILE *fp;
	float **image;
	int lineCount=0, eod;
	int i;
	float dum1,dum2;
	char *geodatFile;
	char line[256];
	fprintf(stderr,"READING XYDEM\n");
	readXYDEMGeoInfo(xyFile, xydem, FALSE);
	/*
	  Open input file
	*/
	fp = fopen(xyFile,"r");
	if(fp == NULL) 
		error("*** getXYDEM: Error opening %s ***\n",xyFile);
	/*
	  Read input file
	*/ 
	image = (float **)mallocImage(xydem->ySize, xydem->xSize);
	for(i=0; i < xydem->ySize; i++) {
		freadBS(image[i],sizeof(float),xydem->xSize, fp, FLOAT32FLAG);        
	}

	xydem->z = image;
	fclose(fp);
} 



void readXYVel(xyVEL *xyvel, char *velFile)
{
	FILE *fp;
	float **image;
	int lineCount=0, eod;
	int i;
	char *geodatFile;
	char line[256];
	char vxfile[1024],vyfile[1024];
	/*
	  Compute geodat file name
	*/
	geodatFile = (char *) malloc((size_t)256); 
	geodatFile[0] ='\0'; 
	geodatFile = strcat(geodatFile,velFile);   /* filename */
	geodatFile = strcat(geodatFile,".vx.geodat");
	fp = openInputFile(geodatFile);
	if(fp == NULL) 
		error("*** getXYVEL: Error opening %s ***\n",velFile);
	/*
	  Read parameters
	*/
	lineCount=getDataString(fp,lineCount,line,&eod); /* Skip # 2 line */
	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%i %i\n",&(xyvel->xSize),&(xyvel->ySize));

	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%lf %lf\n",&(xyvel->deltaX),&(xyvel->deltaY));
	xyvel->deltaX *= MTOKM;  xyvel->deltaY *= MTOKM;

	lineCount=getDataString(fp,lineCount,line,&eod);
	sscanf(line,"%lf %lf\n",&(xyvel->x0),&(xyvel->y0));
	fclose(fp);
	/*
	  Open input file
	*/
	vxfile[0]='\0';
	vyfile[0]='\0';
	strcat(vxfile,velFile);
	strcat(vxfile,".vx");
	strcat(vyfile,velFile);
	strcat(vyfile,".vy");
	/*
	  Read vx file
	*/
	fp = fopen(vxfile,"r");
	if(fp == NULL) 
		error("*** getXYVEL: Error opening %s ***\n",velFile);
	image = (float **)mallocImage(xyvel->ySize,xyvel->xSize);
	for(i=0; i < xyvel->ySize; i++) 
		freadBS(image[i],sizeof(float),xyvel->xSize,fp,FLOAT32FLAG);
	xyvel->vx = image;
	fclose(fp);
	/*
	  Read vy file
	*/
	fp = fopen(vyfile,"r");
	if(fp == NULL) 
		error("*** getXYVEL: Error opening %s ***\n",velFile);
	image = (float **)mallocImage(xyvel->ySize,xyvel->xSize);
	for(i=0; i < xyvel->ySize; i++) 
		freadBS(image[i],sizeof(float),xyvel->xSize,fp,FLOAT32FLAG);
	xyvel->vy = image;
	fclose(fp);
} 
