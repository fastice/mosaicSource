#include "stdio.h"
#include"string.h"
#include "math.h"
#include <stdlib.h>
#include "mosaicSource/common/common.h"
#include "mosaic3d.h"

static void mallocFile(char phase[512],char ***phaseFiles, int i) {
	int j;
	if(strlen(phase) > 512) {
		fprintf(stderr,"%s\n",phase);
		error("file name to long %s",phase);
	}
	(*phaseFiles)[i] = (char *)malloc(strlen(phase)+1);
	for(j=0; j < strlen(phase); j++) 
		(*phaseFiles)[i][j] = phase[j];
	(*phaseFiles)[i][j] = '\0';
}
	 
/*
   Process input file for mosaicDEMs
*/
    void  get3DInputFile(char *inputFile,char ***phaseFiles,   char ***geodatFiles, char ***baselineFiles, 
                           outputImageStructure *outputImage,    float **nDays, float **weights, int *nFiles)
{
    FILE *fp;
    double xo,yo,xs,ys,deltaX,deltaY;
    float weight,nDay;
    int lineCount, eod, i;
    char phase[512],geodat[512],baseline[512];
    char line[1024];
/*
   Open file for input
*/
    fp = openInputFile(inputFile);
/*
   Input nr,na,nlr,nla
*/
    lineCount=getDataString(fp,lineCount,line,&eod);
    if( sscanf(line,"%lf%lf%lf%lf%lf%lf",&xo,&yo,&xs,&ys,&deltaX,&deltaY) != 6)
        error("%s  %i of %s",
            "get3dInputFile -- Missing image parameters at line:",
             lineCount,inputFile); 
     
    outputImage->xSize = (int)(xs/deltaX + 0.5);
    outputImage->ySize = (int)(ys/deltaY + 0.5);
    outputImage->deltaX = deltaX * KMTOM;
    outputImage->deltaY = deltaY * KMTOM;
    outputImage->originX = xo * KMTOM;
    outputImage->originY = yo * KMTOM;
/*
   Input ReMajor,ReMinor, Rc, phic, H
*/
    lineCount=getDataString(fp,lineCount,line,&eod);
    if( sscanf(line,"%i",nFiles) != 1)
        error("%s  %i  of %s",
            "get3dInputFile -- Missing geometric parameters at line:",
             lineCount,inputFile); 
/*
   Malloc space for arrays of filenames
*/
     *phaseFiles = (char **) malloc(sizeof(char *) * (*nFiles)); 
     *geodatFiles = (char **) malloc(sizeof(char *) * (*nFiles));
     *baselineFiles = (char **) malloc(sizeof(char *) * (*nFiles));
     *nDays= (float *)malloc(sizeof(float) * (*nFiles));
     *weights= (float *)malloc(sizeof(float) * (*nFiles));
/*
   Input files
*/
     for(i=0; i < *nFiles; i++) {
         lineCount=getDataString(fp,lineCount,line,&eod);
         if( sscanf(line,"%s %s %s %f %f",
            phase,geodat,baseline,&nDay,&weight) != 5) {
             weight=1.0;
             if( sscanf(line,"%s %s %s %f",
                 phase,geodat,baseline,&nDay) != 4)
                error("%s  %i  of %s",
                   "get3DinputFile -- Missing filename at line:",
                    lineCount,inputFile); 
         }
	 mallocFile( phase,phaseFiles,  i);
	 mallocFile( geodat,geodatFiles,  i);
	 mallocFile( baseline,baselineFiles,  i);
	 /*
          (*phaseFiles)[i] = (char *)malloc(strlen(phase)+1);
         for(j=0; j < strlen(phase); j++) 
             (*phaseFiles)[i][j] = phase[j];
         (*phaseFiles)[i][j] = '\0';

         (*geodatFiles)[i] = (char *)malloc(strlen(geodat)+1);
         for(j=0; j < strlen(geodat); j++) 
             (*geodatFiles)[i][j] = geodat[j];
         (*geodatFiles)[i][j] = '\0';

         (*baselineFiles)[i] = (char *)malloc(strlen(baseline)+1);
         for(j=0; j < strlen(baseline); j++) 
             (*baselineFiles)[i][j] = baseline[j];
         (*baselineFiles)[i][j] = '\0';
	 */
         (*weights)[i]=weight;
         (*nDays)[i] = nDay;

     }
     fclose(fp);
     return;
}

