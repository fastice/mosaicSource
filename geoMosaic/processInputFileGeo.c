#include "stdio.h"
#include"string.h"
#include "math.h"
#include "mosaicSource/common/common.h"
/*#include "ers1/GeoCodeDEM_p/geocodedem.h"*/
#include "geomosaic.h"
#include <stdlib.h>

/*
  Process input file for mosaicDEMs
*/
void  processInputFileGeo(char *inputFile,char ***insarDEMFiles, char ***demInputFiles,  
			  outputImageStructure *outputImage, int *nDEMs, float **weights, char ***antPatFiles)
{
	FILE *fp;
	double xo,yo,xs,ys,deltaX,deltaY;
	int lineCount, eod;
	int i,j;
	float w;
	char insarDEM[1024],demInput[1024],antPat[1024] ;
	char line[512];
	/*
	  Open file for input
	*/
	fp = openInputFile(inputFile);
	/*
	  Input nr,na,nlr,nla
	*/
	lineCount=getDataString(fp,lineCount,line,&eod);
	if( sscanf(line,"%lf%lf%lf%lf%lf%lf",&xo,&yo,&xs,&ys,&deltaX,&deltaY) != 6)
		error("%s  %i of %s",  "processInputFile -- Missing image parameters at line:",  lineCount,inputFile); 
     
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
	if( sscanf(line,"%i",nDEMs) != 1)
		error("%s  %i  of %s",
		      "process InputFile -- Missing geometric parameters at line:",
		      lineCount,inputFile); 
	/*
	  Malloc space for arrays of filenames
	*/
	*insarDEMFiles = (char **) malloc((size_t)(sizeof(char *) * (*nDEMs)));
	*demInputFiles = (char **) malloc((size_t)(sizeof(char *) * (*nDEMs)));
	*antPatFiles =   (char **) malloc((size_t)(sizeof(char *) * (*nDEMs)));
	*weights = (float *) malloc((size_t)(sizeof(float) * (*nDEMs)));
	/*
	  Input files
	*/
	for(i=0; i < *nDEMs; i++) {
		lineCount=getDataString(fp,lineCount,line,&eod);

		if( sscanf(line,"%s %s %f %s",insarDEM,demInput,&w,antPat) != 4) { 
			antPat[0]='\0';
			if( sscanf(line,"%s %s %f",insarDEM,demInput,&w) != 3) { 
				w=1.0; 
				if( sscanf(line,"%s %s",insarDEM,demInput) != 2)
					error("%s  %i  of %s",    "process InputFile -- Missing filename at line:",  lineCount,inputFile); 
			}
		}
		(*insarDEMFiles)[i] = (char *)malloc((size_t)(strlen(insarDEM)+1));
		for(j=0; j < strlen(insarDEM); j++) 
			(*insarDEMFiles)[i][j] = insarDEM[j];
		(*insarDEMFiles)[i][j] = '\0';

		(*demInputFiles)[i] = (char *)malloc((size_t)(strlen(demInput)+1));
		for(j=0; j < strlen(demInput); j++) 
			(*demInputFiles)[i][j] = demInput[j];
		(*demInputFiles)[i][j] = '\0';
		if(strlen(antPat) > 0 && ( (strstr(antPat,"poly") != NULL) || (strstr(antPat,"alos") != NULL)))  {
			(*antPatFiles)[i] = (char *)malloc((size_t)(strlen(antPat)+1));
			for(j=0; j < strlen(antPat); j++) 
				(*antPatFiles)[i][j] = antPat[j];
			(*antPatFiles)[i][j] = '\0';
		} else (*antPatFiles)[i] = NULL;

		(*weights)[i]=w;
	}
	return;
}

