#include "stdio.h"
#include"string.h"
#include "ctype.h"
#include "math.h"
#include <stdlib.h>
#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"

/*
  parse scene input file for siminsar.
*/
void parseSceneFile(char *sceneFile, sceneStructure *scene)
{
	FILE *fp;
	int i;
	char *datFile, buf[1024],buf1[1024],*errorFile;
	int rO,aO,nr,na;
	float deltaR,deltaA;
	int lineCount,eod,nRead;
	char line[1024];
	parseInputFile(sceneFile,&(scene->I) );
	/*
	  Compute baseline increment size
	*/
	scene->bnStep =(scene->bnEnd - scene->bnStart)/   (double)(scene->I.azimuthSize-1.0);
	scene->bpStep =(scene->bpEnd - scene->bpStart)/   (double)(scene->I.azimuthSize-1.0);

	if( scene->toLLFlag == TRUE) {
		/*
		  Read inputfile
		*/
		datFile=scene->llInput;
		fp = openInputFile(datFile);
		lineCount=getDataString(fp,lineCount,line,&eod);
		nRead=sscanf(line,"%i%i%i%i%f%f",&rO,&aO,&nr,&na,&deltaR,&deltaA);
		if(nRead < 6)
			error("%s  %i of %s","readOffsets -- Missing image parameters at line:",  lineCount,datFile);
		scene->aSize=na;
		scene->rSize=nr;
		scene->aO=aO/scene->I.nAzimuthLooks ;
		scene->rO=rO/scene->I.nRangeLooks ;

		scene->dR=deltaR/scene->I.nRangeLooks;
		scene->dA=deltaA/scene->I.nAzimuthLooks;
		scene->latImage = (double **)  malloc(scene->aSize * sizeof(double *) );
		scene->lonImage = (double **) malloc(scene->aSize * sizeof(double *) );

		for(i=0; i < scene->aSize; i++) {
			scene->lonImage[i] = (double *)malloc( scene->rSize * sizeof(double) );		    
			scene->latImage[i] =  (double *)malloc( scene->rSize * sizeof(double) );
		}
		if(scene->maskFlag==TRUE) {
			scene->image = (float **) malloc(scene->aSize * sizeof(float *) );
			for(i=0; i < scene->aSize; i++) {
				scene->image[i] = (float *)malloc( scene->rSize * sizeof(float) );
			}
		}
		fprintf(stderr,"reading %s\n",scene->llInput);
	} else {
		scene->aSize=scene->I.azimuthSize;
		scene->rSize=scene->I.rangeSize;
		scene->aO=0;
		scene->rO=0;
		scene->dR=1;
		scene->dA=1;
		/* Save lat/lon as */
		if(scene->saveLLFlag == TRUE) {
			scene->latImage = (double **)  malloc(scene->aSize * sizeof(double *) );
			scene->lonImage = (double **) malloc(scene->aSize * sizeof(double *) );		
			for(i=0; i < scene->aSize; i++) {
				scene->lonImage[i] = (double *)malloc( scene->rSize * sizeof(double) );		    
				scene->latImage[i] =  (double *)malloc( scene->rSize * sizeof(double) );
			}
		}
		/*
		  Init space for image 
		*/		
		scene->image = (float **) malloc(scene->aSize * sizeof(float *) );		
		for(i=0; i < scene->aSize; i++) {
			scene->image[i] = (float *)malloc( scene->rSize * sizeof(float) );
		}
	}
	return;
}
