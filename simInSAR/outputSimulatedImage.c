#include "stdio.h"
#include"string.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"

/*
  Output simulated image. Writes two files one for image, and xxx.simdat 
  with image header info
*/
void outputSimulatedImage(sceneStructure scene,char *outputFile,        char *demFile, char *displacementFile)
{   
	FILE *imageFP, *imageDatFP;
	int i,j;
	int pSize;
	char *outputFileDat, *tmp;
	char *buf,buf1[1024];
	/*
	  Open image outputfile
	*/
	if(scene.toLLFlag==FALSE) {
		imageFP = fopen(outputFile,"w");
		if(imageFP == NULL) 
			error("*** outputSimulatedImage: Error opening %s ***\n",outputFile);
		/*
		  Output image data
		*/
		pSize = sizeof(float);
		if(scene.maskFlag == TRUE) buf=(char *) malloc(scene.I.rangeSize*sizeof(char));
		for(i=0; i < scene.I.azimuthSize; i++) {
			if(scene.maskFlag == FALSE) fwriteBS(scene.image[i],scene.I.rangeSize,pSize,imageFP,FLOAT32FLAG);
			else {
				for(j=0; j < scene.I.rangeSize; j++)  buf[j]= (char)scene.image[i][j];
				fwriteBS(buf,scene.I.rangeSize,sizeof(char),imageFP,BYTEFLAG);
			}
		}
		/*  free(buf);*/
		fclose(imageFP);
	} else {
		buf1[0]='\0';
		strcpy(buf1,outputFile);
		strcat(buf1,".lat");
		fprintf(stderr,"writing %s\n",buf1);
		imageFP = fopen(buf1,"w");
		for(i=0; i < scene.aSize; i++) {
			fwriteBS(scene.latImage[i],scene.rSize,sizeof(double),imageFP,FLOAT64FLAG);
		}
		fclose(imageFP);
		    
		buf1[0]='\0';
		strcpy(buf1,outputFile);
		strcat(buf1,".lon");
		fprintf(stderr,"writing %s\n",buf1);
		imageFP = fopen(buf1,"w");		    
		for(i=0; i < scene.aSize; i++) {
			fwriteBS(scene.lonImage[i],scene.rSize,sizeof(double),imageFP,FLOAT64FLAG);		    
		}
		if(scene.maskFlag == TRUE) {
			buf1[0]='\0';
			strcpy(buf1,outputFile);
			strcat(buf1,".mask");
			fprintf(stderr,"writing %s %i %i\n",buf1,scene.rSize,scene.aSize );
			imageFP = fopen(buf1,"w");
			if(imageFP == NULL) 
				error("*** outputSimulatedImage: Error opening %s ***\n",outputFile);
			/*
			  Output image data
			*/
			buf=(char *) malloc(scene.rSize*sizeof(char));
			for(i=0; i < scene.aSize; i++) {
				for(j=0; j < scene.rSize; j++)  buf[j]= (char)scene.image[i][j];
				fwriteBS(buf,scene.rSize,sizeof(char),imageFP,BYTEFLAG);
			}
			/*  free(buf);*/
			fclose(imageFP);
		}
	}
	/*
	  Form header filename by adding .simdat suffix to outputFile
	*/
	tmp = (char *) malloc(1024); 
	tmp[0] =  '\0'; 
	tmp = strcat(tmp,outputFile);   /* Compute input filename */
	tmp = strcat(tmp,".simdat");
	outputFileDat = (char *) malloc(strlen(tmp)+1);
	outputFileDat = strcpy(outputFileDat,tmp);
	free(tmp);
	/*
	  Open image outputfile
	*/
	imageDatFP = fopen(outputFileDat,"w");
	if(imageDatFP == NULL) error("*** outputSimulatedImage: Error opening %s ***\n",outputFileDat);
	/*
	  Output image header info.
	*/
	fprintf(imageDatFP,"# 2\n;\n;  Image size (pixels) nx ny \n;\n");
	fprintf(imageDatFP,"%i  %i\n",scene.I.rangeSize,scene.I.azimuthSize);
	fprintf(imageDatFP,"; Baseline start/end\n");
	fprintf(imageDatFP,"%f  %f\n&\n",scene.bnStart,scene.bnEnd);
	fprintf(imageDatFP,"; demFile :        %s\n",demFile);
	fprintf(imageDatFP,"; displacemtnFile: %s\n",displacementFile); 
	fclose(imageDatFP);
	return;
}
