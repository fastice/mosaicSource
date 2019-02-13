#include "stdio.h"
#include"string.h"
#include "math.h"
#include "stdlib.h"
#include "source/common/common.h"
#include "mosaic3d.h"

float dateWeight(double jd1,double jd2,double jdRange1, double jdRange2)
{
	int caseType;
	float weight;
	float percentCover;
	float dT,jdMid,jdMidRange,tDiff;
	/*
	  Time range check. This check makes sure that the center of the observation not to far skewed from center of
	  output window. The max allowed is what would be the case if the obs window and the output window were the same
	  size. This helps estimage avoid getting too skewed in time.
	*/
	dT=jdRange2-jdRange1; /* time range */
	jdMid=0.5*(jd1+jd2);
	jdMidRange=0.5*(jdRange1+jdRange2);
	tDiff=fabs(jdMid-jdMidRange);
	if(tDiff > dT) return(0.0);
	/* 
	   no real date check, so return 1.0
	*/
	if( jdRange2 > HIGHJD && jdRange1 <= LOWJD )   return(1.0);
	percentCover=(jdRange2-jdRange1)/(jd2-jd1);
	/* 
	   This is a resolution check the obs window can be larger than 60% of the output window - 
	   (e.g., for 30 day allows 48 day but not 64) 
	*/	
	if(percentCover < 0.6) return(0.0);
	/*
	  Now do individual cases
	  () -> (jd1,jd2)
	  [] -> [jdRange1,jdRange2]
	*/
	/*	Case 1 [ () ]  : image in output window	*/
	if( jd1 >= jdRange1 && jd2 <= jdRange2) {
		weight = 1.0;
		caseType=1; /* Error diagnostic */
		/* Case 2 [ ( ] )  : image skewed right	*/
	} else if( jd1 >= jdRange1 && jd1 <= jdRange2  && jd2 >= jdRange2) {
		weight=(jdRange2-jd1)/(jd2-jd1);
		if(percentCover < 1.0) weight *= percentCover; /* further deweight low res data */
		caseType=2; /* Error diagnostic */		
		/* Case 3 ( [ ) ]  : image skewed left	*/		
	} else if( jd1 <= jdRange1 && jd2 >= jdRange1  && jd2 <= jdRange2) {
		weight=(jd2-jdRange1)/(jd2-jd1);
		if(percentCover < 1.0) weight *= percentCover; /* further deweight low res data */
		caseType=3; /* Error diagnostic */		
		/* Case 4 ( [ ] )  : output window	 in image	*/			       
	} else if( jd1 <= jdRange1 && jd2 >= jdRange1 ) {
		weight=percentCover;
		caseType=4; /* Error diagnostic */	       
	} else return(0.0);
  	if(weight < 0 || weight > 1) {error("Error computing weights - case %i %f image %lf %lf range %lf %lf\n",caseType,weight,jd1,jd2,jdRange1,jdRange2);}
	/*fprintf(stderr,"weight -- %f\n",weight);*/
	return( weight);
}	

static void readTideCorrections( char *geodatFile, inputImageStructure *inputImage, FILE *fpLog,int i ) {
	FILE *fp;
	int lineCount,eod;
	char line[1024],tideFile[1024],*tmpS;
	
	/* Form name by looking for tide.correction in the geodat directory */
	strcpy(tideFile,geodatFile);
	tmpS=strrchr(tideFile,'/'); 
	tmpS[0]='\0';
	strcat(tideFile,"/tide.correction");
	fp=fopen(tideFile,"r");
	lineCount=0;
	if(fp != NULL) {
		lineCount=getDataString(fp,lineCount,line,&eod);
		sscanf(line,"%lf",&(inputImage->tideCorrection));
		fprintf(stderr,"Tide correction : %s\n",tideFile);              
		fprintf(fpLog,"; Tide correction (data take %i) : %s\n",i,tideFile);              
		fclose(fp);
	} else {
		fprintf(fpLog,"; Tide correction (data take %i) : %s\n",i,"none");              
		inputImage->tideCorrection=0.0;
	}
	fprintf(stderr,"Tide correction value: %f\n",inputImage->tideCorrection);
	fprintf(fpLog,"; Tide correction value : %f\n",inputImage->tideCorrection);
	/*
	  Read spatially varying tide correction if it exists in a "tide.difference" file
	*/	
	strcpy(tideFile,geodatFile);
	tmpS=strrchr(tideFile,'/');	tmpS[0]='\0';
	strcat(tideFile,"/tide.difference");
	fp=fopen(tideFile,"r");
	lineCount=0;
	if(fp != NULL) {
		inputImage->tideDiffFlag =TRUE;
		readXYDEM(tideFile,&(inputImage->tideDiff));
		fprintf(fpLog,"; Tide difference (data take %i) : %s\n",i,tideFile);              
		fclose(fp);
	} else {
		inputImage->tideDiffFlag = FALSE;
		fprintf(fpLog,"; Tide difference (data take %i) : %s\n",i,"none");              
	}
}

/*
  Add a new image and set of paramets to the list, saving the head if need be. Also save max dimensions. 
*/
static void addToList(char *phaseFile,vhParams *dumParams, inputImageStructure *inputImage ,inputImageStructure **ADList,
		      inputImageStructure **ADListHead, vhParams **ADParams, vhParams **ADParamsHead ,  int *maxR,int *maxA,int *nAD ) {

	/*Start of list */
	fprintf(stderr,"1\n");	
	if(*ADList== NULL) {
		*ADListHead=inputImage;
		*ADList = inputImage;
		*ADParams=dumParams;   *ADParamsHead=dumParams;
	} else {
		(*ADList)->next=inputImage; (*ADList)=(*ADList)->next;
		(*ADParams)->next=dumParams;    (*ADParams)=(*ADParams)->next;
	}
	(*ADList)->next=NULL;
	(*ADParams)->next=NULL;		
	(*nAD)++;
	if( strstr(phaseFile,"nophase") == NULL) { 
		*maxR=max(*maxR,inputImage->rangeSize);
		*maxA=max(*maxA,inputImage->azimuthSize);
		fprintf(stderr,"maxR,A %i %i\n",*maxR,*maxA);
	}

}

/*
   Malloc the rowpointers for the image and point them to the shared buff.
 */
static void setupADImageBuffers(inputImageStructure **images,float *buf ) {
	inputImageStructure *image;
	int i;
	for(image=*images; image != NULL; image=image->next) {
		/* Malloc row pointers */
		image->image= (void **) malloc(image->azimuthSize * sizeof(float *));
		/* Set pointers into mem pool */
		for(i=0; i < image->azimuthSize; i++) {
			image->image[i]= &(buf[i*image->rangeSize]);
		}
	}
}

/*
  set up 3d
*/
void  setup3D(int nFiles,char **phaseFiles, char **geodatFiles,  char **baselineFiles,char **offsetFiles,char **azParamsFiles,
	      char **rOffsetFiles,char **rParamsFiles, float *nDays, float *weights, int *crossFlags, inputImageStructure **ascImages, inputImageStructure **descImages,
	      vhParams **ascParams, vhParams **descParams,  int *nAsc, int *nDesc,  int offsetFlag, int rOffsetFlag, int threeDOffFlag, FILE *fpLog,
	      outputImageStructure *outputImage)
{
	extern float *AImageBuffer, *DImageBuffer;
	FILE *fp;	
	vhParams *dumParams, *vhD, *vhA;
	inputImageStructure *inputImage,*desc,*asc;
	double julDay1,julDay2;
	int maxRa,maxAa,maxRd,maxAd;
	float weight, tmpWeight;
	int i;
	fprintf(fpLog,";\n; Entering setUp3D\n;\n");
	/*	  Initializations for maxes	*/
	maxRa=0; maxAa=0; maxRd=0; maxAd=0;
	/*	  Malloc input images	*/
	inputImage = (inputImageStructure *) malloc( sizeof(inputImageStructure)*nFiles);
	/*
	  Parse input file and input data for each image.
	*/
	asc=NULL; desc=NULL;
	*ascImages=NULL; *descImages=NULL;
	*nAsc=0; *nDesc=0;
	for(i=0; i < nFiles; i++) {
		inputImage[i].stateFlag = TRUE;
		inputImage[i].llInit=FALSE;
		/*  Parse the Geodat file	*/
		parseInputFile(geodatFiles[i], &(inputImage[i]));
		julDay1=inputImage[i].julDay;
		julDay2=inputImage[i].julDay+nDays[i];
		/*    test within date range	*/
		inputImage[i].used=FALSE;
		weight=weights[i];
		/*	Compute weight flag based on degree of overlap	*/
		tmpWeight=dateWeight(julDay1,julDay2, outputImage->jd1,outputImage->jd2);
		/* If  full (0.9999) overlap, so use image		*/
		if(tmpWeight > 0.99999 ) inputImage[i].used=TRUE;
		/* otherwise only use image if the timeOverlapFlag set, in which case cascade with other weights */
		else if( outputImage->timeOverlapFlag==TRUE )  {
			weight=weights[i]*tmpWeight;
			inputImage[i].used=TRUE;
		}
		weights[i]=weight;
		inputImage[i].weight = weights[i];
		inputImage[i].crossFlag=crossFlags[i];
		/*
		  Proceed if images is used
		*/
		if( inputImage[i].used == TRUE ) {
			fprintf(stderr,"TTTTTT time weight = %f\n",weight);
			inputImage[i].file=phaseFiles[i];
			dumParams=(vhParams *)malloc(sizeof(vhParams));
			/*  Read tide correction files where they exist */
			readTideCorrections(geodatFiles[i], &(inputImage[i]),fpLog,i );
			/*	  Get baseline  */
			getBaseline(baselineFiles[i],dumParams);
			/*	  Get time info	*/
			fprintf(stderr,"nDays %f\n",nDays[i]);
			dumParams->nDays=nDays[i];
			if(nDays[i] < 0) error("Negative number of days");
			/*
			  Read offsets if necessary - offset flag obsolete - it should always be true
			*/       
			dumParams->offsetFlag=offsetFlag;
			dumParams->rOffsetFlag=rOffsetFlag;
			if(offsetFlag==TRUE || threeDOffFlag==TRUE) {
				dumParams->offsets.file=offsetFiles[i];
				dumParams->offsets.azParamsFile=azParamsFiles[i];
			}
			if(rOffsetFlag==TRUE || threeDOffFlag==TRUE) {
				dumParams->offsets.rParamsFile=rParamsFiles[i];
				dumParams->offsets.rFile=rOffsetFiles[i];
			} else {
				dumParams->offsets.rParamsFile=NULL;
				dumParams->offsets.rFile=NULL;
			}
			dumParams->offsets.bnS=NULL;
			dumParams->offsets.bpS=NULL;
			dumParams->offsets.azInit=FALSE;
			dumParams->offsets.rOffS=0.0;
			dumParams->offsets.deltaB=outputImage->deltaB;
			if(inputImage[i].passType == DESCENDING ) 
				addToList(phaseFiles[i],dumParams,  &(inputImage[i]) ,&desc,   descImages, &vhD, descParams,  &maxRd,&maxAd, nDesc );
			else 
				addToList(phaseFiles[i],dumParams,  &(inputImage[i]) ,&asc,   ascImages, &vhA, ascParams,  &maxRa,&maxAa, nAsc );
		} else {  geodatFiles[i]=NULL;} /* Case where not in time range */
	} 
	fprintf(stderr,"NAsc %i  NDesc %i\n",*nAsc,*nDesc);
	/*
	  Use the same pool of memory for each image. Allocate different 
	  set of pointers for each row.	  Ascending first.
	  Malloc two buffers here using based on asc/desc type. Later in 
	  make3dMosaic the buffers will be re-allocated as needed to hold
	  the two seperate images. 
	*/
	if(maxRa > 0 || maxRd > 0) {
		AImageBuffer= malloc(sizeof(float)*max(maxRa,maxRd)*max(maxAa,maxAd));
		DImageBuffer= malloc(sizeof(float)*max(maxRa,maxRd)*max(maxAa,maxAd));
	}
	/*
	  Setubuffers for the shared memory pool
	 */
	if(*nAsc > 0) setupADImageBuffers(ascImages,AImageBuffer );
	if(*nDesc > 0)  setupADImageBuffers(descImages,DImageBuffer );	

	fprintf(fpLog,";\n; Leaving setUp3D");
	return;
}

