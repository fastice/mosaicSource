#include "stdio.h"
#include"string.h"
#include <math.h>
#include "cRecipes/nrutil.h"
#include <stdlib.h>
#include "mosaicSource/common/common.h"

float *offBuf1,*offBuf2,**offBuf1a,**offBuf1b;
/*
************************ Geocode image InSAR DEM. **************************
*/

/*
  Find region where image exists
*/
void errorsToXY(double er,double ea, double *ex, double *ey, double xyAngle, double hAngle);

static void setBuffer(inputImageStructure *inputImage,float *buf)
{
	int i;
	for(i=0; i < inputImage->azimuthSize; i++) {
		inputImage->image[i]= &(buf[i*inputImage->rangeSize]);
	}
}

/*
  Make mosaic of insar and other dems.
*/
void makeVhMosaic(inputImageStructure *images,vhParams *params,  outputImageStructure *outputImage,float fl)
{   
	extern int HemiSphere;
	extern double Rotation;
	double lat,lon;
	vhParams *currentParams;
	conversionDataStructure *cP;
	inputImageStructure *currentImage;
	double ReH, Re;
	float **vXimage,**vYimage,**vZimage;
	float **scaleX,**scaleY,**scaleZ;
	double tdTemp;
	double twok;
	double dzda,dzdr;      /* Slopes in azimuth and range */
	double range, azimuth; /* range,azimuth coords */
	double Range;          /* Slant Range */
	double phase;          /* Interpolated phase value */
	double x,y,zSp,zWGS84;
	double thetaC, thetaD;  /* Center look angle and deviation from center */
	double ReHfixed,thetaCfixedReH;
	double theta,psi, cotanpsi;           /* Look and inc angle */
	double phiZ;            /* Phase due to topography */
	double scalePhase;      /* Scale factor from phase to range diff */
	double delta;           /* Range difference */
	double hAngle;          /* Heading angle */
	double va,vr;           /* Components of vel in azimuth and range */
	double xyAngle;         /* Angle from north */
	double er,ea,ex,ey;           /* Relative errors */
	double sigmaR,sigmaA;
	double phaseError;
	double tideError; /* Added 05/21/03 */
	double tCenter,tOffCenter,deltaOffCenter;
	float azSLPixSize, dum;
	ShelfMask *shelfMask;    /* Mask with shelf and grounding zone */
	unsigned char sMask;
	int iMin,iMax,jMin,jMax;
	double scX,scY;
	float **vxTmp,**vyTmp,**vzTmp,**fScale,**sxTmp,**syTmp;
	float **errorX,**errorY;
	float da;
	double vx,vy;
	int i,j;
	int count;                /* Current image counter - info only */
	int validData;

	fprintf(stderr,"MOSAICVH\n");
	fprintf(outputImage->fpLog,";\n; Entering makeVhMosaic(.c)\n");
	shelfMask = outputImage->shelfMask;
	/*
	  Pointers to output images
	*/
	setupBuffers(outputImage, &vXimage, &vYimage, &vZimage,&scaleX,&scaleY, &scaleZ,
		     &vxTmp, &vyTmp, &vzTmp,&sxTmp,&syTmp, &fScale,&errorX,&errorY);	
	/*
	  Compute feather scale for existing, then undo prior normalization
	*/
	computeScale((float **)vXimage,fScale, outputImage->ySize,  outputImage->xSize,fl,(float)1.0,(double)(-LARGEINT));
	undoNormalization(outputImage,vXimage,vYimage,vZimage,errorX,errorY, scaleX,scaleY,scaleZ,fScale,FALSE);
	/*
	  Loop over  images
	*/
	currentParams=params;
	tCenter=(outputImage->jd1+outputImage->jd2)*0.5;
	count=1;
	for(currentImage=images; currentImage != NULL;  currentImage=currentImage->next) {
		tOffCenter=currentImage->julDay+currentParams->nDays*0.5;
		deltaOffCenter=tOffCenter-tCenter;
		twok = 4.0* PI / currentImage->par.lambda;
		/*
		  Skip iteration if offsets requested but no offset file
		*/
		if(currentParams->offsetFlag==TRUE && currentParams->offsets.file==NULL) {
			currentParams=currentParams->next;
			fprintf(stderr,"Skipping %s :no offsets given\n",currentImage->file);
			continue;
		}
		fprintf(stderr,"\n\n--- %s\n\n", currentParams->offsets.file);
		if(currentImage->lookDir == LEFT)
			fprintf(stderr,"LEFT %i\n",count);
		else
			fprintf(stderr,"RIGHT %i Weight %f - tOff %f %f %f\n",count, currentImage->weight,deltaOffCenter,tCenter,tOffCenter);
		count++;
		/*
		  Read in image
		*/    
		if(strstr(currentImage->file,"nophase") != NULL) goto ascskip;
		getMosaicInputImage( currentImage);
		/*
		  Read offset file if needed.
		*/
		if(currentParams->offsetFlag==TRUE) {
			readOffsets( &(currentParams->offsets));
			getAzParams(&(currentParams->offsets));
		}
		/*
		  Conversions initialization
		*/
		cP=setupGeoConversions (currentImage, &azSLPixSize, &dum,&Re, &ReH,&thetaC,&ReHfixed,&thetaCfixedReH);
		fprintf(stderr,"RNear %f %f %f %f %f %fn", cP->RNear,thetaC*RTOD,thetaCfixedReH*RTOD,cP->Re,currentImage->par.H,currentImage->par.rn);
		/*
		  Get bounding box of image
		*/
		getRegion(currentImage,&iMin,&iMax,&jMin,&jMax,outputImage);
		/*
		  Now loop over output grid
		*/
		fprintf(stderr,"%i %i %i %i\n",iMin,iMax,jMin,jMax);
		for(i=iMin; i < iMax; i++) {
			if( (i % 100) == 0) fprintf(stderr,"-- %i\n",i);
			y = (outputImage->originY + i*outputImage->deltaY) * MTOKM;
			for(j=jMin; j < jMax; j++) {  
				/*
				  Convert x/y stereographic coords to lat/lon
				*/
				x = (outputImage->originX + j*outputImage->deltaX) * MTOKM;
				xytoll1(x,y,HemiSphere,&lat,&lon,Rotation,outputImage->slat);
				/*
				  Get slope and elevation
				*/
				xyGetZandSlope(lat,lon,x,y,&zSp,&zWGS84,&dzda,&dzdr,cP,  currentParams,currentImage);
				validData=FALSE;
				if(zSp > (MINELEVATION+1) && zSp < 9999.) { /* If valid z ....*/
					/* 
					   Compute phase for topography.
					*/
					llToImageNew(lat,lon,zWGS84,&range,&azimuth,currentImage);
					interpPhaseImage(currentImage,range,azimuth,&phase);
					/*
					  If shelf mask, get mask value
					*/
					sMask=GROUNDED;
					if(shelfMask != NULL) sMask=getShelfMask(shelfMask,x,y);
					if( sMask== NOSOLUTION) {phase = -LARGEINT; };
					/*
					  Process only good phase points 
					*/
					if( phase > -2.0E7 && sMask != GROUNDINGZONE  && ( !(sMask == SHELF && outputImage->noTide==TRUE))) {
						/* Compute look angles,ReH, and Range*/
						geometryInfo(cP, currentImage,azimuth, range,zSp,thetaC, &ReH,&Range,&theta,&thetaD, &psi,zSp);
						/* 
						   Compute phase due to topography and correct phase
						*/
						computePhiZ(&phiZ,zSp,azimuth, currentParams, currentImage, thetaD,Range,ReH,ReHfixed,Re,thetaCfixedReH,&phaseError);
						phase = phase - phiZ;
						/*
						  Shelf correction if mask indicates
						*/
						if(sMask == SHELF) interpTideError(&phaseError, currentImage,currentParams, x,y,psi,twok);
						/* 
						   Convert phase to delta 
						*/
						scalePhase=( 365.25/(double)(currentParams->nDays *twok*sin(psi)) ); 
						delta = phase*scalePhase;
						sigmaR=phaseError * scalePhase;
						/* 
						   Compute angle between range direction and north and xy angle
						*/
						hAngle = computeHeading(lat,lon,0, currentImage,cP);
						computeXYangle(lat,lon,&xyAngle, currentParams->xydem);
						/*	
							Get azimuth component from the offset field.
						*/
						da=interpAzOffset(range,azimuth,&(currentParams->offsets),currentImage,Range,theta ,azSLPixSize);
						sigmaA=interpAzSigma(range,azimuth, &(currentParams->offsets),currentImage,Range,theta ,azSLPixSize);
						/*
						  Note va for left sign flip done in interpAzOffset
						*/
						va = da * (365.25/(double)currentParams->nDays);
						/* 
						   zero slope correction for small vel since vertical vel is < than noise 
						*/
						if(fabs(vr) < 30.0 && fabs(va) < 30.0 ) {
							dzda = 0.0; dzdr = 0.0; 
						}
						cotanpsi=1.0/tan(psi);
						vr = (delta + va*cotanpsi*dzda)/(1.0-cotanpsi*dzdr);
						/*
						  if good date continue
						*/
						if(   da > (-LARGEINT+1)){
							/* 
							   Rotate velocity back to xy coordinates 
							*/               
							rotateFlowDirectionToXY(vr,va,&vx,&vy,xyAngle,hAngle);
							ea= sigmaA * (365.25/(double)currentParams->nDays);
							/* 
							   Weight error by sqrt(2) to avoid double avging when speckle track solution   done (if done).       
							*/
							if(outputImage->rOffsetFlag == TRUE) ea *=1.41421;
							er=sigmaR;         
							/* NOTE THIS RETURNS VARIANCES */           
							errorsToXY(er,ea,&ex,&ey,xyAngle,hAngle);
							scX=1.0/ex;
							scY=1.0/ey;
							vxTmp[i][j] = (float)vx*scX;
							vyTmp[i][j] = (float)vy*scY;
							validData=TRUE;
							if( outputImage->timeOverlapFlag==TRUE ) {
								vzTmp[i][j] =(float)(deltaOffCenter * sqrt(scX*scY));
							} else {
								if(outputImage->vzFlag==VZDEFAULT)	  	vzTmp[i][j] =(float)phase; 
								else if(outputImage->vzFlag==VZHORIZONTAL) vzTmp[i][j] =(float)delta;
								else  if(outputImage->vzFlag==VZVERTICAL)	vzTmp[i][j] = (float)delta*sin(psi)/cos(psi) ;  /* undo h by * sin, then make vert /cos */
							}
							sxTmp[i][j] = scX;
							syTmp[i][j] = scY;
							fScale[i][j] = 1.0;    /* Value for zero feathering */
							currentImage->used=TRUE;
						} 
					} 
				} else { /* End if valid z ...*/
					vxTmp[i][j]=(float)-LARGEINT; fScale[i][j]=0.0;
				} /* end else valid z */
				if(validData==FALSE) {
					vxTmp[i][j]=(float)-LARGEINT; fScale[i][j]=0.0;
				} 
			} /* j loop */
		}  /* i loop */
	ascskip:
		/*
		  Compute scale array for feathering.
		*/    
		if(fl > 0 && (iMax > 0 && jMax > 0) )
			computeScale((float **)vxTmp,fScale, outputImage->ySize,  outputImage->xSize,fl,(float)1.0,(double)(-LARGEINT));
		/*
		  Now sum current result. Falls through if no intersection (iMax&jMax==0)
		*/
		redoNormalization(currentImage->weight,outputImage,iMin, iMax, jMin,jMax, vXimage,vYimage,vZimage,errorX,errorY,
				  scaleX,scaleY,scaleZ,fScale,vxTmp,vyTmp,vzTmp, sxTmp,syTmp,  FALSE);
		/*
		  Update image pointer
		*/
		currentParams=currentParams->next;
	} /* End asc loop */
	/**************************END OF MAIN LOOP ******************************/
	endScale(outputImage,vXimage,vYimage,vZimage,errorX,errorY, scaleX,scaleY,scaleZ,FALSE);	
	fprintf(outputImage->fpLog,";\n; Returning from makeVhMosaic(.c)\n");
	fflush(outputImage->fpLog);
}
