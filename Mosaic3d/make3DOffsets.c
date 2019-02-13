#include "stdio.h"
#include"string.h"
#include <math.h>
#include <stdlib.h>
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "cRecipes/nrutil.h"
#include "source/common/common.h"
#include "mosaic3d.h"
/*
  Compute velocities from range/range offsets data
*/
void make3DOffsets(inputImageStructure *allImages,  vhParams *aParams,  xyDEM *dem,  outputImageStructure *outputImage,float fl,float timeThresh)
{   
	extern int HemiSphere;
	extern double Rotation;
	extern int sepAscDesc;
	inputImageStructure *aOffImage, *dOffImage;    /* Nominal asc/desc image */
	vhParams  *dParams; /* Velocity params */
	conversionDataStructure *aCp,*dCp; /* asc/desc coordinate conversion info */
	ShelfMask *shelfMask;	
	double lat,lon, x,y,zWGS84; /* lat/lon - x,y,z coords */
	double aRange,dRange; /* asc/desc absolute range */
	double arange, drange; /* asc/desc range coordinates in image coords */	
	double aAzimuth, dAzimuth; /* azimuth coordinates */
	double aDelta, dDelta;
	double aReH,dReH, aRe, dRe;  /* Asc/desc Earth radii and radii + alt */	
	double A[2][2], B[2][2];
	double aDemError, dDemError;;
	double aSigmaR,dSigmaR;
	double dzdx,dzdy;     /* Slopes for computing vertical velocity and 3 d solution */
	double aP,dP, aPe,dPe; /* scaled offsets and offset  error */
	double aThetaC, dThetaC,aThetaCfixed,dThetaCfixed,aThetaD, dThetaD; /* Asc/desc center look angle, and dev from center */
	double aTheta, dTheta,aPsi, dPsi;	/* Asc/desc look angle, inc angle */	
	double aSig2Base,dSig2Base;
	double tCenter,tOffCenterA,tOffCenterD,deltaOffCenter;
	double combWeight;
	double ddum1,ddum2;
	double scX,scY;
	double aZSp,dZSp;  /* asc/desc elevations corrected to local sphere */
	double vx,vy,vz; /* velocity solution */	
	float **vXimage,**vYimage,**vZimage,  **errorX,**errorY;; /* velocity and error buffers */
	float **vxTmp,**vyTmp,**vzTmp,**fScale,**sxTmp,**syTmp; /* Temp solutions */	
	float **scaleX,**scaleY,**scaleZ;  /*  scale buffers */
	float dRSLPixSize,aRSLPixSize;
	float dum;
	int iMin,iMax,jMin,jMax; /* range in pixels over which to compute solutions */
	int aa,dd, nTotal;  /* Counters for asc/desc images and total numer of images*/	
	int validData, Aset; /* Flags to indicate a valide solution, and A updates */			
	int i,j;
	unsigned char sMask;
	shelfMask = outputImage->shelfMask;
	fprintf(stderr,"fl = %f\n",(double)fl);
	fprintf(outputImage->fpLog,";\n; Entering make3DOffsets(.c)\n");
	A[0][0]=0; 	A[1][0]=0; 	A[0][1]=0; 	A[1][1]=0;
	/*
	  Pointers to output images
	*/
	setupBuffers(outputImage, &vXimage, &vYimage, &vZimage,&scaleX,&scaleY, &scaleZ,
			  &vxTmp, &vyTmp, &vzTmp,&sxTmp,&syTmp, &fScale,&errorX,&errorY);
	/*
	  Compute feather scale for existing 
	*/
	computeScale((float **)vXimage,fScale, outputImage->ySize,  outputImage->xSize,fl,(float)1.0,(double)(-LARGEINT));
	if(dem->stdLat < 50 || dem->stdLat > 80) error("mosaic3doff invalid slat for dem");
	/*
	  Init array. This undoes the prior normalization so errors are all weighted.
	*/
	undoNormalization(outputImage,vXimage,vYimage,vZimage,errorX,errorY, scaleX,scaleY,scaleZ,fScale,FALSE);
	nTotal=0;
	for(aOffImage=allImages; aOffImage->next != NULL; aOffImage=aOffImage->next) nTotal++;
	/* 
	   MAIN LOOP Loop over ascending images 
	*/
	aa=0;
	fprintf(stderr,"1\n");
	tCenter=(outputImage->jd1+outputImage->jd2)*0.5;	
	for(aOffImage=allImages; aOffImage->next != NULL; aOffImage=aOffImage->next) {
		aa++;
		/*
			Skip if no cross flag set
		*/
		tOffCenterA=aOffImage->julDay+aParams->nDays*0.5;
		/* skip if crossFlag False or weight too small (<5%) */
		if (aOffImage->crossFlag == FALSE || aOffImage->weight < 0.05)
			goto ascskip;
		/*  
		    Check if in output area
		*/			
		getRegion(aOffImage,&iMin,&iMax,&jMin,&jMax,outputImage);
		if(iMin > iMax || jMin > jMax ) goto ascskip;
		/*
		  Setup conversions 
		 */
		aRSLPixSize=aOffImage->rangePixelSize/aOffImage->nRangeLooks;
		readRangeOrRangeOffsets( &(aParams->offsets) ,ASCENDING);
		getRParams(&(aParams->offsets));
		/*
		  Setup conversion parameters
		*/
		aCp=setupGeoConversions (aOffImage, &dum, &aRSLPixSize,&aRe, &aReH,&aThetaC,&ddum1,&aThetaCfixed);
		/*
		      **** SECOND LOOP ***
	      */
		dd=0;
		dParams=aParams->next; /* Start at next element below outer loop */
		for(dOffImage=aOffImage->next; dOffImage != NULL;    dOffImage=dOffImage->next) {
			dd++;
			tOffCenterD=dOffImage->julDay+dParams->nDays*0.5;			
			if( (dOffImage->passType == aOffImage->passType && sepAscDesc == TRUE) || dOffImage->crossFlag==FALSE) goto descskip;
			/*
				Check images close enough in time
			*/
			if(fabs(aOffImage->julDay - dOffImage->julDay) > timeThresh ||  dOffImage->weight < 0.05) goto descskip;    
			fprintf(stderr,"time JD %f %f\n",aOffImage->julDay,dOffImage->julDay);
			/* 
			   Added this 7/31/2015 to skip over images with no overlap 
			*/
			getRegion(dOffImage,&iMin,&iMax,&jMin,&jMax,outputImage);
			if(iMin > iMax || jMin > jMax) goto descskip;  
			/*
			  Get region of  possible intersection - pass if not interect
			*/
			getIntersect(dOffImage,aOffImage,&iMin,&iMax,&jMin,&jMax,outputImage);
			if(iMax==0 && jMax==0)  goto descskip;
			/*
			  Init conversion stuff
			*/
			dCp=setupGeoConversions (dOffImage, &dum, &dRSLPixSize,&dRe, &dReH,&dThetaC,&ddum1,&dThetaCfixed);
			/* 
			   Compute approximate heading by sampling overlap region 
			*/
			computeSceneAlpha(outputImage,aOffImage,dOffImage,aCp,dCp,dem, &iMin,&iMax,&jMin,&jMax);
			if(iMax==0 && jMax==0)  goto descskip;			
			/*
			  Read in descending image if needed (i.e., nozero intersect).
			*/
			readRangeOrRangeOffsets( &(dParams->offsets) ,DESCENDING);
			getRParams(&(dParams->offsets));
			/*
			  Loop over output grid and compute velocities
			*/        
			fprintf(stderr,"---- Asc %i / %i Des %i \n",aa,nTotal,dd);
			Aset=FALSE;
			for(i=iMin; i < iMax; i++) {
				y = (outputImage->originY + i*outputImage->deltaY) * MTOKM;
				for(j=jMin; j < jMax; j++) {  
					/*
					  Convert x/y stereographic coords to lat/lon
					*/
					x = (outputImage->originX + j*outputImage->deltaX) * MTOKM;
					xytoll1(x,y,HemiSphere,&lat,&lon,Rotation,dem->stdLat);
					zWGS84 = getXYHeight(lat,lon,dem,0.0,ELLIPSOIDAL);
					validData=FALSE;
					/* 
					   Process points where elevation is known
					*/
					if(zWGS84 > MINELEVATION) {
						/*  Convert elevations to spherical reference	*/
						aZSp = sphericalElev(zWGS84,lat,aRe);
						dZSp = sphericalElev(zWGS84,lat,dRe);
						/*
						  Compute range azimuth positions
						*/
						llToImageNew(lat,lon,zWGS84,&arange,&aAzimuth,aOffImage);
						geometryInfo(aCp, aOffImage,aAzimuth, arange, aZSp,aThetaCfixed, &aReH,&aRange,&aTheta,&aThetaD,&aPsi,aZSp);					
						llToImageNew(lat,lon,zWGS84,&drange,&dAzimuth,dOffImage);
						geometryInfo(dCp, dOffImage,dAzimuth, drange, dZSp,dThetaCfixed, &dReH,&dRange,&dTheta,&dThetaD,&dPsi,dZSp);
						/*  Interpolate range offsets */	
						aDelta=interpRangeOffset(arange,aAzimuth,   &(aParams->offsets),  aOffImage,aRange,aThetaD ,aRSLPixSize,aTheta,&aDemError);
						dDelta=interpRangeOffset(drange,dAzimuth,   &(dParams->offsets),  dOffImage,dRange,dThetaD ,dRSLPixSize,dTheta,&dDemError);
					}  else {aDelta = -LARGEINT; dDelta=-LARGEINT;} /* End if (z >... */
					/*
					  If shelf mask, get mask value
					*/
					sMask=GROUNDED;
					if(shelfMask != NULL) sMask=getShelfMask(shelfMask,x,y);
					if(sMask== NOSOLUTION) {aDelta = -LARGEINT; dDelta=-LARGEINT;};
					/*
					  If there is valid offsets data from both images then compute velocity
					*/
					if(aDelta > (-LARGEINT+1) &&   dDelta > (-LARGEINT+1) && zWGS84 > MINELEVATION && sMask != GROUNDINGZONE
					   && ( !(sMask == SHELF && outputImage->noTide==TRUE)) ) {
						/*
						  Compute error due to baseline
						*/
						aSig2Base = computeSig2Base( sin(aThetaD),cos(aThetaD), aAzimuth,aOffImage ,&(aParams->offsets));
						dSig2Base = computeSig2Base( sin(dThetaD),cos(dThetaD), dAzimuth,dOffImage ,&(dParams->offsets));

						aSigmaR=interpRangeSigma(arange,aAzimuth, &(aParams->offsets),aOffImage,aRange,aThetaD ,aRSLPixSize);
						dSigmaR=interpRangeSigma(drange,dAzimuth, &(dParams->offsets),dOffImage,dRange,dThetaD ,dRSLPixSize);
						aSigmaR= sqrt(aSigmaR*aSigmaR + aDemError*aDemError +  aSig2Base);
						dSigmaR= sqrt(dSigmaR*dSigmaR + dDemError*dDemError +  dSig2Base);
						/*
						  Tide corrections
						*/  
						if(sMask == SHELF) {
							/* Interp tide errors, set twok (last param) as 1.0 for offsets */
							interpTideError(&aSigmaR, aOffImage,	aParams, x,y, aPsi,1.0);
							interpTideError(&dSigmaR, dOffImage,	dParams, x,y, dPsi,1.0);						
							aDelta -= - aOffImage->tideCorrection * cos(aPsi) * (double)aParams->nDays/365.25;
							dDelta -= - dOffImage->tideCorrection * cos(dPsi)  *	(double)dParams->nDays/365.25;
						} /* ENd if(smask... */
						/*  Update A every set of 10 pixels.  Use Aset to force computation on first calc. */ 						
						if( ((i % 3) == 0 || (j % 3) == 0) || Aset==FALSE ) {
							computeA(lat,lon,x,y,aOffImage,dOffImage,A);
							Aset=TRUE;
						}
						/*
						  Only pursue solution if sufficient difference  in angles for 3d solution
						*/
						if(A[0][0] !=-LARGEINT) {
							/*  Compute B (note B is really C in the TGARS paper	*/
							computeB(x,y,zWGS84,B,&dzdx,&dzdy,aPsi,dPsi,(xyDEM *)dem);
							/*  Scale offsets for velocity computation (scale for m/yr)	*/
							aP = 365.25 * aDelta / ( (double) (aParams->nDays) * sin(aPsi) );
							dP = 365.25 * dDelta / ( (double) (dParams->nDays) * sin(dPsi) );
							aPe = 365.25 * aSigmaR/(aParams->nDays * sin(aPsi));
							dPe = 365.25 * dSigmaR /(dParams->nDays * sin(dPsi));
							/*  Compute velocity */
							computeVxy(aP,dP,aPe,dPe,A,B,&vx,&vy,&scX,&scY);
							/*  Compute vertical velocity	*/
							vz = vx * dzdx + vy * dzdy;
							/*  Update output arrays */
							vxTmp[i][j] = vx*scX;
							vyTmp[i][j] = vy*scY;
							validData=TRUE;
							if( outputImage->timeOverlapFlag==TRUE ) {
								/* For lack of better option, use the average of the two data takes */
								deltaOffCenter = 0.5*(tOffCenterA+tOffCenterD-2.0*tCenter);
								vzTmp[i][j] = (float)(deltaOffCenter * sqrt(scX*scY));
							} else {
								vzTmp[i][j] = vz;
							}
							sxTmp[i][j] = scX; /* This is summing up 1/sigma^2*/
							syTmp[i][j] = scY;
							fScale[i][j] = 1.0;        /* Value for zero feathering */
							vzTmp[i][j]=dzdx;
						}
					}
					if(validData==FALSE) {
						vxTmp[i][j]=(float)-LARGEINT; fScale[i][j]=0.0;
					}
				} /* j loop */
				if( (i % 100) == 0) fprintf(stderr,"-- %i %f %f %f %f  thetas %f %f\n",i,A[0][0],A[0][1],A[1][0],A[1][1], aTheta*RTOD,dTheta*RTOD);								
			}  /* i loop */
			/*
			  Compute scale array for feathering.
			*/    
			if(fl > 0 && (iMax > 0 && jMax > 0) )
				computeScaleLS((float **)vxTmp,fScale, outputImage->ySize, outputImage->xSize,fl,(float)1.0,(double)(-LARGEINT),iMin,iMax,jMin,jMax);
			/*
			  Now sum current result. Falls through if no intersection (iMax&jMax==0)
			*/
			if( outputImage->timeOverlapFlag==TRUE )  {
				combWeight=sqrt(aOffImage->weight * dOffImage->weight );
				fprintf(stderr,"\033[1mComb weight = %lf |Ta-Td| %lf\033[0m\n",combWeight, fabs(aOffImage->julDay - dOffImage->julDay)) ;
			} else combWeight=1.0;

			redoNormalization(combWeight,outputImage,iMin, iMax, jMin,jMax, vXimage,vYimage,vZimage,errorX,errorY,
				  scaleX,scaleY,scaleZ,fScale,vxTmp,vyTmp,vzTmp, sxTmp,syTmp,  FALSE);
		descskip:
			/* Update  nominal descending image pointer */							
			dParams=dParams->next;
		}   /* End desc loop */
	ascskip:
		/* Update  nominal ascending image pointer */					
		aParams=aParams->next;
	} /* End asc loop */
	/**************************END OF MAIN LOOP ******************************/
	fprintf(stderr,"Out of main loop\n");
	/*
	  Adjust scale
	*/
	 endScale(outputImage,vXimage,vYimage,vZimage,errorX,errorY, scaleX,scaleY,scaleZ,FALSE);	
	fprintf(outputImage->fpLog,";\n; Returning from make3DOffsets(.c)\n");
	fflush(outputImage->fpLog);
}













