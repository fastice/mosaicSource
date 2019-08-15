#include "stdio.h"
#include"string.h"
#include <math.h>
#include <stdlib.h>
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "cRecipes/nrutil.h"
#include "mosaicSource/common/common.h"
#include "mosaic3d.h"

static void setBuffer(inputImageStructure *inputImage,float *buf);
static double computePhiZM3d(double *thetaD,  double z,double azimuth,  vhParams *vhParam,inputImageStructure *phaseImage,
			     double Range,double Re,   double ReH,double ReHfixed, double thetaC,double thetaCfixedReH,double *phaseError);
/*
************************ Estimate 3D velocity from phase **************************
*/

/* convert range and time on the lat/lon on the ellipsoid */
static void RTtoLatLon(inputImageStructure *inputImage,double r,double myTime,double *lat,double *lon) {
	int n;
	stateV *sv;
	double xs,ys,zs,vsx,vsy,vsz;
	sv=&(inputImage->sv);
	n=(long int)((myTime - sv->times[1])/(sv->deltaT)+.5);
	n=min(max(0,n-2), sv->nState-NUSESTATE);	
	polintVec(&(sv->times[n]), &(sv->x[n]), &(sv->y[n]),&(sv->z[n]),  &(sv->vx[n]), &(sv->vy[n]),&(sv->vz[n]), myTime, &xs,&ys ,&zs, &vsx,&vsy ,&vsz);
	smlocateZD( xs*MTOKM,  ys*MTOKM,  zs*MTOKM, vsx*MTOKM,vsy*MTOKM, vsz*MTOKM, r*MTOKM,lat,lon,(double)(inputImage->lookDir),0.0);
}

static void printLatLon(	inputImageStructure  *inputImage) {
	int i1,j1;
	double dr,dt,lat,lon;
	dr=(inputImage->par.rf-inputImage->par.rn)/(1);
	dt=(inputImage->azimuthSize* inputImage->nAzimuthLooks/inputImage->par.prf)/(1);
	fprintf(stderr,"%f %f\n",dr,dt);
	for(i1=0; i1 < 2; i1++) 		for(j1=0; j1 < 2; j1++) {
			fprintf(stderr,"%f %f \n",inputImage->par.rc,inputImage->cpAll.sTime);
			RTtoLatLon(inputImage,inputImage->par.rc+i1*dr,inputImage->cpAll.sTime+j1*dt,&lat,&lon); /* lat/lon in first image */
			fprintf(stderr,"%f %f \n",lat,lon);
		}
	RTtoLatLon(inputImage,inputImage->par.rc+0.5*dr,inputImage->cpAll.sTime+0.5*dt,&lat,&lon); /* lat/lon in first image */
	fprintf(stderr,"%f %f \n",lat,lon);	
}

void make3DMosaic(inputImageStructure *ascImages, inputImageStructure *descImages, 
		  vhParams *ascParams, vhParams *descParams,  xyDEM *dem,  outputImageStructure *outputImage,float fl,  int no3d)
{   
	extern int HemiSphere;
	extern double Rotation;
	extern float *AImageBuffer, *DImageBuffer;
	extern int sepAscDesc;
	conversionDataStructure *aCp,*dCp; /* asc/desc coordinate conversion info */
	inputImageStructure  *allImages, *aPhaseImage, *dPhaseImage; /* list of all images, and individual asc/desc images */
	vhParams *aParams, *dParams;
	ShelfMask *shelfMask;	
	double lat,lon, x,y,zWGS84; /* lat/lon - x,y,z coords */
	double aRange,dRange; /* asc/desc absolute range */
	double arange, drange; /* asc/desc range coordinates in image coords */
	double aAzimuth, dAzimuth;
	double aPhase, dPhase;
	double aReH,dReH, aRe, dRe;  /* Asc/desc Earth radii and radii + alt */
	double aReHfixed,aThetaCfixedReH,  dReHfixed,dThetaCfixedReH; /* Fixed geo params */
	double scX,scY; /* X,Y scale factors */	
	double A[2][2], B[2][2]; /* 3d solution matrices */
	double phaseErrorA,phaseErrorD; /* asc/desc phase errors */
	double dzdx,dzdy;     /* Slopes for computing vertical velocity and 3 d solution */
	double aPhiZ, dPhiZ;  /* Phase due to topopgraphy */
	double aP,dP, aPe,dPe; /* scaled phases and phase error */
	double aThetaC, dThetaC,aThetaD, dThetaD; /* Asc/desc center look angle, and dev from center */
	double aTheta, dTheta,aPsi, dPsi;	/* Asc/desc look angle, inc angle */
	double tCenter,tOffCenterA,tOffCenterD,deltaOffCenter; /* Variables for tracking time offsets */
	double combWeight;  /* Weight based on time overlap */
	double aZSp,dZSp;  /* asc/desc elevations corrected to local sphere */
	double twokA,twokD; /* 4pi/lambda */	
	double vx,vy,vz; /* velocity solution */
	float **vXimage,**vYimage,**vZimage,  **errorX,**errorY;; /* velocity and error buffers */
	float **vxTmp,**vyTmp,**vzTmp,**fScale,**sxTmp,**syTmp; /* Temp solutions */	
	float **scaleX,**scaleY,**scaleZ;  /*  scale buffers */
	float dum1,dum2;  /* Placeholder dummys for function calls */
	int validData, Aset; /* Flags to indicate a valide solution, and A updates */		
	int iMin,iMax,jMin,jMax; /* range in pixels over which to compute solutions */
	int aa,dd;  /* Counters for asc/desc images */
	int i,j,i1,j1;
	unsigned char sMask;

	fprintf(outputImage->fpLog,";\n; Entering make3DMosaic(.c)\n");
	if(no3d==TRUE) {
		fprintf(outputImage->fpLog,";\n; no3d flag set, returning;\n; Returning from make3DMosaic(.c)\n");
		return;
	}
	A[0][0]=0; 	A[1][0]=0; 	A[0][1]=0; 	A[1][1]=0;
	shelfMask = outputImage->shelfMask;	
	/*
	  Pointers to output images
	*/
	setupBuffers(outputImage, &vXimage, &vYimage, &vZimage,&scaleX,&scaleY, &scaleZ,
		     &vxTmp, &vyTmp, &vzTmp,&sxTmp,&syTmp, &fScale,&errorX,&errorY);
	if(dem->stdLat < 50 || dem->stdLat > 80) error("mosaic3d invalid slat for dem");
	/*
	  Compute feather scale for existing, and undo normalization
	*/
	computeScale((float **)vXimage,fScale, outputImage->ySize,  outputImage->xSize,fl,(float)1.0,(double)(-LARGEINT));
	undoNormalization(outputImage,vXimage,vYimage,vZimage,errorX,errorY, scaleX,scaleY,scaleZ,fScale,FALSE);
	allImages=ascImages; /* Added 5/30 to avoid using asc/desc */    
	aParams=ascParams;
	/* 
	   MAIN LOOP Loop over ascending images 
	*/
	aa=0;
	tCenter=(outputImage->jd1+outputImage->jd2)*0.5;
	for(aPhaseImage=allImages; aPhaseImage->next != NULL; aPhaseImage=aPhaseImage->next) {
		aa++;
		/* Use for calculating time skew 09/21/17 */
		tOffCenterA=aPhaseImage->julDay+aParams->nDays*0.5;
		/*  
		    Break inner loop if outer image is a nophase
		*/
		fprintf(stderr,"A %s\n",aPhaseImage->file);
		if(strstr(aPhaseImage->file,"nophase") != NULL) goto ascskip;
		/* Added this 7/31/2015 to skip over images with no overlap */
		getRegion(aPhaseImage,&iMin,&iMax,&jMin,&jMax,outputImage);
		/* Skip if no data in range */
		if(iMin > iMax || jMin > jMax) goto ascskip;
		/*  Set buffer, memory channel for sharedmem, and read image		*/		
		setBuffer(aPhaseImage,AImageBuffer);
		getMosaicInputImage(aPhaseImage);
		aPhaseImage->memChan=MEM1;
		twokA=(4.0*PI)/aPhaseImage->par.lambda;
		/*  Setup conversion parameters		*/
		aCp=setupGeoConversions (aPhaseImage, &dum1, &dum2,&aRe, &aReH,&aThetaC,&aReHfixed,&aThetaCfixedReH);
		/* 
		   SECOND LOOP: Loop over descending images:  Changed 05/30/07 to search all images below one from aPhase loop 
		*/   
		dd=0;
		dParams=aParams->next; /* Start at next element below outer loop */
		for(dPhaseImage=aPhaseImage->next; dPhaseImage != NULL;    dPhaseImage=dPhaseImage->next) {
			dd++;
			/* Use for calculating time skew 09/21/17 */			
			tOffCenterD=dPhaseImage->julDay+dParams->nDays*0.5;		
			/*   Check inner image is not a nophase and only proceed with this iteration of the inner loop if so. */
			if(strstr(dPhaseImage->file,"nophase") != NULL ) goto descskip;
			/* 
			   Moved up 9/14/2016 to avoid init ll - had to modify getRegion to use lat/lon from geodat
			*/
			getRegion(dPhaseImage,&iMin,&iMax,&jMin,&jMax,outputImage);

			if(iMin > iMax || jMin > jMax) {
				/* This file has no overlap so, set to nophase - for future loops */
				strncpy(dPhaseImage->file,"nophase",7); dPhaseImage->file[7]='\0';
				fprintf(stderr,"skip %s %i %i\n",dPhaseImage->file,dd,aa);
				goto descskip;
			}
			if(dPhaseImage->passType == aPhaseImage->passType && sepAscDesc == TRUE) goto descskip;    
			/*
			  Get region of  possible intersection;   This provides final iMin/iMax.. to loop over 
			*/
			fprintf(stderr,"%i %i %i %i\n",iMin,iMax,jMin,jMax);
			getIntersect(dPhaseImage,aPhaseImage,&iMin,&iMax,&jMin,&jMax,outputImage);
			if(iMax==0 && jMax==0)  goto descskip;
			/*
			  Init conversion stuff
			*/
			dPhaseImage->memChan=MEM2;
			dCp=setupGeoConversions (dPhaseImage, &dum1, &dum2,&dRe, &dReH,&dThetaC,&dReHfixed,&dThetaCfixedReH);			
			fprintf(stderr,"init dPhaseImage->file %s\n",dPhaseImage->file);
			fprintf(stderr,"dPhaseImage %s %3.0f %5.1f %4i\n",  dPhaseImage->file,dParams->nDays,dPhaseImage->par.lambda,dd );
			fprintf(stderr,"aPhaseImage %s %3.0f %5.1f %4i\n\n",aPhaseImage->file,aParams->nDays,aPhaseImage->par.lambda,aa );
			twokD=(4.0*PI)/dPhaseImage->par.lambda;						
			/*   Compute approximate heading by sampling overlap region  - set iMax,jMax zero if no good solution */
			computeSceneAlpha(outputImage,aPhaseImage,dPhaseImage,aCp,dCp,dem, &iMin,&iMax,&jMin,&jMax);
			if(iMax==0 && jMax==0)  goto descskip; /* no data in range, so skip */			
			/*  Read in descending image if needed (i.e., nozero intersect).	*/
			setBuffer(dPhaseImage,DImageBuffer);
			getMosaicInputImage(dPhaseImage);
			/*
			  Loop over output grid and compute velocities
			*/        
			Aset=FALSE;
			for(i=iMin; i < iMax; i++) {
				if( (i % 100) == 0) fprintf(stderr,"--+ %i %f %f %f %f \n",i,A[0][0],A[0][1],A[1][0],A[1][1]);
				/* y - coordinate */
				y = (outputImage->originY + i*outputImage->deltaY) * MTOKM;
				for(j=jMin; j < jMax; j++) {  
					/*  x-coordinate, then convert x/y stereographic coords to lat/lon	*/
					x = (outputImage->originX + j*outputImage->deltaX) * MTOKM;
					xytoll1(x,y,HemiSphere,&lat,&lon,Rotation,dem->stdLat);
					zWGS84 = getXYHeight(lat,lon,dem,0.0,ELLIPSOIDAL);
					/* Nominal phase error */
					phaseErrorA=PI;phaseErrorD=PI;
					/* 
					   Process points where elevation is known
					*/
					validData=FALSE; /* Assume didn't work until successful */
					if(zWGS84 > MINELEVATION) {
						/*   Convert elevations to spherical reference	*/
						aZSp = sphericalElev(zWGS84,lat,aRe);
						dZSp = sphericalElev(zWGS84,lat,dRe);
						/*
						  Compute range azimuth position
						*/
						llToImageNew(lat,lon,zWGS84,&arange,&aAzimuth,aPhaseImage);
						geometryInfo(aCp, aPhaseImage,aAzimuth, arange, aZSp,aThetaC, &aReH,&aRange,&aTheta,&aThetaD, &aPsi,aZSp);											
						llToImageNew(lat,lon,zWGS84,&drange,&dAzimuth,dPhaseImage);
						geometryInfo(dCp, dPhaseImage,dAzimuth, drange, dZSp,dThetaC, &dReH,&dRange,&dTheta,&dThetaD,&dPsi,dZSp);																	
						/*  Interpolate Phase*/
						interpPhaseImage(aPhaseImage,arange,aAzimuth,&aPhase);
						interpPhaseImage(dPhaseImage,drange,dAzimuth,&dPhase);
						/*  If shelf mask, get mask value */
						sMask=GROUNDED;
						if(shelfMask != NULL) sMask=getShelfMask(shelfMask,x,y);
						if(sMask== NOSOLUTION) {aPhase = -LARGEINT; dPhase=-LARGEINT;};
						/*
						  If there is valid phase data from both images then compute velocity
						*/
						if(aPhase > -LARGEINT &&  dPhase > -LARGEINT && zWGS84 > MINELEVATION && sMask != GROUNDINGZONE
						   && ( !(sMask == SHELF && outputImage->noTide==TRUE)) ) {
							/*
							  Compute phase due to topography.  Everything is looped through
							  and only the pairs where there is a good angular seperation are used. In general, one will be asc and one will be desc, but they could be flipped. 
							  As a consequence, this next step has to look at the flag to see if it should flip the azimuth coordinate when computing the baseline. 
							*/
							aPhiZ = computePhiZM3d(&aThetaD,aZSp,aAzimuth,aParams,aPhaseImage,aRange,aRe,aReH,aReHfixed,aThetaC, aThetaCfixedReH,&phaseErrorA);
							aPhase = aPhase - aPhiZ;
							dPhiZ = computePhiZM3d(&dThetaD,dZSp,dAzimuth,dParams,dPhaseImage,dRange,dRe,dReH,dReHfixed,dThetaC,dThetaCfixedReH,&phaseErrorD); 
							dPhase = dPhase - dPhiZ;
							/*  Tide corrections	*/  
							if(sMask == SHELF) {
								/* update tide correct, and compute phaseImage->tideCorrection */
								interpTideError(&phaseErrorA, aPhaseImage,	aParams, x,y, aPsi,twokA);
								interpTideError(&phaseErrorD, dPhaseImage,	dParams, x,y, dPsi,twokD);
								aPhase -= - aPhaseImage->tideCorrection * cos(aPsi) * twokA* (double)aParams->nDays/365.25;
								dPhase -= - dPhaseImage->tideCorrection * cos(dPsi) * twokD *(double)dParams->nDays/365.25;
							} /* ENd if(smask... */
							/*  Update A every set of 10 pixels.  Use Aset to force computation on first calc. */ 
							if( ((i % 3) == 0 || (j % 3) == 0) || Aset==FALSE ) {
								computeA(lat,lon,x,y,aPhaseImage,dPhaseImage,A);
								Aset=TRUE;
							}
							/*
							  Only pursue solution if sufficient difference  in angles for 3d solution
							*/
							if(A[0][0] !=-LARGEINT) {
								/*  Compute B (note B is really C in the TGARS paper	*/
								computeB(x,y,zWGS84,B,&dzdx,&dzdy,aPsi,dPsi,(xyDEM *)dem);
								/*  Scale phases for velocity computation (scale for m/yr)	*/
								aP = 365.25 * aPhase/(twokA*aParams->nDays * sin(aPsi));
								dP = 365.25 * dPhase/(twokD*dParams->nDays * sin(dPsi));
								aPe = 365.25 * phaseErrorA/(twokA*aParams->nDays * sin(aPsi));
								dPe = 365.25 * phaseErrorD/(twokD*dParams->nDays * sin(dPsi));
								/*  Compute velocity */
								computeVxy(aP,dP,aPe,dPe,A,B,&vx,&vy,&scX,&scY);
								/*  Compute vertical velocity	*/
								vz = vx * dzdx + vy * dzdy;
								/*
								  Update output arrays
								*/
								if(! (scX > -1000. && scX < 1000.) ) error("invalid velocity %f %f %f %f %f %f\n",vx,vy,phaseErrorA,phaseErrorD,aPe,dPe);
								vxTmp[i][j] = vx*scX;
								vyTmp[i][j] = vy*scY;
								if( outputImage->timeOverlapFlag==TRUE ) {
									/* For lack of better option, use the average of the two data takes */
									deltaOffCenter = 0.5*(tOffCenterA+tOffCenterD-2.0*tCenter);
									vzTmp[i][j] = (float)(deltaOffCenter * sqrt(scX*scY));
								} else {  vzTmp[i][j] = dP;	}							
								sxTmp[i][j] = scX; /* This is summing up 1/sigma^2*/
								syTmp[i][j] = scY;
								fScale[i][j] = 1.0;        /* Value for zero feathering */
								aPhaseImage->used=TRUE;
								dPhaseImage->used=TRUE;
								validData=TRUE;
							}
						} 
					}
					if(validData==FALSE) {
						vxTmp[i][j]=(float)-LARGEINT; fScale[i][j]=0.0;
					}
				} /* j loop */
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
				combWeight=sqrt(aPhaseImage->weight * dPhaseImage->weight );
				fprintf(stderr,"\033[1mComb weight = %lf |Ta-Td| %lf\033[0m\n",combWeight, fabs(aPhaseImage->julDay - dPhaseImage->julDay)) ;
			} else combWeight=1.0;
			redoNormalization(combWeight,outputImage,iMin, iMax, jMin,jMax, vXimage,vYimage,vZimage,errorX,errorY,
					  scaleX,scaleY,scaleZ,fScale,vxTmp,vyTmp,vzTmp, sxTmp,syTmp,  FALSE);
			/* ******************************
			   Use end of goto used to skip inner loop for nophase 
			************************************/
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
	fprintf(outputImage->fpLog,";\n; Returning from make3DMosaic(.c)\n");
	fflush(outputImage->fpLog);
}


/*
  Compute phase due to topography
*/
static double computePhiZM3d(double *thetaD, double z,double azimuth,  vhParams *vhParam,inputImageStructure *phaseImage,    double Range,double Re,
			     double ReH,double ReHfixed, double thetaC,double thetaCfixedReH,double *phaseError)                                     
{
	double normAzimuth,imageLength;
	double bn, bp,bSq, delta;
	double theta, thetaDFlat;
	double phiZ;
	double xsq;
	double sig2Base;
	double twok;
	int  i,j;
	double v[7],tmpV[7];
	double sinThetaD, cosThetaD;
	/*
	  Compute baseline
	*/
	twok=4.0*PI/phaseImage->par.lambda;
	imageLength = (double)phaseImage->azimuthSize;
	normAzimuth = (azimuth-0.5*imageLength)/imageLength;
	xsq=normAzimuth*normAzimuth;
	bn = vhParam->Bn + normAzimuth * vhParam->dBn + xsq * vhParam->dBnQ ;
	bp = vhParam->Bp + normAzimuth * vhParam->dBp + xsq * vhParam->dBpQ;
	bSq = bn*bn + bp*bp;
	/* 
	   Compute look angles for nonflat surface
	*/
	theta = acos( (Range*Range + ReH*ReH - pow((Re+z),2.0))/(2.0 * ReH * Range) );
	*thetaD = theta - thetaC;
	sinThetaD=sin(*thetaD);
	cosThetaD=cos(*thetaD);
	/*
	  Vector used to comptue baseline error 
	*/
	v[1]= -twok * sinThetaD;  				v[2]= -twok* cosThetaD ;
	v[3]= -twok * sinThetaD* normAzimuth; 	v[4]= -twok * cosThetaD * normAzimuth;
	v[5]= -twok * sinThetaD* xsq;			v[6]= -twok * cosThetaD * xsq;
	/*
	  C*v
	*/
	for(i=1; i <= 6; i++) {
		tmpV[i]=0;
		for(j=1; j<=6; j++)  tmpV[i] += vhParam->C[i][j] * v[j];
	}
	/*
	  sigma^2= vT * C*v
	*/
	sig2Base=0.0;
	for(j=1; j<=6; j++)  sig2Base += tmpV[j] * v[j];
	/*
	  Assume error more than 1/2 of a fringe, reflects tie point error rather than phase noise
	*/
	*phaseError= sqrt( sig2Base  + min(PI,vhParam->sigma) * min(PI,vhParam->sigma) ); /* note this is returning phase error as sigma */
	/* 
	   This delta uses a varying ReH because it is the computationally correct version
	*/
	delta = sqrt( pow(Range,2.0) -2.0*Range* ( bn*sinThetaD + bp*cosThetaD) + bSq )  - Range;    
	/*
	  Compute thetaD for a flat surface
	  This delta uses a fixed ReH because ultimately it is adding back what was subtracted
	  interferogram phase = phi - phi_flat=(phit + phiv) - phiflat
	  phiz=phithat - phiflat
	  phase - phiz= (phit+phiv)-phiflat - (phithat - phiflat) = phiv + (phit-phithat)
	*/ 
	theta = acos( (Range*Range + ReHfixed*ReHfixed - Re*Re)/( 2.0 * ReHfixed * Range ) );
	thetaDFlat = theta - thetaCfixedReH;
	/* 
	   Substract flat earth phase 
	*/
	delta -=  -bn * sin(thetaDFlat) - bp * cos(thetaDFlat) + bSq * 0.5 / Range;
	phiZ = delta * twok;
	return phiZ;
}





static void setBuffer(inputImageStructure *inputImage,float *buf)
{
	int i;
	for(i=0; i < inputImage->azimuthSize; i++) {
		inputImage->image[i]= &(buf[i*inputImage->rangeSize]);
	}	
}
