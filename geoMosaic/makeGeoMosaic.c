#include "stdio.h"
#include"string.h"
#include <math.h>
#include "source/common/common.h"
#include "geomosaic.h"
#include "cRecipes/nrutil.h"
#include <stdlib.h>

/*
************************ Geocode image InSAR DEM. **************************
*/
int sepAscDesc;
static void readDEMInput(inputImageStructure *inputImage, char *insarDEMFile);
static float interpolatePowerInputImage(inputImageStructure inputImage, double range,double azimuth);
static void freeConversion(inputImageStructure *image);
static float polypat(float theta);
static float polyALOS(float theta);
static void smoothImage( inputImageStructure *inputImage, int smoothL);
static void rsatFinalCal(float **image, int xSize, int ySize);
static void logSigma(float **image,  int xSize, int ySize) ;

/* Avoids using mosaicvh include  - this routine is in mosaicVh_p/makeVhMosaic.c */
void getRegion(inputImageStructure *image, int *iMin,int *iMax, int *jMin,int *jMax,  outputImageStructure *outputImage);
static void mallocTmpBuffers(float ***imageTmp, float ***scaleTmp,  outputImageStructure *outputImage);
static void getGeoMosaicImage(inputImageStructure *inputImage, int *imageDate, int smoothL);
static void applyCorrections(float *value,inputImageStructure *inputImage, double range, double azimuth, double h);
static void geoMosaicScaling(inputImageStructure *inputImage, float **image, float **imageTmp, float **scale, float **scaleTmp, void *dem, outputImageStructure *outputImage, int orbitPriority,int imageDate, int iMin,int iMax, int jMin,int jMax);
static void finalReScale(outputImageStructure *outputImage, float **image, float **scale, int orbitPriority);
static void initImageBuffers(outputImageStructure *outputImage, float **scaleTmp, int orbitPriority);
/*
  Make mosaic of insar and other dems.
*/
void makeGeoMosaic(inputImageStructure *inputImage,   outputImageStructure outputImage,
		   void *dem, int nFiles, int maxR, int maxA,   char **imageFiles,float fl, int smoothL, int orbitPriority)
{   
	extern int HemiSphere;
	extern double Rotation;
	extern int nearestDate;
	extern int hybridZ; 
	extern int rsatFineCal; 
	extern int noPower;
	extern int S1Cal;
	float **scale;

	FILE *fp;
	double x,y;
	double lat,lon;
	double range,azimuth, rangeTmp;
	inputImageStructure tmp;
	double h,hWGS;
	int imageDate;
	double tmp1,tmp2;
	float value; 
	float **image;
	float **tmpImage;
	float **imageTmp, **scaleTmp;
	int i,is,il,s1,s2,l1,l2,i1,j1,k,k1,k2;
	int iMin,iMax,jMin,jMax;
	/*
	  Malloc space for tmp buffers
	*/
	mallocTmpBuffers(&imageTmp, &scaleTmp, &outputImage);
	/* Init image,scale */	
	initImageBuffers(&outputImage,scaleTmp, orbitPriority);
	image = (float **)outputImage.image;
	scale = (float **)outputImage.scale;
	fprintf(stderr,"\norbitPriority %i\n",orbitPriority);
	/*
	  Main mosaicing loop
	*/
	fprintf(stderr,"STARTING MOSAICING\n");
	for(i=0; i < nFiles; i++) {
		/*
		  Read in input image
		*/
		inputImage[i].file=imageFiles[i];
		/*  Get bounding box of image		*/
		getRegion(&(inputImage[i]),&iMin,&iMax,&jMin,&jMax,&outputImage);
		fprintf(stderr,"*** iMin,iMax, jMin,jMax ++ w %i %i %i %i ++ %f\n",iMin,iMax,jMin,jMax,inputImage[i].weight );
		if(inputImage[i].weight < 0.0001) { /* Zero weight so force skip */
			iMin=1; iMax=-1;  jMin=1; jMax=-1; 
		} else {
			getGeoMosaicImage(&(inputImage[i]), &imageDate,  smoothL);
		}
		tmp=(inputImage[i]);
		/*
		  Loop over output grid
		*/
		fprintf(stderr,"%s\n",imageFiles[i]);
		for (i1=iMin; i1 < iMax; i1++) {
			if( (i1 % 100) == 0) fprintf(stderr,"-- %i\n",i1);  
			y = (outputImage.originY + i1*outputImage.deltaY) * MTOKM;
			for(j1=jMin; j1 < jMax; j1++) {            
				/*
				  Convert x/y stereographic coords to lat/lon
				*/
				x = (outputImage.originX + j1*outputImage.deltaX) * MTOKM;
				xytoll1(x,y,HemiSphere,&lat,&lon,Rotation,outputImage.slat);
				/*
				  Get height for given lat/lon  
				*/
				h=getXYHeight(lat,lon,dem,inputImage[i].cpAll.Re,SPHERICAL);
				hWGS=sphericalToWGSElev(h, lat,inputImage[i].cpAll.Re);
				/*
				  Convert lat/lon to image coordinates
				*/       
				llToImageNew(lat,lon,hWGS,&range,&azimuth,&(inputImage[i]));
				/*
				  Interpolate image
				*/
				value=interpolatePowerInputImage(inputImage[i],range,azimuth);
				/* calibration */
				applyCorrections(&value,&(inputImage[i]),range, azimuth,h);
				/* Scaling */
				imageTmp[i1][j1] =value;   /*DEBUG value;*/
				if(value > 0 || (noPower >0 && value > (-LARGEINT+10))) {
					scaleTmp[i1][j1] = 1;                
				} 
			} /* End j1 */
		} /* End i1 */
		/*
		  Compute scale array for feathering.
		*/    
		if(fl > 0 && (iMax > 0 && jMax > 0)  && (iMax > iMin && jMax > jMin) )
			computeScale(imageTmp,scaleTmp, outputImage.ySize, outputImage.xSize,fl, inputImage[i].weight,(float)0.0); 
		/*
		  Now sum current result. Falls through if no intersection (iMax&jMax==0)
		*/         
		geoMosaicScaling(&(inputImage[i]), image,imageTmp,scale, scaleTmp, dem, &outputImage, orbitPriority,imageDate,
				 iMin, iMax, jMin, jMax);
	} /* End for i=0; i < nFiles */
	/*
	  Final rescaling if not nearestDate
	*/
	finalReScale(&outputImage, image, scale, orbitPriority);
	/* Convert values to log db if calibrated  added 11/18/2013 */
	if(rsatFineCal == TRUE) {
		rsatFinalCal( image, outputImage.xSize, outputImage.ySize);		
	}
	if(rsatFineCal == TRUE || S1Cal == TRUE) {
		logSigma( image, outputImage.xSize, outputImage.ySize);		
	}
	fprintf(stderr,"RETURN \n");
	return;
}

static void finalReScale(outputImageStructure *outputImage, float **image, float **scale, int orbitPriority) {
	int i1,j1;
	for (i1=0; i1 < outputImage->ySize; i1++) {
		for(j1=0; j1 < outputImage->xSize; j1++) {
			if(scale[i1][j1] > 0) {
				/* this assumes that dates are larger than 1000, and we won't sum more than 1000*/
				if(scale[i1][j1] < 10000 && orbitPriority < 0) image[i1][j1] /=scale[i1][j1]; 
			}
		}
	}
}

static void geoMosaicScaling(inputImageStructure *inputImage, float **image, float **imageTmp, float **scale, float **scaleTmp, void *dem, outputImageStructure *outputImage, int orbitPriority,int imageDate, int iMin,int iMax, int jMin,int jMax){
	extern int HemiSphere;
	extern double Rotation;	
	extern int hybridZ;
	extern int nearestDate;
	extern int noPower;
	double x,y, hWGS;
	double lat,lon;
	int i1,j1; 
	for(i1=iMin; i1 < iMax; i1++) {
		y = (outputImage->originY + i1*outputImage->deltaY) * MTOKM;
		for(j1=jMin; j1 < jMax; j1++) {
			x = (outputImage->originX + j1*outputImage->deltaX) * MTOKM;
			/* Get elevation */
			if(hybridZ > 0) {
				xytoll1(x,y,HemiSphere,&lat,&lon,Rotation,outputImage->slat);
				hWGS=getXYHeight(lat,lon,dem,inputImage->cpAll.Re,ELLIPSOIDAL);
			} else hWGS=hybridZ-1; /* This will force skip */

			if(imageTmp[i1][j1] > 0 || (noPower > 0 && scaleTmp[i1][j1] > 0) ) { /* Points with valid datat */
				/* case for no nearestDate, no orbitPriority, or hWGS override */
				if((nearestDate <  0  || (nearestDate > 0 && (int)hWGS > hybridZ)) && orbitPriority < 0 ) { /* Summing data or non-nearest date or hybridZ*/
					image[i1][j1] += imageTmp[i1][j1] * scaleTmp[i1][j1]*inputImage->weight;
					scale[i1][j1] += scaleTmp[i1][j1];    
				} else if(orbitPriority > -1 ) { /* orbit priority case */
					if(orbitPriority == ASCENDING) {
						/* either put in data if none already, or if not ascending replace */					
						if( scale[i1][j1] > 0.1 || scale[i1][j1] ==DESCENDING) {							
							image[i1][j1]= imageTmp[i1][j1]*inputImage->weight;
							scale[i1][j1] = inputImage->passType;
						}
					} else if(orbitPriority == DESCENDING) {
						/* either put in data if none already, or if not descending replace */
						if( scale[i1][j1] < -0.1 || scale[i1][j1] == ASCENDING) {
							image[i1][j1]= imageTmp[i1][j1]*inputImage->weight;
							scale[i1][j1] = inputImage->passType;
						}							
					} else error("bad orbit priority flag");
				} else { /* Nearest date, if orbitPriority not st */
					if( (abs(nearestDate -imageDate) < abs(nearestDate-scale[i1][j1])) && scaleTmp[i1][j1] > 0 ) {
						image[i1][j1]= imageTmp[i1][j1]*inputImage->weight;
						scale[i1][j1] = imageDate;
					}
				} /* end else */
			} /* end if imageTmp...*/
			scaleTmp[i1][j1] = 0.0;
		} /* end j1 */
	} /* end i1=iMin sum current ... */
}

static void applyCorrections(float *value, inputImageStructure *inputImage, double range, double azimuth, double h) {
	extern int rsatFineCal;
	extern int S1Cal;
	extern int noPower;
	int index;
	float theta,phi;
	double ReH, aRange,drA;
	if(*value > 0) {
		if( inputImage->rAnt !=NULL && noPower > 0) {
			/* use antenna pattern */
			aRange=range*inputImage->rangePixelSize + inputImage->cpAll.RNear;
			drA=(inputImage->rAnt[1]-inputImage->rAnt[0]);
			index= (aRange-inputImage->rAnt[0]) /drA; 
			index=min(max(0,index),inputImage->patSize);
			*value=*value*inputImage->pAnt[index];
		} else {
			aRange=range*inputImage->rangePixelSize + inputImage->cpAll.RNear;
			ReH=getReH( &(inputImage->cpAll), inputImage, azimuth);
			theta = thetaRReZReH(aRange,(inputImage->cpAll.Re+h),ReH) *RTOD;
			phi=psiRReZReH(aRange,(inputImage->cpAll.Re+h),ReH);
			if( (inputImage->patSize  < 0 && noPower < 0) ) {
				/*   use polynomial function - radarsat or alos specific 
				     At present this only does fine beam RS, but function could easily be modified for other antenna ppatters
				     calibrated  added 11/18/2013
				*/
				if(rsatFineCal == TRUE) {
					*value *= 1.0/polypat(theta);
					/* ADDED 3/4/15 to reference sigma nought  to ellipsoid
					   removed because it made more stripy  *value*= sin(phi)/sin(inputImage->geoData.phic * DTOR);
					   moved to final part           	   *value *=1.0/27.3; 
					   *value-=0.0058;*/
					/* Set noise floor at -30db*/
					/*	      if(*value < 0.001) *value= 0.001;*/
				} else {
					if(inputImage->patSize == RSATFINE) *value *= 1.0/polypat(theta);
					if(inputImage->patSize == ALOS){ *value *= 1.0/polyALOS(theta);}
					/* ADDED 3/4/15 to reference sigma nought  to ellipsoid
					   *value*= sin(phi)/sin(inputImage->geoData.phic * DTOR);
					   Removed because it caused additional striping.
					*/
				} 
			} else if(S1Cal == TRUE) {
				if (inputImage->betaNought > 0) {
					*value *= sin(phi) / pow(inputImage->betaNought,2.0);
				} else {
					fprintf(stderr,"input image %s",inputImage->file);
					error("S1 Cal, but beta Nought < 0 ");
				}
			} else {
				/*
				  This equation is for uncalibrated results - primarily sentinel. The first sin(phi) scales
				  for sigma nought. The sin(phi)**2 appears to cosmetically correct for the incidence angle
				*/
				*value *= sin(phi) * sin(phi)*sin(phi);
				/*	attempt at lambertian backscatter *value *= cos(37.*DTOR)*cos(37.*DTOR)*sin(phi)/(cos(phi)*cos(phi)); */
			}
		}
	}
}


static void getGeoMosaicImage(inputImageStructure *inputImage, int *imageDate,int smoothL) {
	extern int nearestDate;
	float **tmpImage;	
	int extraPad, tmpi;
	int doy[12]={0,31,59,90,120,151,181,212,243,273,304,333};	
	int k1,k2, i;
	initllToImageNew(inputImage); /* Setup conversions */
	*imageDate= inputImage->year*365. + doy[inputImage->month-1] + inputImage->day;
	if(nearestDate > 0) fprintf(stderr,"Image date %i, Nearest Date %i %i %i %i\n",
				    *imageDate,nearestDate,inputImage->year,inputImage->month,inputImage->day);
	getMosaicInputImage( inputImage);
	/* Multilook the image */
	if(smoothL > 0) smoothImage( &(inputImage[i]),smoothL);
	if(inputImage->removePad > 0) {
		extraPad=0;
		tmpi=inputImage->rangeSize * inputImage->nRangeLooks;
		/* 1/23/06 This removes a little extra on wider images (should mostly apply to fine 
		   beam. This should handle any range shifts */
		if(tmpi > 8500) extraPad=120/inputImage->nRangeLooks;
		/* Pad around edges */
		for(k2=0; k2 < inputImage->azimuthSize; k2++)
			for(k1=0; k1 < (inputImage->removePad + extraPad); k1++) {
				tmpImage=(float **)inputImage->image;
				tmpImage[k2][k1]=0.;
				tmpImage[k2][inputImage->rangeSize-k1-1]=0.;
			}
	}
}

static void initImageBuffers(outputImageStructure *outputImage, float **scaleTmp, int orbitPriority) {
	float **image, **scale;
	int i1,j1;
	image = (float **)outputImage->image;
	scale = (float **)outputImage->scale;
	for (i1=0; i1 < outputImage->ySize; i1++) {
		for(j1=0; j1 < outputImage->xSize; j1++) {
			scale[i1][j1]=0.0;
			image[i1][j1]=0.0;
			scaleTmp[i1][j1]=0.0;
			if(orbitPriority > -1) scale[i1][j1]=-1;
		}
	}
}

static void mallocTmpBuffers(float ***imageTmp, float ***scaleTmp,  outputImageStructure *outputImage) {
	size_t bufSize;   
	float *bufI,*bufS;
	int i;
	bufSize = outputImage->xSize * outputImage->ySize *sizeof(float);
	bufI=(float *)malloc((size_t)bufSize);
	bufS=(float *)malloc((size_t)bufSize);
	bufSize=sizeof(float *)*outputImage->ySize;
	*imageTmp=(float **)malloc((size_t)bufSize); 
	*scaleTmp=(float **)malloc((size_t)bufSize);
	for(i=0; i < outputImage->ySize; i++) {
		(*imageTmp)[i] = (float *) &(bufI[i*outputImage->xSize]);
		(*scaleTmp)[i] = (float *) &(bufS[i*outputImage->xSize]);
	}
}


/* final cal for rsat */
static void rsatFinalCal(float **image,  int xSize, int ySize) {
	long int i1,j1;
	for (i1=0; i1 < ySize; i1++) {
		for(j1=0; j1 < xSize; j1++) {
			image[i1][j1] *=1.0/27.3; 
			image[i1][j1] -=0.0058;
		}
	}
}

/* compute 10log(sig) for calibrated results */
static void logSigma(float **image,  int xSize, int ySize) {
	long int i1,j1;
	for (i1=0; i1 < ySize; i1++) {
		for(j1=0; j1 < xSize; j1++) {
			if(image[i1][j1] < 0.001) image[i1][j1]= 0.001;
			image[i1][j1]=10.0*log10(image[i1][j1]);
		}
	}
}



/*
  On 3/21/2013 this just does F1 RADARSAT with a 9th degree polynomial that provides 
  a tight match to the antenna pattern provided by ASF.

  The pattern is in dB, but the gain is return as a look-angle dependent scale factor.
*/
static float polypat(float theta)
{
	double thetaCentered;
	double logGain;
	float gain;
	thetaCentered = (theta -33.399360)/1.880050;
	logGain=0.0;
	logGain +=0.044834  *pow(thetaCentered,9);
	logGain +=-0.048875 *pow(thetaCentered,8);
	logGain +=-0.338039 *pow(thetaCentered,7);
	logGain +=0.414746  *pow(thetaCentered,6);
	logGain +=0.767037  *pow(thetaCentered,5);
	logGain +=-2.244836 *pow(thetaCentered,4);
	logGain +=-1.038396 *pow(thetaCentered,3);
	logGain +=1.074398  *pow(thetaCentered,2);
	logGain +=0.935996  *pow(thetaCentered,1);
	logGain +=1.707515;
	gain=(float)pow(10.0,logGain/10.0);
	return gain;  
}



static float polyALOS(float theta)
{
	double thetaCentered;
	float gain;
	thetaCentered = (theta -34.0);
	gain=0.0;
	gain +=0.00006896  *pow(thetaCentered,4);
	gain +=0.00011602  *pow(thetaCentered,3);
	gain +=-0.05117270 *pow(thetaCentered,2);
	gain +=0.00516289  *pow(thetaCentered,1);
	gain +=1.00243524;
	gain=gain*gain; /* tests indicate this is 2way */
	return gain;
}

 
static float interpolatePowerInputImage(inputImageStructure inputImage,
					double range,double azimuth)
{   
	float result;
	float t,u;
	float p1,p2,p3,p4;
	float **fimage;
	float minvalue;
	int i,j;
	/*
	  If out of input image, return 0
	*/
	if(range < 0.0 || azimuth < 0.0 ||
	   range >= (inputImage.rangeSize) || 
	   azimuth >= (inputImage.azimuthSize) ) return 0.0;
	fimage = (float **)inputImage.image;

	j = (int)range;
	i = (int)azimuth;
	/*
	  Handle border pixels
	*/
	if(j ==  (inputImage.rangeSize-1) ||
	   i == (int) (inputImage.azimuthSize-1) ) {
		return fimage[i][j];
	}
	t = (float)(range - (double)j);
	u = (float)(azimuth - (double)i);
	p1 = fimage[i][j];
	p2 = fimage[i][j+1];
	p3 = fimage[i+1][j+1];
	p4 = fimage[i+1][j];
	/* return largest value if any of the pixels is zero */
	if(p1 == 0 || p2==0 || p3==0 || p4==0) return(max(max(p1,p2),max(p3,p4)));
	minvalue = (float) -LARGEINT * 0.99999;
	if(p1 <= minvalue || p2 <= minvalue || p3 <= minvalue || p4 <= minvalue)
		return max(max(max(p1,p2),p3),p4);
	result = (float)((1.0 - t)*(1.0 - u) * p1 + t * (1.0 - u) * p2 + t * u * p3 +                (1.0 - t) * u* p4);
	return result;
}


static void smoothImage( inputImageStructure *inputImage, int smoothL)
{
	extern float *smoothBuf;
	float **image;
	float *filt;
	float sum, sumF;
	int i,j,k; 
	int hw;
	if(smoothL < 1) return;
	if(smoothL == 1) hw=1; else hw=(int) (smoothL/2);

	filt=(float *)malloc( (size_t) ( (hw*2+1)*sizeof(float) ) );
	filt = &(filt[hw]);
	if(smoothL > 1)  for(i=-hw; i <= hw; i++) filt[i]=1;
	else {filt[-1]=0.5; filt[0]=1.0; filt[1]=0.5;}

	sumF=0;
	for(i=-hw; i <= hw; i++) {
		fprintf(stderr,"filt %f\n",filt[i]);
		sumF+=filt[i];
	}
	
	fprintf(stderr,"hw = %i %f\n",hw,sum);
	if(smoothBuf==NULL) error("memory not allocated for smoothImage\n");

	image=(float **)inputImage->image; 
	/*
	  smooth in range direction - note borders can be irregular so be careful not to smooth across data/zero transition
	*/
	for(i=0; i < inputImage->azimuthSize; i++)  {
		/* Smooth line and save in buf */
		for(j=hw; j < (inputImage->rangeSize-hw); j++) { 
			smoothBuf[j]=0;
			sum=0;
			for(k=-hw; k<= hw; k++)  {
				if(image[i][j+k] > 1.e-6) {
					smoothBuf[j]+=image[i][j+k] * filt[k]; 					
					sum+=filt[k];
				}
			}
			if(sum < sumF ) smoothBuf[j]=0.0;
			else smoothBuf[j] /= sum;
		}
		/* Write result back to main buf */
		for(j=hw; j < inputImage->rangeSize-hw; j++) image[i][j]=smoothBuf[j];
	}
	fprintf(stderr,"-- %i %i \n",i,j);
	/*
	  smooth in azimuth direction
	*/
	for(j=0; j < inputImage->rangeSize; j++)  {
		/* Smooth line and save in buf */
		for(i=hw; i < (inputImage->azimuthSize-hw); i++) { 
			smoothBuf[i]=0;
			sum=0.0;
			for(k=-hw; k<= hw; k++)  {
				if(image[i+k][j]  > 1e-6) {
					sum+=filt[k];
					smoothBuf[i]+=image[i+k][j] * filt[k];
				}
			}
			if(sum < sumF ) smoothBuf[i] = 0.0;
			else smoothBuf[i] /= sum;			
		}
		/* Write result back to main buf */
		for(i=hw; i < inputImage->azimuthSize-hw; i++) image[i][j]=smoothBuf[i];
	}
	fprintf(stderr,"++ %i %i \n",i,j);
}


static void freeConversion(inputImageStructure *image)
{
	int i,j;

	for(i=0; i < image->nAzGrid; i++) {
		for(j=0; j < image->nRGrid; j++) {
			/*         
				   free(image->conversionData[i][j].rNear);
				   free(image->conversionData[i][j].rFar);
				   free(image->conversionData[i][j].ReH);
			*/
			free(image->conversionData[i][j].a2[0]);
			free(image->conversionData[i][j].a2);
		}
		/* free(image->conversionData[i]);*/
	}
	/*
	  free(image->conversionData);
	  free(image->cpAll.rNear);
	  free(image->cpAll.rFar);
	  free(image->cpAll.ReH);
	*/
	free(image->cpAll.a2[0]);
	free(image->cpAll.a2);


}
