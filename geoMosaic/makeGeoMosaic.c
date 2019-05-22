#include "stdio.h"
#include"string.h"
#include <math.h>
#include "mosaicSource/common/common.h"
#include "geomosaic.h"
#include "cRecipes/nrutil.h"
#include <stdlib.h>
#include "time.h"
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
/* from lsmosaic, but only need this */
double xyscale(double latctr,int32 proj);
/*
************************ Geocode image InSAR DEM. **************************
*/
int sepAscDesc;
static void readDEMInput(inputImageStructure *inputImage, char *insarDEMFile);
static float interpolatePowerInputImage(inputImageStructure inputImage, double range,double azimuth);
static float polypat(float theta);
static float polyALOS(float theta);
static void smoothImage( inputImageStructure *inputImage, int smoothL);
static void rsatFinalCal(float **image, int xSize, int ySize);
static void logSigma(float **image,  int xSize, int ySize) ;
static void mallocTmpBuffers(float ***imageTmp, float ***scaleTmp,  outputImageStructure *outputImage,float ***psiBuf, float ***psiBufTmp,float ***gBuf, float ***gBufTmp);
static void getGeoMosaicImage(inputImageStructure *inputImage, int *imageDate, int smoothL);
static float applyCorrections(float *value,inputImageStructure *inputImage, double range, double azimuth, double h);
static void geoMosaicScaling(inputImageStructure *inputImage, float **image, float **imageTmp,float **psiBuf, float **psiBufTmp,float **gBuf, float **gBufTmp,
			     float **scale, float **scaleTmp, void *dem, outputImageStructure *outputImage, int orbitPriority,int imageDate, int iMin,int iMax, int jMin,int jMax);
static void finalReScale(outputImageStructure *outputImage, float **image, float **scale, int orbitPriority);
static void initImageBuffers(outputImageStructure *outputImage, float **scaleTmp, int orbitPriority, float **psiBuf, float **psiBufTmp ,float **gBuf, float **gBufTmp );


static double edgeLength(double x1,double y1,double z1,double x2,double y2,double z2) {
	return sqrt( pow(x1-x2,2) + pow(y1-y2,2) +  pow(z1-z2,2) );
}
static double threeDTriangleArea(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3) {
	double p12,p23,p31,h1, A;
	p12=edgeLength(x1,y1,z1,x2,y2,z2);
	p23=edgeLength(x2,y2,z2,x3,y3,z3);
	p31=edgeLength(x3,y3,z3,x1,y1,z1);
	h1=0.5 * (p12 + p23 + p31);
	A =sqrt(h1 * (h1-p12)*(h1-p23)*(h1-p31));
	return A;
}
static void normVector(double x0, double y0, double z0, double *nx, double *ny,double *nz) {
	double den;
	den=sqrt(pow(x0,2.) +pow(y0,2.) + pow(z0,2.));
	*nx =x0/den;
	*ny =y0/den;
	*nz =z0/den;	
}
/*
  find the normal to three points.
*/
static void threePtNormal(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double *nx,double *ny,double *nz) {
	cross(x2-x1,y2-y1,z2-z1,x3-x1,y3-y1,z3-z1,nx,ny,nz);
	normVector(*nx,*ny,*nz,nx,ny,nz);
}

static void satNorm(double x0, double y0, double z0, double azimuth,inputImageStructure *inputImage, double *nlx, double *nly,double *nlz, double *nsx,double *nsy,double *nsz) {
	conversionDataStructure *cp;
	double vsx,vsy,vsz;
	double xs,ys,zs;
	double myTime;
	double den;
	cp = &(inputImage->cpAll);
	myTime=(azimuth* inputImage->nAzimuthLooks)/cp->prf + cp->sTime;
	getState( myTime,inputImage, &xs,&ys, &zs, &vsx, &vsy, &vsz);
	normVector(x0-xs,y0-ys,z0-zs,nlx,nly,nlz);
	normVector(vsx,vsy,vsz,nsx,nsy,nsz);
}

static void earthNorm(double x0, double y0, double z0, double *nx, double *ny,double *nz) {
	double den;
	normVector(x0,y0,z0,nx,ny,nz);
}


static void rangePlane(double nex,double ney,double nez,double nsx,double nsy, double nsz, double *nrx,double *nry,double *nrz) {
	double den;
	cross(nex,ney,nez,nsx,nsy,nsz,nrx,nry,nrz);
	normVector(*nrx,*nry,*nrz, nrx,nry,nrz);
}

		       
static double threeDArea(double x1[4],double y1[4], double z1[4]) {
	double A1, A2;
	A1=threeDTriangleArea(x1[0],y1[0],z1[0],x1[1],y1[1],z1[1],x1[3],y1[3],z1[3]);
	A2=threeDTriangleArea(x1[1],y1[1],z1[1],x1[2],y1[2],z1[2],x1[3],y1[3],z1[3]);
	return (A1+A2);	
}
static	double nrx,nry,nrz, nsx,nsy,nsz,nlx,nly,nlz;
static void AbAg(double x, double y,double azimuth, double psiE,inputImageStructure *inputImage,  outputImageStructure *outputImage, void *dem,double *Ab, double *Ag, double *shadow, int recycle) {
	extern int HemiSphere;
	extern double Rotation;
	double xc[4],yc[4],z[4];
	double latc[4],lonc[4];
	double xi[4]={-0.5,0.5,0.5,-0.5}; /* Use +/- for a smoother result */
	double yi[4]={-0.5,-0.5,0.5,0.5};
	double x1[4],y1[4],z1[4];
	double x1g[4],y1g[4],z1g[4];
	double x1b[4],y1b[4],z1b[4];			
	double xp[4],yp[4],zp4[4];
	double nex,ney,nez,nex1,ney1,nez1,nex2,ney2,nez2;
	double vdotn,rdotn;
	double x0,y0,z0;
	int i;
	/* Map pixel corners to ECEF */
	x0=0.0; y0=0.0; z0=0.0;
	for(i=0; i < 4; i++) {
		xc[i]=x+xi[i]*outputImage->deltaX *MTOKM;
		yc[i]=y+yi[i]*outputImage->deltaY * MTOKM;
		xytoll1(xc[i],yc[i],HemiSphere,&(latc[i]),&(lonc[i]),Rotation,outputImage->slat);
		z[i]=getXYHeight(latc[i],lonc[i],dem,inputImage->cpAll.Re,ELLIPSOIDAL);
		llToECEF( latc[i],lonc[i],z[i], &(x1[i]), &(y1[i]), &(z1[i]));
		x0 += 0.25*x1[i];   y0 += 0.25*y1[i]; 	z0 += 0.25*z1[i];
	}
	/* Normal to look direction  - toward surface  - gamma plane */
	if(recycle == FALSE)  {
		satNorm(x0,y0,z0,azimuth,inputImage,&nlx,&nly,&nlz,&nsx,&nsy,&nsz);
		/* Normal for plane that contains sat and look vector - beta plane */
		rangePlane(nsx,nsy,nsz, nlx, nly,  nlz, &nrx,&nry,&nrz);
	}
	/* Surface normal */
	if(*shadow > 0) {
		threePtNormal(x1[0],y1[0],z1[0],x1[1],y1[1],z1[1],x1[2],y1[2],z1[2],&nex1,&ney1,&nez1);
		threePtNormal(x1[0],y1[0],z1[0],x1[2],y1[2],z1[2],x1[3],y1[3],z1[3],&nex2,&ney2,&nez2);
		/* Shadow check */	
		*shadow = min( dot(nex1,ney1,nez1,-nlx,-nly,-nlz), dot(nex2,ney2,nez2,-nlx,-nly,-nlz) );
	}	
	/* Project points on to gamma and beta planes */
	for(i=0; i < 4; i++) {
		/*  corner points dotted with line of sight vector */
		vdotn=    dot(x1[i]-x0,y1[i]-y0,z1[i]-z0,nlx,nly,nlz);
		/* Points projected to line of sight plane */
		x1g[i] = x1[i]-vdotn*nlx;   y1g[i] = y1[i]-vdotn*nly;  z1g[i] = z1[i]-vdotn*nlz;
		/* project points on to plane with line of sight and velocity vector */
		rdotn=    dot(x1[i]-x0,y1[i]-y0,z1[i]-z0,nrx,nry,nrz);		
		x1b[i]=x1[i]-rdotn*nrx;	y1b[i]=y1[i]-rdotn*nry;	z1b[i]=z1[i]-rdotn*nrz;
	}
	/* normalization - removed since these a ratioed*/
	*Ag=threeDArea(x1g,y1g,z1g);
        *Ab=threeDArea(x1b,y1b,z1b);
}

static double areaAboutXY(double range,double azimuth,double x,double y,inputImageStructure *inputImage,  outputImageStructure *outputImage,void *dem, double psiE, double value, double *test, int recycle, double *Ab,double *Ag) 
{
	extern int HemiSphere;
	extern double Rotation;
	double  shadow, Ab1,Ag1,rat;
	double x1,y1,lat1,lon1,hWGS,range1,azimuth1;
	double r1,r2,a1,a2;
	int nx=3,ny=2;
	double dx;
	int i,j;
	float tmp;
	shadow=1;
	AbAg(x,y,azimuth,psiE,inputImage,outputImage,dem, Ab,Ag,&shadow,recycle);
	if(shadow < -0.001) return(-1.);
	/*	if((Ab/sin(psiE *DTOR)) > 0.8) return Ab/Ag;*/
	rat=*Ab/ *Ag;
	if(10.*log10(value) < -2. || rat > 0.65 ) return rat;
	/* This part is quite and inefficient way to bin ground samples contributing to the same range bin */
	dx=0.5;	
	r1=range-dx;		r2=range+dx;
	a1=azimuth-dx;	a2=azimuth+dx;
	recycle=TRUE;
	for(i=-nx; i <= nx; i++ ) {
		x1=x+i*outputImage->deltaX *MTOKM;
		for(j=-ny; j <= ny; j++ ) {
			y1=y+j*outputImage->deltaY * MTOKM;
			xytoll1(x1,y1,HemiSphere,&lat1,&lon1,Rotation,outputImage->slat);
			hWGS=getXYHeight(lat1,lon1,dem,inputImage->cpAll.Re,ELLIPSOIDAL);			
			llToImageNew(lat1,lon1,hWGS,&range1,&azimuth1,inputImage);
			if(range1 > r1 && range1 < r2 && azimuth1 > a1 && azimuth1  < a2 && ( !((i ==0) && (j ==0))  )  ) {
				shadow=-1;
				AbAg(x1,y1,azimuth1,psiE,inputImage,outputImage,dem, &Ab1,&Ag1,&shadow,recycle);	
				*Ab+=Ab1;
				*Ag+=Ag1;
			}
		}
	}
	return (*Ab)/(*Ag);
}

static void oversampledxdy(int os, double dx[3],double dy[3]) {
	if(os ==1)  {
		dx[0]=0;	dy[0]=0;
	} else if(os == 2 ) {
		dx[0]=-0.5;	dy[0]=-0.5;
		dx[1]=0.5;	dy[1]=0.5;
	} else if(os=3) {
		dx[0]=-1;		dy[0]=-1;
		dx[1]=0;		dy[1]=0;
		dx[2]=1;		dy[2]=1;				
	} else error("Invalid oversample - should be between 1 and 3");
}

/*
  Make mosaic of insar and other dems.
*/
void makeGeoMosaic(inputImageStructure *inputImage,   outputImageStructure outputImage,
		   void *dem, int nFiles, int maxR, int maxA,   char **imageFiles,float fl, int smoothL,int smoothOut, int orbitPriority, float ***psiData,float ***gamma)
{   
	extern int HemiSphere;
	extern double Rotation;
	extern int nearestDate;
	extern int hybridZ; 
	extern int rsatFineCal; 
	extern int noPower;
	extern int S1Cal;
	double psiE;
	float **scale, **psiBuf, **psiBufTmp,**gBufTmp,**gBuf;
	FILE *fp;
	double x,y, lat,lon;
	double range,azimuth;
	double h,hWGS;
	float value,valueTmp;
	float **imageTmp, **scaleTmp,**image;
	int imageDate;
	int i,j,i1,j1;
	int iMin,iMax,jMin,jMax;
	double RTCscale,Ab,Ag,AbCum,AgCum;
	double test;
	int ios,jos;
	double dx[3],dy[3],invNAvg;
	int recycle, shadow;
	clock_t startTime,lastTime,initTime;	
	oversampledxdy(smoothOut,  dx,dy);
	invNAvg=1.0/(smoothOut*smoothOut);
	/*
	  Malloc space for tmp buffers
	*/
	mallocTmpBuffers(&imageTmp, &scaleTmp, &outputImage,&psiBuf,&psiBufTmp,&gBuf,&gBufTmp);
	*psiData=psiBuf; *gamma=gBuf;
	/* Init image,scale */	
	initImageBuffers(&outputImage,scaleTmp, orbitPriority, psiBuf,psiBufTmp,gBuf,gBufTmp);
	image = (float **)outputImage.image;
	scale = (float **)outputImage.scale;
	fprintf(stderr,"\norbitPriority %i\n",orbitPriority);
	/*
	  Main mosaicing loop
	*/
	fprintf(stderr,"STARTING MOSAICING\n");
	startTime=clock();
	for(i=0; i < nFiles; i++) {
		/*
		  Process input image
		*/
		inputImage[i].file=imageFiles[i];
		/*  Get bounding box of image */
		getRegion(&(inputImage[i]),&iMin,&iMax,&jMin,&jMax,&outputImage);
		fprintf(stderr,"*** iMin,iMax, jMin,jMax ++ w %i %i %i %i ++ %f\n",iMin,iMax,jMin,jMax,inputImage[i].weight );
		if(inputImage[i].weight < 0.0001) { /* Zero weight so force skip */
			iMin=1; iMax=-1;  jMin=1; jMax=-1; 
		} else if(iMax > iMin && jMax > jMin) { /* Read image if inbounds */
			getGeoMosaicImage(&(inputImage[i]), &imageDate,  smoothL);
		}
		/*
		  Loop over output grid
		*/
		fprintf(stderr,"%s\n",imageFiles[i]);
		for (i1=iMin; i1 < iMax; i1++) {
			if( (i1 % 100) == 0) fprintf(stderr,"-- %i\n",i1);  

			for(j1=jMin; j1 < jMax; j1++) {
				value=0;
				psiE=0;
				RTCscale=0.0;
				recycle=FALSE;
				AbCum=0.0; AgCum=0.0;
				shadow=FALSE;
				for(ios=0; ios < smoothOut; ios++) {
					for(jos=0; jos < smoothOut; jos++) {
						y = (outputImage.originY + (i1+dy[ios])*outputImage.deltaY ) * MTOKM;
						x = (outputImage.originX + (j1+dx[jos])*outputImage.deltaX ) * MTOKM ;						
						/*
						  Convert x/y stereographic coords to lat/lon
						*/
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
						valueTmp=interpolatePowerInputImage(inputImage[i],range,azimuth);
						/* bad value so punt here */
						if(valueTmp < 0.0) { break; }
						psiE += applyCorrections(&valueTmp,&(inputImage[i]),range, azimuth,h) * invNAvg;

						if( S1Cal == TRUE) {
							if( areaAboutXY(range,azimuth,x,y,&(inputImage[i]),&outputImage,dem,psiE,valueTmp,&test,recycle, &Ab,&Ag) < -0.001) shadow=TRUE;
							AbCum+=Ab; AgCum+=Ag;
							recycle=TRUE;
						}
						value+=valueTmp*invNAvg;
					}
					/* bad value so punt here */					
					if(valueTmp < 0.0) { break; }					
				}
				if(value > 0 ) {
					psiBufTmp[i1][j1]=psiE;
					/* Note the sin(psiE) undoes the psiE for sigma nought */
					if(shadow == FALSE && S1Cal==TRUE) { 
						gBufTmp[i1][j1]=10.0*log10((AbCum/AgCum)/sin(psiE*DTOR));						
						gBufTmp[i1][j1]=round(gBufTmp[i1][j1]*100.)/100.;
						gBufTmp[i1][j1]=min(max(gBufTmp[i1][j1],-29.9),20.0);
					} else gBufTmp[i1][j1]=-30.0; /* Negative value indicates shadow */
				}  else value=-LARGEINT;	  
				/* Scaling */
				imageTmp[i1][j1] =value;   /*DEBUG value;*/
				if(value > 0 || (noPower >0 && value > inputImage->noData) ) {
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
		geoMosaicScaling(&(inputImage[i]), image,imageTmp,psiBuf, psiBufTmp,gBuf,gBufTmp,scale, scaleTmp, dem, &outputImage, orbitPriority,imageDate,
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
	lastTime=clock();
	fprintf(stderr,"totalTime %f\n",(double)(lastTime-startTime)/CLOCKS_PER_SEC );	
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

static void geoMosaicScaling(inputImageStructure *inputImage, float **image, float **imageTmp,float **psiBuf, float **psiBufTmp,float **gBuf, float **gBufTmp,
			     float **scale, float **scaleTmp, void *dem, outputImageStructure *outputImage, int orbitPriority,int imageDate, int iMin,int iMax, int jMin,int jMax){
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

			if(imageTmp[i1][j1] > 0 || (noPower > 0 && imageTmp[i1][j1] > inputImage->noData) ) { /* Points with valid datat */
				/* case for no nearestDate, no orbitPriority, or hWGS override */
				if((nearestDate <  0  || (nearestDate > 0 && (int)hWGS > hybridZ)) && orbitPriority < 0 ) { /* Summing data or non-nearest date or hybridZ*/
					image[i1][j1] += imageTmp[i1][j1] * scaleTmp[i1][j1]*inputImage->weight;
					scale[i1][j1] += scaleTmp[i1][j1];
					psiBuf[i1][j1] = psiBufTmp[i1][j1];
					gBuf[i1][j1] = gBufTmp[i1][j1];					
				} else if(orbitPriority > -1 ) { /* orbit priority case */
					if(orbitPriority == ASCENDING) {
						/* either put in data if none already, or if not ascending replace */					
						if( scale[i1][j1] > 0.1 || scale[i1][j1] ==DESCENDING) {							
							image[i1][j1]= imageTmp[i1][j1]*inputImage->weight;
							scale[i1][j1] = inputImage->passType;
							psiBuf[i1][j1] = psiBufTmp[i1][j1];
							gBuf[i1][j1] = gBufTmp[i1][j1];							
						}
					} else if(orbitPriority == DESCENDING) {
						/* either put in data if none already, or if not descending replace */
						if( scale[i1][j1] < -0.1 || scale[i1][j1] == ASCENDING) {
							image[i1][j1]= imageTmp[i1][j1]*inputImage->weight;
							scale[i1][j1] = inputImage->passType;
							psiBuf[i1][j1] = psiBufTmp[i1][j1];
							gBuf[i1][j1] = gBufTmp[i1][j1];														
						}							
					} else error("bad orbit priority flag");
				} else { /* Nearest date, if orbitPriority not st */
					if( (abs(nearestDate -imageDate) < abs(nearestDate-scale[i1][j1])) && scaleTmp[i1][j1] > 0 ) {
						image[i1][j1]= imageTmp[i1][j1]*inputImage->weight;
						scale[i1][j1] = imageDate;
						psiBuf[i1][j1] = psiBufTmp[i1][j1];
						gBuf[i1][j1] = gBufTmp[i1][j1];												
					}
				} /* end else */
			} /* end if imageTmp...*/
			scaleTmp[i1][j1] = 0.0;
		} /* end j1 */
	} /* end i1=iMin sum current ... */
}


static float applyCorrections(float *value, inputImageStructure *inputImage, double range, double azimuth, double h) {
	extern int rsatFineCal;
	extern int S1Cal;
	extern int noPower;
	int index;
	float theta,psi;
	double ReH, aRange,drA;
	psi=0.0;
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
			psi=psiRReZReH(aRange,(inputImage->cpAll.Re+h),ReH);
			if( (inputImage->patSize  < 0 && noPower < 0) ) {
				/*   use polynomial function - radarsat or alos specific 
				     At present this only does fine beam RS, but function could easily be modified for other antenna ppatters
				     calibrated  added 11/18/2013
				*/
				if(rsatFineCal == TRUE) {
					*value *= 1.0/polypat(theta);
					/* ADDED 3/4/15 to reference sigma nought  to ellipsoid
					   removed because it made more stripy  *value*= sin(psi)/sin(inputImage->geoData.psic * DTOR);
					   moved to final part           	   *value *=1.0/27.3; 
					   *value-=0.0058;*/
					/* Set noise floor at -30db*/
					/*	      if(*value < 0.001) *value= 0.001;*/
				} else {
					if(inputImage->patSize == RSATFINE) *value *= 1.0/polypat(theta);
					if(inputImage->patSize == ALOS){ *value *= 1.0/polyALOS(theta);}
					/* ADDED 3/4/15 to reference sigma nought  to ellipsoid
					 *value*= sin(psi)/sin(inputImage->geoData.psic * DTOR);
					 Removed because it caused additional striping.
					*/
				} 
			} else if(S1Cal == TRUE) {
				if (inputImage->betaNought > 0) {
					*value *= sin(psi) / pow(inputImage->betaNought,2.0);
				} else {
					fprintf(stderr,"input image %s",inputImage->file);
					error("S1 Cal, but beta Nought < 0 ");
				}
			} else {
				/*
				  This equation is for uncalibrated results - primarily sentinel. The first sin(psi) scales
				  for sigma nought. The sin(psi)**2 appears to cosmetically correct for the incidence angle
				*/
				if(noPower < 0)
					*value *= sin(psi) * sin(psi)*sin(psi);
				/*	attempt at lambertian backscatter *value *= cos(37.*DTOR)*cos(37.*DTOR)*sin(psi)/(cos(psi)*cos(psi)); */
			}
		}
	}
	return (float)(psi* RTOD) ;
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

static void initImageBuffers(outputImageStructure *outputImage, float **scaleTmp, int orbitPriority, float **psiBuf, float **psiBufTmp,float **gBuf, float **gBufTmp )
{
	float **image, **scale;
	int i1,j1;
	image = (float **)outputImage->image;
	scale = (float **)outputImage->scale;

	for (i1=0; i1 < outputImage->ySize; i1++) {
		for(j1=0; j1 < outputImage->xSize; j1++) {
			scale[i1][j1]=0.0;
			image[i1][j1]=0.0;
			scaleTmp[i1][j1]=0.0;
			psiBuf[i1][j1]=0.0;
			psiBufTmp[i1][j1]=0.0;
			gBuf[i1][j1]=-30.0;
			gBufTmp[i1][j1]=-30.0;			
			if(orbitPriority > -1) scale[i1][j1]=-1;
		}
	}
}

static void mallocTmpBuffers(float ***imageTmp, float ***scaleTmp,  outputImageStructure *outputImage,float ***psiBuf, float ***psiBufTmp,float ***gBuf, float ***gBufTmp)
{
	size_t bufSize;   
	float *bufI,*bufS, *bufP, *bufPTmp, *bufG,*bufGTmp;
	int i;
	bufSize = outputImage->xSize * outputImage->ySize *sizeof(float);
	bufI=(float *)malloc((size_t)bufSize);
	if(bufI == NULL) error("Unable to malloc buffer 1 in mallocTmpBuffers\n");		
	bufS=(float *)malloc((size_t)bufSize);
	if(bufS == NULL) error("Unable to malloc buffer 2 in mallocTmpBuffers\n");			
	bufP=(float *)malloc((size_t)bufSize);
	if(bufP == NULL) error("Unable to malloc buffer 3 in mallocTmpBuffers\n");			
	bufPTmp=(float *)malloc((size_t)bufSize);
	if(bufPTmp == NULL) error("Unable to malloc buffer 4 in mallocTmpBuffers\n");
	bufG=(float *)malloc((size_t)bufSize);
	if(bufG == NULL) error("Unable to malloc buffer 5 in mallocTmpBuffers\n");			
	bufGTmp=(float *)malloc((size_t)bufSize);
	if(bufGTmp == NULL) error("Unable to malloc buffer 6 in mallocTmpBuffers\n");
	bufSize=sizeof(float *)*outputImage->ySize;
	*imageTmp=(float **)malloc((size_t)bufSize);
	*scaleTmp=(float **)malloc((size_t)bufSize);
        *psiBuf=(float **)malloc((size_t)bufSize);
        *psiBufTmp=(float **)malloc((size_t)bufSize);
        *gBuf=(float **)malloc((size_t)bufSize);
        *gBufTmp=(float **)malloc((size_t)bufSize);	

	for(i=0; i < outputImage->ySize; i++) {
		(*imageTmp)[i] = (float *) &(bufI[i*outputImage->xSize]);
		(*scaleTmp)[i] = (float *) &(bufS[i*outputImage->xSize]);
		(*psiBuf)[i] = (float *) &(bufP[i*outputImage->xSize]);
		(*psiBufTmp)[i] = (float *) &(bufPTmp[i*outputImage->xSize]);
		(*gBuf)[i] = (float *) &(bufG[i*outputImage->xSize]);
		(*gBufTmp)[i] = (float *) &(bufGTmp[i*outputImage->xSize]);						
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

/* compute 10log(sig) for calibrated results. Clip data to -29.9 */
static void logSigma(float **image,  int xSize, int ySize) {
	long int i1,j1;
	for (i1=0; i1 < ySize; i1++) {
		for(j1=0; j1 < xSize; j1++) {
			/* Limit small values to 0.001001 = -29.92 db , for no data for -30.0 */
			if( image[i1][j1] < 0.001001 ) {
				if( image[i1][j1]  > 1.e-9 ) image[i1][j1]= 0.00102;
				else  image[i1][j1]= 0.0010;
			}
			/* Round to 1/100 of dB - in conversion to tif, it may get rounded further */
			image[i1][j1]=roundf( (10.0*log10(image[i1][j1])) * 100)/100.;
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

 
static float interpolatePowerInputImage(inputImageStructure inputImage,double range,double azimuth)
{   
	float **fimage;
	int i,j;
	/*
	fimage=(float **)inputImage.image;
	j = (int)round(range);
	i = (int)round(azimuth);
	if(range < 0.0 || azimuth < 0.0  ||  j >= inputImage.rangeSize || i >= inputImage.azimuthSize ) return inputImage.noData;
	return fimage[i][j];
	*/
	fimage=(float **)inputImage.image;
	return  bilinearInterp(fimage,range,azimuth,inputImage.rangeSize,inputImage.azimuthSize,inputImage.noData, inputImage.noData);
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


