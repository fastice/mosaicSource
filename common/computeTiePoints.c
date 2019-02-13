#include <math.h>
#include <stdlib.h>
#include "string.h"
#include "source/common/common.h"
/*#include "tiePoints.h"*/
#include "time.h"
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
extern double Rotation;

/*
  Compute location of  tie points for baselines estimation.
*/



void computeTiePoints(inputImageStructure *inputImage,      tiePointsStructure *tiePoints,
		      demStructure dem, int noDEM, char *geodatFile,ShelfMask *shelfMask,int supressOutput)
{   
	FILE *fp;   
	int i,j;
	conversionDataStructure *cP;
	double x,y,zSp,zWGS;
	double xx,yy;
	double lat,lon;
	double azOrigin;
	double num,den,rho,r;
	double RElip;
	double  Re;
	double range,azimuth;
	double deltaZ;
	double vz1;
	int sMask;
	int32 shelfMaskFlag;
	char tideFile[2048],*tmpS;
	ers1Complex **cimage;
	float **fimage;
	clock_t startTime,lastTime,initTime;
	int k=0;
	startTime=clock();
	initllToImageNew( inputImage);
	fprintf(stderr,"Using new \n");
	initTime=clock();
	cP = &(inputImage->cpAll);
	Re = cP->Re;
	shelfMaskFlag = (int32) shelfMask;
	
	if(supressOutput==FALSE) {
		fprintf(stdout,";;\n;; InSAR parameters  Re, RNear, Rc, nlr,nla\n;;\n");  
		fprintf(stdout," %10.3f  %8.3f %8.3f  %i %i %f %f\n",
			cP->Re*MTOKM,cP->RNear*MTOKM,inputImage->par.rc*MTOKM,
			inputImage->nRangeLooks,inputImage->nAzimuthLooks,inputImage->rangePixelSize/inputImage->nRangeLooks, inputImage->par.lambda);
	}
	/*
	  Read spatially varying tide correction
	*/
	tideFile[0]='\0';
	strcpy(tideFile,geodatFile);
	fprintf(stderr,"%s\n",geodatFile);
	tmpS=strrchr(tideFile,'/');
	if(tmpS != NULL) tmpS[0]='\0';
	strcat(tideFile,"/tide.difference");
	fp=fopen(tideFile,"r");
	if(fp != NULL) {
		inputImage->tideDiffFlag =TRUE;
		readXYDEM(tideFile,&(inputImage->tideDiff));
		fprintf(stderr,";; Tide difference : %s\n",tideFile);              
		if(supressOutput==FALSE)  fprintf(stdout,";; Tide difference : %s\n",tideFile);              
	} else {
		inputImage->tideDiffFlag = FALSE;
		if(supressOutput==FALSE)  fprintf(stdout,";; Tide difference : %s\n","none");              
	}
	/*
	  Loop over output points
	*/
	for(i=0; i < tiePoints->npts; i++) {
		lat = tiePoints->lat[i];
		lon = tiePoints->lon[i];
		/* Deal with shelf mask if necessary */
		if(shelfMask != NULL && shelfMaskFlag != -1) {
			lltoxy1(lat,lon,&xx,&yy,Rotation,tiePoints->stdLat);
			sMask=getShelfMask(shelfMask,xx,yy);
			if(sMask == SHELF) {
				deltaZ=interpTideDiff(xx,yy,inputImage->tideDiff);
				if(deltaZ >  MINELEVATION) {
					vz1 = deltaZ/((double)tiePoints->nDays) *365.25;
					tiePoints->vz[i] += vz1;
				}
			}
		}
		/*
		  Get height for given lat/lon, with spherical correction
		*/
		zWGS=tiePoints->z[i];
		if(tiePoints->imageCoords == FALSE ) {
			llToImageNew(lat,lon,zWGS,&range,&azimuth,inputImage);
			RElip = earthRadius(lat*DTOR,EMINOR,EMAJOR)*KMTOM   + tiePoints->z[i];
			/*  Elevation referenced to a sphere */
			zSp = RElip - Re;
		} else {
			/* Obsolete */
			range = lat; /* Lat variable = range */			
			zSp = tiePoints->z[i]; /* Assume r/a/z already spherical */
		}
		tiePoints->z[i] = zSp;
		tiePoints->r[i] = range;
		if(tiePoints->imageCoords == FALSE) {
			tiePoints->a[i] = azimuth;
		} else {
			/* ??? */
			error("check computeTiePoints - for this obsolete line");
			azimuth = inputImage->azimuthSize-lon;
			tiePoints->a[i] = lon;
		}
		azOrigin =  -0.5 * AzimuthPixelSize * inputImage->azimuthSize * inputImage->nAzimuthLooks;
		tiePoints->x[i] = azOrigin + inputImage->nAzimuthLooks * AzimuthPixelSize*azimuth;
	}
	lastTime=clock();
	fprintf(stderr,"initTime %f convertTime %f totalTime %f\n",(double)(initTime-startTime)/CLOCKS_PER_SEC
	,(double)(lastTime-initTime)/CLOCKS_PER_SEC,(double)(lastTime-startTime)/CLOCKS_PER_SEC );
	return;
}
