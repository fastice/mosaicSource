#include "stdio.h"
#include"string.h"
#include <math.h>
#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"
/*
  Function to simulate InSAR image including both terrain and motion effects.

  01/04/00 Modified to handle block updates of Re, ReH. This is kluged a 
  little but since blocks are update by range and rg found uniquely for 
  that range it should be no problem.
*/


static void withHeight(inputImageStructure *inputImage, double *lat, double *lon,double azimuth,double rsl, double zWGS) {
	int n;
	stateV *sv;
	conversionDataStructure *cp;
	 double xs,ys,zs,vsx,vsy,vsz;
	 double myTime;
	cp=&(inputImage->cpAll);
	myTime=(azimuth* inputImage->nAzimuthLooks)/cp->prf + cp->sTime;
	sv=&(inputImage->sv);
	n=(long int)((myTime - sv->times[1])/(sv->deltaT)+.5);
	n=min(max(0,n-2), sv->nState-NUSESTATE);	
	polintVec(&(sv->times[n]), &(sv->x[n]), &(sv->y[n]),&(sv->z[n]),  &(sv->vx[n]), &(sv->vy[n]),&(sv->vz[n]), myTime, &xs,&ys ,&zs, &vsx,&vsy ,&vsz);

	smlocateZD( xs*MTOKM,  ys*MTOKM,  zs*MTOKM, vsx*MTOKM,vsy*MTOKM, vsz*MTOKM, rsl*MTOKM,lat,lon,(double)(inputImage->lookDir),zWGS*MTOKM);	

}



/************************STATIC ROUTINE DECLARATIONS**************************/

unsigned char getShelfMask(ShelfMask *imageMask,double x, double y);
/************************END STATIC ROUTINES DECLARATIONS*********************/


void simInSARimage(sceneStructure *scene, void *dem,xyVEL *xyVel)
{   
	double range;
	int recycle;
	double deltarg;
	double c3, c4;
	double R, Rsquared, phi,theta,thetaD,delta,thetaC, rhoSp,rhoWGS,RCenter, RNear, RWGS,R0,ReWGS;
	double lat,lon;
	double ReH;
	double ii,jj;
	double baselineSquared;
	double toRho, r, rg;
	double deltaToPhase;
	int azIndex,azIndexSave,rIndex,rIndexSave;
	int az;
	double drgStep;
	int iLoop, jLoop;
	double iFloat,jFloat,rgSave;
	double hSp,hWGS;
	double x1,y1;
	double H,Re;
	double rOff,aOff;
	conversionDataStructure *cp,*cpSave;
	demStructure *llDem;
	xyDEM *xyDem;
	int i,j,iIndex;
	int nIter;
	long long nIterTotal;
	float azFix;
	double dRe;
	fprintf(stderr,"SIM INSAR IMAGE\n");
	/*
	  Dem type
	*/
	xyDem = (xyDEM *)dem;		
	/*
	  Init coord conversions
	*/
	initllToImageNew(&(scene->I));
	deltaToPhase = 4.0 * PI /scene->I.par.lambda;
	Re = scene->I.cpAll.Re;
	RCenter = scene->I.cpAll.RCenter;
	ReH=scene->I.cpAll.Re + scene->I.par.H ;
	H=scene->I.par.H ;
	thetaC = acos( (RCenter * RCenter + ReH*ReH - Re*Re )/ (2.0*ReH*RCenter) );
	/* ~1/10 of a pixel  */
	deltarg=(scene->I.rangePixelSize * scene->dR)*0.25/sin(thetaC);
	/* Constants for phase computation */
	c3 =  2.0 * (Re + H);
	c4 = 2.0 * H * Re + pow(H,2.0);
	fprintf(stderr,"Re,thetaC,RNear,RFar, H %f %f %f %f %f %f %f %f %i \n ", Re,thetaC,scene->I.cpAll.RNear,scene->I.cpAll.RFar,H,deltarg,scene->I.rangePixelSize, deltarg,scene->I.nRangeLooks );
	/* Correct offsets relative to pixel size change */
	rOff= 0.5 *(scene->dR - scene->I.rangePixelSize );
	aOff=0.0 ; /* (0.5*(scene->dA - scene->I.azimuthPixelSize ));	*/
	RNear=scene->I.cpAll.RNear;
	/*
	  Loop over output grid. Use units of km.
	*/
	azFix=scene->aO;
	rgSave=-1.0;
	azIndexSave=0;
	nIterTotal=0;
	cp=&(scene->I.cpAll);
	drgStep=deltarg;
	for (iLoop=0; iLoop < scene->aSize; iLoop++) {
		iFloat=iLoop*scene->dA + azFix + aOff;
		i=(int)(iFloat + 0.5);
		/*fprintf(stderr,"iFloat %f\n",iFloat); */
	        iIndex = iLoop;
		ReH=cp->ReH[i];
		rhoSp=rhoRReZReH( RNear,(Re+0), ReH);
		rg = Re * rhoSp -100;    /* Initial ground range on spherical earth */
		range = scene->I.cpAll.RNear + scene->rO*scene->I.rangePixelSize  ;    /* Starting range - gets incremented in loop */
		scene->bn = scene->bnStart + (double)iIndex * scene->bnStep;
		scene->bp = scene->bpStart + (double)iIndex * scene->bpStep;
		if(i == 0) fprintf(stderr,"%f %f \n",scene->bn,scene->bp);
		if(i % 100 == 0) fprintf(stderr,"%i %f %lf \n",i,rg,(double)(nIterTotal/(iLoop * scene->rSize) ));
		/* Loop over range - not this program uses a ground range, rg, on a spherical earth to move across the swath. 
		   The actual ground ranges aren't correct, but the final calculation is  */
		for(jLoop=0; jLoop < scene->rSize; jLoop++) {
			jFloat=scene->rO+jLoop*scene->dR + rOff;
			j=(int)(jFloat + 0.5);
			/* Initialize iteration */
			delta = 0.0;
			nIter=0;
			recycle=FALSE;
			az = i - cp->azOff;
			az= min(max(0,az),cp->azSize-1);				
			while(1) {
				/* For rg on spherical earth, compute lat/lon for the corresponding ellipsoidal earth 
				   and return the corresponding slant range */
				R0=groundRangeToLLNew(rg,iFloat,&lat,&lon,&(scene->I),recycle);
				recycle=TRUE; /* flags to reuse previously computed values that are still valid */
				/* Compute WGS radious at lat/lon */
				ReWGS=earthRadius(lat*DTOR,EMINOR,EMAJOR)*KMTOM;
				/* Comput rho for the ellipsodal radius */
				rhoWGS=rhoRReZReH( R0,(ReWGS+0), ReH);
				/*  Get height for lat/lon referenced to the ellipsoid */
				hWGS= getXYHeight(lat,lon,xyDem,Re, ELLIPSOIDAL);					
				/*  Compute height corrected R    	*/
				R=slantRange( rhoWGS, ReWGS+hWGS,  ReH);
				/* Not the most efficient iteration, but good enough */
				drgStep=(range-R)*0.5;
				if(R > (range-0.1))  break;				
				/* Speed convergence if way out */
				rg += drgStep;    /* Increment ground range */				
				nIter++;
			} /* Endwhile */
			if(jLoop == 0) rgSave=rg;
			nIterTotal +=nIter;
			/*			if(hWGS > 1000 ) fprintf(stderr,"lat lon z az R- %f %f %f %f %f\n",lat,lon,hWGS,(double)az,R);*/
			withHeight(&(scene->I), &lat, &lon, (double)az,range, hWGS);
			/*			if(hWGS > 1000) {
							fprintf(stderr,"lat lon - %f %f\n",lat,lon);
							error("stop");}*/
			/*
			  Get displacement
			*/
			if(scene->maskFlag == TRUE && scene->toLLFlag == FALSE ) {
				lltoxy1(lat,lon,&x1,&y1,xyDem->rot,xyDem->stdLat);
				scene->image[iIndex][jLoop]=getShelfMask(scene->imageMask,x1,y1);
			} else if(scene->heightFlag == TRUE && scene->toLLFlag == FALSE ) {
				scene->image[iIndex][jLoop]=(float)(hWGS);
			} else if(scene->toLLFlag == TRUE) {
				scene->latImage[iIndex][jLoop]=lat;
				scene->lonImage[iIndex][jLoop]=lon;
				if( scene->maskFlag == TRUE ) {
					lltoxy1(lat,lon,&x1,&y1,xyDem->rot,xyDem->stdLat);
					scene->image[iIndex][jLoop]=getShelfMask(scene->imageMask,x1,y1);
				}
			} else {
				/*
				  Use h on single spherical reference for phase
				  even though locally spherical value was used for location
				*/
				hSp= getXYHeight(lat,lon,xyDem,Re, SPHERICAL);
				/* NOT VALIDATED */
				hSp= hSp + Re - scene->I.cpAll.Re;
				Re=scene->I.cpAll.Re;
				/*
				  Compute look angle for height h and range R. 
				  Note I had to make a decision here whether to use ReH or Re + H
				  to compute phases. I went with the less precise Re+H to maintain
				  internal consistency through out the code. This will should have
				  the effect of introducing and erroneous phase ramp, but since the 
				  code is consistent, it should be cancelled via the effective baseline.
				  If I ever choose to update the code to carry Re+h throughout,
				  then this line should replace the theta calculation.
				  theta = acos( ( range*range + ReH*ReH - (Re+h)*(Re+h) )/
				  (2.0*ReH*range) );
				*/
				theta = acos((range*range + c4 - 2.0*hSp * Re - hSp*hSp) / (range*c3));
				thetaD = theta - thetaC;
				/*
				  Add componenent of delta from topography in meters.
				*/
				baselineSquared = pow(scene->bn,2.0) + pow(scene->bp,2.0);
				delta += (sqrt( pow(range,2.0) -2.0*range*(scene->bn*sin(thetaD) + scene->bp*cos(thetaD)) +
						baselineSquared ) - range); 
				/* 
				   Flatten phase image if flatFlag set. Include possible baseline estimation   error.
				*/
				if(scene->flatFlag == TRUE) {
					theta = acos( (range*range + c4)/(range*c3) );
					thetaD = theta - thetaC;
					delta -=  -scene->bn * sin(thetaD) -scene->bp * cos(thetaD) + baselineSquared * 0.5 /range;
				}       
				/*
				  Convert range diff to phase diff.
				*/
				scene->image[iIndex][jLoop] = (float)(delta * deltaToPhase); 
			} /* End else */
			range += scene->I.rangePixelSize * scene->dR;
		}
	}
	fprintf(stderr,"nIt Avg %f\n",(double)nIterTotal/(scene->rSize * scene->aSize));
	fprintf(stderr,"bn bp %f %f \n",scene->bn,scene->bp);

	return;
}




