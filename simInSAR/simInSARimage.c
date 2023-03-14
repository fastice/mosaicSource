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

/*
  Compute flow direction and angle relative to north or south. 
*/
static  void computeXYslope(double lat,double lon, double *dzdx, double *dzdy,  xyDEM xydem)
{   
	extern int HemiSphere;
	double x,y,z,z1x,z1y;
	double dx,dy;

	dx = xydem.deltaX;
	dy = xydem.deltaY;
	lltoxy1(lat,lon,&x,&y,xydem.rot,xydem.stdLat);
	z = interpXYDEM(x,y,xydem);
	z1x = interpXYDEM(x+dx,y,xydem);
	z1y = interpXYDEM(x,y+dy,xydem);	
	if( z <= MINELEVATION || z1x <= MINELEVATION || z1y <= MINELEVATION ||
	    z > 10000.0 || z1x > 10000.0 || z1y > 10000.0 ) {
		*dzdx=0; *dzdy = 0.0;  return;
	}
	/*
	  Compute gradient
	*/
	*dzdx = (z1x-z)/(KMTOM*dx);
	*dzdy = (z1y-z)/(KMTOM*dy);
	return;
}


/************************STATIC ROUTINE DECLARATIONS**************************/

unsigned char getShelfMask(ShelfMask *imageMask,double x, double y);
/************************END STATIC ROUTINES DECLARATIONS*********************/

double interpXYVel(double x, double y, xyVEL *xyVel, double *vx, double *vy)
{
	double xi,yi;
	int i,j, iSize, jSize;
	double t,u, p1, p2, p3, p4;
	double z;
	
	xi = (x - xyVel->x0)/xyVel->deltaX;
	yi = (y - xyVel->y0)/xyVel->deltaY;
	*vx = bilinearInterp(xyVel->vx, xi, yi, xyVel->xSize, xyVel->ySize, -LARGEINT, -LARGEINT);
	*vy = bilinearInterp(xyVel->vy, xi, yi, xyVel->xSize, xyVel->ySize, -LARGEINT, -LARGEINT);
}




void simInSARimage(sceneStructure *scene, void *dem,xyVEL *xyVel)
{
	conversionDataStructure *cp;
	inputImageStructure *inputImage;
	xyDEM *xyDem;
	double R, theta,thetaD,delta,thetaC, rhoSp,rhoWGS,RCenter, RNear, RWGS,R0,ReWGS, ReH, Re;
	double lat,lon;
	double baselineSquared;
	double toRho, r, rg,drgStep, range;
	double deltaToPhase;
	double iFloat,rgSave;
	double hSp,hWGS;
	double x1,y1;
	double rOff,aOff;
	double vx, vy, vr, psi, dzdx, dzdy, dzdr;
	double hAngle, xyAngle;
	long long nIterTotal, nIter;	
	int recycle;
	int iLoop, jLoop, i,iIndex;
	fprintf(stderr,"SIM INSAR IMAGE\n");
	/*	  Dem type	*/
	xyDem = (xyDEM *)dem;		
	/*
	  Init coord conversions
	*/
	inputImage=&(scene->I);	
	initllToImageNew(inputImage);
	deltaToPhase = 4.0 * PI /scene->I.par.lambda;
	Re = inputImage->cpAll.Re;
	RCenter = inputImage->cpAll.RCenter;
	ReH = inputImage->cpAll.Re + inputImage->par.H ;
	thetaC = thetaRReZReH( RCenter, (Re+0), ReH);
	/* Constants for phase computation */
	fprintf(stderr,"Re,thetaC,RNear,RFar, %f %f   %f %f %f %i \n ",
	 Re,thetaC,inputImage->cpAll.RNear,inputImage->cpAll.RFar,inputImage->rangePixelSize,inputImage->nRangeLooks);
	/* Corrections from ml pixels to sl pixels */
	rOff = 0.5*(inputImage->nRangeLooks-1)*inputImage->par.slpR;
	/* First part in slp so need to divide nlooks to get in mlp pixels */
	aOff= 0.5*(inputImage->nAzimuthLooks-1)/inputImage->nAzimuthLooks;
	RNear = inputImage->cpAll.RNear;
	rgSave = -1.0;
	nIterTotal = 0;
	cp = &(inputImage->cpAll);
	/*
	  Loop over output grid.
	*/
	for (iLoop=0; iLoop < scene->aSize; iLoop++) {
		iFloat=iLoop*scene->dA + scene->aO - aOff - cp->azOff;
		i=(int)(iFloat + 0.5);
	    iIndex = iLoop;
		ReH=getReH(cp,&(scene->I), iFloat);
		rhoSp=rhoRReZReH( RNear,(Re+0), ReH);
		rg = Re * rhoSp -100;    /* Initial ground range on spherical earth */
		range = inputImage->cpAll.RNear + scene->rO*inputImage->rangePixelSize - rOff ;    /* Starting range - gets incremented in loop */
		scene->bn = scene->bnStart + (double)iIndex * scene->bnStep;
		scene->bp = scene->bpStart + (double)iIndex * scene->bpStep;
		if(i == 0) fprintf(stderr,"%f %f \n",scene->bn,scene->bp);
		if(i % 100 == 1) fprintf(stderr,"%i %f %lf \n",i,rg,(double)(nIterTotal/(iLoop * scene->rSize) ));
		/* Loop over range - not this program uses a ground range, rg, on a spherical earth to move across the swath. 
		   The actual ground ranges aren't correct, but the final calculation is  */
		for(jLoop=0; jLoop < scene->rSize; jLoop++) {
			/* Initialize iteration */
			delta = 0.0;
			nIter=0;
			recycle=FALSE;
			if( i >= 0 && i < cp->azSize) {
			        /*rangeAzimuthToLL(double range,double iFloat,double *lat, double *lon,sceneStructure *scene) */
				/*nIter=rangeAzimuthToLL(&rg, range, iFloat, rhoSp, ReH,  Re, &lat,&lon,&hWGS,inputImage,xyDem,0.1,0.5 );*/
				while(1) {
					/* For rg on spherical earth, compute lat/lon for the corresponding ellipsoidal earth 
					   and return the corresponding slant range */
					R0=groundRangeToLLNew(rg,iFloat,&lat,&lon,&(scene->I),recycle);
					recycle = TRUE; /* flags to reuse previously computed values that are still valid */
					/* Compute WGS radious at lat/lon */
					ReWGS = earthRadius(lat*DTOR,EMINOR,EMAJOR)*KMTOM;
					/* Comput rho for the ellipsodal radius */
					rhoWGS = rhoRReZReH( R0,(ReWGS+0), ReH);
					/*  Get height for lat/lon referenced to the ellipsoid */
					hWGS = getXYHeight(lat,lon,xyDem,Re, ELLIPSOIDAL);					
					/*  Compute height corrected R    	*/
					R = slantRange( rhoWGS, ReWGS+hWGS,  ReH);
					/* Not the most efficient iteration, but good enough */
					drgStep = (range-R)*0.5;
					if(R > (range-0.1))  break;				
					/* Speed convergence if way out */
					rg += drgStep;    /* Increment ground range */
					/* Avoid infinite loop */
					if(nIter > 8000){fprintf(stderr,"%f %f %f\n",Re,rhoSp,ReH); break; }
					nIter++;
				}  
				if(jLoop == 0) rgSave = rg;
				nIterTotal += nIter;
				withHeight(&(scene->I), &lat, &lon, iFloat,range, hWGS);

				if(scene->maskFlag == TRUE && scene->toLLFlag == FALSE ) {
					lltoxy1(lat, lon, &x1, &y1, xyDem->rot, xyDem->stdLat);
					scene->image[iIndex][jLoop] = getShelfMask(scene->imageMask, x1, y1);
				} else if(scene->heightFlag == TRUE && scene->toLLFlag == FALSE ) {
					scene->image[iIndex][jLoop]= (float)(hWGS);
				} else if(scene->toLLFlag == TRUE || scene->saveLLFlag == TRUE) {
					scene->latImage[iIndex][jLoop]=lat;
					scene->lonImage[iIndex][jLoop]=lon;
					/*fprintf(stderr,"%f %f\n",lat,lon,hWGS); error("stop");*/
					if( scene->maskFlag == TRUE ) {
						lltoxy1(lat,lon,&x1,&y1,xyDem->rot,xyDem->stdLat);
						scene->image[iIndex][jLoop] = getShelfMask(scene->imageMask,x1,y1);
					}
				} else {
					/* Use h on single spherical reference for phase
					  even though locally spherical value was used for location */
					hSp = getXYHeight(lat,lon,xyDem,Re, SPHERICAL);
					/* NOT VALIDATED */
					hSp = hSp + Re - inputImage->cpAll.Re;
					Re = inputImage->cpAll.Re;
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
					theta = thetaRReZReH( range, (Re+hSp), ReH);
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
						theta = thetaRReZReH( range, (Re+0.), ReH);
						thetaD = theta - thetaC;
						delta -= -scene->bn * sin(thetaD) -scene->bp * cos(thetaD) + baselineSquared * 0.5 /range;
					}   
					/*
						Get displacement
					*/  
					if(xyVel->xSize > 0) {
						lltoxy1(lat, lon, &x1, &y1, xyDem->rot,xyDem->stdLat);
						interpXYVel(x1, y1, xyVel, &vx, &vy);
						/* Compute rotation angles */	
						computeXYangle(lat,lon, &xyAngle, *xyDem);
						hAngle = computeHeading(lat, lon, 0., inputImage,cp);
				fprintf(stderr, "%f\n", hAngle * 57.29);
						/* Rotate to ground range */
						vr = vx * cos(hAngle - xyAngle) - vy * sin(hAngle - xyAngle);
						/* Get slopes */
						computeXYslope(lat, lon, &dzdx, &dzdy,  *xyDem);
						/* Rotate to ground range direction */
						dzdr = dzdx * cos(hAngle - xyAngle) - dzdy * sin(hAngle - xyAngle);
						
						psi = psiRReZReH(range ,(Re + hSp), ReH);
						vr = vr * sin(psi) - vr * dzdr * cos(psi);
					}
					/*
					  Convert range diff to phase diff.
					*/
					scene->image[iIndex][jLoop] = (float)(delta * deltaToPhase); 
	scene->image[iIndex][jLoop] = vr * 6./365. * deltaToPhase;	
	
				
				} /* End else */
				range += inputImage->rangePixelSize * scene->dR;
			} else {
				/* No data cases because outsize az range */
				if(jLoop == 10 ) fprintf(stderr,"%i %f %i",iLoop,iFloat,i);
				if(scene->maskFlag == TRUE && scene->toLLFlag == FALSE ) {
					scene->image[iIndex][jLoop]=0;
				} else if(scene->heightFlag == TRUE && scene->toLLFlag == FALSE ) {
					scene->image[iIndex][jLoop]=(float)-9999.0;
				} else if(scene->toLLFlag == TRUE) {
					scene->latImage[iIndex][jLoop]=-LARGEINT;
					scene->lonImage[iIndex][jLoop]=-LARGEINT;
					if( scene->maskFlag == TRUE ) { scene->image[iIndex][jLoop]=0.0;  }
				} else {
					scene->image[iIndex][jLoop]=-2.0e9;
				}
			}				
		} /* end for jLoop ... */
	} /* end for iLoop .... */
	fprintf(stderr,"nIt Avg %f\n",(double)nIterTotal/(scene->rSize * scene->aSize));
	fprintf(stderr,"bn bp %f %f \n",scene->bn,scene->bp);
	return;
}




