#include <math.h> 
#include "azparams.h" 
#include "cRecipes/nrutil.h"
#inlcude <stdlib.h>
/*
   Estimate baseline parameters.
*/


     void constOnlyCoeffs(void *x,int i,double *afunc, int ma);
     void azCoeffs(void *x,int i,double *afunc, int ma);
     void azCoeffsLinear(void *x,int i,double *afunc, int ma);
     void constOnlyCoeffsLinear(void *x,int i,double *afunc, int ma);
    static void getBaselineRates(double *dbcds, double *dbhds, 
                                char *baseFile,double prf,double slPixSize);



    void computeAzParams( tiePointsStructure *tiePoints,
                          inputImageStructure inputImage,char *baseFile )
{   
    double Re, H, RNear, dr,rOffset;
    double *a; /* Solution for params */
    int nParams;
    double theta,  z, r0;
    double **v,**u, chisq;
    double *y, *sig, *w, *chsq;
    double slPixSize;
    int nSingleLook;
    int nData, ma;
    aZmodelValues *x;
    double dbhds,dbcds;
    int i,i1,k,j,npts;
    conversionDataStructure *cP;

    if(tiePoints->constOnlyFlag== TRUE) {
       if(tiePoints->linFlag==FALSE) nParams=1;
       else nParams=2;
    } else {
      if(tiePoints->linFlag==FALSE) nParams=2;
      else nParams=3;
    }
    fprintf(stderr,"nParams = %i\n",nParams);
    a = dvector(1,nParams);
    initllToImage( &inputImage);

    cP = &(inputImage.cpAll);

    Re = cP->Re;
    H = cP->H;
    RNear = cP->RNear;
     fprintf(stderr,"---------------------RNear %f %f \n",RNear,H);
    nSingleLook = inputImage.azimuthSize * inputImage.nAzimuthLooks;
    slPixSize=inputImage.azimuthPixelSize/inputImage.nAzimuthLooks;
    getBaselineRates(&dbcds,&dbhds,baseFile,inputImage.par.prf,slPixSize); 
/*
    if(inputImage.lookDir==LEFT) {
     dbcds *=-1.0; dbhds *=-1.0;
    }
*/
/*
    Range comp params
*/
    dr =  RangePixelSize * inputImage.nRangeLooks;
    rOffset = 0.0;
/*
   Init arrays for least squares fit.
*/
    npts = 0;
    nData = tiePoints->npts;
    for(i=0; i < nData; i++) if(fabs(tiePoints->phase[i]) < 1.0E6) npts++;
    ma = nParams;
    x = (aZmodelValues *) malloc( (npts + 1) * sizeof(aZmodelValues));
    y = dvector(1,npts);
    sig = dvector(1,npts);
    u = dmatrix(1, npts,1,ma);
    v = dmatrix(1,ma,1,ma);
    w = dvector(1,ma);    
/*
   Prep data for solution
*/
    j = 0;
 
    for(i=0; i < nData; i++) {
	if(fabs(tiePoints->phase[i]) < 1.0E6) { /* Use only good points */
	    i1 = j+1;
	    z = tiePoints->z[i];
	    r0 =  RNear + rOffset + tiePoints->r[i] * dr;
	    theta = acos( 
	    ( pow(r0,2.0) + 2.0 * Re * (H - z) + pow(H,2.0) - pow(z,2.0) )/
	    ( 2.0 * (Re + H) * r0 ) );

	    x[i1].r = r0;
	    x[i1].theta = theta;

	    x[i1].x = tiePoints->x[i]/(inputImage.azimuthSize * 
                       inputImage.nAzimuthLooks * AzimuthPixelSize);

            y[i1] = tiePoints->phase[i] + r0*cos(theta) * dbhds;

            if(tiePoints->constOnlyFlag==TRUE) y[i1] -= r0*sin(theta)*dbcds;    
/*
fprintf(stderr,"%i %f %f %f  % f %f\n",i1,tiePoints->phase[i],r0*sin(theta)*dbcds,r0*cos(theta) * dbhds,y[i1],z);
*/

	    sig[i1] = 10.0; 
	    j++;
	} /* End if */
    } /* End for i */

/*
   Solve for parameters
*/
    fprintf(stdout,";\n; Ntiepoints/Ngiven used= %i/%i\n;\n",npts,nData);
    if(tiePoints->constOnlyFlag == TRUE) {
          if(tiePoints->linFlag==FALSE) {
             svdfit((void *)x,y,sig,npts,a,ma,u,v,w, &chisq,&constOnlyCoeffs);
              fprintf(stdout,"%f %e %e %f \n",a[1],dbcds,dbhds,0.);
              fprintf(stderr,"%f %e %e %f\n",a[1],dbcds,dbhds,0.);
          } else {
              svdfit((void *)x,y,sig,npts,a,ma,u,v,w, &chisq,
                  &constOnlyCoeffsLinear);
              fprintf(stdout,"%f %e %e %f \n",a[1],dbcds,dbhds,a[2]);
              fprintf(stderr,"%f %e %e %f\n",a[1],dbcds,dbhds,a[2]);
          }
    } else {
       if(tiePoints->linFlag==FALSE) {
          svdfit((void *)x,y,sig,npts,a,ma,u,v,w, &chisq,&azCoeffs);
          fprintf(stderr,"--  %f %e %e %e  %f\n",
             a[1]*800000.,a[2],dbcds,dbhds,0.);
          fprintf(stdout,"%f %e %e %f\n",a[1]*800000.,a[2],dbhds,0.);
       } else {
          svdfit((void *)x,y,sig,npts,a,ma,u,v,w, &chisq,&azCoeffsLinear);
          fprintf(stderr,"--  %f %e %e %e  %f\n",
             a[1]*800000.,a[2],dbcds,dbhds,a[3]);
          fprintf(stdout,"%f %e %e %f\n",a[1]*800000.,a[2],dbhds,a[3]*800000.);

       }
    }


    return;
}


    void constOnlyCoeffs(void *x,int i,double *afunc, int ma) 
{      
     afunc[1] = 1;
     return;
}

    void constOnlyCoeffsLinear(void *x,int i,double *afunc, int ma) 
{  
     aZmodelValues xx,*xy; 

     xy = (aZmodelValues *)x;
     xx = xy[i];
   
     afunc[1] = 1;
     afunc[2] = xx.x;
     return;
}



    void azCoeffs(void *x,int i,double *afunc, int ma) 
{
     double theta,r;
     aZmodelValues xx,*xy;

     xy = (aZmodelValues *)x;
     xx = xy[i];
  
     theta=xx.theta;
     r=xx.r;

     afunc[1] = 800000.;
     afunc[2] = r*sin(theta);
     return;
}

   void azCoeffsLinear(void *x,int i,double *afunc, int ma) 
{
     double theta,r;
     aZmodelValues xx,*xy;

     xy = (aZmodelValues *)x;
     xx = xy[i];

     theta=xx.theta;
     r=xx.r;

     afunc[1] = 800000.;
     afunc[2] = r*sin(theta);
     afunc[3] = xx.x*800000.;
     return;
}



    static void getBaselineRates(double *dbcds, double *dbhds, 
                                char *baseFile,double prf,double slPixSize)
{
    char line[1024];
    int lineCount, eod;
    double x1,x2,x3;
    FILE *fp;

    fprintf(stderr,"prf,nSingleLook %f %f\n",prf,slPixSize);

    fp = openInputFile(baseFile);
    lineCount=getDataString(fp,lineCount,line,&eod);
    lineCount=getDataString(fp,lineCount,line,&eod);
    if( sscanf(line,"%lf%lf%lf",&x1,&x2,&x3) != 3)
        error("%s  %i of %s",
            "readOffsets -- Missing image parameters at line:",
             lineCount,baseFile); 
    *dbcds=x2/(prf*slPixSize);
    *dbhds=x3/(prf*slPixSize);
    fprintf(stderr,"dbcds,dbhds, %f %f\n", *dbcds,*dbhds);


    close(fp);


}
