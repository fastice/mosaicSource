#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "tiePoints.h"
#include "math.h"

/*
    Add baseline corrections that were removed in the unwrapped image to 
    tiepoint phases.
*/
     void addBaselineCorrections(char *baselineFile, tiePointsStructure *tiePoints, inputImageStructure inputImage)
{
    FILE *fp;
    conversionDataStructure *cP;
    int quadB;
    int n;
    double Bn1, Bn2, Bp1,Bp2,dBn1,dBn2,dBp1,dBp2,BnvC,BpvC;
    double dBnQ1,dBnQ2,dBpQ1,dBpQ2,dBnQv,dBpQv;
    double Bnv, Bpv,dBnv,dBpv,bsq;
    double dum;
    double c1,c2,xAz;
    double theta,thetaC, thetaD,delta;
    double Rc,Re, H,dtop,r0,RNear,rOffset,dr;
    double azLength;
    double omegaA1, omegaA2,omegaA,azRamp;
    int lineCount=0, eod;
    int i;
    double lambda1;
    char line[256];
/*
   Get params
*/
     lambda1=inputImage.par.lambda;
     if(inputImage.isInit != TRUE)   initllToImageNew( &inputImage);
     cP = &(inputImage.cpAll);
     H = inputImage.par.H;
     tiePoints->H = H;
     Re = cP->Re;
     tiePoints->Re = Re;
     RNear = cP->RNear;
     tiePoints->RNear = RNear;
     Rc = (inputImage.par).rc;
     omegaA1 = -9999.0;
     omegaA2 = -9999.0;
     dBnQ1=0.0; dBnQ2=0.0; dBpQ1=0.0; dBpQ2=0.0;
/*
   Input baseline info
*/
    fp = openInputFile(baselineFile);
    lineCount=getDataString(fp,lineCount,line,&eod);
/*
   READ FIRST BASELINE
*/
    if(tiePoints->noRamp == TRUE) {
/* 
    Apply cos(thetad) dependence to ramp 
*/
        if( sscanf(line,"%lf%lf%lf%lf",&Bn1,&Bp1,&dBn1,&omegaA1) != 4)
            error("%s  %i",
            "addBaselineCorrections -- Missing baseline params at line:"
            ,lineCount);
        fprintf(stderr,"Using noRamp option. Cos(thetad) correction applied\n");

    } else if(tiePoints->dBpFlag == TRUE) {
/*
    Use dBn instead of omegaA for post apl stuff
*/
        n = sscanf(line,"%lf%lf%lf%lf%lf%lf",&dum,&dum,&dum,&dum,&dum,&dum);
        quadB=FALSE;
        if(n == 6) quadB=TRUE;

        if(quadB==FALSE) {
           if( sscanf(line,"%lf%lf%lf%lf",&Bn1,&Bp1,&dBn1,&dBp1) != 4)
               error("%s  %i",
              "addBaselineCorrections -- Missing baseline params at line:"
              ,lineCount);
           fprintf(stderr,"Using dBp option. Bp updated along track\n");
        } else {
           sscanf(line,"%lf%lf%lf%lf%lf%lf",
                &Bn1,&Bp1,&dBn1,&dBp1,&dBnQ1,&dBpQ1);
           fprintf(stderr,"Using QUADRATIC Baseline option.\n");
        }

    } else {
        if( sscanf(line,"%lf%lf%lf",&Bn1,&Bp1,&dBn1) != 3)
            error("%s  %i",
            "addBaselineCorrections -- Missing baseline params at line:"
            ,lineCount);
        fprintf(stderr,"Using simple linear azmuth phase ramp, OmegaA\n");
    }

    lineCount=getDataString(fp,lineCount,line,&eod);
/*
    READING SECOND BASELINE
*/
    if(tiePoints->noRamp == TRUE) {
/* 
    Apply cos(thetad) dependence to ramp 
*/
        if( sscanf(line,"%lf%lf%lf%lf",&Bn2,&Bp2,&dBn2,&omegaA2) != 4)
            error("%s  %i",
            "addBaselineCorrections -- Missing baseline params at line:"
            ,lineCount);
        fprintf(stderr,"Using noRamp option. Cos(thetad) correction applied\n");

    } else if(tiePoints->dBpFlag == TRUE) {
/*
    Use dBn instead of omegaA for post apl stuff
*/
        n = sscanf(line,"%lf%lf%lf%lf%lf%lf",&dum,&dum,&dum,&dum,&dum,&dum);
        quadB=FALSE;
        if(n == 6) quadB=TRUE;

        if(quadB==FALSE) {
           if( sscanf(line,"%lf%lf%lf%lf",&Bn2,&Bp2,&dBn2,&dBp2) != 4)
               error("%s  %i",
              "addBaselineCorrections -- Missing baseline params at line:"
              ,lineCount);
           fprintf(stderr,"Using dBp option. Bp updated along track\n");
        } else {
           sscanf(line,"%lf%lf%lf%lf%lf%lf",
                &Bn2,&Bp2,&dBn2,&dBp2,&dBnQ2,&dBpQ2);
           fprintf(stderr,"Using QUADRATIC Baseline option.\n");
        }

    } else {
        if( sscanf(line,"%lf%lf%lf",&Bn2,&Bp2,&dBn2) != 3)
            error("%s  %i",
            "addBaselineCorrections -- Missing baseline params at line:"
            ,lineCount);
        fprintf(stderr,"Using simple linear azmuth phase ramp, OmegaA\n");
    }
/*
    Ouput baseline parms
*/
    if(tiePoints->dBpFlag  == TRUE) fprintf(stderr,"; dBp flag set\n");
    if(tiePoints->noRamp  == TRUE) fprintf(stderr,"; noRamp flag set\n");
    fprintf(stdout,";\n; Number of lines of baseline data\n;\n 3\n");
    fprintf(stdout,";\n; First flattening baseline\n;\n");
    if(tiePoints->dBpFlag == TRUE) {
        fprintf(stdout,"%8.4f  %8.4f  %8.4f %10.7f\n", Bn1, Bp1, dBn1,dBp1);
        fprintf(stdout,";\n; Second flattening baseline\n;\n");
        fprintf(stdout,"%8.4f  %8.4f  %8.4f %10.7f\n", Bn2, Bp2, dBn2,dBp2);   
    } else {
        fprintf(stdout,"%8.4f  %8.4f  %8.4f %10.7f\n", Bn1, Bp1, dBn1,omegaA1);
        fprintf(stdout,";\n; Second flattening baseline\n;\n");
        fprintf(stdout,"%8.4f  %8.4f  %8.4f %10.7f\n", Bn2, Bp2, dBn2,omegaA2);
    }
/*
   Virtual baseline 
*/
    BnvC = Bn1 + Bn2;
    dBnv = dBn1 + dBn2;
    dBnQv = dBnQ1 + dBnQ2;    
    if(tiePoints->dBpFlag  == TRUE) {
         BpvC = Bp1 + Bp2;
         dBpv = dBp1 + dBp2;
         dBpQv = dBpQ1 + dBpQ2;
    }
    Bpv = Bp1 + Bp2;
/* Retain these value for cases where only a subset of the params are est */
    tiePoints->BnCorig=BnvC; tiePoints->BpCorig=BpvC;
    tiePoints->dBnorig=dBnv; tiePoints->dBporig=dBpv;
    tiePoints->dBnQorig=dBnQv; tiePoints->dBpQorig=dBpQv;

    omegaA = omegaA1 + omegaA2;
/*
   Constants for phase computation - note constant ReH used because this is compensating for reflattencomplex
*/
    dtop = (4.0 * PI / lambda1);
    c1 = 2.0 * H * Re + pow(H,2.0);
    c2 = 2.0 * (Re + H);
    thetaC = acos( (pow(Rc,2.0) + c1) / (c2*Rc) );
    tiePoints->thetaC = thetaC;
    azLength = AzimuthPixelSize * inputImage.nAzimuthLooks * inputImage.azimuthSize;
/*
    Range comp params
*/
    dr =  RangePixelSize * inputImage.nRangeLooks;
/*    if(inputImage.passType == DESCENDING) { */
          rOffset = 0.0;
/*        rOffset = (inputImage.rangeSize  - 1) * dr;
          dr = -dr
    } else {
        error("phaseToDelta ASCENDING PASS not implemented yet\n");
    }
*/
/*
    Compute and add corrections
*/

    for(i=0; i < tiePoints->npts; i++) {
/*
    Compute updated baseline and squared baseline
*/
          xAz=(tiePoints->x[i]/azLength);
          Bnv = BnvC + dBnv * xAz + dBnQv*xAz*xAz;
          if(tiePoints->dBpFlag  == TRUE) {
              Bpv = BpvC + dBpv * xAz + dBpQv*xAz*xAz;
          } else 
              azRamp = -(tiePoints->x[i]/azLength) * omegaA *
                       inputImage.nAzimuthLooks * inputImage.azimuthSize;
          bsq = Bnv*Bnv + Bpv*Bpv; 
          tiePoints->bsq[i] = bsq;         
/*
   Estimate range
*/
          r0 =  RNear + rOffset + tiePoints->r[i] * dr;
/*
    Compute theta and thetad for a curved earth
*/
          theta = acos( (pow(r0,2.0) + c1) / (c2*r0) );
          thetaD = (theta - thetaC);
/*
    Compute range and phase difference
*/
          delta = -Bnv * sin(thetaD) - Bpv * cos(thetaD) + bsq/(2.0*r0);
          tiePoints->phase[i] += dtop * delta;
/*
   Added 12/16/94, adds back what was removed in correlate
*/
          if(tiePoints->noRamp == TRUE &&
             (tiePoints->dBpFlag == FALSE)) tiePoints->phase[i] += azRamp;          
/*
 Added this line 12/09/94 for delta^2 correction
*/
          tiePoints->delta[i] = delta;
    }
fprintf(stderr,"leaving addBaselineCorrections.c\n");
 }
