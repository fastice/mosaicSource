#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "rparams.h"
#include "math.h"
/*
    Add baseline corrections that were removed in the unwrapped image to   tiepoint32_t phases.
*/

static void writeOriginalBaselines(tiePointsStructure *tiePoints,
                                   double Bn1, double Bp1, double dBn1, double dBp1,
                                   double Bn2, double Bp2, double dBn2, double dBp2)
{
 /*
    Ouput baseline parms
*/
    if (tiePoints->dBpFlag == TRUE)
        fprintf(stderr, "; dBp flag set\n");
    if (tiePoints->noRamp == TRUE)
        fprintf(stderr, "; noRamp flag set\n");
    fprintf(stdout, ";\n; Number of lines of baseline data\n;\n 3\n");
    fprintf(stdout, ";\n; First flattening baseline\n;\n");
    fprintf(stdout, "%8.4f  %8.4f  %8.4f %10.7f\n", Bn1, Bp1, dBn1, dBp1);
    fprintf(stdout, ";\n; Second flattening baseline\n;\n");
    fprintf(stdout, "%8.4f  %8.4f  %8.4f %10.7f\n", Bn2, Bp2, dBn2, dBp2);
}

void getBaselineFile(char *baselineFile, tiePointsStructure *tiePoints, inputImageStructure inputImage)
{
    FILE *fp;
    conversionDataStructure *cP;
    int32_t quadB;
    int32_t n;
    double Bn1, Bn2, Bp1, Bp2, dBn1, dBn2, dBp1, dBp2, BnvC, BpvC;
    double dBnQ1, dBnQ2, dBpQ1, dBpQ2, dBnQv, dBpQv;
    double Bnv, Bpv, dBnv, dBpv, bsq;
    double dum;
    double c1, c2, xAz;
    double theta, thetaC, thetaD, delta;
    double Rc, Re, H, r0, RNear, rOffset, dr;
    /*    double dtop; */
    double azLength;
    double omegaA1, omegaA2, omegaA, azRamp;
    int32_t lineCount = 0, eod;
    int32_t i;
    char line[256];
    /*
       Get params
    */
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
    dBnQ1 = 0.0;
    dBnQ2 = 0.0;
    dBpQ1 = 0.0;
    dBpQ2 = 0.0;
    /*
       Input baseline info
    */
    fp = fopen(baselineFile, "r");

    if (fp == NULL)
    {
        tiePoints->initWithSV = TRUE;
        writeOriginalBaselines(tiePoints, 0, 0, 0, 0, 0, 0, 0, 0);
        return;
    }
    tiePoints->initWithSV = FALSE;
    lineCount = getDataString(fp, lineCount, line, &eod);
    /*
       READ FIRST BASELINE
    */
    /*
        Use dBn instead of omegaA for post apl stuff
    */
    if (sscanf(line, "%lf%lf%lf%lf", &Bn1, &Bp1, &dBn1, &dBp1) != 4)
        error("%s  %i", "addBaselineCorrections -- Missing baseline params at line:", lineCount);

    lineCount = getDataString(fp, lineCount, line, &eod);
    /*
        READING SECOND BASELINE
    */
    if (sscanf(line, "%lf%lf%lf%lf", &Bn2, &Bp2, &dBn2, &dBp2) != 4)
        error("%s  %i", "addBaselineCorrections -- Missing baseline params at line:", lineCount);

    writeOriginalBaselines(tiePoints, Bn1, Bp1, dBn1, dBp1, Bn2, Bp2, dBn2, dBp2);

    /*
       Virtual baseline
    */
    BnvC = Bn1 + Bn2;
    dBnv = dBn1 + dBn2;
    dBnQv = dBnQ1 + dBnQ2;
    if (tiePoints->dBpFlag == TRUE)
    {
        BpvC = Bp1 + Bp2;
        dBpv = dBp1 + dBp2;
        dBpQv = dBpQ1 + dBpQ2;
    }
    Bpv = Bp1 + Bp2;
    /*
       Retain these value for cases where only a subset of the params are est
    */
    tiePoints->BnCorig = BnvC;
    tiePoints->BpCorig = BpvC;
    tiePoints->dBnorig = dBnv;
    tiePoints->dBporig = dBpv;
    tiePoints->dBnQorig = dBnQv;
    tiePoints->dBpQorig = dBpQv;
    /*
       Constants for phase computation
    */
    thetaC = thetaRReZReH(Rc, (Re + 0), (Re + H));
    tiePoints->thetaC = thetaC;
}
