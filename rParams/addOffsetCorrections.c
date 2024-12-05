#include <math.h>
#include "azparams.h"
#include <stdio.h>
/*
    Compute motion correction.
*/
void addOffsetCorrections(inputImageStructure inputImage,
                          tiePointsStructure *tiePoints)
{
    extern int32_t HemiSphere;
    int32_t i, j;
    conversionDataStructure *cP;
    double x, y, z;
    double lat, lon;
    double azOrigin;
    double num, den, rho1, rho2;
    double range1, range2, azimuth1, azimuth2;
    double va;
    double da, dgr;
    double rot;
    double H, Re, RNear;
    double r1, r2;
    double gRange1, gRange2;
    double hAngle, xyAngle, rotAngle;
    double x1, y1;
    double tSign;
    double deltaT;
    double dlat;
    double hAngle1;
    double phiDisplacement;
    float **fimage;
    /*
       Setup stuff for Shusun Li's conversion routines
    */
    initllToImage(&inputImage);
    cP = &(inputImage.cpAll);
    for (i = 0; i < 6; i++)
    {
        lat = inputImage.latControlPoints[i];
        lon = inputImage.lonControlPoints[i];
    }
    H = cP->H;
    Re = cP->Re;
    RNear = cP->RNear;
    fprintf(stderr, "---------------------RNear %f %f\n", RNear, H);
    deltaT = tiePoints->nDays / 365.25;
    /*
        Loop over tie points
    */
    for (i = 0; i < tiePoints->npts; i++)
    {
        lat = tiePoints->lat[i];
        lon = tiePoints->lon[i];
        /*
            Convert lat/lon to image coordinates
        */
        z = 0.0;
        llToImage(lat, lon, z, &range1, &azimuth1, &inputImage);
        /*
           IRJ changed 10/23/00 (see notes 10/12/00)
           ****if(HemiSphere == SOUTH) dlat = -0.02; else dlat = 0.02;
        */
        if (inputImage->lookDir == LEFT)
            dlat = -0.02;
        else
            dlat = 0.02;
        llToImage(lat + dlat, lon, z, &range2, &azimuth2, &inputImage);
        /*
            Compute azimuth displacement for latitude displacement
        */
        da = (azimuth2 - azimuth1) * inputImage.azimuthPixelSize;
        if (inputImage.passType == DESCENDING)
            da *= -1.0;
        /*
            Compute ground range  displacement for latitude displacement
        */
        r1 = RNear + inputImage.rangePixelSize * range1;
        num = 2.0 * Re * (Re + H + z) + pow(H, 2.0) + pow(z, 2.0) - pow(r1, 2.0);
        den = 2.0 * Re * (Re + H + z) + 2.0 * z * H;
        rho1 = acos(num / den);
        gRange1 = Re * rho1;

        r2 = RNear +  inputImage.rangePixelSize * range2;
        num = 2.0 * Re * (Re + H + z) + pow(H, 2.0) + pow(z, 2.0) - pow(r2, 2.0);
        den = 2.0 * Re * (Re + H + z) + 2.0 * z * H;
        rho2 = acos(num / den);
        gRange2 = Re * rho2;

        dgr = gRange2 - gRange1;
        /*
          Compute track heading
        */
        if (inputImage.lookDir == RIGHT)
            hAngle = atan2(da, dgr);
        else if (inputImage.lookDir == LEFT)
            hAngle = atan2(da, -dgr);
        else
            error("*** computeHeading: invalid lookDir ***");
        /*
           Compute angle for ps coordinats from north
        */
        rot = 45.0;
        if (lat < 0)
            rot = 0.0; /* Use rot=0 for southern latitudes */
        lltoxy1(lat, lon, &x1, &y1, rot, tiePoints->stdLat);
        xyAngle = atan2(-y1, -x1);
        if (HemiSphere == SOUTH)
            xyAngle += PI;
        /*
           Compute angle for CW rotation of ps to radar
        */
        rotAngle = hAngle - xyAngle;
        /*
            Compute vyra by ccw rot from ps to ra
        */
        va = tiePoints->vx[i] * sin(rotAngle) + tiePoints->vy[i] * cos(rotAngle);
        /*
           Compute displacement - meters
        */
        tiePoints->phase[i] -= va * deltaT;
    }

    return;
}

/*
fprintf(stderr,"lat,lon,z -- %f %f %f\n", lat,lon,z);
fprintf(stderr,"dgr,da    -- %f %f\n", dgr,da);
fprintf(stderr,"range1,2  -- %f %f\n", range1,range2);
fprintf(stderr,"grange1,2 -- %f %f \n",gRange1,gRange2);
fprintf(stderr,"azimuth1,2 -- %f %f\n",azimuth1,azimuth2);
fprintf(stderr,"RE,H,rho1,rho2 %f %f %f %f\n",Re,H,rho1*57.29,rho2*57.29);
fprintf(stderr,"hAngle   -- %f\n",hAngle*RTOD);
fprintf(stderr,"xyAngle  -- %f\n",xyAngle*RTOD);
fprintf(stderr,"vx,vy -- %f %f\n",tiePoints->vx[i],tiePoints->vy[i]);
fprintf(stderr,"accross track velocity ---    %f\n",vyra);
fprintf(stderr,"psi --- %f  %f\n",psi*RTOD,theta);
fprintf(stderr,"nDays -- %f\n",tiePoints->nDays);
fprintf(stderr,"\nphiDisplacement --  %f\n",phiDisplacement);
*/
