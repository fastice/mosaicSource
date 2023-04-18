#include "stdio.h"
#include "string.h"
#include <math.h>
#include "mosaicSource/common/common.h"

void smlocateZD(double rxs, double rys, double rzs, double rvsx, double rvsy, double rvsz, double rsl, double *lat, double *lon, double lookDir, double trgalt)
{
    double px, py, pz;         /* pointing vector */
    double plx, ply, plz;      /* pt vector along velocity vector */
    double phx, phy, phz;      /* vector orthogonal to velocity vector */
    double pvx, pvy, pvz;      /* pt vector perpendicular to velocity vector */
    double uvx, uvy, uvz;      /* vertical unit vector */
    double xt, yt, zt;         /* target position vector */
    double elook;              /* effective look angle */
    double delook;             /* effective look angle error */
    double dtx, dty, dtz;      /* target position error */
    double wavl;               /* Wavelength in km */
    double reff;               /* effective earth radius */
    double rt;                 /* earth radius at target */
    double eccsq;              /* eccentricity square */
    double delthrsld = 1.0e-8; /*e ffective look angle convergence threshold */
    double temp, temp1, sinel, cosel;
    double sch; /* Spacecraft altitude */
    /*    double trgalt; */
    int32_t converge;
    int32_t nit, maxIt = 12; /* Number of iteration, and max it number */
    double fd = 0;
    double Rn, latTmp;

    wavl = 0;
    /*
        Calculate the unit pointing vector which is along the velocity vector
           (plx,ply,plz)=(rvsx,rvsy,rvsz)*wvl*fd/2*(|rvsx,rvsy,rvsz|)**2.0d0

        temp= wavl * fd /((rvsx * rvsx + rvsy * rvsy + rvsz*rvsz)*2.0);
        plx = rvsx * temp;
        ply = rvsy * temp;
        plz = rvsz * temp;
    */
    /*
        Calculate the vertical unit vector
        (uvx,uvy,uvz)=(rxs,rys,rzs)/|rxs,rys,rzs|
    */
    temp = sqrt(rxs * rxs + rys * rys + rzs * rzs);
    uvx = rxs / temp;
    uvy = rys / temp;
    uvz = rzs / temp;
    /*
          Calculate the horizontal component of the unit pointing vector
           which is orthogonal to the velocity vector

           (phx,phy,phz)=( (rvsx,rvsy,rvsz) x (uvx,uvy,uvz)/
               (|(rvsx,rvsy,rvsz) x (uvx,uvy,uvz)|) ) *
                sqrt(1-|plx,ply,plz|**2.0d0)
    */
    cross(rvsx, rvsy, rvsz, uvx, uvy, uvz, &phx, &phy, &phz);

    /*    temp1=( plx*plx +ply*ply + plz*plz );*/
    temp = sqrt(phx * phx + phy * phy + phz * phz);
    /*    temp=sqrt(1.- temp1)/temp;*/
    temp = 1.0 / temp;
    phx = phx * temp;
    phy = phy * temp;
    phz = phz * temp;
    /*
        Calculate the vertical component of the unit pointing vector
        which is perpendicular to the velocity vector
         (pvx,pvy,pvz)=(rvsx,rvsy,rvsz)/(|rvsx,rvsy,rvsz|) x (phx,phy,phz)
    */
    cross(rvsx, rvsy, rvsz, phx, phy, phz, &pvx, &pvy, &pvz);
    temp = sqrt(rvsx * rvsx + rvsy * rvsy + rvsz * rvsz);
    pvx = pvx / temp;
    pvy = pvy / temp;
    pvz = pvz / temp;
    /*
        Calculate the effective earth radius

        reff=sqrt((rxs**2+rys**2+rzs**2)/
              ((rxs**2+rys**2)/re**2 +rzs**2/rp**2) )
    */
    reff = sqrt((rxs * rxs + rys * rys + rzs * rzs) / ((rxs * rxs + rys * rys) / (EMAJOR * EMAJOR) + rzs * rzs / (EMINOR * EMINOR)));
    sch = sqrt(rxs * rxs + rys * rys + rzs * rzs) - reff;
    /*
           Find the initial effective look angle
           elook=acos((reff+sch)**2.0d0+rsl**2.0d0-reff**2.0d0)/(2*(reff+sch)*rsl)
    */
    temp = (pow(reff + sch, 2.0) + rsl * rsl - reff * reff) / (2.0 * (reff + sch) * rsl);
    elook = lookDir * acos(temp);

    nit = 0;
    converge = FALSE;
    while (nit < maxIt && converge == FALSE)
    {
        nit = nit + 1;
        /*
               Calculate the initial pointing vector
               (px,py,pz)=(plx,ply,plz)+cos(elook)*(pvx,pvy,pvz)+ sin(elook)*(phx,phy,phz)
        */
        cosel = cos(elook);
        sinel = sin(elook);
        /*
            px= plx + cosel*pvx + sinel*phx;
            py= ply + cosel*pvy + sinel*phy;
            pz= plz + cosel*pvz + sinel*phz;
    */
        px = cosel * pvx + sinel * phx;
        py = cosel * pvy + sinel * phy;
        pz = cosel * pvz + sinel * phz;
        /*
               Calculate the target location
               (tx,ty,tz)=(rxs,rys,rzs)+rsl*(px,py,pz)
        */
        xt = rxs + rsl * px;
        yt = rys + rsl * py;
        zt = rzs + rsl * pz;
        /*
               Calculate the differential position vector
               (dtx,dty,dtz)=rsl*(-sin(elook)*(pvx,pvy,pvz)+rsl*cos(elook)*(phx,phy,phz)
        */
        temp = -sinel * rsl;
        temp1 = cosel * rsl;
        dtx = temp * pvx + temp1 * phx;
        dty = temp * pvy + temp1 * phy;
        dtz = temp * pvz + temp1 * phz;
        /*
            Calculate the error of the effective look angle
              delta elook=(1-(tx**2+ty**2)/(re+trgalt)**2)+tz**2/(rp+trgalt)**2) /
                          (2*((tx*dtx+ty*dty)/(re+h)**2+tz*dtz/(rp+h)**2))
        */
        delook = (1.0 - (xt * xt + yt * yt) / pow(EMAJOR + trgalt, 2.0) - zt * zt / pow(EMINOR + trgalt, 2.0)) /
                 (2.0 * ((xt * dtx + yt * dty) / pow(EMAJOR + trgalt, 2.0) + zt * dtz / pow(EMINOR + trgalt, 2.0)));
        /*
            Check for convergence
        */
        if (fabs(delook) <= delthrsld)
        {
            converge = TRUE;
            xt = xt + delook * dtx;
            yt = yt + delook * dty;
            zt = zt + delook * dtz;

            rt = sqrt(xt * xt + yt * yt + zt * zt);

            *lat = asin(zt / rt); /* geocentric lat */
            *lat = tan(*lat);
            eccsq = 2.0 * F - F * F;
            latTmp = atan(*lat / (1. - eccsq)); /* geodetic lat for 0 elev */
            Rn = earthRadiusCurvatureWGS84(*lat);
            *lat = atan(*lat / (1. - eccsq * (Rn / (Rn + trgalt)))) * RTOD; /* geodetic lat */
            *lon = atan(yt / xt) * RTOD;
            if (xt < 0. && yt >= 0.)
                *lon = 180.0 + *lon;
            if (xt < 0. && yt < 0.)
                *lon = -180.0 + *lon;
            if (*lon < -180.0)
                *lon = *lon + 360.0;
            if (*lon > 180.0)
                *lon = *lon - 360.0;
        }
        else
        { /* If fabs... */
            elook += delook;
        }
    } /* End while nit... */
}
