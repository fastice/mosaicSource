#include "stdio.h"
#include"string.h"
#include <math.h>
#include "source/common/common.h"
/*
  Rotated version.
  Specify rotation in degrees with 0 deg. longitude on negative Y-axis 
  (pointing down, bottom of plot).  Positive rotation rotates 0 deg 
  longitude in CCW direction.
  NSIDC SSMI grids rotated 45. in Arctic and 0. in Antarctic
  Roger C's basin rotated 60. deg.

  Converted from fortran EFD   6.8.90 version to C 03/17/94

*/


    void lltoxy1(double alat,double alon,double *x, double *y, double dlam, 
                double slat)
{ 
/*  
  Converts from geodetic latitude and longitude to Polar
  Stereographic (x,y) coordinates for the Polar regions.
  The standard parallel (lat with no distortion) is 70 deg.
  Equations are from Snyder, J.P. 1982, Map Projections Used
  by the U.S. Geological Survey, Geological Survey Bulletin
  1532, U.S. Gov. Printing Office.  See JPL Tech. Memo.
  3349-85-101 for details.  Sub written by C.S. Morris-
  Apr 29, 1985 (Goddard?).  
  converted to MASSCOMP 7/7/87 DRT
  on SUN      11/89  EAF
  to C        03/94 IRJ

    alat in	latitude, degrees, +N, -S
    alon in     longitude, degrees 0-360(East-West)
    x	 out	x-coordinate in km
    y	 out	y-coordinate in km
    dlam in	rotation, degrees, 0-360
*/
    int i1;
    double t[2];
    double e,e2,re, sn,cm,rho,rlat,tmp;
/*
    Radius of earth (km) -Hughes Ellipsoid
    re=6378.273
*/
    re=6378.137;
/*
    Eccentricity of earth -Hughes Ellipsoid
    e2=0.006693883
    changed to WGS 84 10/14/05
*/
    e2=0.0066943801;
    e=sqrt(e2);
/*
    Standard parallel
    slat=70;
*/

/*
   Test for N or S hemi, set constants as necessary
   For SSM/I grid, Northern Hemisphere  sn=1.
*/
    if(alat < 0) sn = -1.; else sn = 1.;
/*
    Compute x, y
*/
    alat = sn * alat;
    alon = sn * alon;
    if(alat >= 89.995) {
        *x=0.; *y=0.;
    } else {
        rlat=alat;
        for(i1=0; i1 < 2; i1++) {
	    if(i1 == 1) rlat=slat;
	    t[i1] = tan( (PI/4.) - (rlat/(2.*RTOD)) );
            tmp = pow( (1.-e*sin(DTOR*rlat))/(1.+e*sin(DTOR*rlat)),(e/2.0) );
            t[i1] = t[i1]/tmp;
        }
	cm=cos(DTOR*slat)/sqrt(1.-e2 * pow(sin(DTOR*slat),2.0) );
        rho= re * cm * t[0] / t[1];
        *x =  rho * sn * sin( DTOR*(alon+dlam) );
        *y = -rho * sn * cos( DTOR*(alon+dlam) );
    }
/*
    alat = sn * alat;
    alon = sn * alon;
*/
    return;
}
