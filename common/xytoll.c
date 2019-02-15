#include "stdio.h"
#include"string.h"
#include <math.h>
#include "mosaicSource/common/common.h"
/*
 
  Rotated version.
  Specify rotation in degrees with 0 deg. longitude on negative Y-axis 
  (pointing down, bottom of plot).  Positive rotation rotates 0 deg 
  longitude in CCW direction.
  NSIDC SSMI grids rotated 45. in Arctic and 0. in Antarctic
  Roger C's basin rotated 60. deg.

  Converted from fortran EFD   6.8.90 version to C 03/17/94
*/
    void xytoll(double x, double y,int hemi, double *alat,double *alon,
                double dlam)
{ 
/* 

  Converts from Polar Stereographic (x,y) coordinates for the 
  Polar regions to geodetic latitude and longitude.
  The standard parallel (lat with no distortion) is 70 deg.
  Equations are from Snyder, J.P. 1982, Map Projections Used
  by the U.S. Geological Survey, Geological Survey Bulletin
  1532, U.S. Gov. Printing Office.  See JPL Tech. Memo.
  3349-85-101 for details.  Sub written by C.S. Morris-
  Apr 29, 1985 (Goddard?).  
  converted for MASSCOMP 7/7/87 - DRT

  ARGUMENTS
	x       in	x-coordinate in km
	y	in	y-coordinate in km
	alat    out	latitude, degrees, +N, -S
	alon    out	longitude, degrees 0-360(East-West)
	hemi    in	Hemisphere(N or S)
	dlam    in	rotation, degrees, 0-360
*/
    double e,e2,re, slat,sn,rho,t,cm,xpr,ypr,chi,tmp;
/*
    Radius of earth (km) -Hughes Ellipsoid
	re=6378.273
    changed to WGS 84, 10/14/05
*/
      re=6378.137;
/*
    Eccentricity of earth -Hughes Ellipsoid
	e2=0.006693883
*/
	e2=0.0066943801;
	e=sqrt(e2);
/*
    Standard parallel
*/
    slat=70;
/*
    For SSM/I grid,
    Test for N or S hemi, set constants as necessary
*/
    if(hemi == SOUTH) sn=-1.; else sn=1.;
/*
  Compute lat and long
*/
	rho=sqrt( pow(x,2.0)+pow(y,2.0) );
	if(rho <= 0.1) {
          *alon=0.;
	  if(sn <= 0.) *alat=-90.; else *alat=90.;
	} else {
	  cm=cos(DTOR*slat) /sqrt( 1.-e2*( pow(sin(DTOR*slat),2.0) ) );

	  t=tan((PI/4.)-(slat/(2.*RTOD)));
          tmp =  (1.-e*sin(DTOR*slat)) / ( 1.+e*sin(DTOR*slat) );
          t=t/pow( tmp, (e/2.) );
	  t=rho*t/(re*cm);
	  chi=(PI/2.)-2.*atan(t);
	  *alat=chi+
                ((e2/2.)+(5.*pow(e2,2.0)/24.)+(pow(e2,3.0)/12.))*sin(2.0*chi);
          *alat = *alat +
             ((7.*pow(e2,2.0)/48.)+(29.*pow(e2,3.0)/240.))*sin(4.*chi);
          *alat = *alat + (7.*pow(e2,3.0)/120.)*sin(6.*chi);
	  *alat = sn * (*alat) * RTOD;
	  xpr = sn * x;
	  ypr = sn * y;
	  *alon = RTOD * atan2(xpr,-ypr)-sn*dlam;
	  *alon = sn * (*alon);
	  if(*alon < 0.) *alon=*alon+360.;
	  if(*alon > 360.) *alon=*alon-360.;
	}
    return;
}
