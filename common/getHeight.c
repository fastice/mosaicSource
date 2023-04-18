#include "stdio.h"
#include "string.h"
#include <math.h>
#include "mosaicSource/common/common.h"
/*
   Compute height for lat/lon from DEM.
   Modified 2/7/96 to add mode for converting from ellipsoidal to spherical
   heights.

   lat, lon inputs should be in degrees.
   Re in KM
*/
double getHeight(double lat, double lon, demStructure *dem, double Re, int32_t heightFlag)
{
    double h;
    double latGrid, lonGrid;
    double p1, p2, p3, p4;
    double x1, y1, xGrid, yGrid;
    double rDist;
    double t, u;
    int32_t i, j, iSize, jSize;

    if (dem->coordType == LATLON)
    {
        if (lon > 180.)
            lon = lon - 360.;
        lonGrid = (lon - dem->minLon) / dem->deltaLon;
        latGrid = (lat - dem->minLat) / dem->deltaLat;
        j = (int)lonGrid;
        i = (int)latGrid;
        t = lonGrid - (double)j;
        u = latGrid - (double)i;
        jSize = dem->lonSize;
        iSize = dem->latSize;
    }
    else if (dem->coordType == XY)
    {
        lltoxy(lat, lon, &x1, &y1, ROT);
        yGrid = (y1 - dem->minLon) / dem->deltaLon;
        xGrid = (x1 - dem->minLat) / dem->deltaLat;
        j = (int)xGrid;
        i = (int)yGrid;
        t = xGrid - (double)j;
        u = yGrid - (double)i;
        iSize = dem->lonSize;
        jSize = dem->latSize;
    }
    else
        error("getHeight: Invalid coordinate type\n");
    /*
       If location of DEM return border height.
    */
    if (i < 0)
        return 0;
    else if (i + 2 > iSize)
        return 0;
    if (j < 0)
        return 0;
    else if (j + 2 > jSize)
        return 0;

    p1 = (double)dem->demData[i][j];
    p2 = (double)dem->demData[i][j + 1];
    p3 = (double)dem->demData[i + 1][j + 1];
    p4 = (double)dem->demData[i + 1][j];
    /*
        Use bilinear interpolation to compute height.
    */
    h = ((1.0 - t) * (1.0 - u) * p1 + t * (1.0 - u) * p2 +
         t * u * p3 + (1.0 - t) * u * p4);
    if (heightFlag == ELLIPSOIDAL)
        return h;
    /*
       Convert elliptical height to spherical
    */
    /*
       Radial distance from earth center to point
    */
    rDist = h + earthRadius(lat * DTOR, EMINOR, EMAJOR) * KMTOM;
    /*
       Height reference to spherical earth
    */
    h = rDist - Re;

    return h;
}
