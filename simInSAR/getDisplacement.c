#include "stdio.h"
#include"string.h"
#include <math.h>
#include "ers1/common/common.h"
#include "simInSARInclude.h"
/*
   Compute displacement for lat/lon from Displacement map.
*/
    double getDisplacement(double lat, double lon, 
                           displacementStructure *displacements)
{   
/*
    double delta;
    double latGrid,lonGrid;
    double xGrid,yGrid;
    double p1,p2,p3,p4;
    int xSize,ySize,latSize,lonSize;
    double x1,y1;
    double t,u;
    int i,j;

    if(displacements->coordType == LATLON) {
        error("getDisplacement: Lat/Lon type not supported yet\n");
        if(lon > 180.) lon = lon - 360.;
        latGrid = (lon - displacements->minC2)/displacements->deltaC2;
        lonGrid = (lat - displacements->minC1)/displacements->deltaC1;
        j = (int)lonGrid;
        i = (int)latGrid;
        t = lonGrid - (double)j;
        u = latGrid - (double)i;
    } else if(displacements->coordType == XY) {
        lltoxy(lat,lon,&x1,&y1,ROT);
        xGrid = (x1 - displacements->minC1)/displacements->deltaC1;
        yGrid = (y1 - displacements->minC2)/displacements->deltaC2;
        j = (int)xGrid;
        i = (int)yGrid;
        t = xGrid - (double)j;
        u = yGrid - (double)i; 
        if(i < 0 || i+1 > displacements->size2 ||
        j < 0 || j+1 > displacements->size1 ) return 0.0; 
    } else error("getDisplacement: Invalid coord type\n");     
*
   If off location of displacement map, return 0.
*
    p1 = (double)displacements->dR[i][j];
    p2 = (double)displacements->dR[i][j+1];
    p3 = (double)displacements->dR[i+1][j+1];
    p4 = (double)displacements->dR[i+1][j];
*
    Use bilinear interpolation to compute height.
*
    delta = ((1.0 - t)*(1.0 - u) * p1 + t * (1.0 - u) * p2 +
              t * u * p3 +                (1.0 - t) * u* p4);
*/
    return 0;
  
}
