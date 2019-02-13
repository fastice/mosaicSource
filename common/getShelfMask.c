#include "stdio.h"
#include"string.h"
#include "math.h"
#include "common.h"
/*
#include "source/GeoCodeDEM_p/geocodedem.h"
#include "source/computeVH_p/computevh.h"
#include "mosaic3d.h"
*/
/*
   read and return shelf mask value for a given area
*/
    unsigned char getShelfMask(ShelfMask *shelfMask,double x, double y)
{
    int xi,yi;

    xi = (x-shelfMask->x0)/shelfMask->deltaX + 0.5; 
    yi = (y-shelfMask->y0)/shelfMask->deltaY + 0.5;

    if(xi < 0 || yi < 0 || 
        xi >= shelfMask->xSize || yi >= shelfMask->ySize) return(GROUNDED);
    else return(shelfMask->mask[yi][xi]);

} 

/*
fprintf(stderr,"%i %i %f %f %f %lf %lf\n",xi,yi,shelfMask->x0,shelfMask->y0,shelfMask->deltaX,x,y);
*/
