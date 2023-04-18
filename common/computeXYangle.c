#include "stdio.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include "common.h"
/*
    Compute angle relative to north or south - pulls stdlat from dem
*/
void computeXYangle(double lat, double lon, double *xyAngle, xyDEM xydem)
{
    extern int32_t HemiSphere;
    double x, y;

    lltoxy1(lat, lon, &x, &y, xydem.rot, xydem.stdLat);
    *xyAngle = atan2(-y, -x);
    if (HemiSphere == SOUTH)
        *xyAngle += PI;
    return;
}

/*
    Compute angle relative to north or south - need to supply stdLat
*/
void computeXYangleNoDem(double lat, double lon, double *xyAngle, double stdLat)
{
    extern int32_t HemiSphere;
    extern double Rotation;
    double x, y;

    lltoxy1(lat, lon, &x, &y, Rotation, stdLat);
    *xyAngle = atan2(-y, -x);
    if (HemiSphere == SOUTH)
        *xyAngle += PI;
    return;
}
