#include "stdio.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include "common.h"
/*
    Angle between the map grid and true north, differenced against the SAR heading to
    rotate between map (x,y) and radar (range, azimuth) frames.

    All three variants now go through grimpXYAngle(), which returns
    PI/2 + the meridian convergence.  For a polar stereographic that is the original
    atan2(-y,-x), plus PI in the southern hemisphere, evaluated bit for bit; for UTM it
    is the true convergence.  The two are the same formula, not two conventions:
    for a polar aspect the convergence is (lon - lon_0) in the north and -(lon - lon_0)
    in the south, and PI/2 plus that is identical to the atan2 expression.

    Note which projection each variant uses, because they are not interchangeable:
    computeXYangle() takes the DEM's, which is what xyGetZandSlope() wants since it is
    rotating dzdx/dzdy computed on the DEM grid.  Anything rotating a VELOCITY wants the
    output grid instead, and those call sites pass it explicitly.
*/

void computeXYangle(double lat, double lon, double *xyAngle, xyDEM xydem)
{
    /* Angle in the DEM's own projection -- for rotating DEM-grid derivatives. */
    double x, y;

    llToXYProj(lat, lon, &x, &y, &(xydem.proj));
    *xyAngle = grimpXYAngle(lat, lon, x, y, &(xydem.proj));
}

void computeXYangleProj(double lat, double lon, double *xyAngle, const grimpProj *proj)
{
    /* Angle in an explicitly supplied projection -- normally the output grid's. */
    double x, y;

    llToXYProj(lat, lon, &x, &y, proj);
    *xyAngle = grimpXYAngle(lat, lon, x, y, proj);
}

void computeXYangleNoDem(double lat, double lon, double *xyAngle, double stdLat)
{
    /* Legacy entry point: rotation from the global, standard parallel from the caller.
       Kept because the three baseline estimators call it this way. */
    extern int32_t HemiSphere;
    extern double Rotation;
    grimpProj proj = grimpProjFromLegacy(Rotation, stdLat, HemiSphere);
    computeXYangleProj(lat, lon, xyAngle, &proj);
}
