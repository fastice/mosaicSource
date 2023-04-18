#include "stdio.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include "common.h"
/*
    Rotate flow direction vector to xy (PS)  coordinates
*/
void rotateFlowDirectionToXY(double drn, double dan, double *dxn, double *dyn, double xyAngle, double hAngle)
{
    double rotAngle, cosRot, sinRot;
    /*
       Compute angle for CW rotation of ps to radar
    */
    rotAngle = hAngle - xyAngle;
    /*
        Compute range displacement vyra by ccw rot from ps to ra
    */
    sinRot = sin(rotAngle);
    cosRot = cos(rotAngle);
    *dxn = drn * cosRot + dan * sinRot;
    *dyn = -drn * sinRot + dan * cosRot;
}
