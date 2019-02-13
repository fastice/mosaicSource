#include "stdio.h"
#include"string.h"
#include <math.h>
#include <stdlib.h>
#include "common.h"
/*
    Rotate flow direction vector to SAR coordinates
*/
    void rotateFlowDirectionToRA(double dxn,double dyn,  double *dan, double *drn,  double xyAngle, double hAngle)
{
    double rotAngle,cosRot,sinRot;
/*
   Compute angle for CW rotation of ps to radar
*/
    rotAngle = hAngle - xyAngle;
/*
    Compute range displacement vyra by ccw rot from ps to ra
*/      
    sinRot =  sin(rotAngle);
    cosRot =  cos(rotAngle);
    *drn = dxn * cosRot - dyn * sinRot;
    *dan = dxn * sinRot + dyn * cosRot;                  
}
