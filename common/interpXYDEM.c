#include "common.h"

double interpXYDEM(double x, double y, xyDEM xydem)
{
	double xi, yi;
	int32_t i, j, iSize, jSize;
	double t, u, p1, p2, p3, p4;
	double z;

	xi = (x - xydem.x0) / xydem.deltaX;
	yi = (y - xydem.y0) / xydem.deltaY;
	return bilinearInterp(xydem.z, xi, yi, xydem.xSize, xydem.ySize, MINELEVATION, MINELEVATION);
}
