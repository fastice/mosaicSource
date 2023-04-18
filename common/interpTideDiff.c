#include "common.h"

double interpTideDiff(double x, double y, xyDEM xydem)
{
	double xi, yi;
	int32_t i, j, iSize, jSize;
	double t, u, p1, p2, p3, p4;
	double z;

	xi = (x - xydem.x0) / xydem.deltaX;
	yi = (y - xydem.y0) / xydem.deltaY;

	j = (int)xi;
	i = (int)yi;
	t = xi - (double)j;
	u = yi - (double)i;
	iSize = xydem.ySize;
	jSize = xydem.xSize;
	/*
	  If location of DEM return 0.
	*/
	if (i < 0 || i >= iSize || j < 0 || j >= jSize)
		return MINELEVATION;
	/*
	  Border pixels
	*/
	if (i == (iSize - 1) || j == (jSize - 1))
		return xydem.z[i][j];
	/*
	  Set up interpolation
	*/
	p1 = (double)xydem.z[i][j];
	p2 = (double)xydem.z[i][j + 1];
	p3 = (double)xydem.z[i + 1][j + 1];
	p4 = (double)xydem.z[i + 1][j];
	if (p1 <= MINELEVATION || p2 <= MINELEVATION ||
		p3 <= MINELEVATION || p4 <= MINELEVATION)
	{
		return max(max(max(p1, p2), p3), p4);
	}

	/*
	  Use bilinear interpolation to compute height.
	*/
	z = ((1.0 - t) * (1.0 - u) * p1 + t * (1.0 - u) * p2 + t * u * p3 + (1.0 - t) * u * p4);

	return z;
}
