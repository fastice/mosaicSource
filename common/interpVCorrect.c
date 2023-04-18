#include "mosaicSource/common/common.h"

double interpVCorrect(double x, double y, xyDEM *vCorrect)
{
	double xi, yi;
	int32_t i, j, iSize, jSize;
	double t, u, p1, p2, p3, p4;
	double dzdt;

	xi = (x - vCorrect->x0) / vCorrect->deltaX;
	yi = (y - vCorrect->y0) / vCorrect->deltaY;
	dzdt = bilinearInterp(vCorrect->z, xi, yi, vCorrect->xSize, vCorrect->ySize, MINVCORRECT, MINVCORRECT);
	if (dzdt <= MINVCORRECT)
		dzdt = 0.0;
	return dzdt;
}