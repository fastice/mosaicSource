#include "common.h"
/*
  General bilinear interp fimage[azimuth][range] with size nr,na. If any of the four points <=minvalue, then return noData.
*/
float bilinearInterp(float **fimage, double range, double azimuth, int32_t nr, int32_t na, float minvalue, float noData)
{
	int32_t i, j;
	float t, u;
	float p1, p2, p3, p4;
	float result;

	j = (int)range;
	i = (int)azimuth;
	/*
		Return out of range
	*/
	if (range < 0.0 || azimuth < 0.0 || j >= nr || i >= na)
		return noData;
	/*
	  Handle border pixels
	*/
	if ((j == (nr - 1)) && (i >= 0 && i < na))
		return fimage[i][j];
	if ((i == (na - 1)) && (j >= 0 && j < nr))
		return fimage[i][j];
	/*
	  Interpolate interior pixels
	*/
	t = (float)(range - (double)j);
	u = (float)(azimuth - (double)i);
	p1 = fimage[i][j];
	p2 = fimage[i][j + 1];
	p3 = fimage[i + 1][j + 1];
	p4 = fimage[i + 1][j];
	/* if any pixels no data, return nodata - avoids crud on border */
	if (p1 <= minvalue || p2 <= minvalue || p3 <= minvalue || p4 <= minvalue)
		return noData;
	/* final interp */
	result = (float)((1.0 - t) * (1.0 - u) * p1 + t * (1.0 - u) * p2 + t * u * p3 + (1.0 - t) * u * p4);
	return result;
}
