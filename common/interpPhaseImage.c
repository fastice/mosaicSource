#include "common.h"

void interpPhaseImage(inputImageStructure *inputImage, double range, double azimuth, double *phase)
{
	/*
	  If out of input image, return 0
	*/
	*phase = (double)bilinearInterp((float **)inputImage->image, range, azimuth, inputImage->rangeSize, inputImage->azimuthSize, -0.999 * (float)LARGEINT, (float)-LARGEINT);
}
