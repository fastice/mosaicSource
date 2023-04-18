#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "rparams.h"
#include <stdlib.h>
#include <math.h>
/*
  Input range offset and extract phases for tiepoint locations.
*/
void getROffsets(char *phaseFile, tiePointsStructure *tiePoints, inputImageStructure inputImage, Offsets *offsets)
{
	FILE *fp;
	double range, azimuth;
	int32_t i, count;
	/* Get offsets 	*/
	readRangeOffsets(offsets);
	/*
	  Interpolate offsets
	*/
	fprintf(stdout, ";;\n;;Tiepoints row column elevation\n;;\n");
	count = 0;
	for (i = 0; i < tiePoints->npts; i++)
	{
		range = (tiePoints->r[i] * inputImage.nRangeLooks - offsets->rO) / offsets->deltaR;
		azimuth = (tiePoints->a[i] * inputImage.nAzimuthLooks - offsets->aO) / offsets->deltaA;
		tiePoints->phase[i] = bilinearInterp((float **)offsets->dr, range, azimuth,
											 offsets->nr, offsets->na, -0.99 * LARGEINT, (float)-LARGEINT); /* only scale valide values */
		if (tiePoints->phase[i] > -0.98 * LARGEINT)
		{
			tiePoints->phase[i] *= inputImage.rangePixelSize / inputImage.nRangeLooks;
			count++;
		}
		/*  Multiply by -1 for left to get RHS*/
		/*         if(inputImage.lookDir==LEFT) tiePoints->phase[i] *=-1; */
		if (fabs(tiePoints->phase[i]) < LARGEINT / 10 && tiePoints->quiet == FALSE)
			fprintf(stdout, "; %i  %i  %f %f\n", (int)(tiePoints->r[i] + 0.5), (int)(tiePoints->a[i] + 0.5),
					tiePoints->z[i], tiePoints->phase[i]);
	}
	fprintf(stderr, "count %i\n", count);
	fprintf(stdout, ";&\n");
	return;
}
