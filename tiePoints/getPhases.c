#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "tiePoints.h"
#include <stdlib.h>
#include <math.h>

double interpolatePhase(double range, double azimuth, unwrapPhaseStructure phaseImage);
/*
   Input phase image and extract phases for tiepoint locations.
*/
void getPhases(char *phaseFile, tiePointsStructure *tiePoints,
               inputImageStructure inputImage)
{
    FILE *fp;
    unwrapPhaseStructure phaseImage;
    int32_t i;
    /*
       Init image
    */
    phaseImage.phase = (float **)malloc(inputImage.azimuthSize * sizeof(float *));
    for (i = 0; i < inputImage.azimuthSize; i++)
        phaseImage.phase[i] = (float *)malloc(inputImage.rangeSize * sizeof(float));
    phaseImage.rangeSize = inputImage.rangeSize;
    phaseImage.azimuthSize = inputImage.azimuthSize;
    /*
        Input image
    */
    fp = openInputFile(phaseFile);
    for (i = 0; i < inputImage.azimuthSize; i++)
        freadBS(phaseImage.phase[i], inputImage.rangeSize, sizeof(float), fp, FLOAT32FLAG);
    /*
        Interpolate phases.
    */
    fprintf(stdout, ";;\n;;Tiepoints row column elevation\n;;\n");
    for (i = 0; i < tiePoints->npts; i++)
    {

        tiePoints->phase[i] = interpolatePhase(tiePoints->r[i], tiePoints->a[i], phaseImage);
        if (tiePoints->phase[i] > (-LARGEINT + 10))
            fprintf(stdout, "; %i  %i  %f %f\n",
                    (int)(tiePoints->r[i] + 0.5), (int)(tiePoints->a[i] + 0.5), tiePoints->z[i], tiePoints->phase[i]);
    }
    for (i = 0; i < inputImage.azimuthSize; i++)
        free(phaseImage.phase[i]);
    fprintf(stdout, ";&\n");
    return;
}

double interpolatePhase(double range, double azimuth, unwrapPhaseStructure phaseImage)
{
    return (double)bilinearInterp(phaseImage.phase, range, azimuth, phaseImage.rangeSize, phaseImage.azimuthSize,
                                  (float)(-LARGEINT + 10000), (float)(-LARGEINT));
}
