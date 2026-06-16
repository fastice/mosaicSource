#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "tiePoints.h"
#include <stdlib.h>
#include <math.h>
#include "gdal.h"

double interpolatePhase(double range, double azimuth, unwrapPhaseStructure phaseImage);
/*
   Input phase image and extract phases for tiepoint locations.
*/
void getPhases(char *phaseFile, tiePointsStructure *tiePoints,
               inputImageStructure inputImage, int32_t yamlOutput, int32_t verbose)
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
        Input image — use GDAL for .vrt/.tif, raw freadBS otherwise.
    */
    if (strstr(phaseFile, ".vrt") != NULL || strstr(phaseFile, ".tif") != NULL)
    {
        GDALDatasetH hDS = GDALOpen(phaseFile, GA_ReadOnly);
        if (hDS == NULL)
            error("getPhases: GDALOpen failed for %s\n", phaseFile);
        GDALRasterBandH hBand = GDALGetRasterBand(hDS, 1);
        for (i = 0; i < inputImage.azimuthSize; i++)
        {
            if (GDALRasterIO(hBand, GF_Read, 0, i, inputImage.rangeSize, 1,
                             phaseImage.phase[i], inputImage.rangeSize, 1,
                             GDT_Float32, 0, 0) != CE_None)
                error("getPhases: GDALRasterIO failed at line %d of %s\n", i, phaseFile);
        }
        GDALClose(hDS);
    }
    else
    {
        fp = openInputFile(phaseFile);
        for (i = 0; i < inputImage.azimuthSize; i++)
            freadBS(phaseImage.phase[i], inputImage.rangeSize, sizeof(float), fp, FLOAT32FLAG);
    }
    /*
        Interpolate phases.
    */
    if (!yamlOutput && verbose)
        fprintf(stdout, ";;\n;;Tiepoints row column elevation\n;;\n");
    for (i = 0; i < tiePoints->npts; i++)
    {
        tiePoints->phase[i] = interpolatePhase(tiePoints->r[i], tiePoints->a[i], phaseImage);
        if (!yamlOutput && verbose && tiePoints->phase[i] > (-LARGEINT + 10))
            fprintf(stdout, "; %i  %i  %f %f\n",
                    (int)(tiePoints->r[i] + 0.5), (int)(tiePoints->a[i] + 0.5), tiePoints->z[i], tiePoints->phase[i]);
    }
    for (i = 0; i < inputImage.azimuthSize; i++)
        free(phaseImage.phase[i]);
    if (!yamlOutput && verbose)
        fprintf(stdout, ";&\n");
    return;
}

double interpolatePhase(double range, double azimuth, unwrapPhaseStructure phaseImage)
{
    return (double)bilinearInterp(phaseImage.phase, range, azimuth, phaseImage.rangeSize, phaseImage.azimuthSize,
                                  (float)(-LARGEINT + 10000), (float)(-LARGEINT));
}
