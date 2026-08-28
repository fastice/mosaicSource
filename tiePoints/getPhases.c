#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "tiePoints.h"
#include <stdlib.h>
#include <math.h>
#include "gdal.h"

double interpolatePhase(double range, double azimuth, unwrapPhaseStructure phaseImage);
static void readRadarImage(char *file, inputImageStructure inputImage,
                           unwrapPhaseStructure *image, char *caller);
static void freeRadarImage(unwrapPhaseStructure *image);

/*
   Read a full-size radar-grid float image (phase or ionosphere) matching inputImage's
   dimensions. Uses GDAL for .vrt/.tif, raw byte-swapped freadBS otherwise. Factored out
   of getPhases() so getIonosphere() can read its file exactly the same way.
*/
static void readRadarImage(char *file, inputImageStructure inputImage,
                           unwrapPhaseStructure *image, char *caller)
{
    FILE *fp;
    int32_t i;

    image->phase = (float **)malloc(inputImage.azimuthSize * sizeof(float *));
    for (i = 0; i < inputImage.azimuthSize; i++)
    {
        image->phase[i] = (float *)malloc(inputImage.rangeSize * sizeof(float));
    }
    image->rangeSize = inputImage.rangeSize;
    image->azimuthSize = inputImage.azimuthSize;

    if (strstr(file, ".vrt") != NULL || strstr(file, ".tif") != NULL)
    {
        GDALDatasetH hDS = GDALOpen(file, GA_ReadOnly);
        if (hDS == NULL)
            error("%s: GDALOpen failed for %s\n", caller, file);
        if (GDALGetRasterXSize(hDS) != inputImage.rangeSize ||
            GDALGetRasterYSize(hDS) != inputImage.azimuthSize)
            error("%s: %s is %d x %d but geodat gives %d x %d\n", caller, file,
                  GDALGetRasterXSize(hDS), GDALGetRasterYSize(hDS),
                  inputImage.rangeSize, inputImage.azimuthSize);
        GDALRasterBandH hBand = GDALGetRasterBand(hDS, 1);
        for (i = 0; i < inputImage.azimuthSize; i++)
        {
            if (GDALRasterIO(hBand, GF_Read, 0, i, inputImage.rangeSize, 1,
                             image->phase[i], inputImage.rangeSize, 1,
                             GDT_Float32, 0, 0) != CE_None)
                error("%s: GDALRasterIO failed at line %d of %s\n", caller, i, file);
        }
        GDALClose(hDS);
    }
    else
    {
        fp = openInputFile(file);
        for (i = 0; i < inputImage.azimuthSize; i++)
        {
            freadBS(image->phase[i], inputImage.rangeSize, sizeof(float), fp, FLOAT32FLAG);
        }
    }
}

static void freeRadarImage(unwrapPhaseStructure *image)
{
    int32_t i;
    for (i = 0; i < image->azimuthSize; i++)
    {
        free(image->phase[i]);
    }
    free(image->phase);
    image->phase = NULL;
}

/*
   Input phase image and extract phases for tiepoint locations.
*/
void getPhases(char *phaseFile, tiePointsStructure *tiePoints,
               inputImageStructure inputImage, int32_t yamlOutput, int32_t verbose)
{
    unwrapPhaseStructure phaseImage;
    int32_t i;
    /*
        Input image — use GDAL for .vrt/.tif, raw freadBS otherwise.
    */
    readRadarImage(phaseFile, inputImage, &phaseImage, "getPhases");
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
    freeRadarImage(&phaseImage);
    if (!yamlOutput && verbose)
        fprintf(stdout, ";&\n");
    return;
}

/*
   Input an ionospheric phase image (radians, same multilooked grid as the unwrapped phase)
   and extract its value at each tiepoint location into tiePoints->ionPhase.

   The file holds the ionospheric phase itself, in the same sign convention as the phase
   image it accompanies, so consumers SUBTRACT it (phaseCorrected = phase - ion). This is
   ISCE's own convention -- topsApp applies its screen as a*exp(-1.0*J*b) with b =
   topophase.ion (runMergeBursts.py). Note this is the opposite of the range-offset
   ionosphere correction used by rparams/mosaic3d, which is a pre-negated correction that
   gets added.
*/
void getIonosphere(char *ionFile, tiePointsStructure *tiePoints, inputImageStructure inputImage)
{
    unwrapPhaseStructure ionImage;
    int32_t i;

    readRadarImage(ionFile, inputImage, &ionImage, "getIonosphere");
    tiePoints->ionPhase = (double *)malloc(tiePoints->npts * sizeof(double));
    for (i = 0; i < tiePoints->npts; i++)
    {
        tiePoints->ionPhase[i] = interpolatePhase(tiePoints->r[i], tiePoints->a[i], ionImage);
    }
    tiePoints->hasIonosphere = TRUE;
    freeRadarImage(&ionImage);
    return;
}

double interpolatePhase(double range, double azimuth, unwrapPhaseStructure phaseImage)
{
    return (double)bilinearInterp(phaseImage.phase, range, azimuth, phaseImage.rangeSize, phaseImage.azimuthSize,
                                  (float)(-LARGEINT + 10000), (float)(-LARGEINT));
}
