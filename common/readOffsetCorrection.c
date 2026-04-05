#include "stdio.h"
#include "string.h"
#include "math.h"
#include "stdlib.h"
#include "mosaicSource/common/common.h"
#include "gdalIO/gdalIO/grimpgdal.h"

/*
 * Lazily allocate the memory-pool buffers for ionospheric correction data.
 * correctionBuf1/2 hold the float** row-pointer arrays (one entry per azimuth line).
 * Correction1/2 hold the flat float data arrays.
 * Called before first use of each buffer slot; subsequent calls are no-ops.
 */
void mallocCorrectionBuffer(int bufferMode)
{
    extern void *correctionBuf1, *correctionBuf2, *Correction1, *Correction2;

    if (bufferMode == RANGEBUFF && correctionBuf1 == NULL)
    {
        correctionBuf1 = malloc(sizeof(float *) * MAXOFFLENGTH);
        Correction1    = malloc(MAXOFFBUF);
        return;
    }

    if (bufferMode == RANGEUSEAZIMUTHBUFF && correctionBuf2 == NULL)
    {
        correctionBuf2 = malloc(sizeof(float *) * MAXOFFLENGTH);
        Correction2    = malloc(MAXOFFBUF);
        return;
    }
}

/*
 * Read an ionospheric range-offset correction from a VRT file into
 * offsets->rOffCorrection.
 *
 * The VRT metadata must contain:
 *   NumberRangeLooks   -> stored as rOffCorrection.deltaR
 *   NumberAzimuthLooks -> stored as rOffCorrection.deltaA
 *
 * Raster dimensions (xSize, ySize) give nr and na.
 * Origin (r0, a0) is set to deltaR/2, deltaA/2 (centre of first pixel).
 *
 * Data is placed into the pre-allocated pool:
 *   RANGEBUFF          -> correctionBuf1 / Correction1
 *   RANGEUSEAZIMUTHBUFF -> correctionBuf2 / Correction2
 */
void readOffsetCorrection(char *correctionFile, Offsets *offsets, int bufferMode)
{
    extern void *correctionBuf1, *correctionBuf2, *Correction1, *Correction2;
    GDALDatasetH   hDS;
    GDALRasterBandH hBand;
    dictNode *metaData = NULL;
    float **rangeOffsetCorrection;
    float  *data;
    float   deltaR, deltaA;
    int32_t nr, na;
    int     i, status;

    /* Lazy-init the pool for this slot */
    mallocCorrectionBuffer(bufferMode);

    hDS = GDALOpen(correctionFile, GDAL_OF_READONLY);
    if (hDS == NULL)
        error("readOffsetCorrection: cannot open %s\n", correctionFile);

    hBand = GDALGetRasterBand(hDS, 1);
    nr    = GDALGetRasterBandXSize(hBand);
    na    = GDALGetRasterBandYSize(hBand);

    if ((int64_t)nr * na * (int64_t)sizeof(float) > MAXOFFBUF)
        error("readOffsetCorrection: %s exceeds MAXOFFBUF (nr=%i na=%i)\n",
              correctionFile, nr, na);
    if (na > MAXOFFLENGTH)
        error("readOffsetCorrection: %s na=%i exceeds MAXOFFLENGTH=%i\n",
              correctionFile, na, MAXOFFLENGTH);

    readDataSetMetaData(hDS, &metaData);
    deltaR = (float)atof(get_value(metaData, "NumberRangeLooks"));
    deltaA = (float)atof(get_value(metaData, "NumberAzimuthLooks"));

    offsets->rOffCorrection.nr     = nr;
    offsets->rOffCorrection.na     = na;
    offsets->rOffCorrection.rO     = (int32_t)(deltaR / 2);
    offsets->rOffCorrection.aO     = (int32_t)(deltaA / 2);
    offsets->rOffCorrection.deltaR = deltaR;
    offsets->rOffCorrection.deltaA = deltaA;

    if (bufferMode == RANGEBUFF)
    {
        rangeOffsetCorrection = (float **)correctionBuf1;
        data                  = (float *)Correction1;
    }
    else if (bufferMode == RANGEUSEAZIMUTHBUFF)
    {
        rangeOffsetCorrection = (float **)correctionBuf2;
        data                  = (float *)Correction2;
    }
    else
        error("readOffsetCorrection: invalid bufferMode %d\n", bufferMode);

    for (i = 0; i < na; i++)
        rangeOffsetCorrection[i] = &data[i * nr];
    offsets->rOffCorrection.rangeOffsetCorrection = rangeOffsetCorrection;

    status = GDALRasterIO(hBand, GF_Read, 0, 0, nr, na, data,
                          nr, na, GDT_Float32, 0, 0);
    if (status != CE_None)
        error("readOffsetCorrection: GDALRasterIO failed for %s\n", correctionFile);

    /* GDAL replaces VRT nodata pixels with NaN; convert back to -LARGEINT sentinel */
    for (i = 0; i < nr * na; i++)
        if (isnan(data[i])) data[i] = (float)-LARGEINT;

    GDALClose(hDS);
    fprintf(stderr, "readOffsetCorrection: %s  nr=%i na=%i deltaR=%.1f deltaA=%.1f r0=%i a0=%i\n",
            correctionFile, nr, na, deltaR, deltaA,
            offsets->rOffCorrection.rO, offsets->rOffCorrection.aO);
}
