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
 * Read an ionospheric range-offset correction from a VRT/GeoTIFF into
 * offsets->rOffCorrection.
 *
 * The correction file is on the native ROFF grid — the same pixel spacing
 * and origin as the range offsets (r0, a0, deltaR, deltaA in SLC pixels).
 * Those values are copied directly from the already-populated offsets struct
 * rather than re-reading them from file metadata.
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
    float **rangeOffsetCorrection;
    float  *data;
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

    /* Correction is on the ROFF grid — inherit r0/a0/deltaR/deltaA from the
       range offsets (already populated from range.offsets.vrt metadata). */
    offsets->rOffCorrection.nr     = nr;
    offsets->rOffCorrection.na     = na;
    offsets->rOffCorrection.rO     = offsets->rO;
    offsets->rOffCorrection.aO     = offsets->aO;
    offsets->rOffCorrection.deltaR = offsets->deltaR;
    offsets->rOffCorrection.deltaA = offsets->deltaA;

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
            correctionFile, nr, na, offsets->rOffCorrection.deltaR, offsets->rOffCorrection.deltaA,
            offsets->rOffCorrection.rO, offsets->rOffCorrection.aO);
}
