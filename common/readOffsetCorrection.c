#include "stdio.h"
#include "string.h"
#include "math.h"
#include "stdlib.h"
#include "unistd.h"
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
    extern void *correctionBuf1, *correctionBuf2, *correctionBuf3;
    extern void *Correction1, *Correction2, *Correction3;

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

    if (bufferMode == AZIMUTHIONBUFF && correctionBuf3 == NULL)
    {
        correctionBuf3 = malloc(sizeof(float *) * MAXOFFLENGTH);
        Correction3    = malloc(MAXOFFBUF);
        return;
    }
}

/*
 * Load the ionosphere fill-error map that nisargrimpworkflow writes next to a range
 * correction <stem>.ionosphereCorrection.globalFill.offset.vrt as
 * <stem>.ionosphereCorrection.sigma.tif (ionoSigmaVersion sigma-1).  Optional: a missing file,
 * an unrecognised version, or useIonSigma == FALSE leaves corr->ionSigma NULL (no change to
 * the range sigma).  One small buffer per range slot, reused like the correction pool.
 */
static void readIonSigma(char *correctionFile, offsetCorrection *corr, int slot)
{
    extern int32_t useIonSigma;
    static float **rows[2] = {NULL, NULL};
    static float *buf[2] = {NULL, NULL};
    static size_t capRows[2] = {0, 0}, capBuf[2] = {0, 0};
    const char *key = ".ionosphereCorrection.globalFill.offset.vrt";
    char sigmaFile[2048];
    char *p;
    const char *ver, *item;
    GDALDatasetH hDS;
    int32_t nr, na, i;

    corr->ionSigma = NULL;
    if (useIonSigma == FALSE)
    {
        return;
    }
    strncpy(sigmaFile, correctionFile, sizeof(sigmaFile) - 1);
    sigmaFile[sizeof(sigmaFile) - 1] = '\0';
    p = strstr(sigmaFile, key);
    if (p == NULL || strlen(sigmaFile) - strlen(key) + 40 >= sizeof(sigmaFile))
    {
        return;
    }
    strcpy(p, ".ionosphereCorrection.sigma.tif");
    if (access(sigmaFile, F_OK) != 0)
    {
        return;
    }
    hDS = GDALOpen(sigmaFile, GDAL_OF_READONLY);
    if (hDS == NULL)
    {
        return;
    }
    ver = GDALGetMetadataItem(hDS, "ionoSigmaVersion", NULL);
    if (ver == NULL || strcmp(ver, "sigma-1") != 0)
    {
        fprintf(stderr, "readIonSigma: %s has unrecognised ionoSigmaVersion %s -- ignored\n",
                sigmaFile, ver == NULL ? "(none)" : ver);
        GDALClose(hDS);
        return;
    }
    nr = GDALGetRasterXSize(hDS);
    na = GDALGetRasterYSize(hDS);
    if ((size_t)nr * na > capBuf[slot])
    {
        free(buf[slot]);
        buf[slot] = (float *)malloc(sizeof(float) * (size_t)nr * na);
        capBuf[slot] = (size_t)nr * na;
    }
    if ((size_t)na > capRows[slot])
    {
        free(rows[slot]);
        rows[slot] = (float **)malloc(sizeof(float *) * na);
        capRows[slot] = na;
    }
    for (i = 0; i < na; i++)
    {
        rows[slot][i] = &buf[slot][(size_t)i * nr];
    }
    if (GDALRasterIO(GDALGetRasterBand(hDS, 1), GF_Read, 0, 0, nr, na, buf[slot], nr, na,
                     GDT_Float32, 0, 0) != CE_None)
    {
        GDALClose(hDS);
        return;
    }
    item = GDALGetMetadataItem(hDS, "r0", NULL);      corr->sR0 = (item != NULL) ? atof(item) : corr->rO;
    item = GDALGetMetadataItem(hDS, "a0", NULL);      corr->sA0 = (item != NULL) ? atof(item) : corr->aO;
    item = GDALGetMetadataItem(hDS, "deltaR", NULL);  corr->sDeltaR = (item != NULL) ? atof(item) : corr->deltaR;
    item = GDALGetMetadataItem(hDS, "deltaA", NULL);  corr->sDeltaA = (item != NULL) ? atof(item) : corr->deltaA;
    GDALClose(hDS);
    corr->sNr = nr;
    corr->sNa = na;
    corr->ionSigma = rows[slot];
    fprintf(stderr, "Ionosphere fill-error map (%s): %s\n", ver, sigmaFile);
}

/*
 * Read an ionospheric offset correction from a VRT/GeoTIFF into
 * offsets->rOffCorrection (range) or offsets->aOffCorrection (azimuth).
 *
 * The correction file is on the native ROFF grid — the same pixel spacing
 * and origin as the range offsets (r0, a0, deltaR, deltaA in SLC pixels).
 * Those values are copied directly from the already-populated offsets struct
 * rather than re-reading them from file metadata.
 *
 * Data is placed into the pre-allocated pool:
 *   RANGEBUFF           -> correctionBuf1 / Correction1  (range)
 *   RANGEUSEAZIMUTHBUFF -> correctionBuf2 / Correction2  (range, crossing orbit)
 *   AZIMUTHIONBUFF      -> correctionBuf3 / Correction3  (azimuth)
 *
 * The azimuth screen is in SLC azimuth pixels and is ADDED to the azimuth
 * offset, the same sign convention the range screen uses (see the root
 * CLAUDE.md "Ionosphere range correction sign convention").
 */
void readOffsetCorrection(char *correctionFile, Offsets *offsets, int bufferMode)
{
    extern void *correctionBuf1, *correctionBuf2, *correctionBuf3;
    extern void *Correction1, *Correction2, *Correction3;
    offsetCorrection *corr;
    GDALDatasetH   hDS;
    GDALRasterBandH hBand;
    float **offsetCorrectionData;
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

    corr = (bufferMode == AZIMUTHIONBUFF) ? &(offsets->aOffCorrection)
                                          : &(offsets->rOffCorrection);

    /* Correction is on the ROFF grid — inherit r0/a0/deltaR/deltaA from the
       offsets (already populated from the offsets VRT metadata). */
    corr->nr     = nr;
    corr->na     = na;
    corr->rO     = offsets->rO;
    corr->aO     = offsets->aO;
    corr->deltaR = offsets->deltaR;
    corr->deltaA = offsets->deltaA;

    if (bufferMode == RANGEBUFF)
    {
        offsetCorrectionData = (float **)correctionBuf1;
        data                 = (float *)Correction1;
    }
    else if (bufferMode == RANGEUSEAZIMUTHBUFF)
    {
        offsetCorrectionData = (float **)correctionBuf2;
        data                 = (float *)Correction2;
    }
    else if (bufferMode == AZIMUTHIONBUFF)
    {
        offsetCorrectionData = (float **)correctionBuf3;
        data                 = (float *)Correction3;
    }
    else
        error("readOffsetCorrection: invalid bufferMode %d\n", bufferMode);

    for (i = 0; i < na; i++)
        offsetCorrectionData[i] = &data[i * nr];
    if (bufferMode == AZIMUTHIONBUFF)
        corr->azimuthOffsetCorrection = offsetCorrectionData;
    else
        corr->rangeOffsetCorrection = offsetCorrectionData;

    status = GDALRasterIO(hBand, GF_Read, 0, 0, nr, na, data,
                          nr, na, GDT_Float32, 0, 0);
    if (status != CE_None)
        error("readOffsetCorrection: GDALRasterIO failed for %s\n", correctionFile);

    /* GDAL replaces VRT nodata pixels with NaN; convert back to -LARGEINT sentinel */
    for (i = 0; i < nr * na; i++)
        if (isnan(data[i])) data[i] = (float)-LARGEINT;

    GDALClose(hDS);
    if (bufferMode != AZIMUTHIONBUFF)
    {
        readIonSigma(correctionFile, corr, bufferMode == RANGEBUFF ? 0 : 1);
    }
}
