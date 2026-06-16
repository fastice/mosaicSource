/*
 * readComplexAsPower.c
 *
 * Drop-in replacement for getMosaicInputImage that handles complex (SLC) VRT
 * input.  If the VRT data type is complex (GDT_CFloat32 etc.), the data is
 * read into a temporary buffer and converted to power (re^2 + im^2) before
 * storing in inputImage->image.  For real-valued VRTs and non-VRT files the
 * call falls through to getMosaicInputImage unchanged.
 *
 * Integration:
 *   1. In makeGeoMosaic.c, getGeoMosaicImage(), replace:
 *          getMosaicInputImage(inputImage, yMin, yMax);
 *      with:
 *          readComplexAsPower(inputImage, yMin, yMax);
 *
 *   2. Add prototype to geomosaic.h:
 *          void readComplexAsPower(inputImageStructure *inputImage,
 *                                  int32_t yMin, int32_t yMax);
 *
 *   3. Add readComplexAsPower.o to GEOMOSAIC in mosaicSource/Makefile and
 *      to SRCS in geoMosaic/Makefile (same pattern as subpixelRTC.c).
 */

#include <math.h>
#include <stdlib.h>
#include "mosaicSource/common/common.h"
#include "geomosaic.h"
#include "gdal.h"
#include "gdalIO/gdalIO/grimpgdal.h"

void readComplexAsPower(inputImageStructure *inputImage, int32_t yMin, int32_t yMax)
{
    char         vrtBuf[4096], *vrtFile;
    GDALDatasetH hDS;
    GDALRasterBandH hBand;
    int          gdal_type, xSize, ySize, status;
    int32_t      yMin_gdal, yMax_gdal, yMin_valid, yMax_valid, nRows;
    int64_t      i, j, nPix;
    float       *power;   /* points to inputImage->image[0] */
    float       *cplx;    /* temp: interleaved (re, im) for valid rows */

    vrtFile = checkForVrt(inputImage->file, vrtBuf);
    if (vrtFile == NULL) {
        getMosaicInputImage(inputImage, yMin, yMax);
        return;
    }

    /* Peek at the data type before allocating anything */
    hDS = GDALOpen(vrtFile, GDAL_OF_READONLY);
    if (hDS == NULL) {
        getMosaicInputImage(inputImage, yMin, yMax);
        return;
    }
    hBand    = GDALGetRasterBand(hDS, 1);
    gdal_type = GDALGetRasterDataType(hBand);
    xSize    = GDALGetRasterBandXSize(hBand);
    ySize    = GDALGetRasterBandYSize(hBand);

    if (!GDALDataTypeIsComplex(gdal_type)) {
        GDALClose(hDS);
        getMosaicInputImage(inputImage, yMin, yMax);
        return;
    }

    fprintf(stderr, "Complex VRT detected (%s) — converting to power\n", vrtFile);

    /* Convert azimuth bounds from image to VRT (SLC) row coordinates */
    yMin_gdal  = yMin / inputImage->nAzimuthLooks;
    yMax_gdal  = yMax / inputImage->nAzimuthLooks;
    yMin_valid = (yMin_gdal < 0)      ? 0        : yMin_gdal;
    yMax_valid = (yMax_gdal >= ySize) ? ySize - 1 : yMax_gdal;
    nRows      = yMax_valid - yMin_valid + 1;

    /* Initialise the power buffer (multilooked size, not SLC size) to no-data */
    nPix  = (int64_t)inputImage->azimuthSize * inputImage->rangeSize;
    power = (float *)inputImage->image[0];
    for (i = 0; i < nPix; i++)
        power[i] = (float)-LARGEINT;

    /* Temp buffer: only the valid rows, two floats per complex pixel */
    cplx = (float *)malloc((size_t)xSize * nRows * 2 * sizeof(float));
    if (cplx == NULL)
        error("readComplexAsPower: malloc failed (%i rows x %i cols)\n", nRows, xSize);

    /* Read valid rows as CFloat32 regardless of the source complex type */
    status = GDALRasterIO(hBand, GF_Read,
                          0, yMin_valid, xSize, nRows,
                          cplx, xSize, nRows,
                          GDT_CFloat32, 0, 0);
    GDALClose(hDS);
    if (status != CE_None)
        error("readComplexAsPower: GDALRasterIO failed for %s\n", vrtFile);

    /*
     * cplx layout: [re_00, im_00, re_01, im_01, ..., re_(nRows-1)(xSize-1), im_...]
     * Write power into the rows that were read; out-of-range rows stay as -LARGEINT.
     */
    for (i = 0; i < nRows; i++) {
        int64_t out_row = yMin_valid + i;
        for (j = 0; j < xSize; j++) {
            int64_t ci = (i * xSize + j) * 2;
            float re = cplx[ci];
            float im = cplx[ci + 1];
            power[out_row * xSize + j] = re * re + im * im;
        }
    }

    free(cplx);
    fprintf(stderr, "Complex to power done: %i rows x %i cols\n", nRows, xSize);
}
