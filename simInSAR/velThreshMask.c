#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"
#include "gdalIO/gdalIO/grimpgdal.h"

/*
  velThreshMask -- compute a speed-threshold mask from a velocity map.

  For each pixel: mask = 1 if sqrt(vx^2+vy^2) <= velThresh, else 0.
  Nodata pixels (fill sentinel) are set to 0.
  Output is a float32 GeoTIFF in standard upper-left / negative-dy orientation.
  saveAsGeotiff handles the vertical flip from xyVEL's lower-left convention.
*/
void velThreshMask(xyVEL *xyVel, float velThresh, char *outputFile)
{
    int32_t xSize = xyVel->xSize;
    int32_t ySize = xyVel->ySize;
    int32_t i, j;
    float vx, vy, v, mask;
    int32_t nValid = 0, nMasked = 0, nNodata = 0;

    float *buffer = (float *)malloc((size_t)xSize * ySize * sizeof(float));
    if (buffer == NULL)
        error("velThreshMask: failed to allocate mask buffer (%d x %d)\n", xSize, ySize);

    for (i = 0; i < ySize; i++) {
        for (j = 0; j < xSize; j++) {
            vx = xyVel->vx[i][j];
            vy = xyVel->vy[i][j];
            if (vx <= -(float)LARGEINT || vy <= -(float)LARGEINT) {
                mask = 0.0f;
                nNodata++;
            } else {
                v = hypotf(vx, vy);
                if (v > velThresh) {
                    mask = 0.0f;
                    nMasked++;
                } else {
                    mask = 1.0f;
                    nValid++;
                }
            }
            buffer[i * xSize + j] = mask;
        }
    }

    fprintf(stderr, "velThreshMask: %d x %d  velThresh=%.1f m/yr\n", xSize, ySize, velThresh);
    fprintf(stderr, "  valid (speed<=thresh): %d  masked (speed>thresh): %d  nodata: %d\n",
            nValid, nMasked, nNodata);

    /* Geotransform in metres; xyVEL origin is lower-left so y_origin = y0 + ySize*deltaY.
       saveAsGeotiff flips rows vertically before writing. */
    double geoTransform[6];
    geoTransform[0] = xyVel->x0 * KMTOM;
    geoTransform[1] = xyVel->deltaX * KMTOM;
    geoTransform[2] = 0.0;
    geoTransform[3] = (xyVel->y0 + (double)ySize * xyVel->deltaY) * KMTOM;
    geoTransform[4] = 0.0;
    geoTransform[5] = -xyVel->deltaY * KMTOM;

    const char *epsg = getEPSGFromProjectionParams(xyVel->rot, xyVel->stdLat, xyVel->hemisphere);

    saveAsGeotiff(outputFile, buffer, xSize, ySize, geoTransform,
                  epsg, NULL, "GTiff", GDT_Float32, (float)NAN);

    free(buffer);
    fprintf(stderr, "velThreshMask: wrote %s\n", outputFile);
}
