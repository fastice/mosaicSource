#include "common.h"

extern int32_t indentRegionOutput;

/*
  Update azimuthMin/azimuthMax from a single output-grid pixel (i, j).
  Does nothing if the pixel does not project into the SAR image.
*/
static void samplePoint(int i, int j,
                        outputImageStructure *outputImage,
                        inputImageStructure  *currentImage,
                        double maxA, double maxR,
                        float *azimuthMin, float *azimuthMax)
{
    extern int    HemiSphere;
    extern double Rotation;
    double lat, lon, azimuth, range;
    double y = (outputImage->originY + i * outputImage->deltaY) * MTOKM;
    double x = (outputImage->originX + j * outputImage->deltaX) * MTOKM;
    xytoll1(x, y, HemiSphere, &lat, &lon, Rotation, outputImage->slat);
    llToImageNew(lat, lon, 0, &range, &azimuth, currentImage);
    if (azimuth >= 0.0 && azimuth <= maxA && range >= 0.0 && range <= maxR)
    {
        if (azimuth < *azimuthMin) *azimuthMin = (float)azimuth;
        if (azimuth > *azimuthMax) *azimuthMax = (float)azimuth;
    }
}

/*
  Find the azimuth limits (in single-look pixels) of the portion of
  currentImage that intersects the output grid region [imin,imax) x [jmin,jmax).

  Samples both the perimeter of the box (where swath entry/exit points lie
  for a diagonal track) and an interior grid, so that oblique swaths are
  handled correctly.
*/
void getAzimuthBoundsForXYBox(int32_t imin, int32_t imax, int32_t jmin, int32_t jmax,
                               inputImageStructure  *currentImage,
                               outputImageStructure *outputImage,
                               float *azimuthMin, float *azimuthMax)
{
    *azimuthMin =  1e30f;
    *azimuthMax = -1e30f;

    double maxR = (double)(currentImage->rangeSize   - 1);
    double maxA = (double)(currentImage->azimuthSize - 1);

    /* Step size for interior grid — at least 1 to avoid infinite loops */
    int delta = max(1, (int)(12e3 / outputImage->deltaX));
    /* might want to uncomment for debugging later
    fprintf(stderr, "%sgetAzimuthBoundsForXYBox: imin=%i imax=%i jmin=%i jmax=%i delta=%i\n",
            indentRegionOutput ? "\t" : "", imin, imax, jmin, jmax, delta);
    */

    /* --- Sample the perimeter of the box ---
       This is where the swath boundary crosses the output box, so the
       azimuth extremes of the intersection are found here for oblique tracks. */

    /* Top and bottom edges (full column range) */
    for (int j = jmin; j <= jmax; j += delta)
    {
        samplePoint(imin, j, outputImage, currentImage, maxA, maxR, azimuthMin, azimuthMax);
        samplePoint(imax, j, outputImage, currentImage, maxA, maxR, azimuthMin, azimuthMax);
    }
    /* Always hit the corners */
    samplePoint(imin, jmax, outputImage, currentImage, maxA, maxR, azimuthMin, azimuthMax);
    samplePoint(imax, jmax, outputImage, currentImage, maxA, maxR, azimuthMin, azimuthMax);

    /* Left and right edges (full row range) */
    for (int i = imin; i <= imax; i += delta)
    {
        samplePoint(i, jmin, outputImage, currentImage, maxA, maxR, azimuthMin, azimuthMax);
        samplePoint(i, jmax, outputImage, currentImage, maxA, maxR, azimuthMin, azimuthMax);
    }
    /* Always hit the corners */
    samplePoint(imax, jmin, outputImage, currentImage, maxA, maxR, azimuthMin, azimuthMax);
    samplePoint(imax, jmax, outputImage, currentImage, maxA, maxR, azimuthMin, azimuthMax);

    /* --- Sample interior grid ---
       Catches cases where the swath covers the whole box interior. */
    for (int i = imin; i <= imax; i += delta)
        for (int j = jmin; j <= jmax; j += delta)
            samplePoint(i, j, outputImage, currentImage, maxA, maxR, azimuthMin, azimuthMax);

    if (*azimuthMin > *azimuthMax)
    {
        /* No intersection found */
        *azimuthMin = 0.0f;
        *azimuthMax = 0.0f;
        fprintf(stderr, "%sgetAzimuthBoundsForXYBox: no intersection\n", indentRegionOutput ? "\t" : "");
        return;
    }

    /* Add pad and convert multilook -> single-look pixels */
    float pad = (float)(55e3 / currentImage->azimuthPixelSize);
    *azimuthMin = max(*azimuthMin - pad, 0.0f)                              * currentImage->nAzimuthLooks;
    *azimuthMax = min(*azimuthMax + pad, (float)(currentImage->azimuthSize - 1)) * currentImage->nAzimuthLooks;

    /* might want to uncomment for debugging later
    fprintf(stderr, "%sgetAzimuthBoundsForXYBox: azimuthMin=%.1f azimuthMax=%.1f (single-look pixels)\n",
            indentRegionOutput ? "\t" : "", *azimuthMin, *azimuthMax);
    */
}
