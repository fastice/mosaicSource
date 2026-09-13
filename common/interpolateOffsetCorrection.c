#include "mosaicSource/common/common.h"

/*
 * Interpolate the ionospheric range-offset correction at a given SLC
 * range/azimuth coordinate.
 *
 * range, azimuth : SLC pixel coordinates (same system as offsets->rO, deltaR)
 * Returns the interpolated correction value, or noData if out of bounds or
 * any contributing pixel is <= minValue.
 */
float interpolateOffsetIonCorrectionInPixels(offsetCorrection *corr,
                                  double range, double azimuth,
                                  float minValue, float noData)
{
    double r, a;

    r = (range   - corr->rO) / corr->deltaR;
    a = (azimuth - corr->aO) / corr->deltaA;

    return bilinearInterp(corr->rangeOffsetCorrection, r, a,
                          corr->nr, corr->na, minValue, noData);
}

/*
 * Azimuth counterpart: identical sampling, but reads the azimuth slot.
 * Kept as a separate function rather than a flag so the range path stays
 * bit-identical.
 */
float interpolateAzOffsetIonCorrectionInPixels(offsetCorrection *corr,
                                  double range, double azimuth,
                                  float minValue, float noData)
{
    double r, a;

    r = (range   - corr->rO) / corr->deltaR;
    a = (azimuth - corr->aO) / corr->deltaA;

    return bilinearInterp(corr->azimuthOffsetCorrection, r, a,
                          corr->nr, corr->na, minValue, noData);
}

/*
 * Azimuth ionosphere correction at a multi-look (range, azimuth) point, in
 * METRES, ready to be added to the interpolated azimuth offset.
 *
 * Returns 0.0 when no correction is loaded or the screen has no data here, so
 * the caller needs no guard of its own.  Sign convention: ADD, matching the
 * range screen (root CLAUDE.md "Ionosphere range correction sign convention").
 */
double azIonCorrectionMeters(Offsets *offsets, inputImageStructure *myImg,
                             double range, double azimuth, double azSLPixSize)
{
    double rangeSLC, azimuthSLC;
    float c;

    if (offsets->aOffCorrection.azimuthOffsetCorrection == NULL)
        return 0.0;
    computeSLCFromMLCoords(myImg, range, azimuth, &rangeSLC, &azimuthSLC);
    c = interpolateAzOffsetIonCorrectionInPixels(&(offsets->aOffCorrection),
                                                 rangeSLC, azimuthSLC,
                                                 -LARGEINT, 0.0);
    if (c <= -0.98 * LARGEINT)
        return 0.0;
    return (double)c * azSLPixSize;
}
