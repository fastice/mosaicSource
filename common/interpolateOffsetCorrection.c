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
