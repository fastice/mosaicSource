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

/*
 * Variance (m^2 of slant range) of the applied range ionosphere correction from the
 * fill-error map, bilinearly interpolated at an SLC coordinate; 0 when there is no map.
 * Cell j averages offset samples j*n .. j*n+n-1 (n = sDeltaR / deltaR), so its centre sits
 * (n-1)/2 samples past the cell's first sample; outside the map the nearest edge is used.
 */
double rangeIonSigmaVarMeters(offsetCorrection *corr, double rSLC, double aSLC, float rSLPixSize)
{
    double r, a;
    float s;

    if (corr->ionSigma == NULL)
    {
        return 0.0;
    }
    r = (rSLC - corr->sR0) / corr->sDeltaR - 0.5 * (1.0 - corr->deltaR / corr->sDeltaR);
    a = (aSLC - corr->sA0) / corr->sDeltaA - 0.5 * (1.0 - corr->deltaA / corr->sDeltaA);
    r = fmin(fmax(r, 0.0), corr->sNr - 1.001);
    a = fmin(fmax(a, 0.0), corr->sNa - 1.001);
    s = bilinearInterp(corr->ionSigma, r, a, corr->sNr, corr->sNa, -0.9 * LARGEINT, 0.0);
    if (s <= 0.0)
    {
        return 0.0;
    }
    return ((double)s * rSLPixSize) * ((double)s * rSLPixSize);
}
