/*
 * subpixelRTC.c
 *
 * Sub-pixel RTC: samples the input image at a nRangeLooks x nAzimuthLooks
 * grid within the output pixel footprint, accumulating power and area from
 * the same set of sub-pixels so the two are consistent with each other.
 *
 * In the current approach interpolatePowerInputImage takes one sample at the
 * output pixel centre while areaAboutXY searches neighbours for area
 * contributions — the value and the correction describe different things.
 * Here each sub-pixel contributes both its power and its Ab/Ag areas, so a
 * layover zone that spans the output pixel is handled without a neighbour
 * search.
 *
 * ---------------------------------------------------------------------------
 * Integration steps (user performs these in the existing files):
 *
 * 1. makeGeoMosaic.c — remove 'static' from AbAg() and applyCorrections().
 *
 * 2. geomosaic.h — add prototypes:
 *      void  AbAg(double x, double y, double azimuth,
 *                 inputImageStructure *inputImage,
 *                 outputImageStructure *outputImage,
 *                 void *dem, double *Ab, double *Ag,
 *                 double *shadow, int32_t recycle);
 *      float applyCorrections(float *value, inputImageStructure *inputImage,
 *                             double range, double azimuth, double h);
 *      int32_t subPixelGammaRTC(double x, double y,
 *                               inputImageStructure *inputImage,
 *                               outputImageStructure *outputImage,
 *                               void *dem,
 *                               float *power, double *Ab, double *Ag,
 *                               float *psiE);
 *
 * 3. geomosaic.c — add alongside the other global flags:
 *      int32_t useSubPixelRTC = FALSE;
 *    Add to readArgs() and its signature, parsed with e.g. "-subPixelRTC".
 *    Add extern declaration in main():
 *      extern int32_t useSubPixelRTC;
 *    Pass it through to makeGeoMosaic() or leave as global (matches pattern
 *    of S1Cal, noPower, etc.).
 *
 * 4. makeGeoMosaic.c — inside the i1/j1 loop, replace:
 *
 *      value = interpolatePowerInputImage(inputImage[i], range, azimuth);
 *      psiE  = applyCorrections(&value, &(inputImage[i]), range, azimuth, h)
 *              * invNAvg;
 *      if ((S1Cal & TRUE) == TRUE) {
 *          if (areaAboutXY(range, azimuth, x, y, &(inputImage[i]),
 *                          &outputImage, dem, value, &test,
 *                          recycle, &Ab, &Ag) < -0.001)
 *              shadow = TRUE;
 *          AbCum   = Ab;
 *          AgCum   = Ag;
 *          recycle = TRUE;
 *      }
 *
 *    with:
 *
 *      extern int32_t useSubPixelRTC;
 *      if (useSubPixelRTC && (S1Cal & TRUE) == TRUE) {
 *          float   powerVal;
 *          double  AbAcc, AgAcc;
 *          float   psiAcc;
 *          if (subPixelGammaRTC(x, y, &(inputImage[i]), &outputImage, dem,
 *                               &powerVal, &AbAcc, &AgAcc, &psiAcc))
 *              shadow = TRUE;
 *          value = powerVal;
 *          psiE  = psiAcc;
 *          AbCum = AbAcc;
 *          AgCum = AgAcc;
 *      } else {
 *          value = interpolatePowerInputImage(inputImage[i], range, azimuth);
 *          psiE  = applyCorrections(&value, &(inputImage[i]),
 *                                   range, azimuth, h) * invNAvg;
 *          if ((S1Cal & TRUE) == TRUE) {
 *              if (areaAboutXY(range, azimuth, x, y, &(inputImage[i]),
 *                              &outputImage, dem, value, &test,
 *                              recycle, &Ab, &Ag) < -0.001)
 *                  shadow = TRUE;
 *              AbCum   = Ab;
 *              AgCum   = Ag;
 *              recycle = TRUE;
 *          }
 *      }
 *
 * 5. Makefile — add subpixelRTC.o to the object list for geomosaic.
 * ---------------------------------------------------------------------------
 */

#include <math.h>
#include <stdlib.h>
#include "mosaicSource/common/common.h"
#include "geomosaic.h"

/* Must be made non-static in makeGeoMosaic.c — see step 1 above */
extern void  AbAg(double x, double y, double azimuth,
                  inputImageStructure *inputImage,
                  outputImageStructure *outputImage,
                  void *dem, double *Ab, double *Ag,
                  double *shadow, int32_t recycle);
extern float applyCorrections(float *value, inputImageStructure *inputImage,
                              double range, double azimuth, double h);

/*
 * subPixelGammaRTC
 *
 * For the output pixel centred at (x, y) in polar-stereographic km coords,
 * samples the input image on a nRangeLooks x nAzimuthLooks sub-pixel grid
 * covering the full output pixel footprint.  Each sub-pixel contributes its
 * calibrated power and its projected areas Ab (beta plane) and Ag (gamma
 * plane) to running totals.
 *
 * Returns TRUE  — every sub-pixel was in shadow; *power = noData
 *         FALSE — at least one valid sub-pixel found
 *
 * *power  Ab-weighted mean calibrated power (same units / scale as the
 *         value produced by interpolatePowerInputImage + applyCorrections).
 *         Ab-weighting is consistent with the beta-nought reference used by
 *         applyCorrections for S1Cal.
 * *Ab     sum of beta-plane areas over valid sub-pixels
 * *Ag     sum of gamma-plane areas over valid sub-pixels
 * *psiE   mean depression angle in degrees over valid sub-pixels
 */
int32_t subPixelGammaRTC(double x, double y,
                          inputImageStructure *inputImage,
                          outputImageStructure *outputImage,
                          void *dem,
                          float *power, double *Ab, double *Ag, float *psiE)
{
    extern int32_t HemiSphere;
    extern double  Rotation;
    extern int32_t maskLayover;

    outputImageStructure subOutput; /* scaled copy used for AbAg corner offsets */
    double  dxkm, dykm;            /* sub-pixel spacing in km                  */
    double  xk, yl;                /* sub-pixel centre position                */
    double  lat, lon, h, hWGS;
    double  range, azimuth;
    double  Ab1, Ag1, shadow;
    float   val, psi1;
    double  sumPower, sumAb, sumAg, sumPsi;
    int32_t k, l, nr, na, nValid, allShadow, isLayover;

    /*
     * Sub-pixel grid size = output pixel size / input pixel size.
     * rangePixelSize and azimuthPixelSize are already the multilooked slant-
     * range spacings (nLooks * single-look spacing, set in parseInputFile).
     * deltaX/Y are map-projection metres.  The sin(theta) ground-range
     * correction for range is omitted — it requires a DEM lookup — so nr may
     * be slightly over-estimated at steep incidence, which is conservative.
     */
    nr = (int)round(outputImage->deltaX / inputImage->rangePixelSize);
    na = (int)round(outputImage->deltaY / inputImage->azimuthPixelSize);
    nr = (nr < 1) ? 1 : nr;
    na = (na < 1) ? 1 : na;

    /*
     * AbAg computes pixel corner offsets as +/- 0.5 * outputImage->deltaX/Y.
     * Give it a sub-pixel-sized output structure so those corners span only
     * the sub-pixel area (1/nr x 1/na of the full output pixel).
     */
    subOutput        = *outputImage;
    subOutput.deltaX = outputImage->deltaX / nr;
    subOutput.deltaY = outputImage->deltaY / na;

    dxkm = outputImage->deltaX * MTOKM / nr;
    dykm = outputImage->deltaY * MTOKM / na;

    sumPower  = 0.0;
    sumAb     = 0.0;
    sumAg     = 0.0;
    sumPsi    = 0.0;
    nValid    = 0;
    allShadow = TRUE;
    isLayover = FALSE;

    /* Pixel-level layover check: 3 anchor solves give the Jacobian of the
     * map→SLC transform.  If its sign is reversed relative to the expected
     * sign for this orbit direction the whole output pixel is in layover. */
    if (maskLayover) {
        double stepXkm = outputImage->deltaX * MTOKM;
        double stepYkm = outputImage->deltaY * MTOKM;
        double r0, az0, r_px, az_px, r_py, az_py;
        double dRdX, dRdY, dAzdX, dAzdY, Jpix;

        xyToLLProj(x, y, &lat, &lon, &(outputImage->proj));
        h    = getXYHeight(lat, lon, dem, inputImage->cpAll.Re, SPHERICAL);
        hWGS = sphericalToWGSElev(h, lat, inputImage->cpAll.Re);
        llToImageNew(lat, lon, hWGS, &r0, &az0, inputImage);

        xyToLLProj(x + stepXkm, y, &lat, &lon, &(outputImage->proj));
        h    = getXYHeight(lat, lon, dem, inputImage->cpAll.Re, SPHERICAL);
        hWGS = sphericalToWGSElev(h, lat, inputImage->cpAll.Re);
        llToImageNew(lat, lon, hWGS, &r_px, &az_px, inputImage);

        xyToLLProj(x, y + stepYkm, &lat, &lon, &(outputImage->proj));
        h    = getXYHeight(lat, lon, dem, inputImage->cpAll.Re, SPHERICAL);
        hWGS = sphericalToWGSElev(h, lat, inputImage->cpAll.Re);
        llToImageNew(lat, lon, hWGS, &r_py, &az_py, inputImage);

        dRdX  = (r_px  - r0)  / stepXkm;
        dRdY  = (r_py  - r0)  / stepYkm;
        dAzdX = (az_px - az0) / stepXkm;
        dAzdY = (az_py - az0) / stepYkm;
        Jpix  = dRdX * dAzdY - dRdY * dAzdX;

        if (fabs(Jpix) > 1e-6 &&
            (Jpix > 0) != (inputImage->passType == ASCENDING))
            isLayover = TRUE;
    }

    for (k = 0; k < nr; k++) {
        xk = x + (k - 0.5 * (nr - 1)) * dxkm;

        for (l = 0; l < na; l++) {
            yl = y + (l - 0.5 * (na - 1)) * dykm;

            /* Sub-pixel centre -> lat/lon/height */
            xyToLLProj(xk, yl, &lat, &lon, &(outputImage->proj));
            h    = getXYHeight(lat, lon, dem, inputImage->cpAll.Re, SPHERICAL);
            hWGS = sphericalToWGSElev(h, lat, inputImage->cpAll.Re);

            /* lat/lon/height -> SAR image coordinates */
            llToImageNew(lat, lon, hWGS, &range, &azimuth, inputImage);

            /* Sample the input image */
            val = bilinearInterp((float **)inputImage->image,
                                 range, azimuth,
                                 inputImage->rangeSize, inputImage->azimuthSize,
                                 (float)inputImage->noData,
                                 (float)inputImage->noData);
            if (val <= 0)
                continue;

            /* Apply calibration / antenna-pattern correction */
            psi1 = applyCorrections(&val, inputImage, range, azimuth, h);
            if (val <= 0)
                continue;

            /* Compute projected areas; shadow > 0 enables the shadow check */
            shadow = 1.0;
            AbAg(xk, yl, azimuth, inputImage, &subOutput, dem,
                 &Ab1, &Ag1, &shadow, FALSE);

            if (shadow < -0.001)
                continue;   /* sub-pixel in shadow */

            allShadow  = FALSE;
            sumPower  += (double)val * Ab1;  /* Ab-weighted accumulation */
            sumAb     += Ab1;
            sumAg     += Ag1;
            sumPsi    += psi1;
            nValid++;
        }
    }

    if (nValid == 0 || sumAb <= 0.0) {
        *power = (float)inputImage->noData;
        *Ab    = 0.0;
        *Ag    = 0.0;
        *psiE  = 0.0;
        return allShadow;
    }

    *power = (float)(sumPower / sumAb);   /* Ab-weighted mean calibrated power */
    *Ab    = isLayover ? 0.0 : sumAb;
    *Ag    = isLayover ? 0.0 : sumAg;
    *psiE  = (float)(sumPsi / nValid);   /* mean depression angle, degrees     */

    return isLayover;
}
