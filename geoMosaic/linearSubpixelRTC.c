/*
 * linearSubpixelRTC.c
 *
 * Fast sub-pixel RTC via first-order Jacobian approximation.
 *
 * subPixelGammaRTC calls llToImageNew (a Newton solver) for every sub-pixel —
 * nr*na times per output pixel.  subPixelGammaRTCLinear replaces that with
 * three anchor solves (center, center+Δx, center+Δy) and a bilinear Taylor
 * expansion, so the inner loop needs only bilinearInterp + applyCorrections +
 * AbAg.  AbAg itself still calls into the DEM, but the dominant Newton cost
 * is gone.
 *
 * Accuracy: the Jacobian is exact to first order in sub-pixel displacement.
 * At typical NISAR output spacings (200–500 m) the range/azimuth mapping is
 * nearly linear within one pixel, so errors are well below a pixel in the
 * image domain.
 *
 * The finite-difference step uses the full output pixel size (not the
 * sub-pixel size) for a better-conditioned gradient estimate.
 *
 * ---------------------------------------------------------------------------
 * Integration steps (user performs these in the existing files):
 *
 * 1. geomosaic.h — add prototype:
 *      int32_t subPixelGammaRTCLinear(double x, double y,
 *                                     inputImageStructure *inputImage,
 *                                     outputImageStructure *outputImage,
 *                                     void *dem,
 *                                     float *power, double *Ab, double *Ag,
 *                                     float *psiE);
 *
 * 2. makeGeoMosaic.c — in the useSubPixelRTC block, replace:
 *
 *      if (subPixelGammaRTC(x, y, myImg, &outputImage, dem,
 *                           &powerVal, &AbAcc, &AgAcc, &psiAcc))
 *          shadow = TRUE;
 *
 *    with:
 *
 *      extern int32_t linearSubPixelRTC;
 *      if ((linearSubPixelRTC
 *               ? subPixelGammaRTCLinear(x, y, myImg, &outputImage, dem,
 *                                        &powerVal, &AbAcc, &AgAcc, &psiAcc)
 *               : subPixelGammaRTC(x, y, myImg, &outputImage, dem,
 *                                  &powerVal, &AbAcc, &AgAcc, &psiAcc)))
 *          shadow = TRUE;
 *
 * 3. Makefile — add linearSubpixelRTC.o to the SRCS list for geomosaic.
 * ---------------------------------------------------------------------------
 */

#include <math.h>
#include <stdlib.h>
#include "mosaicSource/common/common.h"
#include "geomosaic.h"

extern void  AbAg(double x, double y, double azimuth,
                  inputImageStructure *inputImage,
                  outputImageStructure *outputImage,
                  void *dem, double *Ab, double *Ag,
                  double *shadow, int32_t recycle);
extern float applyCorrections(float *value, inputImageStructure *inputImage,
                              double range, double azimuth, double h);

/*
 * anchor — run the full coordinate chain for one map point and return
 * range and azimuth.  Height is returned so the caller can reuse it
 * for applyCorrections.
 */
static void anchor(double xkm, double ykm,
                   int32_t hemi, double rotation, double slat,
                   inputImageStructure *inputImage, void *dem,
                   double *range, double *azimuth, double *h_out)
{
    double lat, lon, h, hWGS;
    xytoll1(xkm, ykm, hemi, &lat, &lon, rotation, slat);
    h    = getXYHeight(lat, lon, dem, inputImage->cpAll.Re, SPHERICAL);
    hWGS = sphericalToWGSElev(h, lat, inputImage->cpAll.Re);
    llToImageNew(lat, lon, hWGS, range, azimuth, inputImage);
    *h_out = h;
}

/*
 * subPixelGammaRTCLinear
 *
 * Drop-in replacement for subPixelGammaRTC that uses three anchor solves
 * and a first-order Taylor expansion to approximate range/azimuth at every
 * sub-pixel centre, skipping the Newton solver in the inner loop.
 *
 * Return value and output arguments are identical to subPixelGammaRTC.
 */
int32_t subPixelGammaRTCLinear(double x, double y,
                                inputImageStructure *inputImage,
                                outputImageStructure *outputImage,
                                void *dem,
                                float *power, double *Ab, double *Ag,
                                float *psiE)
{
    extern int32_t HemiSphere;
    extern double  Rotation;
    extern int32_t maskLayover;

    outputImageStructure subOutput;
    double  dxkm, dykm;
    double  xk, yl;
    double  Ab1, Ag1, shadow;
    float   val, psi1;
    double  sumPower, sumAb, sumAg, sumPsi;
    int32_t k, l, nr, na, nValid, allShadow, abAgRecycle;

    /* Anchor results */
    double r0,   az0,  h0;   /* center             */
    double r_px, az_px, h_px; /* center + full Δx  */
    double r_py, az_py, h_py; /* center + full Δy  */

    /* Jacobian coefficients */
    double dRdX, dRdY, dAzdX, dAzdY;

    /* Full output pixel size in km — used as the finite-difference step. */
    double stepX = outputImage->deltaX * MTOKM;
    double stepY = outputImage->deltaY * MTOKM;

    nr = (int)round(outputImage->deltaX / inputImage->rangePixelSize);
    na = (int)round(outputImage->deltaY / inputImage->azimuthPixelSize);
    nr = (nr < 1) ? 1 : nr;
    na = (na < 1) ? 1 : na;

    subOutput        = *outputImage;
    subOutput.deltaX = outputImage->deltaX / nr;
    subOutput.deltaY = outputImage->deltaY / na;

    dxkm = stepX / nr;
    dykm = stepY / na;

    /* Three anchor solves. */
    anchor(x,        y,        HemiSphere, Rotation, outputImage->slat,
           inputImage, dem, &r0,   &az0,   &h0);
    anchor(x + stepX, y,        HemiSphere, Rotation, outputImage->slat,
           inputImage, dem, &r_px, &az_px, &h_px);
    anchor(x,        y + stepY, HemiSphere, Rotation, outputImage->slat,
           inputImage, dem, &r_py, &az_py, &h_py);

    dRdX  = (r_px  - r0)  / stepX;
    dRdY  = (r_py  - r0)  / stepY;
    dAzdX = (az_px - az0) / stepX;
    dAzdY = (az_py - az0) / stepY;

    int32_t isLayover = FALSE;
    if (maskLayover) {
        double Jpix = dRdX * dAzdY - dRdY * dAzdX;
        if (fabs(Jpix) > 1e-6 &&
            (Jpix > 0) != (inputImage->passType == ASCENDING))
            isLayover = TRUE;
    }

    sumPower  = 0.0;
    sumAb     = 0.0;
    sumAg     = 0.0;
    sumPsi    = 0.0;
    nValid      = 0;
    allShadow   = TRUE;
    abAgRecycle = FALSE;

    for (k = 0; k < nr; k++) {
        double dxk = (k - 0.5 * (nr - 1)) * dxkm;
        xk = x + dxk;

        for (l = 0; l < na; l++) {
            double dyl = (l - 0.5 * (na - 1)) * dykm;
            double range, azimuth;
            yl = y + dyl;

            /* First-order Taylor expansion — no Newton solve needed. */
            range   = r0   + dRdX * dxk + dRdY * dyl;
            azimuth = az0  + dAzdX * dxk + dAzdY * dyl;

            val = bilinearInterp((float **)inputImage->image,
                                 range, azimuth,
                                 inputImage->rangeSize, inputImage->azimuthSize,
                                 (float)inputImage->noData,
                                 (float)inputImage->noData);
            if (val <= 0)
                continue;

            /* h0 (center height) is used for all sub-pixels; the variation
             * across one output pixel is negligible for the depression angle. */
            psi1 = applyCorrections(&val, inputImage, range, azimuth, h0);
            if (val <= 0)
                continue;

            shadow = 1.0;
            AbAg(xk, yl, azimuth, inputImage, &subOutput, dem,
                 &Ab1, &Ag1, &shadow, abAgRecycle);
            abAgRecycle = TRUE;

            if (shadow < -0.001)
                continue;

            allShadow  = FALSE;
            sumPower  += (double)val * Ab1;
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

    *power = (float)(sumPower / sumAb);
    *Ab    = isLayover ? 0.0 : sumAb;
    *Ag    = isLayover ? 0.0 : sumAg;
    *psiE  = (float)(sumPsi / nValid);

    return isLayover;
}
