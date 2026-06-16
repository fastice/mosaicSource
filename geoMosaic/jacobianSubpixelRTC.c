/*
 * jacobianSubpixelRTC.c
 *
 * Sub-pixel RTC with per-sub-pixel Jacobian weighting.
 *
 * Motivation
 * ----------
 * When geomosaic reads an SLC directly (no pre-multilooking), the sub-pixel
 * accumulation IS the multilooking step.  In foreshortened terrain, multiple
 * map-space sub-pixels (xk, yl) project to the same or adjacent SLC range
 * bins.  Without Jacobian weighting those duplicate samples are counted
 * equally, inflating sumAb and biasing the area ratio.
 *
 * The Jacobian |J| = |∂r/∂x · ∂az/∂y − ∂r/∂y · ∂az/∂x| measures how much
 * unique radar-space area each map-space sub-pixel represents.  Weighting by
 * A_β · |J| makes the accumulation equivalent to GCOV's area-weighted SLC
 * multilooking in radar geometry:
 *   - in foreshortening: many sub-pixels share the same SLC bin; each has
 *     small |J|; their sum Σ|J| ≈ 1 unique radar-pixel worth of area.
 *   - in layover: J < 0 (range decreases as elevation increases up the slope).
 *
 * Layover treatment
 * -----------------
 * If the majority of valid sub-pixels have J < 0 the pixel is in a pure
 * layover zone: the terrain correction is unreliable.  The function then
 * returns TRUE (shadow-equivalent) so the caller suppresses the γ₀
 * correction.  σ₀ (imageTmp) is unaffected.  For mixed pixels (some J < 0,
 * some J > 0 — a foreshortening edge), the |J| weighting handles them
 * continuously and the function returns FALSE with valid accumulated values.
 *
 * Algorithm
 * ---------
 * Pass 1 — fill an nr × na grid of (r, az, h) using full Newton solves.
 *
 * Reference sign — compute J at the grid centre to determine which sign is
 *           "normal" for this orbit direction (ascending vs descending) and
 *           coordinate convention.  Layover is the OPPOSITE sign.
 *
 * Single AbAg call — compute Ab_c, Ag_c at the grid centre only.  All
 *           sub-pixels scale their areas from this single call:
 *             Ab_kl = Ab_c × |J_kl| / |J_c|
 *             Ag_kl = Ag_c × |n_kl| / |n_c|
 *           where n_kl = sqrt(1 + (∂h/∂x)² + (∂h/∂y)²) is the terrain
 *           normal magnitude estimated from the stored hGrid.  This
 *           eliminates O(nr × na) DEM lookups: four AbAg corner-DEM
 *           evaluations are replaced by finite differences in hGrid.
 *
 * Pass 2 — estimate J and n at each sub-pixel by central finite differences
 *           of the stored grids; accumulate with Ab_kl weight.
 *
 * ---------------------------------------------------------------------------
 * Integration steps (user performs these in the existing files):
 *
 * 1. geomosaic.h — add prototype:
 *      int32_t subPixelGammaRTCJacobian(double x, double y,
 *                                        inputImageStructure *inputImage,
 *                                        outputImageStructure *outputImage,
 *                                        void *dem,
 *                                        float *power, double *Ab, double *Ag,
 *                                        float *psiE);
 *
 * 2. geomosaic.c — add alongside the other global flags:
 *      int32_t jacobianSubPixelRTC = FALSE;
 *    Parse "-jacobianSubPixelRTC" to set both useSubPixelRTC and
 *    jacobianSubPixelRTC = TRUE.  Add to usage string.
 *
 * 3. makeGeoMosaic.c — extend the useSubPixelRTC dispatch block:
 *
 *      extern int32_t jacobianSubPixelRTC;
 *      extern int32_t linearSubPixelRTC;
 *      if ((jacobianSubPixelRTC
 *               ? subPixelGammaRTCJacobian(x, y, myImg, &outputImage, dem,
 *                                          &powerVal, &AbAcc, &AgAcc, &psiAcc)
 *           : linearSubPixelRTC
 *               ? subPixelGammaRTCLinear(x, y, myImg, &outputImage, dem,
 *                                        &powerVal, &AbAcc, &AgAcc, &psiAcc)
 *               : subPixelGammaRTC(x, y, myImg, &outputImage, dem,
 *                                  &powerVal, &AbAcc, &AgAcc, &psiAcc)))
 *          shadow = TRUE;
 *
 * 4. Makefile (geoMosaic/) — add jacobianSubpixelRTC.c to SRCS.
 *    Makefile (mosaicSource/) — add jacobianSubpixelRTC.o to GEOMOSAIC.
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
 * subPixelGammaRTCJacobian
 *
 * Drop-in replacement for subPixelGammaRTC that weights each sub-pixel's
 * contribution by Ab_kl = Ab_c × |J_kl|/|J_c|, derived from a single AbAg
 * call at the grid centre rather than one call per sub-pixel.
 *
 * Returns TRUE  — all sub-pixels shadowed, centre pixel shadowed, OR
 *                 majority of valid sub-pixels have J opposite to the centre
 *                 sign (pure layover; γ₀ correction suppressed).
 *         FALSE — at least one valid non-layover sub-pixel; *power, *Ab, *Ag,
 *                 *psiE are set.
 */
int32_t subPixelGammaRTCJacobian(double x, double y,
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
    float   val, psi1;
    double  sumPower, sumAb, sumAg, sumPsi, sumJabs;
    int32_t k, l, nr, na, nValid, nLayover, allShadow, refSign;

    /* Centre AbAg result and terrain-normal reference */
    double  Ab_c, Ag_c, shadow_c, absJ_c, n_c;

    double stepX = outputImage->deltaX * MTOKM;
    double stepY = outputImage->deltaY * MTOKM;

    double *rGrid, *azGrid, *hGrid;

    nr = (int)round(outputImage->deltaX / inputImage->rangePixelSize);
    na = (int)round(outputImage->deltaY / inputImage->azimuthPixelSize);
    nr = (nr < 1) ? 1 : nr;
    na = (na < 1) ? 1 : na;

    subOutput        = *outputImage;
    subOutput.deltaX = outputImage->deltaX / nr;
    subOutput.deltaY = outputImage->deltaY / na;

    dxkm = stepX / nr;
    dykm = stepY / na;

    rGrid  = (double *)malloc((size_t)(nr * na) * sizeof(double));
    azGrid = (double *)malloc((size_t)(nr * na) * sizeof(double));
    hGrid  = (double *)malloc((size_t)(nr * na) * sizeof(double));
    if (!rGrid || !azGrid || !hGrid) {
        free(rGrid); free(azGrid); free(hGrid);
        *power = (float)inputImage->noData;
        *Ab = 0.0; *Ag = 0.0; *psiE = 0.0;
        return TRUE;
    }

    /* ------------------------------------------------------------------
     * Pass 1: full Newton solve at every sub-pixel; store r, az, h.
     * ------------------------------------------------------------------ */
    for (k = 0; k < nr; k++) {
        xk = x + (k - 0.5 * (nr - 1)) * dxkm;
        for (l = 0; l < na; l++) {
            double lat, lon, h, hWGS, range, azimuth;
            yl = y + (l - 0.5 * (na - 1)) * dykm;
            xytoll1(xk, yl, HemiSphere, &lat, &lon, Rotation,
                    outputImage->slat);
            h    = getXYHeight(lat, lon, dem, inputImage->cpAll.Re, SPHERICAL);
            hWGS = sphericalToWGSElev(h, lat, inputImage->cpAll.Re);
            llToImageNew(lat, lon, hWGS, &range, &azimuth, inputImage);
            rGrid [k * na + l] = range;
            azGrid[k * na + l] = azimuth;
            hGrid [k * na + l] = h;
        }
    }

    /* ------------------------------------------------------------------
     * Reference sign: compute J at the grid centre to determine which
     * sign is "normal" for this orbit direction (ascending vs descending)
     * and coordinate convention.  Layover is the OPPOSITE sign.
     * J < 0 for descending right-looking in GrIMP (y northward);
     * J > 0 for ascending.  Using the centre avoids hardcoding either.
     *
     * absJ_c is retained for the per-sub-pixel Ab scaling below.
     * ------------------------------------------------------------------ */
    refSign = 0;
    absJ_c  = 0.0;
    {
        int kc = nr / 2, lc = na / 2;
        int kcp = (kc < nr-1) ? kc+1 : kc, kcm = (kc > 0) ? kc-1 : kc;
        int lcp = (lc < na-1) ? lc+1 : lc, lcm = (lc > 0) ? lc-1 : lc;
        double drc_x  = (kcp > kcm) ? (rGrid [kcp*na+lc] - rGrid [kcm*na+lc]) /
                                       ((kcp-kcm) * dxkm) : 0.0;
        double dazc_x = (kcp > kcm) ? (azGrid[kcp*na+lc] - azGrid[kcm*na+lc]) /
                                       ((kcp-kcm) * dxkm) : 0.0;
        double drc_y  = (lcp > lcm) ? (rGrid [kc*na+lcp] - rGrid [kc*na+lcm]) /
                                       ((lcp-lcm) * dykm) : 0.0;
        double dazc_y = (lcp > lcm) ? (azGrid[kc*na+lcp] - azGrid[kc*na+lcm]) /
                                       ((lcp-lcm) * dykm) : 0.0;
        double Jc = drc_x * dazc_y - drc_y * dazc_x;
        absJ_c = fabs(Jc);
        if      (Jc >  1e-6) refSign =  1;
        else if (Jc < -1e-6) refSign = -1;
        /* refSign == 0: degenerate grid (nr=na=1); layover detection skipped */
    }

    /* ------------------------------------------------------------------
     * Single AbAg call at the grid centre.  All sub-pixels scale their
     * projected areas from Ab_c and Ag_c rather than calling AbAg again.
     *   Ab_kl = Ab_c × |J_kl| / |J_c|
     *   Ag_kl = Ag_c × |n_kl| / |n_c|
     * where n_kl = sqrt(1 + (∂h/∂x)² + (∂h/∂y)²) is the terrain-normal
     * magnitude at sub-pixel (k,l), estimated from hGrid finite differences.
     * ------------------------------------------------------------------ */
    {
        int    kc  = nr / 2,  lc  = na / 2;
        int    kcp = (kc < nr-1) ? kc+1 : kc, kcm = (kc > 0) ? kc-1 : kc;
        int    lcp = (lc < na-1) ? lc+1 : lc, lcm = (lc > 0) ? lc-1 : lc;
        double xk_c = x + (kc - 0.5 * (nr - 1)) * dxkm;
        double yl_c = y + (lc - 0.5 * (na - 1)) * dykm;
        double dhx  = (kcp > kcm) ? (hGrid[kcp*na+lc] - hGrid[kcm*na+lc]) /
                                     ((kcp-kcm) * dxkm) : 0.0;
        double dhy  = (lcp > lcm) ? (hGrid[kc*na+lcp] - hGrid[kc*na+lcm]) /
                                     ((lcp-lcm) * dykm) : 0.0;
        n_c      = sqrt(1.0 + dhx*dhx + dhy*dhy);
        Ab_c     = 0.0;
        Ag_c     = 0.0;
        shadow_c = 1.0;
        AbAg(xk_c, yl_c, azGrid[kc*na+lc], inputImage, &subOutput, dem,
             &Ab_c, &Ag_c, &shadow_c, FALSE);
    }

    if (shadow_c < -0.001 || Ab_c <= 0.0 || Ag_c <= 0.0) {
        free(rGrid); free(azGrid); free(hGrid);
        *power = (float)inputImage->noData;
        *Ab = 0.0; *Ag = 0.0; *psiE = 0.0;
        return TRUE;
    }

    /* ------------------------------------------------------------------
     * Pass 2: accumulate with Ab_kl = Ab_c × |J_kl|/|J_c| weighting.
     * ------------------------------------------------------------------ */
    sumPower  = 0.0;
    sumAb     = 0.0;
    sumAg     = 0.0;
    sumPsi    = 0.0;
    sumJabs   = 0.0;
    nValid    = 0;
    nLayover  = 0;
    allShadow = TRUE;

    for (k = 0; k < nr; k++) {
        xk = x + (k - 0.5 * (nr - 1)) * dxkm;

        for (l = 0; l < na; l++) {
            double range, azimuth, h0;
            double drx, dry, dazx, dazy, dh_x, dh_y;
            double J, absJ, n_kl, Ab1, Ag1, wt;
            int kp, km, lp, lm;

            yl = y + (l - 0.5 * (na - 1)) * dykm;

            range   = rGrid [k * na + l];
            azimuth = azGrid[k * na + l];
            h0      = hGrid [k * na + l];

            /* Jacobian and terrain normal via central differences. */
            kp = (k < nr - 1) ? k + 1 : k;
            km = (k > 0)      ? k - 1 : k;
            lp = (l < na - 1) ? l + 1 : l;
            lm = (l > 0)      ? l - 1 : l;

            drx  = (kp > km) ? (rGrid [kp*na+l] - rGrid [km*na+l]) /
                                ((kp - km) * dxkm) : 0.0;
            dazx = (kp > km) ? (azGrid[kp*na+l] - azGrid[km*na+l]) /
                                ((kp - km) * dxkm) : 0.0;
            dry  = (lp > lm) ? (rGrid [k*na+lp] - rGrid [k*na+lm]) /
                                ((lp - lm) * dykm) : 0.0;
            dazy = (lp > lm) ? (azGrid[k*na+lp] - azGrid[k*na+lm]) /
                                ((lp - lm) * dykm) : 0.0;

            dh_x = (kp > km) ? (hGrid[kp*na+l] - hGrid[km*na+l]) /
                                ((kp - km) * dxkm) : 0.0;
            dh_y = (lp > lm) ? (hGrid[k*na+lp] - hGrid[k*na+lm]) /
                                ((lp - lm) * dykm) : 0.0;

            J     = drx * dazy - dry * dazx;
            absJ  = fabs(J);
            n_kl  = sqrt(1.0 + dh_x*dh_x + dh_y*dh_y);

            val = bilinearInterp((float **)inputImage->image,
                                 range, azimuth,
                                 inputImage->rangeSize, inputImage->azimuthSize,
                                 (float)inputImage->noData,
                                 (float)inputImage->noData);
            if (val <= 0)
                continue;

            psi1 = applyCorrections(&val, inputImage, range, azimuth, h0);
            if (val <= 0)
                continue;

            /* Scale centre areas to this sub-pixel; no AbAg call needed. */
            Ab1 = (absJ_c > 1e-10) ? Ab_c * absJ / absJ_c : Ab_c;
            Ag1 = (n_c    > 1e-10) ? Ag_c * n_kl / n_c    : Ag_c;
            wt  = Ab1;

            allShadow  = FALSE;
            sumPower  += (double)val * wt;
            sumAb     += wt;
            sumAg     += Ag1;
            sumPsi    += (double)psi1 * absJ;
            sumJabs   += absJ;
            nValid++;
            /* Layover: J has opposite sign to the reference (centre pixel).
             * refSign==0 means the grid is degenerate; skip layover counting. */
            if (refSign != 0 && ((J > 0) != (refSign > 0)))
                nLayover++;
        }
    }

    free(rGrid);
    free(azGrid);
    free(hGrid);

    /* Pure layover: majority of valid sub-pixels have J with opposite sign
     * to the reference (centre pixel).
     * Return TRUE so the caller suppresses the γ₀ correction (shadow
     * treatment).  σ₀ in imageTmp is unaffected by the caller's shadow flag. */
    if (nValid > 0 && (nLayover > nValid / 2 || (maskLayover && nLayover > 0))) {
        *power = (float)(sumPower / sumAb);
        *Ab    = 0.0;
        *Ag    = 0.0;
        *psiE  = (float)(sumJabs > 0.0 ? sumPsi / sumJabs : 0.0);
        return TRUE;
    }

    if (nValid == 0 || sumAb <= 0.0) {
        *power = (float)inputImage->noData;
        *Ab    = 0.0;
        *Ag    = 0.0;
        *psiE  = 0.0;
        return allShadow;
    }

    *power = (float)(sumPower / sumAb);
    *Ab    = sumAb;
    *Ag    = sumAg;
    *psiE  = (float)(sumJabs > 0.0 ? sumPsi / sumJabs : 0.0);

    return FALSE;
}
