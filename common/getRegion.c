#include "common.h"
#include <unistd.h>

/* Set by make3DOffsets.c (and similarly-structured callers) while inside the
   crossing-pair inner loop, so this output visually groups with the rest of
   that loop's tab-indented prints. Defaults to FALSE for all other callers. */
int32_t indentRegionOutput = FALSE;

/* Ignore any embedded VRT dataset mask band on offset inputs (readOffsets.c);
   default off, meaning masks are honored by default when present. Defined
   once here (not per-program) since common/readOffsets.c is linked into
   every program in mosaicSource/ -- each program's own CLI parser sets this
   via its own -noMask flag. */
int32_t noMask = FALSE;
/*  Skip the azimuth-offset raster read entirely.  Set by mosaic3d ONLY when -noAzimuthRows is
    given AND the hopper is running, i.e. when nothing downstream can consume offsets->da.  Lives
    here rather than in mosaic3d.c because readOffsets.c is linked into every binary, and an
    extern to a mosaic3d-only symbol would break the link for rparams/azparams/tiepoints. */
int32_t skipAzimuthOffsets = FALSE;
/* Range-offset accuracy term (see rangeAccuracyVar, common/interpOffsets.c).
   Defined HERE rather than in mosaic3d.c, for the same reason noMask is: every
   program that links interpOffsets.c/readOffsets.c needs the symbol, and a
   definition in the mosaic3d main object would leave rparams/azparams with an
   undefined reference (and produce a DT_TEXTREL in the PIE link). */
int32_t rSigmaResidual = TRUE;
/* Azimuth-offset accuracy term (see azimuthAccuracyVar, common/interpOffsets.c).
   The exact analogue of rSigmaResidual above, and NEW: offsets.sigmaAresidual has
   always been parsed but only ever used as a skip/threshold gate, never in the
   error budget -- so the azimuth sigma carried only local matching noise plus the
   azparams baseline covariance, with nothing standing in for long-wavelength
   error.  Consumed by speckleTrackMosaicJoint.c; -noASigmaResidual disables. */
int32_t aSigmaResidual = TRUE;
/* -pairCountLegacy: use the pre-2026-08-26 over-count formula, f = (nOuter +
   nPairs)/2, which treated the PAIR count as if it were the distinct count of
   inner images.  Retained as a control so the counting change can be separated
   from the compounding fix (which is not optional and applies either way). */
int32_t pairCountLegacy = FALSE;
double rSigmaConst = 0.0;
/* Correlation parameter used by inflatePairOverCount() (see
   common/scalingFunctions.c for the model).  BOTH default to 0.6.

   WHAT IT IS.  rho is the fraction of the per-pixel variance that is COMMON to
   every crossing pair and therefore never averages down.  It stands in for
   correlation the cheap over-count factor cannot represent -- principally
   between pairs that SHARE AN IMAGE, since a frame's baseline residual attaches
   to that frame and every pair it joins inherits it.  The exact treatment is
   f = sum(d_i^2)/(2P) over contributing images, but that needs a per-pixel record
   of WHICH images contributed; rho is the fixed-cost stand-in.

   WHY 0.6.  Raising the crossing time threshold adds PAIRS while the contributing
   IMAGES stay fixed (image count moves 0.5% from T=12 to T=37 while pair count
   moves 2.5x), so it adds no information and the measured error is flat -- and is,
   to a few percent across T = 12/37/10000.  A correctly scaled formal error must
   be flat too.  Measured on 406845 px of stable ground against a Sentinel-1
   reference, over 12 product-threshold cells RUN AT rho = 0.6 ITSELF (not
   interpolated), with k = std(difference)/RMS(e), for T = 6/12/37/10000:

       rho=0    k = 1.3 - 3.1, drifting ~2x with threshold   -- excluded
       rho=0.6  phase   0.97 0.99 0.99 0.98
                both    1.01 1.05 1.06 1.06
                offsets 0.77 0.80 0.82 0.83
       rho=1    k = 0.66 - 0.88, overshoots                  -- excluded

   0.6 makes every cell threshold-flat, which is what rho can legitimately fix.
   The residual offsets level (~0.80, flat) is NOT a threshold effect and no
   single rho lifts it without breaking the other two; offsets-only is not a
   shipped NISAR product.

   INDEPENDENT CONFIRMATION.  A sandbox build that consumes each image at most
   once per pixel -- a matching, which has NO over-counting by construction and a
   completely different pairing graph -- calibrates to k = 1 at rho = 0.584
   (phase) and 0.501 (offsets), using the matching's own factor f = rho*m +
   (1-rho).  Two different criteria on two different graphs both land near 0.6.
   Note this places 0.6 at the TOP of the supported 0.50-0.87 range, not its
   centre.  (An earlier version of that experiment ran the matching at f = 1 and
   reported k = 2.4; f = 1 is the rho = 0 case, so it measured the assumption.)

   WHAT IT IS NOT FOR.  rho sets the total variance; it does not fix the spatial
   distribution of the error, which is wrong in a way no scalar addresses -- see
   Documents/mosaic3d.md "Interpreting the formal errors".  Keep the two values
   equal: unequal ones drive f_phase/f_offsets past the >2 + S_p/S_o threshold at
   which the combined error can exceed one of its own inputs (measured at 2.49%
   of pixels with 0.5/0.0).

   Overridable with -rhoPhase / -rhoOffsets. */
double rhoPhase = 0.6;
double rhoOffsets = 0.6;

/*
  Clip a convex polygon against a single half-plane  nx*x + ny*y >= d.
  Input polygon has nIn vertices; output written to xOut/yOut.
  Returns number of output vertices (0 = entirely outside).
*/
int clipByHalfPlane(double *xIn, double *yIn, int nIn,
                     double *xOut, double *yOut,
                     double nx, double ny, double d)
{
    int nOut = 0;
    for (int i = 0; i < nIn; i++)
    {
        int j = (i + 1) % nIn;
        double di = nx * xIn[i] + ny * yIn[i] - d;
        double dj = nx * xIn[j] + ny * yIn[j] - d;
        if (di >= 0.0)
        {
            xOut[nOut] = xIn[i];
            yOut[nOut] = yIn[i];
            nOut++;
        }
        if ((di >= 0.0) != (dj >= 0.0))
        {
            /* Edge crosses boundary — compute and emit intersection point */
            double t = di / (di - dj);
            xOut[nOut] = xIn[i] + t * (xIn[j] - xIn[i]);
            yOut[nOut] = yIn[i] + t * (yIn[j] - yIn[i]);
            nOut++;
        }
    }
    return nOut;
}

/*
  Project the swath's 4 corner control points into polar-stereographic XY (metres),
  angle-sort them into a simple polygon, then clip against the (padded) output image
  rectangle. Shared by getRegion() (which reduces the result to an axis-aligned bbox)
  and getRowBounds() (which needs the polygon itself for per-row intersection).

  xOut/yOut must be sized for at least 8 vertices (4-gon clipped against a rectangle's
  4 half-planes can produce at most 8). Returns vertex count (0 = no overlap).
*/
static int clipSwathToOutputRect(inputImageStructure *image, outputImageStructure *outputImage,
                                  double *xOut, double *yOut)
{
    extern double Rotation;
    double xSwath[4], ySwath[4];   /* swath corners in metres */
    double xA[8], yA[8];           /* work buffers for polygon clipping */
    double xB[8], yB[8];
    double xa1, ya1;
    double xImgMin, xImgMax, yImgMin, yImgMax;
    double pad;
    int n, i;

    /* --- 1. Get swath corner points in metres -------------------------------- */
    /* Control points: [1]=ll, [2]=lr, [3]=ul, [4]=ur.
       Use order ll→lr→ur→ul (indices 1,2,4,3) to form a proper convex CW polygon. */
    {
        int cpIdx[4] = {1, 2, 4, 3};
        for (i = 0; i < 4; i++)
        {
            lltoxy1(image->latControlPoints[cpIdx[i]], image->lonControlPoints[cpIdx[i]],
                    &xa1, &ya1, Rotation, outputImage->slat);
            xSwath[i] = xa1 * KMTOM;
            ySwath[i] = ya1 * KMTOM;
        }
    }
    /* Sort corners CCW by angle from centroid — same fix as getIntersect().
       The fixed {1,2,4,3} index order produces a self-intersecting quadrilateral
       for polar orbits spanning large longitude ranges, which breaks S-H clipping. */
    {
        double cx = 0.0, cy = 0.0, ang[4], tx, ty;
        int ii, jj, imin;
        for (i = 0; i < 4; i++) { cx += xSwath[i]; cy += ySwath[i]; }
        cx *= 0.25; cy *= 0.25;
        for (ii = 0; ii < 4; ii++) ang[ii] = atan2(ySwath[ii] - cy, xSwath[ii] - cx);
        for (ii = 0; ii < 3; ii++) {
            imin = ii;
            for (jj = ii + 1; jj < 4; jj++) if (ang[jj] < ang[imin]) imin = jj;
            if (imin != ii) {
                tx = xSwath[ii]; xSwath[ii] = xSwath[imin]; xSwath[imin] = tx;
                ty = ySwath[ii]; ySwath[ii] = ySwath[imin]; ySwath[imin] = ty;
                tx = ang[ii];    ang[ii]    = ang[imin];    ang[imin]    = tx;
            }
        }
    }

    /* --- 2. Output image rectangle in metres --------------------------------- */
    /* Pad scales with along-track length: at 6000 km the pad is 40 km, minimum 15 km. */
    {
        double midBotX = (xSwath[0] + xSwath[1]) * 0.5;
        double midBotY = (ySwath[0] + ySwath[1]) * 0.5;
        double midTopX = (xSwath[2] + xSwath[3]) * 0.5;
        double midTopY = (ySwath[2] + ySwath[3]) * 0.5;
        double dx = midTopX - midBotX, dy = midTopY - midBotY;
        double trackLength = sqrt(dx * dx + dy * dy);
        pad = max(15000.0, trackLength * (40e3 / 6000e3));
    }
    xImgMin = outputImage->originX - pad;
    xImgMax = outputImage->originX + outputImage->xSize * outputImage->deltaX + pad;
    yImgMin = outputImage->originY - pad;
    yImgMax = outputImage->originY + outputImage->ySize * outputImage->deltaY + pad;

    /* --- 3. Clip swath polygon against all four sides of the output rectangle -
             Sutherland-Hodgman: clip against each half-plane in turn.
             Half-planes: x >= xImgMin, x <= xImgMax,
                          y >= yImgMin, y <= yImgMax               */
    for (i = 0; i < 4; i++) { xA[i] = xSwath[i]; yA[i] = ySwath[i]; }
    n = 4;

    n = clipByHalfPlane(xA, yA, n, xB, yB,  1.0,  0.0, xImgMin);  /* x >= xImgMin */
    if (n == 0) return 0;
    n = clipByHalfPlane(xB, yB, n, xA, yA, -1.0,  0.0, -xImgMax); /* x <= xImgMax */
    if (n == 0) return 0;
    n = clipByHalfPlane(xA, yA, n, xB, yB,  0.0,  1.0, yImgMin);  /* y >= yImgMin */
    if (n == 0) return 0;
    n = clipByHalfPlane(xB, yB, n, xA, yA,  0.0, -1.0, -yImgMax); /* y <= yImgMax */
    if (n == 0) return 0;

    for (i = 0; i < n; i++) { xOut[i] = xA[i]; yOut[i] = yA[i]; }
    return n;
}

/*
  Find region where image footprint overlaps the output grid.

  Projects the 4 corner control points of the SAR swath into polar-
  stereographic XY, clips the resulting (possibly rotated) quadrilateral
  against the axis-aligned output image rectangle, and converts the
  bounding box of the clipped polygon to pixel index bounds.

  Returns iMin==iMax==jMin==jMax==0 when there is no real overlap.
*/
int32_t getRegion(inputImageStructure *image, int32_t *iMin, int32_t *iMax, int32_t *jMin, int32_t *jMax,
               outputImageStructure *outputImage)
{
    double xPoly[8], yPoly[8];
    double minX, maxX, minY, maxY;
    int n, i;

    /* --- 0. File existence check -------------------------------------------- */
    if (image->file != NULL && strstr(image->file, "nophase") == NULL) {
        char vrtBuf[4096];
        snprintf(vrtBuf, sizeof(vrtBuf), "%s.vrt", image->file);
        if (access(image->file, F_OK) != 0 && access(vrtBuf, F_OK) != 0) {
            fprintf(stderr, "getRegion: file not found %s\n", image->file);
            *iMin = 0; *iMax = 0; *jMin = 0; *jMax = 0;
            return 0;
        }
    }

    n = clipSwathToOutputRect(image, outputImage, xPoly, yPoly);
    if (n == 0) goto no_overlap;

    /* --- 4. Bounding box of clipped polygon ---------------------------------- */
    minX = xPoly[0]; maxX = xPoly[0];
    minY = yPoly[0]; maxY = yPoly[0];
    for (i = 1; i < n; i++)
    {
        if (xPoly[i] < minX) minX = xPoly[i];
        if (xPoly[i] > maxX) maxX = xPoly[i];
        if (yPoly[i] < minY) minY = yPoly[i];
        if (yPoly[i] > maxY) maxY = yPoly[i];
    }

    /* --- 5. Convert to output pixel indices ---------------------------------- */
    /* Subtract 1 from min / add 1 to max to avoid truncation cutting the edge row/col */
    {
        int pad2 = max(2, 10e3/outputImage->deltaX); /* Add a little extra padding to catch any rounding issues at the edges */
        *iMin = (int)((minY - outputImage->originY) / outputImage->deltaY) - pad2;
        *jMin = (int)((minX - outputImage->originX) / outputImage->deltaX) - pad2;
        *iMax = (int)((maxY - outputImage->originY) / outputImage->deltaY) + pad2;
        *jMax = (int)((maxX - outputImage->originX) / outputImage->deltaX) + pad2;
    }
    *iMin = max(*iMin, 0);
    *jMin = max(*jMin, 0);
    *iMax = min(outputImage->ySize, *iMax);
    *jMax = min(outputImage->xSize, *jMax);

    fprintf(stderr, "%sgetRegion: clipped bounds x [%.1f %.1f] y [%.1f %.1f] -> i [%i %i] j [%i %i] %f\n",
            indentRegionOutput ? "\t" : "", minX, maxX, minY, maxY, *iMin, *iMax, *jMin, *jMax, outputImage->deltaX);
    return 1;

no_overlap:
    *iMin = 0; *iMax = 0;
    *jMin = 0; *jMax = 0;
    fprintf(stderr, "%sgetRegion: no overlap\n", indentRegionOutput ? "\t" : "");
    return 0;
}

/*
  Tighter per-row column bound within an already-computed [iMin,iMax) row range
  (from getRegion()). For a long, narrow, diagonal swath, getRegion()'s single
  axis-aligned bbox can be much wider than the true footprint at any given row
  (up to ~L²/2 vs the true L·W footprint area for a track of length L, width
  W<<L). This re-clips the same swath polygon getRegion() uses and, for each
  row, intersects it against that row's horizontal line to get a tight column
  range — cheap O(n) edge math, no per-pixel geocoding.

  jRowMin[i]/jRowMax[i] are filled for i in [iMin,iMax) (caller allocates,
  sized to at least outputImage->ySize). Rows with no polygon intersection
  (should only happen from degenerate geometry) fall back to the full
  [jMinFallback,jMaxFallback) range from getRegion() — never narrower than
  necessary, only ever a missed optimization, not a correctness risk.
*/
void getRowBounds(inputImageStructure *image, outputImageStructure *outputImage,
                   int32_t iMin, int32_t iMax, int32_t jMinFallback, int32_t jMaxFallback,
                   int32_t *jRowMin, int32_t *jRowMax)
{
    /* Must approximate checkLL()'s (llToImageNew.c) 50km EUCLIDEAN
       distance-to-boundary tolerance, not just an x-pad on each row's own
       polygon crossing. A row-local x-pad alone misses points near a
       shallow/near-parallel edge (e.g. the swath's along-track end caps):
       such a point can be within 50km of that edge in true distance while
       being far more than 50km away in same-row x-crossing terms. Fixed by
       dilating in the row (y) direction too: row i's bound is the union of
       every row i' within padRows' raw (unpadded) crossings, each then
       widened by pad2 in x. This is a safe (only ever over-inclusive, never
       under-) approximation of the true 50km disk buffer -- any point within
       50km of the polygon boundary has a nearest boundary point at most 50km
       away in EACH of x and y separately, and that boundary point's own row
       falls within padRows of row i by construction. Confirmed necessary by
       a real A/B diff: the row-only-pad version still dropped ~280k valid
       pixels (vs ~2.5M with no pad-matching at all) before this fix. */
    double xPoly[8], yPoly[8];
    int n = clipSwathToOutputRect(image, outputImage, xPoly, yPoly);
    int pad2 = max(2, 50.0e3 / outputImage->deltaX);
    int padRows = max(2, (int)(50.0e3 / outputImage->deltaY));
    int32_t i;
    double *rawMinX, *rawMaxX; /* unpadded per-row polygon x-crossing, indexed [iMin,iMax) */

    if (n == 0)
    {
        for (i = iMin; i < iMax; i++) { jRowMin[i] = jMinFallback; jRowMax[i] = jMaxFallback; }
        return;
    }

    rawMinX = (double *)malloc((size_t)(iMax - iMin) * sizeof(double));
    rawMaxX = (double *)malloc((size_t)(iMax - iMin) * sizeof(double));
    if (rawMinX == NULL || rawMaxX == NULL)
        error("getRowBounds: malloc failed for rawMinX/rawMaxX\n");

    for (i = iMin; i < iMax; i++)
    {
        double y = outputImage->originY + (double)i * outputImage->deltaY; /* metres, matches xPoly/yPoly */
        double xMinRow = 1.0e30, xMaxRow = -1.0e30;
        int k;

        for (k = 0; k < n; k++)
        {
            int k2 = (k + 1) % n;
            double y1 = yPoly[k], y2 = yPoly[k2];
            if (y1 == y2) continue; /* horizontal edge: contributes nothing an adjacent edge doesn't already cover */
            if ((y1 <= y && y2 >= y) || (y1 >= y && y2 <= y))
            {
                double t = (y - y1) / (y2 - y1);
                double x = xPoly[k] + t * (xPoly[k2] - xPoly[k]);
                if (x < xMinRow) xMinRow = x;
                if (x > xMaxRow) xMaxRow = x;
            }
        }
        rawMinX[i - iMin] = xMinRow; /* xMinRow stays 1.0e30 (i.e. "no crossing") if none found */
        rawMaxX[i - iMin] = xMaxRow;
    }

    for (i = iMin; i < iMax; i++)
    {
        double xMinRow = 1.0e30, xMaxRow = -1.0e30;
        int32_t wLo = max(iMin, i - padRows);
        int32_t wHi = min(iMax - 1, i + padRows);
        int32_t w;

        for (w = wLo; w <= wHi; w++)
        {
            if (rawMinX[w - iMin] < xMinRow) xMinRow = rawMinX[w - iMin];
            if (rawMaxX[w - iMin] > xMaxRow) xMaxRow = rawMaxX[w - iMin];
        }

        if (xMinRow > xMaxRow) /* no row in the window had a crossing */
        {
            jRowMin[i] = jMinFallback;
            jRowMax[i] = jMaxFallback;
            continue;
        }

        {
            int32_t jr0 = (int32_t)((xMinRow - outputImage->originX) / outputImage->deltaX) - pad2;
            int32_t jr1 = (int32_t)((xMaxRow - outputImage->originX) / outputImage->deltaX) + pad2;
            jr0 = max(jr0, jMinFallback);
            jr1 = min(jr1, jMaxFallback);
            if (jr0 > jr1)
            {
                /* degenerate after clamping — fall back rather than risk losing real data */
                jr0 = jMinFallback;
                jr1 = jMaxFallback;
            }
            jRowMin[i] = jr0;
            jRowMax[i] = jr1;
        }
    }

    free(rawMinX);
    free(rawMaxX);
}
