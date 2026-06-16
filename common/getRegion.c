#include "common.h"

/*
  Clip a convex polygon against a single half-plane  nx*x + ny*y >= d.
  Input polygon has nIn vertices; output written to xOut/yOut.
  Returns number of output vertices (0 = entirely outside).
*/
static int clipByHalfPlane(double *xIn, double *yIn, int nIn,
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
  Find region where image footprint overlaps the output grid.

  Projects the 4 corner control points of the SAR swath into polar-
  stereographic XY, clips the resulting (possibly rotated) quadrilateral
  against the axis-aligned output image rectangle, and converts the
  bounding box of the clipped polygon to pixel index bounds.

  Returns iMin==iMax==jMin==jMax==0 when there is no real overlap.
*/
void getRegion(inputImageStructure *image, int32_t *iMin, int32_t *iMax, int32_t *jMin, int32_t *jMax,
               outputImageStructure *outputImage)
{
    extern double Rotation;
    double xSwath[4], ySwath[4];   /* swath corners in metres */
    double xA[8], yA[8];           /* work buffers for polygon clipping */
    double xB[8], yB[8];
    double xa1, ya1;
    double xImgMin, xImgMax, yImgMin, yImgMax;
    double minX, maxX, minY, maxY;
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
    /* --- 2. Output image rectangle in metres --------------------------------- */
    pad = 15000.0;
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
    if (n == 0) goto no_overlap;
    n = clipByHalfPlane(xB, yB, n, xA, yA, -1.0,  0.0, -xImgMax); /* x <= xImgMax */
    if (n == 0) goto no_overlap;
    n = clipByHalfPlane(xA, yA, n, xB, yB,  0.0,  1.0, yImgMin);  /* y >= yImgMin */
    if (n == 0) goto no_overlap;
    n = clipByHalfPlane(xB, yB, n, xA, yA,  0.0, -1.0, -yImgMax); /* y <= yImgMax */
    if (n == 0) goto no_overlap;

    /* --- 4. Bounding box of clipped polygon ---------------------------------- */
    minX = xA[0]; maxX = xA[0];
    minY = yA[0]; maxY = yA[0];
    for (i = 1; i < n; i++)
    {
        if (xA[i] < minX) minX = xA[i];
        if (xA[i] > maxX) maxX = xA[i];
        if (yA[i] < minY) minY = yA[i];
        if (yA[i] > maxY) maxY = yA[i];
    }

    /* --- 5. Convert to output pixel indices ---------------------------------- */
    /* Subtract 1 from min / add 1 to max to avoid truncation cutting the edge row/col */
	int pad2 = max(500, 10e3/outputImage->deltaX); /* Add a little extra padding to catch any rounding issues at the edges */
    *iMin = (int)((minY - outputImage->originY) / outputImage->deltaY) - pad2;
    *jMin = (int)((minX - outputImage->originX) / outputImage->deltaX) - pad2;
    *iMax = (int)((maxY - outputImage->originY) / outputImage->deltaY) + pad2;
    *jMax = (int)((maxX - outputImage->originX) / outputImage->deltaX) + pad2;
    *iMin = max(*iMin, 0);
    *jMin = max(*jMin, 0);
    *iMax = min(outputImage->ySize, *iMax);
    *jMax = min(outputImage->xSize, *jMax);

    fprintf(stderr, "getRegion: clipped bounds x [%.1f %.1f] y [%.1f %.1f] -> i [%i %i] j [%i %i] %f\n",
            minX, maxX, minY, maxY, *iMin, *iMax, *jMin, *jMax, outputImage->deltaX);
    return;

no_overlap:
    *iMin = 0; *iMax = 0;
    *jMin = 0; *jMax = 0;
    fprintf(stderr, "getRegion: no overlap\n");
}
