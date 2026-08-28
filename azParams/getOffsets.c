#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "azparams.h"
#include <stdlib.h>
#include "math.h"
/*
   Input azimuth offsets  image and extract phases for tiepoint locations.
*/
void getOffsets(char *phaseFile, tiePointsStructure *tiePoints, inputImageStructure inputImage, Offsets *offsets, int32_t skipLoad)
{
   FILE *fp;
   double range, azimuth;
   int32_t i, count = 0;
   /*
      Init image. skipLoad reuses offsets->da already read from disk (the
      -runFile multi-run path reads azimuth.offsets once, then re-interpolates
      the tiepoints below every run). The interpolation loop always runs.
   */
   if (!skipLoad)
      readAzimuthOffsets(offsets);
   offsets->azInit = FALSE;
   /*
       Interpolate offsets
   */
   /* 2026-06-17: gate header/footer on !quiet so yaml mode (which runs with -quiet) gets clean stdout */
   if (tiePoints->quiet == FALSE)
      fprintf(stdout, ";;\n;;Tiepoints row column elevation\n;;\n");
   for (i = 0; i < tiePoints->npts; i++)
   {
      range = (tiePoints->r[i] * inputImage.nRangeLooks - offsets->rO) / offsets->deltaR;
      azimuth = (tiePoints->a[i] * inputImage.nAzimuthLooks - offsets->aO) / offsets->deltaA;
      /* only scale valide values */
      tiePoints->phase[i] = bilinearInterp((float **)offsets->da, range, azimuth, offsets->nr, offsets->na, -0.99 * LARGEINT, (float)-LARGEINT);
      if (tiePoints->phase[i] > -0.98 * LARGEINT)
      {
         tiePoints->phase[i] *= inputImage.azimuthPixelSize / inputImage.nAzimuthLooks;
         count++;
      }
      /*    Multiply by -1 for left to get RHS????	     */
      /* if(inputImage.lookDir==LEFT) tiePoints->phase[i] *=-1; */
      if (fabs(tiePoints->phase[i]) < (LARGEINT) && tiePoints->quiet == FALSE)
         fprintf(stdout, "; %i  %i  %f %f\n", (int)(tiePoints->r[i] + 0.5), (int)(tiePoints->a[i] + 0.5),
                 tiePoints->z[i], (float)tiePoints->phase[i]);
   }
   fprintf(stderr, "count %i\n", count);
   if (tiePoints->quiet == FALSE)
      fprintf(stdout, ";&\n");
   /*    inputImage->nRangeLooks=tiePoints->deltaR;*/
   return;
}
