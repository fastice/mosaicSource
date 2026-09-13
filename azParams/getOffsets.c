#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "azparams.h"
#include <stdlib.h>
#include "math.h"
#include <unistd.h>
#include "gdalIO/gdalIO/grimpgdal.h"
/*
   Input azimuth offsets  image and extract phases for tiepoint locations.

   noIonosphere == TRUE suppresses the azimuth ionosphere correction even when
   the VRT names one.  Mirrors getROffsets() on the range side.
*/
void getOffsets(char *phaseFile, tiePointsStructure *tiePoints, inputImageStructure inputImage, Offsets *offsets, int32_t noIonosphere, int32_t skipLoad)
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
   {
      /*
        If the correction is wanted and nothing has pre-filled correctionFile,
        peek at the azimuth offsets VRT now so checkForAzimuthIonosphereCorrection
        (inside readAzimuthOffsets) sees a non-empty name and loads the raster.
      */
      if (!noIonosphere && offsets->aOffCorrection.correctionFile[0] == '\0')
      {
         char vrtBuf[2048];
         char *vrtFile = NULL;
         if (has_suffix(offsets->file, ".vrt"))
         {
            vrtFile = offsets->file;
         }
         else
         {
            snprintf(vrtBuf, sizeof(vrtBuf), "%s.vrt", offsets->file);
            if (access(vrtBuf, F_OK) == 0)
               vrtFile = vrtBuf;
         }
         if (vrtFile != NULL)
         {
            GDALDatasetH hDS = GDALOpen(vrtFile, GDAL_OF_READONLY);
            if (hDS != NULL)
            {
               dictNode *metaData = NULL;
               readDataSetMetaData(hDS, &metaData);
               char *ionName = get_value(metaData, "ionosphereAzimuthOffsetCorrection");
               if (ionName != NULL)
               {
                  strncpy(offsets->aOffCorrection.correctionFile, ionName,
                          sizeof(offsets->aOffCorrection.correctionFile) - 1);
                  offsets->aOffCorrection.correctionFile[sizeof(offsets->aOffCorrection.correctionFile) - 1] = '\0';
                  fprintf(stderr, "getOffsets: found azimuth ionosphere correction in VRT: %s\n",
                          offsets->aOffCorrection.correctionFile);
               }
               free_dictionary(metaData);
               GDALClose(hDS);
            }
         }
      }
      readAzimuthOffsets(offsets);
   }
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
      if (offsets->aOffCorrection.azimuthOffsetCorrection != NULL && noIonosphere == FALSE)
      {
         /* SLC pixel coordinates of the tiepoint (not the offset-image index) */
         double rSLC = tiePoints->r[i] * inputImage.nRangeLooks;
         double aSLC = tiePoints->a[i] * inputImage.nAzimuthLooks;
         float ionoCorr = interpolateAzOffsetIonCorrectionInPixels(&offsets->aOffCorrection,
                                                                 rSLC, aSLC,
                                                                 -0.99 * LARGEINT, (float)-LARGEINT);
         /* Correction is in SLC azimuth pixels, same units as da; ADD it, the
            same sign convention the range screen uses. */
         if (ionoCorr > -0.98 * LARGEINT && tiePoints->phase[i] > -0.98 * LARGEINT)
            tiePoints->phase[i] += ionoCorr;
      }
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
