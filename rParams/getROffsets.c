#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "rparams.h"
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "gdalIO/gdalIO/grimpgdal.h"
/*
  Input range offset and extract phases for tiepoint locations.
*/
void getROffsets(char *phaseFile, tiePointsStructure *tiePoints, inputImageStructure inputImage, Offsets *offsets, int32_t noIonosphere)
{
	FILE *fp;
	double range, azimuth;
	int32_t i, count;

	/*
	  If ion correction is requested and no baseline has pre-filled correctionFile,
	  peek at the range offset VRT metadata now so checkForIonosphereCorrection (called
	  inside readRangeOffsets) will see a non-empty correctionFile and load the correction.
	*/
	if (!noIonosphere && offsets->rOffCorrection.correctionFile[0] == '\0')
	{
		char vrtBuf[2048];
		char *vrtFile = NULL;
		/* Try explicit .vrt suffix first, then append .vrt */
		if (has_suffix(phaseFile, ".vrt"))
		{
			vrtFile = phaseFile;
		}
		else
		{
			snprintf(vrtBuf, sizeof(vrtBuf), "%s.vrt", phaseFile);
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
				char *ionName = get_value(metaData, "ionosphereRangeOffsetCorrection");
				if (ionName != NULL)
				{
					strncpy(offsets->rOffCorrection.correctionFile, ionName,
					        sizeof(offsets->rOffCorrection.correctionFile) - 1);
					offsets->rOffCorrection.correctionFile[sizeof(offsets->rOffCorrection.correctionFile) - 1] = '\0';
					fprintf(stderr, "getROffsets: found ionosphere correction in VRT: %s\n",
					        offsets->rOffCorrection.correctionFile);
				}
				free_dictionary(metaData);
				GDALClose(hDS);
			}
		}
	}

	/* Get offsets (checkForIonosphereCorrection inside will load correction if correctionFile set) */
	readRangeOffsets(offsets, FALSE, 0.0f, (float)LARGEINT);
	fprintf(stderr, "	OFFSETS READ\n\n");
	/*
	  Interpolate offsets
	*/
	fprintf(stdout, ";;\n;;Tiepoints row column elevation\n;;\n");
	count = 0;
	for (i = 0; i < tiePoints->npts; i++)
	{
		range = (tiePoints->r[i] * inputImage.nRangeLooks - offsets->rO) / offsets->deltaR;
		azimuth = (tiePoints->a[i] * inputImage.nAzimuthLooks - offsets->aO) / offsets->deltaA;
	
		tiePoints->phase[i] = bilinearInterp((float **)offsets->dr, range, azimuth,
											 offsets->nr, offsets->na, -0.99 * LARGEINT, (float)-LARGEINT); /* only scale valid values */
		if (offsets->rOffCorrection.rangeOffsetCorrection != NULL && noIonosphere == FALSE)
		{
			/* SLC pixel coordinates of the tiepoint (not the offset-image index) */
			double rSLC = tiePoints->r[i] * inputImage.nRangeLooks;
			double aSLC = tiePoints->a[i] * inputImage.nAzimuthLooks;
			float ionoCorr = interpolateOffsetIonCorrectionInPixels(&offsets->rOffCorrection,
																	rSLC, aSLC,
														 			-0.99 * LARGEINT, (float)-LARGEINT);
			//fprintf(stderr, "ionoCorr %f for r %f a %f\n", ionoCorr, range, azimuth);
			if (ionoCorr > -0.98 * LARGEINT)
			{
				//fprintf(stderr, "ionoCorr %f phase %f\n", ionoCorr, tiePoints->phase[i]);
				// Ion correction is in pixels consistent with offsets
				tiePoints->phase[i] += ionoCorr;
			}
		}
		if (tiePoints->phase[i] > -0.98 * LARGEINT)
		{
			tiePoints->phase[i] *= inputImage.rangePixelSize / inputImage.nRangeLooks;
			count++;
		}
		/*  Multiply by -1 for left to get RHS*/
		//        if(inputImage.lookDir==LEFT) tiePoints->phase[i] *=-1; 
		if (fabs(tiePoints->phase[i]) < LARGEINT / 10 && tiePoints->quiet == FALSE)
			fprintf(stdout, "; %i  %i  %f %f\n", (int)(tiePoints->r[i] + 0.5), (int)(tiePoints->a[i] + 0.5),
					tiePoints->z[i], tiePoints->phase[i]);
	}
	fprintf(stderr, "count %i\n", count);
	fprintf(stdout, ";&\n");
	return;
}
