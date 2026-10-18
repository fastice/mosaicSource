#include <stdio.h>

typedef struct landSatImageType
{
	matchResult matches;
	lsFit fitResult;
	char *maskFile;
	double weight;
	struct landSatImageType *next;
	int32_t nImages;
} landSatImage;

/*
  Per-pair constants for evaluating one Landsat pair's velocity at output pixels.  Shared by the
  Landsat round (makeLandSatMosaic) and the hopper's -landsatRows, so both use one implementation.
*/
typedef struct lsPixelSetupType
{
	double x0poly, y0poly;			 /* centre of the fit polynomial (m) */
	double vscalex, vscaley;		 /* pixels -> m/yr */
	double varXGlobal, varYGlobal;	 /* tie-residual + quantization variance (pixels^2) */
	double pX[3], pY[3];			 /* polynomial corrections (pixels) */
	double tOffCenter;				 /* pair mid-time (julian day) */
} lsPixelSetup;

landSatImage *parseLSInputs(char *inputFile, landSatImage *LSImages, double jd1, double jd2, int32_t timeOverlapFlag, double *minJD, double *maxJD);
void makeLandSatMosaic(landSatImage *LSImages, outputImageStructure *outputImage, float fl);
int32_t landSatPrepareImage(landSatImage *currentImage, outputImageStructure *outputImage, lsPixelSetup *c,
							int32_t *iMin, int32_t *iMax, int32_t *jMin, int32_t *jMax);
int32_t landSatVelocityAt(landSatImage *currentImage, lsPixelSetup *c, double x, double y,
						  double *vx, double *vy, double *varX, double *varY);
