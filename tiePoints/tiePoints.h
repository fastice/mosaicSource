/*#include "mosaicSource/ers1Code/ers1.h"*/
/*#include "mosaicSource/GeoCode_p/geocode.h"*/

/* ionosphereMode values -- same convention as rparams (see rParams/rparams.c) */
#define ION_AUTO 0	/* default when -ionosphere is given: fit both ways, keep the better */
#define ION_NONE 1	/* -noIonosphere: never apply the correction */
#define ION_FORCE 2 /* -forceIonosphere: always apply it */
/* Fractional sigma improvement the corrected fit must achieve before it is used */
#define IONSIGMAMARGIN 0.05

/*
   Input phase image and extract phases for tiepoint locations.
*/
void getPhases(char *phaseFile, tiePointsStructure *tiePoints, inputImageStructure inputImage, int32_t yamlOutput, int32_t verbose);
/*
   Input ionospheric phase image (radians, same grid as the unwrapped phase) and extract
   values for tiepoint locations into tiePoints->ionPhase.
*/
void getIonosphere(char *ionFile, tiePointsStructure *tiePoints, inputImageStructure inputImage);
/*
    Add baseline corrections that were removed in the unwrapped image to
    tiepoint phases.
*/
void addBaselineCorrections(char *baselineFile, tiePointsStructure *tiePoints, inputImageStructure inputImage, int32_t yamlOutput, int32_t verbose);
/*
   Estimate baseline parameters.
*/
void computeBaseline(tiePointsStructure *tiePoints, inputImageStructure inputImage, int32_t yamlOutput, int32_t verbose, char *debugFile,
					 int32_t ionosphereMode, double ionSigmaMargin, char *ionosphereFile);
/*
   Compute motion corrections. Writes into phaseOut (parallel array to tiePoints->phase);
   applySquint selects whether the squint(r,a) heading correction is applied.
*/
void addMotionCorrections(inputImageStructure inputImage, tiePointsStructure *tiePoints,
						  int32_t applySquint, double *phaseOut, int32_t verbose);
