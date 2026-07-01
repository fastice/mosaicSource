/*#include "mosaicSource/ers1Code/ers1.h"*/
/*#include "mosaicSource/GeoCode_p/geocode.h"*/

/*
   Input phase image and extract phases for tiepoint locations.
*/
void getPhases(char *phaseFile, tiePointsStructure *tiePoints, inputImageStructure inputImage, int32_t yamlOutput, int32_t verbose);
/*
    Add baseline corrections that were removed in the unwrapped image to
    tiepoint phases.
*/
void addBaselineCorrections(char *baselineFile, tiePointsStructure *tiePoints, inputImageStructure inputImage, int32_t yamlOutput, int32_t verbose);
/*
   Estimate baseline parameters.
*/
void computeBaseline(tiePointsStructure *tiePoints, inputImageStructure inputImage, int32_t yamlOutput, int32_t verbose);
/*
   Compute motion corrections. Writes into phaseOut (parallel array to tiePoints->phase);
   applySquint selects whether the squint(r,a) heading correction is applied.
*/
void addMotionCorrections(inputImageStructure inputImage, tiePointsStructure *tiePoints,
						  int32_t applySquint, double *phaseOut, int32_t verbose);
