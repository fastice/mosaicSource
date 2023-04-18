/*#include "mosaicSource/ers1Code/ers1.h"*/
/*#include "mosaicSource/GeoCode_p/geocode.h"*/

/*
   Input phase image and extract phases for tiepoint locations.
*/
void getPhases(char *phaseFile, tiePointsStructure *tiePoints, inputImageStructure inputImage);
/*
    Add baseline corrections that were removed in the unwrapped image to
    tiepoint phases.
*/
void addBaselineCorrections(char *baselineFile, tiePointsStructure *tiePoints, inputImageStructure inputImage);
/*
   Estimate baseline parameters.
*/
void computeBaseline(tiePointsStructure *tiePoints, inputImageStructure inputImage);
/*
   Compute motion corrections
*/
void addMotionCorrections(inputImageStructure inputImage, tiePointsStructure *tiePoints);
