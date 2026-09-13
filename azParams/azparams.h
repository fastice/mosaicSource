#include "mosaicSource/common/writeTieResidualsGpkg.h"

void getOffsets(char *phaseFile, tiePointsStructure *tiePoints, inputImageStructure inputImage, Offsets *offsets, int32_t noIonosphere, int32_t skipLoad);

void addOffsetCorrections(inputImageStructure *inputImage, tiePointsStructure *tiePoints);

double computeAzParams(tiePointsStructure *tiePoints, inputImageStructure *inputImage, char *baseFile, Offsets *offsets, int32_t yamlOutput, char *debugFile);
