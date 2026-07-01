
void getROffsets(char *phaseFile, tiePointsStructure *tiePoints, inputImageStructure inputImage, Offsets *offsets, int32_t noIonosphere, int32_t skipLoad);

void addOffsetCorrections(inputImageStructure inputImage, tiePointsStructure *tiePoints);

double computeRParams(tiePointsStructure *tiePoints, inputImageStructure inputImage, char *baseFile, Offsets *offsets, int32_t yamlOutput);

void getBaselineFile(char *baselineFile, tiePointsStructure *tiePoints, inputImageStructure inputImage);

void addVelCorrections(inputImageStructure *inputImage, tiePointsStructure *tiePoints);
