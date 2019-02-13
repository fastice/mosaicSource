
  void getROffsets( char *phaseFile, tiePointsStructure *tiePoints,  inputImageStructure inputImage,Offsets *offsets);

   void addOffsetCorrections(inputImageStructure inputImage,  tiePointsStructure *tiePoints);
void computeRParams( tiePointsStructure *tiePoints, inputImageStructure inputImage, char *baseFile ,	Offsets *offsets);
     void getBaselineFile(char *baselineFile,tiePointsStructure *tiePoints, inputImageStructure inputImage);
  void addVelCorrections(inputImageStructure *inputImage,   tiePointsStructure *tiePoints);


  

