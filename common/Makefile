CC =		gcc
CFLAGS  =	 $(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS=  addIrregData.o  bilinearInterp.o  computeHeading.o computePhiZ.o  computeScale.o computeTiePoints.o computeXYangle.o   \
	getBaseline.o  getDataStringSpecial.o getHeight.o  getIrregData.o getMVhInputFile.o  getRegion.o getShelfMask.o getXYHeight.o  \
	initMatrix.o  initRoutines.o  interpOffsets.o  interpPhaseImage.o  interpTideDiff.o interpXYDEM.o \
	llToImageNew.o lltoxy.o lltoxy1.o outputGeocodedImage.o parseInputFile.o  parseIrregFile.o polintVec.o \
	readOffsets.o readShelf.o   readTiePoints.o  readXYDEM.o rotateFlowDirectionToRA.o rotateFlowDirectionToXY.o \
	scalingFunctions.o smlocateZD.o svBase.o xyGetZandSlope.o  xytoll1.o  xytoll.o


.KEEP_STATE:
all:	$(OBJS)

	mv *.o $(MACHTYPE)-$(OSTYPE)
