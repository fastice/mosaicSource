CC =		gcc
MACHTYPE =	$(shell uname -m)
OSTYPE =	$(shell uname -s)
OBJDIR =	$(MACHTYPE)-$(OSTYPE)

$(shell mkdir -p $(OBJDIR))

CFLAGS  =	$(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH) -MMD -MP

SRCS =	addIrregData.c bilinearInterp.c buffers.c computeHeading.c computePhiZ.c \
	computeScale.c computeTiePoints.c computeXYangle.c \
	earthRadiusFunctions.c geojsonCode.c getAzimuthBoundsForXYBox.c getBaseline.c \
	getDataStringSpecial.c getHeight.c getIrregData.c \
	getMVhInputFile.c getRegion.c getShelfMask.c getXYHeight.c groundRangeToLLNew.c julianDay.c \
	initMatrix.c initRoutines.c interpOffsets.c interpPhaseImage.c interpTideDiff.c interpVCorrect.c interpXYDEM.c interpolateOffsetCorrection.c \
	llToImageNew.c lltoxy.c lltoxy1.c outputGeocodedImage.c parseInputFile.c parseIrregFile.c parseSLCvrt.c polintVec.c \
	rangeAzimuthToLL.c readOffsets.c readOffsetCorrection.c readOldPar.c readShelf.c readTiePoints.c readXYDEM.c rotateFlowDirectionToRA.c rotateFlowDirectionToXY.c \
	scalingFunctions.c smlocateZD.c svBase.c vectorFunc.c xyGetZandSlope.c xytoll1.c xytoll.c

OBJS =	$(SRCS:%.c=$(OBJDIR)/%.o)

all:	$(OBJS)

-include $(OBJS:.o=.d)

.KEEP_STATE:


$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $< -o $@
