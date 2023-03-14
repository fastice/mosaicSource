C =		gcc
ROOTDIR =	/Users/ian
PROGDIR =       $(ROOTDIR)/progs/GIT
INCLUDEPATH =	$(ROOTDIR)/progs/GIT
BINDIR =	$(IHOME)/bin/$(MACHTYPE)
#
CFLAGS =	'-O3 -m32 -I$(INCLUDEPATH) $(COMPILEFLAGS)'
CCFLAGS =  '-O3 -m32 -D$(MACHTYPE) $(COMPILEFLAGS) '
#-Wunused-variable'

CCFLAGS1= -O3 -no-pie
# uncomment to debug
#CFLAGS =	'-g -m32 -I$(INCLUDEPATH) $(COMPILEFLAGS)'
#CCFLAGS =  '-g -m32 -D$(MACHTYPE) $(COMPILEFLAGS)'
#CCFLAGS1= '-g'

COMMON=	common/$(MACHTYPE)-$(OSTYPE)/addIrregData.o \
	                common/$(MACHTYPE)-$(OSTYPE)/bilinearInterp.o \
			common/$(MACHTYPE)-$(OSTYPE)/computeHeading.o \
	                common/$(MACHTYPE)-$(OSTYPE)/computePhiZ.o \
			common/$(MACHTYPE)-$(OSTYPE)/computeScale.o \
			common/$(MACHTYPE)-$(OSTYPE)/computeTiePoints.o \
	                common/$(MACHTYPE)-$(OSTYPE)/computeXYangle.o \
	                common/$(MACHTYPE)-$(OSTYPE)/earthRadiusFunctions.o \
			common/$(MACHTYPE)-$(OSTYPE)/getDataStringSpecial.o \
			common/$(MACHTYPE)-$(OSTYPE)/getBaseline.o \
			common/$(MACHTYPE)-$(OSTYPE)/getHeight.o \
	                common/$(MACHTYPE)-$(OSTYPE)/getIrregData.o \
	                common/$(MACHTYPE)-$(OSTYPE)/getMVhInputFile.o \
			common/$(MACHTYPE)-$(OSTYPE)/getRegion.o \
	                common/$(MACHTYPE)-$(OSTYPE)/getShelfMask.o \
			common/$(MACHTYPE)-$(OSTYPE)/getXYHeight.o \
			common/$(MACHTYPE)-$(OSTYPE)/groundRangeToLLNew.o \
			common/$(MACHTYPE)-$(OSTYPE)/initMatrix.o \
	                common/$(MACHTYPE)-$(OSTYPE)/initRoutines.o \
	                common/$(MACHTYPE)-$(OSTYPE)/interpOffsets.o \
			common/$(MACHTYPE)-$(OSTYPE)/interpPhaseImage.o \
	                common/$(MACHTYPE)-$(OSTYPE)/interpTideDiff.o \
	                common/$(MACHTYPE)-$(OSTYPE)/interpVCorrect.o \
			common/$(MACHTYPE)-$(OSTYPE)/interpXYDEM.o \
			common/$(MACHTYPE)-$(OSTYPE)/julianDay.o \
	                common/$(MACHTYPE)-$(OSTYPE)/llToImageNew.o \
	                common/$(MACHTYPE)-$(OSTYPE)/lltoxy.o \
	                common/$(MACHTYPE)-$(OSTYPE)/lltoxy1.o \
	                common/$(MACHTYPE)-$(OSTYPE)/outputGeocodedImage.o \
			common/$(MACHTYPE)-$(OSTYPE)/parseInputFile.o \
	                common/$(MACHTYPE)-$(OSTYPE)/parseIrregFile.o \
			common/$(MACHTYPE)-$(OSTYPE)/polintVec.o \
	                common/$(MACHTYPE)-$(OSTYPE)/rangeAzimuthToLL.o \
	                common/$(MACHTYPE)-$(OSTYPE)/readOffsets.o \
			common/$(MACHTYPE)-$(OSTYPE)/readOldPar.o \
			common/$(MACHTYPE)-$(OSTYPE)/readShelf.o \
			common/$(MACHTYPE)-$(OSTYPE)/readTiePoints.o \
	                common/$(MACHTYPE)-$(OSTYPE)/readXYDEM.o \
			common/$(MACHTYPE)-$(OSTYPE)/rotateFlowDirectionToRA.o \
			common/$(MACHTYPE)-$(OSTYPE)/rotateFlowDirectionToXY.o \
	                common/$(MACHTYPE)-$(OSTYPE)/scalingFunctions.o \
			common/$(MACHTYPE)-$(OSTYPE)/smlocateZD.o \
			common/$(MACHTYPE)-$(OSTYPE)/svBase.o \
			common/$(MACHTYPE)-$(OSTYPE)/vectorFunc.o \
			common/$(MACHTYPE)-$(OSTYPE)/xyGetZandSlope.o \
			common/$(MACHTYPE)-$(OSTYPE)/xytoll1.o \
			common/$(MACHTYPE)-$(OSTYPE)/xytoll.o

STANDARD =	$(PROGDIR)/clib/$(MACHTYPE)-$(OSTYPE)/standard.o

TRIANGLE =	$(PROGDIR)/triangle/$(MACHTYPE)-$(OSTYPE)/triangle.o

RECIPES  =	$(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/polint.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/nrutil.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/ratint.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/four1.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svdfit.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svdcmp.o \
		$(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svdvar.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svbksb.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/pythag.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/hunt.o


LANDSATCODE =		landsatMosaic/$(MACHTYPE)-$(OSTYPE)/readLSOffsets.o \
					landsatMosaic/$(MACHTYPE)-$(OSTYPE)/parseLSInputs.o \
					landsatMosaic/$(MACHTYPE)-$(OSTYPE)/interpLSTies.o \
					landsatMosaic/$(MACHTYPE)-$(OSTYPE)/xyscale.o \
					landsatMosaic/$(MACHTYPE)-$(OSTYPE)/makeLandSatMosaic.o


TARGETS = mosaic3d siminsar rparams azparams coarsereg tiepoints lltora getlocc geomosaic

all: $(TARGETS)

#********************************************************************************
#********************************** mosaic3d ************************************
#********************************************************************************

MOSAIC3D1 =	Mosaic3d/$(MACHTYPE)-$(OSTYPE)/get3DInputFile.o \
		Mosaic3d/$(MACHTYPE)-$(OSTYPE)/make3DMosaic.o \
	        Mosaic3d/$(MACHTYPE)-$(OSTYPE)/make3DOffsets.o \
                Mosaic3d/$(MACHTYPE)-$(OSTYPE)/makeVhMosaic.o \
                Mosaic3d/$(MACHTYPE)-$(OSTYPE)/setup3D.o \
                Mosaic3d/$(MACHTYPE)-$(OSTYPE)/speckleTrackMosaic.o \
                Mosaic3d/$(MACHTYPE)-$(OSTYPE)/writeTieFile.o

MOSAIC3DDIRS =	Mosaic3d common  landsatMosaic $(PROGDIR)/clib

mosaic3d:	
	@for i in ${MOSAIC3DDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		gcc -m32 $(CCFLAGS1) \
                $(MOSAIC3D1)   $(COMMON)   $(STANDARD) $(RECIPES)  $(TRIANGLE) $(LANDSATCODE) \
                -lm  -o $(BINDIR)/mosaic3d Mosaic3d/$(MACHTYPE)-$(OSTYPE)/mosaic3d.o


#********************************************************************************
#********************************** siminsar *************************************
#********************************************************************************
SIMINSAR =	simInSAR/$(MACHTYPE)-$(OSTYPE)/parseSceneFile.o \
                simInSAR/$(MACHTYPE)-$(OSTYPE)/simInSARimage.o \
                simInSAR/$(MACHTYPE)-$(OSTYPE)/outputSimulatedImage.o \
                simInSAR/$(MACHTYPE)-$(OSTYPE)/getDisplacement.o \
                simInSAR/$(MACHTYPE)-$(OSTYPE)/getSlantRangeDEM.o

SIMINSARDIRS =	simInSAR  common   $(PROGDIR)/clib

siminsar:	
	@for i in ${SIMINSARDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		gcc -m32 $(CCFLAGS1) \
                simInSAR/$(MACHTYPE)-$(OSTYPE)/siminsar.o $(SIMINSAR) $(COMMON) $(STANDARD) $(RECIPES) $(TRIANGLE) -lm  -o $(BINDIR)/siminsar


#********************************************************************************
#********************************** rparams *************************************
#********************************************************************************
RPARAMS =       rParams/$(MACHTYPE)-$(OSTYPE)/getROffsets.o \
                rParams/$(MACHTYPE)-$(OSTYPE)/computeRParams.o \
                rParams/$(MACHTYPE)-$(OSTYPE)/getBaselineFile.o \
                rParams/$(MACHTYPE)-$(OSTYPE)/addVelCorrections.o

RPARAMSDIRS =	 rParams  $(PROGDIR)/cRecipes common

rparams:	
	@for i in ${RPARAMSDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		gcc -m32 $(CCFLAGS1) \
                rParams/$(MACHTYPE)-$(OSTYPE)/rparams.o $(RPARAMS)  $(COMMON) $(STANDARD) $(RECIPES) $(TRIANGLE)  \
                -lm  -o $(BINDIR)/rparams

#********************************************************************************
#********************************** azparams *************************************
#********************************************************************************

AZPARAMS =	azParams/$(MACHTYPE)-$(OSTYPE)/addOffsetCorrections.o \
                azParams/$(MACHTYPE)-$(OSTYPE)/computeAzparams.o \
		azParams/$(MACHTYPE)-$(OSTYPE)/getOffsets.o

AZPARAMSDIRS =	 azParams common  $(PROGDIR)/cRecipes 

azparams:	
	@for i in ${AZPARAMSDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		gcc -m32 $(CCFLAGS1) \
                azParams/$(MACHTYPE)-$(OSTYPE)/azparams.o $(AZPARAMS)  $(COMMON) $(STANDARD) $(RECIPES) $(TRIANGLE)  \
                -lm  -o $(BINDIR)/azparams

#********************************************************************************
#********************************** coarsereg ************************************
#********************************************************************************

COARSEREG =  coarseReg/$(MACHTYPE)-$(OSTYPE)/coarseRegister.o

COARSEREGDIRS =	coarseReg  $(PROGDIR)/cRecipes  common

coarsereg:	
	@for i in ${COARSEREGDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		gcc -m32 $(CCFLAGS1) \
                coarseReg/$(MACHTYPE)-$(OSTYPE)/coarsereg.o $(COARSEREG)  $(COMMON) $(STANDARD) $(RECIPES) $(TRIANGLE)  \
                -lm  -o $(BINDIR)/coarsereg


#********************************************************************************
#********************************** tiePoints *************************************
#********************************************************************************

TIEPOINTS =	tiePoints/$(MACHTYPE)-$(OSTYPE)/getPhases.o \
                tiePoints/$(MACHTYPE)-$(OSTYPE)/addBaselineCorrections.o \
                tiePoints/$(MACHTYPE)-$(OSTYPE)/addMotionCorrections.o \
                tiePoints/$(MACHTYPE)-$(OSTYPE)/computeBaseline.o

TIEPOINTSDIRS =	 tiePoints  $(PROGDIR)/cRecipes common

tiepoints:	
	@for i in ${TIEPOINTSDIRS}; do \
		( 	echo "<<< Descending in directory: $$i -- >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		gcc -m32 $(CCFLAGS1) \
                tiePoints/$(MACHTYPE)-$(OSTYPE)/tiepoints.o $(TIEPOINTS) $(STANDARD) $(TRIANGLE) $(RECIPES)  $(COMMON)  \
		-lm  -o $(BINDIR)/tiepoints

#********************************************************************************
#********************************** lltora *************************************
#********************************************************************************

LLTORADIRS =	LLtoRA tiePoints  $(PROGDIR)/cRecipes  common

lltora:	
	@for i in ${LLTORADIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH)  PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
	        gcc -m32 LLtoRA/$(MACHTYPE)-$(OSTYPE)/lltora.o   $(STANDARD) $(RECIPES)  $(COMMON) $(TRIANGLE) \
                      -lm  -o $(BINDIR)/lltora

#********************************************************************************
#********************************** getlocc *************************************
#********************************************************************************

GETLOCC =	getLocC/$(MACHTYPE)-$(OSTYPE)/centerLL.o getLocC/$(MACHTYPE)-$(OSTYPE)/correctTime.o \
	getLocC/$(MACHTYPE)-$(OSTYPE)/glatlon.o

GETLOCCDIRS =	getLocC common

getlocc:
	@for i in ${GETLOCCDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		gcc -m32 $(CCFLAGS1) \
                getLocC/$(MACHTYPE)-$(OSTYPE)/getlocc.o $(GETLOCC) $(COMMON) $(TRIANGLE) $(STANDARD) $(RECIPES) \
                 -lm -o  $(BINDIR)/getlocc


#********************************************************************************
#********************************** geomosaic *************************************
#********************************************************************************

GEOMOSAIC =	geoMosaic/$(MACHTYPE)-$(OSTYPE)/makeGeoMosaic.o \
                geoMosaic/$(MACHTYPE)-$(OSTYPE)/processInputFileGeo.o

GEOMOSAICDIRS =	geoMosaic common landsatMosaic

geomosaic:	
	@for i in ${GEOMOSAICDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		gcc -m32 $(CCFLAGS1)  \
		 $(GEOMOSAIC) $(COMMON)   $(STANDARD) $(RECIPES)  $(TRIANGLE) \
		landsatMosaic/$(MACHTYPE)-$(OSTYPE)/xyscale.o  \
                -lm  -o $(BINDIR)/geomosaic geoMosaic/$(MACHTYPE)-$(OSTYPE)/geomosaic.o
