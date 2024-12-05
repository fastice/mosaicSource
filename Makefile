#
# Modify this section to point to where stuff is.
# Current names are for a specific file system.
# ROOTDIR: root directory for code (e.g. /Users/username directory). Likely should change for linux
# PROGDIR: location for top source code directory (default ROOTDIR/progs/GIT64)
# BINHOME: root directory for binaries
# BINNAME: archetecture dependent basename for bin dir
# BINDIR: directory for binaries (default BINHOME/bin/BINNAME) (will create if doesn't exist)
# INCLUDEPATH: include path (default PROGDIR anything else could cause a problem)
# Various directors can be overridden with environment variable or from make command
# make BINHOME=/a/path/to/binhome
#
# Base directory for code
USER =	$(shell id -u -n)
MACHTYPE = $(shell uname -m)
OSTYPE = $(shelf uname -s)
#
# Default rootdir
ifneq ($(ROOTDIR)),)
	ROOTDIR =	$(dir $(CURDIR))
endif

$(info  "USER is $(USER)")

$(info ROOTDIR="$(ROOTDIR)")
# Default root for source code
ifneq ($(PROGDIR)),)
#	PROGDIR =       $(ROOTDIR)/progs/GIT64
	PROGDIR =       $(dir $(CURDIR))
endif
$(info PROGDIR ="$(PROGDIR)")
#
# Default location root for compiled programs
ifneq ($(BINHOME)),)
	BINHOME =		~$(USER)
endif
$(info BINHOME="$(BINHOME)")
#
# For historical reasons, can compile with 32-bit memory model using MEM=-m32
# In almost all cases, should be compiled as 64bit.
ifneq ($(MEM),-m32)
	BINNAME=	$(MACHTYPE)
	FFTDIR = $(MACHTYPE)-$(OSTYPE)
else
	BINNAME =	i386
	FFTDIR = i386-$(OSTYPE)
endif
$(info "Machtype $(MACHTYPE)")
$(info "BINNAME = $(BINNAME)")
#
# Default binary directory
ifneq ($(BINDIR)),)
	BINDIR =	$(BINHOME)/bin/$(BINNAME)
endif
$(info "BINDIR = $(BINDIR)")
#
# Create bin dir if it doesn't exist
$(shell mkdir -p $(BINDIR))
#
#
# Default include path
ifneq ($(INCLUDEPATH)),)
	INCLUDEPATH =	$(PROGDIR)
endif
$(info INCLUDEPATH ="$(INCLUDEPATH)")
#
# Compiler stuff
#
C =		gcc
#
CFLAGS =	'-O3 $(MEM) -I$(INCLUDEPATH) $(COMPILEFLAGS)'
CCFLAGS =  '-O3 $(MEM) $(COMPILEFLAGS) '
GDAL = -lgdal -lcurl  -lsqlite3 -llzma -lpoppler -lopenjp2 -lssh2 -llcms2
#-Wunused-variable'
#
CCFLAGS1= -O3 
#-no-pie
# uncomment to debug
#CFLAGS =	'-g $(MEM) -I$(INCLUDEPATH) $(COMPILEFLAGS)'
#CCFLAGS =  '-g $(MEM) -D$(MACHTYPE) $(COMPILEFLAGS)'
#CCFLAGS1= '-g'

ifneq ($(OSTYPE),darwin)
	NOPIE =	-no-pie
endif
$(info NOPIE ="$(NOPIE)")

COMMON=	common/$(MACHTYPE)-$(OSTYPE)/addIrregData.o \
	    common/$(MACHTYPE)-$(OSTYPE)/bilinearInterp.o \
		common/$(MACHTYPE)-$(OSTYPE)/computeHeading.o \
	    common/$(MACHTYPE)-$(OSTYPE)/computePhiZ.o \
		common/$(MACHTYPE)-$(OSTYPE)/computeScale.o \
		common/$(MACHTYPE)-$(OSTYPE)/computeTiePoints.o \
	    common/$(MACHTYPE)-$(OSTYPE)/computeXYangle.o \
	    common/$(MACHTYPE)-$(OSTYPE)/earthRadiusFunctions.o \
		common/$(MACHTYPE)-$(OSTYPE)/geojsonCode.o \
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

GDALIO = 	$(PROGDIR)/gdalIO/gdalIO/$(MACHTYPE)-$(OSTYPE)/gdalIO.o \
			$(PROGDIR)/gdalIO/gdalIO/$(MACHTYPE)-$(OSTYPE)/dictionaryCode.o

TARGETS = mosaic3d siminsar rparams azparams coarsereg tiepoints lltora getlocc geomosaic

all: $(TARGETS)

TESTGEO =	testGeo/$(MACHTYPE)-$(OSTYPE)/testgeo.o
TESTGEODIRS =	testGeo $(PROGDIR)/gdalIO/gdalIO  $(PROGDIR)/clib  $(PROGDIR)/mosaicSource/common

testgeo:	
	@for i in ${TESTGEODIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0; \
			cd $(PROGDIR); \
		); done
		g++ $(MEM) $(CCFLAGS1) $(NOPIE) \
                $(TESTGEO) $(GDALIO) $(COMMON) $(STANDARD) $(RECIPES)  $(TRIANGLE)  \
				-lm  -lgdal -o $(BINDIR)/testgeo  -L/usr/lib 	

OFFSETVRT =	offsetVRT/$(MACHTYPE)-$(OSTYPE)/offsetVRT.o
OFFSETVRTDIRS =	offsetVRT $(PROGDIR)/gdalIO/gdalIO  $(PROGDIR)/clib  $(PROGDIR)/mosaicSource/common

offsetvrt:	
	@for i in ${OFFSETVRTDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0; \
			cd $(PROGDIR); \
		); done
		g++ $(MEM) $(CCFLAGS1) $(NOPIE) \
                $(OFFSETVRT) $(GDALIO) $(COMMON) $(STANDARD) $(RECIPES)  $(TRIANGLE)  \
                -lm  -lgdal -o $(BINDIR)/offsetvrt  -L/usr/lib 	
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

MOSAIC3DDIRS =	Mosaic3d common  $(PROGDIR)/gdalIO/gdalIO landsatMosaic $(PROGDIR)/triangle $(PROGDIR)/clib  $(PROGDIR)/cRecipes

mosaic3d:	
	@for i in ${MOSAIC3DDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		g++ $(MEM) $(CCFLAGS1) \
                $(MOSAIC3D1)   $(COMMON) $(STANDARD) $(RECIPES)  $(TRIANGLE) $(LANDSATCODE) $(GDALIO) \
                -lm $(GDAL) -o $(BINDIR)/mosaic3d Mosaic3d/$(MACHTYPE)-$(OSTYPE)/mosaic3d.o


#********************************************************************************
#********************************** siminsar *************************************
#********************************************************************************
SIMINSAR =	simInSAR/$(MACHTYPE)-$(OSTYPE)/parseSceneFile.o \
                simInSAR/$(MACHTYPE)-$(OSTYPE)/simInSARimage.o \
                simInSAR/$(MACHTYPE)-$(OSTYPE)/outputSimulatedImage.o \
                simInSAR/$(MACHTYPE)-$(OSTYPE)/getDisplacement.o \
                simInSAR/$(MACHTYPE)-$(OSTYPE)/getSlantRangeDEM.o

SIMINSARDIRS =	simInSAR  common  $(PROGDIR)/gdalIO/gdalIO  $(PROGDIR)/clib $(PROGDIR)/triangle $(PROGDIR)/clib  $(PROGDIR)/cRecipes

siminsar:	
	@for i in ${SIMINSARDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		g++ $(MEM) $(CCFLAGS1) \
                simInSAR/$(MACHTYPE)-$(OSTYPE)/siminsar.o $(SIMINSAR) $(COMMON)  $(TRIANGLE) $(STANDARD) $(RECIPES) $(GDALIO) \
				-lm -lgdal -o $(BINDIR)/siminsar


#********************************************************************************
#********************************** rparams *************************************
#********************************************************************************
RPARAMS =       rParams/$(MACHTYPE)-$(OSTYPE)/getROffsets.o \
                rParams/$(MACHTYPE)-$(OSTYPE)/computeRParams.o \
                rParams/$(MACHTYPE)-$(OSTYPE)/getBaselineFile.o \
                rParams/$(MACHTYPE)-$(OSTYPE)/addVelCorrections.o

RPARAMSDIRS =	 rParams common $(PROGDIR)/gdalIO/gdalIO $(PROGDIR)/triangle $(PROGDIR)/clib  $(PROGDIR)/cRecipes

rparams:	
	@for i in ${RPARAMSDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		g++ $(MEM) $(CCFLAGS1) \
                rParams/$(MACHTYPE)-$(OSTYPE)/rparams.o $(RPARAMS)  $(COMMON) $(STANDARD) $(RECIPES) $(TRIANGLE) $(GDALIO)  \
                -lm  -lgdal -o $(BINDIR)/rparams

#********************************************************************************
#********************************** azparams *************************************
#********************************************************************************

AZPARAMS =	azParams/$(MACHTYPE)-$(OSTYPE)/addOffsetCorrections.o \
            azParams/$(MACHTYPE)-$(OSTYPE)/computeAzparams.o \
			azParams/$(MACHTYPE)-$(OSTYPE)/getOffsets.o

AZPARAMSDIRS =	 azParams common $(PROGDIR)/gdalIO/gdalIO $(PROGDIR)/cRecipes $(PROGDIR)/triangle $(PROGDIR)/clib  $(PROGDIR)/cRecipes  

azparams:	
	@for i in ${AZPARAMSDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		g++ $(MEM) $(CCFLAGS1) \
                azParams/$(MACHTYPE)-$(OSTYPE)/azparams.o $(AZPARAMS)  $(COMMON) $(STANDARD) $(RECIPES) $(TRIANGLE) $(GDALIO) \
                -lm  -lgdal -o $(BINDIR)/azparams
#********************************************************************************
#********************************** computeBaseline *************************************
#********************************************************************************

BASELINEDIRS =	 computeBaseline common $(PROGDIR)/gdalIO/gdalIO $(PROGDIR)/triangle $(PROGDIR)/clib  $(PROGDIR)/cRecipes

computebaseline:	
	@for i in ${BASELINEDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		g++ $(MEM) $(CCFLAGS1) \
                computeBaseline/$(MACHTYPE)-$(OSTYPE)/computebaseline.o   $(COMMON) $(STANDARD) $(RECIPES) $(TRIANGLE) $(GDALIO)  \
                -lm  -lgdal -o $(BINDIR)/computebaseline

#********************************************************************************
#********************************** coarsereg ************************************
#********************************************************************************

COARSEREG =  coarseReg/$(MACHTYPE)-$(OSTYPE)/coarseRegister.o

COARSEREGDIRS =	coarseReg  $(PROGDIR)/cRecipes  common $(PROGDIR)/triangle $(PROGDIR)/clib  $(PROGDIR)/cRecipes

coarsereg:	
	@for i in ${COARSEREGDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		g++ $(MEM) $(CCFLAGS1) \
                coarseReg/$(MACHTYPE)-$(OSTYPE)/coarsereg.o $(COARSEREG)  $(COMMON) $(STANDARD) $(RECIPES) $(TRIANGLE) $(GDALIO)  \
                -lm -lgdal -o $(BINDIR)/coarsereg


#********************************************************************************
#********************************** tiePoints *************************************
#********************************************************************************

TIEPOINTS =	tiePoints/$(MACHTYPE)-$(OSTYPE)/getPhases.o \
                tiePoints/$(MACHTYPE)-$(OSTYPE)/addBaselineCorrections.o \
                tiePoints/$(MACHTYPE)-$(OSTYPE)/addMotionCorrections.o \
                tiePoints/$(MACHTYPE)-$(OSTYPE)/computeBaseline.o

TIEPOINTSDIRS =	 tiePoints  $(PROGDIR)/cRecipes common $(PROGDIR)/triangle $(PROGDIR)/clib  $(PROGDIR)/cRecipes

tiepoints:	
	@for i in ${TIEPOINTSDIRS}; do \
		( 	echo "<<< Descending in directory: $$i -- >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		g++ $(MEM) $(CCFLAGS1) \
                tiePoints/$(MACHTYPE)-$(OSTYPE)/tiepoints.o $(TIEPOINTS) $(STANDARD) $(TRIANGLE) $(RECIPES)  $(COMMON) $(GDALIO) \
		-lm -lgdal -o $(BINDIR)/tiepoints

#********************************************************************************
#********************************** lltora *************************************
#********************************************************************************

LLTORADIRS =	LLtoRA tiePoints  $(PROGDIR)/cRecipes  common $(PROGDIR)/triangle $(PROGDIR)/clib  $(PROGDIR)/cRecipes

lltora:	
	@for i in ${LLTORADIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH)  PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
	        g++ $(MEM) LLtoRA/$(MACHTYPE)-$(OSTYPE)/lltora.o   $(STANDARD) $(RECIPES)  $(COMMON) $(TRIANGLE) $(GDALIO) \
                      -lm -lgdal -o $(BINDIR)/lltora

#********************************************************************************
#********************************** getlocc *************************************
#********************************************************************************

GETLOCC =	getLocC/$(MACHTYPE)-$(OSTYPE)/centerLL.o \
			getLocC/$(MACHTYPE)-$(OSTYPE)/correctTime.o \
			getLocC/$(MACHTYPE)-$(OSTYPE)/glatlon.o \

GETLOCCDIRS =	getLocC common $(PROGDIR)/triangle $(PROGDIR)/clib  $(PROGDIR)/cRecipes

getlocc:
	@for i in ${GETLOCCDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		g++ $(MEM) $(CCFLAGS1) \
                getLocC/$(MACHTYPE)-$(OSTYPE)/getlocc.o $(GETLOCC) $(COMMON) $(TRIANGLE) $(STANDARD) $(RECIPES) $(GDALIO) \
                 -lm -lgdal -o  $(BINDIR)/getlocc


#********************************************************************************
#********************************** geomosaic *************************************
#********************************************************************************

GEOMOSAIC =	geoMosaic/$(MACHTYPE)-$(OSTYPE)/makeGeoMosaic.o \
                geoMosaic/$(MACHTYPE)-$(OSTYPE)/processInputFileGeo.o

GEOMOSAICDIRS =	geoMosaic common landsatMosaic $(PROGDIR)/triangle $(PROGDIR)/clib  $(PROGDIR)/cRecipes

geomosaic:	
	@for i in ${GEOMOSAICDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/mosaicSource; \
		); done
		g++ $(MEM) $(CCFLAGS1)  \
		 $(GEOMOSAIC) $(COMMON)   $(STANDARD) $(RECIPES)  $(TRIANGLE) $(GDALIO) \
		landsatMosaic/$(MACHTYPE)-$(OSTYPE)/xyscale.o  \
                -lm -lgdal -o $(BINDIR)/geomosaic geoMosaic/$(MACHTYPE)-$(OSTYPE)/geomosaic.o
