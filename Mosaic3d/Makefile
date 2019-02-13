CC =		gcc
CFLAGS  =	 $(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS=  get3DInputFile.o make3DOffsets.o make3DMosaic.o  makeVhMosaic.o  mosaic3d.o  setup3D.o speckleTrackMosaic.o  writeTieFile.o 

.KEEP_STATE:
all:	$(OBJS)

	mv *.o $(MACHTYPE)-$(OSTYPE)
