CC =		gcc
CFLAGS  =	 $(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS=  get3DInputFile.o make3DOffsets.o make3DMosaic.o  makeVhMosaic.o  mosaic3d.o  setup3D.o speckleTrackMosaic.o  writeTieFile.o 

MACHTYPE = $(shell uname -m)
OSTYPE = $(shell uname -s)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

.KEEP_STATE:
all:	$(OBJS)

	mv *.o $(MACHTYPE)-$(OSTYPE)
