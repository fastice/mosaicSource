CC =		gcc
CFLAGS  =	$(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS= siminsar.o parseSceneFile.o  simInSARimage.o   outputSimulatedImage.o getDisplacement.o  getSlantRangeDEM.o 
.KEEP_STATE:
all:	$(OBJS)





	mv *.o $(MACHTYPE)-$(OSTYPE)
