CC =		gcc
CFLAGS  =	 $(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS= addOffsetCorrections.o azparams.o computeAzparams.o getOffsets.o   

.KEEP_STATE:
all:	$(OBJS)

	mv *.o $(MACHTYPE)-$(OSTYPE)
