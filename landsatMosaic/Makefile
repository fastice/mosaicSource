CC =		gcc
CFLAGS  =	 $(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS=	parseLSInputs.o readLSOffsets.o makeLandSatMosaic.o interpLSTies.o xyscale.o

MACHTYPE = $(shell uname -m)
OSTYPE = $(shell uname -s)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

.KEEP_STATE:
all:	$(OBJS)
	mv *.o $(MACHTYPE)-$(OSTYPE)
