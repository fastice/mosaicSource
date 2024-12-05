CC =		gcc
CFLAGS  =	 $(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS= rparams.o getROffsets.o computeRParams.o getBaselineFile.o \
      addVelCorrections.o

MACHTYPE = $(shell uname -m)
OSTYPE = $(shell uname -s)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

.KEEP_STATE:
all:	$(OBJS)

	mv *.o $(MACHTYPE)-$(OSTYPE)
