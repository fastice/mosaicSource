CC =		gcc
CFLAGS  =	 $(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS= addBaselineCorrections.o addMotionCorrections.o computeBaseline.o \
	getPhases.o tiepoints.o  


.KEEP_STATE:
all:	$(OBJS)

	mv *.o $(MACHTYPE)-$(OSTYPE)
