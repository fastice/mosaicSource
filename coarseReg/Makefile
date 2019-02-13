CC =		gcc
CFLAGS  =	 $(FLAGS) -c -I$(INCLUDEPATH)

OBJS= coarseRegister.o coarsereg.o

.KEEP_STATE:
all:	$(OBJS)

	mv *.o $(MACHTYPE)-$(OSTYPE)
