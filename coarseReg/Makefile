CC =		gcc
CFLAGS  =	 $(FLAGS) -c -I$(INCLUDEPATH)

OBJS= coarseRegister.o coarsereg.o


MACHTYPE = $(shell uname -m)
OSTYPE = $(shelf uname -s)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

.KEEP_STATE:
all:	$(OBJS)

	mv *.o $(MACHTYPE)-$(OSTYPE)
