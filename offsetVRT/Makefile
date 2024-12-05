CC =		gcc
CFLAGS =	$(FLAGS) -c -I$(INCLUDEPATH)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

OBJS =	offsetVRT.o

MACHTYPE = $(shell uname -m)
OSTYPE = $(shelf uname -s)

all:	$(OBJS)
	mv *.o $(MACHTYPE)-$(OSTYPE)
.KEEP_STATE:
