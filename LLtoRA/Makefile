CC =		gcc
CFLAGS  =	$(FLAGS)  -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS= lltora.o

MACHTYPE = $(shell uname -m)
OSTYPE = $(shelf uname -s)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

.KEEP_STATE:
all:	$(OBJS)

	mv *.o $(MACHTYPE)-$(OSTYPE)
