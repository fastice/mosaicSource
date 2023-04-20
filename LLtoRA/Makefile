CC =		gcc
CFLAGS  =	$(FLAGS)  -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS= lltora.o

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))

.KEEP_STATE:
all:	$(OBJS)

	mv *.o $(MACHTYPE)-$(OSTYPE)
