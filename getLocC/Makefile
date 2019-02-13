CC =		gcc
CFLAGS  =	$(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH)

OBJS=  centerLL.o correctTime.o  getlocc.o  glatlon.o   readOldPar.o

.KEEP_STATE:
all:	$(OBJS)

	mv *.o $(MACHTYPE)-$(OSTYPE)
