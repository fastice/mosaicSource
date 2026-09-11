CC =		gcc
MACHTYPE =	$(shell uname -m)
OSTYPE =	$(shell uname -s)
OBJDIR =	$(MACHTYPE)-$(OSTYPE)

$(shell mkdir -p $(OBJDIR))

CFLAGS  =	$(FLAGS) -c -DPAF=$(PAF) -I$(INCLUDEPATH) -MMD -MP -fopenmp

SRCS =	get3DInputFile.c make3DOffsets.c make3DOffsetsJoint.c make3DMosaic.c make3DMosaicJoint.c makeVhMosaic.c \
	mosaic3d.c mosaicHopper.c mosaicHopper3D.c mosaicTrue3D.c mosaicTrue3DPhase.c obsDump.c setup3D.c speckleTrackMosaic.c writeTieFile.c

OBJS =	$(SRCS:%.c=$(OBJDIR)/%.o)

all:	$(OBJS)

-include $(OBJS:.o=.d)

.KEEP_STATE:


$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $< -o $@
