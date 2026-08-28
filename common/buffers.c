#include "mosaicSource/common/common.h"

float *smoothBuf;
float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
float *AIonBuffer = NULL, *DIonBuffer = NULL; /* Ionospheric phase pools, parallel to the image buffers */
char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2, *SEBuf; /* Buffers for offset and azimuth parameter interpolation, and special culling cases. */
void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4, *offSEBuffSpace;
void *lBuf1, *lBuf2, *lBuf3, *lBuf4, *lSEBuf;
void *correctionBuf1=NULL, *correctionBuf2=NULL, *Correction1=NULL, *Correction2=NULL; /* Buffers for ionospheric correction */
