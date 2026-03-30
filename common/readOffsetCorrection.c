#include "stdio.h"
#include "string.h"
#include "math.h"
#include "stdlib.h"
#include "mosaicSource/common/common.h"

// Malloc the space if hasn't been malloced yet. This is called from readOffsets after the offsets 
// have been read and the pass type is known. So this will only be called if there are offsets to read,
// and it will only be called once per pass type (e.g., if there are only ASC offsets, then it will only be called once for ASC).

void mallocCorrectionBuffer(int bufferMode)
{
    extern void *correctionBuf1, *correctionBuf2, *Correction1, *Correction2;
    
    if(bufferMode == RANGEBUFF && correctionBuf1 == NULL)
    {
        correctionBuf1 = malloc(sizeof(float *) * MAXOFFLENGTH);
        Correction1 = malloc(MAXOFFBUF);
        return;
    } 

    if(bufferMode == RANGEUSEAZIMUTHBUFF && correctionBuf2 == NULL)
    {
        correctionBuf2 = malloc(sizeof(float *) * MAXOFFLENGTH);
        Correction2 = malloc(MAXOFFBUF);
        return;
    } 
    return;
}

void readOffsetCorrection(char *correctionFile, Offsets *offsets, int bufferMode)
{
    // On first use malloc the correction buffers, otherwise will do nothing.
    return;
    /*
    mallocCorrectionBuffer(bufferMode);

    if(bufferMode == RANGEBUFF)
    {

        mapBuffer(nCols, nLines, (float **)Correction, NULL, (float *)ACorrectionBuf, NULL);
        return;
    } 

    if(bufferMode == AZIMUTHBUFF && DCorrectionBuf != NULL)
    {
        mapBuffer(nCols, nLines, (float **)DCorrection, NULL, (float *)DCorrectionBuf, NULL);
        return;
    } */
    return;
}