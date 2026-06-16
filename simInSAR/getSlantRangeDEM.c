#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"
#include <stdlib.h>
/*
   Input slant range dem.
*/
void getSlantRangeDEM(char *demFile, demStructure *dem,
                      sceneStructure scene)
{
    FILE *fp;
    int32_t i;
    /*
       Init image
    */
    fprintf(stderr, "begin getSlantRangeDEM\n");
    dem->demData = mallocImage(scene.I.azimuthSize, scene.I.rangeSize);
    dem->lonSize = scene.I.rangeSize;
    dem->latSize = scene.I.azimuthSize;
    /*
        Input image
    */
    fp = openInputFile(demFile);
    for (i = 0; i < scene.I.azimuthSize; i++)
        freadBS(dem->demData[i], scene.I.rangeSize, sizeof(float), fp, FLOAT32FLAG);
    fprintf(stderr, "end getSlantRangeDEM\n");

    return;
}
