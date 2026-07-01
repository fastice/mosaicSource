#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#undef MNST
#include "math.h"
//#include "mosaicSource/common/geojsonCode.h"
#include "mosaicSource/common/common.h"
#define STR_BUFFER_SIZE 1024
#define STR_BUFF(fmt, ...) ({ \
    char* __buf = (char*)calloc(STR_BUFFER_SIZE, sizeof(char)); \
    snprintf(__buf, STR_BUFFER_SIZE, fmt, ##__VA_ARGS__); \
    __buf; \
})


int32_t RangeSize = RANGESIZE;				/* Range size of complex image */
int32_t AzimuthSize = AZIMUTHSIZE;			/* Azimuth size of complex image */
int32_t BufferSize = BUFFERSIZE;			/* Size of nonoverlap region of the buffer */
int32_t BufferLines = 512;					/* # of lines of nonoverlap in buffer */
int32_t sepAscDesc = TRUE;
int32_t HemiSphere = NORTH;
double Rotation = 45.;
double SLat = -91.;
int llConserveMem = 999;

float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
void *lBuf1, *lBuf2, *lBuf3, *lBuf4;
int main(int argc, char *argv[])
{
inputImageStructure testImage1;
inputImageStructure testImage2;
 GDALAllRegister();

//parseGeojson("geodat10x2.geojson", &testImage1);
parseInputFile("geodat10x2.in", &testImage1);
FILE *fp1 = fopen("geojsonResult", "w");
dumpParsedInput(&testImage1, fp1);

//parseInputFile("geodat.geojson", &testImage2);
//FILE *fp2 = fopen("geodatResult", "w");
//dumpParsedInput(&testImage2, fp2);
}