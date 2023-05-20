#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <gdal.h>
#include <gdal_vrt.h>
#include <cpl_conv.h>
#include <stdint.h>
#include "mosaicSource/common/common.h"
#include "gdalIO/gdalIO/grimpgdal.h"

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
double RangePixelSize = RANGEPIXELSIZE;		/* Range PixelSize */
double AzimuthPixelSize = AZIMUTHPIXELSIZE; /* Azimuth PixelSize */
int32_t sepAscDesc = TRUE;
int32_t HemiSphere = NORTH;
double Rotation = 45.;
double SLat = -91.;
int llConserveMem = 999;

float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
void *lBuf1, *lBuf2, *lBuf3, *lBuf4;

void getArgs(int argc, char *argv[], char *filename[], char **datfile, char **vrtfile, int32_t *xSize, int32_t *ySize, int32_t *byteSwap, int32_t *nBands)
{
    int i;
    int c;
    int option_index = 0;
    *byteSwap = FALSE; *xSize = -1; *ySize = -1;
    while (1)
    {   
        static struct option long_options[] = {
            {"xSize", required_argument, 0, 'x'},
            {"ySize", required_argument, 0, 'y'},
            {"byteSwap", no_argument, 0, 'b'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}};

        c = getopt_long(argc, argv, "xyb:h", long_options, &option_index);
        if (c == -1)
            break;
    fprintf(stderr,"%i %c %i %i\n", c,c, option_index, optind);
        switch (c)
        {
        case 0:
            if (long_options[option_index].flag != 0)
                break;
            printf("option %s", long_options[option_index].name);
            if (optarg)
                printf(" with arg %s", optarg);
            printf("\n");
            break;
        case 'x':
            *xSize = atoi(optarg);
            break;
        case 'y':
            *ySize = atoi(optarg);
            break;
        case 'b':
            *byteSwap = TRUE;
            break;
        case 'h':
            printf("Usage: %s [--xSize] xSize [--ySize] ySize [--byteSwap] [--help] offsetFile\n", argv[0]);
            exit(EXIT_SUCCESS);
        case '?':
            break;
        default:
            printf("Unknown option %c\n", c);
            exit(EXIT_FAILURE);
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "Expected non-optional argument after options\n");
        exit(EXIT_FAILURE);
    }
    *nBands = 0;
    for(i=optind; i < argc; i++) {
        filename[*nBands] = argv[i];
        fprintf(stderr, "arg %i %s\n", *nBands, filename[*nBands]);
        (*nBands)++;
    }
    *datfile = appendSuffix(*filename, ".dat", malloc(strlen(filename[0])+5));
    *vrtfile = appendSuffix(*filename, ".vrt", malloc(strlen(filename[0])+5));
}

int main(int argc, char **argv) {
    dictNode *metaData = NULL;
    int32_t xSize, ySize, byteSwap, i;
    char *filenames[3], *datfile, *vrtfile;
    int nBands;
    Offsets offsets;
    double geoTransform[6] = {-0.5, 1., 0., -0.5, 0., 1.};
    GDALAllRegister();
    getArgs(argc, argv, filenames, &datfile, &vrtfile, &xSize, &ySize, &byteSwap, &nBands);
    fprintf(stderr, "nBands %i\n", nBands);
    for(i=0; i < nBands; i++) { 
        fprintf(stderr,"filename %s xSize %i ySize %i byteSwap %i\n", filenames[i], xSize, ySize, byteSwap);
    }
    readOffsetParams(datfile, &offsets, FALSE);
    insert_node(&metaData, "r0", STR_BUFF("%i", offsets.rO));
	insert_node(&metaData, "a0", STR_BUFF("%i", offsets.aO));
    //insert_node(&metaData, "nr", STR_BUFF("%i", offsets.nr));
	//insert_node(&metaData, "na", STR_BUFF("%i", offsets.na));
	insert_node(&metaData, "deltaA", STR_BUFF("%f", offsets.deltaA));
	insert_node(&metaData, "deltaR", STR_BUFF("%f", offsets.deltaR));
	insert_node(&metaData, "sigmaStreaks", STR_BUFF("%f", offsets.sigmaStreaks));
	insert_node(&metaData, "sigmaRange", STR_BUFF("%f", offsets.sigmaRange));
    insert_node(&metaData, "geo1", STR_BUFF("%s", offsets.geo1));
	insert_node(&metaData, "geo2", STR_BUFF("%s", offsets.geo2));
    printDictionary(metaData);
	/*
	   read geodat files if the exist
	 */
  
    makeVRT(vrtfile, offsets.nr, offsets.na, GDT_Float32, filenames, nBands, geoTransform, byteSwap, metaData);
}