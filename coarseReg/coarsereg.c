#include "stdio.h"
#include"string.h"
#include "mosaicSource/common/common.h"
#include <sys/types.h>
#include <sys/types.h>
#include <sys/time.h>
#include <math.h>

/*
  This program performs registration using corner point locations
  in geodat files.

  There are two output options. The default is the range varying shifts
  necessary to register two images. The range width and pixel size
  are read from the geodat file.
  
  The -intial option only outputs the average offset. This is 
  used for intial registration prior to the more accurate registration
  via cross correlation. The outputs are the fractional range/azimuth
  pixel offsets. Although the values can be fractional, they yield
  integer single-look pixel values when multiplied by the number of looks.

*/


static void readArgs(int argc,char *argv[], char **inputFile1,char **inputFile2);
static void usage();

/* 
   Global variable definitions
*/
int RangeSize=RANGESIZE;       /* Range size of complex image */
int AzimuthSize=AZIMUTHSIZE;   /* Azimuth size of complex image */
int BufferSize=BUFFERSIZE;     /* Size of nonoverlap region of the buffer */
int BufferLines = 512;         /* # of lines of nonoverlap in buffer */
double RangePixelSize=RANGEPIXELSIZE; /* Range PixelSize */
double AzimuthPixelSize=AZIMUTHPIXELSIZE; /* Azimuth PixelSize */
int HemiSphere=NORTH;
double Rotation=45.;
double SLat=-91.0;
char *Abuf1,*Abuf2,*Dbuf1,*Dbuf2;
int llConserveMem=999; /* Kluge to maintain backwards compat 9/13/06 */

float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
void *offBufSpace1,*offBufSpace2,*offBufSpace3,*offBufSpace4;
void *lBuf1,*lBuf2,*lBuf3,*lBuf4;

void coarseRegister(inputImageStructure *inputImage1,inputImageStructure *inputImage2);


void main(int argc, char *argv[])
{   
	inputImageStructure inputImage;       /* input image */
	char *outputFile,*inputFile1, *inputFile2;
	inputImageStructure inputImage1,inputImage2;
	int  passType;
	int i,j;                               /* LCV */
	/* 
	   Read command line args and compute filenames
	*/
	readArgs(argc,argv,&inputFile1,&inputFile2);
	/*
	  Parse input file and input data for each image.
	*/  
	inputImage1.passType =DESCENDING; /* default */
	inputImage1.stateFlag = TRUE;
	parseInputFile(inputFile1, &(inputImage1));

	inputImage2.passType =DESCENDING; /* default */
	inputImage2.stateFlag = TRUE;
	parseInputFile(inputFile2, &(inputImage2));
	/*
	  Compute Offsets
	*/
	coarseRegister(&inputImage1,&inputImage2);
	return; 
}


static void readArgs(int argc,char *argv[],  char **inputFile1,char **inputFile2)
{
	int filenameArg;
	char *argString;
	int i,n;
	if(argc != 3) usage();
	*inputFile1 = argv[argc-2];
	*inputFile2 = argv[argc-1];
	return;
}

 
static void usage()
{ 
	error("\n%s\n\n%s\n\n%s\n",
	      "Compute initial (coarse) registration offsets, output to stdout",
	      "Usage:", 
	      " coarsereg geodat1 geodat2\n");
}
