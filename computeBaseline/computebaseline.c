#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mosaicSource/common/common.h"
//#include "rparams.h"
#include <sys/types.h>
#include <sys/time.h>
#include "gdalIO/gdalIO/grimpgdal.h"
/*
  Estimate range params using tiepoints.

  This program uses some of the routines for geocode,
  which means there is alot of unused junk to initialize everything correctly.
*/

static void readArgs(int32_t argc, char *argv[], char **geodatFile1, char **geodatFile2, char **baselineFile);
int32_t llConserveMem = 999; /* NO mem conserve Kluge to maintain backwards compat 9/13/06 */
static void usage();

/*
   Global variables definitions
*/
int32_t RangeSize = RANGESIZE;				/* Range size of complex image */
int32_t AzimuthSize = AZIMUTHSIZE;			/* Azimuth size of complex image */
int32_t BufferSize = BUFFERSIZE;			/* Size of nonoverlap region of the buffer */
int32_t BufferLines = 512;					/* # of lines of nonoverlap in buffer */
int32_t HemiSphere = NORTH;
double Rotation = 45.;
double SLat = -91.0;
int32_t sepAscDesc = TRUE;

int main(int argc, char *argv[])
{
	extern char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
	FILE *fp;
	demStructure dem;
	tiePointsStructure tiePoints;
	double nDays, stdLat;
	inputImageStructure inputImage1, inputImage2; /* input image */
	outputImageStructure outputImage;
	char *geodatFile1, *geodatFile2, *baselineFile;
	double dt1t2, ReH, Re, thetaC, azTime1, azTime2;
	double bn1, bp1, bn2, bp2, dbp, dbn, bn, bp;
	double bTCN[3];
	/*
	   Read command line args and compute filenames
	*/
	readArgs(argc, argv, &geodatFile1, &geodatFile2, &baselineFile);
	GDALAllRegister();
	/*
	  Parse input file
	*/
	
	parseInputFile(geodatFile1, &inputImage1);
	initllToImageNew(&inputImage1);
	parseInputFile(geodatFile2, &inputImage2);
	initllToImageNew(&inputImage2);
	// Compute theta
	Re = inputImage1.cpAll.Re;
	ReH = getReH(&(inputImage1.cpAll), &inputImage1, (inputImage1.azimuthSize) / 2);
	// Theta C uses H from geodat, which is consistent with that used to define baseline
	thetaC = thetaRReZReH(inputImage1.cpAll.RCenter, (Re + 0), (ReH));
	// Time offset bdtween images
	dt1t2 = inputImage1.cpAll.sTime - inputImage2.cpAll.sTime;
	// First and last azTtimes
	azTime1 = inputImage1.cpAll.sTime + 0 * inputImage1.nAzimuthLooks / inputImage1.par.prf;
	azTime2 = inputImage1.cpAll.sTime + inputImage1.azimuthSize * inputImage1.nAzimuthLooks / inputImage1.par.prf;
	// Compute baselines in TCN and covert to bnorm, bperp
	svBaseTCN(azTime1, dt1t2, &(inputImage1.sv), &(inputImage2.sv), bTCN);
	svBnBp(azTime1, thetaC, dt1t2, &(inputImage1.sv), &(inputImage2.sv), &bn1, &bp1, inputImage1.lookDir);
	svBaseTCN(azTime2, dt1t2, &(inputImage1.sv), &(inputImage2.sv), bTCN);
	svBnBp(azTime2, thetaC, dt1t2, &(inputImage1.sv), &(inputImage2.sv), &bn2, &bp2, inputImage1.lookDir);
	
	fprintf(stderr, "--- %f %f %f %f\n",bn1, bp1, bn2, bp2);
	bn = (bn1 + bn2) * 0.5;
	bp = (bp1 + bp2) * 0.5;
	dbn = bn2 - bn1;
	dbp = bp2 - bp1;
	fprintf(stderr, "Baseline %f %f %f %f\n", bn, bp, dbn, dbp);
	fp = fopen(baselineFile, "w");
	fprintf(fp, "%f %f %f 0.0 %f\n", bn, bp, dbn, dbp);
	fclose(fp);
}

static void readArgs(int32_t argc, char *argv[], char **geodatFile1, char **geodatFile2, char **baselineFile)
{
	if (argc != 4)
		usage();
	*geodatFile1 = argv[argc - 3];
	*geodatFile2 = argv[argc - 2];
	*baselineFile = argv[argc - 1];
}

static void usage()
{
	error(
		"\n\n%s\n%s\n%s\n",
		"Compute baseline from two geodats",
		"Usage:",
		" computeBaseline geodat1 geodat2 baselinefile \n");
}
