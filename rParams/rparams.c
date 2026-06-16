#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "rparams.h"
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>
#include "gdalIO/gdalIO/grimpgdal.h"
//#include "mosaicSource/common/common.h"
/*
  Estimate range params using tiepoints.

  This program uses some of the routines for geocode,
  which means there is alot of unused junk to initialize everything correctly.

  Ionosphere correction:
    When estimateIonosphere.py has been run, SetupNISAR.py stamps the key
      ionosphereRangeOffsetCorrection = <basename>
    on range.offsets.vrt.  getROffsets peeks at that VRT metadata before
    calling readRangeOffsets; checkForIonosphereCorrection (inside
    readRangeOffsets) then loads the named correction file — a GeoTIFF on
    the native ROFF grid in SLC pixels — into offsets.rOffCorrection.
    The default mode (ION_AUTO) runs with and without the correction and
    emits whichever baseline solution has the lower residual sigma; the
    chosen mode is recorded in the output baseline file via the
    ;* offsetCorrectionFile line so mosaic3d can re-apply the same decision.
*/

/* ionosphereMode values */
#define ION_AUTO    0   /* default: run with and without, pick best sigma */
#define ION_NONE    1   /* -noIonosphere: never apply correction */
#define ION_FORCE   2   /* -forceIonosphere: always apply if file exists */

static void readArgs(int32_t argc, char *argv[], char **geodatFile, char **tiePointFile, char **offsetFile,
					 char **baselineFile, tiePointsStructure *tiepoints, char **shelfMaskFile, int32_t *ionosphereMode);
static void setMapProjectionForHemisphere(tiePointsStructure *tiePoints);


int32_t llConserveMem = 999; /* NO mem conserve Kluge to maintain backwards compat 9/13/06 */


static void usage();
//#define MAXOFFBUF 72000000
//#define MAXOFFLENGTH 30000
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
	FILE *tiePointFp;
	demStructure dem;
	tiePointsStructure tiePoints;
	double nDays, stdLat;
	inputImageStructure inputImage, inputImage2; /* input image */
	outputImageStructure outputImage;
	char *demFile, *geodatFile, *tiePointFile, *offsetFile, *baselineFile;
	double deltaR, deltaA;
	char *outputFile;
	char *shelfMaskFile;
	Offsets offsets;
	memset(&offsets, 0, sizeof(offsets));
	int32_t imageFlag, passType, noDEM, noRamp, dBpFlag, timeReverseFlag, ionosphereMode;
	int32_t bufferSize;
	int32_t imageCoords;
	int32_t linFlag;
	int32_t i, j; /* LCV */
	Abuf1 = NULL;
	Abuf2 = NULL;
	Dbuf1 = NULL;
	Dbuf2 = NULL;
	GDALAllRegister();
	/* Used for pointers rows to above buffer space */
	lBuf3 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	lBuf4 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	offBufSpace3 = (void *)malloc(MAXOFFBUF);
	offBufSpace4 = (void *)malloc(MAXOFFBUF);
	/*
	   Read command line args and compute filenames
	*/
	readArgs(argc, argv, &geodatFile, &tiePointFile, &offsetFile, &baselineFile, &tiePoints, &shelfMaskFile, &ionosphereMode);
	/*
	  Parse input file
	*/
	noDEM = TRUE;
	offsets.rFile = offsetFile;
	parseInputFile(geodatFile, &inputImage);
	/*
	  Input tiepoints
	*/
	tiePointFp = openInputFile(tiePointFile);
	tiePoints.motionFlag = TRUE;
	readTiePoints(tiePointFp, &tiePoints, noDEM);
	/*
	  Determine hemisphere
	*/
	setTiePointsMapProjectionForHemisphere(&tiePoints);
	/*
	  Call to set up stuff, don't really use the outputimage
	*/
	outputImage.noMem = TRUE;
	outputImage.originX = LARGEINT;
	outputImage.originY = LARGEINT;
	outputImage.xSize = -1; 
	outputImage.ySize = -1;
	initOutputImage(&outputImage, inputImage);
	/*
	  Read shelf mask if needed
	*/
	outputImage.fpLog = stderr;
	if (shelfMaskFile != NULL)
	{
		fprintf(stderr, "shelfMask %s\n", shelfMaskFile);
		readShelf(&outputImage, shelfMaskFile);
	}
	else
		outputImage.shelfMask = NULL;
	/*
	  Compute image coords and get z from dem if necessary
	*/
	computeTiePoints(&inputImage, &tiePoints, dem, noDEM, geodatFile, outputImage.shelfMask, FALSE);
	/*
	  Probe: run getBaselineFile and getROffsets with stdout suppressed to determine
	  whether ion correction is available and to set up offsets/tiePoints state.
	  computeTiePoints (above) has already written the InSAR params header to stdout;
	  flush that now so it reaches the real output before any redirect.
	*/
	fflush(stdout);
	{
		int devnull_fd = open("/dev/null", O_WRONLY);
		int probe_saved = dup(STDOUT_FILENO);
		dup2(devnull_fd, STDOUT_FILENO);
		close(devnull_fd);

		//getBaselineFile(baselineFile, &tiePoints, inputImage);
		fprintf(stderr, "------- %s\n", offsets.rFile);
		getROffsets(offsetFile, &tiePoints, inputImage, &offsets, ionosphereMode == ION_NONE);
		fprintf(stderr, "OFFSETS READ\n\n");

		dup2(probe_saved, STDOUT_FILENO);
		close(probe_saved);
	}
	fprintf(stderr, "%s %s %i \n", offsets.geo1, offsets.geo2, (int)tiePoints.deltaB);
	if (offsets.geo1 != NULL && offsets.geo2 != NULL &&
		( (tiePoints.deltaB != DELTABNONE) || (tiePoints.initWithSV == TRUE)))
	{
		parseInputFile(offsets.geo2, &inputImage2);
		fprintf(stderr, "inputImage2 %s\n", offsets.geo2);
		//getBaselineFile(baselineFile, &tiePoints, inputImage);
		initllToImageNew(&inputImage2);
		memcpy(&(offsets.sv2), &(inputImage2.sv), sizeof(inputImage2.sv));
		offsets.dt1t2 = inputImage.cpAll.sTime - inputImage2.cpAll.sTime;
		fprintf(stderr, "times %f %f %f\n", inputImage.cpAll.sTime, inputImage2.cpAll.sTime, offsets.dt1t2);
		svOffsets(&inputImage, &inputImage2, &offsets, &(tiePoints.cnstR), &(tiePoints.cnstA));
	} else if(tiePoints.deltaB != DELTABNONE)
		error("SV baselines but geodats not specified in .dat file ");
	if (tiePoints.deltaB == DELTABNONE)
	{
		tiePoints.cnstA = 0.0;
		tiePoints.cnstR = 0.0;
	}
	fprintf(stderr, "Rg/Az offsets %10.5f %10.5f\n", tiePoints.cnstR, tiePoints.cnstA);

	/*
	  Estimate baseline solution.
	  For each branch: re-run getBaselineFile+getROffsets (so their output is included in
	  the final file alongside the computeRParams solution), then addVelCorrections, then solve.
	  For ION_AUTO: capture each run to a tmpfile and emit the lower-sigma result.
	*/
//fprintf(stderr, "ionosphereMode %i %s\n", ionosphereMode, offsets.rOffCorrection.correctionFile);
//error("h");
	if (offsets.rOffCorrection.rangeOffsetCorrection != NULL && ionosphereMode == ION_FORCE)
	{
		/* Single run with ion correction */
		getBaselineFile(baselineFile, &tiePoints, inputImage);
		getROffsets(offsetFile, &tiePoints, inputImage, &offsets, FALSE);
		addVelCorrections(&inputImage, &tiePoints);
		computeRParams(&tiePoints, inputImage, baselineFile, &offsets);
	}
	else if (offsets.rOffCorrection.rangeOffsetCorrection != NULL && ionosphereMode == ION_AUTO)
	{
		/* Run with and without ion correction; emit the lower-sigma result */
		char tmp1[] = "/tmp/rparams_ion_XXXXXX";
		char tmp2[] = "/tmp/rparams_noion_XXXXXX";
		int fd1 = mkstemp(tmp1);
		int fd2 = mkstemp(tmp2);
		if (fd1 < 0 || fd2 < 0) error("rparams: mkstemp failed\n");
		int saved_stdout = dup(STDOUT_FILENO);

		/* --- Run 1: with ionosphere correction --- */
		dup2(fd1, STDOUT_FILENO);
		getBaselineFile(baselineFile, &tiePoints, inputImage);
		getROffsets(offsetFile, &tiePoints, inputImage, &offsets, FALSE);
		addVelCorrections(&inputImage, &tiePoints);
		double sigma_ion = computeRParams(&tiePoints, inputImage, baselineFile, &offsets);
		fflush(stdout);

		/* --- Run 2: without ionosphere correction --- */
		dup2(fd2, STDOUT_FILENO);
		getBaselineFile(baselineFile, &tiePoints, inputImage);
		getROffsets(offsetFile, &tiePoints, inputImage, &offsets, TRUE);
		addVelCorrections(&inputImage, &tiePoints);
		/* clear correctionFile so it is not written to the no-ion output */
		offsets.rOffCorrection.correctionFile[0] = '\0';
		double sigma_noion = computeRParams(&tiePoints, inputImage, baselineFile, &offsets);
		fflush(stdout);

		/* --- Restore stdout and emit the winner --- */
		dup2(saved_stdout, STDOUT_FILENO);
		close(saved_stdout);
		close(fd1);
		close(fd2);

		int use_ion = (sigma_ion <= sigma_noion);
		fprintf(stderr, "sigma with ion correction: %f  without: %f  -- using %s\n",
		        sigma_ion, sigma_noion, use_ion ? "with ion" : "without ion");

		FILE *winner = fopen(use_ion ? tmp1 : tmp2, "r");
		char buf[4096];
		size_t n;
		while ((n = fread(buf, 1, sizeof(buf), winner)) > 0)
			fwrite(buf, 1, n, stdout);
		fclose(winner);
		fprintf(stdout, "; sigma with ion correction: %f  without: %f  -- using %s\n",
		        sigma_ion, sigma_noion, use_ion ? "with ion" : "without ion");
		unlink(tmp1);
		unlink(tmp2);
	}
	else
	{
		/* ION_NONE, or no correction file available */
		getBaselineFile(baselineFile, &tiePoints, inputImage);
		getROffsets(offsetFile, &tiePoints, inputImage, &offsets, ionosphereMode == ION_NONE);
		addVelCorrections(&inputImage, &tiePoints);
		computeRParams(&tiePoints, inputImage, baselineFile, &offsets);
	}
}

static void usage()
{
	fprintf(stderr,
		"\nCompute parameters to calibrate range offsets\n"
		"Usage:\n"
		"  rparams [options] geodatFile tiepointsFile offsetFile baselineFile\n"
		"\nOptions:\n"
		"  -nDays <days>      Temporal baseline in days (default: 24)\n"
		"  -shelfMask <file>  Shelf mask file for tidal corrections\n"
		"  -constOnly         Estimate only the constant term\n"
		"  -deltaBQ           Estimate quadratic correction to state vector baseline\n"
		"  -deltaBC           Estimate constant correction to Bp component of baseline\n"
		"  -quadB             Estimate quadratic baseline terms\n"
		"  -bnbpOnly          Estimate only bn and bp\n"
		"  -bnbpdBpOnly       Estimate only bn, bp, and dBp\n"
		"  -bpdBpOnly         Estimate only bp and dBp\n"
		"  -noIonosphere      Do not apply ionosphere correction if one exists\n"
		"  -forceIonosphere   Always apply ionosphere correction if file exists\n"
		"  -quiet             Don't echo tiepoints to solution\n"
		"\nPositional arguments (required, in order):\n"
		"  geodatFile         Geodat parameter file\n"
		"  tiepointsFile      Tiepoint location file (lat,lon,z,vx,vy,vz)\n"
		"  offsetFile         Range offset file (offsetFile.dat must also exist)\n"
		"  baselineFile       CW state vector baseline file\n");
	exit(1);
}

static void readArgs(int32_t argc, char *argv[], char **geodatFile, char **tiePointFile,
					 char **offsetFile, char **baselineFile, tiePointsStructure *tiePoints,
					 char **shelfMaskFile, int32_t *ionosphereMode)
{
	int32_t bnbpFlag = FALSE, bpFlag = FALSE, bnbpdBpFlag = FALSE, bpdBpFlag = FALSE;
	int32_t constOnlyFlag = FALSE, quadB = FALSE, deltaB = DELTABNONE;
	double nDays = 24;
	int32_t i;
	*ionosphereMode = ION_AUTO;

	*shelfMaskFile = NULL;
	tiePoints->quiet = FALSE;

	if (argc < 5)
		usage();

	for (i = 1; i < argc - 4; i++)
	{
		if (strcmp(argv[i], "-nDays") == 0)
		{
			if (++i >= argc - 4) usage();
			nDays = atof(argv[i]);
		}
		else if (strcmp(argv[i], "-shelfMask") == 0 || strcmp(argv[i], "-shelfMaskFile") == 0)
		{
			if (++i >= argc - 4) usage();
			*shelfMaskFile = argv[i];
		}
		else if (strcmp(argv[i], "-bnbpOnly") == 0)
		{
			if (bpdBpFlag || bnbpdBpFlag)
				error("bnbpOnly, bpdBpOnly, bnbpdBpOnly are mutually exclusive");
			bnbpFlag = TRUE;
		}
		else if (strcmp(argv[i], "-bnbpdBpOnly") == 0)
		{
			if (bpdBpFlag || bnbpFlag)
				error("bnbpOnly, bpdBpOnly, bnbpdBpOnly are mutually exclusive");
			bnbpdBpFlag = TRUE;
		}
		else if (strcmp(argv[i], "-bpdBpOnly") == 0)
		{
			if (bnbpFlag || bnbpdBpFlag)
				error("bnbpOnly, bpdBpOnly, bnbpdBpOnly are mutually exclusive");
			bpdBpFlag = TRUE;
		}
		else if (strcmp(argv[i], "-quadB") == 0)
			quadB = TRUE;
		else if (strcmp(argv[i], "-constOnly") == 0)
			constOnlyFlag = TRUE;
		else if (strcmp(argv[i], "-deltaBQ") == 0)
			deltaB = DELTABQUAD;
		else if (strcmp(argv[i], "-deltaBC") == 0)
			deltaB = DELTABCONST;
		else if (strcmp(argv[i], "-noIonosphere") == 0)
			*ionosphereMode = ION_NONE;
		else if (strcmp(argv[i], "-forceIonosphere") == 0)
			*ionosphereMode = ION_FORCE;
		else if (strcmp(argv[i], "-quiet") == 0)
			tiePoints->quiet = TRUE;
		else
		{
			fprintf(stderr, "Unknown option: %s\n", argv[i]);
			usage();
		}
	}

	*geodatFile   = argv[argc - 4];
	*tiePointFile = argv[argc - 3];
	*offsetFile   = argv[argc - 2];
	*baselineFile = argv[argc - 1];

	tiePoints->constOnlyFlag = constOnlyFlag;
	tiePoints->linFlag       = TRUE;
	tiePoints->quadB         = quadB;
	tiePoints->dBpFlag       = TRUE;
	tiePoints->bnbpFlag      = bnbpFlag;
	tiePoints->bpFlag        = bpFlag;
	tiePoints->bnbpdBpFlag   = bnbpdBpFlag;
	tiePoints->bpdBpFlag     = bpdBpFlag;
	tiePoints->nDays         = nDays;
	tiePoints->vrFlag        = FALSE;
	tiePoints->deltaB        = deltaB;

	if (tiePoints->constOnlyFlag)
		fprintf(stderr, "\n(****Constant only fit*****\n");
	if (tiePoints->linFlag)
		fprintf(stderr, "\n(****Including linear term fit*****\n");
}
