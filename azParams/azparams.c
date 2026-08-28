#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "mosaicSource/common/common.h"
#include "azparams.h"
#include <sys/types.h>
#include <sys/time.h>
#include "math.h"
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
//#include "gdalIO/gdalIO/grimpgdal.h"
/*
  Estimate baseline using tiepoints.

  This program uses some of the routines for geocode,
  which means there is alot of unused junk to initialize everything correctly.
*/

static void readArgs(int32_t argc, char *argv[], char **geodatFile, char **tiePointFile, char **offsetFile, char **baselineFile, tiePointsStructure *tiePoints, int32_t *yamlOutput, int32_t *debugFlag, char **outputFile, char **runFile);
static void usage();
static const char *azparamsModeString(tiePointsStructure *tiePoints);

/* ---------------------------------------------------------------
   Multi-run (-runFile) support -- mirrors rParams/rparams.c. Lets the tie
   workflow run the three azimuth fits (linear / SV-const / SV-linear) from a
   single invocation, reading azimuth.offsets once instead of three times.
   --------------------------------------------------------------- */
#define AZPARAMS_MAX_RUNS 64
#define AZPARAMS_PATH_LEN 4096

/* Per-run fit mode -> tiePoints flags (mirrors the 3 old azparams calls). */
#define AZMODE_LINEAR 0   /* -linear            (az.est)          */
#define AZMODE_SVCONST 1  /* -useSV -constOnly  (az.est.const)    */
#define AZMODE_SVLINEAR 2 /* -useSV -linear     (az.est.svlinear) */

typedef struct {
	int32_t mode;
	char tiefile[AZPARAMS_PATH_LEN];
	char outfile[AZPARAMS_PATH_LEN];
	char oldfile[AZPARAMS_PATH_LEN]; /* optional old file to unlink after write */
} azRunSpec_t;

static void readAzRunFile(const char *specFile, azRunSpec_t *runs, int *nRuns)
{
	FILE *fp = fopen(specFile, "r");
	if (fp == NULL) error("azparams: cannot open -runFile %s", specFile);
	*nRuns = 0;
	char line[AZPARAMS_PATH_LEN * 5];
	while (fgets(line, sizeof(line), fp))
	{
		char *nl = strchr(line, '\n');
		if (nl) *nl = '\0';
		char *p = line;
		while (*p == ' ' || *p == '\t') p++;
		if (*p == '#' || *p == '\0') continue;
		if (*nRuns >= AZPARAMS_MAX_RUNS)
			error("azparams: runFile exceeds %d entries", AZPARAMS_MAX_RUNS);
		azRunSpec_t *r = &runs[*nRuns];
		r->oldfile[0] = '\0';
		char token[64];
		int n = sscanf(p, "%63s %4095s %4095s %4095s", token, r->tiefile, r->outfile, r->oldfile);
		if (n < 3) error("azparams: bad runFile line: %s", line);
		if      (strcmp(token, "LINEAR")   == 0) r->mode = AZMODE_LINEAR;
		else if (strcmp(token, "SVCONST")  == 0) r->mode = AZMODE_SVCONST;
		else if (strcmp(token, "SVLINEAR") == 0) r->mode = AZMODE_SVLINEAR;
		else error("azparams: unknown mode token '%s' in runFile", token);
		(*nRuns)++;
	}
	fclose(fp);
	if (*nRuns == 0) error("azparams: runFile %s has no valid entries", specFile);
}

/* Apply a run's fit mode to the tiePoints flags. */
static void setAzMode(tiePointsStructure *tiePoints, int32_t mode)
{
	switch (mode)
	{
	case AZMODE_LINEAR:
		tiePoints->linFlag = TRUE;  tiePoints->constOnlyFlag = FALSE; tiePoints->deltaB = DELTABNONE;  break;
	case AZMODE_SVCONST:
		tiePoints->linFlag = FALSE; tiePoints->constOnlyFlag = TRUE;  tiePoints->deltaB = DELTABCONST; break;
	case AZMODE_SVLINEAR:
		tiePoints->linFlag = TRUE;  tiePoints->constOnlyFlag = FALSE; tiePoints->deltaB = DELTABCONST; break;
	}
}

/* Ignore any embedded VRT dataset mask band on offset inputs; default off
   (mask honored when present). Defined once in common/getRegion.c since
   common/readOffsets.c is linked into every program in mosaicSource/. */
extern int32_t noMask;


int32_t llConserveMem = 999; /* NO mem conserve Kluge to maintain backwards compat 9/13/06 */
int32_t sepAscDesc;

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
//float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 not use only for mosaic3d compatability */

int main(int argc, char *argv[])
{
	extern char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
	
	FILE *tiePointFp;
	demStructure dem;
	tiePointsStructure tiePoints;
	double nDays, stdLat;
	inputImageStructure inputImage; /* input image */
	outputImageStructure outputImage;
	char *demFile, *geodatFile, *tiePointFile, *offsetFile, *baselineFile;
	char *outputFile, *runFile;
	Offsets offsets;
	memset(&offsets, 0, sizeof(offsets));
	int32_t imageFlag, passType, noDEM, noRamp, dBpFlag, timeReverseFlag;
	int32_t bufferSize;
	int32_t imageCoords;
	int32_t constOnlyFlag, linFlag;
	int32_t yamlOutput, debugFlag;
	int32_t i, j; /* LCV */
	Abuf1 = NULL;
	Abuf2 = NULL;
	Dbuf1 = NULL;
	Dbuf2 = NULL;
	/* Used for pointers rows to above buffer space */
	lBuf1 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	lBuf2 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	offBufSpace1 = (void *)malloc(MAXOFFBUF);
	offBufSpace2 = (void *)malloc(MAXOFFBUF);
	GDALAllRegister();
	/*
	   Read command line args and compute filenames
	*/
	outputFile = NULL;
	runFile = NULL;
	readArgs(argc, argv, &geodatFile, &tiePointFile, &offsetFile, &baselineFile, &tiePoints, &yamlOutput, &debugFlag, &outputFile, &runFile);
	/*
	  Parse input file
	*/
	inputImage.stateFlag = TRUE;
	noDEM = TRUE;
	offsets.file = offsetFile;
	parseInputFile(geodatFile, &inputImage);
	/* ================================================================
	   MULTI-RUN PATH (-runFile): run the 3 azimuth fits from one command,
	   reading azimuth.offsets once (skipLoad). Mirrors rParams/rparams.c.
	   ================================================================ */
	if (runFile != NULL)
	{
		static azRunSpec_t runs[AZPARAMS_MAX_RUNS];
		int nRuns;
		int offsetsLoaded = 0;
		readAzRunFile(runFile, runs, &nRuns);
		outputImage.noMem = TRUE;
		outputImage.originX = LARGEINT;
		outputImage.originY = LARGEINT;
		initOutputImage(&outputImage, inputImage);
		/* Each run reproduces the single-run sequence exactly (fresh tiepoint
		   read + geocode, with the fit mode set first), so the only thing shared
		   is the azimuth.offsets read: getOffsets loads it on the first run and
		   reuses offsets->da (skipLoad) thereafter. This keeps every run
		   bit-identical to the old separate azparams invocations. */
		for (i = 0; i < nRuns; i++)
		{
			setAzMode(&tiePoints, runs[i].mode);
			tiePointFp = openInputFile(runs[i].tiefile);
			tiePoints.motionFlag = TRUE;
			readTiePoints(tiePointFp, &tiePoints, noDEM);
			fclose(tiePointFp);
			setTiePointsMapProjectionForHemisphere(&tiePoints);
			tiePoints.imageCoords = FALSE;
			/* Redirect stdout to this run's output file. With -quiet the geocode
			   and tiepoint listing are suppressed, so the file holds only the
			   computeAzParams result -- same as the old "azparams ... > az.est". */
			int out_fd = open(runs[i].outfile, O_WRONLY | O_CREAT | O_TRUNC, 0666);
			if (out_fd < 0) error("azparams: cannot open output file %s", runs[i].outfile);
			int run_saved = dup(STDOUT_FILENO);
			dup2(out_fd, STDOUT_FILENO);
			close(out_fd);
			computeTiePoints(&inputImage, &tiePoints, dem, noDEM, geodatFile, NULL, tiePoints.quiet);
			getOffsets(offsetFile, &tiePoints, inputImage, &offsets, offsetsLoaded);
			offsetsLoaded = 1;
			tiePoints.initWithSV = FALSE;
			if (offsets.geo1 != NULL && offsets.geo2 != NULL)
			{
				tiePoints.initWithSV = TRUE;
				svInitAzParams(&inputImage, &offsets);
			}
			if (tiePoints.deltaB == DELTABNONE)
			{
				tiePoints.cnstA = 0.0;
				tiePoints.cnstR = 0.0;
			}
			addOffsetCorrections(&inputImage, &tiePoints);
			char runDebugFileBuf[AZPARAMS_PATH_LEN + 32];
			char *runDebugFile = NULL;
			if (debugFlag)
			{
				snprintf(runDebugFileBuf, sizeof(runDebugFileBuf), "%s.residuals.gpkg", runs[i].outfile);
				runDebugFile = runDebugFileBuf;
			}
			computeAzParams(&tiePoints, &inputImage, baselineFile, &offsets, yamlOutput, runDebugFile);
			fflush(stdout);
			dup2(run_saved, STDOUT_FILENO);
			close(run_saved);
			if (runs[i].oldfile[0] != '\0')
				if (unlink(runs[i].oldfile) != 0 && errno != ENOENT)
					fprintf(stderr, "azparams: warning: could not remove %s: %s\n",
							runs[i].oldfile, strerror(errno));
		}
		return 0;
	}
	/*
	  Set inputs
	*/
	if (tiePoints.constOnlyFlag == TRUE)
		fprintf(stderr, "\n(****Constant only fit*****\n");
	if (tiePoints.linFlag == TRUE)
		fprintf(stderr, "\n(****Including linear term fit*****\n");
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
	tiePoints.imageCoords = FALSE;
	/*
	  Call to set up stuff, don't really use the outputimage
	*/
	outputImage.noMem = TRUE;
	outputImage.originX = LARGEINT;
	outputImage.originY = LARGEINT;
	initOutputImage(&outputImage, inputImage);
	/*
	  Compute image coords and get z from dem if necessary
	*/
	computeTiePoints(&inputImage, &tiePoints, dem, noDEM, geodatFile, NULL, tiePoints.quiet);
	/*
	  Extract phases from phase file.
	*/
	getOffsets(offsetFile, &tiePoints, inputImage, &offsets, 0);
	fprintf(stderr, "%s %s %i \n", offsets.geo1, offsets.geo2, (int)tiePoints.deltaB);
	tiePoints.initWithSV = FALSE;
	if (offsets.geo1 != NULL && offsets.geo2 != NULL)
	{
		tiePoints.initWithSV = TRUE;
		svInitAzParams(&inputImage, &offsets);
		fprintf(stderr, "sv fit %f %f %f %f\n", offsets.azFit[0], offsets.azFit[1], offsets.azFit[2], offsets.azFit[3]);
	}
	if(tiePoints.deltaB == DELTABNONE) 
	{
		tiePoints.cnstA = 0.0;
		tiePoints.cnstR = 0.0;
	}
	/*
	  Add baseline corrections for previously removed baselines
	*/
	addOffsetCorrections(&inputImage, &tiePoints);
	/*
	  Output results for checking to sterr

	*/
	if (tiePoints.quiet == FALSE)
	{
		for (i = 0; i < tiePoints.npts; i++)
			if (fabs(tiePoints.phase[i]) < 200000)
				fprintf(stderr, "%8.1f %8.1f %8.1f ---  %7.2f %7.2f --- %f ---- %f\n",
						tiePoints.x[i], tiePoints.y[i], tiePoints.z[i], tiePoints.r[i], tiePoints.a[i], tiePoints.phase[i], tiePoints.vyra[i]);
		fprintf(stderr, "\n");
	}
	/*
	  Estimate baseline solution and output to stdout (or -outputFile, if given)
	*/
	char *debugFile = NULL;
	char debugFileBuf[4096];
	if (debugFlag)
	{
		if (outputFile != NULL)
			snprintf(debugFileBuf, sizeof(debugFileBuf), "%s.residuals.gpkg", outputFile);
		else
			snprintf(debugFileBuf, sizeof(debugFileBuf), "azparams.%s.gpkg", azparamsModeString(&tiePoints));
		debugFile = debugFileBuf;
	}

	FILE *outFp = NULL;
	int savedStdoutForFile = -1;
	if (outputFile != NULL)
	{
		outFp = fopen(outputFile, "w");
		if (outFp == NULL)
			error("azparams: cannot open -outputFile %s", outputFile);
		savedStdoutForFile = dup(STDOUT_FILENO);
		dup2(fileno(outFp), STDOUT_FILENO);
	}

	computeAzParams(&tiePoints, &inputImage, baselineFile, &offsets, yamlOutput, debugFile);

	if (outFp != NULL)
	{
		fflush(stdout);
		dup2(savedStdoutForFile, STDOUT_FILENO);
		close(savedStdoutForFile);
		fclose(outFp);
	}
}

static const char *azparamsModeString(tiePointsStructure *tiePoints)
{
	if (tiePoints->constOnlyFlag == TRUE)
		return tiePoints->linFlag ? "constOnlyLinear" : "constOnly";
	return tiePoints->linFlag ? "linear" : "default";
}

static void readArgs(int argc, char *argv[], char **geodatFile, char **tiePointFile, char **offsetFile, char **baselineFile, tiePointsStructure *tiePoints, int32_t *yamlOutput, int32_t *debugFlag, char **outputFile, char **runFile)
{
	int32_t filenameArg;
	char *argString;
	double nDays = -1;
	int32_t linFlag = FALSE, constOnlyFlag = FALSE;
	int32_t deltaB = DELTABNONE;
	int32_t i, n, nPos;
	*runFile = NULL;
	/* First pass: -runFile drops the tiepoint-file positional (per-run in the
	   spec), so positionals become geodat / offsets / baseline (3 not 4). */
	for (i = 1; i < argc; i++)
		if (strcmp(argv[i], "-runFile") == 0 && i + 1 < argc)
		{
			*runFile = argv[i + 1];
			break;
		}
	nPos = (*runFile != NULL) ? 3 : 4;
	if (argc < nPos + 2 || argc > 25)
		usage(); /* Check number of args */
	n = argc - nPos - 1;
	tiePoints->quiet = FALSE;
	*yamlOutput = FALSE;
	*debugFlag = FALSE;
	*outputFile = NULL;
	for (i = 1; i <= n; i++)
	{
		argString = strchr(argv[i], '-');
		if (strstr(argString, "runFile") != NULL)
		{
			i++; /* value already captured in the first pass */
		}
		else if (strstr(argString, "constOnly") != NULL)
		{
			constOnlyFlag = TRUE;
		}
		else if (strstr(argString, "yaml") != NULL)
		{
			/* 2026-06-17: yaml output flag */
			*yamlOutput = TRUE;
		}
		else if (strstr(argString, "linear") != NULL)
		{
			linFlag = TRUE;
		}
		else if (strstr(argString, "nDays") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &nDays);
			i++;
		}
		else if (strstr(argString, "useSV") != NULL)
		{
			deltaB = DELTABCONST;
		}
		else if (strstr(argString, "quiet") != NULL)
		{
			tiePoints->quiet = TRUE;
		}
		else if (strstr(argString, "outputFile") != NULL)
		{
			*outputFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "debug") != NULL)
			*debugFlag = TRUE;
		else if (strstr(argString, "noMask") != NULL)
			noMask = TRUE;
		else if (i != n)
			usage();
	}
	if (*runFile != NULL && *outputFile != NULL)
		error("azparams: -runFile and -outputFile are mutually exclusive");
	if (nDays < 0)
		error("nDays not specified in azparams");
	if (*runFile != NULL)
	{
		*geodatFile = argv[argc - 3];
		*tiePointFile = NULL; /* per-run in the spec */
		*offsetFile = argv[argc - 2];
		*baselineFile = argv[argc - 1];
	}
	else
	{
		*geodatFile = argv[argc - 4];
		*tiePointFile = argv[argc - 3];
		*offsetFile = argv[argc - 2];
		*baselineFile = argv[argc - 1];
	}
	tiePoints->deltaB = deltaB;
	tiePoints->linFlag = linFlag;
	tiePoints->constOnlyFlag = constOnlyFlag;
	tiePoints->nDays = nDays;
	return;
}

static void usage()
{
	/* Single concatenated literal (no %s indirection) -- same reasoning as
	   tiepoints.c's usage(): avoids the format/argument-count coupling that
	   caused a real off-by-one bug there. */
	error(
		"\n\n"
		"Compute parameters to calibrate azimuth offsets\n"
		"Usage:\n"
		" azparams  -constOnly -linear -useSV -quiet -yaml -noMask -nDays nDays geodatFile tiepointsFile offsetFile "
		"baselineFile\n"
		"where\n"
		"  constOnly    = do not estimate baseline dependent terms, can be combined with linear fit\n"
		"  linear       = add linear along track fit to either constOnly or the baseline parameter solution\n"
		"  useSV     = as any of the four other options except it adds a correction after an SV determined offset removed\n"
		"  quiet          = don't echo tiepoints to solution\n"
		"  yaml           = write YAML output instead of legacy semicolon format\n"
		"  noMask       = ignore any embedded VRT dataset mask band on offsetFile; default off (mask honored when present)\n"
		"  debug        = write every tiepoint used in the fit, plus its residual, to a\n"
		"                 GeoPackage (<outputFile>.residuals.gpkg, or azparams.<mode>.gpkg\n"
		"                 if -outputFile not given)\n"
		"  outputFile   = write the solution to <path> instead of stdout; also names the\n"
		"                 -debug GeoPackage\n"
		"  nDays        = temporals basline in days (default=24)\n"
		"  geoDatFile   = geo param file\n"
		"  tiepointFile = Tiepoint location file in (lat,lon,z,vx,vy,vz)\n"
		"  offsetFile  =  azimuth offset file (offsetFile.dat must also exist)\n"
		"  baselineFile = cw state vector baseline file \n");
}
