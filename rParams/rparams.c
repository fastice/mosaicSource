#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mosaicSource/common/common.h"
#include "rparams.h"
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
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
					 char **baselineFile, tiePointsStructure *tiepoints, char **shelfMaskFile, int32_t *ionosphereMode,
					 char **runFile, int32_t *yamlOutput);
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

/* ---------------------------------------------------------------
   Multi-run (-runFile) support
   --------------------------------------------------------------- */
#define RPARAMS_MAX_RUNS 64
#define RPARAMS_PATH_LEN 4096

typedef struct {
	int32_t deltaB;
	char tiefile[RPARAMS_PATH_LEN];
	char outfile[RPARAMS_PATH_LEN];
	char oldfile[RPARAMS_PATH_LEN];   /* optional old file to unlink after write */
} runSpec_t;

static void readRunFile(const char *specFile, runSpec_t *runs, int *nRuns)
{
	FILE *fp = fopen(specFile, "r");
	if (fp == NULL) error("rparams: cannot open -runFile %s", specFile);
	*nRuns = 0;
	char line[RPARAMS_PATH_LEN * 5];
	while (fgets(line, sizeof(line), fp)) {
		char *nl = strchr(line, '\n'); if (nl) *nl = '\0';
		char *p = line; while (*p == ' ' || *p == '\t') p++;
		if (*p == '#' || *p == '\0') continue;
		if (*nRuns >= RPARAMS_MAX_RUNS)
			error("rparams: runFile exceeds %d entries", RPARAMS_MAX_RUNS);
		runSpec_t *r = &runs[*nRuns];
		r->oldfile[0] = '\0';
		char token[64];
		int n = sscanf(p, "%63s %4095s %4095s %4095s", token, r->tiefile, r->outfile, r->oldfile);
		if (n < 3)
			error("rparams: bad runFile line: %s", line);
		if      (strcmp(token, "NONE")         == 0) r->deltaB = DELTABNONE;
		else if (strcmp(token, "DELTABCONST")  == 0) r->deltaB = DELTABCONST;
		else if (strcmp(token, "DELTABQUAD")   == 0) r->deltaB = DELTABQUAD;
		else error("rparams: unknown deltaB token '%s' in runFile", token);
		(*nRuns)++;
	}
	fclose(fp);
	if (*nRuns == 0) error("rparams: runFile %s has no valid entries", specFile);
}

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
	int32_t yamlOutput = 0;
	int32_t i, j; /* LCV */
	char *runFile = NULL;
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
	readArgs(argc, argv, &geodatFile, &tiePointFile, &offsetFile, &baselineFile, &tiePoints, &shelfMaskFile, &ionosphereMode, &runFile, &yamlOutput);
	/*
	  Parse input file
	*/
	noDEM = TRUE;
	offsets.rFile = offsetFile;
	parseInputFile(geodatFile, &inputImage);
	/*
	  Call to set up outputimage (noMem — no pixel buffers allocated)
	*/
	outputImage.noMem = TRUE;
	outputImage.originX = LARGEINT;
	outputImage.originY = LARGEINT;
	outputImage.xSize = -1;
	outputImage.ySize = -1;
	initOutputImage(&outputImage, inputImage);
	outputImage.fpLog = stderr;
	if (shelfMaskFile != NULL)
	{
		fprintf(stderr, "shelfMask %s\n", shelfMaskFile);
		readShelf(&outputImage, shelfMaskFile);
	}
	else
		outputImage.shelfMask = NULL;

	/* ================================================================
	   SINGLE-RUN PATH  (legacy: no -runFile)
	   ================================================================ */
	if (runFile == NULL)
	{
		tiePointFp = openInputFile(tiePointFile);
		tiePoints.motionFlag = TRUE;
		readTiePoints(tiePointFp, &tiePoints, noDEM);
		setTiePointsMapProjectionForHemisphere(&tiePoints);
		if (yamlOutput) {
			/* suppress tiepoint header to stdout in yaml mode */
			int sup_fd = open("/dev/null", O_WRONLY);
			int sup_saved = dup(STDOUT_FILENO);
			dup2(sup_fd, STDOUT_FILENO);
			close(sup_fd);
			computeTiePoints(&inputImage, &tiePoints, dem, noDEM, geodatFile, outputImage.shelfMask, FALSE);
			fflush(stdout);
			dup2(sup_saved, STDOUT_FILENO);
			close(sup_saved);
		} else {
			computeTiePoints(&inputImage, &tiePoints, dem, noDEM, geodatFile, outputImage.shelfMask, FALSE);
		}
		fflush(stdout);
		{
			int devnull_fd = open("/dev/null", O_WRONLY);
			int probe_saved = dup(STDOUT_FILENO);
			dup2(devnull_fd, STDOUT_FILENO);
			close(devnull_fd);
			fprintf(stderr, "------- %s\n", offsets.rFile);
			getROffsets(offsetFile, &tiePoints, inputImage, &offsets, ionosphereMode == ION_NONE, 0);
			dup2(probe_saved, STDOUT_FILENO);
			close(probe_saved);
		}
		fprintf(stderr, "%s %s %i \n", offsets.geo1, offsets.geo2, (int)tiePoints.deltaB);
		if (offsets.geo1 != NULL && offsets.geo2 != NULL &&
			((tiePoints.deltaB != DELTABNONE) || (tiePoints.initWithSV == TRUE)))
		{
			parseInputFile(offsets.geo2, &inputImage2);
			fprintf(stderr, "inputImage2 %s\n", offsets.geo2);
			initllToImageNew(&inputImage2);
			memcpy(&(offsets.sv2), &(inputImage2.sv), sizeof(inputImage2.sv));
			offsets.dt1t2 = inputImage.cpAll.sTime - inputImage2.cpAll.sTime;
			fprintf(stderr, "times %f %f %f\n", inputImage.cpAll.sTime, inputImage2.cpAll.sTime, offsets.dt1t2);
			svOffsets(&inputImage, &inputImage2, &offsets, &(tiePoints.cnstR), &(tiePoints.cnstA));
		} else if (tiePoints.deltaB != DELTABNONE)
			error("SV baselines but geodats not specified in .dat file ");
		if (tiePoints.deltaB == DELTABNONE)
		{
			tiePoints.cnstA = 0.0;
			tiePoints.cnstR = 0.0;
		}
		fprintf(stderr, "Rg/Az offsets %10.5f %10.5f\n", tiePoints.cnstR, tiePoints.cnstA);

		if (offsets.rOffCorrection.rangeOffsetCorrection != NULL && ionosphereMode == ION_FORCE)
		{
			getBaselineFile(baselineFile, &tiePoints, inputImage);
			getROffsets(offsetFile, &tiePoints, inputImage, &offsets, FALSE, 0);
			addVelCorrections(&inputImage, &tiePoints);
			computeRParams(&tiePoints, inputImage, baselineFile, &offsets, yamlOutput);
		}
		else if (offsets.rOffCorrection.rangeOffsetCorrection != NULL && ionosphereMode == ION_AUTO)
		{
			char tmp1[] = "/tmp/rparams_ion_XXXXXX";
			char tmp2[] = "/tmp/rparams_noion_XXXXXX";
			int fd1 = mkstemp(tmp1);
			int fd2 = mkstemp(tmp2);
			if (fd1 < 0 || fd2 < 0) error("rparams: mkstemp failed\n");
			int saved_stdout = dup(STDOUT_FILENO);

			dup2(fd1, STDOUT_FILENO);
			getBaselineFile(baselineFile, &tiePoints, inputImage);
			getROffsets(offsetFile, &tiePoints, inputImage, &offsets, FALSE, 0);
			addVelCorrections(&inputImage, &tiePoints);
			double sigma_ion = computeRParams(&tiePoints, inputImage, baselineFile, &offsets, yamlOutput);
			fflush(stdout);

			dup2(fd2, STDOUT_FILENO);
			getBaselineFile(baselineFile, &tiePoints, inputImage);
			getROffsets(offsetFile, &tiePoints, inputImage, &offsets, TRUE, 0);
			addVelCorrections(&inputImage, &tiePoints);
			offsets.rOffCorrection.correctionFile[0] = '\0';
			double sigma_noion = computeRParams(&tiePoints, inputImage, baselineFile, &offsets, yamlOutput);
			fflush(stdout);

			dup2(saved_stdout, STDOUT_FILENO);
			close(saved_stdout);
			close(fd1);
			close(fd2);

			/* sigma<0 means that attempt found no solution (see fewPoints() in
			   computeRParams.c) -- treat it as worse than any real fit rather than
			   letting a negative number win a naive numeric comparison. */
			int use_ion = (sigma_noion < 0) ? 1 : (sigma_ion >= 0 && sigma_ion <= sigma_noion);
			fprintf(stderr, "sigma with ion correction: %f  without: %f  -- using %s\n",
			        sigma_ion, sigma_noion, use_ion ? "with ion" : "without ion");
			FILE *winner = fopen(use_ion ? tmp1 : tmp2, "r");
			char buf[4096];
			size_t n;
			while ((n = fread(buf, 1, sizeof(buf), winner)) > 0)
				fwrite(buf, 1, n, stdout);
			fclose(winner);
			if (yamlOutput)
				fprintf(stdout, "sigmaWithIonCorrection: %f\nsigmaWithoutIonCorrection: %f\nusingIon: %s\n",
				        sigma_ion, sigma_noion, use_ion ? "True" : "False");
			else
				fprintf(stdout, "; sigma with ion correction: %f  without: %f  -- using %s\n",
				        sigma_ion, sigma_noion, use_ion ? "with ion" : "without ion");
			unlink(tmp1);
			unlink(tmp2);
		}
		else
		{
			getBaselineFile(baselineFile, &tiePoints, inputImage);
			getROffsets(offsetFile, &tiePoints, inputImage, &offsets, ionosphereMode == ION_NONE, 0);
			addVelCorrections(&inputImage, &tiePoints);
			computeRParams(&tiePoints, inputImage, baselineFile, &offsets, yamlOutput);
		}
		return 0;
	}

	/* ================================================================
	   MULTI-RUN PATH  (-runFile)
	   ================================================================ */
	static runSpec_t runs[RPARAMS_MAX_RUNS];
	int nRuns;
	readRunFile(runFile, runs, &nRuns);

	/*
	  Load first tiefile so computeTiePoints has real tiepoints, then capture
	  the header it writes to stdout (image params + tide line) into a buffer
	  so we can prepend it to every per-run output file.
	*/
	tiePointFp = openInputFile(runs[0].tiefile);
	tiePoints.motionFlag = TRUE;
	readTiePoints(tiePointFp, &tiePoints, noDEM);
	fclose(tiePointFp);
	setTiePointsMapProjectionForHemisphere(&tiePoints);

	char hdr_tmp[] = "/tmp/rparams_hdr_XXXXXX";
	int hdr_fd = mkstemp(hdr_tmp);
	if (hdr_fd < 0) error("rparams: mkstemp for header failed");
	{
		int pre_hdr = dup(STDOUT_FILENO);
		dup2(hdr_fd, STDOUT_FILENO);
		close(hdr_fd);
		computeTiePoints(&inputImage, &tiePoints, dem, noDEM, geodatFile, outputImage.shelfMask, FALSE);
		fflush(stdout);
		dup2(pre_hdr, STDOUT_FILENO);
		close(pre_hdr);
	}
	FILE *hfp = fopen(hdr_tmp, "r");
	fseek(hfp, 0, SEEK_END);
	long hdrLen = ftell(hfp);
	rewind(hfp);
	char *hdrBuf = (char *)malloc(hdrLen + 1);
	if (fread(hdrBuf, 1, hdrLen, hfp) != (size_t)hdrLen)
		error("rparams: failed to read header tmpfile");
	fclose(hfp);
	unlink(hdr_tmp);

	/* Probe: load offsets->dr from disk and detect ionosphere correction file */
	fflush(stdout);
	{
		int devnull_fd = open("/dev/null", O_WRONLY);
		int probe_saved = dup(STDOUT_FILENO);
		dup2(devnull_fd, STDOUT_FILENO);
		close(devnull_fd);
		fprintf(stderr, "------- %s\n", offsets.rFile);
		getROffsets(offsetFile, &tiePoints, inputImage, &offsets, ionosphereMode == ION_NONE, 0);
		fflush(stdout); /* flush probe output to /dev/null before restoring */
		dup2(probe_saved, STDOUT_FILENO);
		close(probe_saved);
	}

	/* SV offsets: needed by any DELTAB run; compute once with shared geometry */
	{
		int anyDeltaB = 0;
		for (i = 0; i < nRuns; i++)
			if (runs[i].deltaB != DELTABNONE) { anyDeltaB = 1; break; }
		fprintf(stderr, "%s %s %i \n", offsets.geo1, offsets.geo2, anyDeltaB);
		if (offsets.geo1 != NULL && offsets.geo2 != NULL && anyDeltaB)
		{
			parseInputFile(offsets.geo2, &inputImage2);
			fprintf(stderr, "inputImage2 %s\n", offsets.geo2);
			initllToImageNew(&inputImage2);
			memcpy(&(offsets.sv2), &(inputImage2.sv), sizeof(inputImage2.sv));
			offsets.dt1t2 = inputImage.cpAll.sTime - inputImage2.cpAll.sTime;
			fprintf(stderr, "times %f %f %f\n", inputImage.cpAll.sTime, inputImage2.cpAll.sTime, offsets.dt1t2);
			svOffsets(&inputImage, &inputImage2, &offsets, &(tiePoints.cnstR), &(tiePoints.cnstA));
		}
		else if (anyDeltaB && (offsets.geo1 == NULL || offsets.geo2 == NULL))
			error("SV baselines but geodats not specified in .dat file");
		else
		{
			tiePoints.cnstA = 0.0;
			tiePoints.cnstR = 0.0;
		}
	}
	fprintf(stderr, "Rg/Az offsets %10.5f %10.5f\n", tiePoints.cnstR, tiePoints.cnstA);

	/* Per-run loop: reload tiepoints when tiefile changes, redirect stdout per run */
	char currentTiefile[RPARAMS_PATH_LEN];
	strncpy(currentTiefile, runs[0].tiefile, sizeof(currentTiefile) - 1);
	currentTiefile[sizeof(currentTiefile) - 1] = '\0';

	for (i = 0; i < nRuns; i++)
	{
		/* Reload tiepoints when tiefile changes */
		if (strcmp(runs[i].tiefile, currentTiefile) != 0)
		{
			strncpy(currentTiefile, runs[i].tiefile, sizeof(currentTiefile) - 1);
			currentTiefile[sizeof(currentTiefile) - 1] = '\0';
			tiePointFp = openInputFile(currentTiefile);
			tiePoints.motionFlag = TRUE;
			readTiePoints(tiePointFp, &tiePoints, noDEM);
			fclose(tiePointFp);
			setTiePointsMapProjectionForHemisphere(&tiePoints);
			/* Geocode only — header already captured, suppress stdout output */
			{
				int sup_fd = open("/dev/null", O_WRONLY);
				int sup_saved = dup(STDOUT_FILENO);
				dup2(sup_fd, STDOUT_FILENO);
				close(sup_fd);
				computeTiePoints(&inputImage, &tiePoints, dem, noDEM, geodatFile, outputImage.shelfMask, TRUE);
				dup2(sup_saved, STDOUT_FILENO);
				close(sup_saved);
			}
		}

		/* Redirect stdout to this run's output file */
		int out_fd = open(runs[i].outfile, O_WRONLY | O_CREAT | O_TRUNC, 0666);
		if (out_fd < 0) error("rparams: cannot open output file %s", runs[i].outfile);
		int run_saved = dup(STDOUT_FILENO);
		dup2(out_fd, STDOUT_FILENO);
		close(out_fd);

		/* Write the shared header (skipped in yaml mode — header is not valid yaml) */
		if (!yamlOutput) {
			fwrite(hdrBuf, 1, hdrLen, stdout);
			fflush(stdout);
		}

		/* Set deltaB for this run; save correctionFile (ION_AUTO clears it) */
		tiePoints.deltaB = runs[i].deltaB;
		char savedCorrFile[2048];
		strncpy(savedCorrFile, offsets.rOffCorrection.correctionFile, sizeof(savedCorrFile) - 1);
		savedCorrFile[sizeof(savedCorrFile) - 1] = '\0';

		/* Run estimation — skipLoad=1 reuses offsets->dr loaded by probe */
		if (offsets.rOffCorrection.rangeOffsetCorrection != NULL && ionosphereMode == ION_FORCE)
		{
			getBaselineFile(baselineFile, &tiePoints, inputImage);
			getROffsets(offsetFile, &tiePoints, inputImage, &offsets, FALSE, 1);
			addVelCorrections(&inputImage, &tiePoints);
			computeRParams(&tiePoints, inputImage, baselineFile, &offsets, yamlOutput);
		}
		else if (offsets.rOffCorrection.rangeOffsetCorrection != NULL && ionosphereMode == ION_AUTO)
		{
			char tmp1[] = "/tmp/rparams_ion_XXXXXX";
			char tmp2[] = "/tmp/rparams_noion_XXXXXX";
			int fd1 = mkstemp(tmp1);
			int fd2 = mkstemp(tmp2);
			if (fd1 < 0 || fd2 < 0) error("rparams: mkstemp failed\n");
			int saved_stdout = dup(STDOUT_FILENO);

			dup2(fd1, STDOUT_FILENO);
			getBaselineFile(baselineFile, &tiePoints, inputImage);
			getROffsets(offsetFile, &tiePoints, inputImage, &offsets, FALSE, 1);
			addVelCorrections(&inputImage, &tiePoints);
			double sigma_ion = computeRParams(&tiePoints, inputImage, baselineFile, &offsets, yamlOutput);
			fflush(stdout);

			dup2(fd2, STDOUT_FILENO);
			getBaselineFile(baselineFile, &tiePoints, inputImage);
			getROffsets(offsetFile, &tiePoints, inputImage, &offsets, TRUE, 1);
			addVelCorrections(&inputImage, &tiePoints);
			offsets.rOffCorrection.correctionFile[0] = '\0';
			double sigma_noion = computeRParams(&tiePoints, inputImage, baselineFile, &offsets, yamlOutput);
			fflush(stdout);

			dup2(saved_stdout, STDOUT_FILENO);
			close(saved_stdout);
			close(fd1);
			close(fd2);

			/* sigma<0 means that attempt found no solution (see fewPoints() in
			   computeRParams.c) -- treat it as worse than any real fit rather than
			   letting a negative number win a naive numeric comparison. */
			int use_ion = (sigma_noion < 0) ? 1 : (sigma_ion >= 0 && sigma_ion <= sigma_noion);
			fprintf(stderr, "sigma with ion correction: %f  without: %f  -- using %s\n",
			        sigma_ion, sigma_noion, use_ion ? "with ion" : "without ion");
			FILE *winner = fopen(use_ion ? tmp1 : tmp2, "r");
			char buf[4096];
			size_t n;
			while ((n = fread(buf, 1, sizeof(buf), winner)) > 0)
				fwrite(buf, 1, n, stdout);
			fclose(winner);
			if (yamlOutput)
				fprintf(stdout, "sigmaWithIonCorrection: %f\nsigmaWithoutIonCorrection: %f\nusingIon: %s\n",
				        sigma_ion, sigma_noion, use_ion ? "True" : "False");
			else
				fprintf(stdout, "; sigma with ion correction: %f  without: %f  -- using %s\n",
				        sigma_ion, sigma_noion, use_ion ? "with ion" : "without ion");
			unlink(tmp1);
			unlink(tmp2);
		}
		else
		{
			getBaselineFile(baselineFile, &tiePoints, inputImage);
			getROffsets(offsetFile, &tiePoints, inputImage, &offsets, ionosphereMode == ION_NONE, 1);
			addVelCorrections(&inputImage, &tiePoints);
			computeRParams(&tiePoints, inputImage, baselineFile, &offsets, yamlOutput);
		}

		/* Restore correctionFile (may have been cleared by ION_AUTO no-ion branch) */
		strncpy(offsets.rOffCorrection.correctionFile, savedCorrFile,
		        sizeof(offsets.rOffCorrection.correctionFile) - 1);
		offsets.rOffCorrection.correctionFile[sizeof(offsets.rOffCorrection.correctionFile) - 1] = '\0';

		fflush(stdout);
		dup2(run_saved, STDOUT_FILENO);
		close(run_saved);

		/* Remove old non-yaml file now that yaml version is written */
		if (runs[i].oldfile[0] != '\0') {
			if (unlink(runs[i].oldfile) != 0 && errno != ENOENT)
				fprintf(stderr, "rparams: warning: could not remove %s: %s\n",
				        runs[i].oldfile, strerror(errno));
		}
	}
	free(hdrBuf);
	return 0;
}

static void usage()
{
	fprintf(stderr,
		"\nCompute parameters to calibrate range offsets\n"
		"Usage (single run):\n"
		"  rparams [options] geodatFile tiepointsFile offsetFile baselineFile\n"
		"Usage (multi-run):\n"
		"  rparams [options] -runFile specFile geodatFile offsetFile baselineFile\n"
		"\nOptions:\n"
		"  -nDays <days>      Temporal baseline in days (default: 24)\n"
		"  -shelfMask <file>  Shelf mask file for tidal corrections\n"
		"  -runFile <file>    Run-spec file for multi-run mode (replaces tiepointsFile)\n"
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
		"  -yaml              Write baseline output in YAML format\n"
		"\nPositional arguments (single-run):\n"
		"  geodatFile         Geodat parameter file\n"
		"  tiepointsFile      Tiepoint location file (lat,lon,z,vx,vy,vz)\n"
		"  offsetFile         Range offset file (offsetFile.dat must also exist)\n"
		"  baselineFile       CW state vector baseline file\n"
		"\nPositional arguments (multi-run, -runFile given):\n"
		"  geodatFile         Geodat parameter file\n"
		"  offsetFile         Range offset file\n"
		"  baselineFile       CW state vector baseline file\n"
		"\nRun-spec file format (one run per line):\n"
		"  NONE|DELTABCONST|DELTABQUAD  /abs/path/tiefile  outfile  [oldfile]\n"
		"  oldfile: optional old (non-yaml) file to remove after writing outfile\n");
	exit(1);
}

static void readArgs(int32_t argc, char *argv[], char **geodatFile, char **tiePointFile,
					 char **offsetFile, char **baselineFile, tiePointsStructure *tiePoints,
					 char **shelfMaskFile, int32_t *ionosphereMode, char **runFile, int32_t *yamlOutput)
{
	int32_t bnbpFlag = FALSE, bpFlag = FALSE, bnbpdBpFlag = FALSE, bpdBpFlag = FALSE;
	int32_t constOnlyFlag = FALSE, quadB = FALSE, deltaB = DELTABNONE;
	double nDays = 24;
	int32_t i, nPos;
	*ionosphereMode = ION_AUTO;
	*shelfMaskFile = NULL;
	*runFile = NULL;
	*yamlOutput = 0;
	tiePoints->quiet = FALSE;

	/* First pass: detect -runFile so we know how many positional args to expect */
	for (i = 1; i < argc - 1; i++)
	{
		if (strcmp(argv[i], "-runFile") == 0)
		{
			*runFile = argv[i + 1];
			break;
		}
	}
	/* multi-run: 3 positionals (geodat offsetFile baselineFile)
	   single-run: 4 positionals (geodat tiefile offsetFile baselineFile) */
	nPos = (*runFile != NULL) ? 3 : 4;
	if (argc < nPos + 1)
		usage();

	for (i = 1; i < argc - nPos; i++)
	{
		if (strcmp(argv[i], "-nDays") == 0)
		{
			if (++i >= argc - nPos) usage();
			nDays = atof(argv[i]);
		}
		else if (strcmp(argv[i], "-shelfMask") == 0 || strcmp(argv[i], "-shelfMaskFile") == 0)
		{
			if (++i >= argc - nPos) usage();
			*shelfMaskFile = argv[i];
		}
		else if (strcmp(argv[i], "-runFile") == 0)
		{
			if (++i >= argc - nPos) usage();
			/* already captured in first pass */
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
		else if (strcmp(argv[i], "-yaml") == 0)
			*yamlOutput = 1;
		else
		{
			fprintf(stderr, "Unknown option: %s\n", argv[i]);
			usage();
		}
	}

	if (*runFile != NULL)
	{
		/* multi-run: 3 positionals — no tiefile on command line */
		*geodatFile   = argv[argc - 3];
		*tiePointFile = NULL;
		*offsetFile   = argv[argc - 2];
		*baselineFile = argv[argc - 1];
	}
	else
	{
		*geodatFile   = argv[argc - 4];
		*tiePointFile = argv[argc - 3];
		*offsetFile   = argv[argc - 2];
		*baselineFile = argv[argc - 1];
	}

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
