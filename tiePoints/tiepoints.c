#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "mosaicSource/common/common.h"
#include "tiePoints.h"
#include <sys/types.h>
#include <sys/time.h>
#include <math.h>
#include <unistd.h>
int32_t sepAscDesc;
/*
  Estimate baseline using tiepoints.

  This program uses some of the routines for geocode,
  which means there is alot of unused junk to initialize everything correctly.
*/

static void readArgs(int argc, char *argv[], int32_t *imageFlag, int32_t *passType,
					 int32_t *noDEM, int32_t *noRamp, char **demFile, char **inputFile,
					 char **tiePointFile, char **phaseFile,
					 char **baselineFile, int32_t *dBpFlag,
					 int32_t *imageCoords, int32_t *motionFlag, int32_t *timeReverseFlag,
					 double *nDays, int32_t *quadB, double *stdLat,
					 int32_t *bnbpFlag, int32_t *bpFlag, int32_t *bnbpdBpFlag,
					 int32_t *bpdBpFlag, int32_t *vrFlag, char **shelfMaskFile,
					 int32_t *yamlOutput, int32_t *verbose, int32_t *debugFlag, char **outputFile,
					 char **ionosphereFile, int32_t *ionosphereMode, double *ionSigmaMargin);

static void usage();

int32_t llConserveMem = 999; /* mem conserve Kluge to maintain backwards compat 9/13/06 */
/*
   Global variables definitions
*/
int32_t RangeSize = RANGESIZE;				/* Range size of complex image */
int32_t AzimuthSize = AZIMUTHSIZE;			/* Azimuth size of complex image */
int32_t BufferSize = BUFFERSIZE;			/* Size of nonoverlap region of the buffer */
int32_t BufferLines = 512;					/* # of lines of nonoverlap in buffer */
char *shelfMaskFile;
int32_t HemiSphere = NORTH;
double Rotation = 45.;
double SLat = -91.0;

int main(int argc, char *argv[])
{
	FILE *tiePointFp;
	demStructure dem;
	tiePointsStructure tiePoints;
	double nDays, stdLat;
	inputImageStructure inputImage; /* input image */
	outputImageStructure outputImage;
	char *demFile, *inputFile, *tiePointFile, *phaseFile, *baselineFile;
	char *outputFile;
	int32_t imageFlag, passType, noDEM, noRamp, dBpFlag, timeReverseFlag;
	int32_t bufferSize;
	int32_t imageCoords;
	int32_t motionFlag, quadB, bnbpFlag, bpFlag, bnbpdBpFlag, bpdBpFlag, vrFlag;
	int32_t yamlOutput, verbose, debugFlag;
	char *ionosphereFile;
	int32_t ionosphereMode;
	double ionSigmaMargin;
	int32_t i, j; /* LCV */
	/*
	   Read command line args and compute filenames
	*/
	GDALAllRegister();
	outputFile = NULL;
	readArgs(argc, argv, &imageFlag, &passType, &noDEM, &noRamp, &demFile,
			 &inputFile, &tiePointFile, &phaseFile, &baselineFile,
			 &dBpFlag, &imageCoords, &motionFlag, &timeReverseFlag, &nDays, &quadB,
			 &stdLat, &bnbpFlag, &bpFlag, &bnbpdBpFlag, &bpdBpFlag, &vrFlag, &shelfMaskFile,
			 &yamlOutput, &verbose, &debugFlag, &outputFile,
			 &ionosphereFile, &ionosphereMode, &ionSigmaMargin);
	/*
	  Parse input file
	*/

	inputImage.imageType = imageFlag;
	inputImage.passType = passType;
	inputImage.stateFlag = TRUE;
	parseInputFile(inputFile, &inputImage);

	tiePoints.motionFlag = motionFlag;
	tiePoints.timeReverseFlag = timeReverseFlag;
	tiePoints.nDays = nDays;
	tiePoints.quadB = quadB;
	tiePoints.bnbpFlag = bnbpFlag;
	tiePoints.bpFlag = bpFlag;
	tiePoints.bnbpdBpFlag = bnbpdBpFlag;
	tiePoints.bpdBpFlag = bpdBpFlag;
	tiePoints.stdLat = stdLat;
	tiePoints.vrFlag = vrFlag;

	/*
	  Input tiepoints
	*/
	tiePointFp = openInputFile(tiePointFile);
	readTiePoints(tiePointFp, &tiePoints, noDEM);

	if (tiePoints.lat[0] < 0)
	{
		HemiSphere = SOUTH;
		tiePoints.stdLat = 71;
		Rotation = 0.0;
	}

	{
		const char *passStr = (passType == ASCENDING) ? "Ascending" : "Descending";
		const char *fitStr  = quadB       ? "Quadratic" :
		                      bnbpFlag    ? "bn,bp"     :
		                      bpFlag      ? "bpOnly"    :
		                      bnbpdBpFlag ? "bn,bp,dBp" :
		                      bpdBpFlag   ? "bp,dBp"    : "Linear";
		const char *hemStr  = (HemiSphere == SOUTH) ? "South" : "North";
		fprintf(stderr, "%s %s %s fit nDays=%.6g%s\n",
		        passStr, hemStr, fitStr, nDays,
		        timeReverseFlag ? " timeReverse" : "");
	}

	tiePoints.noRamp = noRamp;
	tiePoints.dBpFlag = dBpFlag;
	tiePoints.imageCoords = imageCoords;
	if (verbose) {
		if (tiePoints.noRamp == TRUE)
			fprintf(stderr, "NO RAMP\n");
		else
			fprintf(stderr, "RAMP\n");
	}
	/*
	  Call to set up stuff, don't really use the outputimage
	*/
	outputImage.noMem = TRUE;
	initOutputImage(&outputImage, inputImage);
	outputImage.fpLog = stderr;
	/*
	  Read shelf mask if needed
	*/
	if (shelfMaskFile != NULL)
	{
		readShelf(&outputImage, shelfMaskFile);
	}
	else
		outputImage.shelfMask = NULL;
	/*
	  Compute image coords and get z from dem if necessary
	*/
	computeTiePoints(&inputImage, &tiePoints, dem, noDEM, inputFile, outputImage.shelfMask, yamlOutput);
	/*
	  Extract phases from phase file.
	*/
	getPhases(phaseFile, &tiePoints, inputImage, yamlOutput, verbose);
	/*
	  Extract the ionospheric phase at the same tiepoints, if one was given. Sampling it
	  here is equivalent to sampling it after the corrections below, since every downstream
	  correction (addBaselineCorrections, addMotionCorrections) is a purely additive
	  perturbation of phase[] that does not depend on its value.
	*/
	if (ionosphereFile != NULL && ionosphereMode != ION_NONE)
	{
		getIonosphere(ionosphereFile, &tiePoints, inputImage);
	}
	/*
	  Add baseline corrections for previously removed baselines
	*/
	addBaselineCorrections(baselineFile, &tiePoints, inputImage, yamlOutput, verbose);
	/*
	  Motion corrections. Always compute the unsquinted solution; additionally compute the
	  squinted solution when the input geodat carries squint coefficients, so computeBaseline()
	  can fit and emit both (see mosaicSource/CLAUDE.md "Squint"). The squinted pass runs first
	  so that tiePoints.phase/vyra end up holding the unsquinted (default) values for the
	  verbose dump below, matching prior behavior.
	*/
	tiePoints.hasSquintSolution = FALSE;
	tiePoints.phaseSquint = NULL;
	if (motionFlag == TRUE || vrFlag == TRUE)
	{
		if (verbose) fprintf(stderr, "Before motion corrections\n");
		if (inputImage.hasSquintPolynomial)
		{
			tiePoints.phaseSquint = (double *)malloc(tiePoints.npts * sizeof(double));
			addMotionCorrections(inputImage, &tiePoints, TRUE, tiePoints.phaseSquint, verbose);
			tiePoints.hasSquintSolution = TRUE;
		}
		{
			double *phaseNoSquint = (double *)malloc(tiePoints.npts * sizeof(double));
			addMotionCorrections(inputImage, &tiePoints, FALSE, phaseNoSquint, verbose);
			memcpy(tiePoints.phase, phaseNoSquint, tiePoints.npts * sizeof(double));
			free(phaseNoSquint);
		}
		if (verbose) fprintf(stderr, "After motion corrections\n");
	}
	/*
	  Output results for checking to sterr
	*/
	if (verbose)
	{
		for (i = 0; i < tiePoints.npts; i++)
			if (fabs(tiePoints.phase[i]) < 200000)
				fprintf(stderr, "%8.1f %8.1f %8.1f ---  %7.2f %7.2f --- %f ---- %f\n",
						tiePoints.x[i],
						tiePoints.y[i], tiePoints.z[i], tiePoints.r[i], tiePoints.a[i],
						tiePoints.phase[i], tiePoints.vyra[i]);
		fprintf(stderr, "\n");
	}

	/*
	  Estimate baseline solution and output to stdout (or -outputFile, if given)
	*/
	char *debugFile = NULL;
	char debugFileBuf[4096];
	if (debugFlag)
	{
		const char *modeStr = quadB       ? "quadB" :
		                      bnbpFlag    ? "bnbp" :
		                      bpFlag      ? "bpOnly" :
		                      bnbpdBpFlag ? "bnbpdBp" :
		                      bpdBpFlag   ? "bpdBp" : "default";
		if (outputFile != NULL)
			snprintf(debugFileBuf, sizeof(debugFileBuf), "%s.residuals.gpkg", outputFile);
		else
			snprintf(debugFileBuf, sizeof(debugFileBuf), "tiepoints.%s.gpkg", modeStr);
		debugFile = debugFileBuf;
	}

	FILE *outFp = NULL;
	int savedStdoutForFile = -1;
	if (outputFile != NULL)
	{
		outFp = fopen(outputFile, "w");
		if (outFp == NULL)
			error("tiepoints: cannot open -outputFile %s", outputFile);
		savedStdoutForFile = dup(STDOUT_FILENO);
		dup2(fileno(outFp), STDOUT_FILENO);
	}

	computeBaseline(&tiePoints, inputImage, yamlOutput, verbose, debugFile,
					ionosphereMode, ionSigmaMargin, ionosphereFile);

	if (outFp != NULL)
	{
		fflush(stdout);
		dup2(savedStdoutForFile, STDOUT_FILENO);
		close(savedStdoutForFile);
		fclose(outFp);
	}
}

static void readArgs(int argc, char *argv[], int32_t *imageFlag, int32_t *passType,
					 int32_t *noDEM, int32_t *noRamp, char **demFile, char **inputFile,
					 char **tiePointFile, char **phaseFile, char **baselineFile, int32_t *dBpFlag,
					 int32_t *imageCoords, int32_t *motionFlag, int32_t *timeReverseFlag,
					 double *nDays, int32_t *quadB, double *stdLat, int32_t *bnbpFlag, int32_t *bpFlag, int32_t *bnbpdBpFlag,
					 int32_t *bpdBpFlag, int32_t *vrFlag, char **shelfMaskFile, int32_t *yamlOutput, int32_t *verbose,
					 int32_t *debugFlag, char **outputFile,
					 char **ionosphereFile, int32_t *ionosphereMode, double *ionSigmaMargin)
{
	int32_t filenameArg;
	char *argString;
	int32_t i, n;

	if (argc < 5 || argc > 28)
		usage(); /* Check number of args */

	*imageFlag = DESCENDING; /* Default */
	*passType = -1;
	n = argc - 5;
	*noDEM = TRUE;
	*noRamp = FALSE;
	*bpFlag = FALSE;
	*bnbpFlag = FALSE;
	*bpdBpFlag = FALSE;
	*bnbpdBpFlag = FALSE;
	*dBpFlag = FALSE;
	*imageCoords = FALSE;
	*motionFlag = FALSE;
	*timeReverseFlag = FALSE;
	*shelfMaskFile = NULL;
	*nDays = 3.0;
	*stdLat = 70.0;
	*quadB = FALSE;
	*vrFlag = FALSE;
	*yamlOutput = FALSE;
	*verbose = FALSE;
	*debugFlag = FALSE;
	*outputFile = NULL;
	*ionosphereFile = NULL;
	*ionosphereMode = ION_AUTO;
	*ionSigmaMargin = IONSIGMAMARGIN;
	for (i = 1; i <= n; i++)
	{
		argString = strchr(argv[i], '-');
		if (strstr(argString, "descending") != NULL && *passType < 0)
			*passType = DESCENDING;
		else if (strstr(argString, "ascending") != NULL && *passType < 0)
			*passType = ASCENDING;
		else if (strstr(argString, "noDEM") != NULL)
			*noDEM = TRUE;
		else if (strstr(argString, "noRamp") != NULL)
			*noRamp = TRUE;
		else if (strstr(argString, "shelfMask") != NULL)
		{
			*shelfMaskFile = argv[i + 1];
			i++;
		}
		/* Ionosphere flags -- matched longest-first, since the loop uses strstr, not strcmp */
		else if (strstr(argString, "noIonosphere") != NULL)
		{
			*ionosphereMode = ION_NONE;
		}
		else if (strstr(argString, "forceIonosphere") != NULL)
		{
			*ionosphereMode = ION_FORCE;
		}
		else if (strstr(argString, "ionSigmaMargin") != NULL)
		{
			sscanf(argv[i + 1], "%lf", ionSigmaMargin);
			i++;
		}
		else if (strstr(argString, "ionosphere") != NULL)
		{
			*ionosphereFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "quadB") != NULL)
		{
			*quadB = TRUE;
			if (*bpFlag == TRUE || *bnbpFlag == TRUE ||
				*bpdBpFlag == TRUE || *bnbpdBpFlag == TRUE)
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -"
					  "incompatable");
		}
		else if (strstr(argString, "bnbpOnly") != NULL)
		{
			*bnbpFlag = TRUE;
			if (*bpFlag == TRUE || *quadB == TRUE ||
				*bpdBpFlag == TRUE || *bnbpdBpFlag == TRUE)
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -"
					  "incompatable");
		}
		else if (strstr(argString, "bpOnly") != NULL)
		{
			*bpFlag = TRUE;
			if (*quadB == TRUE || *bnbpFlag == TRUE ||
				*bpdBpFlag == TRUE || *bnbpdBpFlag == TRUE)
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -"
					  "incompatable");
		}
		else if (strstr(argString, "bnbpdBpOnly") != NULL)
		{
			*bnbpdBpFlag = TRUE;
			if (*bpFlag == TRUE || *quadB == TRUE ||
				*bpdBpFlag == TRUE || *bnbpFlag == TRUE)
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -"
					  "incompatable");
		}
		else if (strstr(argString, "bpdBpOnly") != NULL)
		{
			*bpdBpFlag = TRUE;
			if (*bpFlag == TRUE || *quadB == TRUE ||
				*bnbpFlag == TRUE || *bnbpdBpFlag == TRUE)
				error("bnbpOnly,bnOnly,bndBn,bnbp,dBnOnly,quadB -"
					  "incompatable");
		}
		else if (strstr(argString, "imageCoords") != NULL)
			*imageCoords = TRUE;
		else if (strstr(argString, "vr") != NULL)
			*vrFlag = TRUE;
		else if (strstr(argString, "dBp") != NULL)
			*dBpFlag = TRUE;
		else if (strstr(argString, "center") != NULL)
			{ if (*verbose) fprintf(stderr, "ignoring obsolete center flag\n"); }
		else if (strstr(argString, "yaml") != NULL)
			*yamlOutput = TRUE;
		else if (strstr(argString, "verbose") != NULL)
			*verbose = TRUE;
		else if (strstr(argString, "motion") != NULL)
			*motionFlag = TRUE;
		else if (strstr(argString, "useSquint") != NULL)
			fprintf(stderr, "note: tiepoints -useSquint is deprecated and ignored -- "
					"tiepoints now always computes both the unsquinted and (when squint "
					"data is available) squinted baseline solutions in -yaml output; "
					"mosaic3d -useSquint selects which one to use.\n");
		else if (strstr(argString, "timeReverse") != NULL)
			*timeReverseFlag = TRUE;
		else if (strstr(argString, "nDays") != NULL)
		{
			sscanf(argv[i + 1], "%lf", nDays);
			i++;
		}
		else if (strstr(argString, "stdLat") != NULL)
		{
			sscanf(argv[i + 1], "%lf", stdLat);
			i++;
		}
		else if (strstr(argString, "outputFile") != NULL)
		{
			*outputFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "debug") != NULL)
			*debugFlag = TRUE;
		else if (i != n)
			usage();
	}
	if (*dBpFlag == TRUE && *noRamp == TRUE)
	{
		fprintf(stderr, "*** -dBp incompatible with -noRamp ***");
		usage();
	}
	else if (*imageCoords == TRUE && *motionFlag == TRUE)
	{
		fprintf(stderr, "*** motionFlag incompatible with imageCoords ***");
		usage();
	}
	if (*passType < 0)
		*passType = DESCENDING;
	*imageFlag = PHASE;
	*demFile = NULL;
	if (*noDEM == FALSE)
		*demFile = argv[argc - 5];
	*inputFile = argv[argc - 4];
	*tiePointFile = argv[argc - 3];
	*phaseFile = argv[argc - 2];
	*baselineFile = argv[argc - 1];
	return;
}

static void usage()
{
	/* Single concatenated literal (no %s indirection) -- this used to be a
	   two-part format string with one %s per argument line, but the argument
	   count had drifted one ahead of the placeholder count (harmless in C,
	   since extra varargs are simply ignored, but it silently dropped the
	   last line, "baselineFile = ..."). Plain concatenation removes that
	   whole class of bug since there's nothing left to miscount. */
	error(
		"\n\n"
		"Compute baseline for given tiepoints.\n"
		"At least 4 tiepoints must be given in the tiepoints file\n"
		"Output is to stdout (or -outputFile, if given)\n"
		"\n"
		"Usage:\n"
		" tiepoints  -vr -imageCoords -motion  -dBp -passType -noRamp -noDEM \\\n"
		"            -shelfMask shelfMask -timeReverse -nDays nDays -quadB  \\\n"
		"            -ionosphere ionFile -noIonosphere -forceIonosphere       \\\n"
		"            -ionSigmaMargin frac                                     \\\n"
		"            dem geoInput tiepointsFile uwPhase baselineFile\n"
		"\n"
		"where\n"
		"  bpOnly      estimate only bp\n"
		"  bnbpOnly    estimate only bn,bp\n"
		"  bpdBpOnly   estimate only bp,dBp\n"
		"  bnbpdBpOnly estimate only bn,bp,dBp\n"
		"  imageCoords  = tiepoints in imageCoords instead of lat/lon\n"
		"  passType     = DESCENDING (default) or ASCENDING\n"
		"  noRamp       = if set, treat apply cos(thetad) depedence to omegaA\n"
		"  dBp          = use deltaBp in place of omegaA\n"
		"  noDEM        = set flag if no DEM used (obsolete - DEM NEVER USED)\n"
		"  shelfMask    = shelfMask file to indicate use tidal corrections\n"
		"  timeReverse  = flag to reverse time when order of orbits switched\n"
		"  nDays        = temporals basline in days (default=3)\n"
		"  quadB        = flag for quadratic fit\n"
		"  dem          = (OMIT if dem not used) dem file for tiepoint elevations\n"
		"  motion       = moving tiepoints, tiepoint file lat,lon,z,vx,vy,vz (m/yr)\n"
		"  vr           = moving tiepoints, tiepoint file lat,lon,z,vr,vz (m/yr)\n"
		"  useSquint    = deprecated, ignored -- both solutions are now always computed\n"
		"                 (in -yaml output) when squintCoefficients are in geoInput\n"
		"  debug        = write every tiepoint used in the fit, plus its residual, to a\n"
		"                 GeoPackage (<outputFile>.residuals.gpkg, or tiepoints.<mode>.gpkg\n"
		"                 if -outputFile not given). Writes both residuals_noSquint and\n"
		"                 residuals_squint layers when squintCoefficients are in geoInput\n"
		"  outputFile   = write the solution to <path> instead of stdout; also names the\n"
		"                 -debug GeoPackage\n"
		"  ionosphere   = ionospheric phase file (radians, same multilooked grid as uwPhase).\n"
		"                 When given, the baseline is fit both with and without the\n"
		"                 correction and the better fit is kept; the choice and the file\n"
		"                 name are recorded in the output so mosaic3d applies the same\n"
		"                 correction to the phase. Omit for the uncorrected fit (default)\n"
		"  noIonosphere      = never apply the correction, even if -ionosphere is given\n"
		"  forceIonosphere   = always apply the correction, without comparing sigmas\n"
		"  ionSigmaMargin f  = fractional sigma improvement the corrected fit must achieve\n"
		"                 before it is used (default 0.05)\n"
		"\n"
		"  geoInputFile = geo params file\n"
		"  tiepointFile = Tiepoint location file in (lat,lon) or (lat,lon,z)\n"
		"  uwPhaseFile  = unwrapped phase image\n"
		"  baselineFile = baseline used to unwrap phases\n");
}
