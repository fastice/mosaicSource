#include <stdio.h>
#include "string.h"
#include "mosaicSource/common/common.h"
#include "mosaic3d.h"
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
#include "mosaicSource/landsatMosaic/landSatMosaic.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include "ogr_srs_api.h"
/*
  Mosaic several insar dems with altimetry dem.


  5/30/07 - Major revision for the 3D part. Modified so that decision to
  derive velocity from two images is base on the track heading rather than asc/desc
  pairs. This required significant changes to the mem conservation scheme. Specifically,
  instead of saving memory seperately for ASC and DESC images, the memory now needs to
  be designated by first or second image in the pair (which will shift). So a memChan
  variable was added to inputImage. If set, this will override the asc/desc allocation.
  So this gets explicitly set before each image is started in the two loops. Likewise
  the image buffers are still pre-allocated, but two global pointers (AImageBuffer and DImageBuffer)
  are used to carry around the locations of the blocks of memory. Prior to reading each image, the
  buffer pointers are set as needed with explicit routine (setBuffer).
*/
static void computeDateRange(char **date1, char **date2, outputImageStructure *outputImage, inputImageStructure *images, vhParams *params, double minJDLS, double maxJDLS);
static void removeOutOfBounds(outputImageStructure *outputImage, inputImageStructure **ascImages, vhParams **ascParams, int32_t *nImages);
static void findOutBounds(outputImageStructure *outputImage, inputImageStructure *ascImages, inputImageStructure *descImages,
						  landSatImage *LSImages, int32_t *autoSize, int32_t writeBlank);
static void mallocOutputImage(outputImageStructure *outputImage);
static void readArgs(int32_t argc, char *argv[], mosaicArgs *args,
					 referenceVelocity *refVel, outputImageStructure *outputImage);
static void write3Doutput(outputImageStructure outputImage, char *outFileBase);
static void write3DTiffOutput(outputImageStructure outputImage, char *outFileBase, char *driverType, const char *epsg, char *date1,  char *date2);
static void removeStaleOutputs(char *outFileBase, int32_t writingTiff);
static int32_t writeMetaFile(inputImageStructure *image, outputImageStructure *outputImage, vhParams *params, char *outFileBase,
							 char *demFile, int32_t writeBlank);
void caldat(int32_t julian, int32_t *mm, int32_t *id, int32_t *iyyy);
static void readReferenceVelMosaic(referenceVelocity *refVel, outputImageStructure *outputImage);
static void processMosaicDate(outputImageStructure *outputImage, char *date1, char *date2);
static void usage();
static void logInputs3d(outputImageStructure *outputImage, char *outFileBase, char *inputFile, char *demFile, char *irregFile,
						char *shelfMaskFile, char *extraTieFile,
						char *tideFile, char *verticalCorrectionFile, float fl, int32_t statsFlag, int32_t threeDOffFlag,
						double tieThresh, referenceVelocity *refVel, mosaicArgs *args);
static void logInputFiles3d(outputImageStructure *outputImage, char **geodatFiles, char **phaseFiles, char **baselineFiles,
							char **rOffsetFiles, char **rParamsFiles, char **offsetFiles, char **azParamsFiles, float *nDays, float *weights,
							int32_t *crossFlags, int32_t nFiles, int32_t offsetFlag);
static void get3DProj(inputImageStructure *ascImages, inputImageStructure *descImages, int32_t nAsc, int32_t nDesc,
					  int32_t northFlag, outputImageStructure *outputImage);
static void malloc3DBuffers();
static void init3DImages(outputImageStructure *outputImage, referenceVelocity *refVel);
static void consolodateLists(inputImageStructure **images, vhParams **params, inputImageStructure *ascImages,
							 inputImageStructure *descImages, vhParams *ascParams, vhParams *descParams, int32_t nAsc, int32_t nDesc);
/* moved to mosaic3d.h
   #define MAXOFFBUF 30000000
   #define MAXOFFLENGTH 30000
*/
/*
   Global variables definitions
*/
int32_t RangeSize = RANGESIZE;				/* Range size of complex image */
int32_t AzimuthSize = AZIMUTHSIZE;			/* Azimuth size of complex image */
int32_t BufferSize = BUFFERSIZE;			/* Size of nonoverlap region of the buffer */
int32_t BufferLines = 512;					/* # of lines of nonoverlap in buffer */
int32_t sepAscDesc = TRUE;
int32_t HemiSphere = NORTH;
double Rotation = 45.;
double SLat = -91.;

int32_t llConserveMem = 1234; /* Kluge to maintain backwards compat 9/13/06 */
int32_t useSquint = FALSE; /* apply squint(r,a) heading correction (phase/make3DMosaic.c only) */
/* Joint (normal-equations) crossing-orbit solvers are the DEFAULT as of 2026-08-29.  They loop
   once per PRODUCT rather than once per pair, so each measurement enters the solution exactly
   once and the n_A x n_D over-count cannot arise.  Measured: same velocity, 12-13% more coverage
   on NISAR, 18.6-71x faster.  -legacyPairPhase / -legacyPairRange restore the originals.
   NOTE the reported ex/ey are then Cov = N^-1, which propagates measurement noise only and is
   optimistic against an external reference (1-sigma coverage ~43% rather than 68%); the missing
   term is common-mode error at the pixel and is not modellable from within.  See
   Documents/crossingOrbitRedundancy.md and mosaicSource/CLAUDE.md. */
int32_t legacyPairPhase = FALSE;
int32_t legacyPairRange = FALSE;
double jointMaxSigma = JOINTMAXSIGMAPHASEDEF;
double jointMaxSigmaRange = JOINTMAXSIGMARANGEDEF;
double jointErrScale = 1.0;
int32_t noAzimuthRows = FALSE;
/* True 3-component solver (mosaicTrue3D.c).  Offsets only, no surface-parallel constraint. */
int32_t true3D = FALSE;
int32_t true3DProject = FALSE;
double true3DMaxSigma = JOINTMAXSIGMARANGEDEF;
int32_t true3DDiag = FALSE;
int32_t true3DPhase = FALSE;
/*  DEFAULT SOLVER since 2026-09.  A bare mosaic3d run uses the 2D hopper: one per-pixel
    normal-equation solve carrying phase, range-offset and azimuth-offset rows, each weighted by
    its own sigma, replacing the legacy four-round pipeline (crossing phase -> crossing range
    offsets -> vh -> speckle).  -legacyCode selects that pipeline instead.  Legacy round-selection
    flags are TRANSLATED into row switches rather than ignored -- see translateLegacyFlags(). */
int32_t hopper = TRUE;
int32_t hopper3D = FALSE;
/*  -legacyCode: the pre-2026-09 four-round pipeline.  Also implied by -stats and by
    -legacyPairPhase/-legacyPairRange, none of which have a hopper equivalent. */
int32_t legacyCode = FALSE;
/*  DEFAULT IS 2D EVERYWHERE (-1), deliberately.  Measured on three NISAR sectors, the
    unconstrained 3-component solve costs horizontal accuracy against the multi-year S1
    reference in EVERY case -- +15.6%, +70.7%, +144.4% -- and the penalty is largest where the
    2D solution is best, because dropping the surface-parallel constraint can only add variance
    to an already-good answer.  So 3D is strictly opt-in: pass a value >= 0 to enable it.
    See Documents/hopper3DPlan.md. */
double hopper3DMaxSigma = -1.0;
int32_t noPhaseRows = FALSE;
int32_t noRangeRows = FALSE;
/*  -sigmaAThreshVel X : drop azimuth rows whose azparams residual, converted to m/yr
    (sigmaAresidual * 365.25/nDays), exceeds X.  The existing -sigmaAThresh is in METRES over the
    pair interval and is left alone.  Default -1 = inert. */
double sigmaAThreshVel = 150.0;
/*  -gateNEff : in the output rejection gate, replace the raw row count with the weighted
    effective count (sum w)^2 / sum(w^2).  The gate estimates per-measurement sigma as
    sigma_worst*sqrt(n), which assumes comparable weights; azimuth rows sit ~4 orders of magnitude
    below phase, so they pad n without informing lambda_min and inflate the gate ~1.7x.
    DEFAULT ON since 2026-09-02; -noGateNEff restores the raw count.  Wired in BOTH hoppers, which
    it must be: -hopper3D -hopper3DMaxSigma -1 is documented to reproduce -hopper, and gating them
    differently would break that. */
int32_t gateNEff = TRUE;
/*  -gateAbsolute : gate on the worst-direction formal sigma ALONE, dropping the sqrt(n) factor,
    so the test reads "do not publish a velocity whose formal 1-sigma exceeds X m/yr".
    The default sigma_worst*sqrt(nEff) form estimates a PER-MEASUREMENT sigma, but sigma_worst
    also carries the geometric dilution: where the look directions are clustered (Antarctic
    coastal rim, offsets-only pixels -- phase and range share the LOS direction) the statistic is
    inflated by geometry rather than by noise, and tightens further as n grows.  Measured on the
    Antarctic multi-year mosaic: the default rejects 5.1M pixels the pair solvers kept, on the
    ENTIRE coastal margin, whose chi-square is 0.23 (i.e. the measurements agree) and whose
    formal error is 9 m/yr.  Greenland is unaffected because its rows are better calibrated
    (chi2 ~0.6 vs ~0.2), which is exactly why a single sqrt(n) threshold does not port between
    archives and an absolute one does.  DEFAULT OFF; -noGateAbsolute restores it.  Wired in BOTH
    hoppers for the same reason as gateNEff. */
int32_t gateAbsolute = FALSE;
/*  -gateSpeedFrac F : speed-aware form of the absolute cap, effective cap = max(X, F*|v|)
    where X is -jointMaxSigma.  See mosaic3d.h.

    DEFAULT 0.03 since 2026-09; pass -gateSpeedFrac 0 to disable.  A fixed cap in m/yr is a
    tightening constraint as speed rises -- 50 m/yr is 0.7% of a 7 km/yr velocity but 50% of a
    100 m/yr one -- so an absolute gate rejects fast ice for having the large ABSOLUTE error that
    fast ice necessarily has.  Measured on Greenland, retention above 5 km/yr: 35.6% (n-normalised
    default), 92.2% (-gateAbsolute 50), 99.8% (+ this term), with NO change below 500 m/yr.

    Note max(): the term can only ever RAISE the cap, so enabling it can add pixels but never
    remove them.  That is what makes a default change safe -- no pixel that passed before fails
    now. */
double gateSpeedFrac = 0.03;
/*  -maxChi2 X : reject a solved pixel whose REDUCED chi-square exceeds X.  A blunder screen,
    not a precision test -- chi2 asks whether the rows agree with EACH OTHER, so it catches the
    case the formal sigma cannot: a confidently wrong solve (one frame with an unwrapping error,
    say) whose reported error is small.  Measured on the Antarctic -gateAbsolute mosaic: the 583
    pixels that changed by >1000 m/yr against the interpolated product carry median chi2 9e5
    while good pixels sit at 0.38, so X=100 removes 81 % of them at 0.00 % cost in Antarctica and
    0.80 % in Greenland (X=1000: 80 % / 0.09 %).  Keep X LOOSE.  A tight cut is wrong: chi2 also
    carries real temporal variability, and 22 % of Greenland's valid pixels exceed 2.
    DEFAULT -1 = off.  Wired in BOTH hoppers. */
double maxChi2 = -1.0;
/* -obsDump <pointsFile>: per-observation dump at a short list of lat/lon points, for the
   GPS forward model (Documents/gpsForwardModelPlan.md).  NULL = inert. */
char *obsDumpFile = NULL;
/* -noErrorGate: keep pixels whose velocity is valid but whose formal error is not.
   Default FALSE, i.e. such pixels ARE removed.  See toSigma(). */
int32_t noErrorGate = FALSE;
/* Inflate crossing-pair errors by (n_A+n_D)/2 to account for the N_A x N_D pairs not being
   independent (see inflatePairOverCount, common/scalingFunctions.c).  Default ON, matching
   the behaviour of the GNSS-validated products; -noPairOverCount disables it.  Without any
   inflation the crossing-orbit errors come out several-fold BELOW the observed GNSS scatter
   (measured: median ex 0.136 vs 2.803 for the same sector of the validated product), so the
   correction is needed; it is only its magnitude that was wrong.  See mosaicSource/CLAUDE.md. */
int32_t pairOverCount = TRUE;
/* Include the rparams tie-point fit residual (offsets.sigmaRresidual, metres of
   slant range) in the range-offset error budget, in quadrature with the local
   matching sigma.  This is the exact analogue of the phase path's
   min(6*PI, vhParam->sigma) term (make3DMosaic.c), which has always been there.
   Without it the offsets budget contains ONLY interpRangeSigma -- the .sr band,
   a local neighbourhood scatter computed after a local plane is removed -- so it
   describes matching noise and is structurally blind to long-wavelength error
   (ionosphere, orbit ramps).  Measured against Sentinel-1 over stable ground,
   that made the offsets formal error ~6x too small AND smaller than the phase
   one, inverting the inverse-variance weighting so crossing offsets dragged the
   combined solution.  Default ON; -noRSigmaResidual restores the old budget. */
extern int32_t rSigmaResidual;
extern int32_t pairCountLegacy;  /* defined in common/getRegion.c */
/* -rSigmaConst X : add a FIXED X metres of slant range to every offsets frame's
   error budget, in quadrature (see rangeAccuracyVar, common/interpOffsets.c).
   Prefer this to -rSigmaResidual: it raises the offsets budget relative to
   phase -- the miscalibration that corrupts the combined product -- without
   disturbing the relative weighting between offset frames, which the per-frame
   residual does at a measured cost to velocity accuracy. */
extern double rSigmaConst;
/* -rhoPhase X / -rhoOffsets X : correlation parameter driving
   inflatePairOverCount()'s factor
       f = rho*nPairs + (1 - rho)*(nOuter + nInner)/2
   (see common/scalingFunctions.c, and common/getRegion.c for the full rationale).
   BOTH default to 0.5, set from the measured threshold dependence: holding the
   images fixed and raising the crossing threshold 3.3x leaves the measured
   accuracy flat (MAD ratio 0.99), so the formal error should be flat too --
   rho = 0 gives 0.52-0.56, rho = 0.5 gives 0.91-0.95.  It is standing in for
   correlation between pairs that share an image, which the approximate
   (n_A+n_D)/2 count misses once a time threshold makes the pairing graph
   non-complete.  Note rho = 0.5 with n_A = n_D reproduces the legacy
   -pairCountLegacy formula exactly (verified on real data to 0.01%). */
extern double rhoPhase;
extern double rhoOffsets;
extern int32_t noMask; /* ignore any embedded VRT dataset mask band on offset inputs; defined once in common/getRegion.c since every program shares readOffsets.c */

int main(int argc, char *argv[])
{
	extern char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
	extern void *lBuf1, *lBuf2, *lBuf3, *lBuf4;
	extern void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
	FILE *fp;
	xyDEM dem, verticalCorrection;					  /* Main dem, flow dir dem */
	referenceVelocity refVel;						  /* Reference velocity used for clipping and initializing */
	landSatImage *LSImages;							  /* List of landsat images to include */
	outputImageStructure outputImage;				  /* Output image */
	vhParams *ascParams, *descParams, *params, *tmpP; /* Param lists */
	irregularData *irregDat, *irregTmp;
	inputImageStructure *ascImages=NULL, *descImages=NULL, *images=NULL, *tmp; /* Lists of asc/desc and all  images */
	double tieThresh;											/* Limiting value used for tiepoints */
	mosaicArgs args;
	char **phaseFiles, **geodatFiles, **baselineFiles;
	float *weights, *nDays;
	double minLSJD = HIGHJD, maxLSJD = 0; // Initial values
	char **azParamsFiles, **offsetFiles; /* Az paramter and offest files */
	char **rOffsetFiles, **rParamsFiles; /* Range offset file */
	int32_t i, j;								  /* LCV */
	int32_t haveData;
	int32_t nAsc, nDesc, nFiles;			  /* Number of ascending/descending and all files */
	int32_t offsetFlag = TRUE;				  /* Flag to indicate do both offset solution where needed, always true old option removed */
	int32_t *crossFlags;
	int32_t autoSize, count;
	const char *epsg=NULL;

	GDALAllRegister();
	if (getenv("OMP_NUM_THREADS") == NULL)
		omp_set_num_threads(4);
	/*
	   Read command line args and compute filenames
	*/
	readArgs(argc, argv, &args, &refVel, &outputImage);
	outputImage.sigmaAThresh = args.sigmaAThresh;
	/* Added August 2021 to set projection parameters from DEM */
	readXYDEMGeoInfo(args.demFile, &dem, TRUE);

	/* Removed no offset flag version */
	processMosaicDate(&outputImage, args.date1, args.date2);
	/* Write to log */
	logInputs3d(&outputImage, args.outFileBase, args.inputFile, args.demFile, args.irregFile, args.shelfMaskFile,
				args.extraTieFile, args.tideFile, args.verticalCorrectionFile, args.fl, args.statsFlag,
				args.threeDOffFlag, args.tieThresh, &refVel, &args);
	/*
	  read inputfile
	*/
	getMVhInputFile(args.inputFile, &phaseFiles, &geodatFiles, &baselineFiles, &offsetFiles, &azParamsFiles, &rOffsetFiles, &rParamsFiles,
					&outputImage, &nDays, &weights, &crossFlags, &nFiles, offsetFlag, outputImage.rOffsetFlag, args.threeDOffFlag);
	if (outputImage.rOffsetFlag == TRUE)
		fprintf(stderr, "*** USING RANGE/AZIMUTH OFFSETS in solution *** \n");
	fprintf(outputImage.fpLog, ";\n; **** SOURCE DATA **** \n");
	/*
	  Parse inputfiles and set everything up.
	*/
	setup3D(nFiles, phaseFiles, geodatFiles, baselineFiles, offsetFiles, azParamsFiles, rOffsetFiles, rParamsFiles, nDays, weights, crossFlags,
			&ascImages, &descImages, &ascParams, &descParams, &nAsc, &nDesc, offsetFlag, outputImage.rOffsetFlag, args.threeDOffFlag,
			outputImage.fpLog, &outputImage);
	logInputFiles3d(&outputImage, geodatFiles, phaseFiles, baselineFiles, rOffsetFiles, rParamsFiles, offsetFiles, azParamsFiles,
					nDays, weights, crossFlags, nFiles, offsetFlag);

	fprintf(stderr, "%s %s\n", args.date1, args.date2);

	/*
	  Determine hemisphere
	*/
	get3DProj(ascImages, descImages, nAsc, nDesc, args.north, &outputImage);
	outputImage.slat = dem.stdLat;
	/* Process landsat images */
	LSImages = NULL;
	if (args.landSatFile != NULL)
	{
		LSImages = parseLSInputs(args.landSatFile, LSImages, outputImage.jd1, outputImage.jd2, outputImage.timeOverlapFlag, &minLSJD, &maxLSJD);
		fprintf(stderr, "min/max JD from Landsat %f  %f\n", minLSJD, maxLSJD);
	}
	/* Find bounding box */
	fprintf(stderr, "nAsc/nDesc %i %i\n", nAsc, nDesc);
	findOutBounds(&outputImage, ascImages, descImages, LSImages, &autoSize, args.writeBlank);
	fprintf(stderr, "xSize=%d, ySize=%d\n",outputImage.xSize, outputImage.ySize );
	/* Remove images that are outside output area */
	removeOutOfBounds(&outputImage, &ascImages, &ascParams, &nAsc);
	removeOutOfBounds(&outputImage, &descImages, &descParams, &nDesc);
	/* Allocate image buffers sized to surviving images only */
	allocateOffsetBuffers(ascImages, descImages);
	/*
	  Read shelf mask
	*/
	if (args.shelfMaskFile != NULL)
	{
		fprintf(outputImage.fpLog, "; Starting reading Shelf mask file\n");
		fflush(outputImage.fpLog);
		readShelf(&outputImage, args.shelfMaskFile);
		fprintf(outputImage.fpLog, "; Finished reading Shelf mask file\n");
		fflush(outputImage.fpLog);
	}
	else
	{
		outputImage.shelfMask = NULL;
	}

	if (args.verticalCorrectionFile != NULL)
	{
		fprintf(outputImage.fpLog, "; Starting reading vertical correction file\n");
		fflush(outputImage.fpLog);
		readXYDEM(args.verticalCorrectionFile, &verticalCorrection);
		outputImage.verticalCorrection = &verticalCorrection;
		fprintf(outputImage.fpLog, "; Finished reading vertical correction file\n");
		fflush(outputImage.fpLog);
	}
	else
	{
		outputImage.verticalCorrection = NULL;
	}
	/* Unconditional: only make3DMosaicJoint() ever sets this, and a stack-allocated
	   outputImage would otherwise leave it as garbage and trigger a bogus .nobs write. */
	outputImage.jointNObs = NULL;
	outputImage.jointChi2 = NULL;
	/* Only mosaicTrue3D() ever sets these; NULL is what tells the writer below to skip the
	   ".vz3d"/".ez" bands. */
	outputImage.vZ3D = NULL;
	outputImage.errorZ = NULL;
	outputImage.hopperMode = NULL;
	/*
	   read ref vel file for error clipping
	*/
	if (refVel.velFile != NULL && (nAsc + nDesc) > 0)
	{
		readReferenceVelMosaic(&refVel, &outputImage);
	}
	/*
	  Init output image memory
	*/
	/* Single contributing image, nothing else touching the output grid: the weighted
	   multi-image accumulation buffers/math are provably a no-op (see
	   mosaicSource/CLAUDE.md "single-image fast path"), so mallocOutputImage() and
	   speckleTrackMosaic() can skip them. Deliberately narrow/auto-detected (no CLI
	   flag) rather than covering every mode. outputImage.no3d must be TRUE too:
	   make3DMosaic() is called unconditionally below whenever (nAsc+nDesc)>0 (it is
	   NOT gated by noVhFlag) and only skips touching the accumulation buffers via its
	   own early "if (no3d == TRUE) return;" -- without that, it would call
	   setupBuffers()/computeScale() on the buffers this path leaves NULL. */
	/*  legacyCode is REQUIRED: this path leaves the accumulation buffers NULL and relies on the
	    solver's own "if (no3d == TRUE) return;" to avoid touching them.  Under the hopper -no3d
	    no longer short-circuits the solver, so without this gate the hopper would dereference
	    those NULL buffers. */
	outputImage.singleImageFastPath =
		(legacyCode == TRUE &&
		 outputImage.rOffsetFlag == TRUE && (nAsc + nDesc) == 1 &&
		 args.threeDOffFlag == FALSE && outputImage.noVhFlag == TRUE &&
		 outputImage.no3d == TRUE && true3D == FALSE &&
		 args.landSatFile == NULL && args.irregFile == NULL &&
		 outputImage.makeTies == FALSE && args.statsFlag == FALSE &&
		 refVel.initMapFlag == FALSE && outputImage.timeOverlapFlag == FALSE);
	malloc3DBuffers();
	mallocOutputImage(&outputImage);
	/* Consolodate lists */
	consolodateLists(&images, &params, ascImages, descImages, ascParams, descParams, nAsc, nDesc);
	computeDateRange(&args.date1, &args.date2, &outputImage, images, params, minLSJD, maxLSJD);
	tmpP = params;
	for (tmp = images; tmp != NULL; tmp = tmp->next, tmpP = tmpP->next)
	{
		fprintf(stderr, "%s %s\n", tmpP->offsets.file, tmp->file);
		if (strstr(tmp->file, "nophase") != NULL)
		{
			fprintf(stderr, "Skipping %s :no phase given\n", tmp->file);
			continue;
		}
		fprintf(stderr, "--\n");
	}
	//  Init and input DEM
	fprintf(outputImage.fpLog, ";\n; About to Read XYDEM %s\n", args.demFile);
	/* readXYDEM() ("Read a full XY DEM") loads the entire DEM file regardless of how
	   much of it the output grid actually needs -- for a circumpolar DEM (e.g. the
	   full ~20741x20741 Antarctic Copernicus DEM) that can be several times larger
	   than a single-track output grid actually requires. Crop to the output grid's
	   own extent instead, via the readXYDEMcrop() this file already delegates to
	   (readXYDEMcrop.c:readXYImageGDALCropped()/getCropBounds() already support
	   real cropping; readXYDEM() just always passes the 0,0,0,0 "no crop" sentinel).
	   Units are km, matching getCropBounds()'s internal convention (outputImage's
	   own origin/size/spacing are metres). Padded by demPadKm for the DEM
	   interpolation stencil (xyGetZandSlope.c/interpXYDEM.c) and any small mismatch
	   between the DEM's own projection parameters and the output grid's -- cheap
	   insurance given the crop is already ~10x smaller than the full DEM. */
	if ((nAsc + nDesc) == 0 && LSImages == NULL)
	{
		/* Tiepoint-only run with no input frames -- e.g. a bedrock/-extraTies
		   tie that just echoes the extra tie points (writeExtraTies() in
		   writeTieFile.c). There is no meaningful output grid: findOutBounds()
		   hits its "no data case", then autosizes from zero images, leaving a
		   garbage/overflowed extent (xSize can be INT_MIN). Cropping the DEM to
		   that grid throws "Area requested does not fit in DEM" and aborts the
		   tie generation (this is what broke bedrock tiefiles when the crop was
		   added -- readXYDEM used to load the whole DEM). writeExtraTies still
		   needs the DEM (getXYHeight for z, dem->stdLat for the projection) and
		   the extra ties span the whole region, so read the FULL DEM here. */
		readXYDEM(args.demFile, &dem);
	}
	else
	{
		double demPadKm = 10.0;
		double demXmin = outputImage.originX * MTOKM - demPadKm;
		double demXmax = (outputImage.originX + outputImage.xSize * outputImage.deltaX) * MTOKM + demPadKm;
		double demYmin = outputImage.originY * MTOKM - demPadKm;
		double demYmax = (outputImage.originY + outputImage.ySize * outputImage.deltaY) * MTOKM + demPadKm;
		readXYDEMcrop(args.demFile, &dem, demXmin, demXmax, demYmin, demYmax);
	}
	for (tmpP = params; tmpP != NULL; tmpP = tmpP->next)
		tmpP->xydem = dem;
	fprintf(outputImage.fpLog, ";\n; Returned from readXYDEM\n");
	fflush(outputImage.fpLog);
	/*  Name the solver and the row set actually in use.  The default changed 2026-09 from the
	    four-round pipeline to the hopper, so an unannotated log would leave no trace of which
	    one produced a given product. */
	{
		const char *solverName = (legacyCode == TRUE) ? "legacy 4-round pipeline"
							   : (hopper3D == TRUE)
									 ? ((hopper3DMaxSigma >= 0.0) ? "hopper3D (3D enabled)"
																  : "hopper3D (forced 2D)")
									 : "hopper2D";
		int32_t k;
		for (k = 0; k < 2; k++)
		{
			FILE *fp = (k == 0) ? stderr : outputImage.fpLog;
			fprintf(fp, "; SOLVER: %s", solverName);
			if (legacyCode == FALSE)
			{
				fprintf(fp, "  rows:%s%s%s",
						noPhaseRows ? "" : " phase", noRangeRows ? "" : " range",
						noAzimuthRows ? "" : " azimuth");
				if (noPhaseRows && noRangeRows && noAzimuthRows)
					fprintf(fp, " NONE -- this will produce an empty mosaic");
			}
			fprintf(fp, "\n");
		}
		fflush(outputImage.fpLog);
	}
	/*  Per-observation dump.  Must come AFTER findOutBounds/readXYDEM: it maps lat/lon to
	    output-grid cells, so it needs the final grid AND dem.stdLat (outputImage.slat is
	    flagged "not fully implemented" in geocode.h).  Output is per-sector-named, since
	    makemosaic runs one mosaic3d per sector and a fixed name would collide. */
	if (obsDumpFile != NULL)
	{
		char obsDumpOut[2048];
		snprintf(obsDumpOut, sizeof(obsDumpOut), "%s.obsDump", args.outFileBase);
		obsDumpInit(obsDumpFile, obsDumpOut, &outputImage, dem.stdLat);
	}
	// Init values
	/* init3DImages() unconditionally touches image/image2/image3/scale/scale2/scale3,
	   which mallocOutputImage() left NULL for singleImageFastPath -- mallocOutputImage()
	   already did the equivalent sentinel init for that path's own buffers. */
	if (outputImage.singleImageFastPath == FALSE)
		init3DImages(&outputImage, &refVel);
	/*
	  Step 0: Switched to first map since to accomdate discard of large dt.
	  *******************************START Landsat mosaics******************************
	  */
	if (args.landSatFile != NULL && LSImages != NULL)
	{
		if (args.statsFlag == TRUE)
			error("Landsat  incompatible with stats flag, which is for speckle tracked offsets only\n");
		makeLandSatMosaic(LSImages, &outputImage, args.fl);
	}
	/*
	  Step 0: Mosaic using ascending and descending data where possible.
	*/
	if ((nAsc + nDesc) > 0)
	{
		/* Both routines take the same args and both early-return on no3d BEFORE touching the
		   accumulation buffers, so the singleImageFastPath precondition above (which requires
		   no3d == TRUE) holds identically for either. */
		/*  no3d is passed as FALSE to the hoppers: under the hopper it is a legacy round
		    selector that has already been translated into noPhaseRows, so letting it reach the
		    solver's early return would emit an EMPTY mosaic (the pre-2026-09 behaviour, and a
		    silent one -- a blank product looks identical to a failed run). */
		if (hopper3D == TRUE)
		{
			mosaicHopper3D(images, descImages, params, descParams, &dem, &outputImage, args.fl,
						   FALSE, args.timeThreshPhase, &refVel);
		}
		else if (hopper == TRUE)
		{
			mosaicHopper(images, descImages, params, descParams, &dem, &outputImage, args.fl,
						 FALSE, args.timeThreshPhase, &refVel);
		}
		else if (true3DPhase == TRUE)
		{
			mosaicTrue3DPhase(images, descImages, params, descParams, &dem, &outputImage, args.fl,
							  outputImage.no3d, args.timeThreshPhase);
		}
		else if (legacyPairPhase == TRUE)
		{
			make3DMosaic(images, descImages, params, descParams, &dem, &outputImage, args.fl, outputImage.no3d, args.timeThreshPhase);
		}
		else
		{
			make3DMosaicJoint(images, descImages, params, descParams, &dem, &outputImage, args.fl, outputImage.no3d, args.timeThreshPhase);
		}
	}
	/*
	  Step 1:
	*/
	/*  A hopper round already consumed the range offsets as rows, so the crossing-offsets round
	    must not also run -- that would enter the same measurements twice.  Step 3 (speckle) is
	    gated the same way further down; step 0 dispatches to the hopper instead of the crossing
	    phase solver, so phase is covered by construction. */
	if ((nAsc + nDesc) > 0 && args.threeDOffFlag == TRUE && legacyCode == TRUE)
	{
		if (legacyPairRange == TRUE)
		{
			make3DOffsets(images, params, &dem, &outputImage, args.fl, args.timeThresh);
		}
		else
		{
			make3DOffsetsJoint(images, params, &dem, &outputImage, args.fl, args.timeThresh);
		}
		fprintf(stderr, "End of 3d offsets %i\n", args.threeDOffFlag);
	}
	/*
	   Step 2: Make phase/az offset velocity
	*/
	if (outputImage.noVhFlag == FALSE && (nAsc + nDesc) > 0 && legacyCode == TRUE)
	{
		makeVhMosaic(images, params, &outputImage, args.fl);
	}
	/* Same misattribution as the speckleTrackMosaic branch below -- split for the
	   same reason. Existing message text left verbatim, typo included, so anything
	   grepping for it keeps working. */
	else if (outputImage.noVhFlag == TRUE)
	{
		fprintf(outputImage.fpLog, ";\n; NoVh flag set, not Enterng makeVhMosaic\n;\n");
	}
	else
	{
		fprintf(outputImage.fpLog, ";\n; No ascending/descending images for this piece, not Enterng makeVhMosaic\n;\n");
	}
	/*
	  Step 3: Include fully speckle-tracked data
	*/
	/*  -hopper already consumed the range and azimuth offsets as rows in its own solve, so the
	    speckle round must NOT also run -- that would enter the same measurements twice.
	    -rOffsets is still required with -hopper, because it is what makes setup3D parse the
	    offsets files in the first place. */
	if (outputImage.rOffsetFlag == TRUE && (nAsc + nDesc) > 0 && legacyCode == TRUE)
	{
		if (true3D == TRUE)
		{
			mosaicTrue3D(images, params, &outputImage, args.fl, &refVel);
		}
		else
		{
			speckleTrackMosaic(images, params, &outputImage, args.fl, &refVel, args.statsFlag);
		}
	}
	/* Two distinct reasons to skip, and the message has to say which. Reporting the
	   flag unconditionally is wrong whenever the flag is set but the sector simply
	   has no images -- a reader then sees "rOffset flag False" in a log whose own
	   header block says "rOffset Flag : 1", which is a contradiction that invites
	   the conclusion that speckle tracking never ran anywhere in the mosaic. */
	else if (outputImage.rOffsetFlag == FALSE)
	{
		fprintf(outputImage.fpLog, ";\n; rOffset flag False, not Entering speckleTrackMosaic\n;\n");
	}
	else
	{
		fprintf(outputImage.fpLog, ";\n; No ascending/descending images for this piece, not Entering speckleTrackMosaic\n;\n");
	}
	/********************************END Landsat mosaics******************************	*/
	/*
	  Step 6: Include irregularly interpolated data, if specified.
	*/
	if (args.irregFile != NULL)
	{
		fprintf(stderr, "*** incorporating irregularly gridded data from file %s :  ***\n\n", args.irregFile);
		irregDat = NULL;
		parseIrregFile(args.irregFile, &irregDat);
		for (irregTmp = irregDat; irregTmp != NULL; irregTmp = irregTmp->next)
		{
			fprintf(stderr, "|%s|\n", irregTmp->file);
			irregTmp->maxLength = 15;
			irregTmp->maxArea = 75.;
		}
		getIrregData(irregDat);
		addIrregData(irregDat, &outputImage, args.fl);
	}
	/*
	   write meta file
	*/
	haveData = writeMetaFile(images, &outputImage, params, args.outFileBase, args.demFile, args.writeBlank);
	if (args.landSatFile != NULL)
		haveData = TRUE;
	/*
	  Output result
	*/
	if (outputImage.makeTies == FALSE)
	{
		if (haveData == TRUE || refVel.initMapFlag == TRUE)
		{
			/* Drop any stale output of the OTHER format so a dir that
			   switches between binary and GeoTIFF does not keep leftover
			   files of the format no longer being written. */
			removeStaleOutputs(args.outFileBase, (args.COG == TRUE || args.GTiff == TRUE));
			if(args.COG == FALSE && args.GTiff == FALSE) {
				write3Doutput(outputImage, args.outFileBase);
			}
			else
			{
				char *driverType;
				if(args.COG == TRUE) driverType = "COG"; else driverType = "GTiff";
				write3DTiffOutput(outputImage, args.outFileBase, driverType, epsg, args.date1, args.date2);
			}

			remove("NoOutput_noDataInRange");
		}
		else
		{
			fp = fopen("NoOutput_noDataInRange", "w");
			error("No output because no data in range, check dates \n");
		}
	}
	else
	{
		writeTieFile(&outputImage, &dem, &verticalCorrection, args.outFileBase, args.tieThresh, args.extraTieFile, args.tideFile, autoSize);
	}
	obsDumpClose(); /* no-op unless -obsDump was given */
}

static void computeDateRange(char **date1, char **date2, outputImageStructure *outputImage, inputImageStructure *images,
							 vhParams *params, double minJDLS, double maxJDLS)
{
	double jd1 = minJDLS, jd2 = maxJDLS;
	int year, month, day, hour, minute, second;
	/* Already supplied on the command line - nothing to compute. */
	if (*date1 != NULL && *date2 != NULL)
	{
		return;
	}
	/* A Landsat-only mosaic has no SAR images, but parseLSInputs has still
	   filled in the Landsat julian range, which is what the dates come from.
	   Returning here left date1/date2 NULL, and write3DTiffOutput then passed
	   NULL to insert_node -> strdup(NULL) -> segfault under -GTiff. */
	if (images == NULL && !(maxJDLS > 0. && minJDLS < HIGHJD))
	{
		return;
	}
	if(images != NULL) {
		vhParams *currentParams;
		inputImageStructure *currentImage;
		for (currentImage = images, currentParams = params;
			 currentImage != NULL;
			 currentImage = currentImage->next, currentParams = currentParams->next)
		{
			// fprintf(stderr, "jd %lf", currentImage->julDay);
			jd1 = min(jd1, currentImage->julDay);
			jd2 = max(jd2, currentImage->julDay + currentParams->nDays);
		}
	}
	/* else: Landsat only, so jd1/jd2 keep the Landsat range they started with */
	jd_to_date_and_time(jd1, &year, &month, &day, &hour, &minute, &second);
	if(*date1 == NULL) { *date1 = malloc(11); sprintf(*date1, "%4d-%02d-%02d", year, month, day);}
	jd_to_date_and_time(jd2, &year, &month, &day, &hour, &minute, &second);
	if(*date2 == NULL) { *date2 = malloc(11); sprintf(*date2, "%4d-%02d-%02d", year, month, day);}
	fprintf(stderr, "%s %s", *date1, *date2);
}

static void removeOutOfBounds(outputImageStructure *outputImage, inputImageStructure **ascImages, vhParams **ascParams, int32_t *nImages)
{
	/* Go through lists of images and remove ones that are outside the output area */
	int32_t iMin, iMax, jMin, jMax;
	int32_t count, countAdd = 0, countRemove = 0;
	/* Skip empty list */
	if (*ascImages == NULL || *nImages == 0)
		return;
	inputImageStructure *newList, *newListHead, *tmp;
	vhParams *ptmp, *newParams, *newParamsHead;
	/* Init list */
	newList = NULL;
	newParams = NULL;
	newListHead = NULL;
	newParamsHead = NULL;
	count = 0;
	for (tmp = *ascImages, ptmp = *ascParams; tmp != NULL; tmp = tmp->next, ptmp = ptmp->next)
	{
		/* Keep only images with overlap and accessible files */
		if (getRegion(tmp, &iMin, &iMax, &jMin, &jMax, outputImage))
		{
			if (newList == NULL)
			{
				newList = tmp;
				newListHead = tmp;
				newParams = ptmp;
				newParamsHead = ptmp;
			}
			else
			{
				newList->next = tmp;
				newList = tmp;
				newParams->next = ptmp;
				newParams = ptmp;
			}
			countAdd++;
		}
		else
		{
			*nImages -= 1;
			countRemove++;
		}
	}
	/* Terminate lists */
	if (newList != NULL)
	{
		newList->next = NULL;
		newParams->next = NULL;
	}
	*ascImages = newListHead;
	*ascParams = newParamsHead;
}

static void consolodateLists(inputImageStructure **images, vhParams **params, inputImageStructure *ascImages, inputImageStructure *descImages,
							 vhParams *ascParams, vhParams *descParams, int32_t nAsc, int32_t nDesc)
{
	inputImageStructure *tmp; /* Tmp list for looping */
	vhParams *tmpP;
	fprintf(stderr, "**** nAsc %i, nDesc %i\n", nAsc, nDesc);
	if (nAsc == 0 && nDesc > 0)
	{
		*images = descImages;
		*params = descParams;
	}
	else if (nAsc > 0 && nDesc == 0)
	{
		*images = ascImages;
		*params = ascParams;
	}
	else if (nAsc > 0 && nDesc > 0)
	{
		*images = ascImages;
		/* Skip to end of list */
		for (tmp = ascImages; tmp->next != NULL; tmp = tmp->next)
			;
		/* append */
		tmp->next = descImages;
		*params = ascParams;
		for (tmpP = ascParams; tmpP->next != NULL; tmpP = tmpP->next)
			;
		tmpP->next = descParams;
	}
	else
	{
		*images = NULL;
		*params = NULL;
	}
}

static void init3DImages(outputImageStructure *outputImage, referenceVelocity *refVel)
{
	float **vXimage, **vYimage, **vZimage;
	float **scaleX, **scaleY, **scaleZ;
	float **errorX, **errorY;
	float vxPt, vyPt, exPt, eyPt;
	double x, y;
	int32_t i, j;

	vXimage = (float **)outputImage->image;
	vYimage = (float **)outputImage->image2;
	vZimage = (float **)outputImage->image3;
	errorX = (float **)outputImage->errorX;
	errorY = (float **)outputImage->errorY;
	scaleX = (float **)outputImage->scale;
	scaleY = (float **)outputImage->scale2;
	scaleZ = (float **)outputImage->scale3;
	for (i = 0; i < outputImage->ySize; i++)
	{
		y = (outputImage->originY + i * outputImage->deltaY) * MTOKM;
		for (j = 0; j < outputImage->xSize; j++)
		{
			x = (outputImage->originX + j * outputImage->deltaX) * MTOKM;
			if (refVel->initMapFlag == TRUE)
			{
				/*	fprintf(stderr,"%i %i ",i,j); */
				if (outputImage->timeOverlapFlag == TRUE)
					error("can't use reference map and timeoverlap flag at the same time");
				refVelInterp(x, y, refVel, &vxPt, &vyPt, &exPt, &eyPt);
				if (vxPt > -2e6 && vxPt < 2e6)
				{
					scaleX[i][j] = 1.0 / (exPt * exPt);
					vXimage[i][j] = vxPt;
					scaleY[i][j] = 1.0 / (eyPt * eyPt);
					vYimage[i][j] = vyPt;
					scaleZ[i][j] = 0.0;
					vZimage[i][j] = 0.0;
					errorX[i][j] = exPt * exPt;
					errorY[i][j] = eyPt * eyPt;
				}
			}
			else
			{
				errorX[i][j] = 0.0;
				errorY[i][j] = 0.0;
				scaleX[i][j] = 0.0;
				vXimage[i][j] = -2.e9;
				scaleY[i][j] = 0.0;
				vYimage[i][j] = -2.e9;
				scaleZ[i][j] = 0.0;
				vZimage[i][j] = -2.e9;
			}
		}
	}
}

static void malloc3DBuffers()
{
	extern char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2, *SEBuf;
	extern void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4, *offSEBuffSpace;
	extern void *lBuf1, *lBuf2, *lBuf3, *lBuf4, *lSEBuf;

	Abuf1 = malloc(MAXADBUF);
	Abuf2 = malloc(MAXADBUF2);
	Dbuf1 = malloc(MAXADBUF);
	Dbuf2 = malloc(MAXADBUF2);
	SEBuf = malloc(max(MAXADBUF, MAXADBUF2));
	/*
	  Init buffer memory
	*/
	offBufSpace1 = (void *)malloc(MAXOFFBUF);
	offBufSpace2 = (void *)malloc(MAXOFFBUF);
	offBufSpace3 = (void *)malloc(MAXOFFBUF);
	offBufSpace4 = (void *)malloc(MAXOFFBUF);
	offSEBuffSpace = (void *)malloc(MAXOFFBUF);

	/* Used for pointers rows to above buffer space */
	lBuf1 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	lBuf2 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	lBuf3 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	lBuf4 = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	lSEBuf = (void *)malloc(sizeof(float *) * MAXOFFLENGTH);
	if(lSEBuf == NULL || lBuf1 == NULL || lBuf2 == NULL || lBuf3 == NULL || lBuf4 == NULL || offBufSpace1 == NULL || offBufSpace2 == NULL || offBufSpace3 == NULL || offBufSpace4 == NULL ||
	   Abuf1 == NULL || Abuf2 == NULL || Dbuf1 == NULL || Dbuf2 == NULL || SEBuf == NULL)
		error("Error allocating 3D buffers");
}

static void get3DProj(inputImageStructure *ascImages, inputImageStructure *descImages, int32_t nAsc, int32_t nDesc, int32_t northFlag, outputImageStructure *outputImage)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	float lat;

	fprintf(outputImage->fpLog, "; Determining Hemisphere\n;\n");
	if (nAsc != 0)
		lat = ascImages->latControlPoints[1];
	else if (nDesc != 0)
		lat = descImages->latControlPoints[1];
	else
	{
		if (northFlag == FALSE)
		{
			fprintf(stderr, "*** No Asc or Desc images specified - assuming south ***");
			lat = -80.;
		}
		else
		{
			fprintf(outputImage->fpLog, "; *** No Asc or Desc images specified - using north flag ***");
			lat = 80.;
		}
	}
	/* modifed august 2021 to use dem to set projection */
	if (lat < 0)
	{
		HemiSphere = SOUTH; /* Rotation=0.0;  outputImage->slat=71.0;*/
		fprintf(stderr, "\n ***** Southern Hemisphere  ******\n");
		fprintf(outputImage->fpLog, "; ***** Southern Hemisphere  ******\n;\n");
	}
	else
	{
		/* outputImage->slat=70.0;*/
		HemiSphere = NORTH; /* Rotation=0.0;  outputImage->slat=71.0;*/
		fprintf(stderr, "\n ***** Northern Hemisphere  ******\n");
		fprintf(outputImage->fpLog, "; ***** Northern Hemisphere  ******\n;\n");
	}
	fflush(outputImage->fpLog);
}

static void logInputFiles3d(outputImageStructure *outputImage, char **geodatFiles, char **phaseFiles, char **baselineFiles, char **rOffsetFiles,
							char **rParamsFiles, char **offsetFiles, char **azParamsFiles, float *nDays, float *weights, int32_t *crossFlags, int32_t nFiles, int32_t offsetFlag)
{
	int32_t count = 0, i;
	/* loop through and log each file */
	for (i = 0; i < nFiles; i++)
	{
		if (geodatFiles[i] != NULL)
		{
			count++;
			fprintf(outputImage->fpLog, ";\n; Data Take %i\n", count);
			fprintf(outputImage->fpLog, "; phaseFile     : %s\n; Geodat File   : %s\n; Baseline File : %s\n", phaseFiles[i], geodatFiles[i], baselineFiles[i]);
			fprintf(outputImage->fpLog, "; nDays         : %f", nDays[i]);
			fprintf(outputImage->fpLog, "; Weight        : %f\n", weights[i]);
			fprintf(outputImage->fpLog, "; CrossFlag        : %i\n", crossFlags[i]);
			if (offsetFlag == TRUE && offsetFiles[i] != NULL && azParamsFiles[i] != NULL)
				fprintf(outputImage->fpLog, "; Az offsetFile : %s\n; Az Params     : %s\n", offsetFiles[i], azParamsFiles[i]);
			if (outputImage->rOffsetFlag == TRUE && rOffsetFiles[i] != NULL)
				fprintf(outputImage->fpLog, "; Range Offsets : %s\n; Rbaseline     : %s\n", rOffsetFiles[i], rParamsFiles[i]);
			else
				fprintf(outputImage->fpLog, "; Range Offsets : %s\n; Rbaseline     : %s\n", "none", "none");
		}
	}
	fflush(outputImage->fpLog);
}

static void logInputs3d(outputImageStructure *outputImage, char *outFileBase, char *inputFile, char *demFile, char *irregFile, char *shelfMaskFile,
						char *extraTieFile, char *tideFile, char *verticalCorrectionFile, float fl, int32_t statsFlag,
						int32_t threeDOffFlag, double tieThresh, referenceVelocity *refVel, mosaicArgs *args)
{
	char *logFile;
	logFile = (char *)malloc(strlen(outFileBase) + 16);
	logFile[0] = '\0';
	logFile = strcat(logFile, outFileBase);
	logFile = strcat(logFile, ".log");
	outputImage->fpLog = fopen(logFile, "w");
	if (outputImage->fpLog == NULL)
		error("Could not open log file %s", logFile);

	fprintf(outputImage->fpLog, "; Mosaic3d log file\n");
	fprintf(outputImage->fpLog, ";\n; ***Command Line Parameters***\n;\n");
	fprintf(outputImage->fpLog, "; inputFile        : %s\n", inputFile);
	fprintf(outputImage->fpLog, "; demFile          : %s\n", demFile);
	if (irregFile != NULL)
		fprintf(outputImage->fpLog, "; irregFile        : %s\n", irregFile);
	else
		fprintf(outputImage->fpLog, "; irregFile        : %s\n", "None");
	if (shelfMaskFile != NULL)
		fprintf(outputImage->fpLog, "; shelfMaskFile    : %s\n", shelfMaskFile);
	else
		fprintf(outputImage->fpLog, "; shelfMaskFile        : %s\n", "None");
	if (outputImage->makeTies == TRUE && extraTieFile != NULL)
		fprintf(outputImage->fpLog, "; extraTieFile     : %s\n", extraTieFile);
	else
		fprintf(outputImage->fpLog, "; extraTieFile     : %s\n", "None");
	if (tideFile != NULL)
		fprintf(outputImage->fpLog, "; tideFile         : %s\n", tideFile);
	else
		fprintf(outputImage->fpLog, "; tideFile         : %s\n", "None");
	if (verticalCorrectionFile != NULL)
		fprintf(outputImage->fpLog, "; verticalCorrection         : %s\n", verticalCorrectionFile);
	else
		fprintf(outputImage->fpLog, "; verticalCorrection         : %s\n", "None");
	fprintf(outputImage->fpLog, "; vertCorrSuffix   : %s\n", outputImage->verticalCorrectionSuffix != NULL ? outputImage->verticalCorrectionSuffix : "None");
	fprintf(outputImage->fpLog, "; landSatFile      : %s\n", args->landSatFile != NULL ? args->landSatFile : "None");
	fprintf(outputImage->fpLog, "; date1            : %s\n", args->date1 != NULL ? args->date1 : "None");
	fprintf(outputImage->fpLog, "; date2            : %s\n", args->date2 != NULL ? args->date2 : "None");
	fprintf(outputImage->fpLog, "; OutputFile base  : %s\n", outFileBase);
	fprintf(outputImage->fpLog, "; Feather length   : %f\n", fl);
	fprintf(outputImage->fpLog, "; sigmaAThresh     : %f\n", outputImage->sigmaAThresh);
	fprintf(outputImage->fpLog, "; timeThresh       : %f\n", args->timeThresh);
	fprintf(outputImage->fpLog, "; timeThreshPhase  : %f\n", args->timeThreshPhase);
	fprintf(outputImage->fpLog, "; NoVh Flag        : %i\n", outputImage->noVhFlag);
	fprintf(outputImage->fpLog, "; rOffset Flag     : %i\n", outputImage->rOffsetFlag);
	fprintf(outputImage->fpLog, "; No3D Flag        : %i\n", outputImage->no3d);
	fprintf(outputImage->fpLog, "; Stats Flag        : %i\n", statsFlag);
	fprintf(outputImage->fpLog, "; Maketies Flag    : %i\n", outputImage->makeTies);
	fprintf(outputImage->fpLog, "; No Tide  Flag    : %i\n", outputImage->noTide);
	fprintf(outputImage->fpLog, "; ThreeD Off Flag    : %i\n", threeDOffFlag);
	fprintf(outputImage->fpLog, "; TimeOverlapFlag    : %i\n", outputImage->timeOverlapFlag);
	fprintf(outputImage->fpLog, "; DeltaB    : %i\n", outputImage->deltaB);
	fprintf(outputImage->fpLog, "; useSquint Flag   : %i\n", useSquint);
	fprintf(outputImage->fpLog, "; noMask Flag      : %i\n", noMask);
	{
		/* Which crossing-orbit solver ran, and its gates.  0 = joint (the default); when joint
		   ran, the rhoPhase/rhoOffsets/pairOverCount factors below do not apply to that round,
		   since a joint solve forms no pairs and cannot over-count. */
		extern int32_t legacyPairPhase, legacyPairRange;
		extern double jointMaxSigma, jointMaxSigmaRange, jointErrScale;
		fprintf(outputImage->fpLog, "; legacyPairPhase  : %i\n", legacyPairPhase);
		fprintf(outputImage->fpLog, "; legacyPairRange  : %i\n", legacyPairRange);
		/*  One cap acts under the hopper (jointMaxSigma, gating the single system that holds
		    phase, range and azimuth rows), so log that alone.  The per-round split is only
		    meaningful when the legacy rounds actually run, so it is logged only then --
		    printing both unconditionally is what made the Phase/Range naming look like a
		    per-observable choice. */
		fprintf(outputImage->fpLog, "; jointMaxSigma    : %f\n", jointMaxSigma);
		if (legacyCode == TRUE)
		{
			fprintf(outputImage->fpLog, "; jointMaxSigRange : %f  (legacy crossing-offsets round)\n",
					jointMaxSigmaRange);
		}
		fprintf(outputImage->fpLog, "; jointErrScale    : %f\n", jointErrScale);
		{
			extern int32_t gateNEff, gateAbsolute;
			extern double maxChi2;
			fprintf(outputImage->fpLog, "; gateNEff         : %i\n", gateNEff);
			fprintf(outputImage->fpLog, "; gateAbsolute     : %i\n", gateAbsolute);
			fprintf(outputImage->fpLog, "; maxChi2          : %f\n", maxChi2);
			fprintf(outputImage->fpLog, "; gateSpeedFrac    : %f\n", gateSpeedFrac);
		}
		{
			extern int32_t noAzimuthRows, aSigmaResidual;
			fprintf(outputImage->fpLog, "; noAzimuthRows    : %i\n", noAzimuthRows);
			fprintf(outputImage->fpLog, "; aSigmaResidual   : %i\n", aSigmaResidual);
		}
		{
			extern int32_t true3D, true3DProject;
			extern double true3DMaxSigma;
			fprintf(outputImage->fpLog, "; true3D           : %i\n", true3D);
			fprintf(outputImage->fpLog, "; true3DProject    : %i\n", true3DProject);
			fprintf(outputImage->fpLog, "; true3DMaxSigma   : %f\n", true3DMaxSigma);
		}
	}
	{
		extern int32_t pairOverCount, rSigmaResidual, pairCountLegacy;
		extern double rSigmaConst, rhoPhase, rhoOffsets;
		fprintf(outputImage->fpLog, "; pairOverCount    : %i\n", pairOverCount);
		fprintf(outputImage->fpLog, "; rSigmaResidual   : %i\n", rSigmaResidual);
		fprintf(outputImage->fpLog, "; rSigmaConst      : %f\n", rSigmaConst);
		fprintf(outputImage->fpLog, "; pairCountLegacy  : %i\n", pairCountLegacy);
		fprintf(outputImage->fpLog, "; rhoPhase         : %f\n", rhoPhase);
		fprintf(outputImage->fpLog, "; rhoOffsets       : %f\n", rhoOffsets);
	}
	fprintf(outputImage->fpLog, "; sepAscDesc Flag  : %i\n", sepAscDesc);
	fprintf(outputImage->fpLog, "; north Flag       : %i\n", args->north);
	fprintf(outputImage->fpLog, "; writeBlank Flag  : %i\n", args->writeBlank);
	fprintf(outputImage->fpLog, "; GTiff Flag       : %i\n", args->GTiff);
	fprintf(outputImage->fpLog, "; COG Flag         : %i\n", args->COG);
	fprintf(outputImage->fpLog, "; outputRA Flag    : %i\n", outputImage->outputRAFlag);
	fprintf(outputImage->fpLog, "; vzFlag           : %i\n", outputImage->vzFlag);
	fprintf(outputImage->fpLog, "; ompThreads       : %d\n", omp_get_max_threads());
	if (refVel->velFile != NULL)
	{
		fprintf(outputImage->fpLog, "; Reference Velocity    : %s\n", refVel->velFile);
		fprintf(outputImage->fpLog, "; ClipFlag    : %f\n", refVel->clipThresh);
		fprintf(outputImage->fpLog, "Include Ref Vel in mosaic %i\n", refVel->initMapFlag);
		fprintf(stderr, "Include Ref Vel in mosaic %i\n", refVel->initMapFlag);
		fprintf(outputImage->fpLog, "Compute errors mode using ref Vel %i\n", refVel->initMapFlag);
	}
	if (outputImage->makeTies == TRUE)
		fprintf(outputImage->fpLog, "; tieThresh        : %lf\n", tieThresh);
	fflush(outputImage->fpLog);
}

static void toSigma(outputImageStructure outputImage)
{
	/*  Variance -> sigma.
	    NOTE: a pixel with a valid velocity and NO valid error is NOT an error state.  The
	    workflow (mosaicworkflow/setupquarters.py, interpMosiacs) gap-fills .vx/.vy but
	    deliberately does NOT interpolate .ex/.ey, so the absence of an error IS the flag that
	    a pixel was interpolated rather than measured.  Do not gate the velocity on the error
	    here -- it would discard exactly the pixels the workflow intends to ship as filled.
	    mosaic3d's own output is already self-consistent (both bands no-data together);
	    the asymmetry appears downstream, by design. */
	for (int i = 0; i < outputImage.ySize; i++)
	{
		for (int j = 0; j < outputImage.xSize; j++)
		{
			if (outputImage.errorX[i][j] > 0.0)
			{
				outputImage.errorX[i][j] = (float)sqrt((double)outputImage.errorX[i][j]);
			}
			if (outputImage.errorY[i][j] > 0.0)
			{
				outputImage.errorY[i][j] = (float)sqrt((double)outputImage.errorY[i][j]);
			}
			if (outputImage.errorZ != NULL && outputImage.errorZ[i][j] > 0.0)
			{
				outputImage.errorZ[i][j] = (float)sqrt((double)outputImage.errorZ[i][j]);
			}
		}
	}
}

static void filterDT(outputImageStructure outputImage)
{
	/* Added 0.49 to ensure last day of month gets included. */
		double maxdT = (outputImage.jd2 - outputImage.jd1 + 1.0) * 0.5 + 0.49; /* Remove .49 ? */
		fprintf(stderr, "maxdT %f %i\n", maxdT, outputImage.timeOverlapFlag);
		if (maxdT < 0 || maxdT > 20000)
			fprintf(stderr, "invalid dT (<0 or > 20000) in writing ouput");
		/* recast pointers */
		float **vx = (float **)outputImage.image;
		float **vy = (float **)outputImage.image2;
		float **dT = (float **)outputImage.image3;
		for (int i = 0; i < outputImage.ySize; i++)
		{
			for (int j = 0; j < outputImage.xSize; j++)
			{
				/* Ensure only data in date range written */
				if (dT[i][j] > maxdT || dT[i][j] < (-maxdT))
				{
					vx[i][j] = -2.e9;
					vy[i][j] = -2.e9;
					dT[i][j] = -2.e9;
					outputImage.errorX[i][j] = -2.e9;
					outputImage.errorY[i][j] = -2.e9;
				}
			}
		}
}

/* Write a multi-band VRT pointing to separate big-endian float32 binary files. */
static void writeBinaryVRT(const char *vrtFile, const char **srcFiles,
                            const char **bandDescriptions, int nBands,
                            int xSize, int ySize, double *geoTransform,
                            const char *epsg, float noDataValue)
{
	GDALDriverH driver = GDALGetDriverByName("VRT");
	GDALDatasetH ds = GDALCreate(driver, vrtFile, xSize, ySize, 0, GDT_Unknown, NULL);
	if (ds == NULL) {
		fprintf(stderr, "writeBinaryVRT: failed to create %s\n", vrtFile);
		return;
	}
	GDALSetGeoTransform(ds, geoTransform);
	if (epsg != NULL) {
		OGRSpatialReferenceH srs = OSRNewSpatialReference(NULL);
		if (OSRImportFromEPSG(srs, atoi(epsg)) == OGRERR_NONE) {
			char *wkt = NULL;
			OSRExportToWkt(srs, &wkt);
			GDALSetProjection(ds, wkt);
			CPLFree(wkt);
		}
		OSRDestroySpatialReference(srs);
	}
	for (int i = 0; i < nBands; i++) {
		char fileOpt[2048];
		sprintf(fileOpt, "SourceFilename=%s", extract_filename((char *)srcFiles[i]));
		char *options[] = {fileOpt, "relativeToVRT=1", "subclass=VRTRawRasterBand",
		                   "BYTEORDER=MSB", NULL};
		GDALAddBand(ds, GDT_Float32, options);
		GDALRasterBandH band = GDALGetRasterBand(ds, i + 1);
		GDALSetDescription(band, bandDescriptions[i]);
		GDALSetRasterNoDataValue(band, noDataValue);
	}
	GDALClose(ds);
}

/* Write VRT sidecars alongside the flat-binary output files. */
static void write3DFlatVRTs(outputImageStructure outputImage, char *outFileBase,
                             const char *outFileVx, const char *outFileVy,
                             const char *outFileVz,
                             const char *outFileEx, const char *outFileEy)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	double geoTransform[6];
	const float noDataValue = -2.0e9;
	char *vrtFile;
	const char *vxDesc, *vyDesc, *vzDesc, *exDesc, *eyDesc;

	computeGeoTransform(geoTransform, outputImage.originX, outputImage.originY,
	                    outputImage.xSize, outputImage.ySize,
	                    outputImage.deltaX, outputImage.deltaY);
	const char *epsg = getEPSGFromProjectionParams(Rotation, SLat, HemiSphere);

	if (outputImage.outputRAFlag) {
		vxDesc = "vr"; vyDesc = "va";
		exDesc = "er"; eyDesc = "ea";
	} else {
		vxDesc = "vx"; vyDesc = "vy";
		exDesc = "ex"; eyDesc = "ey";
	}
	/*  With -timeOverlap the third plane is NOT a vertical velocity: it is "dT", the
	    PRECISION-WEIGHTED MEAN DATE OFFSET of the contributing data from the product's
	    nominal centre date, in days, signed (negative = data skewed early).  It is a
	    weighted FIRST MOMENT, not an interval, a duration, or a span -- an annual product
	    whose data all falls in March and one split evenly between January and June can
	    report the same dT.  Weighted by sqrt(scX*scY), i.e. by how much each measurement
	    actually contributed to the velocity, so it dates the ANSWER rather than the
	    acquisitions.  Also load-bearing, not merely informational: filterDT() below nulls
	    vx/vy/ex/ey wherever |dT| exceeds the half-window.  */
	vzDesc = outputImage.timeOverlapFlag ? "dT" : "vz";

	/* stem.vrt: vx + vy as two named bands */
	vrtFile = appendSuffix(outFileBase, ".vrt", (char *)malloc(strlen(outFileBase) + 5));
	const char *velFiles[] = {outFileVx, outFileVy};
	const char *velDescs[] = {vxDesc, vyDesc};
	writeBinaryVRT(vrtFile, velFiles, velDescs, 2,
	               outputImage.xSize, outputImage.ySize, geoTransform, epsg, noDataValue);
	free(vrtFile);

	/* stem.vz.vrt or stem.dT.vrt */
	const char *vzSuffix = outputImage.timeOverlapFlag ? ".dT.vrt" : ".vz.vrt";
	vrtFile = appendSuffix(outFileBase, (char *)vzSuffix, (char *)malloc(strlen(outFileBase) + 8));
	const char *vzFiles[] = {outFileVz};
	const char *vzDescs[] = {vzDesc};
	writeBinaryVRT(vrtFile, vzFiles, vzDescs, 1,
	               outputImage.xSize, outputImage.ySize, geoTransform, epsg, noDataValue);
	free(vrtFile);

	/* stem.err.vrt: ex + ey as two named bands */
	vrtFile = appendSuffix(outFileBase, ".err.vrt", (char *)malloc(strlen(outFileBase) + 9));
	const char *errFiles[] = {outFileEx, outFileEy};
	const char *errDescs[] = {exDesc, eyDesc};
	writeBinaryVRT(vrtFile, errFiles, errDescs, 2,
	               outputImage.xSize, outputImage.ySize, geoTransform, epsg, noDataValue);
	free(vrtFile);
}

/* Remove stale output files of the OTHER format (binary <-> GeoTIFF) so a
   velocity directory that is re-run with a different output format does not keep
   leftover files of the format no longer being written. Called just before the
   write, so removing the shared ".vrt" here is safe -- the chosen writer rewrites
   it. remove() on a missing file is a harmless no-op, so the stem list can be a
   superset covering the XY/RA (vx/vy vs vr/va) and vz/dT variants. */
static void removeStaleOutputs(char *outFileBase, int32_t writingTiff)
{
	const char *stems[] = {"vx", "vy", "vr", "va", "vz", "dT", "ex", "ey", "er", "ea"};
	const int32_t nStems = (int32_t)(sizeof(stems) / sizeof(stems[0]));
	char path[LINEMAX];
	int32_t i;
	if (writingTiff)
	{
		/* Writing .tif: drop binary flat files, their .geodat sidecars, and the
		   binary VRT sidecars. */
		const char *vrts[] = {".vrt", ".vz.vrt", ".dT.vrt", ".err.vrt"};
		for (i = 0; i < nStems; i++)
		{
			snprintf(path, sizeof(path), "%s.%s", outFileBase, stems[i]);
			remove(path);
			snprintf(path, sizeof(path), "%s.%s.geodat", outFileBase, stems[i]);
			remove(path);
		}
		for (i = 0; i < (int32_t)(sizeof(vrts) / sizeof(vrts[0])); i++)
		{
			snprintf(path, sizeof(path), "%s%s", outFileBase, vrts[i]);
			remove(path);
		}
	}
	else
	{
		/* Writing binary: drop the GeoTIFF counterparts (write3Doutput itself
		   rewrites the shared .vrt and the binary .vz.vrt/.err.vrt). */
		for (i = 0; i < nStems; i++)
		{
			snprintf(path, sizeof(path), "%s.%s.tif", outFileBase, stems[i]);
			remove(path);
		}
	}
}

static void write3Doutput(outputImageStructure outputImage, char *outFileBase)
{
	char *outFileVx, *outFileVy, *outFileVz; /* Output files */
	char *outFileEx, *outFileEy;
	/* Output file names velocity */
	if (outputImage.outputRAFlag)
	{
		outFileVx = appendSuffix(outFileBase, ".vr", (char *)malloc(strlen(outFileBase) + 4));
		outFileVy = appendSuffix(outFileBase, ".va", (char *)malloc(strlen(outFileBase) + 4));
	}
	else
	{
		outFileVx = appendSuffix(outFileBase, ".vx", (char *)malloc(strlen(outFileBase) + 4));
		outFileVy = appendSuffix(outFileBase, ".vy", (char *)malloc(strlen(outFileBase) + 4));
	}
	if (outputImage.timeOverlapFlag == FALSE)
	{
		outFileVz = appendSuffix(outFileBase, ".vz", (char *)malloc(strlen(outFileBase) + 4));
	}
	else
	{
		outFileVz = appendSuffix(outFileBase, ".dT", (char *)malloc(strlen(outFileBase) + 4));
	}
	if (outputImage.outputRAFlag)
	{
		outFileEx = appendSuffix(outFileBase, ".er", (char *)malloc(strlen(outFileBase) + 4));
		outFileEy = appendSuffix(outFileBase, ".ea", (char *)malloc(strlen(outFileBase) + 4));
	}
	else
	{
		outFileEx = appendSuffix(outFileBase, ".ex", (char *)malloc(strlen(outFileBase) + 4));
		outFileEy = appendSuffix(outFileBase, ".ey", (char *)malloc(strlen(outFileBase) + 4));
	}
	// convert error variance to sigma
	toSigma(outputImage);
	//  Fileter by dT if timeOverlap
	if (outputImage.timeOverlapFlag == TRUE)
	{
		filterDT(outputImage);
	}
	outputGeocodedImage(outputImage, outFileVx);
	free(outputImage.image[0]);
	free(outputImage.image);
	outputImage.image = outputImage.image2;
	outputGeocodedImage(outputImage, outFileVy);
	free(outputImage.image[0]);
	free(outputImage.image);
	outputImage.image = outputImage.image3;
	outputGeocodedImage(outputImage, outFileVz);
	free(outputImage.image[0]);
	free(outputImage.image);
	outputImage.image = (void **)outputImage.errorX;
	outputGeocodedImage(outputImage, outFileEx);
	free(outputImage.image[0]);
	free(outputImage.image);
	outputImage.image = (void **)outputImage.errorY;
	outputGeocodedImage(outputImage, outFileEy);
	free(outputImage.image[0]);
	free(outputImage.image);
	/* True-3D vertical velocity and its formal error.  Written only when mosaicTrue3D() ran;
	   deliberately separate bands so the shared vz/dT plane keeps its existing meaning. */
	if (outputImage.vZ3D != NULL)
	{
		char *outFileVz3D = appendSuffix(outFileBase, ".vz3d", (char *)malloc(strlen(outFileBase) + 7));
		outputImage.image = (void **)outputImage.vZ3D;
		outputGeocodedImage(outputImage, outFileVz3D);
		free(outputImage.image[0]);
		free(outputImage.image);
		outputImage.vZ3D = NULL;
	}
	if (outputImage.errorZ != NULL)
	{
		char *outFileEz = appendSuffix(outFileBase, ".ez", (char *)malloc(strlen(outFileBase) + 4));
		outputImage.image = (void **)outputImage.errorZ;
		outputGeocodedImage(outputImage, outFileEz);
		free(outputImage.image[0]);
		free(outputImage.image);
		outputImage.errorZ = NULL;
	}
	/* mosaicHopper3D only: which branch each pixel took (3 = unconstrained 3D, 2 = projected). */
	if (outputImage.hopperMode != NULL)
	{
		char *outFileMode = appendSuffix(outFileBase, ".mode", (char *)malloc(strlen(outFileBase) + 7));
		outputImage.image = (void **)outputImage.hopperMode;
		outputGeocodedImage(outputImage, outFileMode);
		free(outputImage.image[0]);
		free(outputImage.image);
		outputImage.hopperMode = NULL;
	}
	/* Diagnostic band: per-pixel measurement count from the joint solver.
	   Only written when make3DMosaicJoint() actually ran and handed the plane over. */
	if (outputImage.jointNObs != NULL)
	{
		char *outFileNObs = appendSuffix(outFileBase, ".nobs", (char *)malloc(strlen(outFileBase) + 7));
		outputImage.image = (void **)outputImage.jointNObs;
		outputGeocodedImage(outputImage, outFileNObs);
		free(outputImage.image[0]);
		free(outputImage.image);
		outputImage.jointNObs = NULL;
	}
	if (outputImage.jointChi2 != NULL)
	{
		char *outFileChi2 = appendSuffix(outFileBase, ".chi2", (char *)malloc(strlen(outFileBase) + 7));
		outputImage.image = (void **)outputImage.jointChi2;
		outputGeocodedImage(outputImage, outFileChi2);
		free(outputImage.image[0]);
		free(outputImage.image);
		outputImage.jointChi2 = NULL;
	}
	write3DFlatVRTs(outputImage, outFileBase, outFileVx, outFileVy, outFileVz, outFileEx, outFileEy);
}



/*
  Write one diagnostic plane as its own GeoTIFF and free it.

  SEPARATE FILES, not extra bands on the main .vrt: downstream code (mosaicworkflow, the merge
  step, validateMosaic) expects that VRT to be exactly the 5-band vx/vy/vz/ex/ey velocity product,
  and appending bands would break every consumer.  This mirrors the binary path, which writes
  .vz3d/.ez/.nobs/.chi2/.mode as sibling files.

  NULL plane means the solver that owns it never ran, which is what tells us not to write.
  outputImage is passed by value in both writers, so the caller's copy is untouched -- the plane
  is freed here exactly once, same as the binary path does.
*/
static void saveDiagBandTiff(float **plane, char *outFileBase, const char *suffix,
							 double *geoTransform, const char *epsg, dictNode *meta,
							 char *driverType, int32_t xSize, int32_t ySize)
{
	char *outFile;
	float *buf;
	int32_t i;
	if (plane == NULL)
	{
		return;
	}
	/*  MUST copy into a contiguous buffer.  The velocity planes come from mallocOutputImage(),
	    which allocates ONE block per plane and points the rows into it -- which is why
	    saveAsGeotiff(..., outputImage.image[0], ...) is correct for them.  These diagnostic
	    planes come from mallocImage() (common/initRoutines.c:865), which mallocs EVERY ROW
	    SEPARATELY, so plane[0] is a single row, not the image.  Passing it directly reads
	    xSize*ySize floats off the end of one row -- undefined behaviour that happened to look
	    plausible for .vz3d/.ez/.mode/.chi2 (adjacent mallocs land adjacent often enough) and
	    produced obvious garbage for .nobs (median 0, p99 6.5e20).  Do not "optimise" this copy
	    away without first changing how the planes are allocated. */
	buf = (float *)malloc((size_t)xSize * (size_t)ySize * sizeof(float));
	if (buf == NULL)
	{
		error("saveDiagBandTiff: malloc failed for %s\n", suffix);
	}
	for (i = 0; i < ySize; i++)
	{
		memcpy(buf + (size_t)i * (size_t)xSize, plane[i], (size_t)xSize * sizeof(float));
	}
	outFile = appendSuffix(outFileBase, (char *)suffix,
						   (char *)malloc(strlen(outFileBase) + strlen(suffix) + 2));
	saveAsGeotiff(outFile, buf, xSize, ySize, geoTransform, epsg, meta,
				  driverType, GDT_Float32, -2.0e9);
	free(buf);
	/*  Free every row, not just row 0 -- these were malloc'd individually.  (The velocity band
	    writers above free only row 0 plus the pointer array, which is correct for their
	    contiguous allocation.) */
	for (i = 0; i < ySize; i++)
	{
		free(plane[i]);
	}
	free(plane);
}

static void write3DTiffOutput(outputImageStructure outputImage, char *outFileBase, char *driverType, const char *epsg, char *date1,  char *date2)
{
	char *outFileVx, *outFileVy, *outFileVz; /* Output files */
	char *outFileEx, *outFileEy;
	double geoTransform[6];
	float noDataValues[5] = {-2.0e9, -2.0e9, -2.0e9, -2.0e9, -2.0e9};
	/* Output file names velocity */
	if (outputImage.outputRAFlag)
	{
		outFileVx = appendSuffix(outFileBase, ".vr.tif", (char *)malloc(strlen(outFileBase) + 8));
		outFileVy = appendSuffix(outFileBase, ".va.tif", (char *)malloc(strlen(outFileBase) + 8));
	}
	else
	{
		outFileVx = appendSuffix(outFileBase, ".vx.tif", (char *)malloc(strlen(outFileBase) + 8));
		outFileVy = appendSuffix(outFileBase, ".vy.tif", (char *)malloc(strlen(outFileBase) + 8));
	}
	if (outputImage.timeOverlapFlag == FALSE)
	{
		outFileVz = appendSuffix(outFileBase, ".vz.tif", (char *)malloc(strlen(outFileBase) + 8));
	}
	else
	{
		outFileVz = appendSuffix(outFileBase, ".dT.tif", (char *)malloc(strlen(outFileBase) + 8));
	}
	if (outputImage.outputRAFlag)
	{
		outFileEx = appendSuffix(outFileBase, ".er.tif", (char *)malloc(strlen(outFileBase) + 8));
		outFileEy = appendSuffix(outFileBase, ".ea.tif", (char *)malloc(strlen(outFileBase) + 8));
	}
	else
	{
		outFileEx = appendSuffix(outFileBase, ".ex.tif", (char *)malloc(strlen(outFileBase) + 8));
		outFileEy = appendSuffix(outFileBase, ".ey.tif", (char *)malloc(strlen(outFileBase) + 8));
	}
	// convert error variance to sigma
	toSigma(outputImage);
	//  Fileter by dT if timeOverlap
	if (outputImage.timeOverlapFlag == TRUE)
	{
		filterDT(outputImage);
	}
    // Compute Geotransform
	computeGeoTransform(geoTransform, outputImage.originX, outputImage.originY, outputImage.xSize,
					    outputImage.ySize, outputImage.deltaX, outputImage.deltaY);
	// Set up meta data
	dictNode *summaryMetaData = NULL;
	char *timeStamp = timeStampMeta();
	insert_node(&summaryMetaData, "FirstDate", date1);
	insert_node(&summaryMetaData, "LastDate", date2);
	insert_node(&summaryMetaData, "CreationTime", timeStamp);
	// Get epsg code
	if(epsg == NULL)
	{
		epsg = getEPSGFromProjectionParams(Rotation, SLat, HemiSphere); 
	}	
	// Save files as tiffs
	// Vx
	saveAsGeotiff(outFileVx, (float *)outputImage.image[0], outputImage.xSize,
				outputImage.ySize, geoTransform, epsg, summaryMetaData, driverType, GDT_Float32, noDataValues[0]);
	free(outputImage.image[0]);
	free(outputImage.image);
	// Vy
	saveAsGeotiff(outFileVy, (float *)outputImage.image2[0], outputImage.xSize,
				 outputImage.ySize, geoTransform, epsg, summaryMetaData, driverType, GDT_Float32, noDataValues[1]);
	free(outputImage.image2[0]);
	free(outputImage.image2);
	// Vz
	saveAsGeotiff(outFileVz, (float *)outputImage.image3[0], outputImage.xSize,
				 outputImage.ySize, geoTransform, epsg, summaryMetaData, driverType, GDT_Float32, noDataValues[2]);
	free(outputImage.image3[0]);
	free(outputImage.image3);
	// Ex
	saveAsGeotiff(outFileEx, (float *)outputImage.errorX[0], outputImage.xSize,
				 outputImage.ySize, geoTransform, epsg, summaryMetaData, driverType, GDT_Float32, noDataValues[3]);
	free(outputImage.errorX[0]);
	free(outputImage.errorX);
	// Ey
	saveAsGeotiff(outFileEy, (float *)outputImage.errorY[0], outputImage.xSize,
				 outputImage.ySize, geoTransform, epsg, summaryMetaData, driverType, GDT_Float32, noDataValues[4]);
	free(outputImage.errorY[0]);
	free(outputImage.errorY);	
	//
	// Create VRT
	char *vrtFile = appendSuffix(outFileBase, ".vrt", (char *)malloc(strlen(outFileBase) + 5));
	const char *bands[] = {outFileVx, outFileVy, outFileVz, outFileEx, outFileEy};
	makeTiffVRT(vrtFile, bands, 5, noDataValues, summaryMetaData);
	/*
	  Diagnostic bands from the 3D and joint solvers.  Each is NULL unless the solver that owns it
	  actually ran, so a plain phase or offsets run writes none of them and is byte-identical to
	  before.  Until 2026-08-30 these were written on the binary path only, so every -GTiff/-COG
	  run -- i.e. all of production -- silently dropped them, which made the mosaicHopper3D 2D/3D
	  split unmeasurable in exactly the products people look at.
	*/
	saveDiagBandTiff(outputImage.vZ3D, outFileBase, ".vz3d.tif", geoTransform, epsg,
					 summaryMetaData, driverType, outputImage.xSize, outputImage.ySize);
	saveDiagBandTiff(outputImage.errorZ, outFileBase, ".ez.tif", geoTransform, epsg,
					 summaryMetaData, driverType, outputImage.xSize, outputImage.ySize);
	saveDiagBandTiff(outputImage.hopperMode, outFileBase, ".mode.tif", geoTransform, epsg,
					 summaryMetaData, driverType, outputImage.xSize, outputImage.ySize);
	saveDiagBandTiff(outputImage.jointNObs, outFileBase, ".nobs.tif", geoTransform, epsg,
					 summaryMetaData, driverType, outputImage.xSize, outputImage.ySize);
	saveDiagBandTiff(outputImage.jointChi2, outFileBase, ".chi2.tif", geoTransform, epsg,
					 summaryMetaData, driverType, outputImage.xSize, outputImage.ySize);
}

static void processMosaicDate(outputImageStructure *outputImage, char *date1, char *date2)
{
	int32_t m1, m2, d1, d2, y1, y2;

	fprintf(stderr, "date1,date2 %s %s\n", date1, date2);
	if (date1 != NULL)
	{
		if (sscanf(date1, "%2d-%2d-%4d", &m1, &d1, &y1) != 3)
			error("invalid date %s\n", date1);
		outputImage->jd1 = juldayDouble(m1, d1, y1);
		fprintf(stderr, "%i %i %i %lf\n", m1, d1, y1, outputImage->jd1);
	}
	else
		outputImage->jd1 = 0.0;

	if (date2 != NULL)
	{
		if (sscanf(date2, "%2d-%2d-%4d", &m2, &d2, &y2) != 3)
			error("invalid date %s\n", date2);
		outputImage->jd2 = (double)juldayDouble(m2, d2, y2) + 0.9999999; /* add 0.999 to ensure end of day */
		fprintf(stderr, "%i %i %i %lf\n", m2, d2, y2, outputImage->jd2);
	}
	else
		outputImage->jd2 = HIGHJD + 1000.; /* way past present */

	if ((date1 != NULL || date2 != NULL) && (outputImage->jd2 < outputImage->jd1))
		error("date 1 follows date2\n");
	fprintf(stderr, "date1,date2 %f %f\n", outputImage->jd1, outputImage->jd2);
}

static int32_t writeMetaFile(inputImageStructure *image, outputImageStructure *outputImage, vhParams *params, char *outFileBase, char *demFile, int32_t writeBlank)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	FILE *fpMeta;
	char *metaFile;
	char sensor[128];
	char *tmp;
	int32_t mm, day, year;
	double dday;
	time_t tproc;
	struct tm *tS, tTmp;
	double x, y, lat, lon;
	size_t sMax = 128;
	char *timeString;
	int32_t hour, minute, second;
	inputImageStructure *imageTmp;
	int32_t haveData;
	/*
	  If single pair result, make a meta file with the date
	*/
	metaFile = (char *)malloc(LINEMAX);
	metaFile[0] = '\0';
	metaFile = strcat(metaFile, outFileBase);
	metaFile = strcat(metaFile, ".meta");
	fpMeta = fopen(metaFile, "w");
	timeString = (char *)malloc(sizeof(char) * 128);
	haveData = writeBlank; /* If writeBlank this will force have data to true */
	for (imageTmp = image; imageTmp != NULL; imageTmp = imageTmp->next)
	{
		if (imageTmp->used == TRUE)
		{ /* Only output images that were used */
			haveData = TRUE;
			/*
				Julian date for central estimate
			*/
			fprintf(fpMeta, "Central Julian Date (CE) for Pair = %10.3lf\n", imageTmp->julDay + params->nDays / 2);
			/*
			  Date for first image
			*/
			caldat((int32_t)(imageTmp->julDay), &mm, &day, &year);
			fprintf(stderr, "mm,day,year %i %i %i %f\n", mm, day, year, imageTmp->julDay);
			/* note uses month as 0 to 11 */
			tTmp.tm_year = year - 1900;
			tTmp.tm_mon = mm - 1;
			tTmp.tm_mday = day;
			strftime(timeString, sMax, "%b:%d:%Y", &tTmp);
			fprintf(fpMeta, "First Image Date (MM:DD:YYYY) = %s\n", timeString);
			/*
			  Date for first second image
			*/
			caldat((int32_t)(imageTmp->julDay + params->nDays), &mm, &day, &year);
			tTmp.tm_year = year - 1900;
			tTmp.tm_mon = mm - 1;
			tTmp.tm_mday = day;
			strftime(timeString, sMax, "%b:%d:%Y", &tTmp);
			fprintf(fpMeta, "Second Image Date (MM:DD:YYYY) = %s\n", timeString);
			/*
			  Nominal time
			*/
			dday = (imageTmp->julDay - (int32_t)imageTmp->julDay);
			hour = (int)(dday * 24);
			minute = (int)((dday * 24 - hour) * 60);
			second = (int)(dday * 86400 - hour * 3600 - minute * 60 + .5);
			tTmp.tm_sec = second;
			tTmp.tm_hour = hour;
			tTmp.tm_min = minute;
			strftime(timeString, sMax, "%H:%M:%S", &tTmp);
			fprintf(fpMeta, "Nominal Time for Pair (HH:MM:SS) = %s\n", timeString);
			/*
			  Sensor
			*/
			sensor[0] = '\0';
			tmp = imageTmp->par.label;
			/*
			  Extra sensor from geodat title. Only implemented now for TSX/TDX
			*/
			if (strstr(tmp, "TDX") != NULL)
				strcat(sensor, "TSX/TDX");
			else if (strstr(imageTmp->par.label, "TSX") != NULL)
				strcat(sensor, "TSX/TDX");
			else
				strcat(sensor, "not_specfied");
			fprintf(fpMeta, "Sensor = %s\n", sensor);
		}

	} /* End for image loop */
	x = outputImage->originX + 0.5 * ((outputImage->xSize - 1) * outputImage->deltaX);
	x *= 0.001;
	y = outputImage->originY + 0.5 * ((outputImage->ySize - 1) * outputImage->deltaY);
	y *= 0.001;
	xytoll1(x, y, HemiSphere, &lat, &lon, Rotation, outputImage->slat);
	if (lon > 180)
		lon -= 360;
	fprintf(fpMeta, "Product Center Latitude  = %10.5lf\n", lat);
	fprintf(fpMeta, "Product Center Longitude = %10.5lf\n", lon);
	/*
		DEM
	*/
	if (strstr(demFile, "gimp1") != NULL)
	{
		fprintf(fpMeta, "DEM version = GIMP DEM V1\n");
	}
	else
	{
		if (strstr(demFile, "gimp2") != NULL)
			fprintf(fpMeta, "DEM version = GIMP DEM V2\n");
		else
			fprintf(fpMeta, "DEM version = Custom\n");
	}
	/*
	  time
	*/
	tproc = time(&tproc);
	tS = localtime(&tproc);
	strftime(timeString, sMax, "%b-%d-%Y-%H:%M:%S\n", tS);
	fprintf(fpMeta, "Production Date/Time = %s\n", timeString);
	if (haveData == FALSE)
		remove(metaFile);
	return (haveData);
}

static void readArgs(int32_t argc, char *argv[], mosaicArgs *args,
					 referenceVelocity *refVel, outputImageStructure *outputImage)
{
	extern int32_t sepAscDesc;
	char *argString;
	float tmp;
	int32_t i, n;
	int32_t noVhFlag, no3d, rOffsetFlag, vzFlag, noTide, timeOverlapFlag;
	/*  Which legacy round-selection flags were EXPLICITLY given.  The translation below needs
	    "was it passed", not "what is its value" -- with none passed the derivation must not fire
	    at all, or a bare run would silently lose its range rows (range is ON only when -3dOff or
	    -rOffsets is present, whereas the modern default is all three row types). */
	int32_t gaveNo3d = FALSE, gaveNoVh = FALSE, gave3dOff = FALSE, gaveROffsets = FALSE;
	int32_t gaveRowFlag = FALSE; /* an explicit -no*Rows always beats a derived value */
	int32_t deltaB;
	char *verticalCorrectionSuffix;
	int32_t iceOnly;
	int32_t flipSquint;

	/* Raised 30 -> 44 (2026-08-26) when -rSigmaConst, which takes a VALUE and so
	   costs two argv slots, pushed a routine production command over the limit --
	   it failed by printing the usage text, which reads as an unrecognised flag
	   rather than an arg-count overflow.  Every optional flag added over the years
	   has eaten into this headroom; the message had also gone stale at 27.
	   Raised again 44 -> 52 for -rhoPhase/-rhoOffsets, two more value-taking
	   flags at two argv slots each.  Raised 52 -> 62 for the joint-solver flags:
	   -legacyPairPhase/-legacyPairRange (one each), -jointMaxSigma/-jointMaxSigmaRange
	   and -jointErrScale (two each), plus a slot of headroom.  Raised 62 -> 70 for the joint
	   speckle solver: -noAzimuthRows, -noASigmaResidual (one each), the retired-but-accepted
	   -speckleTrackJoint (one) and -jointMaxSigmaSpeckle (two), plus headroom.  Raised 70 -> 76 for -true3D /
	   -true3DProject (one each) and -true3DMaxSigma (two), plus headroom.  Raised 76 -> 80 for
	   -hopper3D (one) and -hopper3DMaxSigma (two), plus headroom. */
	if (argc < 4 || argc > 84)
	{
		fprintf(stderr, "Arg count out of range (max 80): %i\n", argc);
		usage(); /* Check number of args */
	}
	n = argc - 4;
	/*
	  Defaults
	*/
	refVel->clipFlag = FALSE;
	refVel->clipThresh = 100000;
	refVel->initMapFlag = FALSE;
	refVel->velFile = NULL;
	args->date1 = NULL;
	args->date2 = NULL;
	noTide = FALSE;
	outputImage->makeTies = FALSE;
	args->irregFile = NULL;
	rOffsetFlag = FALSE;
	args->shelfMaskFile = NULL;
	args->extraTieFile = NULL;
	args->landSatFile = NULL;
	noVhFlag = FALSE;
	timeOverlapFlag = FALSE;
	no3d = FALSE;
	args->threeDOffFlag = FALSE;
	args->tieThresh = 100.0; /* Limit on velocity for tie points */
	args->fl = 0.0;
	args->tideFile = NULL;
	args->north = FALSE;
	args->timeThresh = 12;
	args->timeThreshPhase = 548; /* Allow pairing with in 1.5 years */
	args->sigmaAThresh = 1000.0;
	args->writeBlank = FALSE;
	vzFlag = VZDEFAULT;
	args->statsFlag = FALSE;
	args->verticalCorrectionFile = NULL;
	verticalCorrectionSuffix = NULL;
	iceOnly = FALSE;
	flipSquint = FALSE;
	args->COG = FALSE;
	args->GTiff = FALSE;
	args->outputRAFlag = FALSE;
	/* Added this flag to sort ignore crossing orbits or like asc/desc types May 6 2014 */
	sepAscDesc = TRUE;
	deltaB = DELTABNONE;
	/*
	  Parse command line args
	*/
	for (i = 1; i <= n; i++)
	{
		argString = strchr(argv[i], '-');
		if (strstr(argString, "center") != NULL)
			fprintf(stderr, "Ignoring obsolete center flag\n");
		else if (strstr(argString, "xyDEM") != NULL)
			fprintf(stderr, "xyDEM flag obsolete - xydem is the default");
		else if (strstr(argString, "writeBlank") != NULL)
			args->writeBlank = TRUE;
		/* Tested early on purpose.  This parser dispatches on strstr, so a flag is silently
		   swallowed by any SHORTER flag string contained in it that is tested first.  Longest
		   first within each family: -jointMaxSigma/-jointMaxSigmaRange share the prefix
		   -jointMaxSigma, and -legacyPairPhase/-legacyPairRange share -legacyPair. */
		else if (strstr(argString, "legacyPairPhase") != NULL)
		{
			legacyPairPhase = TRUE;
		}
		else if (strstr(argString, "legacyPairRange") != NULL)
		{
			legacyPairRange = TRUE;
		}
		/*  Legacy alias: sets the primary cap only, leaving jointMaxSigmaRange at its own
		    value -- that is what preserves pre-2026-09 behaviour for templates that pass it.
		    MUST precede the bare -jointMaxSigma, which is a substring of it. */
		else if (strstr(argString, "jointMaxSigmaPhase") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &jointMaxSigma);
			i++;
		}
		else if (strstr(argString, "jointMaxSigmaRange") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &jointMaxSigmaRange);
			i++;
		}
		/*  RETIRED 2026-09: -jointMaxSigmaSpeckle went with speckleTrackMosaicJoint.  Accepted
		    and ignored (it takes a value, so the value must still be consumed) rather than
		    rejected, so an old script keeps running.  MUST precede the bare -jointMaxSigma
		    below, which is a substring of it. */
		else if (strstr(argString, "jointMaxSigmaSpeckle") != NULL)
		{
			fprintf(stderr, "; note: -jointMaxSigmaSpeckle is retired and ignored "
							"(the joint speckle solver was removed); use -jointMaxSigma\n");
			i++;
		}
		/*  -jointMaxSigma: the canonical name.  The Phase/Range split is a legacy-pipeline
		    artifact -- under the hopper there is ONE system holding phase, range and azimuth
		    rows and only jointMaxSigma was ever read, so naming it "Phase" misdescribes
		    what is gated.  This sets both, which is the same thing under the hopper and the
		    obvious meaning under -legacyCode.  Tested AFTER the longer names above: strstr
		    dispatch means the bare prefix would otherwise swallow all of them. */
		else if (strstr(argString, "jointMaxSigma") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &jointMaxSigma);
			jointMaxSigmaRange = jointMaxSigma;
			i++;
		}
		else if (strstr(argString, "jointErrScale") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &jointErrScale);
			i++;
		}
		/*  RETIRED 2026-09: speckleTrackMosaicJoint was removed.  Its purpose was to keep the
		    full 2x2 covariance that speckleTrackMosaic isotropises, and the hopper does that
		    by construction while also pooling rows across images.  Accepted and ignored so an
		    old invocation still runs -- it now simply gets the hopper. */
		else if (strstr(argString, "speckleTrackJoint") != NULL)
		{
			fprintf(stderr, "; note: -speckleTrackJoint is retired and ignored -- the hopper "
							"supersedes it (use -noPhaseRows for a speckle-only solve)\n");
		}
		else if (strstr(argString, "noAzimuthRows") != NULL)
		{
			noAzimuthRows = TRUE;
			gaveRowFlag = TRUE;
		}
		/* Not swallowed by the rSigmaResidual tests further down: those require a lowercase
		   'r' immediately before "SigmaResidual", and this has an 'A' there. */
		else if (strstr(argString, "noASigmaResidual") != NULL)
		{
			aSigmaResidual = FALSE;
		}
		/* Longest first: "true3DMaxSigma" and "true3DProject" both contain "true3D". */
		else if (strstr(argString, "true3DMaxSigma") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &true3DMaxSigma);
			i++;
		}
		else if (strstr(argString, "noErrorGate") != NULL)
		{
			noErrorGate = TRUE;
		}
		else if (strstr(argString, "sigmaAThreshVel") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &sigmaAThreshVel);
			i++;
		}
		else if (strstr(argString, "maxChi2") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &maxChi2);
			i++;
		}
		else if (strstr(argString, "noGateAbsolute") != NULL)
		{
			gateAbsolute = FALSE;
		}
		else if (strstr(argString, "gateSpeedFrac") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &gateSpeedFrac);
			i++;
		}
		else if (strstr(argString, "gateAbsolute") != NULL)
		{
			gateAbsolute = TRUE;
		}
		else if (strstr(argString, "noGateNEff") != NULL)
		{
			gateNEff = FALSE;
		}
		else if (strstr(argString, "gateNEff") != NULL)
		{
			gateNEff = TRUE; /* accepted for symmetry; already the default */
		}
		else if (strstr(argString, "legacyCode") != NULL)
		{
			legacyCode = TRUE;
			hopper = FALSE;
		}
		else if (strstr(argString, "noPhaseRows") != NULL)
		{
			noPhaseRows = TRUE;
			gaveRowFlag = TRUE;
		}
		else if (strstr(argString, "noRangeRows") != NULL)
		{
			noRangeRows = TRUE;
			gaveRowFlag = TRUE;
		}
		else if (strstr(argString, "obsDump") != NULL)
		{
			obsDumpFile = argv[i + 1];
			i++;
		}
		/*  LONGEST FIRST.  "hopper3DMaxSigma" and "hopper3D" both contain "hopper", and this
		    parser dispatches on strstr -- so a bare "hopper" test placed first would swallow
		    -hopper3D (silently running the 2D module) and would leave -hopper3DMaxSigma's value
		    token to fall through to usage().  Same trap the root CLAUDE.md records for
		    -rhoOffsets vs offsets. */
		else if (strstr(argString, "hopper3DMaxSigma") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &hopper3DMaxSigma);
			i++;
		}
		else if (strstr(argString, "hopper3D") != NULL)
		{
			hopper3D = TRUE;
		}
		else if (strstr(argString, "hopper") != NULL)
		{
			hopper = TRUE;
		}
		else if (strstr(argString, "true3DPhase") != NULL)
		{
			true3DPhase = TRUE;
		}
		else if (strstr(argString, "true3DDiag") != NULL)
		{
			true3DDiag = TRUE;
			true3D = TRUE;
		}
		else if (strstr(argString, "true3DProject") != NULL)
		{
			true3DProject = TRUE;
			true3D = TRUE;
		}
		else if (strstr(argString, "true3D") != NULL)
		{
			true3D = TRUE;
		}
		/* The rho flags MUST be tested before "rOffsets"/"offsets" below.  This
		   parser dispatches on strstr, so a longer flag containing a shorter
		   one is silently swallowed by the shorter test.  "-rhoOffsets" escapes
		   "offsets" only because of the capital O, which is far too fragile a
		   thing to depend on. */
		else if (strstr(argString, "rhoPhase") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &rhoPhase);
			i++;
		}
		else if (strstr(argString, "rhoOffsets") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &rhoOffsets);
			i++;
		}
		else if (strstr(argString, "rOffsets") != NULL)
		{
			rOffsetFlag = TRUE;
			gaveROffsets = TRUE;
		}
		else if (strstr(argString, "offsets") != NULL)
			fprintf(stderr, "ignoring obsolet offsets flag - always enabled\n");
		else if (strstr(argString, "-initMap") != NULL)
			refVel->initMapFlag = TRUE;
		else if (strstr(argString, "makeTies") != NULL)
			outputImage->makeTies = TRUE;
		else if (strstr(argString, "timeOverlap") != NULL)
			timeOverlapFlag = TRUE;
		else if (strstr(argString, "stats") != NULL)
			args->statsFlag = TRUE;
		else if (strstr(argString, "COG") != NULL)
			args->COG = TRUE;
		else if (strstr(argString, "GTiff") != NULL)
			args->GTiff = TRUE;
		else if (strstr(argString, "vzFlag") != NULL)
		{
			if (sscanf(argv[i + 1], "%i\n", &vzFlag) != 1)
				usage();
			i++;
		}
		else if (strstr(argString, "sigmaAThresh") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &args->sigmaAThresh);
			i++;
		}
		else if (strstr(argString, "tieThresh") != NULL)
		{
			if (sscanf(argv[i + 1], "%lf\n", &args->tieThresh) != 1)
				usage();
			i++;
		}
		else if (strstr(argString, "extraTies") != NULL)
		{
			args->extraTieFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "tideFile") != NULL)
		{
			args->tideFile = argv[i + 1];
			i++;
		} // Make sure this goes before shorter verticalCorrection
		else if (strstr(argString, "verticalCorrectionSuffix") != NULL)
		{
			verticalCorrectionSuffix = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "iceOnly") != NULL)
		{
			iceOnly = TRUE;
		}
		/* Must be tested BEFORE "useSquint": strstr would match the shorter flag
		   inside neither, but keep them adjacent so the ordering stays obvious. */
		else if (strstr(argString, "flipSquint") != NULL)
		{
			extern int32_t flipSquintSign;
			flipSquintSign = TRUE;
			flipSquint = TRUE;
		}
		else if (strstr(argString, "verticalCorrection") != NULL)
		{
			args->verticalCorrectionFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "useSquint") != NULL)
		{
			useSquint = TRUE;
		}
		/* longest-first: -noPairOverCount must be tested before -pairOverCount */
		else if (strstr(argString, "noPairOverCount") != NULL)
		{
			pairOverCount = FALSE;
		}
		else if (strstr(argString, "pairOverCount") != NULL)
		{
			pairOverCount = TRUE;
		}
		else if (strstr(argString, "pairCountLegacy") != NULL)
		{
			pairCountLegacy = TRUE;
		}
		/* longest-first, as above */
		else if (strstr(argString, "noRSigmaResidual") != NULL)
		{
			rSigmaResidual = FALSE;
		}
		else if (strstr(argString, "rSigmaResidual") != NULL)
		{
			rSigmaResidual = TRUE;
		}
		else if (strstr(argString, "rSigmaConst") != NULL)
		{
			sscanf(argv[i + 1], "%lf", &rSigmaConst);
			i++;
		}
		else if (strstr(argString, "noMask") != NULL)
		{
			noMask = TRUE;
		}
		else if (strstr(argString, "fl") != NULL)
		{
			sscanf(argv[i + 1], "%f", &args->fl);
			i++;
		}
		else if (strstr(argString, "timeThresh") != NULL)
		{
			sscanf(argv[i + 1], "%f", &args->timeThresh);
			i++;
		}
		else if (strstr(argString, "timePhaseThresh") != NULL)
		{
			sscanf(argv[i + 1], "%f", &args->timeThreshPhase);
			i++;
		}
		else if (strstr(argString, "irreg") != NULL)
		{
			args->irregFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "noVh") != NULL)
		{
			noVhFlag = TRUE;
			gaveNoVh = TRUE;
		}
		else if (strstr(argString, "no3d") != NULL)
		{
			no3d = TRUE;
			gaveNo3d = TRUE;
		}
		else if (strstr(argString, "3dOff") != NULL)
		{
			args->threeDOffFlag = TRUE;
			gave3dOff = TRUE;
		}
		else if (strstr(argString, "noSepAscDesc") != NULL)
		{
			sepAscDesc = FALSE;
		}
		else if (strstr(argString, "noTide") != NULL)
		{
			noTide = TRUE;
		}
		else if (strstr(argString, "SVAlongTrack") != NULL)
		{
			deltaB = DELTABQUAD;
		}
		else if (strstr(argString, "SVConst") != NULL)
		{
			deltaB = DELTABCONST;
		}
		else if (strstr(argString, "outputRA") != NULL)
		{
			args->outputRAFlag = TRUE;
		}
		else if (strstr(argString, "north") != NULL)
		{
			args->north = TRUE;
		}
		else if (strstr(argString, "landSat") != NULL)
		{
			args->landSatFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "refVel") != NULL)
		{
			refVel->velFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "clipThresh") != NULL)
		{
			sscanf(argv[i + 1], "%f", &tmp);
			i++;
			refVel->clipFlag = TRUE;
			refVel->clipThresh = tmp;
		}
		else if (strstr(argString, "shelfMask") != NULL)
		{
			args->shelfMaskFile = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "date1") != NULL)
		{
			args->date1 = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "date2") != NULL)
		{
			args->date2 = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "ompThreads") != NULL)
		{
			int32_t nThreads = 0;
			if (i + 1 < argc && argv[i + 1][0] != '-' && argv[i + 1][0] != '\0')
			{
				sscanf(argv[i + 1], "%d", &nThreads);
				i++;
			}
			if (nThreads > 0)
			{
				omp_set_num_threads(nThreads);
				fprintf(stderr, "\033[1;3;34mompThreads set to %d\033[0m\n", nThreads);
			}
			else
				fprintf(stderr, "\033[1;3;34mompThreads using default (%d)\033[0m\n", omp_get_max_threads());
		}
		else
			usage();
	}
	if (refVel->initMapFlag == TRUE && refVel->velFile == NULL)
		error("initMap set but no reference velocity provided\n");
	if (args->statsFlag == TRUE)
	{
		fprintf(stderr, "setting flags to force speckletracking only since in computeError or stats mode");
		no3d = TRUE;
		args->threeDOffFlag = FALSE;
		noVhFlag = TRUE;
		rOffsetFlag = TRUE;
	}
	if (args->statsFlag == TRUE)
	{
		printf("stats and timeOverlap flags incompatible, setting timeOverlap flag to False");
		timeOverlapFlag = FALSE;
	}
	/*  ---- Solver selection, and translation of the legacy round-selection flags -------------
	    -stats is a speckle-round OUTPUT mode (unweighted mean/sigma/count) with no hopper
	    equivalent, and the legacyPair* solvers are legacy by definition, so all three force the
	    old pipeline rather than being silently ignored.  setupquarters.py:230 can append
	    -legacyPairPhase/-legacyPairRange onto ANY template at runtime, including hopper ones, so
	    that combination is reachable even though it appears in no file on disk. */
	if (args->statsFlag == TRUE || legacyPairPhase == TRUE || legacyPairRange == TRUE)
	{
		if (legacyCode == FALSE)
		{
			fprintf(stderr, "; note: -stats/-legacyPair* have no hopper equivalent -- using -legacyCode\n");
		}
		legacyCode = TRUE;
		hopper = FALSE;
	}
	if (legacyCode == TRUE)
	{
		hopper = FALSE;
		hopper3D = FALSE;
	}
	else
	{
		/*  Translate the legacy flags into row switches so an existing pair template runs
		    unchanged and measures what it always measured.  The three conditions come from which
		    rounds consume which observable:
		        phase   -- crossing round (off under -no3d) OR vh round (off under -noVh)
		        range   -- crossing-offsets round (-3dOff) OR speckle round (-rOffsets)
		        azimuth -- vh round (on unless -noVh) OR speckle round (-rOffsets)
		    setup3D.c:319 states the phase disjunction itself.  Verified against all eight flag
		    combinations in production use; see Documents/mosaic3d.md. */
		if (gaveNo3d || gaveNoVh || gave3dOff || gaveROffsets)
		{
			int32_t wantPhase = !(no3d == TRUE && noVhFlag == TRUE);
			int32_t wantRange = (args->threeDOffFlag == TRUE || rOffsetFlag == TRUE);
			int32_t wantAz = (noVhFlag == FALSE || rOffsetFlag == TRUE);
			char derived[256];
			derived[0] = '\0';
			if (wantPhase == FALSE && noPhaseRows == FALSE)
			{
				noPhaseRows = TRUE;
				strcat(derived, " -noPhaseRows");
			}
			if (wantRange == FALSE && noRangeRows == FALSE)
			{
				noRangeRows = TRUE;
				strcat(derived, " -noRangeRows");
			}
			if (wantAz == FALSE && noAzimuthRows == FALSE)
			{
				noAzimuthRows = TRUE;
				strcat(derived, " -noAzimuthRows");
			}
			fprintf(stderr, "; DEPRECATED legacy flags:%s%s%s%s\n",
					gaveNo3d ? " -no3d" : "", gaveNoVh ? " -noVh" : "",
					gave3dOff ? " -3dOff" : "", gaveROffsets ? " -rOffsets" : "");
			fprintf(stderr, "; translated for the hopper as:%s\n",
					derived[0] != '\0' ? derived : " (all row types kept)");
			if (gaveRowFlag == TRUE)
			{
				fprintf(stderr, "; explicit -no*Rows flags were also given and take precedence\n");
			}
			fprintf(stderr, "; pass the -no*Rows flags directly to silence this, or -legacyCode "
							"for the old pipeline\n");
		}
		/*  setup3D only parses the range-offset files when -rOffsets or -3dOff is set
		    (setup3D.c:382-400).  The hopper consumes them as rows, so they must always be
		    loaded; row selection is then handled by noRangeRows above, not by whether the file
		    was read.  Set AFTER the translation, which needs the as-given value. */
		rOffsetFlag = TRUE;
	}
	/*  Skip the azimuth-offset RASTER READ when nothing will consume it.  Requires BOTH
	    -noAzimuthRows (so the hopper suppresses the row) and the hopper actually running:
	    under -legacyCode, makeVhMosaic reads offsets->da gated only on its own offsetFlag and
	    would consume an unread buffer.  Pure I/O saving; the row gating is unchanged. */
	{
		extern int32_t skipAzimuthOffsets;
		skipAzimuthOffsets = (noAzimuthRows == TRUE && legacyCode == FALSE);
		if (skipAzimuthOffsets == TRUE)
		{
			fprintf(stderr, "; -noAzimuthRows with the hopper: skipping the azimuth offset read\n");
		}
	}
	if (args->COG == TRUE && args->GTiff == TRUE)
	{
		error("Select COG or GTiff but not both");
	}
	/* Must have az offsets to do range offsets */
	args->inputFile = argv[argc - 3];
	args->demFile = argv[argc - 2];
	args->outFileBase = argv[argc - 1];
	outputImage->noVhFlag = noVhFlag;
	outputImage->no3d = no3d;
	outputImage->rOffsetFlag = rOffsetFlag;
	outputImage->noTide = noTide;
	outputImage->vzFlag = vzFlag;
	outputImage->deltaB = deltaB;
	outputImage->outputRAFlag = args->outputRAFlag;
	outputImage->timeOverlapFlag = timeOverlapFlag;
	outputImage->verticalCorrectionSuffix = verticalCorrectionSuffix;
	outputImage->iceOnly = iceOnly;
	outputImage->flipSquint = flipSquint;
	return;
}

static void usage()
{
	error("\033[1m\n\n%s\n\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\033[0m\n\n\n",
		  "mosaic3d: mosaic phase and speckle data to create a velocity mosaic",
		  "Usage:",
		  " mosaic3d -north -GTiff -COG -writeBlank -makeTies -tieThresh -extraTies extraTieFile -date1 MM-DD-YYYY -date2 MM-DD-YYYY -timeOverlap -tideFile tideFile \\",
		  " \t-verticalCorrection vcFile -verticalCorrectionSuffix suffix -no3d -3dOff -deltaBQ -deltaBC -noVh -landSat landSatList  -refVel refVelFile -initMap -clipThresh clipThresh -lsClip clipVal  \\",
		  " \t-shelfMask shelfMask -fl fl  -rOffsets -timeThresh timeThresh -irregFile irregFile -vzFlag flag -stats -noSepAscDesc -noTide -ompThreads N \\",
		  " \tinputFile demFile outPutImage\n",
		  " where :",
		  "\tnorth =\t\t\t Force northern hemisphere",
		  "\tGTiff =\t\t\t Save to geotiff files",
		  "\tCOG =\t\t\t Save to Cloud Optimized Geotiff",
		  "\twriteBlank =\t\t\t Force writing of results with no data",
		  "\tmakeTies =\t\t Flag to produce tiepoint file instead of usual output ",
		  "\ttieThresh =\t\t Only use values less than tieThresh for tiepoints",
		  "\textraTies =\t\t File with additional ties to include with data file",
		  "\tdate1,date2 =\t\t Time range for data (include only interval in range (date1,date2). format -> MM-DD-YYYY",
		  "\ttimeOverlap =\t\t Default is False so only images full in (date1,date2) are included. \n\t\t\t\t If set,  images that overlap (date1,date2) are included and weighted accordingly",
		  "\ttideFile =\t\t Tidal correction used for creating tiepoints",
		  "\tverticalCorrection =\t\t Vertical velocity correction (sub/emerg vel) in m/yr",
		  "\tverticalCorrectionSuffix =\t\t Suffix to add to phase, baseline, and rBaseline files for  velocity correction (sub/emerg vel) in m/yr",
		  "\tno3d =\t\t\t No crossing orbit solution with phase",
		  "\t3dOff =\t\t\t Crossing orbit solution with offsets",
		  "\tSVAlongTrack =\t\tIf available, use along track correction to  state vector solution",
		  "\tSVConst =\t\t If available, use const only correction on top of state vector solution",
		  "\tnoVh =\t\t\t Use only the crossing orbit solution",
		  "\tlandSat =\t\t File containing list of landsat offsets to include in mosaic",
		  "\trefVel =\t\t Velocity mosaic used to clip errors ",
		  "\tinitMap =\t\t Interpolate refVel as startin point for mosaic ",
		  "\tclipThresh =\t\t Use with refVel to clip differences > clipThresh for slow moving regions (< 100 m/yr).  Honoured by the hopper solvers as well as the legacy speckle round",
		  "\tlsClip (not implemented)=\t\t Clip landsat mosaic using refVelFile ",
		  "\tshelfMask =\t\t Shelf mask file (includes nodata)",
		  "\tfl =\t\t\t Feather length",
		  "\toffsets =\t\t Use offset azimuth offset data for second component ",
		  "\trOffsets =\t\t Use offset data for both components where needed ",
		  "\ttimeThresh =\t\t For 3D crossing offset solution, only use pairs within timeThresh days (def=12)",
		  "\ttimePhaseThresh =\t\t For 3D crossing phase solution, only use pairs within timePhaseThresh days (def=548 days: 1.5 years)",
		  "\tirregFile =\t\t File with list of irregular input data",
		  "\tvzFlag =\t\t Used for changing output in vertical channel 0 for default (vz correction), 1 horizontal 1/sin(psi), 2 for 1/vertical cos(psi), 3 flag for LOS scaled to m/yr, 4 inc angle",
		  "\tstats =\t\t compute unweight mean vx, and vy and standard dev of the inputs (ex,ey) and number of points (vz) - works only for speckleTrack",
		  "\toutputRA =\t\t Output range (vr/er) and azimuth (va/ea) components instead of rotating to map vx/vy",
		  "\tnoSepAscDesc =\t\t For 3d ignore asc/desc and use heading, default use asc/desc",
		  "\tnoTide =\t\t Don't compute value if on shelf",
	  "\tompThreads N =\t\t Set OpenMP thread count (default 4, or OMP_NUM_THREADS if set in shell)",
		  "\tinputFile =\t\t File with input params, dem and geodat filenames",
		  "\tdemFile =\t\t File with nonInsar dem",
		  "\tuseSquint =\t\t Apply per-image squint(r,a) heading correction (phase/crossing-orbit solution only); default off",
		  "\tnoMask =\t\t Ignore any embedded VRT dataset mask band on offset inputs; default off (mask honored when present)",
		  "\tnoPairOverCount =\t\t Disable the (n_A+n_D)/2 crossing-orbit error inflation for pair non-independence; default is enabled",
		  "\trSigmaConst X =\t\t Diagnostic: add a FIXED X metres instead of each frame's own residual; default 0 = off",
		  "\tpairCountLegacy =\t\t Control: use the pre-2026-08 over-count formula (nOuter+nPairs)/2 instead of (nOuter+nPairs/nOuter)/2; default off",
		  "\trhoPhase X =\t\t Pair-correlation parameter for the PHASE crossing over-count factor; default 0.6",
		  "\trhoOffsets X =\t\t Pair-correlation parameter for the crossing-OFFSETS over-count factor; default 0.6",
		  "\tnoRSigmaResidual =\t\t Omit the rparams tie-point fit residual from the range-offset budget, leaving only the local .sr matching sigma; default is to include it",
		  "\tlegacyPairPhase =\t\t Use the ORIGINAL pairwise crossing-PHASE solver; default is the joint (normal-equations) solver",
		  "\tlegacyPairRange =\t\t Use the ORIGINAL pairwise crossing-RANGE solver; default is the joint (normal-equations) solver",
		  "\tjointMaxSigma X =\t\t Reject pixels with sigmaWorst*sqrt(n) > X m/yr (effective per-measurement sigma).  Under the hopper this gates the ONE system holding phase, range and azimuth rows; default 35, 0 = off",
		  "\tjointMaxSigmaPhase/Range X =\t Legacy per-round overrides of -jointMaxSigma (crossing-phase / crossing-range rounds).  Only distinguishable under -legacyCode; defaults 35 / 100",
		  "\tjointErrScale X =\t\t Scale the joint solvers reported sigma by X (caller-supplied 1-sigma calibration); default 1 = off",
		  "\tnoAzimuthRows =\t\t Drop every azimuth row (range-only solve); the control for what the azimuth offsets contribute; default off",
		  "\tspeckleTrackJoint, jointMaxSigmaSpeckle X =\t RETIRED, accepted and ignored (superseded by the hopper)",
		  "\tnoASigmaResidual =\t\t Omit the azparams tie-point fit residual from the azimuth-offset budget; default is to include it",
		  "\ttrue3D =\t\t Solve vx,vy,vz from range+azimuth offsets with NO surface-parallel constraint (offsets only); writes .vz3d/.ez; default off",
		  "\ttrue3DProject =\t\t Reduction test for -true3D: force the surface-parallel projection of the same 3x3 accumulator; not a production path",
		  "\tgateSpeedFrac =\t\t Make the -jointMaxSigma cap speed-aware: effective cap = max(X, F*|v|), so fast ice is not rejected for having proportionally small error.  DEFAULT 0.03; 0 = off",
		  "\tmaxChi2 =\t\t Reject solved pixels whose reduced chi-square exceeds X -- a blunder screen; keep it loose (100-1000).  Default -1 = off",
		  "\tgateAbsolute =\t\t Gate on the formal sigma alone (no sqrt(n) factor), i.e. reject only pixels whose 1-sigma exceeds -jointMaxSigma; default off",
		  "\tnoErrorGate =\t\t Keep pixels with a valid velocity but an invalid formal error; default is to remove them",
		  "\tnoPhaseRows / noRangeRows / noAzimuthRows =\t\t Drop a row type from the hopper solve; default is all three",
		  "\tsigmaAThreshVel =\t\t Drop azimuth rows whose azparams sigma exceeds X m/yr (sigmaAresidual scaled by 365.25/nDays).  DEFAULT 150; <0 disables",
		  "\tnoGateNEff =\t\t Use the RAW row count in the output rejection gate instead of the weighted effective count.  The weighted count is the DEFAULT",
		  "\tlegacyCode =\t\t Use the pre-2026-09 four-round pipeline instead of the hopper (the hopper is now the DEFAULT).  Implied by -stats and -legacyPair*",
		  "\tobsDump =\t\t File of 'lat lon [name]' points; dumps every contributing observation (a_i, w_i, d_i) at those pixels to <outFileBase>.obsDump.  Requires -hopper or -hopper3D (instrumented in those two solvers)",
		  "\thopper =\t\t Put PHASE, RANGE OFFSETS and AZIMUTH OFFSETS in ONE per-frame-weighted solve instead of averaging separate rounds; default off",
		  "\ttrue3DPhase =\t\t Solve vx,vy,vz from PHASE with no surface-parallel constraint; different ionosphere path from -true3D; default off",
		  "\ttrue3DMaxSigma X =\t\t Reject true3D pixels with sigmaWorst*sqrt(n) > X m/yr; default 100, 0 = off",
		  "\toutputImage =\t\t Root of output image (e.g., mosaicOffsets)");
}

/*
  Estimate output area for image based on location of inputs
*/
static void findOutBounds(outputImageStructure *outputImage, inputImageStructure *ascImages, inputImageStructure *descImages, landSatImage *LSImages, int32_t *autoSize, int32_t writeBlank)
{
	double minX, maxX, minY, maxY, x, y;
	double minXLS, maxXLS, minYLS, maxYLS;
	inputImageStructure *tmp; /* Tmp list for looping */
	landSatImage *LStmp;
	int32_t j;
	/*
	  No data case for tiepoints: no input frames (e.g. a bedrock/-extraTies tie
	  that only echoes external tie points). There is no meaningful output grid.
	  Keep xSize/ySize at 0 -- the mosaic pixel loops then run zero iterations,
	  contributing nothing (writeTieFile still echoes the -extraTies) -- but do
	  NOT fall through to the autosize block below: autosizing from zero images
	  leaves min/max at their 1e30 sentinels, and (int32_t)1e30 is out of range
	  and overflows to a negative xSize/ySize. Skip it explicitly here.
	*/
	if (ascImages == NULL && descImages == NULL && LSImages == NULL && writeBlank == FALSE)
	{
		outputImage->xSize = (int)0;
		outputImage->ySize = (int)0;
		outputImage->originX = 0;
		outputImage->originY = 0;
		outputImage->deltaX = 1.;
		outputImage->deltaY = 1.;
		/* autoSize=TRUE tells writeTieFile() to output ALL tie points (huge
		   +/-1e9 bounds) instead of clipping to this degenerate output grid --
		   the extraTies span the whole region and must not be filtered out.
		   The original code got this by falling through to the autosize block
		   below (which sets *autoSize=TRUE); we skip that block to avoid its
		   integer overflow, so set the same flag explicitly here. */
		*autoSize = TRUE;
	}
	/* If ouput size zero, auto size */
	else if (outputImage->xSize == 0 || outputImage->ySize == 0)
	{
		minXLS = (double)LARGEINT;
		maxXLS = -(double)LARGEINT;
		minYLS = (double)LARGEINT;
		maxYLS = -(double)LARGEINT;
		if (LSImages != NULL)
		{
			for (LStmp = LSImages; LStmp != NULL; LStmp = LStmp->next)
			{
				minXLS = min(minXLS, LStmp->matches.x0 * MTOKM);
				maxXLS = max(maxXLS, (LStmp->matches.x0 + LStmp->matches.stepX * LStmp->matches.dx * (LStmp->matches.nx - 1)) * MTOKM);
				minYLS = min(minYLS, LStmp->matches.y0 * MTOKM);
				maxYLS = max(maxYLS, (LStmp->matches.y0 + LStmp->matches.stepY * LStmp->matches.dy * (LStmp->matches.ny - 1)) * MTOKM);
			}
			fprintf(stderr, "LS Bounds  %lf %lf %lf %lf\n", minXLS, maxXLS, minYLS, maxYLS);
		}
		/*
		  Find output bounds
		*/
		fprintf(outputImage->fpLog, ";\n; Entering findOutBounds (mosaic3d.c)\n;\n");
		minX = 1.e30;
		minY = 1.0e30;
		maxX = -1.e30;
		maxY = -1.0e30;
		for (tmp = ascImages; tmp != NULL; tmp = tmp->next)
		{
			fprintf(stderr, "A %f %f %f %f\n", tmp->minX, tmp->minY, tmp->maxX, tmp->maxY);
			minX = min(tmp->minX, minX);
			minY = min(tmp->minY, minY);
			maxX = max(tmp->maxX, maxX);
			maxY = max(tmp->maxY, maxY);
			/*
			for(j=1; j < 5; j++) {
				lltoxy1(tmp->latControlPoints[j],tmp->lonControlPoints[j],&x,&y,
					Rotation,outputImage->slat);
				minX=min(x,minX);minY=min(y,minY);maxX=max(x,maxX); maxY=max(y,maxY);
			}*/
		}
		for (tmp = descImages; tmp != NULL; tmp = tmp->next)
		{
			fprintf(stderr, "D %f %f %f %f\n", tmp->minX, tmp->minY, tmp->maxX, tmp->maxY);
			minX = min(tmp->minX, minX);
			minY = min(tmp->minY, minY);
			maxX = max(tmp->maxX, maxX);
			maxY = max(tmp->maxY, maxY);
			/*
			for(j=1; j < 5; j++) {
				lltoxy1(tmp->latControlPoints[j],tmp->lonControlPoints[j],&x,&y,
					Rotation,outputImage->slat);
				minX=min(x,minX);minY=min(y,minY);maxX=max(x,maxX); maxY=max(y,maxY);
				fprintf(stderr,"%f %f\n", x,y);
				fprintf(stderr,"%f %f %f %f ..\n", tmp->minX, tmp->maxX, tmp->minY, tmp->maxY);
			} */
		}
		if (LSImages != NULL)
		{
			minX = min(minX, minXLS);
			maxX = max(maxX, maxXLS);
			minY = min(minY, minYLS);
			maxY = max(maxY, maxYLS);
		}
		minX = (double)((int32_t)minX - 3);
		maxX = (double)((int32_t)maxX + 3);
		minY = (double)((int32_t)minY - 3);
		maxY = (double)((int32_t)maxY + 3);
		fprintf(stderr, "minX,maxX,minY,maxY %f %f %f %f %i %i", minX, maxX, minY, maxY, outputImage->xSize, outputImage->ySize);
		*autoSize = TRUE;
		fprintf(stderr, "\n\n*** AUTOSIZING REGION *** *\n\n");
		fprintf(outputImage->fpLog, ";   *** AUTOSIZING REGION ***\n;\n");
		fprintf(stderr, "minX,maxX,minY,maxY %f %f %f %f %i %i %f\n", minX, maxX, minY, maxY, outputImage->xSize, outputImage->ySize, outputImage->deltaX);
		outputImage->xSize = (int)((maxX - minX) / (outputImage->deltaX * MTOKM) + 0.5);
		outputImage->ySize = (int)((maxY - minY) / (outputImage->deltaY * MTOKM) + 0.5);
		fprintf(stderr, "minX,maxX,minY,maxY %f %f %f %f %i %i\n", minX, maxX, minY, maxY, outputImage->xSize, outputImage->ySize);
		outputImage->originX = minX * KMTOM;
		outputImage->originY = minY * KMTOM;
	}
	else
	{
		fprintf(stderr, "\n;   *** USER SPECIFIED SIZE ***\n;\n");
		*autoSize = FALSE;
	}

	fprintf(stderr, "dimensions %i %i %f %f %f %f\n\n", outputImage->xSize, outputImage->ySize,
			outputImage->originX, outputImage->originY, outputImage->deltaX, outputImage->deltaY);
	fprintf(outputImage->fpLog, "; Size       : %i %i\n; Origin     : %f %f \n; Spacing    : %f %f\n;\n",
			outputImage->xSize, outputImage->ySize, outputImage->originX, outputImage->originY, outputImage->deltaX, outputImage->deltaY);
	fprintf(outputImage->fpLog, ";\n; Leaving findOutBounds (mosaic3d.c)\n;\n");
	fflush(outputImage->fpLog);
}

static void mallocOutputImage(outputImageStructure *outputImage)
{
	float *buf1, *buf2, *buf3, *buf1s, *buf2s, *buf3s;
	float *bufex, *bufey;
	float **vXimage, **vYimage, **vZimage;
	float **scaleX, **scaleY, **scaleZ;
	int32_t bufSize, i, j, k;
	float *bufx, *bufy, *bufz, *bufs, *bufsx, *bufsy;
	fprintf(outputImage->fpLog, ";\n; Entering  mallocOutputImage (from mosaic3d)\n");
	outputImage->imageType = POWER;
	bufSize = outputImage->ySize * sizeof(float *);
	fprintf(stderr, "%i\n", bufSize);

	if (outputImage->singleImageFastPath == TRUE)
	{
		/* One contributing image, nothing else touching the grid: the weighted
		   accumulation buffers (image/image2/image3/scale/scale2/scale3) and the
		   per-image weighting scratch (sxTmp/syTmp/fScale) are provably unnecessary
		   -- see the singleImageFastPath comment in geocode.h and speckleTrackMosaic.c.
		   Leave them NULL entirely instead of allocating ~700MB-scale buffers per
		   grid dimension that would never be read. */
		outputImage->image = NULL;
		outputImage->image2 = NULL;
		outputImage->image3 = NULL;
		outputImage->scale = NULL;
		outputImage->scale2 = NULL;
		outputImage->scale3 = NULL;
		outputImage->sxTmp = NULL;
		outputImage->syTmp = NULL;
		outputImage->fScale = NULL;
		outputImage->vxTmp = (float **)malloc(bufSize);
		outputImage->vyTmp = (float **)malloc(bufSize);
		outputImage->vzTmp = (float **)malloc(bufSize);
		outputImage->errorX = (float **)malloc(bufSize);
		outputImage->errorY = (float **)malloc(bufSize);

		bufSize = outputImage->xSize * outputImage->ySize * sizeof(float);
		bufx = (float *)malloc(bufSize);
		bufy = (float *)malloc(bufSize);
		bufz = (float *)malloc(bufSize);
		bufex = (float *)malloc(bufSize);
		bufey = (float *)malloc(bufSize);

		for (i = 0; i < outputImage->ySize; i++)
		{
			outputImage->vxTmp[i] = (float *)&(bufx[i * outputImage->xSize]);
			outputImage->vyTmp[i] = (float *)&(bufy[i * outputImage->xSize]);
			outputImage->vzTmp[i] = (float *)&(bufz[i * outputImage->xSize]);
			outputImage->errorX[i] = (float *)&(bufex[i * outputImage->xSize]);
			outputImage->errorY[i] = (float *)&(bufey[i * outputImage->xSize]);
		}

		/* Sentinel-fill up front: speckleTrackMosaic()'s fast path only visits rows
		   within its per-row-tightened footprint bound, so pixels outside it (and
		   never-touched rows entirely) must already read back as "no data" once
		   image/image2/image3 are aliased onto these buffers. */
		for (j = 0; j < outputImage->ySize; j++)
		{
			for (k = 0; k < outputImage->xSize; k++)
			{
				outputImage->vxTmp[j][k] = -LARGEINT;
				outputImage->vyTmp[j][k] = -LARGEINT;
				outputImage->vzTmp[j][k] = -LARGEINT;
				outputImage->errorX[j][k] = -LARGEINT;
				outputImage->errorY[j][k] = -LARGEINT;
			}
		}
		fprintf(outputImage->fpLog, "; Returned from  mallocOutputImage (singleImageFastPath)\n");
		fflush(outputImage->fpLog);
		return;
	}

	outputImage->image = (void **)malloc(bufSize);
	outputImage->image2 = (void **)malloc(bufSize);
	outputImage->image3 = (void **)malloc(bufSize);
	outputImage->scale = (float **)malloc(bufSize);
	outputImage->scale2 = (float **)malloc(bufSize);
	outputImage->scale3 = (float **)malloc(bufSize);
	outputImage->vxTmp = (float **)malloc(bufSize);
	outputImage->sxTmp = (float **)malloc(bufSize);
	outputImage->vyTmp = (float **)malloc(bufSize);
	outputImage->syTmp = (float **)malloc(bufSize);
	outputImage->vzTmp = (float **)malloc(bufSize);
	outputImage->fScale = (float **)malloc(bufSize);
	outputImage->errorX = (float **)malloc(bufSize);
	outputImage->errorY = (float **)malloc(bufSize);

	bufSize = outputImage->xSize * outputImage->ySize * sizeof(float);
	buf1 = (float *)malloc(bufSize);
	buf2 = (float *)malloc(bufSize);
	buf3 = (float *)malloc(bufSize);
	buf1s = (float *)malloc(bufSize);
	buf2s = (float *)malloc(bufSize);
	buf3s = (float *)malloc(bufSize);
	bufx = (float *)malloc(bufSize);
	bufy = (float *)malloc(bufSize);
	bufz = (float *)malloc(bufSize);
	bufs = (float *)malloc(bufSize);
	bufsx = (float *)malloc(bufSize);
	bufsy = (float *)malloc(bufSize);
	bufex = (float *)malloc(bufSize);
	bufey = (float *)malloc(bufSize);

	for (i = 0; i < outputImage->ySize; i++)
	{
		outputImage->image[i] = (void *)&(buf1[i * outputImage->xSize]);
		outputImage->image2[i] = (void *)&(buf2[i * outputImage->xSize]);
		outputImage->image3[i] = (void *)&(buf3[i * outputImage->xSize]);
		outputImage->scale[i] = (float *)&(buf1s[i * outputImage->xSize]);
		outputImage->scale2[i] = (float *)&(buf2s[i * outputImage->xSize]);
		outputImage->scale3[i] = (float *)&(buf3s[i * outputImage->xSize]);
		outputImage->vxTmp[i] = (float *)&(bufx[i * outputImage->xSize]);
		outputImage->vyTmp[i] = (float *)&(bufy[i * outputImage->xSize]);
		outputImage->vzTmp[i] = (float *)&(bufz[i * outputImage->xSize]);
		outputImage->sxTmp[i] = (float *)&(bufsx[i * outputImage->xSize]);
		outputImage->syTmp[i] = (float *)&(bufsy[i * outputImage->xSize]);
		outputImage->fScale[i] = (float *)&(bufs[i * outputImage->xSize]);
		outputImage->errorX[i] = (float *)&(bufex[i * outputImage->xSize]);
		outputImage->errorY[i] = (float *)&(bufey[i * outputImage->xSize]);
	}

	vXimage = (float **)outputImage->image;
	vYimage = (float **)outputImage->image2;
	vZimage = (float **)outputImage->image3;
	scaleX = (float **)outputImage->scale;
	scaleY = (float **)outputImage->scale2;
	scaleZ = (float **)outputImage->scale3;

	for (j = 0; j < outputImage->ySize; j++)
	{
		for (k = 0; k < outputImage->xSize; k++)
		{
			scaleX[j][k] = 0.0;
			vXimage[j][k] = -LARGEINT;
			scaleY[j][k] = 0.0;
			vYimage[j][k] = -LARGEINT;
			scaleZ[j][k] = 0.0;
			vZimage[j][k] = -LARGEINT;
		}
	}
	fprintf(outputImage->fpLog, "; Returned from  mallocOutputImage\n");
	fflush(outputImage->fpLog);
}

static void readReferenceVelMosaic(referenceVelocity *refVel, outputImageStructure *outputImage)
{
	FILE *fp;
	char *geodatFile, *vxFile, *vyFile, *exFile, *eyFile;
	double maxX, maxY;
	char line[1500];
	uint32_t lineCount = 0;
	int32_t eod;
	float *tmp, *tmp1;
	float dum1, dum2;
	uint32_t nx, ny;
	double x0, y0;
	uint32_t xoff, yoff, tail; /* offset into velocity file in samples */
	uint32_t i;
	/*
	  geodat file name
	*/
	fprintf(stderr, "**** Start reading reference velocity file ****\n");
	geodatFile = (char *)malloc(strlen(refVel->velFile) + 11);
	geodatFile[0] = '\0';
	geodatFile = strcpy(geodatFile, refVel->velFile);
	geodatFile = strcat(geodatFile, ".vx.geodat");
	fprintf(stderr, "vel geodat file %s\n", geodatFile);
	/*
	   velocity file names
	*/
	vxFile = (char *)malloc(strlen(refVel->velFile) + 4);
	vxFile[0] = '\0';
	vxFile = strcpy(vxFile, refVel->velFile);
	vxFile = strcat(vxFile, ".vx");
	vyFile = (char *)malloc(strlen(refVel->velFile) + 4);
	vyFile[0] = '\0';
	vyFile = strcpy(vyFile, refVel->velFile);
	vyFile = strcat(vyFile, ".vy");
	fprintf(stderr, "vx,vy file %s %s\n", vxFile, vyFile);
	/*
	  Open geodat file
	*/
	fp = openInputFile(geodatFile);
	if (fp == NULL)
		error("*** readREfVelFile: Error opening %s ***\n", geodatFile);
	/*
	  Read parameters
	*/
	lineCount = getDataString(fp, lineCount, line, &eod); /* Skip # 2 line */
	lineCount = getDataString(fp, lineCount, line, &eod);
	sscanf(line, "%f %f\n", &dum1, &dum2); /* read as float in case fp value */
	nx = (int)dum1;
	ny = (int)dum2;

	lineCount = getDataString(fp, lineCount, line, &eod);
	sscanf(line, "%lf %lf\n", &(refVel->dx), &(refVel->dy));

	lineCount = getDataString(fp, lineCount, line, &eod);
	sscanf(line, "%lf %lf\n", &(x0), &(y0));
	x0 *= KMTOM;
	y0 *= KMTOM;
	fclose(fp);
	fprintf(stderr, "%i %i \n %f %f \n %f %f \n", nx, ny, refVel->dx, refVel->dy, x0, y0);
	refVel->x0 = max(x0, outputImage->originX);
	refVel->y0 = max(y0, outputImage->originY);
	fprintf(stderr, "output grid %lf %lf \n", outputImage->originX, outputImage->originY);
	maxX = min(x0 + refVel->dx * (nx - 1), outputImage->originX + outputImage->deltaX * (outputImage->xSize - 1));
	maxY = min(y0 + refVel->dy * (ny - 1), outputImage->originY + outputImage->deltaY * (outputImage->ySize - 1));
	refVel->nx = (maxX - refVel->x0) / (refVel->dx) + 1;
	refVel->ny = (maxY - refVel->y0) / (refVel->dy) + 1;
	fprintf(stderr, "%f %f  %i %i \n", refVel->x0, refVel->y0, refVel->nx, refVel->ny);
	/* compute file offsets */
	xoff = (uint32_t)((refVel->x0 - x0) / refVel->dx + 0.5);
	yoff = (uint32_t)((refVel->y0 - y0) / refVel->dy + 0.5);
	fprintf(stderr, "xoff,yoff %i %i\n", xoff, yoff);
	/*
	  Malloc array
	*/
	refVel->vx = (float **)malloc(refVel->ny * sizeof(float *));
	tmp = (float *)malloc(refVel->nx * refVel->ny * sizeof(float));
	/*
	  Open vx file
	*/
	fp = openInputFile(vxFile);
	fseek(fp, (nx * yoff * sizeof(float)), SEEK_SET);
	tail = max(nx - (xoff + refVel->nx), 0);
	for (i = 0; i < refVel->ny; i++)
	{
		tmp1 = &(tmp[i * refVel->nx]);
		fseek(fp, (xoff * sizeof(float)), SEEK_CUR);
		freadBS(tmp1, sizeof(float), refVel->nx, fp, FLOAT32FLAG);
		fseek(fp, (tail * sizeof(float)), SEEK_CUR);
		refVel->vx[i] = tmp1;
	}
	fclose(fp);
	/*
	  Malloc array
	*/
	refVel->vy = (float **)malloc(refVel->ny * sizeof(float *));
	tmp = (float *)malloc(refVel->nx * refVel->ny * sizeof(float));
	/*
	  Open vy file
	*/
	fp = openInputFile(vyFile);
	fseek(fp, (nx * yoff * sizeof(float)), SEEK_SET);
	for (i = 0; i < refVel->ny; i++)
	{
		tmp1 = &(tmp[i * refVel->nx]);
		fseek(fp, (xoff * sizeof(float)), SEEK_CUR);
		freadBS(tmp1, sizeof(float), refVel->nx, fp, FLOAT32FLAG);
		fseek(fp, (tail * sizeof(float)), SEEK_CUR);
		refVel->vy[i] = tmp1;
	}
	fclose(fp);

	/* Error files  if needed */
	if (refVel->initMapFlag == TRUE)
	{
		/*
		   velocity file names
		*/
		exFile = (char *)malloc(strlen(refVel->velFile) + 4);
		exFile[0] = '\0';
		exFile = strcpy(exFile, refVel->velFile);
		exFile = strcat(exFile, ".ex");
		eyFile = (char *)malloc(strlen(refVel->velFile) + 4);
		eyFile[0] = '\0';
		eyFile = strcpy(eyFile, refVel->velFile);
		eyFile = strcat(eyFile, ".ey");
		/*
		  Malloc array
		*/
		refVel->ex = (float **)malloc(refVel->ny * sizeof(float *));
		tmp = (float *)malloc(refVel->nx * refVel->ny * sizeof(float));
		/*
		  Open ex file
		*/
		fp = openInputFile(exFile);
		fseek(fp, (nx * yoff * sizeof(float)), SEEK_SET);
		tail = max(nx - (xoff + refVel->nx), 0);
		for (i = 0; i < refVel->ny; i++)
		{
			tmp1 = &(tmp[i * refVel->nx]);
			fseek(fp, (xoff * sizeof(float)), SEEK_CUR);
			freadBS(tmp1, sizeof(float), refVel->nx, fp, FLOAT32FLAG);
			fseek(fp, (tail * sizeof(float)), SEEK_CUR);
			refVel->ex[i] = tmp1;
		}
		fclose(fp);
		/*
		  Malloc array
		*/
		refVel->ey = (float **)malloc(refVel->ny * sizeof(float *));
		tmp = (float *)malloc(refVel->nx * refVel->ny * sizeof(float));
		/*
		  Open ey file
		*/
		fp = openInputFile(eyFile);
		fseek(fp, (nx * yoff * sizeof(float)), SEEK_SET);
		tail = max(nx - (xoff + refVel->nx), 0);
		for (i = 0; i < refVel->ny; i++)
		{
			tmp1 = &(tmp[i * refVel->nx]);
			fseek(fp, (xoff * sizeof(float)), SEEK_CUR);
			freadBS(tmp1, sizeof(float), refVel->nx, fp, FLOAT32FLAG);
			fseek(fp, (tail * sizeof(float)), SEEK_CUR);
			refVel->ey[i] = tmp1;
		}
		fclose(fp);
		fprintf(stderr, "**** reading error map ****\n");
	}
	fprintf(stderr, "**** End reading reference velocity file ****\n");
}

#define IGREG 2299161

void caldat(int32_t julian, int32_t *mm, int32_t *id, int32_t *iyyy)
{
	int32_t ja, jalpha, jb, jc, jd, je;

	if (julian >= IGREG)
	{
		jalpha = (int32_t)(((float)(julian - 1867216) - 0.25) / 36524.25);
		ja = julian + 1 + jalpha - (int32_t)(0.25 * jalpha);
	}
	else
		ja = julian;
	jb = ja + 1524;
	jc = (int32_t)(6680.0 + ((float)(jb - 2439870) - 122.1) / 365.25);
	jd = (int32_t)(365 * jc + (0.25 * jc));
	je = (int32_t)((jb - jd) / 30.6001);
	*id = jb - jd - (int32_t)(30.6001 * je);
	*mm = je - 1;
	if (*mm > 12)
		*mm -= 12;
	*iyyy = jc - 4715;
	if (*mm > 2)
		--(*iyyy);
	if (*iyyy <= 0)
		--(*iyyy);
}
#undef IGREG
