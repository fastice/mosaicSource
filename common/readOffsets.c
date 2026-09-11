#include "stdio.h"
#include "string.h"
#include <stdlib.h>
#include "common.h"
#include <libgen.h>
#include <unistd.h>
#include <time.h>
#include "gdalIO/gdalIO/grimpgdal.h"

extern int32_t indentRegionOutput;
extern int32_t noMask; /* ignore any embedded VRT dataset mask band on offset inputs */
extern int32_t skipAzimuthOffsets; /* skip the azimuth raster read; see getRegion.c */

/* Running grand totals of I/O time across the whole run, accumulated by
   readRangeOrRangeOffsets() and getMosaicInputImage() respectively; printed
   once by make3DOffsets.c/make3DMosaic.c at the end of their main loops. */
double totalOffsetsIOTime = 0.0;
double totalPhaseIOTime = 0.0;

static double now(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}


#define AZONLY 30
#define RGANDAZ 31
#define RGONLY 32
#define AZFORRANGE 33

static void readCov(FILE *fp, int32_t n, double C[7][7], double *sigmaResidual, char *line);
static void readOffsetFile(float **data, int32_t nr, int32_t na, char *offsetFile);
static void mapBuffer(int32_t nr, int32_t na, float **d, float **s,  float *buffSpaceD, float *buffSpaceS);
static void initOffsetBuffers(Offsets *offsets, int32_t mode);
static char *mergePath(char *file1, char *path);
static char *RgOffsetsParamName(char *rParamsFile, char *newFile, int32_t deltaB, char *verticalCorrectionSuffix);
static char *AzOffsetsParamName(char *aParamsFile, char *newFile, int32_t deltaB);
static int is_not_final_range(const char *vrtFile);
static int get_se_correction_path(const char *vrtFile, char *outPath, size_t outSize);

static char *checkForOffsetsVrt(char *filename, char *vrtBuff) 
{
	char noExtensionName[2048], *vrtFile; 
	int32_t noExtensionNameLength; 
	vrtFile = checkForVrt(filename, vrtBuff);
	// This will return the vrt if found (e.g., range.offsets.vrt)
	if(vrtFile != NULL) return vrtFile;
	// Else see if for xxx.yyy.da (.dr, .sa, .sx) there is an xxx.yyy.vrt
	noExtensionName[0] = '\0';
	noExtensionNameLength = strlen(filename)-3;
	strncpy(noExtensionName, filename, noExtensionNameLength);
	noExtensionName[noExtensionNameLength] = '\0';
	return checkForVrt(noExtensionName, vrtBuff);
}

static GDALRasterBandH getBandAndMeta(GDALDatasetH hDS, Offsets *offsets, int32_t band, char *path);
/*
   Read the offset data and paramter files
 */
void readOffsetDataAndParams(Offsets *offsets, float azimuthMin, float azimuthMax)
{	
	getRParams(offsets);
	fprintf(stderr, "Reading ionosphere correction file: {%s}\n", offsets->rOffCorrection.correctionFile);
	
	readBothOffsets(offsets, azimuthMin, azimuthMax);
	getAzParams(offsets);
	fprintf(stderr, "Offsets and parameters read\n");
	if (offsets->deltaB != DELTABNONE && offsets->geo2 == NULL)
		error("offsets deltaB set but no second geodat for %s\n", offsets->rFile);
	//fprintf(stderr,"\nDEBUG: geo1/2 %s %s\n", offsets->geo1, offsets->geo2);
}

static void mapBuffer(int32_t nr, int32_t na, float **d, float **s, float *buffSpaceD, float *buffSpaceS)
{
	int i;
	for (i = 0; i < na; i++)
	{
		d[i] = &(buffSpaceD[i * nr]);
		s[i] = &(buffSpaceS[i * nr]);
	}
}

static void initOffsetBuffers(Offsets *offsets, int32_t mode)
{
	extern void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4, *offSEBuffSpace;
	extern void *lBuf1, *lBuf2, *lBuf3, *lBuf4, *lSEBuf;
	float *fBuf1, *fBuf2, *fBuf3, *fBuf4;
	int32_t i, nr, na;
	nr = offsets->nr;
	na = offsets->na;
	if (nr * na * 4 > MAXOFFBUF)
		error("Offsets Buffer Size exceeds MAXOFFBUF of %i; nr %i na %i mode %i\n %s\n", MAXOFFBUF, nr, na, offsets->file);

	switch(mode) {
		// If RANGANDAZ Fall through to both
		case RGONLY: 
			offsets->dr = (float **)lBuf3;
			offsets->sr = (float **)lBuf4;
			fBuf3 = (float *)offBufSpace3;
			fBuf4 = (float *)offBufSpace4;
			// fprintf(stderr, "RANGe BUFFERS\n");
			mapBuffer(nr, na, offsets->dr, offsets->sr, fBuf3, fBuf4);
		 	break;
		case AZONLY:
			offsets->da = (float **)lBuf1;
			offsets->sa = (float **)lBuf2;
			fBuf1 = (float *)offBufSpace1;
			fBuf2 = (float *)offBufSpace2;
			// fprintf(stderr, "AZIMUTH BUFFERS\n");
			mapBuffer(nr, na, offsets->da, offsets->sa, fBuf1, fBuf2);
		    break;
		case AZFORRANGE:
			offsets->dr = (float **)lBuf1;
			offsets->sr = (float **)lBuf2;
			fBuf3 = (float *)offBufSpace1;
			fBuf4 = (float *)offBufSpace2;
			// fprintf(stderr, "RANGE FLIP BUFFERS\n");
			mapBuffer(nr, na, offsets->dr, offsets->sr, fBuf3, fBuf4);
			break;
		case SEBUF:
			offsets->SECorrection = (float **)lSEBuf;
			fBuf1 = (float *)offSEBuffSpace;
			//fprintf(stderr, "Special Solid Earth BUFFERS\n");
			mapBuffer(nr, na, offsets->SECorrection, offsets->SECorrection, fBuf1, fBuf1);
			break;		
		default:
		error("initOffsets invalid buffer code %i", mode);
			break;
	}
}

static char *mergePath(char *file1, char *path)
{
	char fileTmp[1024], *tmp, *merged;
	fileTmp[0] = '\0';
	if (path != NULL)
	{
		tmp = strcat(fileTmp, path);
		tmp = strcat(fileTmp, "/");
	}
	tmp = strcat(fileTmp, file1);
	merged = (char *)calloc(strlen(fileTmp) + 1, sizeof(char));
	merged[0] = '\0';
	return (strcat(merged, fileTmp));
}
/*
  Read the params, use merge=TRUE to update geodat paths, FALSE to leave unchanged.
*/
void readOffsetParams(char *datFile, Offsets *offsets, int32_t merge)
{
	FILE *fp;
	int32_t rO, aO, nr, na;
	float deltaA, deltaR;
	double sigmaS, sigmaR;
	char line[1024], *tmp;
	char file1[512], file2[512], *path, buf[2048];
	int32_t lineCount = 0, eod;
	int32_t nRead;
	/* See if vrt exits */
	buf[0] = '\0';
	/* Read params file */
	fp = openInputFile(datFile);
	lineCount = getDataString(fp, lineCount, line, &eod);
	nRead = sscanf(line, "%i%i%i%i%f%f%lf%lf", &rO, &aO, &nr, &na, &deltaR, &deltaA, &sigmaS, &sigmaR);

	if (nRead != 6 && nRead != 7 && nRead != 8)
		error("%s  %i of %s", "readOffsets -- Missing image parameters at line:", lineCount, datFile);
	else if (nRead == 6)
	{
		//fprintf(stderr, "**** WARNING-MISSING SIGMA STREAKS for %s\n", datFile);
		sigmaS = 0.0;
		sigmaR = 0.0;
	}
	else if (nRead == 7)
		sigmaR = 0.0;
	/* load param in structure */
	offsets->rO = rO;
	offsets->aO = aO;
	offsets->deltaA = deltaA;
	offsets->deltaR = deltaR;
	offsets->sigmaStreaks = sigmaS;
	offsets->sigmaRange = sigmaR;
	fprintf(stderr, "SigmaStreaks/Range = %f %f %i\n", sigmaS, sigmaR, nRead);
	offsets->nr = nr;
	offsets->na = na;
	// fprintf(stderr, "nr na %i %i\n", nr, na);
	/*
	   read geodat files if the exist
	 */
	tmp = fgets(line, 1024, fp);
	offsets->geo1 = NULL;
	offsets->geo2 = NULL;
	if (tmp != NULL)
	{
		sscanf(line, "%s %s %s", file1, file2, datFile);
		if (merge == TRUE)
			path = dirname(strcpy(buf, datFile));
		else
			path = NULL;
		offsets->geo1 = mergePath(file1, path);
		offsets->geo2 = mergePath(file2, path);
		//fprintf(stderr, "Found geo1 & geo2. %s %s \n", offsets->geo1, offsets->geo2);
	}
	fclose(fp);
}

/*
	read a single offset file
*/
static void readOffsetFile(float **data, int32_t nr, int32_t na, char *offsetFile)
{
	FILE *fp;
	int32_t i;
	/*fprintf(stderr,"--- Reading offset file %s",offsetFile);*/
	fp = openInputFile(offsetFile);
	for (i = 0; i < na; i++)
		freadBS(data[i], sizeof(float), nr, fp, FLOAT32FLAG);
	fclose(fp);
	/* fprintf(stderr,"- done - \n");	*/
}

static void readRParamsYaml(FILE *fp, Offsets *offsets)
{
	char line[512];
	int inC = 0, ci = 0;
	/* zero all output fields */
	offsets->bn = 0.0; offsets->bp = 0.0;
	offsets->dBn = 0.0; offsets->dBp = 0.0;
	offsets->rConst = 0.0;
	offsets->dBnQ = 0.0; offsets->dBpQ = 0.0;
	offsets->sigmaRresidual = 0.0;
	offsets->rOffCorrection.correctionFile[0] = '\0';
	for (int i = 1; i <= 6; i++)
		for (int j = 1; j <= 6; j++)
			offsets->Cr[i][j] = 0.0;

	while (fgets(line, sizeof(line), fp)) {
		if      (sscanf(line, "sigma: %lf",  &offsets->sigmaRresidual) == 1) { inC = 0; }
		else if (sscanf(line, "Bn: %lf",     &offsets->bn)             == 1) { inC = 0; }
		else if (sscanf(line, "Bp: %lf",     &offsets->bp)             == 1) { inC = 0; }
		else if (sscanf(line, "dBn: %lf",    &offsets->dBn)            == 1) { inC = 0; }
		else if (sscanf(line, "dBp: %lf",    &offsets->dBp)            == 1) { inC = 0; }
		else if (strncmp(line, "cnst:", 5) == 0) {
			double cnstTmp = 0.0;
			sscanf(line + 5, "%lf", &cnstTmp);
			if (offsets->deltaB == DELTABNONE)
				offsets->rConst = cnstTmp;
			inC = 0;
		}
		else if (sscanf(line, "dBnQ: %lf",   &offsets->dBnQ)           == 1) { inC = 0; }
		else if (sscanf(line, "dBpQ: %lf",   &offsets->dBpQ)           == 1) { inC = 0; }
		else if (strncmp(line, "C:", 2) == 0) { inC = 1; ci = 0; }
		else if (inC && strstr(line, "- [") && ci < 6) {
			char *p = strstr(line, "[");
			if (p)
				sscanf(p + 1, "%lf, %lf, %lf, %lf, %lf, %lf",
				       &offsets->Cr[ci+1][1], &offsets->Cr[ci+1][2],
				       &offsets->Cr[ci+1][3], &offsets->Cr[ci+1][4],
				       &offsets->Cr[ci+1][5], &offsets->Cr[ci+1][6]);
			ci++;
		}
		else if (strstr(line, "offsetCorrectionFile:") != NULL) {
			char *p = strstr(line, "offsetCorrectionFile:") + strlen("offsetCorrectionFile:");
			while (*p == ' ' || *p == '\t') p++;
			if (strncmp(p, "nil", 3) != 0) {
				strncpy(offsets->rOffCorrection.correctionFile, p,
				        sizeof(offsets->rOffCorrection.correctionFile) - 1);
				offsets->rOffCorrection.correctionFile[sizeof(offsets->rOffCorrection.correctionFile) - 1] = '\0';
				int32_t len = strlen(offsets->rOffCorrection.correctionFile);
				while (len > 0 && (offsets->rOffCorrection.correctionFile[len-1] == '\n' ||
				                   offsets->rOffCorrection.correctionFile[len-1] == '\r' ||
				                   offsets->rOffCorrection.correctionFile[len-1] == ' '))
					offsets->rOffCorrection.correctionFile[--len] = '\0';
			}
			inC = 0;
		}
		else { inC = 0; }
	}
	fprintf(stderr, "%srange sigma*sqrt(X2/n) = %lf (m)\n", indentRegionOutput ? "\t" : "", offsets->sigmaRresidual);
}

static char *RgOffsetsParamName(char *rParamsFile, char *newFile, int32_t deltaB, char *verticalCorrectionSuffix)
{
	char *suffix[3] = {"", ".deltabp", ".quad"};
	char *yamlSuffix[3] = {".yaml", ".deltabp.yaml", ".deltabquad.yaml"};
	char tmpSuffix[64], fullSuffix[64];
	size_t blen;

	if (deltaB > DELTABQUAD || deltaB < DELTABNONE)
		error("invalide deltaB flag %i", deltaB);

	/* 2026-06-24: rParamsFile is written once by the Python orchestration layer
	   with a fixed name (typically the .deltabp.yaml form), independent of which
	   deltaB mode mosaic3d ends up actually requesting at runtime -- so the mode
	   "embedded" in the given name can't be trusted. Strip any existing mode infix
	   (.deltabp/.deltabquad) plus the .yaml extension to recover the canonical base
	   (e.g. "rBaseline"), then re-derive the filename for the CURRENT deltaB. This
	   was the bug behind mosaic3d silently reading the DELTABCONST fit even when
	   run without -SVConst (deltaB=NONE), since the old code just returned
	   whatever name it was given unchanged. getRParams()'s caller falls back to
	   the as-given rParamsFile if this derived name doesn't exist on disk, so
	   older frames missing a mode-specific variant still work as before. */
	blen = strlen(rParamsFile);
	if (blen > 5 && strcmp(rParamsFile + blen - 5, ".yaml") == 0) {
		char *yamlInfix[3] = {"", ".deltabp", ".deltabquad"};
		char base[2048];
		size_t baseLen = blen - 5; /* strip ".yaml" */
		int k;
		strncpy(base, rParamsFile, baseLen);
		base[baseLen] = '\0';
		for (k = 1; k <= 2; k++) {
			size_t ilen = strlen(yamlInfix[k]);
			if (baseLen > ilen && strcmp(base + baseLen - ilen, yamlInfix[k]) == 0) {
				base[baseLen - ilen] = '\0';
				break;
			}
		}
		strcpy(newFile, base);
		strcat(newFile, yamlInfix[deltaB]);
		if (verticalCorrectionSuffix != NULL && deltaB != DELTABNONE) {
			/* e.g. rBaseline.deltabp.yaml -> rBaseline.deltabp.480.yaml */
			strcat(newFile, ".");
			strncat(newFile, verticalCorrectionSuffix, 32);
		}
		strcat(newFile, ".yaml");
		return newFile;
	}

	/* If rParamsFile doesn't end in .yaml, check whether we should build a yaml name.
	   We do this when the deltaB-specific yaml file would exist but the old file might not.
	   The caller (getRParams) handles the actual open/fallback logic; here we just construct
	   the canonical name for the requested deltaB mode. Use old suffix logic. */
	fullSuffix[0] = '\0';
	if(verticalCorrectionSuffix != NULL)
	{
		appendSuffix(suffix[deltaB], ".", tmpSuffix);
		appendSuffix(tmpSuffix, verticalCorrectionSuffix, fullSuffix);
		return (appendSuffix(rParamsFile, fullSuffix, newFile));
	}
	return (appendSuffix(rParamsFile, suffix[deltaB], newFile));
}

void getRParams(Offsets *offsets)
{
	FILE *fp;
	double bn, bp, dBn, dBp, rConst, dBnQ, dBpQ;
	int32_t nBaselines;
	double dum;
	int32_t ci;
	int32_t lineCount = 0, eod, special;
	int32_t i, j;
	char line[256], *tmp, paramFile[2048];
	//error("STOP HERE..");
	if((strlen(offsets->rParamsFile)) > 2040) 
	{
		error("filname: %s exceeds 1800 characters");		
	}	
	//error("STOP HERE");
	//  Input parm info
	RgOffsetsParamName(offsets->rParamsFile, paramFile, offsets->deltaB, offsets->verticalCorrectionSuffix);

	/* Case 1: paramFile explicitly ends in .yaml → must use yaml reader */
	{
		size_t plen = strlen(paramFile);
		if (plen > 5 && strcmp(paramFile + plen - 5, ".yaml") == 0) {
			fp = fopen(paramFile, "r");
			if (fp == NULL && strcmp(paramFile, offsets->rParamsFile) != 0) {
				/* The correctly mode-derived name doesn't exist -- e.g. an older
				   frame only ever processed for one deltaB mode. Fall back to the
				   name exactly as given (RgOffsetsParamName()'s pre-fix behavior)
				   rather than erroring just because the re-derived name is missing. */
				fprintf(stderr, "getRParams: %s not found, falling back to %s\n", paramFile, offsets->rParamsFile);
				fp = fopen(offsets->rParamsFile, "r");
			}
			if (fp == NULL)
				error("getRParams: cannot open yaml baseline file %s (also tried %s)", paramFile, offsets->rParamsFile);
			readRParamsYaml(fp, offsets);
			fclose(fp);
			return;
		}
	}

	/* Case 2: try old text format */
	fp = fopen(paramFile, "r");
	if (fp == NULL)
	{
		/* Case 3: auto-detect yaml fallback — try two names:
		     3a: paramFile + ".yaml" (e.g. rBaseline.deltabp.yaml, compatible with old suffix)
		     3b: rParamsFile + canonical yaml infix (e.g. rBaseline.deltabp.yaml per plan) */
		{
			const char *canonInfix[3] = {".yaml", ".deltabp.yaml", ".deltabquad.yaml"};
			int dBcano = (offsets->deltaB >= 0 && offsets->deltaB <= 2) ? offsets->deltaB : 0;
			char yamlName[2200];

			/* 3a: old suffix + .yaml */
			snprintf(yamlName, sizeof(yamlName), "%s.yaml", paramFile);
			fp = fopen(yamlName, "r");
			if (fp == NULL) {
				/* 3b: base name + canonical infix (with optional vcSuffix) */
				if (offsets->verticalCorrectionSuffix != NULL && dBcano != DELTABNONE) {
					/* e.g. rBaseline.deltabp.480.yaml */
					char base[2048];
					snprintf(base, sizeof(base), "%s%s", offsets->rParamsFile, canonInfix[dBcano]);
					size_t blen2 = strlen(base);
					/* insert vcSuffix before .yaml */
					snprintf(yamlName, sizeof(yamlName), "%.*s.%s.yaml",
					         (int)(blen2 - 5), base, offsets->verticalCorrectionSuffix);
				} else {
					snprintf(yamlName, sizeof(yamlName), "%s%s",
					         offsets->rParamsFile, canonInfix[dBcano]);
				}
				fp = fopen(yamlName, "r");
			}
			if (fp != NULL) {
				fprintf(stderr, "getRParams: using yaml fallback %s\n", yamlName);
				readRParamsYaml(fp, offsets);
				fclose(fp);
				return;
			}
		}
		/* Case 4: revert to DELTABNONE base file (existing behavior) */
		fp = openInputFile(offsets->rParamsFile);
		offsets->deltaB = DELTABNONE;
	}
	/*fprintf(stderr,"deltaB %i\n",offsets->deltaB);*/
	for (i = 1; i <= 6; i++)
		for (j = 1; j <= 6; j++)
			offsets->Cr[i][j] = 0.0;
	/*
	  Skip past initial data lines
	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%i", &nBaselines) != 1)
		error("getRparams -- Missing baseline params at line %i: of %s", lineCount, paramFile);
	if (nBaselines < 1 || nBaselines > 3)
		error("getRParams -- invalid number of baselines at line %i of %s\n", lineCount, paramFile);

	if (nBaselines > 1)
		lineCount = getDataString(fp, lineCount, line, &eod);
	if (nBaselines > 2)
		lineCount = getDataString(fp, lineCount, line, &eod);
	/*
	  Read covariance matrix if there is one.
	*/
	readCov(fp, 6, offsets->Cr, &(offsets->sigmaRresidual), line);
	/* for(i=1; i <=6; i++) fprintf(stderr,"%le %le %le %le %le %le \n",
		(offsets->Cr[i][1]),(offsets->Cr[i][2]),(offsets->Cr[i][3]),(offsets->Cr[i][4]),(offsets->Cr[i][5]),(offsets->Cr[i][6]));*/
	fprintf(stderr, "%srange sigma*sqrt(X2/n) = %lf (m)\n", indentRegionOutput ? "\t" : "", offsets->sigmaRresidual);
	/*
	  Input baseline estimated with tiepoints.
	*/
	dBpQ = 0.0;
	dBnQ = 0.0;
	rConst = LARGEINT;
	if (sscanf(line, "%lf%lf%lf%lf%lf%lf%lf", &bn, &bp, &dBn, &dBp, &rConst, &dBnQ, &dBpQ) != 7)
	{
		if (sscanf(line, "%lf%lf%lf%lf%lf", &bn, &bp, &dBn, &dBp, &rConst) != 5)
			error("\n\ngetRoffsets:Invalid range offset baseline file\nFile: %s\nLine: %s", offsets->rParamsFile, line);
	}
	offsets->bn = bn;
	offsets->bp = bp;
	offsets->dBn = dBn;
	offsets->dBp = dBp;
	/* This parameter will get calculated in the SV basline init routine if its used, so only set for computed baseline */
	//offsets->rConst = 0.; this was overwriting earlier values
	if (offsets->deltaB == DELTABNONE)
	{
		if(rConst < LARGEINT)
		{
			offsets->rConst = rConst;
		}
		else
		{
			error("getRparams: rCoonst not initialized");
		}
	}
	offsets->dBnQ = dBnQ;
	offsets->dBpQ = dBpQ;

	/*
	  Scan all remaining lines (including after '&') for a ;* offsetCorrectionFile entry.
	  Format: ;* offsetCorrectionFile <path>
	  Only stores the path — caller is responsible for loading the data.
	*/
	offsets->rOffCorrection.correctionFile[0] = '\0';
	{
		char sc_line[256], *tmp2;
		while (fgets(sc_line, sizeof(sc_line), fp) != NULL)
		{
			/* Must start with ';' and contain '*' to be a special line */
			if (sc_line[0] != ';') continue;
			if (strchr(sc_line, '*') == NULL) continue;
			tmp2 = strstr(sc_line, "offsetCorrectionFile");
			if (tmp2 != NULL)
			{
				tmp2 += strlen("offsetCorrectionFile");
				while (*tmp2 == ' ' || *tmp2 == '\t') tmp2++;
				strncpy(offsets->rOffCorrection.correctionFile, tmp2,
				        sizeof(offsets->rOffCorrection.correctionFile) - 1);
				offsets->rOffCorrection.correctionFile[sizeof(offsets->rOffCorrection.correctionFile) - 1] = '\0';
				/* strip trailing whitespace/newline */
				int32_t len = strlen(offsets->rOffCorrection.correctionFile);
				while (len > 0 && (offsets->rOffCorrection.correctionFile[len-1] == '\n' ||
				                   offsets->rOffCorrection.correctionFile[len-1] == '\r' ||
				                   offsets->rOffCorrection.correctionFile[len-1] == ' '))
					offsets->rOffCorrection.correctionFile[--len] = '\0';
				fprintf(stderr, "getRParams: found offsetCorrectionFile %s\n",
				        offsets->rOffCorrection.correctionFile);
				break; /* found it, no need to continue */
			}
		}
	}

	//fprintf(stderr, "bn %f %f %f bp %f %f %f off %f\n", offsets->bn, offsets->dBn,
	//	offsets->dBnQ, offsets->bp, offsets->dBp, offsets->dBpQ, offsets->rConst);
	//error("STOP %s\n", paramFile);
	fclose(fp);
}


/* 2026-06-17: parse yaml output written by azparams -yaml */
static void readAzParamsYaml(FILE *fp, Offsets *offsets)
{
	char line[512];
	int inC = 0, ci = 0;
	int i, j;

	offsets->sigmaAresidual = 0.0;
	offsets->c1 = 0.0;
	offsets->dbcds = 0.0;
	offsets->dbhds = 0.0;
	offsets->doffdx = 0.0;
	for (i = 1; i <= 4; i++)
		for (j = 1; j <= 4; j++)
			offsets->Ca[i][j] = 0.0;

	while (fgets(line, sizeof(line), fp)) {
		if      (sscanf(line, "sigma: %lf",  &offsets->sigmaAresidual) == 1) { inC = 0; }
		else if (sscanf(line, "cnst: %lf",   &offsets->c1)             == 1) { inC = 0; }
		else if (sscanf(line, "dbcds: %lf",  &offsets->dbcds)          == 1) { inC = 0; }
		else if (sscanf(line, "dbhds: %lf",  &offsets->dbhds)          == 1) { inC = 0; }
		else if (sscanf(line, "doffdx: %lf", &offsets->doffdx)         == 1) { inC = 0; }
		else if (strncmp(line, "C:", 2) == 0) { inC = 1; ci = 0; }
		else if (inC && strstr(line, "- [") && ci < 4) {
			char *p = strstr(line, "[");
			if (p)
				sscanf(p + 1, "%lf, %lf, %lf, %lf",
				       &offsets->Ca[ci+1][1], &offsets->Ca[ci+1][2],
				       &offsets->Ca[ci+1][3], &offsets->Ca[ci+1][4]);
			ci++;
		}
		else { inC = 0; }
	}
	fprintf(stderr, "azimuth sigma*sqrt(X2/n) = %lf (m)\n", offsets->sigmaAresidual);
}

static char *AzOffsetsParamName(char *aParamsFile, char *newFile, int32_t deltaB)
{
	char *suffix[3] = {"", ".const", ".svlinear"};
	size_t blen;

	if (deltaB > DELTABQUAD || deltaB < DELTABNONE)
		error("invalid deltaB flag %i", deltaB);

	/* 2026-06-17: if aParamsFile ends in .yaml, insert mode suffix before .yaml */
	blen = strlen(aParamsFile);
	if (blen > 5 && strcmp(aParamsFile + blen - 5, ".yaml") == 0) {
		strncpy(newFile, aParamsFile, blen - 5);
		newFile[blen - 5] = '\0';
		strcat(newFile, suffix[deltaB]);
		strcat(newFile, ".yaml");
		return newFile;
	}

	return (appendSuffix(aParamsFile, suffix[deltaB], newFile));
}

/*
  Input azimuth parameter info. Modifed 3/1/16 to read in covariance matrix
*/
void getAzParams(Offsets *offsets)
{
	FILE *fp;
	int32_t lineCount = 0, eod;
	int32_t i, j;
	char line[256], paramFile[1024];

	AzOffsetsParamName(offsets->azParamsFile, paramFile, offsets->deltaB);

	/* 2026-06-17: Case 1: paramFile explicitly ends in .yaml → must use yaml reader */
	{
		size_t plen = strlen(paramFile);
		if (plen > 5 && strcmp(paramFile + plen - 5, ".yaml") == 0) {
			fp = fopen(paramFile, "r");
			if (fp == NULL)
				error("getAzParams: cannot open yaml az params file %s", paramFile);
			readAzParamsYaml(fp, offsets);
			fclose(fp);
			return;
		}
	}

	/* Case 2: try old text format */
	fp = fopen(paramFile, "r");
	if (fp == NULL)
	{
		/* Case 3: auto-detect yaml fallback — try paramFile + ".yaml" */
		{
			char yamlName[1200];
			snprintf(yamlName, sizeof(yamlName), "%s.yaml", paramFile);
			fp = fopen(yamlName, "r");
			if (fp != NULL) {
				fprintf(stderr, "getAzParams: using yaml fallback %s\n", yamlName);
				readAzParamsYaml(fp, offsets);
				fclose(fp);
				return;
			}
		}
		/* Case 4: revert to DELTABNONE base file */
		fp = openInputFile(offsets->azParamsFile);
		offsets->deltaB = DELTABNONE;
		if (fp == NULL)
			error("getAzParams: Error opening %s\n", offsets->azParamsFile);
	}
	/*
	  Skip past initial data lines
	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	/*
	  Read covariance matrix if there is one.
	*/
	for (i = 1; i <= 4; i++)
		for (j = 1; j <= 4; j++)
			offsets->Ca[i][j] = 0.0;
	readCov(fp, 4, offsets->Ca, &(offsets->sigmaAresidual), line);
	fprintf(stderr, "azimuth sigma*sqrt(X2/n) = %lf (m)\n", offsets->sigmaAresidual);
	/* read params from fit */
	if (sscanf(line, "%lf%lf%lf%lf", &(offsets->c1), &(offsets->dbcds), &(offsets->dbhds), &(offsets->doffdx)) != 4)
	{
		error("getAzParams  -- Missing baseline params at line: %s %i\n %s", paramFile, lineCount, line);
	}
	fclose(fp);
}

void initOffParams(Offsets *offsets){
	offsets->nr = 0;
	offsets->na = 0;
	offsets->rO = 0;
	offsets->aO = 0;
	offsets->deltaR = 0;
	offsets->deltaA =0;
	offsets->sigmaStreaks = 0.0;
	offsets->sigmaRange = 0.0;
	offsets->geo1 = NULL;
	offsets->geo2 = NULL;
	offsets->rOffCorrection.rangeOffsetCorrection = NULL;
	offsets->rOffCorrection.correctionFile[0] = '\0';
	/* 2026-06-24: bnS/bpS are lazily computed by svInitBnBp() (see svInterpBnBp()'s
	   "if bnS==NULL" check), indexed according to aO/deltaA/na at the time they were
	   computed. When the same Offsets struct gets reused across multiple different
	   pairings within one mosaic3d run (e.g. make3DOffsets.c trying several candidate
	   crossing-orbit images against the same descending frame), aO/deltaA/na get
	   correctly refreshed here on every call, but bnS/bpS previously didn't -- so the
	   lazy-init only ever fired on the first pairing a given image was visited with,
	   and every later pairing silently reused a baseline array indexed for the wrong
	   azimuth window. Freeing them here forces a fresh, correctly-indexed
	   recomputation for every pairing. */
	free(offsets->bnS);
	free(offsets->bpS);
	offsets->bnS = NULL;
	offsets->bpS = NULL;
}

static GDALRasterBandH getBandAndMeta(GDALDatasetH hDS, Offsets *offsets, int32_t band, char *path)
{
	dictNode *metaData = NULL;
	GDALRasterBandH hBand;
	float tmp;
	// Get meta data
	hBand = GDALGetRasterBand(hDS, band);
	
	readDataSetMetaData(hDS, &metaData);
	// Write to offsets
	offsets->nr = GDALGetRasterBandXSize(hBand);
	offsets->na = GDALGetRasterBandYSize(hBand);
	offsets->rO = atoi(get_value(metaData, "r0"));
	offsets->aO = atoi(get_value(metaData, "a0"));
	// fprintf(stderr, "R0,A)\n");
	offsets->deltaR = atof(get_value(metaData, "deltaR"));
	offsets->deltaA = atof(get_value(metaData, "deltaA"));
	// Check sigmas not already set
	if( (tmp=atof(get_value(metaData, "sigmaStreaks"))) > 0.0) offsets->sigmaStreaks = tmp;
	if( (tmp=atof(get_value(metaData, "sigmaRange"))) > 0.0) offsets->sigmaRange = tmp;
	// fprintf(stderr,"PATH %s\n", path);
	offsets->geo1 = mergePath(get_value(metaData, "geo1"), path);
	offsets->geo2 = mergePath(get_value(metaData, "geo2"), path);
	fprintf(stderr, "%sGEO2 %s\n", indentRegionOutput ? "\t" : "", offsets->geo2);
	// Get Band
	return hBand;
}


static void mapBandDescriptionsToBandNumbers(GDALDatasetH hDS, int32_t bandNumbers[5])
{	
	/* 
	return band numbers for in bandnumbers as 
	bandNumbers[1] band number for azimuth offsets
	bandNumbers[2] band number for range offsets
	bandNumbers[3] band number for azimuth errors
	bandNumbers[4] band number for range errors
	*/
	const char *description;
	GDALRasterBandH hBand;
	for(int i=1; i <= 4; i++) bandNumbers[i] = 0;
	int32_t nBands = GDALGetRasterCount(hDS);
	// fprintf(stderr, "Number of bands in offset file %i\n", nBands);
	if(nBands > 4) error("mapBandDescriptionsToBandNumbers: to many (%i) bands for offset file\n", nBands);
	for(int i=1; i <= nBands; i++) {
		hBand = GDALGetRasterBand(hDS, i);
		// fprintf(stderr, "Band %i\n", i);
		description = GDALGetMetadataItem(hBand, "Description", NULL);
		if(description == NULL)
		{
			description = GDALGetDescription(hBand);
		}
		if (description == NULL || strlen(description) == 0) error("No description found for this band.\n");
		//fprintf(stderr, "AfterBand %i %s\n", i, description);
		// fprintf(stderr, "Description %s %i\n", description, i);
		if(strstr(description, "AzimuthOffsets") != NULL) bandNumbers[1] = i;
		else if(strstr(description, "AzimuthSigma")  != NULL) bandNumbers[3] = i;
		else if(strstr(description, "RangeOffsets") != NULL) bandNumbers[2] = i;
		else if(strstr(description, "RangeSigma") != NULL) bandNumbers[4] = i;
		else error("mapBandDescriptionsToBandNumbers: invalid band name (%s) for offset file %s\n", description);
		//fprintf(stderr, "%i, %i %i %i %i\n", i, bandNumbers[1], bandNumbers[2], bandNumbers[3], bandNumbers[4]);
	}
	// fprintf(stderr, "%i %i %i %i\n", bandNumbers[1], bandNumbers[2], bandNumbers[3], bandNumbers[4]);
}

int32_t readGDALOffsets(GDALDatasetH hDS, Offsets *offsets, int bufferMode, float azimuthMin, float azimuthMax)
{
	int32_t status;
	float *data;
	char *path, buf[2048];
	int32_t bandNumbers[5];
	GDALRasterBandH hBand;
	int32_t iAzMin, iAzMax, nRows, k;
	int32_t maskApplied = FALSE;

	mapBandDescriptionsToBandNumbers(hDS, bandNumbers);
	// Handle various buffer cases
	buf[0] = '\0';
	switch (bufferMode)
	{
	case AZIMUTHBUFF:
		// fprintf(stderr, "AZIMUTH BUFF\n");
		path = dirname(strcpy(buf, offsets->file));
		hBand = getBandAndMeta(hDS, offsets, bandNumbers[1], path);
		initOffsetBuffers(offsets, AZONLY);
		data = offsets->da[0];
		break;
	case RANGEBUFF:
		fprintf(stderr, "%sRANGE BUFF  %s\n", indentRegionOutput ? "\t" : "", offsets->rFile);
		path = dirname(strcpy(buf, offsets->rFile));
		hBand = getBandAndMeta(hDS, offsets, bandNumbers[2], path);
		initOffsetBuffers(offsets, RGONLY);
		data = offsets->dr[0];
		break;
	case RANGEUSEAZIMUTHBUFF:
		fprintf(stderr, "%sRANGE FLIP BUFF %s\n", indentRegionOutput ? "\t" : "", offsets->rFile);
		path = dirname(strcpy(buf, offsets->rFile));
		hBand = getBandAndMeta(hDS, offsets, bandNumbers[2], path);
		initOffsetBuffers(offsets, AZFORRANGE);
		data = offsets->dr[0];
		break;
	case AZIMUTHERRORBUFF:
		// fprintf(stderr, "AZIMUTH ERROR BUFF\n")
		hBand = GDALGetRasterBand(hDS, bandNumbers[3]);
		data = offsets->sa[0];
		break;
	case RANGEERRORBUFF:
		// fprintf(stderr, "RANGE ERROR BUFF\n");
		hBand = GDALGetRasterBand(hDS, bandNumbers[4]);
		data = offsets->sr[0];
		break;
	default:
		error("Invalide code readGDALoffstes");
	}

	/* Clamp azimuth range to valid rows */
	azimuthMin = (azimuthMin - offsets->aO) / offsets->deltaA;
	azimuthMax = (azimuthMax - offsets->aO) / offsets->deltaA;
	iAzMin = max(0, (int32_t)azimuthMin);
	iAzMax = min(offsets->na - 1, (int32_t)azimuthMax);
	nRows  = iAzMax - iAzMin + 1;

	if (iAzMin == 0 && iAzMax == offsets->na - 1)
	{
		/* Full image read */
		fprintf(stderr, "%s\033[32mFull read of %i rows and %i columns\033[0m\n", indentRegionOutput ? "\t" : "", offsets->na, offsets->nr);
		status = GDALRasterIO(hBand, GF_Read, 0, 0, offsets->nr, offsets->na, data,
							  offsets->nr, offsets->na, GDT_Float32, 0, 0);
	}
	else
	{
		/* Partial read: fill entire buffer with -LARGEINT, then read desired rows
		   in-place so row indices match those of a full read */
		for (k = 0; k < offsets->nr * offsets->na; k++)
			data[k] = (float)-LARGEINT;
		fprintf(stderr, "%s\033[34mPartial read of rows %i to %i (of %i) and %i columns\033[0m\n", indentRegionOutput ? "\t" : "", iAzMin, iAzMax, offsets->na, offsets->nr);
		status = GDALRasterIO(hBand, GF_Read,
							  0, iAzMin, offsets->nr, nRows,
							  data + iAzMin * offsets->nr, offsets->nr, nRows,
							  GDT_Float32, 0, 0);
	}

	if (status != CE_None)
	{
		error("readGDALOffsets: GDALRasterIO failed\n");
	}

	/* Convert any NaN pixels (source nodata) to -LARGEINT sentinel */
	for (k = 0; k < offsets->nr * offsets->na; k++)
	{
		if (isnan(data[k])) data[k] = (float)-LARGEINT;
	}

	/* Apply an embedded VRT dataset mask band (e.g., range.offsets.good.tif/
	   azimuth.offsets.good.tif from autocleanNISAR.py), if present, to the
	   primary value bands only -- downstream code gates validity purely on
	   dr[i][j]/da[i][j] > -LARGEINT, so the sigma bands don't need masking. */
	if (!noMask && (bufferMode == AZIMUTHBUFF || bufferMode == RANGEBUFF || bufferMode == RANGEUSEAZIMUTHBUFF))
	{
		int32_t maskFlags = GDALGetMaskFlags(hBand);
		if (!(maskFlags & GMF_ALL_VALID))
		{
			GDALRasterBandH hMaskBand = GDALGetMaskBand(hBand);
			GByte *maskData = (GByte *)malloc(offsets->nr * offsets->na * sizeof(GByte));
			if (maskData == NULL)
			{
				error("readGDALOffsets: failed to allocate mask buffer\n");
			}
			if (iAzMin == 0 && iAzMax == offsets->na - 1)
			{
				status = GDALRasterIO(hMaskBand, GF_Read, 0, 0, offsets->nr, offsets->na, maskData,
									  offsets->nr, offsets->na, GDT_Byte, 0, 0);
			}
			else
			{
				for (k = 0; k < offsets->nr * offsets->na; k++)
				{
					maskData[k] = 1;
				}
				status = GDALRasterIO(hMaskBand, GF_Read,
									  0, iAzMin, offsets->nr, nRows,
									  maskData + iAzMin * offsets->nr, offsets->nr, nRows,
									  GDT_Byte, 0, 0);
			}
			if (status != CE_None)
			{
				error("readGDALOffsets: GDALRasterIO failed reading mask band\n");
			}
			for (k = 0; k < offsets->nr * offsets->na; k++)
			{
				if (maskData[k] == 0) data[k] = (float)-LARGEINT;
			}
			free(maskData);
			maskApplied = TRUE;
		}
	}
	return maskApplied;
}

/*
 This combines funtionality of historical readOffsets and readAzimuthOffsets.
*/
void readOffsetsOptionalErrors(Offsets *offsets, int32_t includeErrors, float azimuthMin, float azimuthMax)
{
	char *datFile, buf[1024], bufa[2048], vrtBuffer[2048], *vrtFile;
	char *eFileA, *file;
	GDALDatasetH hDS;

	//fprintf(stderr, "]n\noffsets file %s\n\n", offsets->file);
	if (has_suffix(offsets->file, ".vrt") == TRUE)
		vrtFile = strcpy(vrtBuffer, offsets->file);
	else
		vrtFile = checkForOffsetsVrt(offsets->file, vrtBuffer);
	fprintf(stderr, "vrtFile %s\n", vrtFile);
	//vrtFile = checkForOffsetsVrt(offsets->file, vrtBuffer);
	if (vrtFile != NULL)
	{	// Zero params — but preserve any correctionFile pre-set by caller
		char savedCorrFile2[2048];
		strncpy(savedCorrFile2, offsets->rOffCorrection.correctionFile, sizeof(savedCorrFile2) - 1);
		savedCorrFile2[sizeof(savedCorrFile2) - 1] = '\0';
		initOffParams(offsets);
		if (savedCorrFile2[0] != '\0')
			strncpy(offsets->rOffCorrection.correctionFile, savedCorrFile2,
			        sizeof(offsets->rOffCorrection.correctionFile));
		//fprintf(stderr, "OPENING VRT %s\n", vrtFile);
		// Open data set
		hDS = GDALOpen(vrtFile, GDAL_OF_READONLY);
		/*  -noAzimuthRows under the hopper: the buffers are still allocated (initOffParams/
		    readGDALOffsets set the geometry), but the rasters are not read.  Safe only because
		    the same flag forces useAzRow FALSE, so nothing reads offsets->da afterwards.  This is
		    pure I/O saving -- on a NISAR frame the azimuth band is the same size as the range one. */
		if (skipAzimuthOffsets == FALSE)
		{
			// Read azimuthg offsets and errors
			readGDALOffsets(hDS, offsets, AZIMUTHBUFF, azimuthMin, azimuthMax);
			if(includeErrors == TRUE)
				readGDALOffsets(hDS, offsets, AZIMUTHERRORBUFF, azimuthMin, azimuthMax);
		}
		GDALClose(hDS);
	}
	else
	{
		datFile = appendSuffix(offsets->file, ".dat", buf);
		//fprintf(stderr, "OPENING DAT %s\n", datFile);
		// Read the data file
		readOffsetParams(datFile, offsets, TRUE);
		// Init memory
		initOffsetBuffers(offsets, AZONLY);
		// Read files
		if (skipAzimuthOffsets == FALSE)
		{
			readOffsetFile(offsets->da, offsets->nr, offsets->na, offsets->file);
			if(includeErrors == TRUE) {
				eFileA = appendSuffix(offsets->file, ".sa", bufa);
				readOffsetFile(offsets->sa, offsets->nr, offsets->na, eFileA);
			}
		}
	}
}

/*
   Read azimuth offsets with errors
 */
void readOffsets(Offsets *offsets) {
	readOffsetsOptionalErrors(offsets, TRUE, 0.0f, (float)LARGEINT);
}

/*
   Read azimuth offsets only
 */
void readAzimuthOffsets(Offsets *offsets) {
	readOffsetsOptionalErrors(offsets, FALSE, 0.0f, (float)LARGEINT);
}


static void checkForIonosphereCorrection(GDALDatasetH hDS, Offsets *offsets, int bufferMode)
{
	dictNode *metaData = NULL;
	char ionospherePath[2048];
	char tmp[2048];

	/* Only load if correctionFile was pre-filled (either by getRParams reading a baseline
	   that recorded it, or by getROffsets peeking at the VRT before readRangeOffsets).
	   This preserves consistency: if the baseline was estimated without correction,
	   we don't retroactively apply one. */
	//error("XXX: checkForIonosphereCorrection not fully implemented yet %s\n", offsets->rOffCorrection.correctionFile);
	if (offsets->rOffCorrection.correctionFile[0] == '\0')
		return;

	readDataSetMetaData(hDS, &metaData);
	char *ionosphereCorrection = get_value(metaData, "ionosphereRangeOffsetCorrection");
	fprintf(stderr, "%sReading Ionosphere correction file: %s\n", indentRegionOutput ? "\t" : "", ionosphereCorrection);
	
	/* Baseline specified a file but VRT metadata has none */
	if (ionosphereCorrection == NULL)
		error("checkForIonosphereCorrection: baseline file specifies offsetCorrectionFile '%s' "
		      "but VRT has no ionosphereRangeOffsetCorrection entry\n",
		      offsets->rOffCorrection.correctionFile);

	/* Compare basenames only */
	strncpy(tmp, offsets->rOffCorrection.correctionFile, sizeof(tmp) - 1);
	char *baselineBase = basename(tmp);
	if (strcmp(baselineBase, ionosphereCorrection) != 0)
		error("checkForIonosphereCorrection: baseline offsetCorrectionFile basename '%s' "
		      "does not match VRT ionosphereRangeOffsetCorrection '%s'\n",
		      baselineBase, ionosphereCorrection);

	/* Names match — build full path and load */
	strncpy(tmp, offsets->rFile, sizeof(tmp) - 1);
	snprintf(ionospherePath, sizeof(ionospherePath), "%s/%s", dirname(tmp), ionosphereCorrection);
	if (access(ionospherePath, F_OK) != 0)
		error("checkForIonosphereCorrection: correction file not found: %s\n", ionospherePath);
	//error("STOP: ionospheric correction loading not implemented yet\n");
	readOffsetCorrection(ionospherePath, offsets, bufferMode);
}

/*
   Read range  offsets (for now no sigma)
 */
void readRangeOffsets(Offsets *offsets, int32_t includeErrors, float azimuthMin, float azimuthMax)
{
	char *datFile, buf[2048], bufd[2048], vrtBuffer[2048], *vrtFile, SEPath[2048];
	char *eFileR, *file;
	float *data;
	GDALDatasetH hDS;
	/*
	  Read inputfile
	*/ 
	if (has_suffix(offsets->rFile, ".vrt") == TRUE)
		vrtFile = strcpy(vrtBuffer, offsets->rFile);
	else
		vrtFile = checkForOffsetsVrt(offsets->rFile, vrtBuffer);
	fprintf(stderr, "vrtFile %s\n", vrtFile);
	if (vrtFile != NULL)
	{	// Zero parameters — but preserve any correctionFile pre-set by caller
		char savedCorrFile[2048];
		strncpy(savedCorrFile, offsets->rOffCorrection.correctionFile, sizeof(savedCorrFile) - 1);
		savedCorrFile[sizeof(savedCorrFile) - 1] = '\0';
		initOffParams(offsets);
		if (savedCorrFile[0] != '\0')
			strncpy(offsets->rOffCorrection.correctionFile, savedCorrFile,
			        sizeof(offsets->rOffCorrection.correctionFile));
		fprintf(stderr, "OPENING VRT %s\n", vrtFile);
		// Open data set
		hDS = GDALOpen(vrtFile, GDAL_OF_READONLY);
		fprintf(stderr, "VRT opened\n");
		// Read data and close
		readGDALOffsets(hDS, offsets, RANGEBUFF, azimuthMin, azimuthMax);
		fprintf(stderr, "Range offsets read from VRT xxx\n");
		if(includeErrors == TRUE)
			readGDALOffsets(hDS, offsets, RANGEERRORBUFF, azimuthMin, azimuthMax);
		// Check for ionospheric correction
		fprintf(stderr, "Checking for ionospheric correction\n");
		checkForIonosphereCorrection(hDS, offsets, RANGEBUFF);
		//fprintf(stderr, "%s\n", offsets->rFile);
		fprintf(stderr, "Range offsets read from VRT\n");
		// Check  if intermediate product that needs SE correction
		// This is kluge, which should only be invoked if creating velocity_nocull product from intermediate offset products.
		if(is_not_final_range(vrtFile) && get_se_correction_path(offsets->geo1, SEPath, sizeof(SEPath))) {
			fprintf(stderr, "\033[1;34mWARNING: applying SE correction to %s which is an intermediate file.\n\033[0m", vrtFile);
			initOffsetBuffers(offsets, SEBUF);
			data =  offsets->SECorrection[0];
			int xSize, ySize, dataType;
			dictNode *metaDictionary = NULL;
			readRasterVRT(SEPath, 1, &xSize, &ySize,  &dataType, &metaDictionary, data, 0, (int32_t)LARGEINT);
			// Apply SE correction to range offsets
			for(int i=0; i < ySize; i++)
			{
				for(int j=0; j < xSize; j++)
				{
					if(offsets->dr[i][j] > -LARGEINT)
						offsets->dr[i][j] -= offsets->SECorrection[i][j];
				}
			}
		} 
		GDALClose(hDS);
	}
	else
	{
		// Zero parameters
		initOffParams(offsets);
		// read offset param
		datFile = appendSuffix(offsets->rFile, ".dat", buf);
		//fprintf(stderr, "OPENING DAT %s\n", datFile);
		readOffsetParams(datFile, offsets, TRUE);
		// setup buffers
		initOffsetBuffers(offsets, RGONLY);
		// Read files
		readOffsetFile(offsets->dr, offsets->nr, offsets->na, offsets->rFile);
		if(includeErrors == TRUE) {
			eFileR = appendSuffix(offsets->rFile, ".sr", bufd);
			readOffsetFile(offsets->sr, offsets->nr, offsets->na, eFileR);
		}
	}
}

/*
	read offsets and error files
*/
void readBothOffsets(Offsets *offsets, float azimuthMin, float azimuthMax)
{
	char *datFile, buf[1024], bufa[1024], bufd[1024], bufvrt[2048];
	char *eFileA, *eFileR;
	char *file;
	int32_t i;
	/*
	  Read azimuth offsets followed by range offsets
	*/
	readOffsetsOptionalErrors(offsets, TRUE, azimuthMin, azimuthMax);
	readRangeOffsets(offsets, TRUE, azimuthMin, azimuthMax);
	fprintf(stderr, "SIGMA FINAL %f %f", offsets->sigmaStreaks, offsets->sigmaRange);
}


/*
 * Returns 1 if:
 *   - filename is exactly "range.offsets.fast"
 *   OR
 *   - filename contains ".interp.dr"
 * Otherwise returns 0.
 */
static int is_not_final_range(const char *vrtFile)
{
    const char *filename;
    const char *p;

    if (vrtFile == NULL)
        return 0;

    /* Find last '/' */
    filename = strrchr(vrtFile, '/');

    /* Also check for Windows '\' */
    p = strrchr(vrtFile, '\\');
    if (p != NULL && (filename == NULL || p > filename))
        filename = p;

    /* Move past separator if found */
    if (filename != NULL)
        filename++;
    else
        filename = vrtFile;
    /* Substring match */
	//fprintf(stderr, "Checking if %s is an intermediate range offset file\n", filename);
    if (strstr(filename, "range.offsets.fast") != NULL)
        return 1;

    /* Substring match */
    if (strstr(filename, "cull.interp") != NULL)
        return 1;

    return 0;
}


/*
 * If offsets.SECorrection.vrt exists in same directory as vrtFile:
 *   - writes full path into outPath
 *   - returns 1
 * Otherwise:
 *   - returns 0
 */
int get_se_correction_path(const char *vrtFile,
                           char *outPath,
                           size_t outSize)
{
    const char *filename;
    const char *p;
    size_t dir_len;
    FILE *fp;

    if (!vrtFile || !outPath || outSize == 0)
        return 0;

    /* Find last path separator */
    filename = strrchr(vrtFile, '/');
    p = strrchr(vrtFile, '\\');
    if (p && (!filename || p > filename))
        filename = p;

    if (filename)
        dir_len = filename - vrtFile + 1;  /* include separator */
    else
        dir_len = 0;

    if (dir_len + strlen("offsets.SECorrection.vrt") + 1 > outSize)
        return 0;

    /* Copy directory */
    if (dir_len > 0)
    {
        strncpy(outPath, vrtFile, dir_len);
        outPath[dir_len] = '\0';
    }
    else
    {
        outPath[0] = '\0';
    }

    /* Append filename */
    strcat(outPath, "offsets.SECorrection.vrt");

    /* Check existence */
    fp = fopen(outPath, "r");
    if (fp)
    {
        fclose(fp);
        return 1;
    }

    return 0;
}


/*
  This reads the range offsets, but uses the azimuth  offsets buffer for the asc and the range for the descending
*/
void readRangeOrRangeOffsets(Offsets *offsets, int32_t orbitType, float azimuthMin, float azimuthMax)
{
	char *datFile, buf[2048], bufd[2048], vrtBuffer[2048], *vrtFile;
	char *eFileR, *file;
	GDALDatasetH hDS;
	int bufferMode;
	// Zero parameters — but preserve any correctionFile pre-set by caller
	char savedCorrFile3[2048];
	strncpy(savedCorrFile3, offsets->rOffCorrection.correctionFile, sizeof(savedCorrFile3) - 1);
	savedCorrFile3[sizeof(savedCorrFile3) - 1] = '\0';
	initOffParams(offsets);
	if (savedCorrFile3[0] != '\0')
		strncpy(offsets->rOffCorrection.correctionFile, savedCorrFile3,
		        sizeof(offsets->rOffCorrection.correctionFile));
	if (has_suffix(offsets->rFile, ".vrt") == TRUE)
		vrtFile = strcpy(vrtBuffer, offsets->rFile);
	else
		vrtFile = checkForOffsetsVrt(offsets->rFile, vrtBuffer);
	//fprintf(stderr, "vrtFile %s\n", vrtFile);
	//vrtFile = checkForOffsetsVrt(offsets->rFile, vrtBuffer);
	if (vrtFile != NULL)
	{
		//fprintf(stderr, "OPENING VRT %s\n", vrtFile);
		// Open data set
		double tBlockStart = now();
		double t0 = tBlockStart;
		hDS = GDALOpen(vrtFile, GDAL_OF_READONLY);
		fprintf(stderr, "%sGDALOpen range vrt time: %.3f s\n", indentRegionOutput ? "\t" : "", now() - t0);
		// Read data and close
		t0 = now();
		int32_t maskApplied;
		if (orbitType == ASCENDING) maskApplied = readGDALOffsets(hDS, offsets, RANGEUSEAZIMUTHBUFF, azimuthMin, azimuthMax);
		else maskApplied = readGDALOffsets(hDS, offsets, RANGEBUFF, azimuthMin, azimuthMax);
		readGDALOffsets(hDS, offsets, RANGEERRORBUFF, azimuthMin, azimuthMax);
		fprintf(stderr, "%sRead range+sigma%s bands time: %.3f s\n", indentRegionOutput ? "\t" : "",
				maskApplied ? "+mask" : "", now() - t0);
		t0 = now();
		checkForIonosphereCorrection(hDS, offsets,
									orbitType == ASCENDING ? RANGEUSEAZIMUTHBUFF : RANGEBUFF);
		fprintf(stderr, "%sIonosphere correction check/read time: %.3f s\n", indentRegionOutput ? "\t" : "", now() - t0);
		fprintf(stderr, "%sTotal range read time: %.3f s\n", indentRegionOutput ? "\t" : "", now() - tBlockStart);
		totalOffsetsIOTime += now() - tBlockStart;
		GDALClose(hDS);
	}
	else
	{
		double tBlockStart = now();
		// read offset param
		datFile = appendSuffix(offsets->rFile, ".dat", buf);
		//fprintf(stderr, "OPENING DAT %s\n", datFile);
		readOffsetParams(datFile, offsets, TRUE);
		// setup buffers
		if (orbitType == ASCENDING)
			initOffsetBuffers(offsets, AZFORRANGE);
		else
			initOffsetBuffers(offsets, RGONLY);
		// Read files
		readOffsetFile(offsets->dr, offsets->nr, offsets->na, offsets->rFile);
		eFileR = appendSuffix(offsets->rFile, ".sr", bufd);
		readOffsetFile(offsets->sr, offsets->nr, offsets->na, eFileR);
		fprintf(stderr, "%sTotal range read time: %.3f s\n", indentRegionOutput ? "\t" : "", now() - tBlockStart);
		totalOffsetsIOTime += now() - tBlockStart;
	}
}

/*
  Input phase or power image for geocode
*/
void getMosaicInputImage(inputImageStructure *inputImage, int32_t yMin, int32_t yMax)
{
	FILE *fp;
	float **fimage;
	float *imageLine;
	char vrtBuf[4092], *vrtFile;
	int32_t i, j;
	dictNode *metaOut = NULL;
  	int xSize, ySize, dataType, status;
	double tFuncStart = now();
	/*
	  Open image
	*/
	vrtFile = checkForVrt(inputImage->file, vrtBuf);
	fimage = (float **)inputImage->image;
	if(vrtFile != NULL)
	{
		imageLine = inputImage->image[0];
		fprintf(stderr, "%sReading %s (BAND 1)\n", indentRegionOutput ? "\t" : "", vrtFile);
		double t0 = now();
		readRasterVRT(vrtFile, 1, &xSize, &ySize, &dataType, &metaOut, imageLine,
					  yMin / inputImage->nAzimuthLooks, yMax / inputImage->nAzimuthLooks);
		fprintf(stderr, "%sGDAL read time: %.3f s\n", indentRegionOutput ? "\t" : "", now() - t0);
		fprintf(stderr, "%sazimuthSize %i %i rangeSize %i %i\n", indentRegionOutput ? "\t" : "", inputImage->azimuthSize, ySize,
			inputImage->rangeSize, xSize);
		fprintf(stderr, "%sVRT read\n", indentRegionOutput ? "\t" : "");
		fprintf(stderr, "%sTotal phase read time: %.3f s\n", indentRegionOutput ? "\t" : "", now() - tFuncStart);
		totalPhaseIOTime += now() - tFuncStart;
		return;
	} else
	{
		fprintf(stderr, "%sReading %s (NO VRT FILE)\n", indentRegionOutput ? "\t" : "", inputImage->file);
		if (strstr(inputImage->file, "nophase") == NULL)
		{
			fp = fopen(inputImage->file, "r");
			if (fp == NULL)
				error("*** getPhaseOrPowerImage: Error opening %s ***\n",
					inputImage->file);
		}
		else
		{
			fp = NULL;
		}
		imageLine = fimage[0];
		if (fp != NULL)
		{
			freadBS(imageLine, sizeof(float), inputImage->rangeSize * inputImage->azimuthSize, fp, FLOAT32FLAG);
		}
		else
		{ /* nophase case */
			for (i = 0; i < inputImage->azimuthSize; i++)
				for (j = 0; j < inputImage->rangeSize; j++)
					fimage[i][j] = -2.0e9;
			if (fp != NULL)
			fclose(fp);
		}
		
	}
	fprintf(stderr, "%scompleted \n", indentRegionOutput ? "\t" : "");
	fprintf(stderr, "%sTotal phase read time: %.3f s\n", indentRegionOutput ? "\t" : "", now() - tFuncStart);
	totalPhaseIOTime += now() - tFuncStart;
	return;
}


/*
   Read the ionospheric phase image (radians) that the baseline file named for this image,
   into the pool pointed to by buf, honoring the same partial-read row window as
   getMosaicInputImage(). No-op when this image has no correction.

   Unlike getMosaicInputImage(), this opens the file with GDALOpen directly rather than
   going through checkForVrt(), so a plain .tif works without a sidecar .vrt wrapper --
   the same approach tiepoints uses in tiePoints/getPhases.c.
*/
void getIonospherePhaseImage(inputImageStructure *inputImage, float *buf, int32_t yMin, int32_t yMax)
{
	GDALDatasetH hDS;
	GDALRasterBandH hBand;
	int32_t i, iYMin, iYMax, nRows;
	int64_t k, nPix;
	double tFuncStart = now();

	if (inputImage->ionospherePhase == NULL || inputImage->ionospherePhaseFile[0] == '\0')
	{
		return;
	}
	/* Repoint the row pointers at the shared pool, mirroring setBuffer() in make3DMosaic.c.
	   A NULL buf means "keep the pool setupADImageBuffers already assigned" -- what
	   makeVhMosaic.c needs, since it never calls setBuffer() either. */
	if (buf != NULL)
	{
		for (i = 0; i < inputImage->azimuthSize; i++)
		{
			inputImage->ionospherePhase[i] = &(buf[(int64_t)i * inputImage->rangeSize]);
		}
	}
	else
	{
		buf = inputImage->ionospherePhase[0];
	}
	hDS = GDALOpen(inputImage->ionospherePhaseFile, GDAL_OF_READONLY);
	if (hDS == NULL)
		error("getIonospherePhaseImage: cannot open %s\n", inputImage->ionospherePhaseFile);
	if (GDALGetRasterXSize(hDS) != inputImage->rangeSize ||
		GDALGetRasterYSize(hDS) != inputImage->azimuthSize)
		error("getIonospherePhaseImage: %s is %d x %d but the phase image is %d x %d\n",
			  inputImage->ionospherePhaseFile, GDALGetRasterXSize(hDS), GDALGetRasterYSize(hDS),
			  inputImage->rangeSize, inputImage->azimuthSize);
	hBand = GDALGetRasterBand(hDS, 1);
	if (hBand == NULL)
		error("getIonospherePhaseImage: cannot get band 1 of %s\n", inputImage->ionospherePhaseFile);
	/* Fill with the nodata sentinel, then read only the requested rows in place, so
	   unread rows interpolate as invalid rather than as stale data. */
	nPix = (int64_t)inputImage->rangeSize * inputImage->azimuthSize;
	for (k = 0; k < nPix; k++)
	{
		buf[k] = (float)-LARGEINT;
	}
	iYMin = yMin / inputImage->nAzimuthLooks;
	iYMax = yMax / inputImage->nAzimuthLooks;
	if (iYMin < 0)
	{
		iYMin = 0;
	}
	if (iYMax >= inputImage->azimuthSize)
	{
		iYMax = inputImage->azimuthSize - 1;
	}
	nRows = iYMax - iYMin + 1;
	if (nRows > 0)
	{
		if (GDALRasterIO(hBand, GF_Read, 0, iYMin, inputImage->rangeSize, nRows,
						 &(buf[(int64_t)iYMin * inputImage->rangeSize]),
						 inputImage->rangeSize, nRows, GDT_Float32, 0, 0) != CE_None)
			error("getIonospherePhaseImage: GDALRasterIO failed for %s\n", inputImage->ionospherePhaseFile);
	}
	GDALClose(hDS);
	fprintf(stderr, "%sRead ionosphere phase %s (rows %i to %i): %.3f s\n",
			indentRegionOutput ? "\t" : "", inputImage->ionospherePhaseFile, iYMin, iYMax,
			now() - tFuncStart);
	totalPhaseIOTime += now() - tFuncStart;
}

/*
   Read 4x4 or 6x6 cov matrix for the az/rg params file format
*/
static void readCov(FILE *fp, int32_t n, double C[7][7], double *sigmaResidual, char *line)
{
	int32_t lineCount = 0, eod;
	int32_t count;
	int32_t special, ci;
	char *tmp;
	double fdum1, fdum2, fdum3, fdum4, fdum5, fdum6;
	count = 0;
	special = TRUE;
	while (special == TRUE)
	{
		lineCount = getDataStringSpecial(fp, lineCount, line, &eod, '*', &special);
		ci = 0;
		tmp = strstr(line, "C_1");
		if (tmp != NULL)
		{
			tmp += 3;
			ci = 1;
		}
		if (tmp == NULL)
		{
			tmp = strstr(line, "C_2");
			if (tmp != NULL)
			{
				tmp += 3;
				ci = 2;
			}
		}
		if (tmp == NULL)
		{
			tmp = strstr(line, "C_3");
			if (tmp != NULL)
			{
				tmp += 3;
				ci = 3;
			}
		}
		if (tmp == NULL)
		{
			tmp = strstr(line, "C_4");
			if (tmp != NULL)
			{
				tmp += 3;
				ci = 4;
			}
		}
		if (tmp == NULL)
		{
			tmp = strstr(line, "C_5");
			if (tmp != NULL)
			{
				tmp += 3;
				ci = 5;
			}
		}
		if (tmp == NULL)
		{
			tmp = strstr(line, "C_6");
			if (tmp != NULL)
			{
				tmp += 3;
				ci = 6;
			}
		}
		if (ci > 0)
		{
			fdum1 = 0;
			fdum2 = 0;
			fdum3 = 0.0;
			fdum4 = 0.0;
			fdum5 = 0.0;
			fdum6 = 0.0;
			if (n == 4)
				sscanf(tmp, "%lf%lf%lf%lf", &fdum1, &fdum2, &fdum3, &fdum4);
			else
				sscanf(tmp, "%lf%lf%lf%lf%lf%lf", &fdum1, &fdum2, &fdum3, &fdum4, &fdum5, &fdum6);
			C[ci][1] = fdum1;
			C[ci][2] = fdum2;
			C[ci][3] = fdum3;
			C[ci][4] = fdum4;
			if (n == 6)
			{
				C[ci][5] = fdum5;
				C[ci][6] = fdum6;
			}
		}
		else
		{
			tmp = strstr(line, "sigma*sqrt(X2/n)=");
			if (tmp != NULL)
			{
				tmp += strlen("sigma*sqrt(X2/n)=");
				sscanf(tmp, "%lf", sigmaResidual);
			}
		}
		/* Make sure no infintite loop if file wrong */
		count++;
		if (count > 8)
			error("Problem reading covariance matrix");
	}
}
