#include "stdio.h"
#include "string.h"
#include "common.h"
#include "math.h"
#include <libgen.h>
#include <unistd.h>

/* One labeled solution block parsed from a tiepoints -yaml baseline file (noSquint:/squint:) */
typedef struct
{
	double Bn, Bp, dBn, dBp, dBnQ, dBpQ, sigma;
	double C[7][7];
	int32_t seen;
} yamlBaselineBlock;

/*
  Store the ionospheric phase file named by a tiepoints baseline file into
  params->ionospherePhaseFile. A relative name is resolved against the baseline file's own
  directory, so the two travel together; an absolute path is used as given. "nil" (the
  no-correction marker tiepoints writes) leaves the field empty.

  The check that the file exists happens here rather than at read time so a mis-plumbed
  path fails immediately with a clear message, matching checkForIonosphereCorrection() on
  the range-offset side (readOffsets.c).
*/
static void setIonospherePhaseFile(vhParams *params, char *ionName, char *baselineFile)
{
	char tmp[2048];

	while (*ionName == ' ' || *ionName == '\t')
	{
		ionName++;
	}
	int32_t len = strlen(ionName);
	while (len > 0 && (ionName[len - 1] == '\n' || ionName[len - 1] == '\r' || ionName[len - 1] == ' '))
	{
		ionName[--len] = '\0';
	}
	if (len == 0 || strncmp(ionName, "nil", 3) == 0)
	{
		params->ionospherePhaseFile[0] = '\0';
		return;
	}
	if (ionName[0] == '/')
	{
		strncpy(params->ionospherePhaseFile, ionName, sizeof(params->ionospherePhaseFile) - 1);
	}
	else
	{
		strncpy(tmp, baselineFile, sizeof(tmp) - 1);
		tmp[sizeof(tmp) - 1] = '\0';
		snprintf(params->ionospherePhaseFile, sizeof(params->ionospherePhaseFile),
				 "%s/%s", dirname(tmp), ionName);
	}
	params->ionospherePhaseFile[sizeof(params->ionospherePhaseFile) - 1] = '\0';
	if (access(params->ionospherePhaseFile, F_OK) != 0)
		error("getBaseline: %s names ionosphere correction %s, which does not exist\n",
			  baselineFile, params->ionospherePhaseFile);
	fprintf(stderr, "getBaseline: found ionosphereCorrectionFile %s\n", params->ionospherePhaseFile);
}

/*
  Input baseline info. useSquint selects which labeled YAML block (noSquint/squint) to use
  when the baseline file has both -- see mosaicSource/CLAUDE.md "Squint". Ignored for
  non-YAML baseline files, which predate the squint feature and carry only one solution.
*/
void getBaseline(char *baselineFile, vhParams *params, int32_t noPhase, int32_t useSquint)
{
	FILE *fp;
	double BnEst, BpEst, dBnEst, dBpEst, dBpEstQ, dBnEstQ;
	int32_t nBaselines;
	char *tmp;
	double dum;
	double fdum1, fdum2, fdum3, fdum4, fdum5, fdum6;
	int32_t lineCount = 0, eod, special;
	int32_t i, j;
	int32_t ci;
	char line[256];
	/*
	  Input parm info
	*/
	for (i = 1; i <= 6; i++)
		for (j = 1; j <= 6; j++)
			params->C[i][j] = 0;
	params->sigma = PI; /* Default value */
	params->applyFlatEarth = 0;
	params->ionospherePhaseFile[0] = '\0';
	/*
	  YAML baseline file (extension .yaml): flat-earth correction, no topo term.
	*/
	{
		size_t blen = strlen(baselineFile);
		if (blen > 5 && strcmp(baselineFile + blen - 5, ".yaml") == 0)
		{
			FILE *yfp = fopen(baselineFile, "r");
			if (yfp == NULL)
				error("getBaseline: cannot open YAML file %s", baselineFile);
			yamlBaselineBlock noSquint, squint;
			int32_t hasSquintSolution = 0;
			memset(&noSquint, 0, sizeof(noSquint));
			memset(&squint, 0, sizeof(squint));
			/* section: 0 = top-level (applyFlatEarth/hasSquintSolution/nDays), 1 = noSquint:, 2 = squint: */
			int32_t section = 0, inC = 0, ci = 0;
			yamlBaselineBlock *blk = NULL;
			while (fgets(line, 256, yfp))
			{
				/* Top-level (0-indent) keys, checked before whitespace-stripping so they
				   can't be confused with an indented child key of the same name. */
				if (line[0] != ' ' && line[0] != '\t')
				{
					if (strncmp(line, "noSquint:", 9) == 0)
						{ section = 1; blk = &noSquint; inC = 0; continue; }
					else if (strncmp(line, "squint:", 7) == 0)
						{ section = 2; blk = &squint; inC = 0; continue; }
					else if (strncmp(line, "applyFlatEarth: true", 20) == 0)
						{ params->applyFlatEarth = 1; continue; }
					else if (strncmp(line, "hasSquintSolution: true", 23) == 0)
						{ hasSquintSolution = 1; continue; }
					else if (sscanf(line, "nDays: %lf", &params->nDays) == 1)
						{ continue; }
					else if (strncmp(line, "ionosphereCorrectionFile:", 25) == 0)
						{ setIonospherePhaseFile(params, line + 25, baselineFile); continue; }
					/* Legacy flat (pre-squint, single-solution) files have no noSquint:/
					   squint: headers at all -- their Bn:/Bp:/etc. keys sit at 0-indent.
					   Route those into the noSquint block directly. */
					section = 1; blk = &noSquint;
				}
				char *p = line;
				while (*p == ' ' || *p == '\t') p++;
				if      (sscanf(p, "sigma: %lf", &blk->sigma) == 1) { inC = 0; }
				else if (sscanf(p, "Bn: %lf",     &blk->Bn)     == 1) { inC = 0; }
				else if (sscanf(p, "Bp: %lf",     &blk->Bp)     == 1) { inC = 0; }
				else if (sscanf(p, "dBn: %lf",    &blk->dBn)    == 1) { inC = 0; }
				else if (sscanf(p, "dBp: %lf",    &blk->dBp)    == 1) { inC = 0; }
				else if (sscanf(p, "dBnQ: %lf",   &blk->dBnQ)   == 1) { inC = 0; }
				else if (sscanf(p, "dBpQ: %lf",   &blk->dBpQ)   == 1) { inC = 0; }
				else if (strncmp(p, "C:", 2) == 0)
					{ inC = 1; ci = 0; blk->seen = 1; }
				else if (inC && strstr(p, "- [") != NULL && ci < 6)
				{
					char *pb = strstr(p, "[");
					if (pb)
						sscanf(pb + 1, "%lf, %lf, %lf, %lf, %lf, %lf",
							   &blk->C[ci+1][1], &blk->C[ci+1][2],
							   &blk->C[ci+1][3], &blk->C[ci+1][4],
							   &blk->C[ci+1][5], &blk->C[ci+1][6]);
					ci++;
				}
				else { inC = 0; }
			}
			fclose(yfp);

			yamlBaselineBlock *sel = (useSquint && hasSquintSolution) ? &squint : &noSquint;
			params->Bn = sel->Bn;
			params->Bp = sel->Bp;
			params->dBn = sel->dBn;
			params->dBp = sel->dBp;
			params->dBnQ = sel->dBnQ;
			params->dBpQ = sel->dBpQ;
			params->sigma = sel->sigma;
			for (i = 1; i <= 6; i++)
				for (j = 1; j <= 6; j++)
					params->C[i][j] = sel->C[i][j];
			return;
		}
	}
	if (strstr(baselineFile, "nobaseline") != NULL || noPhase == TRUE)
	{
		/* No baseline used for this solution so return */
		params->Bn = 0;
		params->Bp = 0;
		params->dBn = 0;
		params->dBp = 0;
		params->dBnQ = 0;
		params->dBpQ = 0;
		return;
	}
	fp = openInputFile(baselineFile);
	/*
	  Skip past initial data lines
	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%i", &nBaselines) != 1)
		error("getBaseline -- Missing baseline params at line: %i of %s", lineCount, baselineFile);
	if (nBaselines < 1 || nBaselines > 3)
		error("getBaseline -- invalid number of baselines at line %i of %s\n", lineCount, baselineFile);
	if (nBaselines > 1)
		lineCount = getDataString(fp, lineCount, line, &eod);
	if (nBaselines > 2)
		lineCount = getDataString(fp, lineCount, line, &eod);
	/*
	  Read covariance matrix if there is one.
	*/
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
			sscanf(tmp, "%lf%lf%lf%lf%lf%lf", &fdum1, &fdum2, &fdum3, &fdum4, &fdum5, &fdum6);
			params->C[ci][1] = fdum1;
			params->C[ci][2] = fdum2;
			params->C[ci][3] = fdum3;
			params->C[ci][4] = fdum4;
			params->C[ci][5] = fdum5;
			params->C[ci][6] = fdum6;
		}
		else
		{
			tmp = strstr(line, "sigma*sqrt(X2/n)=");
			if (tmp != NULL)
			{
				tmp += strlen("sigma*sqrt(X2/n)=");
				sscanf(tmp, "%lf", &fdum1);
				params->sigma = fdum1;
			}
		}
	}
	//fprintf(stderr, "sigma*sqrt(X2/n) = %lf\n", params->sigma);
	/*
	  Input baseline estimated with tiepoints.
	*/
	if (sscanf(line, "%lf%lf%lf%lf%lf%lf", &dum, &dum, &dum, &dum, &dum, &dum) == 6)
	{
		fprintf(stderr, "\n***BASELINE FIT WAS QUADRATIC***\n");
		sscanf(line, "%lf%lf%lf%lf%lf%lf", &BnEst, &BpEst, &dBnEst, &dBpEst, &dBnEstQ, &dBpEstQ);
	}
	else
	{
		if (sscanf(line, "%lf%lf%lf%lf", &BnEst, &BpEst, &dBnEst, &dBpEst) != 4)
			error("%s %i", "getBaseline  -- Missing baseline params at line:", lineCount);
		dBnEstQ = 0.0;
		dBpEstQ = 0.0;
	}
	/* populate baseline params structure */
	params->Bn = BnEst;
	params->Bp = BpEst;
	params->dBn = dBnEst;
	params->dBp = dBpEst;
	params->dBnQ = dBnEstQ;
	params->dBpQ = dBnEstQ;
	/*
	  Scan the remaining lines (including those after '&') for a ;* ionosphereCorrectionFile
	  entry, the legacy-text counterpart of the YAML key above. Same scan pattern getRParams()
	  uses for ;* offsetCorrectionFile (readOffsets.c).
	*/
	{
		char scLine[256], *tmp2;
		while (fgets(scLine, sizeof(scLine), fp) != NULL)
		{
			/* Must start with ';' and contain '*' to be a special line */
			if (scLine[0] != ';')
			{
				continue;
			}
			if (strchr(scLine, '*') == NULL)
			{
				continue;
			}
			tmp2 = strstr(scLine, "ionosphereCorrectionFile");
			if (tmp2 != NULL)
			{
				setIonospherePhaseFile(params, tmp2 + strlen("ionosphereCorrectionFile"), baselineFile);
				break;
			}
		}
	}
	fclose(fp);
}
