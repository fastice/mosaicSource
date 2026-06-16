#include "stdio.h"
#include "string.h"
#include "common.h"
#include "math.h"

/*
  Input baseline info.
*/
void getBaseline(char *baselineFile, vhParams *params, int32_t noPhase)
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
			params->Bn = 0.0; params->Bp = 0.0;
			params->dBn = 0.0; params->dBp = 0.0;
			params->dBnQ = 0.0; params->dBpQ = 0.0;
			int inC = 0, ci = 0;
			while (fgets(line, 256, yfp))
			{
				if      (sscanf(line, "nDays: %lf",  &params->nDays)  == 1) { inC = 0; }
				else if (sscanf(line, "sigma: %lf",  &params->sigma)  == 1) { inC = 0; }
				else if (sscanf(line, "Bn: %lf",     &params->Bn)     == 1) { inC = 0; }
				else if (sscanf(line, "Bp: %lf",     &params->Bp)     == 1) { inC = 0; }
				else if (sscanf(line, "dBn: %lf",    &params->dBn)    == 1) { inC = 0; }
				else if (sscanf(line, "dBp: %lf",    &params->dBp)    == 1) { inC = 0; }
				else if (sscanf(line, "dBnQ: %lf",   &params->dBnQ)   == 1) { inC = 0; }
				else if (sscanf(line, "dBpQ: %lf",   &params->dBpQ)   == 1) { inC = 0; }
				else if (strncmp(line, "applyFlatEarth: true", 20) == 0)
					{ params->applyFlatEarth = 1; inC = 0; }
				else if (strncmp(line, "C:", 2) == 0)
					{ inC = 1; ci = 0; }
				else if (inC && strstr(line, "- [") != NULL && ci < 6)
				{
					char *p = strstr(line, "[");
					if (p)
						sscanf(p + 1, "%lf, %lf, %lf, %lf, %lf, %lf",
							   &params->C[ci+1][1], &params->C[ci+1][2],
							   &params->C[ci+1][3], &params->C[ci+1][4],
							   &params->C[ci+1][5], &params->C[ci+1][6]);
					ci++;
				}
				else { inC = 0; }
			}
			fclose(yfp);
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
	fclose(fp);
}
