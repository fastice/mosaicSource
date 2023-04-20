#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "mosaicSource/common/common.h"
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
#include "mosaicSource/landsatMosaic/landSatMosaic.h"
static void parseFitFile(lsFit *fitDat);
char *lineBufLS, *keyBuf, *valueBuf; /* Input line buffer */
/*
  Parse a file with list of lansat images for inclusion in mosaic3d.
*/
float dateWeight(double jd1, double jd2, double jdRange1, double jdRange2);

landSatImage *parseLSInputs(char *inputFile, landSatImage *LSImages, double jd1, double jd2, int timeOverlapFlag)
{
	extern char *lineBufLS, *keyBuf, *valueBuf; /* Input line buffer */
	FILE *fp;
	int32_t notdone, lineCount, lineLength;
	landSatImage *tmpImage, *firstImage, *prevImage;
	char *line, *lineSave;
	char *matchFile, *fitFile, *maskFile;
	float tmpWeight;
	lineBufLS = NULL;
	keyBuf = NULL;
	valueBuf = NULL;
	keyBuf = (char *)malloc(sizeof(char) * (LINEMAX + 1));
	valueBuf = (char *)malloc(sizeof(char) * (LINEMAX + 1));
	lineBufLS = (char *)malloc(sizeof(char) * (LINEMAX + 1));

	line = (char *)malloc(sizeof(char) * (LINEMAX + 1));
	lineSave = line;
	fprintf(stderr, "Enter ParseLSInputs -- %s \n", inputFile);
	/*
	   Open Input File
	*/
	fp = openInputFile(inputFile);
	/*
	  Loop through lines
	*/
	notdone = TRUE;
	firstImage = NULL;
	while (notdone == TRUE)
	{ /* Loop to read lines */
		line = lineSave;
		lineLength = (int32_t)fgetline(fp, line, LINEMAX); /* Read line */
		lineCount++;
		if (strchr(line, ENDDATA) != NULL || lineLength == 0)
		{
			notdone = FALSE; /* End of data, set exit flag */
		}
		else if (strchr(line, COMMENT) == NULL)
		{						  /* If not comment, parse */
			line = stripWS(line); /* Strip any white space */
			matchFile = strtok(line, " ");
			fitFile = strtok(NULL, " ");
			maskFile = strtok(NULL, " ");
			/*fprintf(stderr,"%s \n %s\n %s\n",matchFile,fitFile,maskFile);*/
			/*
			   Malloc image structure
			*/
			tmpImage = (landSatImage *)malloc(sizeof(landSatImage));
			tmpImage->next = NULL;
			if (LSImages == NULL)
			{
				firstImage = tmpImage;
				prevImage = NULL;
				LSImages = tmpImage;
				firstImage->nImages = 0;
			}
			else
			{
				LSImages->next = tmpImage;
				prevImage = LSImages;
				LSImages = tmpImage;
			}
			/*
			  Parse mask file
			*/
			if (maskFile != NULL)
			{
				LSImages->maskFile = (char *)malloc(sizeof(char) * (strlen(maskFile) + 1));
				LSImages->maskFile[0] = '\0';
				strcpy(LSImages->maskFile, maskFile);
			}
			else
			{
				LSImages->maskFile = NULL;
			}
			/*
			  Parse match data
			*/
			LSImages->fitResult.matchFile = (char *)malloc(sizeof(char) * (strlen(matchFile) + 1));
			LSImages->fitResult.matchFile[0] = '\0';
			strcpy(LSImages->fitResult.matchFile, matchFile);
			readLSOffsets(&(LSImages->fitResult), &(LSImages->matches), FALSE, maskFile); /* False forces skip reading data */
			/*
			  Check if in allowed time range
			*/
			tmpWeight = dateWeight(LSImages->matches.jdEarly, LSImages->matches.jdLate, jd1, jd2);
			LSImages->weight = 0.0;
			/* Case 1 overlapflag - accept anthing with better than 5% overlap, or strictly in range > 99% overlap */
			/*fprintf(stderr,"tmpW %f %i\n",tmpWeight,timeOverlapFlag);	*/
			if ((timeOverlapFlag == TRUE && tmpWeight > 0.05) || (timeOverlapFlag == FALSE && tmpWeight > 0.99))
			{
				/* Image is in range so populate with fit and other data */
				firstImage->nImages++;
				LSImages->weight = tmpWeight;
				/*
				  Parse fit data
				*/
				if (fitFile == NULL)
				{
					LSImages->fitResult.fitFile = (char *)malloc(sizeof(char) * (strlen(line) + 1 + 4));
					LSImages->fitResult.fitFile[0] = '\0';
					strcpy(LSImages->fitResult.fitFile, matchFile);
					strcat(LSImages->fitResult.fitFile, ".fit");
				}
				else
				{
					LSImages->fitResult.fitFile = (char *)malloc(sizeof(char) * (strlen(fitFile) + 1));
					LSImages->fitResult.fitFile[0] = '\0';
					strcpy(LSImages->fitResult.fitFile, fitFile);
				}
				parseFitFile(&(LSImages->fitResult));
			}
			else
			{
				/* Image not used (not in date range) so back up, effectively discarding this image */
				LSImages = prevImage;
			}
		}
	}
	fclose(fp);
	if (firstImage != NULL)
		fprintf(stderr, "Exit ParseLSInputs %i\n", firstImage->nImages);
	else
		fprintf(stderr, "Exit ParseLSInputs NULL\n");
	free(line);
	return (firstImage);
}

static void parseFitFile(lsFit *fitDat)
{
	int32_t notdone, i, j;						/* Loop flag */
	int32_t linelength, lineCount;				/* Input line length */
	extern char *lineBufLS, *keyBuf, *valueBuf; /* Input line buffer */
	char *line, *keyword, *value;
	FILE *fp;
	double fdum1, fdum2, fdum3;

	keyword = (char *)keyBuf;
	value = (char *)valueBuf;
	line = (char *)lineBufLS;
	/*
	  Open file
	*/
	fp = openInputFile(fitDat->fitFile);
	/*
	  Loop through lines
	*/
	notdone = TRUE;
	lineCount = 0;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fitDat->Cx[i][j] = LARGEINT;
			fitDat->Cy[i][j] = LARGEINT;
		}
	}
	for (i = 0; i < MAXFITPARAM; i++)
	{
		fitDat->pX[i] = NODATA;
		fitDat->pY[i] = NODATA;
	}
	while (notdone == TRUE && lineCount < 30)
	{											  /* Loop to read lines */
		linelength = fgetline(fp, line, LINEMAX); /* Read line */
		lineCount++;
		if (strchr(line, ENDDATA) != NULL)
			notdone = FALSE; /* End of data, set exit flag */
		else if (strchr(line, COMMENT) == NULL)
		{ /* If not comment, parse */
			if (strchr(line, '=') != NULL)
			{
				parseKeyValue(line, keyword, value);
				if (strstr(keyword, "Xfit") != NULL)
				{
					if (sscanf(value, "%lf%lf%lf\n", &fdum1, &fdum2, &fdum3) != 3)
						error("parse fit X fit");
					fitDat->pX[0] = fdum1;
					fitDat->pX[1] = fdum2;
					fitDat->pX[2] = fdum3;
				}
				if (strstr(keyword, "Yfit") != NULL)
				{
					if (sscanf(value, "%lf%lf%lf\n", &fdum1, &fdum2, &fdum3) != 3)
						error("parse fit Y fit");
					fitDat->pY[0] = fdum1;
					fitDat->pY[1] = fdum2;
					fitDat->pY[2] = fdum3;
				}
				/*
				  Read residuals in units of meters
				*/
				if (strstr(keyword, "X_residual") != NULL)
				{
					if (sscanf(value, "%lf%lf%lf\n", &fdum1, &fdum2, &fdum3) != 3)
					{
						fitDat->sigmaXRes = NODATA;
					}
					else
						fitDat->sigmaXRes = fdum3;
				}
				if (strstr(keyword, "Y_residual") != NULL)
				{
					if (sscanf(value, "%lf%lf%lf\n", &fdum1, &fdum2, &fdum3) != 3)
					{
						fitDat->sigmaYRes = NODATA;
					}
					else
						fitDat->sigmaYRes = fdum3;
				}

				if (strstr(keyword, "CX_11_12_13") != NULL)
				{
					if (sscanf(value, "%lf%lf%lf\n", &fdum1, &fdum2, &fdum3) == 3)
					{
						fitDat->Cx[0][0] = fdum1;
						fitDat->Cx[0][1] = fdum2;
						fitDat->Cx[0][2] = fdum3;
					}
				}

				if (strstr(keyword, "CX_21_22_23") != NULL)
				{
					if (sscanf(value, "%lf%lf%lf\n", &fdum1, &fdum2, &fdum3) == 3)
					{
						fitDat->Cx[1][0] = fdum1;
						fitDat->Cx[1][1] = fdum2;
						fitDat->Cx[1][2] = fdum3;
					}
				}

				if (strstr(keyword, "CX_31_32_33") != NULL)
				{
					if (sscanf(value, "%lf%lf%lf\n", &fdum1, &fdum2, &fdum3) == 3)
					{
						fitDat->Cx[2][0] = fdum1;
						fitDat->Cx[2][1] = fdum2;
						fitDat->Cx[2][2] = fdum3;
					}
				}
				if (strstr(keyword, "CY_11_12_13") != NULL)
				{
					if (sscanf(value, "%lf%lf%lf\n", &fdum1, &fdum2, &fdum3) == 3)
					{
						fitDat->Cy[0][0] = fdum1;
						fitDat->Cy[0][1] = fdum2;
						fitDat->Cy[0][2] = fdum3;
					}
				}

				if (strstr(keyword, "CY_21_22_23") != NULL)
				{
					if (sscanf(value, "%lf%lf%lf\n", &fdum1, &fdum2, &fdum3) == 3)
					{
						fitDat->Cy[1][0] = fdum1;
						fitDat->Cy[1][1] = fdum2;
						fitDat->Cy[1][2] = fdum3;
					}
				}

				if (strstr(keyword, "CY_31_32_33") != NULL)
				{
					if (sscanf(value, "%lf%lf%lf\n", &fdum1, &fdum2, &fdum3) == 3)
					{
						fitDat->Cy[2][0] = fdum1;
						fitDat->Cy[2][1] = fdum2;
						fitDat->Cy[2][2] = fdum3;
					}
				}
			} /* end if line=... */
		}	  /* E nd else */
	}		  /* End while */
	if (lineCount >= 30)
		error("> 30 lines in fit file, likely error %s \n", fitDat->fitFile);

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			/*fprintf(stderr,"%f %f \n",fitDat->Cx[i][j],fitDat->Cy[i][j]);*/
			if (fitDat->Cx[i][j] == LARGEINT || fitDat->Cy[i][j] == LARGEINT)
				error("-- Missing Covarance Elements--");
		}
	}

	for (i = 0; i < MAXFITPARAM; i++)
	{
		if (fitDat->pX[i] == NODATA || fitDat->pY[i] == NODATA)
			error("Error reading fit");
	}
	/*fprintf(stderr,"Xfit_const_coeffX_coeff %le %le %le\n",fitDat->pX[0],fitDat->pX[1],fitDat->pX[2]);
	  fprintf(stderr,"Yfit_const_coeffX_coeff %le %le %le\n",fitDat->pY[0],fitDat->pY[1],fitDat->pY[2]);*/
	/*	free(keyword);free(value);	free(line);*/
	fclose(fp);
}
