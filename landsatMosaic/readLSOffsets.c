#include "stdio.h"
#include "string.h"
#include <stdlib.h>
#include "geotiff/xtiffio.h"   /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "clib/standard.h"
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"

#define MAXTIEPOINTS 500000
static void parseLSOffsetMeta(char *datFile, lsFit *fitDat, matchResult *matches);
float **LSreadFloatImage(char *file, int32 nx, int32 ny);
uint8 **LSreadByteImage(char *file, int32 nx, int32 ny);
/*
  Read tiepoint file for tiepoints.
*/

void readLSOffsets(lsFit *fitDat, matchResult *matches, int32 readData, char *maskFile)
{
	char *offsetDatFile, *offXFile, *offYFile, *typeFile;
	uint8 **mask;
	int32 i, j;
	/*
	  Parse Offset meta data
	*/
	offsetDatFile = (char *)malloc(strlen(fitDat->matchFile) + 5);
	offsetDatFile[0] = '\0';
	offsetDatFile = strcpy(offsetDatFile, fitDat->matchFile);
	offsetDatFile = strcat(offsetDatFile, ".dat");
	parseLSOffsetMeta(offsetDatFile, fitDat, matches);
	/*
	fprintf(stderr,"------------ MATCH PARAMS --------- \n");
	fprintf(stderr,"x0,y0 %lf %lf\n",matches->x0,matches->y0);
	fprintf(stderr,"dx,dy %lf %lf\n",matches->dx,matches->dy);
	fprintf(stderr,"dx,dy %lf %lf\n",matches->jdEarly,matches->jdLate);
	fprintf(stderr,"nx,ny %i %i\n",matches->nx,matches->ny);
	fprintf(stderr,"stepX,stepY %u %u\n",matches->stepX,matches->stepY);
	fprintf(stderr,"EPSG %i\n",fitDat->proj);
	fprintf(stderr,"deltaT %f\n",fitDat->deltaT);
	fprintf(stderr,"------------ END MATCH PARAMS --------- \n");*/
	/*
	  Read data X data
	*/
	if (readData == FALSE)
		return; /* This option is so file read up front in mosaicker doesn't leak memory */
	fprintf(stderr, "Read -  X offsets \n");
	offXFile = (char *)malloc(strlen(fitDat->matchFile) + 4);
	offXFile[0] = '\0';
	offXFile = strcpy(offXFile, fitDat->matchFile);
	offXFile = strcat(offXFile, ".dx");
	matches->X = LSreadFloatImage(offXFile, matches->nx, matches->ny);
	/*
	  Read data Y data
	*/
	fprintf(stderr, "Read - Y offsets \n");
	offYFile = (char *)malloc(strlen(fitDat->matchFile) + 4);
	offYFile[0] = '\0';
	offYFile = strcpy(offYFile, fitDat->matchFile);
	offYFile = strcat(offYFile, ".dy");
	matches->Y = LSreadFloatImage(offYFile, matches->nx, matches->ny);
	/*
	  Read match type data
	*/
	fprintf(stderr, "Read - Y offsets \n");
	typeFile = (char *)malloc(strlen(fitDat->matchFile) + 7);
	typeFile[0] = '\0';
	typeFile = strcpy(typeFile, fitDat->matchFile);
	typeFile = strcat(typeFile, ".mtype");
	matches->type = LSreadByteImage(typeFile, matches->nx, matches->ny);
	/*
	 Mask data
	*/
	if (maskFile != NULL)
	{
		fprintf(stderr, "masking data\n");
		mask = LSreadByteImage(maskFile, matches->nx, matches->ny);
		for (i = 0; i < matches->ny; i++)
			for (j = 0; j < matches->nx; j++)
				if (mask[i][j] == 0)
				{
					matches->Y[i][j] = NODATA;
					matches->X[i][j] = NODATA;
					matches->type[i][j] = 0; /* Not sure if this actually gets used */
				}
	}
}

/***************************LSreadFloatImage*************************************
input a floating point image file, with size nx by ny
***********************************************************************************/
float **LSreadFloatImage(char *file, int32 nx, int32 ny)
{
	float **image;
	float *tmp;
	int32 i;
	FILE *fp;
	/*
	   open file
	*/

	image = (float **)malloc(sizeof(float *) * ny);
	tmp = (float *)malloc(sizeof(float) * (size_t)nx * (size_t)ny);
	for (i = 0; i < ny; i++)
		image[i] = &(tmp[i * nx]);
	/*
	   open file
	*/
	fp = openInputFile(file);
	/*
	   read file
	*/
	for (i = 0; i < ny; i++)
	{
		freadBS(image[i], sizeof(float), (size_t)(nx), fp, FLOAT32FLAG);
	}
	fclose(fp);
	return (image);
}

uint8 **LSreadByteImage(char *file, int32 nx, int32 ny)
{
	uint8 **image;
	uint8 *tmp;
	int32 i;
	size_t ret;
	FILE *fp;
	/*
	   open file
	*/
	image = (uint8 **)malloc(sizeof(uint8 *) * ny);
	tmp = (uint8 *)malloc(sizeof(uint8) * (size_t)nx * (size_t)ny);
	for (i = 0; i < ny; i++)
		image[i] = &(tmp[i * nx]);
	/*
	   open file
	*/
	fp = openInputFile(file);
	/*
	   read file
	*/
	for (i = 0; i < ny; i++)
	{
		ret = fread(image[i], sizeof(uint8), (size_t)(nx), fp);
	}
	fclose(fp);
	return (image);
}

static void parseLSOffsetMeta(char *datFile, lsFit *fitDat, matchResult *matches)
{
	int32 notdone, idum, i;		 /* Loop flag */
	int32 linelength, lineCount; /* Input line length */
	char *line;					 /* Input line buffer */
	char *keyword, *value;
	uint32 udum;
	FILE *fp;
	double fdum;
	keyword = (char *)malloc(sizeof(char) * LINEMAX + 1);
	value = (char *)malloc(sizeof(char) * LINEMAX + 1);
	fprintf(stderr, "%s\n", datFile);
	line = (char *)malloc(sizeof(char) * (LINEMAX + 1));
	/*
	  Open file
	*/
	fp = openInputFile(datFile);
	/*
	  Loop through lines
	*/
	notdone = TRUE;
	matches->x0 = NODATA;
	matches->y0 = NODATA;
	matches->dx = NODATA;
	matches->dy = NODATA;
	matches->nx = 0;
	matches->ny = 0;
	fitDat->proj = -1;
	fitDat->deltaT = NODATA;

	while (notdone == TRUE)
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
				if (strstr(keyword, "fileEarly") != NULL)
				{
					for (i = 0; i < strlen(value); i++)
						if (value[i] == '\n' || value[i] == '\r')
							value[i] = '\0'; /* strip carriage return */
					matches->fileEarly = (char *)malloc((strlen(value) + 1) * sizeof(char));
					strcpy(matches->fileEarly, value);
				}
				if (strstr(keyword, "fileLate") != NULL)
				{
					for (i = 0; i < strlen(value); i++)
						if (value[i] == '\n' || value[i] == '\r')
							value[i] = '\0'; /* strip carriage return */
					matches->fileLate = (char *)malloc((strlen(value) + 1) * sizeof(char));
					strcpy(matches->fileLate, value);
				}
				if (strstr(keyword, "x0") != NULL)
				{
					if (sscanf(value, "%lf\n", &fdum) != 1)
						error("parse LS Offsets x0");
					matches->x0 = fdum;
				}
				if (strstr(keyword, "y0") != NULL)
				{
					if (sscanf(value, "%lf\n", &fdum) != 1)
						error("parse LS Offsets y0");
					matches->y0 = fdum;
				}
				if (strstr(keyword, "dx") != NULL)
				{
					if (sscanf(value, "%lf\n", &fdum) != 1)
						error("parse LS Offsets dx");
					matches->dx = fdum;
				}
				if (strstr(keyword, "dy") != NULL)
				{
					if (sscanf(value, "%lf\n", &fdum) != 1)
						error("parse LS Offsets dy");
					matches->dy = fdum;
				}
				if (strstr(keyword, "earlyImageJD") != NULL)
				{
					if (sscanf(value, "%lf\n", &fdum) != 1)
						error("parse LS Offsets earlyImageJD");
					matches->jdEarly = fdum;
				}
				if (strstr(keyword, "lateImageJD") != NULL)
				{
					if (sscanf(value, "%lf\n", &fdum) != 1)
						error("parse LS Offsets lateImageJD");
					matches->jdLate = fdum;
				}

				if (strstr(keyword, "nx") != NULL)
				{
					if (sscanf(value, "%i\n", &idum) != 1)
						error("parse LS Offsets nx");
					matches->nx = idum;
				}
				if (strstr(keyword, "ny") != NULL)
				{
					if (sscanf(value, "%i\n", &idum) != 1)
						error("parse LS Offsets nx");
					matches->ny = idum;
				}
				if (strstr(keyword, "stepX") != NULL)
				{
					if (sscanf(value, "%u\n", &udum) != 1)
						error("parse LS Offsets nx");
					matches->stepX = udum;
				}
				if (strstr(keyword, "stepY") != NULL)
				{
					if (sscanf(value, "%u\n", &udum) != 1)
						error("parse LS Offsets nx");
					matches->stepY = udum;
				}
				if (strstr(keyword, "EPSG") != NULL)
				{
					if (sscanf(value, "%i\n", &idum) != 1)
						error("parse LS Offsets proj: %s", value);
					fitDat->proj = idum;
				}
				if (strstr(keyword, "slowFlag") != NULL)
				{
					if (sscanf(value, "%i\n", &idum) != 1)
						error("parse LS Offsets slowFlag");
					fitDat->slowFlag = idum;
				}
				if (strstr(keyword, "IntervalBetweenImages") != NULL)
				{
					if (sscanf(value, "%lf\n", &fdum) != 1)
						error("parse LS Offsets deltaT");
					fitDat->deltaT = fdum;
				}

				if (strstr(keyword, "Mean_sigmaX") != NULL)
				{
					if (sscanf(value, "%lf\n", &fdum) != 1)
					{
						matches->meanSigmaX = NODATA;
					}
					matches->meanSigmaX = fdum;
				}
				if (strstr(keyword, "Mean_sigmaY") != NULL)
				{
					if (sscanf(value, "%lf\n", &fdum) != 1)
					{
						matches->meanSigmaY = NODATA;
					}
					matches->meanSigmaY = fdum;
				}
			}
		} /* End else */
	}	  /* End while */
	if (matches->x0 < (NODATA + 1))
		error("ParseLSOffsetMeta: Invalid x0\n");
	if (matches->y0 < (NODATA + 1))
		error("ParseLSOffsetMeta: Invalid y0\n");
	if (matches->dx < (NODATA + 1))
		error("ParseLSOffsetMeta: Invalid dx\n");
	if (matches->dy < (NODATA + 1))
		error("ParseLSOffsetMeta: Invalid dy\n");
	if (matches->nx <= 0)
		error("ParseLSOffsetMeta: Invalid x0\n");
	if (matches->nx <= 0)
		error("ParseLSOffsetMeta: Invalid y0\n");
	if (fitDat->deltaT < (NODATA + 1))
		error("ParseLSOffsetMeta: Invalid deltaT\n");
	if (fitDat->proj <= 0)
		error("ParseLSOffsetMeta: Invalid proj %i\n", fitDat->proj);
	fclose(fp);
}

void parseKeyValue(char *line, char *keyword, char *value)
{
	char tmp1[LINEMAX];
	uint32 i, j;
	/*
	  Process keyword
	*/
	i = 0;
	while (line[i] != '=' && i < LINEMAX)
	{
		tmp1[i] = line[i];
		i++;
	};
	if (i == LINEMAX)
		error("\nparseKeyValue: Invalid keyword pair\n");
	tmp1[i] = '\0';
	strcpy(keyword, tmp1);
	/*
	  Process value
	*/
	j = 0;
	i++;
	while (line[i] != '\0' && i < LINEMAX)
	{
		tmp1[j] = line[i];
		i++;
		j++;
	};
	if (i == LINEMAX)
		error("\nparseKeyValue: Invalid keyword pair\n");
	tmp1[i] = '\0';
	strcpy(value, tmp1);
}

/*
  Strip leading and trailing spaces and tabs from a string
*/

char *stripWS(char *string)
{
	char *tmp1;
	while (*string == ' ' || *string == '\t' || *string == '\n')
		string++;
	tmp1 = string + strlen(string) - 1;

	while (*tmp1 == ' ' || *tmp1 == '\t' || *tmp1 == '\n')
	{
		/*fprintf(stderr,"=%s\n",tmp1);*/
		*tmp1 = '\0';
		tmp1--;
	}
	return (string);
}
