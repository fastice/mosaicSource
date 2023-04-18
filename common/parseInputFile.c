#include "stdio.h"
#include "ctype.h"
#include <string.h>
#include <stdlib.h>
#include "math.h"
#include "common.h"

/*
   Parse a geodatXxY.in file - parseInputFile - with numerous supporting routines to parse data
 */

/* Parse data string and populate image date fields */
static void parseImageDate(inputImageStructure *inputImage, char *line, char *inputFile)
{
	char *tmp, *tmp1;
	int32_t m, done;
	int32_t mm, dd, yyyy;
	fprintf(stderr, "%s %s", line, inputFile);
	if ((tmp = strchr(line, ':')) == NULL)
	{
		error("--invalid date in %s %s\n", inputFile, tmp);
	}
	else
		tmp++;
	m = 0;
	done = 0;
	while (done < 1)
	{
		if (tmp[m] == '1' || tmp[m] == '2' || tmp[m] == '3' || tmp[m] == '4' || tmp[m] == '5' ||
			tmp[m] == '6' || tmp[m] == '7' || tmp[m] == '8' || tmp[m] == '9' || m >= strlen(tmp))
			done = 1;
		else if (tmp[m] == '0')
			tmp[m] = ' ';
		m++;
	}
	/*fprintf(stderr,"%s\n",tmp);*/
	sscanf(tmp, "%i\n", &dd);
	if (tmp1 = strstr(tmp, "JAN"), tmp1 != NULL)
		mm = 1;
	else if (tmp1 = strstr(tmp, "FEB"), tmp1 != NULL)
		mm = 2;
	else if (tmp1 = strstr(tmp, "MAR"), tmp1 != NULL)
		mm = 3;
	else if (tmp1 = strstr(tmp, "APR"), tmp1 != NULL)
		mm = 4;
	else if (tmp1 = strstr(tmp, "MAY"), tmp1 != NULL)
		mm = 5;
	else if (tmp1 = strstr(tmp, "JUN"), tmp1 != NULL)
		mm = 6;
	else if (tmp1 = strstr(tmp, "JUL"), tmp1 != NULL)
		mm = 7;
	else if (tmp1 = strstr(tmp, "AUG"), tmp1 != NULL)
		mm = 8;
	else if (tmp1 = strstr(tmp, "SEP"), tmp1 != NULL)
		mm = 9;
	else if (tmp1 = strstr(tmp, "OCT"), tmp1 != NULL)
		mm = 10;
	else if (tmp1 = strstr(tmp, "NOV"), tmp1 != NULL)
		mm = 11;
	else if (tmp1 = strstr(tmp, "DEC"), tmp1 != NULL)
		mm = 12;
	else
		error("++invalid date in %s %s\n", inputFile, tmp);
	sscanf((tmp1 + 3), "%i\n", &yyyy);
	inputImage->year = yyyy;
	inputImage->month = mm;
	inputImage->day = dd;
	/* fprintf(stderr,"mm,dd,yyyy %i %i %i\n",mm,dd,yyyy);*/
}

static void parseTime(inputImageStructure *inputImage, char *line, int32_t eod)
{
	if (sscanf(line, "%i %i %lf\n", &(inputImage->par.hr), &(inputImage->par.min),
			   &(inputImage->par.sec)) != 3)
		error("Invalid image start time \n %s\n", line);
	if (eod == TRUE)
		error("missing data after time");
}

static void parsePRF(inputImageStructure *inputImage, char *line, int32_t eod)
{
	if (sscanf(line, "%lf\n", &(inputImage->par.prf)) != 1)
		error("Invalid prf");
	if (eod == TRUE)
		error("missing data after prf");
}

static void parseLambda(inputImageStructure *inputImage, char *line, int32_t eod)
{
	if (sscanf(line, "%lf\n", &(inputImage->par.lambda)) != 1)
		error("Invalid wavelength");
	if (eod == TRUE)
		error("missing data after lamda");
}

static int32_t parseStateVectors(FILE *fp, int32_t lineCount, inputImageStructure *inputImage, char *inputFile)
{
	stateV *sv;
	char line[512];
	int32_t eod, i;

	sv = &(inputImage->sv);
	/*
	  If state flag set read state vector information and initialize
	  par structure
	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%i", &(sv->nState)) != 1)
		error("missing number of state vectors");
	if (sv->nState < MNST || sv->nState > MXST)
		error("\nInvalid number of state vectors - %i\n", sv->nState);

	lineCount = getDataString(fp, lineCount, line, &eod); /* Get next line */
	if (eod == TRUE)
		error("missing state vector info at line %i\n", lineCount);
	sscanf(line, "%lf", &(sv->t0));

	lineCount = getDataString(fp, lineCount, line, &eod); /* Get next line */
	if (eod == TRUE)
		error("missing state vector info at line %i\n", lineCount);
	sscanf(line, "%lf", &(sv->deltaT));
	/*
	  Read state vector info - NOTE USING 1-nState indexing for polint
	*/
	for (i = 1; i <= sv->nState; i++)
	{
		/* Position */
		lineCount = getDataString(fp, lineCount, line, &eod); /* Get next line */
		if (sscanf(line, "%le %le%le", &(sv->x[i]), &(sv->y[i]), &(sv->z[i])) != 3)
			error("Invalid state x,y,z %i %s\n", lineCount, inputFile);
		/* Velocity */
		lineCount = getDataString(fp, lineCount, line, &eod); /* Get next line */
		if (sscanf(line, "%le %le%le", &(sv->vx[i]), &(sv->vy[i]), &(sv->vz[i])) != 3)
			error("Invalid vx,vy,vz at line %i %s\n", lineCount, inputFile);
		sv->times[i] = sv->t0 + (i - 1) * sv->deltaT;
	}
	return (lineCount);
}

static void parseDeltaT(inputImageStructure *inputImage, int32_t julDayInt, char *line)
{
	SARData *par;
	double deltaTCorrect;

	par = &(inputImage->par);
	deltaTCorrect = 0.0;
	deltaTCorrect += (inputImage->nAzimuthLooks - 1.0) * 0.5 / inputImage->par.prf;
	if (strstr(line, "deltaT") != NULL)
		sscanf((line + 6), "%lf\n", &deltaTCorrect);
	par->sec += deltaTCorrect;
	if (par->sec > 60.0)
	{
		par->sec -= 60.0;
		par->min += 1;
	}
	if (par->sec < .0)
	{
		par->sec += 60;
		par->min -= 1;
	}
	if (par->min > 59)
	{
		par->min -= 60;
		par->hr += 1;
	}
	if (par->min < 0)
	{
		par->min += 60;
		par->hr -= 1;
	}
	/*fprintf(stderr,"***DeltaTCorrect* %f %i %i %f\n",deltaTCorrect,par->hr,par->min,par->sec);    */
	inputImage->julDay = (double)julDayInt + (par->hr * 3600. + par->min * 60. + par->sec) / 86400.;
	/*fprintf(stderr,"julday %li %lf\n",julDayInt,inputImage->julDay);*/
}

static void parseSize(inputImageStructure *inputImage, char *line, char *inputFile)
{
	int32_t nr, na, nlr, nla;

	if (sscanf(line, "%i%i%i%i", &nr, &na, &nlr, &nla) != 4)
		error("%s  %i %s %s", "parseInputFile -- Missing image size parameters", "File: ", inputFile);
	inputImage->rangeSize = nr;
	inputImage->azimuthSize = na;
	inputImage->nRangeLooks = nlr;
	inputImage->nAzimuthLooks = nla;
	/*fprintf(stderr,"%i %i %i %i %i\n",nr,na,nlr,nla,sizeof(nla));*/
}

static void parseGeoInfo(inputImageStructure *inputImage, char *line, char *inputFile)
{
	double ReMajor, ReMinor, Rc, phic, H, rangeError;
	/*
	   Try to read range error if available else try no range error
	*/
	if (sscanf(line, "%lf%lf%lf%lf%lf%lf", &ReMajor, &ReMinor, &Rc, &phic, &H, &rangeError) != 6)
	{
		rangeError = 0.0;
		if (sscanf(line, "%lf%lf%lf%lf%lf", &ReMajor, &ReMinor, &Rc, &phic, &H) != 5)
			error("%s  %s %i %s", inputFile, "parseInputFile -- Missing geometric parameters at line:\n", line);
	}
	inputImage->par.rc = Rc * KMTOM - rangeError;
	inputImage->par.H = H * KMTOM;
	inputImage->par.ReMajor = ReMajor;
	inputImage->par.ReMinor = ReMinor;
	/* This is applied here, so it should not be applied elsewhere */
	inputImage->par.rangeError = rangeError;
	/* fprintf(stderr,"*** Range Error = %f\n",rangeError);*/
}

static int32_t parseControlPoints(FILE *fp, int32_t lineCount, inputImageStructure *inputImage, char *inputFile)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	extern double SLat;
	double lat, lon, x, y, sLat;
	int32_t eod, i;
	char line[512];
	/*
	  Init corner point arrays (allocate 4 for corner 1 for center and 1 for additional point).
	*/
	inputImage->latControlPoints = (double *)malloc((size_t)(NCONTROLPOINTS * sizeof(double)));
	inputImage->lonControlPoints = (double *)malloc((size_t)(NCONTROLPOINTS * sizeof(double)));
	/*
	  Input corner points
	*/
	if (SLat < -90.)
	{ /* SLat not already set, then guess */
		if (HemiSphere == NORTH)
			sLat = 70.0;
		else
			sLat = 71.0;
	}
	sLat = SLat;
	inputImage->minX = 1e30;
	inputImage->maxX = -1.e30;
	inputImage->minY = 1e30;
	inputImage->maxY = -1.e30;
	for (i = 1; i < 5; i++)
	{
		lineCount = getDataString(fp, lineCount, line, &eod);
		if (sscanf(line, "%lf %lf", &lat, &lon) != 2)
			error("%s %s %i %s %i", inputFile, "parseInputFile -- Missing parameters for control point: ", i, "at line: ", lineCount);
		inputImage->latControlPoints[i] = lat;
		inputImage->lonControlPoints[i] = lon;
		lltoxy1(lat, lon, &x, &y, Rotation, sLat);
		inputImage->minX = min(inputImage->minX, x);
		inputImage->maxX = max(inputImage->maxX, x);
		inputImage->minY = min(inputImage->minY, y);
		inputImage->maxY = max(inputImage->maxY, y);
	}
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%lf %lf", &lat, &lon) != 2)
		error("%s %i %s %i", "parseInputFile -- Missing parameters for control point: ", i, "at line: ", lineCount);
	inputImage->latControlPoints[0] = lat;
	inputImage->lonControlPoints[0] = lon;
	return (lineCount);
}

static void parsePixelSize(inputImageStructure *inputImage, char *line, char *inputFile)
{
	if (sscanf(line, "%lf %lf", &RangePixelSize, &AzimuthPixelSize) != 2)
		error("%s %s %s", inputFile, "parseInputFile -- Missing pixel size data at line : \n ", line);
	inputImage->azimuthPixelSize = (double)inputImage->nAzimuthLooks * AzimuthPixelSize;
	inputImage->rangePixelSize = (double)inputImage->nRangeLooks * RangePixelSize;
	inputImage->par.slpA = AzimuthPixelSize;
	inputImage->par.slpR = RangePixelSize;
	/*fprintf(stderr,"Range/Azimuth Pixel sizes %f %f\n",RangePixelSize,AzimuthPixelSize);*/
}

static void parseAscDesc(inputImageStructure *inputImage, char *line, char *inputFile)
{
	int32_t i;
	for (i = 0; i < 256 && line[i] != '\0'; i++)
		line[i] = (char)tolower(line[i]);
	/* These warnings are likely obsolete but left in anyway just in case */
	if (strstr(line, "descending") != NULL)
	{
		inputImage->passType = DESCENDING;
		fprintf(stderr, "*** Descending - ");
	}
	else if (strstr(line, "ascending") != NULL)
	{
		inputImage->passType = ASCENDING;
		fprintf(stderr, "*** Ascending - ");
	}
	else
		error("%s %s", inputFile, "parseInputFile: missing ascending or descending flag\n");
}

static void parseLookDir(inputImageStructure *inputImage, char *line, char *inputFile)
{
	int32_t i;
	for (i = 0; i < 256 && line[i] != '\0'; i++)
		line[i] = (char)tolower(line[i]);
	if (strstr(line, "left") != NULL)
	{
		inputImage->lookDir = LEFT;
		fprintf(stderr, " LEFT LOOKING DATA ****\n");
	}
	else if (strstr(line, "right") != NULL)
	{
		inputImage->lookDir = RIGHT;
		fprintf(stderr, " RIGHT LOOKING DATA ***\n");
	}
	else
		error("%s %s", inputFile, "no left/right in geodata file \n");
}

void parseInputFile(char *inputFile, inputImageStructure *inputImage)
{
	FILE *fp;
	int32_t lineCount, eod;
	SARData *par;
	stateV *sv;
	int32_t julDayInt;
	int32_t i;
	char line[1024];
	lineCount = 0;
	/* pull out pointers, for shorthand */
	par = &(inputImage->par);
	sv = &(inputImage->sv);
	/*	  Open file for input	*/
	fp = openInputFile(inputFile);
	/*	  Parse label	*/
	line[0] = '\0';
	fgetline(fp, line, 512); /* Read line */
	par->label = (char *)malloc(strlen(line) + 1);
	par->label[0] = '\0';
	strcat(par->label, line);
	fprintf(stderr, "\n----\npar->label %s", par->label);
	/*	  Parse date	*/
	fgetline(fp, line, 512); /* Read line */
	parseImageDate(inputImage, line, inputFile);
	julDayInt = julday(inputImage->month, inputImage->day, inputImage->year);
	lineCount = 2;
	/*	  Input nr,na,nlr,nla	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	parseSize(inputImage, line, inputFile);
	/*	  Input ReMajor,ReMinor, Rc, phic, H	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	parseGeoInfo(inputImage, line, inputFile);
	/*   Read corner points	*/
	lineCount = parseControlPoints(fp, lineCount, inputImage, inputFile);
	/*	  Compute pixels sizes fpr multilook images	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	parsePixelSize(inputImage, line, inputFile);
	/*	  parse descending/ascending flag	*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	parseAscDesc(inputImage, line, inputFile);
	/* 	   parse look dir	 */
	lineCount = getDataString(fp, lineCount, line, &eod); /* Get next line */
	parseLookDir(inputImage, line, inputFile);
	lineCount = getDataString(fp, lineCount, line, &eod); /* Get next line */
	/* This is an obsolete hold over */
	if (eod == TRUE || strstr(line, "state") == NULL)
		error("no state flag present in geodat");
	/* Get corrected time */
	lineCount = getDataString(fp, lineCount, line, &eod);
	parseTime(inputImage, line, eod);
	/* PRF */
	lineCount = getDataString(fp, lineCount, line, &eod);
	parsePRF(inputImage, line, eod);
	/* Wavelength */
	lineCount = getDataString(fp, lineCount, line, &eod);
	parseLambda(inputImage, line, eod);
	/* State vectors */
	lineCount = parseStateVectors(fp, lineCount, inputImage, inputFile);
	/*fprintf(stderr,"SV ts, n: %f %i\n",sv->t0,sv->nState);*/
	/* Get delta T and correct Time*/
	lineCount = getDataString(fp, lineCount, line, &eod);
	/* Parse deltaT if one exists */
	parseDeltaT(inputImage, julDayInt, line);
	inputImage->par.rn = inputImage->par.rc - inputImage->rangePixelSize * (inputImage->rangeSize - 1) * 0.5;
	inputImage->par.rf = inputImage->par.rc + inputImage->rangePixelSize * (inputImage->rangeSize - 1) * 0.5;
	/* end add */
	if (inputImage->lookDir == LEFT)
		par->lookDir = -1.0;
	else
		par->lookDir = 1.0;
	/*	  close up and exit        */
	fclose(fp);
	return;
}
