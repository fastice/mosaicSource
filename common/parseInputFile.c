#include "stdio.h"
#include "ctype.h"
#include <string.h>
#include <stdlib.h>
#include "math.h"
#include "mosaicSource/common/common.h"
//#include "mosaicSource/common/geojsonCode.h"
/*
   Parse a geodatXxY.in file - parseInputFile - with numerous supporting routines to parse data
 */

/* Parse data string and populate image date fields */
static void parseImageDate(inputImageStructure *inputImage, char *line, char *inputFile)
{
	char *tmp, *tmp1;
	int32_t m, done;
	int32_t mm, dd, yyyy;
	// fprintf(stderr, "%s %s", line, inputFile);
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

static void parseTime(inputImageStructure *inputImage, char *line, int32_t eod, char *format)
{
	if (sscanf(line, format, &(inputImage->par.hr), &(inputImage->par.min),
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


static int32_t parseStateVectorsGeojson(OGRFeatureH myFeature, inputImageStructure *inputImage)
{
	stateV *sv;
	int32_t nPos, nVel, i;
	const double *position, *velocity;
	sv = &(inputImage->sv);
	/*
	  If state flag set read state vector information and initialize
	  par structure
	*/
	sv->nState = OGR_F_GetFieldAsInteger(myFeature, OGR_F_GetFieldIndex(myFeature, "NumberOfStateVectors"));
	sv->deltaT = OGR_F_GetFieldAsInteger(myFeature, OGR_F_GetFieldIndex(myFeature, "StateVectorInterval"));
	if (sv->nState < MNST || sv->nState > MXST)
		error("\nInvalid number of state vectors - %i\n", sv->nState);
	sv->t0 = OGR_F_GetFieldAsDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "TimeOfFirstStateVector"));
	/*
	  Read state vector info - NOTE USING 1-nState indexing for polint
	*/
	// fprintf(stderr, "nState %i\n", sv->nState);
	for(i=1; i <= sv->nState; i++)
	{
		position = OGR_F_GetFieldAsDoubleList(myFeature, OGR_F_GetFieldIndex(myFeature, svTag(i, "Pos")), &nPos);
		velocity = OGR_F_GetFieldAsDoubleList(myFeature, OGR_F_GetFieldIndex(myFeature, svTag(i, "Vel")), &nVel);
		if(nPos != 3 || nVel != 3)
			error("parseStateVectorsGeojson: invalid number of position or velocity values %i %i\n", nPos, nVel);
		sv->x[i] = position[0];
		sv->y[i] = position[1];
		sv->z[i] = position[2];
		sv->vx[i] = velocity[0];
		sv->vy[i] = velocity[1];
		sv->vz[i] = velocity[2];
		// fprintf(stderr, "%i %s %s %i %lf %lf\n", i, svTag(i, "Pos"), svTag(i, "Vel"), sv->nState, sv->x[i], sv->vz[i]);
		sv->times[i] = sv->t0 + (i - 1) * sv->deltaT;
	}
}

static void updateDeltaT(inputImageStructure *inputImage, int32_t julDayInt, double deltaTCorrect) {

	SARData *par = &(inputImage->par);
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

static void parseDeltaT(inputImageStructure *inputImage, int32_t julDayInt, char *line)
{
	double deltaTCorrect;
	
	deltaTCorrect = 0.0;
	// Correction from SLC to ML coordinates
	deltaTCorrect += (inputImage->nAzimuthLooks - 1.0) * 0.5 / inputImage->par.prf;
	// Override value
	if (strstr(line, "deltaT") != NULL)
		sscanf((line + 6), "%lf\n", &deltaTCorrect);
	updateDeltaT(inputImage, julDayInt, deltaTCorrect);
}

static void parseDeltaTGeojson(OGRFeatureH myFeature, inputImageStructure *inputImage, int32_t julDayInt)
{
	double deltaTCorrect = 0.0;

	deltaTCorrect = OGR_F_GetFieldAsDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "deltaT"));
	// Correction from SLC to ML coordinates when no deltaT (<1e-6)
	if(fabs(deltaTCorrect) < 1e-6)
		deltaTCorrect += (inputImage->nAzimuthLooks - 1.0) * 0.5 / inputImage->par.prf;
	// Override value
	updateDeltaT(inputImage, julDayInt, deltaTCorrect);
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

static void computeControlPointsXY(inputImageStructure *inputImage)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	extern double SLat;
	double lat, lon, x, y, sLat;
	int i;
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
	else
		sLat = SLat;
	inputImage->minX = 1e30;
	inputImage->maxX = -1.e30;
	inputImage->minY = 1e30;
	inputImage->maxY = -1.e30;
	for (i = 1; i < 5; i++)
	{
		lat = inputImage->latControlPoints[i];
		lon = inputImage->lonControlPoints[i];
		lltoxy1(lat, lon, &x, &y, Rotation, sLat);
		// fprintf(stderr, "%lf %lf %lf %lf %lf\n", lat, lon, x, y, sLat);
		inputImage->minX = min(inputImage->minX, x);
		inputImage->maxX = max(inputImage->maxX, x);
		inputImage->minY = min(inputImage->minY, y);
		inputImage->maxY = max(inputImage->maxY, y);
	}
	// fprintf(stderr, "xMin, xMax, yMin, yMax %lf %lf %lf %lf\n",
	//		inputImage->minX, inputImage->maxX, inputImage->minY, inputImage->maxY);
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
	for (i = 1; i < 5; i++)
	{
		lineCount = getDataString(fp, lineCount, line, &eod);
		if (sscanf(line, "%lf %lf", &lat, &lon) != 2)
			error("%s %s %i %s %i", inputFile, "parseInputFile -- Missing parameters for control point: ", i, "at line: ", lineCount);
		inputImage->latControlPoints[i] = lat;
		inputImage->lonControlPoints[i] = lon;
	}
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%lf %lf", &lat, &lon) != 2)
		error("%s %i %s %i", "parseInputFile -- Missing parameters for control point: ", i, "at line: ", lineCount);
	inputImage->latControlPoints[0] = lat;
	inputImage->lonControlPoints[0] = lon;
	computeControlPointsXY(inputImage);
	return lineCount;
}

static void parseControlPointsGeoJson(OGRFeatureH myFeature, inputImageStructure *inputImage)
{
	double lat[5], lon[5];
	int32_t i, nItems;
	// Space for lat/long
	inputImage->latControlPoints = (double *)malloc((size_t)(NCONTROLPOINTS * sizeof(double)));
	inputImage->lonControlPoints = (double *)malloc((size_t)(NCONTROLPOINTS * sizeof(double)));
	// Get the geometry
	OGRGeometryH myGeometry = OGR_F_GetGeometryRef(myFeature);
	OGRGeometryH myRing = OGR_G_GetGeometryRef(myGeometry, 0);
	OGR_G_GetPoints(myRing, lat, sizeof(double), lon, sizeof(double), NULL, 0);
	// Remap points
	int32_t index[5] = {0, 0, 3, 1, 2};
	for (i = 1; i < 5; i++)
	{
		inputImage->latControlPoints[i] = lat[index[i]];
		inputImage->lonControlPoints[i] = lon[index[i]];
		fprintf(stderr, "%f %f\n",lat[index[i]],lon[index[i]]);
	}
	double const *ll = OGR_F_GetFieldAsDoubleList(myFeature, OGR_F_GetFieldIndex(myFeature, "CenterLatLon"), &nItems);
	inputImage->latControlPoints[0] = ll[0];
	inputImage->lonControlPoints[0] = ll[1];
	computeControlPointsXY(inputImage);
}

static void parsePixelSize(inputImageStructure *inputImage, char *line, char *inputFile)
{
	double azimuthPixelSize, rangePixelSize;
	if (sscanf(line, "%lf %lf", &rangePixelSize, &azimuthPixelSize) != 2)
		error("%s %s %s", inputFile, "parseInputFile -- Missing pixel size data at line : \n ", line);
	inputImage->azimuthPixelSize = (double)inputImage->nAzimuthLooks * azimuthPixelSize;
	inputImage->rangePixelSize = (double)inputImage->nRangeLooks * rangePixelSize;
	inputImage->par.slpA = azimuthPixelSize;
	inputImage->par.slpR = rangePixelSize;
}

static int determineAscDesc(inputImageStructure *inputImage, const char *line)
{
	int32_t lineLen = strlen(line);
	// Make case insensitive
	char lineLower[lineLen + 1];
	for (int i = 0; i <= lineLen; i++)
		lineLower[i] = (char)tolower(line[i]);
	/* These warnings are likely obsolete but left in anyway just in case */
	if (strstr(lineLower, "descending") != NULL)
	{
		inputImage->passType = DESCENDING;
		fprintf(stderr, "*** Descending - ");
		return TRUE;
	}
	else if (strstr(lineLower, "ascending") != NULL)
	{
		inputImage->passType = ASCENDING;
		fprintf(stderr, "*** Ascending - ");
		return TRUE;
	}
	return FALSE;
}

static void parseAscDesc(inputImageStructure *inputImage, char *line, char *inputFile)
{
	if (determineAscDesc(inputImage, line) == FALSE)
		error("%s %s", inputFile, "parseInputFile: missing ascending or descending flag\n");
}
static int determineLookDir(inputImageStructure *inputImage, const char *line)
{
	int32_t lineLen = strlen(line);
	// Make case insensitive
	char lineLower[lineLen + 1];
	for (int i = 0; i <= lineLen; i++)
		lineLower[i] = (char)tolower(line[i]);
	if (strstr(lineLower, "left") != NULL)
	{
		inputImage->lookDir = LEFT;
		fprintf(stderr, "LEFT LOOKING DATA ****\n");
	}
	else if (strstr(lineLower, "right") != NULL)
	{
		inputImage->lookDir = RIGHT;
		fprintf(stderr, "RIGHT LOOKING DATA ***\n");
		
	} else return FALSE;
	// Set the value for par
	if (inputImage->lookDir == LEFT)
		inputImage->par.lookDir = -1.0;
	else
		inputImage->par.lookDir = 1.0;
	return TRUE;
	
}

static void parseLookDir(inputImageStructure *inputImage, char *line, char *inputFile)
{
	if (determineLookDir(inputImage, line) == FALSE) {
		error("%s %s \n %s", inputFile, "parseLookDir: missing left/right flag\n", line);
	}
}

static char *checkGeojson(char *input)
{
	char *geojson;
	fprintf(stderr, "%s", input);
	if (strstr(input, ".geojson") != NULL)
	{
		geojson = input;
	}
	else if (strstr(input, ".in") != NULL)
	{
		size_t len = strlen(input);
		geojson = malloc(len + 6);
		geojson[0] = '\0';
		strncpy(geojson, input, len - 3);
		geojson[len - 3] = '\0';
		strcat(geojson, ".geojson");
	}
	if (access(geojson, F_OK) == 0) return geojson;
	return NULL;
}


void parseInputFile(char *inputFile, inputImageStructure *inputImage)
{
	FILE *fp;
	int32_t lineCount, eod;
	SARData *par;
	stateV *sv;
	int32_t julDayInt;
	int32_t i;
	char *geojson;
	char line[1024];
	lineCount = 0;
	/* pull out pointers, for shorthand */
	par = &(inputImage->par);
	sv = &(inputImage->sv);
	geojson = checkGeojson(inputFile);
	// If a geojson version exist use that
	if(geojson != NULL) {
		fprintf(stderr, "\nParsing %s\n", geojson);
		parseGeojson(geojson, inputImage);
		//dumpParsedInput(inputImage, stderr);
		return;
	}
	/*	  Open file for input	*/
	fp = openInputFile(inputFile);
	/*	  Parse label	*/
	line[0] = '\0';
	fgetline(fp, line, 512); /* Read line */
	par->label = (char *)malloc(strlen(line) );
	for(i=14; i < strlen(line)-1; i++) par->label[i-14] = line[i];
	par->label[i-14] = '\0';
	// strcat(par->label, line);
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
	parseTime(inputImage, line, eod, "%i %i %lf");
	fprintf(stderr, "Time %i %i %f\n", inputImage->par.hr, inputImage->par.min, inputImage->par.sec);	
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
	/*	  close up and exit        */
	fclose(fp);
	//dumpParsedInput(inputImage, stderr);
	return;
}


void parseGeojson(char *inputFile, inputImageStructure *inputImage)
{
	SARData *par;
	stateV *sv;
	float dum1, dum2, dum3, dum4, dum5;
	int32_t idum1, idum2, idum3, idum4, idum5;
	int32_t julDayInt;
	// Check file exists
	if (access(inputFile, F_OK) != 0)
		error("parseGeojson: File does not exist %s", inputFile);
	// Get driver
	GDALDriverH myDriver = (GDALDriverH) GDALGetDriverByName("geojson");
	// OPen data sets
	GDALDatasetH myData = (GDALDatasetH) GDALOpenEx(inputFile, GDAL_OF_VECTOR, NULL, NULL, NULL);
	// Get the layer
	OGRLayerH myLayer = OGR_DS_GetLayer(myData, 0);
	// Create a feature
	OGRFeatureH myFeature = OGR_L_GetFeature(myLayer, 0);
	// Pull out par and sv structures to streamline code
	par = &(inputImage->par);
	sv = &(inputImage->sv);
	// Get label
	const char *imageName = OGR_F_GetFieldAsString(myFeature, OGR_F_GetFieldIndex(myFeature, "ImageName"));
	par->label = strcpy(calloc(sizeof(char), strlen(imageName) + 1), imageName);
	// Get date
	OGR_F_GetFieldAsDateTimeEx(myFeature, OGR_F_GetFieldIndex(myFeature, "Date"),
							   &(inputImage->year), &(inputImage->month), &(inputImage->day), &idum1, &idum2, &dum1, &idum3);
	// Compute julian data
	julDayInt = julday(inputImage->month, inputImage->day, inputImage->year);
	//	Input nr,na,nlr,nla
	inputImage->rangeSize = OGR_F_GetFieldAsInteger(myFeature, OGR_F_GetFieldIndex(myFeature, "MLRangeSize"));
	inputImage->azimuthSize = OGR_F_GetFieldAsInteger(myFeature, OGR_F_GetFieldIndex(myFeature, "MLAzimuthSize"));
	inputImage->nRangeLooks = OGR_F_GetFieldAsInteger(myFeature, OGR_F_GetFieldIndex(myFeature, "NumberRangeLooks"));
	inputImage->nAzimuthLooks = OGR_F_GetFieldAsInteger(myFeature, OGR_F_GetFieldIndex(myFeature, "NumberAzimuthLooks"));
	//	  Input ReMajor,ReMinor, Rc, phic, H	]
	double rangeError = OGR_F_GetFieldAsDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "RangeErrorCorrection"));
	inputImage->par.rc = OGR_F_GetFieldAsDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "MLCenterRange")) - rangeError;
	inputImage->par.H = OGR_F_GetFieldAsDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "SpaceCraftAltitude"));
	inputImage->par.ReMajor = OGR_F_GetFieldAsDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "EarthRadiusMajor")) * MTOKM;
	inputImage->par.ReMinor = OGR_F_GetFieldAsDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "EarthRadiusMinor")) * MTOKM;
	inputImage->par.rangeError = rangeError;
	// Parse Control points
	parseControlPointsGeoJson(myFeature, inputImage);
	// Sizes for ingle and multilook images
	inputImage->par.slpR = OGR_F_GetFieldAsDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "SLCRangePixelSize"));
	inputImage->par.slpA = OGR_F_GetFieldAsDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "SLCAzimuthPixelSize"));
	inputImage->rangePixelSize = (double)inputImage->nRangeLooks * inputImage->par.slpR;
	inputImage->azimuthPixelSize = (double)inputImage->nAzimuthLooks * inputImage->par.slpA;
	//
	if (determineAscDesc(inputImage,
						 OGR_F_GetFieldAsString(myFeature, OGR_F_GetFieldIndex(myFeature, "PassType"))) == FALSE)
		error("parseGeojson: could not parse TrackDirection");
	if (determineLookDir(inputImage,
						 OGR_F_GetFieldAsString(myFeature, OGR_F_GetFieldIndex(myFeature, "LookDirection"))) == FALSE)
		error("parseGeojson: could not parse LookDirection");
	//  Get corrected time
	const char* timeString = OGR_F_GetFieldAsString(myFeature, OGR_F_GetFieldIndex(myFeature, "CorrectedTime"));
	parseTime(inputImage, (char *) timeString, FALSE, "%2u %2u %lf");
	fprintf(stderr, "Time %i %i %f\n", inputImage->par.hr, inputImage->par.min, inputImage->par.sec);
	error("STOP");
	// PRF
	inputImage->par.prf = OGR_F_GetFieldAsDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "PRF"));
	// Wavelength
	inputImage->par.lambda = OGR_F_GetFieldAsDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "Wavelength"));
	// State vectors
	parseStateVectorsGeojson(myFeature, inputImage);
	// Parse deltaT and correct time.
	parseDeltaTGeojson(myFeature, inputImage, julDayInt);
	inputImage->par.rn = inputImage->par.rc - inputImage->rangePixelSize * (inputImage->rangeSize - 1) * 0.5;
	inputImage->par.rf = inputImage->par.rc + inputImage->rangePixelSize * (inputImage->rangeSize - 1) * 0.5;
	// Clean up
	OGRReleaseDataSource(myData);

}


void dumpParsedInput(inputImageStructure *inputImage, FILE *fp) {
	SARData *par;
	int i;
	stateV *sv;
	par = &(inputImage->par);
	sv = &(inputImage->sv);
	// Get label
	fprintf(fp, "Image Name %s\n", par->label);
	fprintf(fp, "Date %i %i %i\n", inputImage->year, inputImage->month, inputImage->day);
	fprintf(fp, "nr,na,nlr,nla  %i %i %i %i\n", inputImage->rangeSize, inputImage->azimuthSize,
			inputImage->nRangeLooks, inputImage->nAzimuthLooks);
	//	  Input ReMajor,ReMinor, Rc, phic, H	]
	fprintf(fp, "r, H, ReMa, ReM, rangeErr %lf %lf %lf %lf %lf\n",
			inputImage->par.rc, inputImage->par.H, inputImage->par.ReMajor, inputImage->par.ReMinor, inputImage->par.rangeError);
	// Parse Control points
	for(i=0;  i < 5; i++) {
		fprintf(fp,"%lf %lf\n", inputImage->latControlPoints[i], inputImage->lonControlPoints[i]);
	}
	// Sizes for ingle and multilook images
	fprintf(fp, "SL delta R/A %lf %lf ML delta R/A %lf %lf\n", inputImage->par.slpR , inputImage->par.slpA,
		inputImage->rangePixelSize, inputImage->azimuthPixelSize);
	//
	fprintf(fp, "Lookdir %i %f\n", inputImage->lookDir, par->lookDir);
	fprintf(fp, "Pass type %i\n", inputImage->passType);
	//  Get corrected time
	fprintf(fp, "time %i %i %12.8lf julday %lf\n", inputImage->par.hr, inputImage->par.min, inputImage->par.sec, inputImage->julDay);
	// PRF
	fprintf(fp, "PRF, lambda %lf %lf\n", inputImage->par.prf, inputImage->par.lambda);
	// State vectors
	fprintf(fp, "%i %lf %lf\n", sv->nState, sv->deltaT, sv->t0);
	for(i=1; i <=sv->nState; i++) {
		fprintf(fp, "%i %lf %lf %lf   %lf %lf %lf\n", i, sv->x[i], sv->y[i], sv->z[i], sv->vx[i], sv->vy[i], sv->vz[i]);
	}
}
