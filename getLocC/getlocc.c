#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#undef MNST
#include "math.h"
//#include "mosaicSource/common/geojsonCode.h"
#include "mosaicSource/common/common.h"
/*
  Create geodat file from .par file from clw processor.
*/

static void readArgs(int argc, char *argv[], int32_t *nlr, int32_t *nla, int32_t *noffset, double *squint, char **parFile,
					 char **geodatFile, int32_t *passType, double *lookDir, int32_t *squintTime, int32_t *ersFix,
					 int32_t *surveyFlag, double *lambda);
static void usage();
static void echoSarD(FILE *fp, SARData *sarD, stateV *sv, int32_t nlr, int32_t nla);
/* These routines are largely fixed to this program and not used elsewhere - so not bothering with include file*/

void correctTime(SARData *sarD, double *squint, double noffset, double *tskew, double *toffset, int32_t squintTime);
void centerLL(SARData *sarD, stateV *sv, int32_t nla, double *lat, double *lon, double deltaT);
void glatlon(SARData *sarD, stateV *sv, int32_t nlr, int32_t nla, double ***lat1, double ***lon1, int32_t ma, int32_t mr, double deltaT);
void formatGeoJSON(char *myFile, int indent);
static void writeCoordinates(OGRFeatureH myFeature,
							 double **lat, double **lon, double latc, double lonc,
							 int32_t passType, FILE *fp);

static void writeTimeOffsets(OGRFeatureH myFeature, SARData *sarD,
							 int32_t noffset, int32_t squintTime, double squint, int32_t passType, FILE *fp);

static void writePassType(OGRFeatureH myFeature, int32_t lookDir, int32_t passType, FILE *fp);
static void computeSARLookGeom(OGRFeatureH myFeature, SARData *sarD, double latc, FILE *fp);
static void echoSarDgeojson(OGRFeatureH myFeature, SARData *sarD, stateV *sv, int32_t nlr, int32_t nla);
static void writeTimeAndWavelength(OGRFeatureH myFeature, SARData *sarD, FILE *fp);
static void writeStateVectors(OGRFeatureH myFeature, stateV *sv, FILE *fp);
/*
   Global variables definitions (NOT USED, NEED FOR LINKING ERS CODE
*/
int32_t RangeSize = RANGESIZE;				/* Range size of complex image */
int32_t AzimuthSize = AZIMUTHSIZE;			/* Azimuth size of complex image */
int32_t BufferSize = BUFFERSIZE;			/* Size of nonoverlap region of the buffer */
int32_t BufferLines = 512;					/* # of lines of nonoverlap in buffer */
int32_t HemiSphere = NORTH;
double Rotation = 45.;
double SLat = -91.0;
void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
void *lBuf1, *lBuf2, *lBuf3, *lBuf4;
char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
int32_t llConserveMem = 999; /* Kluge to maintain backwards compat 9/13/06 */

// Hard code wkt to avoid proj lib path issues.
const char *wkt = "GEOGCRS[\"WGS 84\", DATUM[\"World Geodetic System 1984\", \
        ELLIPSOID[\"WGS 84\",6378137,298.257223563, \
            LENGTHUNIT[\"metre\",1]]], PRIMEM[\"Greenwich\",0, \
        ANGLEUNIT[\"degree\",0.0174532925199433]],\
    CS[ellipsoidal,2], AXIS[\"geodetic latitude (Lat)\",north,\
            ORDER[1], ANGLEUNIT[\"degree\",0.0174532925199433]],\
        AXIS[\"geodetic longitude (Lon)\",east,\
            ORDER[2],  ANGLEUNIT[\"degree\",0.0174532925199433]],\
    USAGE[SCOPE[\"Horizontal component of 3D system.\"],AREA[\"World.\"],BBOX[-90,-180,90,180]],\
    ID[\"EPSG\",4326]]";

int main(int argc, char *argv[])
{
	FILE *fp;
	int32_t nlr, nla; /* Number of range and azimuth looks */
	int32_t noffset;
	char *parFile, *geodatFile; /* Input and output files */
	char *geojsonFile;
	int32_t passType; /* 0 -> descending, 1 -> ascending */
	double tdum1 = 0., tdum2 = 0., tdum3 = 0.;
	double **lat, **lon;
	double latc, lonc;
	double lookDir;
	double lambda;
	double squint = 0.;
	int32_t ersFix;
	int32_t squintTime, surveyFlag;
	int32_t i;
	int32_t nMLR, nMLA;
	double deltaRgML;
	double mldT;
	SARData sarD;
	stateV sv;
	// Parse the arguments
	readArgs(argc, argv, &nlr, &nla, &noffset, &squint, &parFile, &geodatFile, &passType, &lookDir,
			 &squintTime, &ersFix, &surveyFlag, &lambda);

	geojsonFile = calloc(sizeof(char), strlen(geodatFile) + 5);
	geojsonFile['\0'];
	strncpy(geojsonFile, geodatFile, strlen(geodatFile) - 3);
	strcat(geojsonFile, ".geojson");
	fprintf(stderr, "geosjsonFile: %s geodatFile: %s\n", geojsonFile, geodatFile);
	sarD.lambda = lambda;
	/* Read old style gamm cw par file - no keywords */
	readOldPar(parFile, &sarD, &sv);
	nMLR = sarD.nSlpR / nlr; /* Mulitlook size */
	nMLA = sarD.nSlpA / nla;
	deltaRgML = nlr * sarD.slpR;
	/*
	  Make sure range convention obeyed.
	*/
	sarD.lookDir = lookDir;
	/* Sort of obsolete, but just in case */
	if (surveyFlag == TRUE)
		sarD.prf /= 8.;
	/*
	  Make sure ranges centered on multilook pixels rather than sl pixels
	*/
	sarD.rn += (nlr - 1) * sarD.slpR * 0.5;
	sarD.rf = sarD.rn + (nMLR - 1) * deltaRgML;
	sarD.rc = (sarD.rn + sarD.rf) * 0.5;
	/* Nearly obsolete, keeping just in case */
	if (ersFix == TRUE)
	{
		fprintf(stderr, "*** APPLYING ERS FIX *** \n");
		sarD.rn -= 54.0;
		sarD.rc -= 54.0;
		sarD.rf -= 54.0;
	}
	/* Start writing geodat file */
	GDALAllRegister();
	OGRDataSourceH myDS = getGeojsonDataSet(geojsonFile);
	OGRFeatureDefnH geoTemplate = createFeatureDef(sv.nState);
	OGRFeatureH myFeature = OGR_F_Create(geoTemplate);
	OGRSpatialReferenceH mySRS = OSRNewSpatialReference(wkt);
	// Replaced this with hard code
	//OGRErr srsErr = OSRImportFromEPSG(mySRS, 4326);
	OGRLayerH myLayer = OGR_DS_CreateLayer(myDS, "", mySRS, wkbPolygon, NULL);
	// Write the basic data that appears near the front of the old format
	echoSarDgeojson(myFeature, &sarD, &sv, nlr, nla);
	// Open old type geodate
	fp = fopen(geodatFile, "w");
	// Write the SAR data
	echoSarD(fp, &sarD, &sv, nlr, nla);
	/*
	  Correct starting time for squint and offset
	*/
	writeTimeOffsets(myFeature, &sarD, noffset, squintTime, squint, passType, fp);
	// Range/azimuth sizes
	fprintf(fp, ";\n; rangesize,azimuthsize,nrangelooks,nazimuthlooks\n;\n");
	fprintf(fp, "%i  %i  %i  %i\n", nMLR, nMLA, nlr, nla);
	OGR_F_SetFieldInteger(myFeature, OGR_F_GetFieldIndex(myFeature, "MLRangeSize"), nMLR);
	OGR_F_SetFieldInteger(myFeature, OGR_F_GetFieldIndex(myFeature, "MLAzimuthSize"), nMLA);
	// Compute center ll and SAR geometry
	centerLL(&sarD, &sv, nla, &latc, &lonc, 0.0);
	computeSARLookGeom(myFeature, &sarD, latc, fp);
	/*
	  Note added mldT add, then subtract to be consistent time to glat/lon corrected  for center of multiLook pixel.
	*/
	mldT = 0.5 * (nla - 1);
	correctTime(&sarD, &tdum2, mldT, &tdum3, &tdum1, squintTime);
	glatlon(&sarD, &sv, nlr, nla, &lat, &lon, 3, 3, 0.0);
	correctTime(&sarD, &tdum2, -mldT, &tdum3, &tdum1, squintTime); /* Undo correction so time is SLC not ML */
	writeCoordinates(myFeature, lat, lon, latc, lonc, passType, fp);
	// Write single look resolutions to old format, geojson handled by echo
	fprintf(fp, ";\n; Range/azimuth single look pixel sizes \n;\n");
	fprintf(fp, "%f  %f\n", sarD.slpR, sarD.slpA);
	// Byte order - fix MSB for now
	OGR_F_SetFieldString(myFeature, OGR_F_GetFieldIndex(myFeature, "ByteOrder"), "MSB");
	// Write pass type and lookdir
	writePassType(myFeature, lookDir, passType, fp);
	// Time and wavelength
	writeTimeAndWavelength(myFeature, &sarD, fp);
	// Output the state vectors
	writeStateVectors(myFeature, &sv, fp);
	// Add Feature to layer and close to force write
	int err = OGR_L_CreateFeature(myLayer, myFeature);
	GDALClose(myDS);
	// Reformat the ouput
	formatGeoJSON(geojsonFile, 2);
}

static void readArgs(int argc, char *argv[], int32_t *nlr, int32_t *nla, int32_t *noffset, double *squint, char **parFile,
					 char **geodatFile, int32_t *passType, double *lookDir, int32_t *squintTime, int32_t *ersFix,
					 int32_t *surveyFlag, double *lambda)
{
	int32_t i;
	if (argc < 9 || argc > 13)
		usage(); /* Check number of args */
	*squintTime = FALSE;
	*ersFix = FALSE;
	*surveyFlag = FALSE;
	*lambda = LAMBDAERS1;
	if (argc >= 10)
	{
		for (i = 1; i < argc - 8; i++)
		{
			if (strstr(argv[i], "squintTime") != NULL)
				*squintTime = TRUE;
			else if (strstr(argv[i], "ersFix") != NULL)
				*ersFix = TRUE;
			else if (strstr(argv[i], "survey") != NULL)
				*surveyFlag = TRUE;
			else if (strstr(argv[i], "lambda") != NULL)
			{
				i++;
				sscanf(argv[i], "%lf", lambda);
			}
			else
				usage();
		}
	}
	sscanf(argv[argc - 8], "%i", nlr);
	sscanf(argv[argc - 7], "%i", nla);
	sscanf(argv[argc - 6], "%i", noffset);
	// fprintf(stderr, "%i %i %i\n", *nlr, *nla, *noffset);
	sscanf(argv[argc - 5], "%lf", squint);
	*parFile = argv[argc - 4];
	*geodatFile = argv[argc - 3];
	sscanf(argv[argc - 2], "%i", passType);
	sscanf(argv[argc - 1], "%lf", lookDir);
	if (*lookDir < 0.0)
		*lookDir = -1.0;
	else
		*lookDir = 1.0;
	return;
}

static void usage()
{
	error("getloc -lambda lambda -ersFix -squintTime nlr nla noffset squintAngle parfile outfile\n\n"
		  " lambda = wavelength\n"
		  " passtype = 0 for desc and 1 for right\n"
		  " squintTime flag to specify squint in time\n"
		  "lookDir = +1.0 for right and -1.0 for left");
}

/* This program writes the largely comment information to the geodat file */
static void echoSarD(FILE *fp, SARData *sarD, stateV *sv, int nlr, int nla)
{
	char *months[12] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
	fprintf(fp, "; Image name: %s\n", sarD->label);
	fprintf(fp, "; Image date: %02d %s %4i\n", sarD->day, months[sarD->month - 1], sarD->year);
	fprintf(fp, "; Image time: %i %i %f\n", sarD->hr, sarD->min, sarD->sec);
	fprintf(fp, "; Nominal center lat,lon: %f %f\n", 0., 0.);
	fprintf(fp, "; track direction: %f\n", 0.);
	fprintf(fp, "; S/C altitude: %f\n", sarD->H);
	fprintf(fp, "; Average height above terrain: %f\n", 0.0);
	fprintf(fp, "; Vel along track: %f\n", 0.0);
	fprintf(fp, "; PRF :   %f\n", sarD->prf);
	fprintf(fp, "; near/cen/far range : %f %f %f\n", sarD->rn, sarD->rc, sarD->rf);
	fprintf(fp, "; Range pixel spacing :   %f\n", sarD->slpR * (double)(nlr));
	fprintf(fp, "; Number of looks (rg,az) :   %i %i\n", nlr, nla);
	fprintf(fp, "; Azimuth pixel spacing :   %f\n", sarD->slpA * (double)nla);
	fprintf(fp, "; Number of pixels (rg,az) :  %i  %i\n", sarD->nSlpR / nlr, sarD->nSlpA / nla);
	fprintf(fp, "; Number of state vectors :   %i\n", sv->nState);
	fprintf(fp, "; Start time of state vectors :   %f\n", sv->t0);
	fprintf(fp, "; Interval between 2 state vectors :   %f\n", sv->deltaT);
	fprintf(fp, "; Look direction  :   %f\n", sarD->lookDir);
}

static void writeStateVectors(OGRFeatureH myFeature, stateV *sv, FILE *fp)
{
	/*
	  Write state vectors
	*/
	double position[3], velocity[3];
	// Parameter used to fix timing errors used Radarsat and ALOS, code con
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "deltaT"), 0.0);
	// State vector time info
	OGR_F_SetFieldInteger(myFeature, OGR_F_GetFieldIndex(myFeature, "NumberOfStateVectors"), sv->nState);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "TimeOfFirstStateVector"), sv->t0);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "StateVectorInterval"), sv->deltaT);
	fprintf(fp, "; number of state vectors\n%i\n", sv->nState);
	fprintf(fp, "; time of first vector \n%f\n", sv->t0);
	fprintf(fp, "; state vector interval \n%f\n", sv->deltaT);
	fprintf(fp, "; state vectors \n");
	// Write the state vectors
	for (int32_t i = 1; i <= sv->nState; i++)
	{
		fprintf(fp, "%.9e  %.9e  %.9e\n", sv->x[i], sv->y[i], sv->z[i]);
		fprintf(fp, "%.9e  %.9e  %.9e\n", sv->vx[i], sv->vy[i], sv->vz[i]);
		position[0] = sv->x[i];
		position[1] = sv->y[i];
		position[2] = sv->z[i];
		velocity[0] = sv->vx[i];
		velocity[1] = sv->vy[i];
		velocity[2] = sv->vz[i];
		OGR_F_SetFieldDoubleList(myFeature, OGR_F_GetFieldIndex(myFeature, svTag(i, "Pos")), 3, position);
		OGR_F_SetFieldDoubleList(myFeature, OGR_F_GetFieldIndex(myFeature, svTag(i, "Vel")), 3, velocity);
	}
}

static void writeTimeAndWavelength(OGRFeatureH myFeature, SARData *sarD, FILE *fp)
{
	char *timeString;
	fprintf(fp, ";\n; Flag to indicate state vectors and associated data\n;\n");
	fprintf(fp, "state\n");
	fprintf(fp, "; time after squint and skew corrections\n %i %i %11.7lf\n", sarD->hr, sarD->min, sarD->sec);
	fprintf(fp, "; prf \n%f\n", sarD->prf);
	fprintf(fp, "; wavelength\n%f\n", sarD->lambda);
	timeString = calloc(sizeof(char), 40);
	sprintf(timeString, "%02d %02d %011.8lf", sarD->hr, sarD->min, sarD->sec);
	OGR_F_SetFieldString(myFeature, OGR_F_GetFieldIndex(myFeature, "CorrectedTime"), timeString);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "Wavelength"), sarD->lambda);
}

static void computeSARLookGeom(OGRFeatureH myFeature, SARData *sarD, double latc, FILE *fp)
{
	double theta, phic, rearth;
	rearth = earthRadius(latc * DTOR, EMINOR, EMAJOR) * KMTOM;
	theta = sarD->rc * sarD->rc + 2.0 * sarD->H * rearth + sarD->H * sarD->H;
	theta = theta / (2.0 * sarD->rc * (rearth + sarD->H));
	theta = acos(theta);
	fprintf(stdout, "theta/latc %f %f\n", theta * RTOD, latc);
	phic = ((rearth + sarD->H) / rearth) * sin(theta);
	phic = asin(phic) * RTOD;
	/* Ouput info */
	fprintf(fp, ";\n; ReMajor, ReMinor, Rc, phic, h\n;\n");
	fprintf(fp, "%12.6f %12.6f %15.10f %8.4f %13.7f\n", EMAJOR, EMINOR, sarD->rc * MTOKM, phic, sarD->H * MTOKM);

	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "EarthRadiusMajor"), EMAJOR * KMTOM);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "EarthRadiusMinor"), EMINOR * KMTOM);
	// OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "MLRangeCenter"), sarD->rc * MTOKM);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "MLIncidenceCenter"), phic);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "SpaceCraftAltitude"), sarD->H);
}

static void writeTimeOffsets(OGRFeatureH myFeature, SARData *sarD,
							 int32_t noffset, int32_t squintTime, double squint, int32_t passType, FILE *fp)
{
	double toffset, tskew;
	char *passTypeString;
	//
	// Old format stuf
	correctTime(sarD, &squint, (double)noffset, &tskew, &toffset, squintTime);
	fprintf(fp, "; Offset of first recordin complex image (s) : %f\n", toffset);
	fprintf(fp, "; Skew offset (s), squint (deg) : %f  %lf\n", tskew, squint);
	/* Indicate pass type - asc/desc */
	if (passType == DESCENDING)
		fprintf(fp, ";\n; Descending Pass \n");
	else if (passType == ASCENDING)
		fprintf(fp, ";\n; Ascending Pass \n");
	else
		error("getlocc: Invalid Pass Type");
	//
	// New format
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "TimeToFirstSLCSample"), toffset);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "SkewOffset"), tskew);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "Squint"), squint);
}

static void writePassType(OGRFeatureH myFeature, int32_t lookDir, int32_t passType, FILE *fp)
{
	char *passTypeString, *lookDirString;

	if (passType == DESCENDING)
	{
		passTypeString = "Descending";
		fprintf(fp, ";\ndescending\n");
	}
	else
	{
		passTypeString = "Ascending";
		fprintf(fp, ";\nascending\n");
	}
	/* Look dir */
	fprintf(fp, ";\n; Look direction\n;\n");
	if (lookDir < 0)
	{
		fprintf(fp, "left\n");
		lookDirString = "Left";
	}
	else
	{
		fprintf(fp, "right\n");
		lookDirString = "Right";
	}
	// Save geojson fields
	OGR_F_SetFieldString(myFeature, OGR_F_GetFieldIndex(myFeature, "PassType"), passTypeString);
	OGR_F_SetFieldString(myFeature, OGR_F_GetFieldIndex(myFeature, "LookDirection"), lookDirString);
}

static void writeCoordinates(OGRFeatureH myFeature,
							 double **lat, double **lon, double latc, double lonc,
							 int32_t passType, FILE *fp)
{
	OGRGeometryH myGeometry;
	// Order corner points for ll,lr.. convention
	// Note for the geojson this will order the points ll, lr, ur, ul
	fprintf(fp, ";\n; ll,lr,ul,ur\n;\n");
	switch (passType)
	{
	case DESCENDING:
	{
		fprintf(fp, "%15.10lf %15.10lf\n", lat[2][2], lon[2][2]);
		fprintf(fp, "%15.10lf %15.10lf\n", lat[2][0], lon[2][0]);
		fprintf(fp, "%15.10lf %15.10lf\n", lat[0][2], lon[0][2]);
		fprintf(fp, "%15.10lf %15.10lf\n", lat[0][0], lon[0][0]);
		// ll, ul, ur, lr
		double lat1[4] = {lat[2][2], lat[0][2], lat[0][0], lat[2][0]};
		double lon1[4] = {lon[2][2], lon[0][2], lon[0][0], lon[2][0]};
		myGeometry = createGeometry(lat1, lon1);
		break;
	}
	case ASCENDING:
	{
		fprintf(fp, "%15.10lf %15.10lf\n", lat[0][0], lon[0][0]);
		fprintf(fp, "%15.10lf %15.10lf\n", lat[0][2], lon[0][2]);
		fprintf(fp, "%15.10lf %15.10lf\n", lat[2][0], lon[2][0]);
		fprintf(fp, "%15.10lf %15.10lf\n", lat[2][2], lon[2][2]);
		// ll, ul, ur, lr
		double lat1[4] = {lat[0][0], lat[2][0], lat[2][2], lat[0][2]};
		double lon1[4] = {lon[0][0], lon[2][0], lon[2][2], lon[0][2]};
		// double lat1[4] = {lat[0][2], lat[0][0], lat[2][0], lat[2][2]};
		//  double lon1[4] = {lon[0][2], lon[0][0], lon[2][0], lon[2][2]};
		myGeometry = createGeometry(lat1, lon1);
		break;
	}
	default:
		error("*** INVALID PASS TYPE ***\n");
	}
	double ll[2] = {latc, lonc};
	fprintf(fp, "%15.10lf %15.10lf\n", latc, lonc);
	OGR_F_SetGeometryDirectly(myFeature, myGeometry);
	OGR_F_SetFieldDoubleList(myFeature, OGR_F_GetFieldIndex(myFeature, "CenterLatLon"), 2, ll);
}

void formatGeoJSON(char *myFile, int indent)
{
	// Much longer than string,
	int maxChar = 50000;
	char geojsonString[maxChar], *tmp;
	char formattedGeoJSON[maxChar];
	FILE *fp;
	// Create a new string to hold the formatted GeoJSON.
	// char *formattedGeoJSON = malloc(strlen(geojsonString) + indent * 2 + 1);
	fp = fopen(myFile, "r");
	tmp = fgets(geojsonString, maxChar, fp);
	fclose(fp);
	// Reformat the string with line breaks and indentation
	for (int j = 0; j < maxChar; j++)
		formattedGeoJSON[j] = ' ';
	// Iterate over the characters in the GeoJSON string, adding line breaks and indents as needed.
	int i = 0;
	int k = 0;
	int skip = 0;
	int braces = 0;
	while (geojsonString[i] != '\0')
	{
		if (geojsonString[i] == '[')
			skip++;
		else if (geojsonString[i] == ']')
			skip--;
		if (geojsonString[i] == '{')
		{
			braces++;
			// Add a line break and indent.
			formattedGeoJSON[k++] = geojsonString[i];
			formattedGeoJSON[k++] = '\n';
			for (int j = 0; j < braces * indent - 1; j++)
				formattedGeoJSON[k++] = ' ';
		}
		else if (geojsonString[i] == ',')
		{
			// Add a comma and a space.
			formattedGeoJSON[k++] = ',';
			if (skip == 0)
			{
				formattedGeoJSON[k++] = '\n';
				for (int j = 0; j < braces * indent - 1; j++)
					formattedGeoJSON[k++] = ' ';
			}
		}
		else if (geojsonString[i] == '}')
		{
			// Add a line break and indent.
			formattedGeoJSON[k++] = '\n';
			braces--;
			if (braces > 0)
			{
				for (int j = 0; j < braces * indent - 1; j++)
					formattedGeoJSON[k++] = ' ';
			}
			formattedGeoJSON[k++] = '}';
		}
		else if (geojsonString[i] == ' ')
		{
			if (skip == 0)
				formattedGeoJSON[k++] = ' ';
		}
		else
		{
			// Just copy the character.
			formattedGeoJSON[k++] = geojsonString[i];
		}
		i++;
	}

	// Terminate the string.
	formattedGeoJSON[k] = '\0';
	// Print the formatted GeoJSON.
	// Overwrite original
	fp = fopen(myFile, "w");
	fprintf(fp, "%s", formattedGeoJSON);
	fclose(fp);
	// Free the memory allocated for the formatted GeoJSON.
}

static void echoSarDgeojson(OGRFeatureH myFeature, SARData *sarD, stateV *sv, int32_t nlr, int32_t nla)
{
	OGR_F_SetFieldString(myFeature, OGR_F_GetFieldIndex(myFeature, "ImageName"), sarD->label);
	OGR_F_SetFieldDateTimeEx(myFeature, OGR_F_GetFieldIndex(myFeature, "Date"), sarD->year, sarD->month, sarD->day, 0, 0, 0., 0);
	float sec = (float)sarD->sec;
	OGR_F_SetFieldDateTimeEx(myFeature, OGR_F_GetFieldIndex(myFeature, "NominalTime"), 0, 0, 0, sarD->hr, sarD->min, sec, 0);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "PRF"), sarD->prf);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "MLNearRange"), sarD->rn);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "MLCenterRange"), sarD->rc);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "MLFarRange"), sarD->rf);
	// Keep for future need
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "RangeErrorCorrection"), 0.);

	OGR_F_SetFieldInteger(myFeature, OGR_F_GetFieldIndex(myFeature, "NumberRangeLooks"), nlr);
	OGR_F_SetFieldInteger(myFeature, OGR_F_GetFieldIndex(myFeature, "NumberAzimuthLooks"), nla);

	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "SLCRangePixelSize"), sarD->slpR);
	OGR_F_SetFieldDouble(myFeature, OGR_F_GetFieldIndex(myFeature, "SLCAzimuthPixelSize"), sarD->slpA);
}
