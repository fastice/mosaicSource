#include "stdio.h"
#include <stdlib.h>
#include <string.h>
#include "mosaicSource/common/common.h"
#include "math.h"
#include "clib/standard.h"
#include <sys/types.h>

#define DEM 1
#define VELOCITY 2

//void readXYDEMGeoInfo(char *xyFile, void *xydem, int32_t resetProjection);
static void readXYGeodatFile(char *xyFile, double *x0, double *y0, double *deltaX, double *deltaY, 
							int32_t *xSize, int32_t *ySize, double *rot, int32_t *hemisphere, double *stdLat);
static float **readXYImageCrop(char *xyFile, void *xyImage, double xmin, double xmax, double  ymin, double ymax, int type);
static float **readXYImageGDALCropped(char *xyFile, void *xyImage, double xMin, double xMax,
						   		double yMin, double yMax, int32_t type);
static void getCropBounds(GDALDatasetH hDataset, double xMin, double xMax, double yMin, double yMax,
				   int *xOffset, int *yOffset, void *obj, int32_t *flip, int32_t type);
static void readXYGeoInfo(char *xyFile, void *xyImage, int32_t resetProjection, int type);
static void readXYGeoInfoGDAL(char *xyFile, void *xyImage, int type);
static void setProjection(void *obj, int type, double rot, int hemisphere, double stdLat);
static void getProjection(void *obj, int type, double *rot, int *hemisphere, double *stdLat);
static void readXYProjInfoGDAL(char *xyFile, void *obj,  int type);
static void getXYSize(void *xyImage, int32_t type, int *xSize, int *ySize);
static void setGeom(void *obj, int type, double x0, double y0, double deltaX, double deltaY, int32_t xSize, int32_t ySize);
static void getGeom(void *obj, int type, double *x0, double *y0, double *deltaX, double *deltaY, int32_t *xSize, int32_t *ySize);
static int32_t getEPSG(GDALDatasetH hDataset);
static void getImageSize(GDALDatasetH hDataset, int32_t *xSize, int32_t *ySize);
static void getGeoTransform(GDALDatasetH hDataset, double *geoTransform);
static double parseWKT(OGRSpatialReferenceH hSRS, const char *projection, const char *param);
static char *strcasestr1(const char *haystack, const char *needle);

void readXYDEMGeoInfo(char *xyFile, xyDEM *xydem, int32_t resetProjection)
{
	// Read the geometric information for a DEM
	readXYGeoInfo(xyFile, xydem, resetProjection, DEM);
}

void readXYDEM(char *xyFile, xyDEM *xydem) {
	// Read a full XY DEM
	fprintf(stderr, "READING DEM\n");
	// All zeros indicates no cropping.
	readXYDEMcrop(xyFile, xydem, 0., 0., 0., 0.);
}

void readXYDEMcrop(char *xyFile, xyDEM *xyDEM, double xmin, double xmax, double ymin, double ymax)
{   // Read a cropped area from a DEM
	xyDEM->z = readXYImageCrop(xyFile, xyDEM,  xmin,  xmax, ymin, ymax, DEM);
	if(xyDEM->z == NULL)
		error("Could readXYDEMcrop %s\n", xyFile);
}

void readXYVel(xyVEL *xyvel, char *velFile)
{	// read a full velocity map
	readXYCropVel(xyvel, velFile, 0., 0., 0., 0. );
}

void readXYCropVel(xyVEL *xyvel, char *velFile, double xmin, double xmax, double ymin, double ymax)
{   // read a cropped velocity map
	char *vxFile, *vyFile;
	if(has_extension(velFile, ".tif") == TRUE || has_extension(velFile, ".vrt")) {
		//Option 1 use a wild card *
		vxFile = replace_wildcard(velFile, "*", "vx");
		vyFile = replace_wildcard(velFile, "*", "vy");
		// Option to use vv in place of vx, vy
		vxFile = replace_wildcard(velFile, "vv", "vx");
		vyFile = replace_wildcard(velFile, "vv", "vy");
		if(vxFile == NULL || vyFile == NULL) 
			error("velocity file name %s has not '*' wildcard for vx and vy location");
	} else {
		vxFile = appendSuffix(velFile, ".vx",
							   (char *)malloc((size_t)strlen(velFile) + 4));
		vyFile = appendSuffix(velFile, ".vx",
							   (char *)malloc((size_t)strlen(velFile) + 4));
	}
	// Check file exists. TRUE forces abort if file does not exist
	fprintf(stderr, "Files: %s %s\n", vxFile, vyFile);
	error("STOP");
	fileExists(vxFile, TRUE);
	fileExists(vyFile, TRUE);
	
	xyvel->vx = readXYImageCrop(vxFile, xyvel,  xmin,  xmax, ymin, ymax, VELOCITY);
	xyvel->vy = readXYImageCrop(vyFile, xyvel,  xmin,  xmax, ymin, ymax, VELOCITY);
}

static void readXYGeodatFile(char *xyFile, double *x0, double *y0, double *deltaX, double *deltaY, 
							int32_t *xSize, int32_t *ySize, double *rot, int32_t *hemisphere, double *stdLat)
{
	extern int32_t HemiSphere;
	extern double Rotation;
	extern double SLat;
	int32_t lineCount = 0, eod;
	float dum1, dum2;
	char line[256];
	// Defaults
	*rot = Rotation;
	*hemisphere = HemiSphere;
	*stdLat = SLat;
	char *geodatFile = appendSuffix(xyFile, ".geodat",
							(char *)malloc((size_t)strlen(xyFile) + 8));
	fprintf(stderr, "Geodat %s\n", geodatFile);
	FILE *fp = openInputFile(geodatFile);
	if (fp == NULL)
		error("*** getXYDEM: Error opening %s ***\n", xyFile);
	// Read parameters
	lineCount = getDataString(fp, lineCount, line, &eod); /* Skip # 2 line */
	lineCount = getDataString(fp, lineCount, line, &eod);
	sscanf(line, "%f %f\n", &dum1, &dum2);
	*xSize = (int)dum1;
	*ySize = (int)dum2;
	lineCount = getDataString(fp, lineCount, line, &eod);
	fprintf(stderr, "%s\n", line);
	sscanf(line, "%lf %lf\n", deltaX, deltaY);
	*deltaX *= MTOKM;
	*deltaY *= MTOKM;
	lineCount = getDataString(fp, lineCount, line, &eod);
	sscanf(line, "%lf %lf\n", x0, y0);
	fprintf(stderr, "%s\n", line);
	lineCount = getDataString(fp, lineCount, line, &eod);
	fprintf(stderr, "%s\n", line);
	/* Hemisphere */
	if (eod != TRUE)
	{
		fprintf(stderr, "Hemispher\n");
		if (strstr(line, "SOUTH"))
			HemiSphere = SOUTH;
		else
			HemiSphere = NORTH;
		lineCount = getDataString(fp, lineCount, line, &eod);
		*hemisphere = HemiSphere;
	}
	/* Rotation */
	if (eod != TRUE)
	{
		fprintf(stderr, "Rot\n");
		sscanf(line, "%lf", rot);
		lineCount = getDataString(fp, lineCount, line, &eod);
	}
	/* Std lat */
	if (eod != TRUE)
	{
		fprintf(stderr, "Stdlat\n");
		sscanf(line, "%lf", stdLat);
	}
	fclose(fp);
}


static float **readXYImageCrop(char *xyFile, void *xyImage, double xmin, double xmax, double  ymin, double ymax, int type)
{  // Read a cropped area from an XY image binary, tiff, or VRT. 
	extern int32_t HemiSphere;
	extern double Rotation;
	extern double SLat;
	double x0DEM, y0DEM, deltaX, deltaY;
	double x0Crop, y0Crop;
	double stdLat, rot;
	int32_t hemisphere;
	int32_t xSizeDEM, ySizeDEM, xSizeCrop, ySizeCrop;
	float **image;
	int32_t i, j;
	double xmaxDEM, ymaxDEM, xminDEM, yminDEM;
	int32_t x1, x2, y1, y2;
	off_t offset1;
	fprintf(stderr, "READING *** cropped *** %s\n", xyFile);
	if (has_extension(xyFile, ".tif") == TRUE || has_extension(xyFile, ".vrt")) {
		return readXYImageGDALCropped(xyFile, xyImage,  xmin, xmax, ymin, ymax, type);
	}
	// Read geometric info
	readXYGeodatFile(xyFile, &x0DEM, &y0DEM, &deltaX, &deltaY, &xSizeDEM, &ySizeDEM, &rot, &hemisphere, &stdLat);
	//
	//  Now compute parameters for cropped dem
	xminDEM = x0DEM; /* lower left corner defines min extent of input dem */
	yminDEM = y0DEM;
	xmaxDEM = xminDEM + deltaX * (float)(xSizeDEM - 1); /* upper right corner defines max extent */
	ymaxDEM = yminDEM + deltaY * (float)(ySizeDEM - 1);
	// If all zero, no cropping so set to dem extent
	fprintf(stderr, "input dem limits %f %f %f %f\n", xminDEM, xmaxDEM, yminDEM, ymaxDEM);
	fprintf(stderr, "requested dem limits %f %f %f %f\n", xmin, xmax, ymin, ymax);
    if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0) {
		xmin = xminDEM;
		xmax = xmaxDEM;
		ymin = yminDEM;
		ymax = ymaxDEM;
	} else {
		xmin = max(xminDEM, xmin);
		xmax = min(xmaxDEM, xmax);
		ymin = max(yminDEM, ymin);
		ymax = min(ymaxDEM, ymax);
	}
	fprintf(stderr, "output dem limits %f %f %f %f\n", xmin, xmax, ymin, ymax);
	xSizeCrop = (int)((xmax - xmin)/deltaX) + 1;
	x1 = (int) ((xmin - xminDEM)/ deltaX);
	x2 = x1 + xSizeCrop - 1;
	ySizeCrop = (int)((ymax - ymin)/deltaY) + 1;
	y1 = (int) ((ymin - yminDEM)/ deltaY);
	y2 = y1 + ySizeCrop - 1;
	x0Crop = x0DEM + x1 * deltaX;
	y0Crop = y0DEM + y1 * deltaY;

	fprintf(stderr, "pixel limits %i %i %i %i\n", x1, x2, y1, y2);
	fprintf(stderr, "Dem Size %i %i\n", xSizeCrop, ySizeCrop);
	fprintf(stderr, "Dem Res %f %f\n", deltaX, deltaY);
	fprintf(stderr, "Dem Origin %f %f\n", x0Crop, y0Crop);
	//
	// Open image file file
	FILE *fp = fopen(xyFile, "r");
	if (fp == NULL)
		error("*** getXYDEM: Error opening %s ***\n", xyFile);
	// Read Image file
	image = (float **)mallocImage(ySizeCrop, xSizeCrop);
	for (i = y1, j = 0; i <= y2; i++, j++)
	{
		offset1 = ((off_t) xSizeDEM * (off_t)i + (off_t)x1) * (off_t)sizeof(float);
		fseeko(fp, offset1, SEEK_SET);
		freadBS(image[j], sizeof(float), xSizeCrop, fp, FLOAT32FLAG);
	}
	setProjection(xyImage, type, rot, hemisphere, stdLat);
	setGeom(xyImage, type, x0Crop, y0Crop, deltaX, deltaY, xSizeCrop, ySizeCrop);
	fclose(fp);
	return image;
}

static float **readXYImageGDALCropped(char *xyFile, void *xyImage, double xMin, double xMax,
						   		double yMin, double yMax, int32_t type)
// Read a cropped area from the DEM.
{
	extern int32_t HemiSphere;
	int xOffset, yOffset;
	int32_t xSize, ySize;
	float **image;
	int32_t i, j;
	int32_t flip;
	// Projection for the DEM
	readXYProjInfoGDAL(xyFile, xyImage, type);
	// Now open
	GDALDatasetH hDataset = GDALOpen((const char *)xyFile, GA_ReadOnly);
	// Get the crop bounds
	getCropBounds(hDataset, xMin, xMax, yMin, yMax, &xOffset, &yOffset, xyImage, &flip, type);
	getXYSize(xyImage, type, &xSize, &ySize);
	// Get the band (for now assume single band)
	GDALRasterBandH hBand = GDALGetRasterBand(hDataset, 1);
	if (hBand == NULL)
	{
		fprintf(stderr, "Failed to get raster band.\n");
		GDALClose(hDataset);
		return NULL;
	}
	// Create a buffer to hold the cropped data
	fprintf(stderr, "xSize, ySize %i %i\n", xSize, ySize);
	float *buffer = (float *)CPLMalloc(sizeof(float) * xSize * ySize);
	// Set pointers to rows
	image = (float **)malloc((size_t)(ySize * sizeof(float *)));
	// Flip for -dy
	for (i = 0, j = ySize-1; i < ySize; i++, j--)
		if (flip == TRUE)
			image[i] = &(buffer[j * xSize]);
		else
			image[i] = &(buffer[i * xSize]);
	// Read the cropped region
	CPLErr err = GDALRasterIO(hBand, GF_Read, xOffset, yOffset, xSize, ySize,
							  buffer, xSize, ySize, GDT_Float32, 0, 0);
	return image;
}

static void getCropBounds(GDALDatasetH hDataset, double xMin, double xMax, double yMin, double yMax,
				   int *xOffset, int *yOffset, void *obj, int32_t *flip, int32_t type)
// Get the bounds for the DEM and set the origin and size for the DEM. Crops to edge of dem if requested
// area extends off an edge.
{
	double geoTransform[6];
	int xSizeDEM, ySizeDEM;
	double xMinDEM, xMaxDEM, yMinDEM, yMaxDEM, dx, dy;
	int32_t xSize, ySize;
	double x0, y0, deltaX, deltaY;
	int crop = TRUE;
	if (xMin == 0 && xMax == 0 && yMin == 0 && yMax == 0)
		crop = FALSE;
	//
	getImageSize(hDataset, &xSizeDEM, &ySizeDEM);
	getGeoTransform(hDataset, geoTransform);
	// Working internally in km
	for(int i=0; i < 6; i++) geoTransform[i] *= MTOKM;
	dx = geoTransform[1];
	dy = geoTransform[5];
	// X Dem limits changed to pixel centers
	xMinDEM = geoTransform[0] + dx * 0.5;
	xMaxDEM = xMinDEM + dx * (xSizeDEM - 1);
	if (crop == FALSE)
	{
		xMin = xMinDEM;
		xMax = xMaxDEM;
	}
	// Crop to boundaries if needed
	xMin = max(xMin, xMinDEM);
	xMax = min(xMax, xMaxDEM);
	*xOffset = (int) ((xMin - xMinDEM) / dx);
	xSize = (int) ((xMax - xMin) / dx) + 1;
	x0 = (xMinDEM + (*xOffset) * dx);
	//
	// Handle up or lower origin as indicated by sign of dy
	yMinDEM = geoTransform[3] + dy * 0.5;
	yMaxDEM = yMinDEM + dy * (ySizeDEM - 1);
	if (dy < 0)
	{
		// Swap yMinDEM and yMaxDEM for negative dy
		double temp = yMinDEM;
		yMinDEM = yMaxDEM;
		yMaxDEM = temp;
		*flip = TRUE;
	}
	else
	{
		*flip = FALSE;
	}
	// Set yMin and yMax if cropping is disabled
	if (crop == FALSE)
	{
		yMin = yMinDEM;
		yMax = yMaxDEM;
	}
	// Crop to boundaries if needed
	yMin = max(yMin, yMinDEM);
	yMax = min(yMax, yMaxDEM);
	// Compute yOffset, ySize, and y0
	ySize = (int) ((yMax - yMin) / fabs(dy)) + 1;
	*yOffset = (int)((yMin - yMinDEM) / fabs(dy));
	if(*flip == TRUE) *yOffset = (ySizeDEM -1) - (*yOffset + ySize  -1);
	y0 = (yMaxDEM + (*yOffset) * dy + (dy < 0 ? (ySize - 1) * dy : 0)) ;
	// Check area requested fits in DEM
	if ((xMin > xMaxDEM || xMax < xMinDEM || yMin > yMaxDEM || yMax < yMinDEM) && crop == TRUE)
	{
		fprintf(stderr, "%i %i\n", xSize, ySize);
		fprintf(stderr, "%f %f %f %f\n", xMin, xMax, yMin, yMax);
		fprintf(stderr, "%f %f %f %f\n", xMinDEM, xMaxDEM, yMinDEM, yMaxDEM);
		error("Area requested does not fit in DEM");
	}
	deltaX = fabs(dx);
	deltaY = fabs(dy);
	fprintf(stderr, "xMin, xMax, yMin, yMax %f %f %f %f\n", xMin, xMax, yMin, yMax);
	fprintf(stderr, "xMinDEM, xMaxDEM, yMinDEM, yMaxDEM %f %f %f %f\n", xMinDEM, xMaxDEM, yMinDEM, yMaxDEM);
	fprintf(stderr, "xOffset, yOffset, xSize, ySize %i %i %i %i\n", *xOffset, *yOffset, xSize, ySize);
	// Now set results in vel or dem structure
	setGeom(obj, type, x0, y0, deltaX, deltaY, xSize, ySize);
}

static void readXYGeoInfo(char *xyFile, void *xyImage, int32_t resetProjection, int type)
{ // Get the geometry info for and xyDEM or xyVel for binary, tiff, or vrt image.
	extern int32_t HemiSphere;
	extern double Rotation;
	extern double SLat;
	double x0, y0, deltaX, deltaY;
	int32_t xSize, ySize;
	double rot, stdLat;
	int32_t hemisphere;
	
	if (has_extension(xyFile, ".tif") == TRUE || has_extension(xyFile, ".vrt"))
	{	// Read from a tiff or vrt
		readXYGeoInfoGDAL(xyFile, xyImage, type);
		// Geth parameters from the xyImage
		getGeom(xyImage, type, &x0, &y0, &deltaX, &deltaY, &xSize, &ySize);
		getProjection(xyImage, type, &rot, &hemisphere, &stdLat);
	} 
	else 
	{	
		readXYGeodatFile(xyFile, &x0, &y0, &deltaX,&deltaY, &xSize, &ySize, &rot, &hemisphere, &stdLat);
		setProjection(xyImage, type, rot, hemisphere, stdLat);
		setGeom(xyImage, type, x0, y0, deltaX, deltaY, xSize, ySize);
	}
	HemiSphere = hemisphere;
	// Print projection summary
	fprintf(stderr, "Dem Size %i %i\n", xSize, ySize);
	fprintf(stderr, "Dem Res %f %f\n", deltaX, deltaY);
	fprintf(stderr, "Dem Origin %f %f\n", x0, y0);
	if (HemiSphere == SOUTH)
		fprintf(stderr, "SOUTHER HEMISPHERE PROJECTION\n");
	else
		fprintf(stderr, "NORTHERN HEMISPHERE PROJECTION\n");
	fprintf(stderr, "Rotation = %f\n", rot);
	fprintf(stderr, "Standard Lat = %f\n", stdLat);
	// Force projection to match dem
	if (resetProjection == TRUE)
	{
		Rotation = rot;
		SLat = stdLat;
	}
}

void readXYGeoInfoGDAL(char *xyFile, void *xyImage, int32_t type)
{// Read the geometric and projection info from a tiff or vrt and save them for xyDEM or xyVel (type=DEM or VELOCITY)
	int xOffset, yOffset;
	int32_t flip;
	// Projection for the DEM
	readXYProjInfoGDAL(xyFile, xyImage, type);
	// Now open
	GDALDatasetH hDataset = GDALOpen((const char *)xyFile, GA_ReadOnly);
	// Get the crop bounds
	getCropBounds(hDataset, 0, 0, 0, 0, &xOffset, &yOffset, xyImage, &flip, type);
	// Get the band (for now assume single band)
	GDALClose(hDataset);
}

static void getProjection(void *obj, int type, double *rot, int *hemisphere, double *stdLat) 
{// Get the projection parameters for an xyDEM or xyVel as indicated by type with DEM or VELOCITY
	if(type == DEM) {
		xyDEM *xydem = (xyDEM *) obj;
		*rot = xydem->rot;
		*hemisphere = xydem->hemisphere;
		*stdLat = xydem->stdLat;
	} else if(type == VELOCITY)
	{
		xyVEL *xyvel = (xyVEL *) obj;
		*rot = xyvel->rot;
		*hemisphere = xyvel->hemisphere;
		*stdLat = xyvel->stdLat;
	}
}

static void setProjection(void *obj, int type, double rot, int hemisphere, double stdLat) 
{ // Set the projection parameters for an xyDEM or xyVel as indicated by type with DEM or VELOCITY
	if(type == DEM) {
		xyDEM *xydem = (xyDEM *) obj;
		xydem->rot = rot;
		xydem->hemisphere = hemisphere;
		xydem->stdLat = stdLat;
	} else if(type == VELOCITY)
	{
		xyVEL *xyvel = (xyVEL *) obj;
		xyvel->rot = rot;
		xyvel->hemisphere = hemisphere;
		xyvel->stdLat = stdLat;
	}
}

static char *strcasestr1(const char *haystack, const char *needle)
 {
	// because some systems don't have strcasestr
    size_t needle_len = strlen(needle);
    if (needle_len == 0) return (char *)haystack;

    for (const char *p = haystack; *p; ++p) {
        if (tolower((unsigned char)*p) == tolower((unsigned char)*needle)) {
            if (strncasecmp(p, needle, needle_len) == 0) {
                return (char *)p;
            }
        }
    }
    return NULL;
}

static double parseWKT(OGRSpatialReferenceH hSRS, const char *projection, const char *param)
{
	double paramValue;
	char buf[1024];
	int success = 0;
    paramValue = OSRGetProjParm(hSRS, param, 0.0, &success);
	if(success == 1) 
	{
		return paramValue;
	}
    // Backup parse from wkt string directly
	//printf("%s", projection);
	buf[0] = '\0';
	sprintf(buf, "PARAMETER[\"%s\",", param);
    char *paramString = strcasestr1(projection, buf);
	if (paramString != NULL) {
		buf[0] = '\0';
		sprintf(buf, "PARAMETER[\"%s\",%%lf", param);
		sscanf(paramString , buf, &paramValue);
		printf("%s (manual parsing): %f\n", param, paramValue);
	} else {
		printf("%s parameter not found in WKT.\n", param);
		return (double) -LARGEINT;
	}
	return paramValue;
}

static void getPolarStereoParams(GDALDatasetH hDataset, double *standardLat, double *rot)
{

	const char *projectionConst = GDALGetProjectionRef(hDataset);
    if (projectionConst == NULL || strlen(projectionConst) == 0) {
        error("No projection found in the file.\n");
        return;
    }
	char *projection = strdup(projectionConst);
    OGRSpatialReferenceH hSRS = OSRNewSpatialReference(NULL);
    if (OSRImportFromWkt(hSRS, (char **)&projection) != OGRERR_NONE) {
        error("Failed to import spatial reference from WKT.\n");
    }
	*standardLat = parseWKT(hSRS, projectionConst, "latitude_of_origin");
	if(*standardLat < -90) {
		*standardLat = parseWKT(hSRS, projectionConst, "latitude_of_standard_parallel");
	}
	*rot = -1 * parseWKT(hSRS, projectionConst, "central_meridian");
	
	fprintf(stderr, "Rotation %lf\n", *rot);
	fprintf(stderr, "Standard lat %lf\n", *standardLat);
	if(*rot < -180. || *standardLat <= -90) 
	{
		error("Could not parse Rotation and/or Standard lat");
	}
    OSRDestroySpatialReference(hSRS);
}

void readXYProjInfoGDAL(char *xyFile, void *obj,  int type)
// Read and set the projection from a tiff or DEM parameters for with type as  DEM or VELOCITY
{	
	extern double Rotation;
	extern double SLat;
	double standardLat, rot;
	//
	GDALDatasetH hDataset = GDALOpen((const char *)xyFile, GA_ReadOnly);
	// Get epsg
	int32_t epsg = getEPSG(hDataset);
	fprintf(stderr, "EPSG %i\n", epsg);
	// Set project parameters for vrt.
	if(epsg == 3413)
	{
		setProjection(obj, type, 45.0, NORTH, 70.0);
	}
	 else if(epsg == 3031)
	{
		setProjection(obj, type, 0.0, SOUTH, 71.0);
	}
	else
	{
		getPolarStereoParams(hDataset, &standardLat, &rot);
		if(standardLat < 0)
		{
			setProjection(obj, type, rot, SOUTH, standardLat);
		} else 
		{
			setProjection(obj, type, rot, NORTH, standardLat);
		}
	}
	GDALClose(hDataset);
}

static void getXYSize(void *xyImage, int32_t type, int *xSize, int *ySize)
{// set size info for an xyDEM or xyVel structure specified with type as DEM or VELOCITY
	if(type == DEM) 
	{
		*xSize = ((xyDEM *)xyImage)->xSize;	
		*ySize = ((xyDEM *)xyImage)->ySize;	
	}
	else if(type == VELOCITY)
	{
		*xSize = ((xyVEL *)xyImage)->xSize;
		*ySize = ((xyVEL *)xyImage)->ySize;
	}		
	else
			error("readXYGeoInfo invalid type %i", type);	
}

static void setGeom(void *obj, int type, double x0, double y0, double deltaX, double deltaY, int32_t xSize, int32_t ySize) 
{// set geometry params from xyDEM or xyVel structure specified with type as DEM or VELOCITY
	if(type == DEM) {
		xyDEM *xydem= (xyDEM *) obj;
		xydem->x0 = x0;
		xydem->y0 = y0;
		xydem->deltaX = deltaX; 
		xydem->deltaY = deltaY;
		xydem->xSize = xSize;
		xydem->ySize = ySize;
	} else if(type == VELOCITY)
	{
		xyVEL *xyvel = (xyVEL *) obj;
		xyvel->x0 = x0;
		xyvel->y0 = y0;
		xyvel->deltaX = deltaX; 
		xyvel->deltaY = deltaY;
		xyvel->xSize = xSize;
		xyvel->ySize = ySize;
	}
}

static void getGeom(void *obj, int type, double *x0, double *y0, double *deltaX, 
					double *deltaY, int32_t *xSize, int32_t *ySize) 
// get geometry params from xyDEM or xyVel structure specified with type as DEM or VELOCITY
{
	if(type == DEM) {
		xyDEM *xydem= (xyDEM *) obj;
		*x0 = xydem->x0;
		*y0 = xydem->y0;
		*deltaX = xydem->deltaX; 
		*deltaY = xydem->deltaY;
		*xSize = xydem->xSize;
		*ySize = xydem->ySize;
	} else if(type == VELOCITY)
	{
		xyVEL *xyvel = (xyVEL *) obj;
		*x0 = xyvel->x0;
		*y0 = xyvel->y0;
		*deltaX = xyvel->deltaX; 
		*deltaY = xyvel->deltaY;
		*xSize = xyvel->xSize;
		*ySize = xyvel->ySize;
	}
}



static int32_t getEPSG(GDALDatasetH hDataset)
{	// Get epsg from a tiff or vrt
	int32_t epsg;
	// Get the SRS
	const char *projectionConst = GDALGetProjectionRef(hDataset);
	char *projection = strdup(projectionConst);
	OGRSpatialReferenceH hSRS = OSRNewSpatialReference(NULL);

	if (projection != NULL && strlen(projection) > 0)
	{
		// Import the spatial reference from the dataset's projection
		if (OSRImportFromWkt(hSRS, &projection) != OGRERR_NONE)
		{
			fprintf(stderr, "Failed to import spatial reference.\n");
			return 0;
		}
	}
	else
	{
		printf("No projection found in the file.\n");
	}
	// Attempt to identify the EPSG code
	const char *epsgCode = OSRGetAuthorityCode(hSRS, NULL);
	if (epsgCode != NULL)
	{
		printf("EPSG Code: %s\n", epsgCode);
		epsg = atoi(epsgCode); // Convert EPSG code to integer
		printf("EPSG Code as integer: %d\n", epsg);
	}
	else
	{
		printf("Failed to determine the EPSG code.\n");
	}
	OSRDestroySpatialReference(hSRS);
	fprintf(stderr, "Returing epsg %i", epsg);
	return epsg;
}

static void getImageSize(GDALDatasetH hDataset, int32_t *xSize, int32_t *ySize)
// Get tiff image size from a tiff or vrt
{
	*xSize = GDALGetRasterXSize(hDataset); // Get width (number of pixels in X direction)
	*ySize = GDALGetRasterYSize(hDataset); // Get height (number of pixels in Y direction)
	fprintf(stderr, "xSize, ySize %i %i\n", *xSize, *ySize);
}


static void getGeoTransform(GDALDatasetH hDataset, double *geoTransform)
// Get geotransform from a tiff or vrt
{
	if (GDALGetGeoTransform(hDataset, geoTransform) == CE_None)
	{
		printf("GeoTransform:\n");
		printf("  Origin : (%.6f, %.6f)\n", geoTransform[0], geoTransform[3]);
		printf("  Pixel Size: (%.6f, %.6f)\n", geoTransform[1], geoTransform[5]);
		printf("  Rotation: X: %.6f, Y: %.6f\n", geoTransform[2], geoTransform[4]);
	}
	else
	{
		printf("GeoTransform not available for this dataset.\n");
	}
}
