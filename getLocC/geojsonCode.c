#include <stdlib.h>
#include <string.h>
#include "mosaicSource/getLocC/geojsonCode.h"


OGRDataSourceH getGeojsonDataSet(char *geojsonFile, int registerGDAL) {

    // Register gdal if reqested

    if(registerGDAL == TRUE)
        fprintf(stderr, "Registering GDAL\n");
        GDALAllRegister();
    // Get the geosjson driver
    GDALDriverH  myDriver = GDALGetDriverByName("GeoJSONSeq");
    // Remove any existing file since driver won't overwrite
    if (access(geojsonFile, F_OK) == 0)
	{
		remove(geojsonFile);
	}
    GDALDatasetH myData = GDALCreate(myDriver, geojsonFile,  0, 0, 0, GDT_Unknown, NULL);
    return myData;
}

OGRGeometryH createGeometry(double *lat, double *lon) {
    // Create an SRS for the geometry
    OGRSpatialReferenceH mySRS = OSRNewSpatialReference("");
    OGRErr srsErr = OSRImportFromEPSG(mySRS, 4326);
    //   Create the geometry
    OGRGeometryH myGeom = OGR_G_CreateGeometry(wkbPolygon);
    OGR_G_AssignSpatialReference(myGeom, mySRS);
    OGRGeometryH myRing = OGR_G_CreateGeometry(wkbLinearRing);
    OGR_G_SetPoints(myRing, 4, lat, sizeof(double), lon, sizeof(double), NULL, sizeof(double));
    OGR_G_AddGeometryDirectly(myGeom, myRing);
    return myGeom;
}

const char *svTag(int32_t i, char *svType) {
    // Malloc space for stateVectorPositionNNN + \0
    char *myTagSpace = calloc(50, sizeof(char));
    sprintf(myTagSpace, "SV_%s_%i", svType, i);
    const char *myTag = myTagSpace;
    return  myTag;
}


OGRFeatureDefnH createFeatureDef(int32_t nState)
{
    int i;
    OGRFeatureDefnH myFeatureDef = OGR_FD_Create("properties");
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("ImageName", OFTString));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("Date", OFTDate));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("NominalTime", OFTTime));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("TrackDirection", OFTReal));
    //
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("NumberRangeLooks", OFTInteger));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("NumberAzimuthLooks", OFTInteger));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("MLRangeSize", OFTInteger));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("MLAzimuthSize", OFTInteger));
    //
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("VelocityAlongTrack", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("PRF", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("MLNearRange", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("MLCenterRange", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("MLFarRange", OFTReal));
   
    //
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("LookDirection", OFTString));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("PassType", OFTString));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("CenterLatLon", OFTRealList));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("TimeToFirstSLCSample", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("SkewOffset", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("Squint", OFTReal));
    //
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("EarthRadiusMajor", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("EarthRadiusMinor", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("MLIncidenceCenter", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("SpaceCraftAltitude", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("CorrectedTime", OFTTime));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("Wavelength", OFTReal));
    // Single look 
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("RangeSLCPixelSize", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("AzimuthSLCPixelSize", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("deltaT", OFTReal));
    // State vectors
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("NumberOfStateVectors", OFTInteger));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("TimeOfFirstStateVector", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("StateVectorInterval", OFTReal));

    for(i=1; i <= nState; i++) {
        OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create(svTag(i, "Pos"), OFTRealList));
        OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create(svTag(i, "Vel"), OFTRealList));
    }
   
    return myFeatureDef;
}
/*
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

*/
