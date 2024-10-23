#include <stdlib.h>
#include <string.h>
//#include "mosaicSource/common/geojsonCode.h"
#include "mosaicSource/common/common.h"

OGRDataSourceH getGeojsonDataSet(char *geojsonFile) {
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
    // Removed SRS to avoid problems with proj lib
    // Create an SRS for the geometry
    //OGRSpatialReferenceH mySRS = OSRNewSpatialReference("");
    //OGRErr srsErr = OSRImportFromEPSG(mySRS, 4326);
    //   Create the geometry
    OGRGeometryH myGeom = OGR_G_CreateGeometry(wkbPolygon);
    //OGR_G_AssignSpatialReference(myGeom, mySRS);
    OGRGeometryH myRing = OGR_G_CreateGeometry(wkbLinearRing);
    OGR_G_SetPoints(myRing, 4, lat, sizeof(double), lon, sizeof(double), NULL, sizeof(double));
    OGR_G_CloseRings(myRing);
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
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("RangeErrorCorrection", OFTReal));
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
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("CorrectedTime", OFTString));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("Wavelength", OFTReal));
    // Single look 
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("SLCRangePixelSize", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("SLCAzimuthPixelSize", OFTReal));
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("deltaT", OFTReal));
     // Byte Order
    OGR_FD_AddFieldDefn(myFeatureDef, OGR_Fld_Create("ByteOrder", OFTString));
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

