#include <unistd.h>
#include <stdio.h>
#include "mosaicSource/common/writeTieResidualsGpkg.h"

/*
  Write tie points + per-point fit residuals as a GeoPackage point layer for
  -debug output on tiepoints/rparams/azparams. Modeled on the OGR vector-writing
  pattern already established by common/geojsonCode.c + getLocC/getlocc.c, adapted
  for the GPKG driver (a strict-schema format -- fields must be added to the real
  layer via OGR_L_CreateField() before any feature is written, unlike the
  GeoJSONSeq driver getlocc.c targets).
*/

OGRDataSourceH openTieResidualsGpkg(const char *gpkgPath, int32_t append)
{
	GDALDatasetH ds;
	if (!append)
	{
		GDALDriverH driver = GDALGetDriverByName("GPKG");
		if (driver == NULL)
			error("openTieResidualsGpkg: GPKG driver not available");
		if (access(gpkgPath, F_OK) == 0)
			remove(gpkgPath);
		ds = GDALCreate(driver, gpkgPath, 0, 0, 0, GDT_Unknown, NULL);
	}
	else
	{
		ds = GDALOpenEx(gpkgPath, GDAL_OF_VECTOR | GDAL_OF_UPDATE, NULL, NULL, NULL);
	}
	if (ds == NULL)
		error("openTieResidualsGpkg: could not open/create %s", gpkgPath);
	return (OGRDataSourceH)ds;
}

int32_t writeTieResidualsLayer(OGRDataSourceH ds, const char *layerName,
								tiePointsStructure *tiePoints, tieResidualsType *fit,
								const char *residualFieldName, const char *residualDescription,
								double rotation, double stdLat, int32_t hemisphere)
{
	const char *epsg = getEPSGFromProjectionParams(rotation, stdLat, hemisphere);
	OGRSpatialReferenceH srs = OSRNewSpatialReference(NULL);
	OSRImportFromEPSG(srs, atoi(epsg));

	OGRLayerH layer = OGR_DS_CreateLayer((GDALDatasetH)ds, layerName, srs, wkbPoint, NULL);
	OSRDestroySpatialReference(srs);
	if (layer == NULL)
		error("writeTieResidualsLayer: could not create layer %s", layerName);

	OGR_L_CreateField(layer, OGR_Fld_Create("id", OFTInteger), TRUE);
	OGR_L_CreateField(layer, OGR_Fld_Create("lat", OFTReal), TRUE);
	OGR_L_CreateField(layer, OGR_Fld_Create("lon", OFTReal), TRUE);
	OGR_L_CreateField(layer, OGR_Fld_Create("x_km", OFTReal), TRUE);
	OGR_L_CreateField(layer, OGR_Fld_Create("y_km", OFTReal), TRUE);
	OGR_L_CreateField(layer, OGR_Fld_Create("range", OFTReal), TRUE);
	OGR_L_CreateField(layer, OGR_Fld_Create("azimuth", OFTReal), TRUE);
	OGR_L_CreateField(layer, OGR_Fld_Create("z", OFTReal), TRUE);
	OGR_L_CreateField(layer, OGR_Fld_Create("weight", OFTReal), TRUE);
	OGRFieldDefnH residField = OGR_Fld_Create(residualFieldName, OFTReal);
	OGR_Fld_SetAlternativeName(residField, residualDescription);
	OGR_L_CreateField(layer, residField, TRUE);

	OGRFeatureDefnH featureDefn = OGR_L_GetLayerDefn(layer);
	int32_t idIdx = OGR_FD_GetFieldIndex(featureDefn, "id");
	int32_t latIdx = OGR_FD_GetFieldIndex(featureDefn, "lat");
	int32_t lonIdx = OGR_FD_GetFieldIndex(featureDefn, "lon");
	int32_t xkmIdx = OGR_FD_GetFieldIndex(featureDefn, "x_km");
	int32_t ykmIdx = OGR_FD_GetFieldIndex(featureDefn, "y_km");
	int32_t rangeIdx = OGR_FD_GetFieldIndex(featureDefn, "range");
	int32_t azIdx = OGR_FD_GetFieldIndex(featureDefn, "azimuth");
	int32_t zIdx = OGR_FD_GetFieldIndex(featureDefn, "z");
	int32_t weightIdx = OGR_FD_GetFieldIndex(featureDefn, "weight");
	int32_t residIdx = OGR_FD_GetFieldIndex(featureDefn, residualFieldName);

	for (int32_t k = 0; k < fit->npts; k++)
	{
		int32_t i = fit->origIndex[k];
		double lat = tiePoints->lat[i], lon = tiePoints->lon[i];
		double xkm, ykm;
		/* True polar-stereo coords -- NOT tiePoints->x[i]/y[i], which hold a scaled
		   azimuth fit-coordinate (common/computeTiePoints.c), not map coords. */
		lltoxy1(lat, lon, &xkm, &ykm, rotation, stdLat);

		OGRFeatureH feature = OGR_F_Create(featureDefn);
		OGR_F_SetFieldInteger(feature, idIdx, i);
		OGR_F_SetFieldDouble(feature, latIdx, lat);
		OGR_F_SetFieldDouble(feature, lonIdx, lon);
		OGR_F_SetFieldDouble(feature, xkmIdx, xkm);
		OGR_F_SetFieldDouble(feature, ykmIdx, ykm);
		OGR_F_SetFieldDouble(feature, rangeIdx, tiePoints->r[i]);
		OGR_F_SetFieldDouble(feature, azIdx, tiePoints->a[i]);
		OGR_F_SetFieldDouble(feature, zIdx, tiePoints->z[i]);
		OGR_F_SetFieldDouble(feature, weightIdx, tiePoints->weight[i]);
		OGR_F_SetFieldDouble(feature, residIdx, fit->residual[k]);

		OGRGeometryH geom = OGR_G_CreateGeometry(wkbPoint);
		/* Geometry MUST be in the layer SRS's units (metres for EPSG:3413/3031) --
		   NOT lat/lon degrees, NOT the km values used for the x_km/y_km attributes. */
		OGR_G_SetPoint_2D(geom, 0, xkm * 1000.0, ykm * 1000.0);
		OGR_F_SetGeometryDirectly(feature, geom);

		OGR_L_CreateFeature(layer, feature);
		OGR_F_Destroy(feature);
	}
	return TRUE;
}

void closeTieResidualsGpkg(OGRDataSourceH ds)
{
	GDALClose((GDALDatasetH)ds);
}
