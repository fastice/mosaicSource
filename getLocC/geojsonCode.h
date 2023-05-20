#include <ogr_core.h>
#include <ogr_srs_api.h>
#include <cpl_conv.h>
#include <gdal.h>
#include "mosaicSource/common/common.h"


OGRDataSourceH getGeojsonDataSet(char *geojsonFile, int registerGDAL);

OGRGeometryH createGeometry(double *lat, double *lon);
OGRFeatureDefnH createFeatureDef(int32_t nState);
const char *svTag(int32_t i, char *svType);