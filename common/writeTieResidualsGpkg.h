#ifndef WRITETIERESIDUALSGPKG_H
#define WRITETIERESIDUALSGPKG_H
#include <ogr_core.h>
#include <ogr_srs_api.h>
#include <cpl_conv.h>
#include <gdal.h>
#include "mosaicSource/common/common.h"

/* Per-point residuals recomputed after a winning svdfit() call. origIndex[k] is the
   0-based index into tiePointsStructure's lat/lon/r/a/weight/z arrays that residual[k]
   corresponds to. Both arrays are malloc'd by the producer (fitBaseline()/
   computeRParams()/computeAzParams()) and owned by whichever caller eventually writes
   or discards them. */
typedef struct
{
	int32_t npts;
	int32_t *origIndex;
	double *residual;
} tieResidualsType;

OGRDataSourceH openTieResidualsGpkg(const char *gpkgPath, int32_t append);
int32_t writeTieResidualsLayer(OGRDataSourceH ds, const char *layerName,
								tiePointsStructure *tiePoints, tieResidualsType *fit,
								const char *residualFieldName, const char *residualDescription,
								double rotation, double stdLat, int32_t hemisphere);
void closeTieResidualsGpkg(OGRDataSourceH ds);
#endif
