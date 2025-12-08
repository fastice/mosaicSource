#define _XOPEN_SOURCE
#include <time.h>
#include "clib/standard.h"
#include "stdio.h"
#include <stdlib.h>
#include <string.h>
#include "mosaicSource/common/common.h"
#include "gdalIO/gdalIO/grimpgdal.h"

static void flagError(int count, int nReturned, char *varName)
{
    if( count != nReturned) 
        error("Could only parse %i values instead of %i for %s", nReturned, count, varName);
}

static void readVRTState(dictNode *metaData, stateV *sv)
{
    char velKey[19], posKey[19];
    flagError(1, sscanf(get_value(metaData, "NumberOfStateVectors"), "%d",
        &(sv->nState)), "NumberOfStateVectors");

    flagError(1, sscanf(get_value(metaData, "StateVectorInterval"), "%lf",
        &(sv->deltaT)), "StateVectorInterval");
    
    flagError(1, sscanf(get_value(metaData, "TimeOfFirstStateVector"), "%lf",
        &(sv->t0)), "TimeOfFirstStateVector");

    for(int i=1; i <= sv->nState; i++) 
    {
        sprintf(velKey, "SV_Vel_%d", i);
        sprintf(posKey, "SV_Pos_%d", i);
        flagError(3, sscanf(get_value(metaData, velKey), "[%lf, %lf, %lf]",
            &(sv->vx[i]), &(sv->vy[i]), &(sv->vz[i])), velKey);
        flagError(3, sscanf(get_value(metaData, posKey), "[%lf, %lf, %lf]",
            &(sv->x[i]), &(sv->y[i]), &(sv->z[i])), posKey);
    }
}

void parseSLCVrt(char *vrtFile, SARData *sarD, stateV *sv, int32_t *byteOrder)
{
    GDALDatasetH hDS;
    dictNode *metaData = NULL, *metaDataMT=NULL;
    // Read and Process metadata
    hDS = GDALOpen(vrtFile, GDAL_OF_READONLY);
	readDataSetMetaData(hDS, &metaData);
    //printDictionary(metaData);

    //sarD->label = strdup(vrtFile);
    flagError(6, sscanf(get_value(metaData, "datetime"), "%d-%d-%d %d:%d:%lf",
        &(sarD->year), &(sarD->month), &(sarD->day),
        &(sarD->hr), &(sarD->min), &(sarD->sec)), "datetime");
    
    // PRF
    flagError(1, sscanf(get_value(metaData, "PRF"), "%lf", &(sarD->prf)), "PRF");
    // Assume zero Doppler
    for (int i = 0; i < 4; i++) 
    {
        sarD->fd[i] = 0;
    }
    // time to first sample
    flagError(1, sscanf(get_value(metaData, "SLCFirstZeroDopplerTime"), "%lf",
        &(sarD->echoTD)), "SLCFirstZeroDopplerTime");
    // Near, far, and center range
    flagError(1, sscanf(get_value(metaData, "SLCNearRange"), "%lf",
        &(sarD->rn)), "SLCNearRange");
    flagError(1, sscanf(get_value(metaData, "SLCFarRange"), "%lf",
        &(sarD->rf)), " SLCFarRange");
    sarD->rc = 0.5 * (sarD->rn + sarD->rf);
    // Range and Azimuth pixel Size
    flagError(1, sscanf(get_value(metaData, "SLCRangePixelSize"), "%lf",
        &(sarD->slpR)), "SLCRangePixelSize");
    flagError(1, sscanf(get_value(metaData, "SLCAzimuthPixelSize"), "%lf",
        &(sarD->slpA)), "SLCAzimuthPixelSize");
    // Image Size
    flagError(1, sscanf(get_value(metaData, "SLCRangeSize"), "%d",
        &(sarD->nSlpR)), "SLCRangeSize");
    flagError(1, sscanf(get_value(metaData, "SLCAzimuthSize"), "%d",
        &(sarD->nSlpA)), "SLCAzimuthSize");
    // Get state
    readVRTState(metaData, sv);
    // Byte order
    *byteOrder = -1;
    if(strstr("LSB", get_value(metaData, "ByteOrder")) != NULL) 
    {
        *byteOrder = LSB;
    } else if(strstr("MSB", get_value(metaData, "ByteOrder")) != NULL) 
    {
        *byteOrder = MSB;   
    }  
}