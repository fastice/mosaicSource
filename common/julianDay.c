#include "math.h"
#include "string.h"
#include "float.h"
#include "common.h"

int32_t julday(int32_t mm, int32_t id, int32_t iyyy)
{
        int32_t jul;
        int32_t ja, jy = iyyy, jm;

        if (jy < 0)
                ++jy;

        if (mm > 2)
        {
                jm = mm + 1;
        }
        else
        {
                --jy;
                jm = mm + 13;
        }
        /*(15+31*(10+12*1582))*/
        jul = (int32_t)(floor(365.25 * jy) + floor(30.6001 * jm) + id + 1720995);

        if (id + 31L * (mm + 12 * iyyy) >= (15 + 31 * (10 + 12 * 1582)))
        {
                ja = (int)(0.01 * jy);
                jul += 2 - ja + (int)(0.25 * ja);
        }
        return jul;
}

double juldayDouble(int32_t mm, int32_t id, int32_t iyyy)
{
        return ((double)julday(mm, id, iyyy) - 0.5);
}
