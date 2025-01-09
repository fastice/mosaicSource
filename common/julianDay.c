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


void julian_to_gregorian(long jd, int *year, int *month, int *day) {
    long l, n, i, j, k;

    // Correct Julian Day Number for Gregorian Calendar
    l = jd + 68569;
    n = (4 * l) / 146097;
    l = l - (146097 * n + 3) / 4;
    i = (4000 * (l + 1)) / 1461001;
    l = l - (1461 * i) / 4 + 31;
    j = (80 * l) / 2447;
    *day = l - (2447 * j) / 80;
    l = j / 11;
    *month = j + 2 - (12 * l);
    *year = 100 * (n - 49) + i + l;
}

void jd_to_date_and_time(double jd, int *year, int *month, int *day, int *hour, int *minute, int *second) {
    long intPart = (long)jd;  // Integer part of the Julian Day Number
    double fracPart = jd - intPart;  // Fractional part of the Julian Day

    // Convert the integer part of the Julian Day Number to a Gregorian date
    julian_to_gregorian(intPart, year, month, day);

    // Convert the fractional part to time (hours, minutes, seconds)
    double dayFraction = fracPart * 24.0;  // Multiply by 24 to get the time in hours
    *hour = (int)dayFraction;
    dayFraction -= *hour;
    *minute = (int)(dayFraction * 60);
    *second = (int)((dayFraction * 60 - *minute) * 60);
}