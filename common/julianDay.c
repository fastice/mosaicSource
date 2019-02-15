#include "math.h"
#include "string.h"
#include "float.h"
#include "common.h"


long julday(int mm, int id, int iyyy)
{
        long jul;
        int ja,jy=iyyy,jm;

        if (jy < 0) ++jy;

        if (mm > 2) {
                jm=mm+1;
        } else {
                --jy;
                jm=mm+13;
        }
	/*(15+31*(10+12*1582))*/
        jul = (long) (floor(365.25*jy)+floor(30.6001*jm)+id+1720995);

        if ( id+31L*(mm+12*iyyy) >=  (15+31*(10+12*1582))) {
                ja=(int)(0.01 * jy);
                jul += 2-ja+ (int)(0.25*ja);
        }
        return jul;
}

double juldayDouble(int mm, int id, int iyyy) {
	return ((double)julday(mm,id,iyyy) -0.5);
}
