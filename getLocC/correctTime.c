#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "math.h"
#include "mosaicSource/common/common.h"

/*
  Correct starting time for squint and skew
*/
void correctTime(SARData *sarD, double *squint, double noffset, double *tskew, double *toffset, int32_t squintTime)
{
  /* ADDED FLAG 12/20/99 for fixing RS data */
  if (squintTime == FALSE)
  {
    *tskew = (sarD->rc * sin((DTOR * (*squint))));
    *tskew = *tskew / (sarD->slpA * sarD->prf);
  }
  else
  {
    *tskew = *squint;
    *squint = *tskew * (sarD->slpA * sarD->prf) / sarD->rc;
    *squint = asin(*squint) * RTOD;
  }
  *toffset = (double)noffset / sarD->prf;
  /*
     Correct time
  */
  sarD->sec += *tskew + *toffset;

  while (sarD->sec > 60.0)
  {
    sarD->sec -= 60.0;
    sarD->min += 1;
  }
  while (sarD->min > 59)
  {
    sarD->min -= 60;
    sarD->hr += 1;
  }
  if (sarD->hr > 24)
    fprintf(stderr,
            "\n\n**** WARNING TIME CROSSING 24 hour boundary ***\n\n");
}
