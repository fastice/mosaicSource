#include "stdio.h"
#include "string.h"
#include "math.h"
/*
#include "mosaicSource/GeoCodeDEM_p/geocodedem.h"
#include "mosaicSource/computeVH_p/computevh.h"
#include "mosaicvh.h"*/
#include "common.h"
#include <stdlib.h>

/*
   Process input file for mosaicDEMs
*/
void parseIrregFile(char *irregFile, irregularData **irregData)
{
    FILE *fp;
    int32_t lineCount, eod;
    irregularData *tmp, *tmp1;
    int32_t i;
    char line[1024];
    /*
       Open file for input
    */
    fp = openInputFile(irregFile);
    /*
       Input nr,na,nlr,nla
    */
    eod = FALSE;
    while (eod == FALSE)
    {
        line[0] = '\0';
        lineCount = getDataString(fp, lineCount, line, &eod);

        i = strlen(line);
        if (i > 0 && eod == FALSE)
        {
            tmp = (irregularData *)malloc(sizeof(irregularData));
            tmp->file = (char *)malloc((i + 1) * sizeof(char));
            tmp->file = strcpy(tmp->file, line);
            tmp->file[i - 1] = '\0';

            tmp->next = NULL;
            tmp->x = NULL;
            tmp->y = NULL;
            tmp->vx = NULL;
            tmp->vy = NULL;
            tmp->link = NULL;
            if (*irregData == NULL)
            {
                *irregData = tmp;
                tmp1 = tmp;
            }
            else
                tmp1->next = tmp;
            tmp1 = tmp;
        }
    }

    return;
}
