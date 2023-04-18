#include "stdio.h"
#include "string.h"
#include "common.h"
/*
   Read through file returning data strings.
   If end of data eod = TRUE else eod = FALSE.

   This allows a commented line, with a special character to be returned
   e.g., ;*
*/

int32_t getDataStringSpecial(FILE *fp, int32_t lineCount, char *line, int32_t *eod, char special, int32_t *specialFound)
{
	int32_t lineLength;
	int32_t count;
	count = 0;
	*specialFound = FALSE;
	while (TRUE)
	{											  /* Loop to read lines */
		lineLength = fgetline(fp, line, LINEMAX); /* Read line */
		lineCount++;
		count++;
		if (strchr(line, COMMENT) == NULL &&
			strpbrk(line, "0123456789.ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") != NULL)
		{
			*eod = FALSE;
			return lineCount; /* Data */
		}
		else if (strchr(line, special) != NULL && strpbrk(line, "0123456789.ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") != NULL)
		{
			*eod = FALSE;
			*specialFound = TRUE;
			return lineCount; /* Data */
		}
		if (strchr(line, ENDDATA) != NULL && strchr(line, COMMENT) == NULL)
		{ /* End of data */
			*eod = TRUE;
			return lineCount;
		}
		if (count >= 100000)
			error("getDataString: Missing end of data");
	}
}
