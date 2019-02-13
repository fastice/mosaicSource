#include "stdio.h"
#include "string.h"
#include "common.h"
/* 
   Read through file returning data strings.
   If end of data eod = TRUE else eod = FALSE.

   This allows a commented line, with a special character to be returned
   e.g., ;*
*/

int getDataStringSpecial(FILE *fp, int lineCount, char *line,int *eod, char special, int *specialFound)
{
	int lineLength;  
	int count;
	count = 0;
	*specialFound=FALSE;
	while( TRUE ) {  /* Loop to read lines */
		lineLength = fgetline(fp,line,LINEMAX); /* Read line */
		lineCount++;
		count++;
		if( strchr(line,COMMENT) == NULL &&   
		    strpbrk(line, "0123456789.ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") != NULL ) {
			*eod = FALSE;
			return lineCount; /* Data */ 
		} else if(  strchr(line,special)  != NULL  && strpbrk(line, "0123456789.ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")  != NULL ) {
			*eod = FALSE;
			*specialFound=TRUE;
			return lineCount; /* Data */ 
		}
		if( strchr(line,ENDDATA) != NULL &&  strchr(line,COMMENT) == NULL ) { /* End of data */
			*eod = TRUE;
			return lineCount;
		}
		if( count >= 100000) error("getDataString: Missing end of data");
	}
}
