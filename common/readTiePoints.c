#include "stdio.h"
#include "string.h"
#include "mosaicSource/common/common.h"
/*#include "tiePoints.h"*/
#include <stdlib.h>

/*
  Read tiepoint file for tiepoints.
*/
void readTiePoints(FILE *fp, tiePointsStructure *tiePoints, int noDEM)
{
	double lat, lon;
	int linelength;   /* Input line length */
	int notdone;      /* Loop flag */
	double z=0;
	double vx=0.0, vy=0.0,vz=0.0;
	char lineBuffer[LINEMAX+1];
	char *line;       /* Input line buffer */
	int lineCount=0;
	if(noDEM == TRUE) fprintf(stderr,"NO DEM used\n");
	line = lineBuffer;  /* Allocate line buffer */
	notdone = TRUE;
  
	/* Modified 4/02/07 to fix crash with declared in structure */
	tiePoints->bsq = (double *)malloc(sizeof(double)*MAXTIEPOINTS); /* Modified 12/09/94 to vector */
	tiePoints->lat = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->lon = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->x = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->y = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->z = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->r = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->a = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->phase = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->delta = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->vx = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->vy = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->vz = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	tiePoints->vyra = (double *)malloc(sizeof(double)*MAXTIEPOINTS);
	/* End 4/2/7 fix */

	tiePoints->npts=0;
	while( notdone == TRUE ) {                 /* Loop to read lines */
		linelength = fgetline(fp,line,LINEMAX); /* Read line */
		lineCount++;
		if( strchr(line,ENDDATA) != NULL ) 
			notdone = FALSE;                    /* End of data, set exit flag */
		else if( strchr(line,COMMENT) == NULL ) {   /* If not comment, parse */
			if(tiePoints->motionFlag == TRUE) { 
				if(sscanf(line,"%lf%lf%lf%lf%lf%lf",  &lat,&lon,&z,&vx,&vy,&vz) != 6) {/*motion ties*/
					fprintf(stderr,"%lf %lf %lf %lf %lf %lf",lat,lon,z,vx,vy,vz);
					error("%s  %i", "readTiePoints: -- Invalid # of parameters at line:", lineCount);
				}
			} else if(tiePoints->vrFlag == TRUE) { 
				if(sscanf(line,"%lf%lf%lf%lf%lf%lf",  &lat,&lon,&z,&vx,&vy,&vz) != 5) {/*motion ties*/
					fprintf(stderr,"%lf %lf %lf %lf %lf %lf",lat,lon,z,vx,vy,vz);
					error("%s  %i", "readTiePoints: -- Invalid # of parameters at line:",   lineCount);
				}
			} else { /* Regular lat/lon/z points */
				if(sscanf(line,"%lf%lf%lf",&lat,&lon,&z) != 3) {/*z too*/
					fprintf(stderr,"%f %f %f",lat,lon,z);
					error("%s  %i","readTiePoints: -- Invalid # of parameters at line:",  lineCount);
				}
			}
			/* 
			   Assign tiepoints and update counter
			*/
			tiePoints->lat[tiePoints->npts] = lat;
			tiePoints->lon[tiePoints->npts] = lon;
			tiePoints->z[tiePoints->npts] = z;    
			tiePoints->vx[tiePoints->npts] = vx;     
			tiePoints->vy[tiePoints->npts] = vy;     
			tiePoints->vz[tiePoints->npts] = vz;       
			(tiePoints->npts)++;
			if(tiePoints->npts >= MAXTIEPOINTS) /* Too many tiepts ? */
				error("readTiePoints -- MAXTIEPOINTS=%i exceeded", MAXTIEPOINTS);
		} /* End else */
	} /* End while */
}
