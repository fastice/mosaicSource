#include "stdio.h"
#include"string.h"
#include "math.h"
/*#include "source/GeoCodeDEM_p/geocodedem.h"*/
/*
#include "source/computeVH_p/computevh.h"
#include "mosaicvh.h"
*/
#include <stdlib.h>
#include "common.h"
#define VNULL ((void *)NULL)
#define EXIT_FAILURE_1 -1
static int GMT_delaunay (double *x_in, double *y_in, int n, int **link);
static void *GMT_memory (void *prev_addr, size_t nelem, size_t size, char *progname);
/*
   Process input file for mosaicDEMs
*/
    void  getIrregData(irregularData *irregData)
{
    extern int HemiSphere;
    extern double Rotation;
    FILE *fp;
    double slat;
    int lineCount, eod;
    irregularData *currentData;
    int i;
    double lat,lon,vx,vy,x,y,speed,bearing;
    double xa,ya,x2,y2,x2a,y2a,theta,phi;
    char line[2048];
    int nLines;
    slat=70.0;
    if(HemiSphere==SOUTH) slat=71.0;
    fprintf(stderr,"slat = %f\n",slat);
/*
   Open file for input
*/
    for(currentData=irregData;currentData !=NULL;currentData=currentData->next){
/*
   Pass one: determine number of data lines
*/  
        fp = openInputFile(currentData->file);
        eod=FALSE;
        nLines=0;
        lineCount=getDataString(fp,lineCount,line,&eod); /* Skip # line */
        while(eod == FALSE) {
            lineCount=getDataString(fp,lineCount,line,&eod);
            if(strlen(line) > 5) nLines++;
        }  
/*
    Malloc space
*/    
        currentData->x = (double *)malloc(sizeof(double)*nLines);  
        currentData->y = (double *)malloc(sizeof(double)*nLines);  
        currentData->vx = (double *)malloc(sizeof(double)*nLines);  
        currentData->vy = (double *)malloc(sizeof(double)*nLines);
        currentData->nData=nLines;  
/*
   Pass 2: read and convert data
*/
        fprintf(stderr,"%s has %i lines\n",currentData->file,nLines);
        rewind(fp);
        lineCount=getDataString(fp,lineCount,line,&eod); /* Skip # line */
        i=0;
        eod=FALSE;
        while(eod == FALSE) {
            lineCount=getDataString(fp,lineCount,line,&eod);
            if(strlen(line) > 5) {
                if(sscanf(line,"%lf%lf%lf%lf",&lat,&lon,&speed,&bearing) != 4) 
                 error("missing data in %s at line %i\n %s\n",
                    currentData->file,lineCount,line);

                lltoxy1(lat,lon,&x,&y,Rotation,slat);
                currentData->x[i] = x; currentData->y[i]=y;
/*
   Convert to xy velocities
*/
                theta = atan2(y,x);
                x2=(x*cos(theta) + y*sin(theta))*KMTOM;
                y2=(-x*sin(theta) + y*cos(theta))*KMTOM;
                phi=PI;
                if(HemiSphere == SOUTH) phi=0.0;
                bearing=bearing*DTOR;
                x2a=x2 + speed*cos(-(bearing-phi));
                y2a=y2 + speed*sin(-(bearing-phi));
                xa = x2a*cos(theta) - y2a*sin(theta);
                ya = x2a*sin(theta) + y2a*cos(theta);
                vx=xa-x*KMTOM; 
                vy=ya-y*KMTOM;
                currentData->vx[i]=vx; currentData->vy[i]=vy;
/*
   Update counter
*/
                i++;
            }
        } 
/*
  Do Delauny triangulation using GMT routine.
*/ 

        fclose(fp);
        currentData->nTri=GMT_delaunay(currentData->x,currentData->y, nLines,&(currentData->link));
    }
     return;
}

/*
    Routines from GMT.
*/

#define REAL double

#include "triangle/triangle.h"

static int GMT_delaunay (double *x_in, double *y_in, int n, int **link)
{
	/* GMT interface to the triangle package; see above for references.
	 * All that is done is reformatting of parameters and calling the
	 * main triangulate routine.  Thanx to Alain Coat for the tip.
	 */

	int i, j;
	struct triangulateio In, Out, vorOut;

	/* Set everything to 0 and NULL */

	memset ((void *)&In,	 0, sizeof (struct triangulateio));
	memset ((void *)&Out,	 0, sizeof (struct triangulateio));
	memset ((void *)&vorOut, 0, sizeof (struct triangulateio));

	/* Allocate memory for input points */

	In.numberofpoints = n;
	In.pointlist = (double *) GMT_memory ((void *)NULL, (size_t)(2 * n), sizeof (double), "GMT_delaunay");

	/* Copy x,y points to In structure array */

	for (i = j = 0; i < n; i++) {
		In.pointlist[j++] = x_in[i];
		In.pointlist[j++] = y_in[i];
	}

	/* Call Jonathan Shewchuk's triangulate algorithm */

	triangulate ("zIQB", &In, &Out, &vorOut);

	*link = Out.trianglelist;	/* List of node numbers to return via link */

	if (Out.pointlist) free ((void *)Out.pointlist);

	return (Out.numberoftriangles);
}


static void *GMT_memory (void *prev_addr, size_t nelem, size_t size, char *progname)
{
	void *tmp;

	if (nelem == 0) return(VNULL); /* Take care of n = 0 */
	
	if (prev_addr) {
		if ((tmp = realloc ((void *) prev_addr, (size_t)(nelem * size))) == VNULL) {
			fprintf (stderr, "GMT Fatal Error: %s could not reallocate more memory, n = %d\n", progname, nelem);
			exit (EXIT_FAILURE_1);
		}
	}
	else {
		if ((tmp = calloc ((size_t) nelem, (unsigned) size)) == VNULL) {
			fprintf (stderr, "GMT Fatal Error: %s could not allocate memory, n = %d\n", progname, nelem);
			exit (EXIT_FAILURE_1);
		}
	}
	return (tmp);
}
