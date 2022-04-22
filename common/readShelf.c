#include "stdio.h"
#include"string.h"
#include "math.h"
#include <stdlib.h>
#include "mosaicSource/common/common.h" 
/*#include "mosaic3d.h"*/

void readGeodatFile(char *geodatFile, double *x0, double *y0, double *deltaX, double *deltaY, int *sx, int *sy){
    char line[256];
    int lineCount=0, eod;
    FILE *fp;
    float dum1, dum2;
/* 
    Open geodat file
*/ 
    fp = openInputFile(geodatFile);
    if(fp == NULL) 
        error("*** readShelf: Error opening %s ***\n",geodatFile);
/*
   Read parameters
*/
    lineCount=getDataString(fp,lineCount,line,&eod); /* Skip # 2 line */
    lineCount=getDataString(fp,lineCount,line,&eod);
    sscanf(line,"%f %f\n",&dum1,&dum2); /* read as float in case fp value */
    *sx = (int)dum1;
    *sy = (int)dum2; 

    lineCount=getDataString(fp,lineCount,line,&eod);
    sscanf(line,"%lf %lf\n",deltaX, deltaY);
    *deltaX *= MTOKM;  *deltaY *= MTOKM;

    lineCount=getDataString(fp,lineCount,line,&eod);
    sscanf(line,"%lf %lf\n",x0, y0);
    fprintf(stderr,"%s\n",line);
    fclose(fp);
}
/*
   Process input file for mosaicDEMs
*/
    void readShelf(outputImageStructure *outputImage,char *shelfMaskFile)
{
    FILE *fp;
    char *geodatFile;
    ShelfMask *shelfMask;
    
    unsigned char *tmp,*tmp1;
    float dum1,dum2;
    int i;
    fprintf(outputImage->fpLog,";\n; Entering readShelf(.c)\n;\n");
    fflush(outputImage->fpLog);
    shelfMask=(ShelfMask *)malloc(sizeof(ShelfMask));
    outputImage->shelfMask=shelfMask;
/*
     geodat file name
*/
   fprintf(stderr,"ShelfMaskFile %s\n",shelfMaskFile);
   geodatFile = (char *)malloc(strlen(shelfMaskFile)+8);
   geodatFile[0]='\0';
   geodatFile=strcpy(geodatFile,shelfMaskFile);
   geodatFile=strcat(geodatFile,".geodat"); 
   fprintf(stderr,"Shelfmask geodat file %s\n",geodatFile); 
   readGeodatFile(geodatFile, &(shelfMask->x0), &(shelfMask->y0), &(shelfMask->deltaX), 
                    &(shelfMask->deltaY), &(shelfMask->xSize), &(shelfMask->ySize) );
   shelfMask->rot=0;
   shelfMask->hemisphere=SOUTH;
   shelfMask->stdLat = 71.0;

    fprintf(stderr,"%i %i \n %f %f \n %f %f \n %f %i %f \n",
       shelfMask->xSize,shelfMask->ySize,
       shelfMask->deltaX,shelfMask->deltaY,
       shelfMask->x0,shelfMask->y0,
       shelfMask->rot,shelfMask->hemisphere,shelfMask->stdLat);
    /*
       Malloc array
    */
     shelfMask->mask = (unsigned char **)
	  malloc(shelfMask->ySize * sizeof(unsigned char *));
     tmp = (unsigned char *)malloc(shelfMask->xSize *shelfMask->ySize * 
           sizeof(unsigned char));
/* 
    Open shelfFile file
*/
    fp = openInputFile(shelfMaskFile);
    if(fp == NULL) 
        error("*** readShelf: Error opening %s ***\n",shelfMaskFile);
    for(i=0; i < shelfMask->ySize; i++) { 
        tmp1 =&(tmp[i*shelfMask->xSize]);
        freadBS(tmp1,sizeof(unsigned char),shelfMask->xSize,fp,BYTEFLAG);
        shelfMask->mask[i]=tmp1;
    }
    fprintf(outputImage->fpLog,";\n; Leaving readShelf(.c)\n;\n");
    fflush(outputImage->fpLog);
    fclose(fp);
} 


