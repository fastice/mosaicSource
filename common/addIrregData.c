#include "stdio.h"
#include"string.h"
#include <math.h>
/*#include "cRecipes/nrutil.h"*/
#include "common.h"

static void getTriBox(int index,int *iMin,int *iMax,int *jMin,int *jMax, irregularData *data,outputImageStructure *outputImage);
static int inTriangle(double x, double y, int n, irregularData *data);
static int goodTriangle(int n, irregularData *data);
static void linTriInterp(double x, double y, int n, irregularData *data, double *vx, double *vy);

/*
  Make mosaic of insar and other dems.
*/
void addIrregData(irregularData *irregDat, outputImageStructure *outputImage, float fl)
{   
	irregularData *currentData;
	float **vXimage,**vYimage,**vZimage;
	float **scaleX,**scaleY,**scaleZ;
	double scX,scY;
	int *triangles, nTri;
	double x,y;
	double *xp,*yp,*vxp,*vyp;
	float **vxTmp,**vyTmp,**vzTmp,**fScale,**sxTmp,**syTmp;
	int iMin, iMax, jMin,jMax;
	double ex,ey;
	double vx,vy;
	int i,j,n,k;
	/*
	  Pointers to output images
	*/
	vXimage = (float **)outputImage->image;
	vYimage = (float **)outputImage->image2;
	vZimage = (float **)outputImage->image3;
	scaleX = (float **)outputImage->scale;
	scaleY = (float **)outputImage->scale2;
	scaleZ = (float **)outputImage->scale3;
	vxTmp = (float **) outputImage->vxTmp;
	vyTmp = (float **) outputImage->vyTmp;
	vzTmp = (float **) outputImage->vzTmp;
	sxTmp = (float **) outputImage->sxTmp;
	syTmp = (float **) outputImage->syTmp;
	fScale= (float **) outputImage->fScale;
	/*
	  Compute feather scale for existing 
	*/
	if(fl == 0) {
		for (j=0; j < outputImage->ySize; j++) 
			for(k=0; k < outputImage->xSize; k++) fScale[j][k]=1.0;
	} else {
		computeScale((float **)vXimage,fScale, outputImage->ySize,  outputImage->xSize,fl,(float)1.0,(double)(-LARGEINT));
	}
	/*
	  Init array.
	*/
	for (j=0; j < outputImage->ySize; j++) {
		for(k=0; k < outputImage->xSize; k++) {
			vXimage[j][k] *= scaleX[j][k]*fScale[j][k];
			vYimage[j][k] *= scaleY[j][k]*fScale[j][k];
			scaleX[j][k] *= fScale[j][k];
			scaleY[j][k] *= fScale[j][k];
		}
	}
	/*
	  Loop over irregular data sets
	*/
	for(currentData=irregDat;currentData != NULL;currentData=currentData->next){
		triangles=currentData->link;
		nTri=currentData->nTri;
		xp=currentData->x;  yp=currentData->y;
		vxp=currentData->vx; vyp=currentData->vy;
		/* Init arrays */
		for (j=0; j < outputImage->ySize; j++) {
			for(k=0; k < outputImage->xSize; k++) {
				sxTmp[j][k]=0.0; syTmp[j][k]=0.0;
				vxTmp[j][k]=-2.0E9; vyTmp[j][k]=-2.0e9;      
			}
		}

		for(n=0; n < nTri; n++) {
			/*
			  Compute bounding box for triangle
			*/
			getTriBox(n,&iMin,&iMax, &jMin, &jMax, currentData, outputImage);

			/*
			  Loop over points in box
			*/
			for(i=iMin; i <= iMax; i++) {
				y = (outputImage->originY + i*outputImage->deltaY) * MTOKM;
				for(j=jMin; j <= jMax; j++) {  
					x = (outputImage->originX + j*outputImage->deltaX) * MTOKM;
					/*
					  If a point lies in a triangle, interpolate
					*/
					if(inTriangle(x,y,n,currentData) == TRUE) {
						if(goodTriangle(n,currentData) == TRUE) {
							/*
							  Interpolate
							*/
							linTriInterp(x,y,n,currentData,&vx,&vy);
							/*     ************ Improve scales      */
							ex=200.0;
							ey=200.0;
							scX=1.0/(ex*ex);
							scY=1.0/(ey*ey);
							vxTmp[i][j] = (float)vx*scX; 
							vyTmp[i][j] = (float)vy*scY; 
							sxTmp[i][j] = scX;
							syTmp[i][j] = scY;
						}
					} /* End if inTr.. */
				} /* End for(j.. */
			}  /* End for i... */
		}  /* End for n... */
		/*
		  Compute scale array for feathering.
		*/    
		if(fl > 0 )
			computeScale((float **)vxTmp,fScale, outputImage->ySize,  outputImage->xSize,fl,(float)1.0,(double)(-LARGEINT));
		/*
		  Now sum current result. 
		*/         
		for(i=0; i <  outputImage->ySize; i++) {
			for(j=0; j < outputImage->xSize; j++) {
				if(vxTmp[i][j] > (-LARGEINT+1)) {
					vXimage[i][j] += vxTmp[i][j]*fScale[i][j]; 
					scaleX[i][j]  += sxTmp[i][j]*fScale[i][j];
					vYimage[i][j] += vyTmp[i][j]*fScale[i][j]; 
					scaleY[i][j]  += syTmp[i][j]*fScale[i][j];     
				}
			}
		}

	} /* End for currentData .. */
	/**************************END OF MAIN LOOP ******************************/
	/*
	  Adjust scale
	*/
	for (j=0; j < outputImage->ySize; j++) {
		for(k=0; k < outputImage->xSize; k++) {
			if(scaleX[j][k] > 0.0)  vXimage[j][k] /= scaleX[j][k];
			else vXimage[j][k] = -LARGEINT;
			if(scaleY[j][k] > 0.0)  vYimage[j][k] /= scaleY[j][k];
			else vYimage[j][k] = -LARGEINT;
		} 
	}
}

/*
  Linear interp using triangle n.  Adapted from GMT 3.1
*/ 
static void linTriInterp(double x, double y, int n, irregularData *data,  double *vx, double *vy)
{
	double ux[3],uy[3];
	double *xx,*yy, *zvx, *zvy;
	double zkjx, zkjy,zljx,zljy,zjx,zjy,zkx,zky,zlx,zly;
	double ax,ay,bx,by,cx,cy,f;
	double xkj,ykj,xlj,ylj;
	int *link;

	link=data->link;
	xx=data->x; yy=data->y; zvx=data->vx; zvy = data->vy;
	/* 
	   Find equation for the plane as z = ax + by + c
	*/
	n=n*3;		
	ux[0] =  xx[link[n]];  uy[0] =  yy[link[n]];
	zjx = zvx[link[n]];  zjy=zvy[link[n]]; n++;
	ux[1] = xx[link[n]];   uy[1] = yy[link[n]];
	zkx = zvx[link[n]];  zky=zvy[link[n]]; n++;
	ux[2] = xx[link[n]];   uy[2] = yy[link[n]];
	zlx = zvx[link[n]];  zly=zvy[link[n]]; n++;

	xkj = ux[1] - ux[0];  ykj = uy[1] - uy[0];
	xlj = ux[2] - ux[0];  ylj = uy[2] - uy[0];

	zkjx = zkx - zjx; zkjy = zky - zjy;
	zljx = zlx - zjx; zljy = zly - zjy;
			
	f = 1.0 / (xkj * ylj - ykj * xlj);
	ax = -f * (ykj * zljx - zkjx * ylj);   ay = -f * (ykj * zljy - zkjy * ylj);
	bx = -f * (zkjx * xlj - xkj * zljx);   by = -f * (zkjy * xlj - xkj * zljy);
	cx = -ax * ux[1] - bx * uy[1] + zkx;   cy = -ay * ux[1] - by * uy[1] + zky;

  
	*vx = ax*x + bx*y + cx;
	*vy = ay*x + by*y + cy;
}

/*
  Perform various tests to throw out unwanted triangles
*/ 
static int goodTriangle(int n, irregularData *data)
{
	double x1,x2,x3,y1,y2,y3;
	double l1,l2,l3;
	double area;

	x1=data->x[data->link[n*3]]; y1=data->y[data->link[n*3]];
	x2=data->x[data->link[n*3+1]]; y2=data->y[data->link[n*3+1]];
	x3=data->x[data->link[n*3+2]]; y3=data->y[data->link[n*3+2]];
	/*
	  Discard triangle if edge to long
	*/
	l1 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	l2 = sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3));
	l3 = sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));
	if(l1 >= data->maxLength || l2 >= data->maxLength || l3 >= data->maxLength) {
		return(FALSE);
	}
	/*
	  Discard triangle area too large
	*/
	area=(x2*y3-y2*x3) + (x1*y2 - y1*x2) - (x1*y3 - y1*x3);
	area *= 0.5;
	if(area > data->maxArea) return(FALSE);

	return(TRUE);

}

/*
  Check if in triangle
*/
static int inTriangle(double x, double y, int n, irregularData *data)
{
	double vx1,vx2,vx3,vy1,vy2,vy3;
	double xp12,xp23,xp31;

	vx1=data->x[data->link[n*3]] - x; vy1=data->y[data->link[n*3]] - y;
	vx2=data->x[data->link[n*3+1]]-x; vy2=data->y[data->link[n*3+1]] - y;
	vx3=data->x[data->link[n*3+2]]-x; vy3=data->y[data->link[n*3+2]] - y;

	xp12= vx1*vy2 - vy1 * vx2;
	xp23 =vx2*vy3 - vy2 * vx3;
	xp31 =vx3*vy1 - vy3 * vx1;

	if( (xp12 >= 0 && xp23 >= 0.0 && xp31 >= 0.0) ||
	    (xp12 <= 0.&& xp23 <= 0.0 && xp31 <= 0.0) ) return(TRUE);

	return(FALSE);
}


/*
  Find region where image exists
*/
static void getTriBox(int index,int *iMin,int *iMax,int *jMin,int *jMax,  irregularData *data,outputImageStructure *outputImage)
{
	double xa1,ya1;
	double minX,maxX,minY,maxY;
	int i;
	/*
	  Loop through points to find max and min locations
	*/   
	for(i=0; i < 3; i++) {
		xa1 = data->x[data->link[index*3+i]];
		ya1 = data->y[data->link[index*3+i]];
		if(i==0) {minX=xa1; maxX=xa1; minY=ya1; maxY=ya1;}
		else {
			minX=min(minX,xa1); maxX=max(maxX,xa1);
			minY=min(minY,ya1); maxY=max(maxY,ya1);
		}
	}   
	/*
	  Compute i,j min,max with pad
	*/ 
	*iMin=(int)((minY*KMTOM-outputImage->originY)/outputImage->deltaY);
	*jMin=(int)((minX*KMTOM-outputImage->originX)/outputImage->deltaX);
	*iMax=(int)((maxY*KMTOM-outputImage->originY)/outputImage->deltaY);
	*jMax=(int)((maxX*KMTOM-outputImage->originX)/outputImage->deltaX);
	(*iMin)--; (*jMin)--; (*iMax)++; (*jMax)++;

	*iMin=max(*iMin,0); *jMin=max(*jMin,0);
	*iMax=min(outputImage->ySize-1,*iMax); *jMax=min(outputImage->xSize-1,*jMax);
}


