#include "math.h"
#include "common.h"
#include "cRecipes/nrutil.h"

static void intError(double *xa,double x) {
	int j;
	fprintf(stderr,"ERROR POLINTVEC\n");
	for (j=1;j<=NUSESTATE;j++) 
		fprintf(stderr,"%i %f %f\n",j, xa[j], x);                           
	        error("Error in routine polintvec");
}
/* This is based on POLINT, but it works about 4x faster than 6 seperate calls */
void polintVec(double xa[], double y1[],double y2[],double y3[],double y4[],double y5[],double y6[],
		      double x, double *yr1,double *yr2,double *yr3,double *yr4,double *yr5,double *yr6 )
{
	int i,m,ns=1,j, k;
	int na;
	double denInv,den,dif,dift,ho,hp,w;
	double c[7][NUSESTATE+1],d[7][NUSESTATE+1];

	dif=fabs(x-xa[1]);
	for (i=1;i<=NUSESTATE;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[1][i]=y1[i]; c[2][i]=y2[i]; c[3][i]=y3[i]; c[4][i]=y4[i]; c[5][i]=y5[i]; c[6][i]=y6[i];
		d[1][i]=y1[i]; d[2][i]=y2[i]; d[3][i]=y3[i]; d[4][i]=y4[i]; d[5][i]=y5[i]; d[6][i]=y6[i];
	}
	*yr1=y1[ns]; *yr2=y2[ns];  *yr3=y3[ns]; *yr4=y4[ns]; *yr5=y5[ns]; 	*yr6=y6[ns];
	ns--;	
	for (m=1;m<NUSESTATE;m++) {
		for (i=1;i<=NUSESTATE-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			if ( (den=ho-hp) == 0.0) intError(xa,x);
			denInv=1./den;
			for( k=1; k<= 6; k++) {						
				w=c[k][i+1]-d[k][i];
				den=w*denInv;
				d[k][i]=hp*den;
				c[k][i]=ho*den;
			}
		}
		na=ns+1;
		if(2*ns < (NUSESTATE-m)) {
			*yr1+=c[1][na];
			*yr2+=c[2][na];
			*yr3+=c[3][na];
			*yr4+=c[4][na];
			*yr5+=c[5][na];
			*yr6+=c[6][na];								
		} else {
			*yr1+=d[1][ns];
			*yr2+=d[2][ns];
			*yr3+=d[3][ns];
			*yr4+=d[4][ns];
			*yr5+=d[5][ns];
			*yr6+=d[6][ns];
			ns--;						
		}
	}
}
#undef NRANSI
