#include "common.h"

void initFloatMatrix(float **x,long int nr,long int nc,float initValue) {
	long int j,k;
	for(j=0; j < nr; j++)
		for(k=0; k < nc; k++) x[j][k]=initValue;
}

void initDoubleMatrix(double **x,long int nr,long int nc,double initValue) {
	long int j,k;
	for(j=0; j < nr; j++)
		for(k=0; k < nc; k++) x[j][k]=initValue;
}
