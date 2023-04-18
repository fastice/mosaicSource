#include "common.h"

void initFloatMatrix(float **x, int32_t nr, int32_t nc, float initValue)
{
	int32_t j, k;
	for (j = 0; j < nr; j++)
		for (k = 0; k < nc; k++)
			x[j][k] = initValue;
}

void initDoubleMatrix(double **x, int32_t nr, int32_t nc, double initValue)
{
	int32_t j, k;
	for (j = 0; j < nr; j++)
		for (k = 0; k < nc; k++)
			x[j][k] = initValue;
}
