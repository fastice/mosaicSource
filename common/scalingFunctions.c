#include "math.h"
#include "common.h"

/*
   At the start of each mosaicking round, undo the prior normalization prior to starting.
 */
void undoNormalization(outputImageStructure *outputImage, float **vXimage, float **vYimage, float **vZimage, float **errorX, float **errorY, float **scaleX, float **scaleY, float **scaleZ, float **fScale, int32_t statsFlag)
{
	int32_t j, k;
	if (statsFlag == FALSE)
	{
		for (j = 0; j < outputImage->ySize; j++)
		{
			for (k = 0; k < outputImage->xSize; k++)
			{
				vXimage[j][k] *= scaleX[j][k] * fScale[j][k];
				vYimage[j][k] *= scaleY[j][k] * fScale[j][k];
				vZimage[j][k] *= scaleZ[j][k] * fScale[j][k];
				errorX[j][k] *= (scaleX[j][k] * scaleX[j][k]) * fScale[j][k] * fScale[j][k];
				errorY[j][k] *= (scaleY[j][k] * scaleY[j][k]) * fScale[j][k] * fScale[j][k];
				scaleX[j][k] *= fScale[j][k];
				scaleY[j][k] *= fScale[j][k];
				scaleZ[j][k] *= fScale[j][k];
			}
		}
	}
}
/*
  For intermediate products, do scaling.
*/
void redoNormalization(float myWeight, outputImageStructure *outputImage, int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax,
					   float **vXimage, float **vYimage, float **vZimage, float **errorX, float **errorY, float **scaleX, float **scaleY, float **scaleZ, float **fScale,
					   float **vxTmp, float **vyTmp, float **vzTmp, float **sxTmp, float **syTmp, int32_t statsFlag)
{
	float weight;
	int32_t i, j;

	for (i = iMin; i < iMax; i++)
	{
		for (j = jMin; j < jMax; j++)
		{
			if (vxTmp[i][j] > (-LARGEINT + 1))
			{
				/* added weight Dec 2016 - need to verify it workds */
				if (statsFlag == FALSE)
					weight = myWeight;
				else
					weight = 1.0;
				vXimage[i][j] += vxTmp[i][j] * fScale[i][j] * weight;
				scaleX[i][j] += sxTmp[i][j] * fScale[i][j] * weight;
				vYimage[i][j] += vyTmp[i][j] * fScale[i][j] * weight;
				scaleY[i][j] += syTmp[i][j] * fScale[i][j] * weight;
				if (outputImage->timeOverlapFlag == TRUE)
				{
					vZimage[i][j] += vzTmp[i][j] * fScale[i][j] * weight;
					scaleZ[i][j] += sqrt(syTmp[i][j] * sxTmp[i][j]) * fScale[i][j] * weight;
				}
				else
				{
					vZimage[i][j] += vzTmp[i][j]; /*fScale[i][j] * weight; */
					scaleZ[i][j] += 1;
				}
				if (statsFlag == FALSE)
				{
					errorX[i][j] += fScale[i][j] * fScale[i][j] * sxTmp[i][j] * weight * weight;
					errorY[i][j] += fScale[i][j] * fScale[i][j] * syTmp[i][j] * weight * weight;
				}
				else
				{
					/* Sum for variance */
					errorX[i][j] += vxTmp[i][j] * vxTmp[i][j];
					errorY[i][j] += vyTmp[i][j] * vyTmp[i][j];
				}
			}
		}
	}
}

/*
	last scaling after all products have been added for a round.
*/

void endScale(outputImageStructure *outputImage, float **vXimage, float **vYimage, float **vZimage, float **errorX, float **errorY, float **scaleX, float **scaleY, float **scaleZ, int32_t statsFlag)
{
	int32_t j, k;

	for (j = 0; j < outputImage->ySize; j++)
	{
		for (k = 0; k < outputImage->xSize; k++)
		{
			if (scaleX[j][k] > 0.0)
				vXimage[j][k] /= scaleX[j][k];
			else
				vXimage[j][k] = -LARGEINT;
			if (scaleY[j][k] > 0.0)
				vYimage[j][k] /= scaleY[j][k];
			else
				vYimage[j][k] = -LARGEINT;
			if (scaleZ[j][k] > 0.0)
				vZimage[j][k] /= scaleZ[j][k];
			else
				vZimage[j][k] = -LARGEINT;
			if (statsFlag == FALSE)
			{
				if (scaleX[j][k] > 0.0)
					errorX[j][k] /= (scaleX[j][k] * scaleX[j][k]);
				else
					errorX[j][k] = -LARGEINT;
				if (scaleY[j][k] > 0.0)
					errorY[j][k] /= (scaleY[j][k] * scaleY[j][k]);
				else
					errorY[j][k] = -LARGEINT;
			}
			else
			{
				/* Compute sigma^2 - sqrt applied on output*/
				if (scaleX[j][k] > 1.0)
				{
					errorX[j][k] /= scaleX[j][k];
					errorX[j][k] = errorX[j][k] - (vXimage[j][k] * vXimage[j][k]);
					errorY[j][k] /= scaleY[j][k];
					errorY[j][k] = errorY[j][k] - (vYimage[j][k] * vYimage[j][k]);
					vZimage[j][k] = scaleX[j][k];
				}
				else
				{
					errorX[j][k] = 0.0;
					errorY[j][k] = 0.0;
				}
			}
		}
	}
}
