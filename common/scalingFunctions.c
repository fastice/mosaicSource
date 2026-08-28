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
					if (syTmp[i][j] > 0.0f && sxTmp[i][j] > 0.0f)
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
  Correct the crossing-pair error over-counting for one mosaicking round.

  make3DMosaic()/make3DOffsets() accumulate every crossing pair as an independent
  observation, but a pixel covered by n_A ascending and n_D descending images yields
  only n_A + n_D independent measurements, not n_A * n_D.

  How much the naive 1/sum(w) understates the variance depends on HOW CORRELATED the
  per-pair errors are, which is what `rho` parameterises.  Splitting each pair's error
  into a part attached to the contributing images and a part common to every image at
  that pixel, and writing P = n_A*n_D for the pair count:

      f(rho) = rho*P + (1 - rho)*(n_A + n_D)/2

    rho = 0   purely per-image error.  Averaging P pairs built from n_A + n_D images
              gains sqrt((n_A+n_D)/2) less than the naive count claims.  This is the
              right model when random measurement noise dominates.
    rho = 1   error entirely common to every image at the pixel (a DEM error, a slope
              error, an unmodelled vertical rate).  Such a term averages down not at
              all, so the naive estimate is optimistic by the full P.
    rho = 0.5 half the variance non-averaging.  For n_A = n_D = n this equals
              (n + n^2)/2, i.e. EXACTLY the legacy `-pairCountLegacy` formula -- so
              that historical choice was implicitly a rho = 0.5 claim.

  Defaults are set per solution type by the caller (rhoOffsets = 0, rhoPhase = 0.5;
  see mosaic3d.c) because the two observables sit in different regimes: crossing
  offsets are dominated by matching noise, which really does average down, while the
  phase residual is dominated by frame-scale systematics that do not.  Measured
  against the Sentinel-1 reference over stable ground, the rho that reproduces the
  observed standard deviation is ~0.53-0.56 for phase and ~0.08 for offsets.

  The general form above is symmetric in n_A and n_D; the legacy formula
  (nOuter + nPairs)/2 is not, and drifts from it as the two counts diverge.

  Two things this routine is careful about, both of which an earlier in-line version
  of this correction got wrong:

  1) It inflates only THIS round's own contribution to the errorX/errorY
     accumulators (error - error0, error0 being the snapshot taken right after
     undoNormalization), never what earlier rounds already deposited there.  The
     earlier version multiplied the whole buffer after endScale(), so
     undoNormalization() folded the inflation back into the accumulator and the next
     round inflated it a second time -- a phase mosaic's errors were multiplied again
     by the crossing-offsets round's own factor, which is why adding (noisier)
     crossing offsets made the reported errors grow instead of shrink.

  2) It derives n_D from the pair count rather than treating the pair count as n_D.
     nOuter is an exact per-pixel count of distinct outer-loop images (the caller
     dedupes it with a per-image contribution map); nPairs counts contributing pairs.
     With the standard list order (all ascending, then all descending -- see
     consolodateLists() in mosaic3d.c) and sepAscDesc set, the outer-loop
     contributors are exactly the ascending images, so n_D = nPairs / nOuter.  The
     earlier version used nOuter + nPairs, i.e. n_A + n_A*n_D, which for a modest
     3x3 geometry is 12 rather than 6.

  Must be called BEFORE endScale(), while errorX/errorY are still accumulators.
*/
void inflatePairOverCount(outputImageStructure *outputImage, float **errorX, float **errorY,
						  float **errorX0, float **errorY0, float **nOuter, float **nPairs,
						  double rho)
{
	extern int32_t pairCountLegacy;
	double nA, nP, nD, f;
	int32_t j, k;

	/* Guard here rather than at the CLI so every caller is covered, including any
	   that sets the value programmatically.  rho is a variance fraction. */
	if (rho < 0.0)
	{
		rho = 0.0;
	}
	if (rho > 1.0)
	{
		rho = 1.0;
	}
	for (j = 0; j < outputImage->ySize; j++)
	{
		for (k = 0; k < outputImage->xSize; k++)
		{
			nA = (double)nOuter[j][k];
			nP = (double)nPairs[j][k];
			if (nA < 1.0 || nP < 1.0)
			{
				continue;
			}
			if (pairCountLegacy == TRUE)
			{
				/* Legacy: nPairs used directly as if it were the distinct
				   inner-image count.  Retained bit-for-bit so the July-vintage
				   products can be reproduced exactly; it ignores rho, and is
				   the asymmetric relative of f(rho = 0.5). */
				f = 0.5 * (nA + nP);
			}
			else
			{
				nD = nP / nA;
				f = rho * nP + (1.0 - rho) * 0.5 * (nA + nD);
			}
			/* A single pair (n_A = n_D = 1) is already independent -- no inflation. */
			if (f < 1.0)
			{
				f = 1.0;
			}
			errorX[j][k] = errorX0[j][k] + (errorX[j][k] - errorX0[j][k]) * (float)f;
			errorY[j][k] = errorY0[j][k] + (errorY[j][k] - errorY0[j][k]) * (float)f;
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
