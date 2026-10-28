#include <math.h>
/*#include "mosaicSource/GeoCodeDEM_p/geocodedem.h"*/
#include "cRecipes/nrutil.h"
#include "common.h"
/*
	Save kernel memory to recyle on repeated calls
*/
double **rDistSave = NULL;
/*
  Compute scale for feathering.
*/
void computeScale(float **inImage, float **scale, int32_t azimuthSize, int32_t rangeSize, float fl, float weight, double minVal)
{
	extern double **rDistSave;
	double **rDist, rA;					  /* Radial distance kernel */
	int32_t j, k, is, il, s1, s2, l1, l2; /* Loop variabiles */
	float minV;
	/* 	No feathring case - set weight to all 1's	*/
	if (fl == 0)
	{
		initFloatMatrix(scale, azimuthSize, rangeSize, 1.0);
		return;
	}
	/*
	  Set up featherin - radial distance from center of kernel.
	  Allocate space if first time, or recycle previoius
	*/
	if (rDistSave == NULL)
	{
		rDistSave = dmatrix(-fl, fl, -fl, fl);
		fillRadialKernel(rDistSave, fl, weight);
	}
	rDist = rDistSave;
	/*
	  Set initial value for scale array with specified weight value
	*/
	initFloatMatrix(scale, (int32_t)azimuthSize, (int32_t)rangeSize, weight);
	/*
	  Now loop through weights
	*/
	for (j = 0; j < azimuthSize; j++)
	{
		for (k = 0; k < rangeSize; k++)
		{
			/*  Adding the feathering at valid points. */
			if (inImage[j][k] > minVal)
			{
				/* Find edge pixels - check all neighbors and if one is non-valid, its an edge */
				minV = 1.0e30; /* will get reduced if its an edge */
				s1 = max(0, j - 1);
				s2 = min(azimuthSize - 1, j + 1);
				l1 = max(0, k - 1);
				l2 = min(rangeSize - 1, k + 1);
				/* find minV for point and its neigbors */
				for (is = s1; is <= s2; is++)
					for (il = l1; il <= l2; il++)
					{
						minV = FMIN(minV, (float)(inImage[is][il]));
					}
				/*
					if minV <= minVal its a border pixel, so set scale to distance from border.
					 Use minvalue where kernels overlap.
				 */
				if ((minV <= minVal))
				{
					s1 = max(0, j - fl);
					s2 = min(azimuthSize - 1, j + fl);
					l1 = max(0, k - fl);
					l2 = min(rangeSize - 1, k + fl);
					/* set the scale to the rad dist, or a less value of already set */
					for (is = s1; is <= s2; is++)
					{
						for (il = l1; il <= l2; il++)
						{
							scale[is][il] = FMIN(rDist[is - j][il - k], scale[is][il]);
						} /* il */
					}	  /* is */
				}		  /* End if(minV... */
			}			  /* End if(inp... */
		}				  /* k */
	}					  /* j */
}

void fillRadialKernel(double **rDist, float fl, float weight)
{
	int is, il; /* Loop variabiles */
	double rA;
	/*    Compute radial distance for a pre-allocated kernel */
	for (is = 0; is <= fl; is++)
	{
		for (il = 0; il <= fl; il++)
		{
			rA = weight * min(max(sqrt((double)(is * is) + (double)(il * il)), 0.5) / fl, 1);
			rDist[il][is] = rA;
			rDist[il][-is] = rA;
			rDist[-il][is] = rA;
			rDist[-il][-is] = rA;
		} /* end for(il... */
	}	  /* end for(is... */
}

/*
  computeScaleEDT -- the same feather ramp as computeScaleLS, from an exact Euclidean
  distance transform instead of stamping a (2fl+1)^2 kernel at every edge pixel.

  computeScaleLS sets, for every pixel, scale = weight * min(max(d, 0.5)/fl, 1), where d is the
  distance to the nearest EDGE pixel (a valid pixel with an invalid 8-neighbour; neighbours off
  the grid do not count, so a sector border is not an edge).  Stamping costs O(nEdge * fl^2) and
  runs on one thread, which made the hopper3D per-product taper roughly double the run time at
  fl = 100.  The Felzenszwalb-Huttenlocher transform gives the same d exactly in O(region), as two
  separable 1-D passes, each parallel over lines.  For integer fl the result is identical to
  computeScaleLS; for fractional fl (e.g. from -flkm) computeScaleLS's window truncates to
  floor(fl) pixels, so edges between floor(fl) and fl away are seen here and not there.

  Only [iMin,iMax) x [jMin,jMax) of scale is written (computeScaleLS also resets rows iMax..end);
  callers must only read that region.
*/
#define EDTINF 1.0e20

static void edt1d(const double *f, int32_t n, double *d, int32_t *v, double *z)
{
	int32_t q, k = 0;
	double s;
	v[0] = 0;
	z[0] = -EDTINF;
	z[1] = EDTINF;
	for (q = 1; q < n; q++)
	{
		s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k])) / (2.0 * q - 2.0 * v[k]);
		while (s <= z[k])
		{
			k--;
			s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k])) / (2.0 * q - 2.0 * v[k]);
		}
		k++;
		v[k] = q;
		z[k] = s;
		z[k + 1] = EDTINF;
	}
	k = 0;
	for (q = 0; q < n; q++)
	{
		while (z[k + 1] < q)
		{
			k++;
		}
		d[q] = (double)(q - v[k]) * (double)(q - v[k]) + f[v[k]];
	}
}

void computeScaleEDT(float **inImage, float **scale, int32_t azimuthSize, int32_t rangeSize, float fl, float weight,
					 double minVal, int32_t iMin, int32_t iMax, int32_t jMin, int32_t jMax)
{
	int32_t nr = iMax - iMin, nc = jMax - jMin, nMax;
	double *d2;
	int32_t i;
	if (nr <= 0 || nc <= 0)
	{
		return;
	}
	nMax = (nr > nc) ? nr : nc;
	d2 = (double *)malloc(sizeof(double) * (size_t)nr * (size_t)nc);
	if (d2 == NULL)
	{
		error("computeScaleEDT: malloc failed (%i x %i)", nr, nc);
	}
	/* Feature map: 0 at edge pixels, "infinity" elsewhere -- the edge test of computeScaleLS */
#pragma omp parallel for schedule(static)
	for (i = 0; i < nr; i++)
	{
		int32_t j, is, il, ii = i + iMin, edge;
		for (j = 0; j < nc; j++)
		{
			int32_t jj = j + jMin;
			edge = FALSE;
			if (inImage[ii][jj] > minVal)
			{
				for (is = max(0, ii - 1); is <= min(azimuthSize - 1, ii + 1) && edge == FALSE; is++)
				{
					for (il = max(0, jj - 1); il <= min(rangeSize - 1, jj + 1); il++)
					{
						if (inImage[is][il] <= minVal)
						{
							edge = TRUE;
							break;
						}
					}
				}
			}
			d2[(size_t)i * nc + j] = (edge == TRUE) ? 0.0 : EDTINF;
		}
	}
	/* Pass 1: down each column.  Pass 2: along each row.  Per-thread line buffers. */
#pragma omp parallel
	{
		double *f = (double *)malloc(sizeof(double) * nMax);
		double *d = (double *)malloc(sizeof(double) * nMax);
		double *z = (double *)malloc(sizeof(double) * (nMax + 1));
		int32_t *v = (int32_t *)malloc(sizeof(int32_t) * nMax);
		int32_t j, k;
		if (f == NULL || d == NULL || z == NULL || v == NULL)
		{
			error("computeScaleEDT: malloc failed for line buffers");
		}
#pragma omp for schedule(static)
		for (j = 0; j < nc; j++)
		{
			for (k = 0; k < nr; k++)
			{
				f[k] = d2[(size_t)k * nc + j];
			}
			edt1d(f, nr, d, v, z);
			for (k = 0; k < nr; k++)
			{
				d2[(size_t)k * nc + j] = d[k];
			}
		}
#pragma omp for schedule(static)
		for (k = 0; k < nr; k++)
		{
			for (j = 0; j < nc; j++)
			{
				f[j] = d2[(size_t)k * nc + j];
			}
			edt1d(f, nc, d, v, z);
			for (j = 0; j < nc; j++)
			{
				d2[(size_t)k * nc + j] = d[j];
			}
		}
		free(f);
		free(d);
		free(z);
		free(v);
	}
	/* The ramp, with exactly fillRadialKernel's arithmetic */
#pragma omp parallel for schedule(static)
	for (i = 0; i < nr; i++)
	{
		int32_t j;
		double rA;
		for (j = 0; j < nc; j++)
		{
			rA = weight * min(max(sqrt(d2[(size_t)i * nc + j]), 0.5) / fl, 1);
			scale[i + iMin][j + jMin] = (float)FMIN(rA, weight);
		}
	}
	free(d2);
}
