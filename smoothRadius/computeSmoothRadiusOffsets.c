#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "smoothRadius.h"

/*
  Per-pixel smoothing-radius map from the MERGED range and azimuth offsets.

  For each candidate half-width r, dr and da are each box-averaged nIter times (box ->
  triangular -> Gaussian-ish) and compared with the unsmoothed value.  A pixel keeps r while
  BOTH components stay within their own tolerance, and is locked the moment either violates --
  the largest contiguous safe radius, the same rule as
  mosaicSource/simInSAR/computeSmoothRadius.c and simoffsets.computeSmoothRadiusMap().

  Two differences from the simInSAR version, both deliberate:
    - it works on the merged dr/da that only simoffsets has (after the velocity and
      vertical-correction terms are added), and locks on either component;
    - the range and azimuth half-widths are capped independently (maxRadiusR, maxRadiusA)
      rather than being tied by the pixel aspect ratio, so an isotropic sweep is just
      maxRadiusR == maxRadiusA.  Offsets grids are resampled to square pixels, so that is
      the normal setting.

  The sweeps look sequential because of the locking rule, but each radius's pass/fail verdict
  depends only on the unsmoothed input: a pixel's answer is (first r that fails) - 1, which is
  assembled afterwards.  So the sweeps run in parallel and are folded in order.  The serial
  form's early exit ("stop once every pixel is locked") is dropped with no practical loss --
  measured on real frames, some pixels stay unlocked all the way to r = 50, so it never fired.
*/

/*
  Masked box mean with half-widths (wr, wa), zero-padded at the edges.  Separable running
  sums, so cost is independent of window size.  Numerator and denominator use the same
  window, so the window-area normalisation cancels and no-data pixels drop out of both.
  A pixel with no valid neighbour keeps its input value.
*/
static void maskedBoxMean(const float *in, const unsigned char *valid, float *out,
						  int32_t nr, int32_t na, int32_t wr, int32_t wa,
						  double *sumBuf, int32_t *cntBuf)
{
	int32_t i, j, add, rem;
	double sum;
	int32_t cnt;
	/* Pass 1: running sum along range within each azimuth row. */
	for (i = 0; i < na; i++)
	{
		const float *row = in + (size_t)i * nr;
		const unsigned char *vrow = valid + (size_t)i * nr;
		double *srow = sumBuf + (size_t)i * nr;
		int32_t *crow = cntBuf + (size_t)i * nr;
		sum = 0.0;
		cnt = 0;
		for (j = 0; j <= wr && j < nr; j++)
		{
			if (vrow[j]) { sum += row[j]; cnt++; }
		}
		srow[0] = sum;
		crow[0] = cnt;
		for (j = 1; j < nr; j++)
		{
			add = j + wr;
			rem = j - wr - 1;
			if (add < nr && vrow[add]) { sum += row[add]; cnt++; }
			if (rem >= 0 && vrow[rem]) { sum -= row[rem]; cnt--; }
			srow[j] = sum;
			crow[j] = cnt;
		}
	}
	/* Pass 2: running sum along azimuth over the range-summed buffers. */
	for (j = 0; j < nr; j++)
	{
		sum = 0.0;
		cnt = 0;
		for (i = 0; i <= wa && i < na; i++)
		{
			sum += sumBuf[(size_t)i * nr + j];
			cnt += cntBuf[(size_t)i * nr + j];
		}
		out[j] = (cnt > 0) ? (float)(sum / cnt) : in[j];
		for (i = 1; i < na; i++)
		{
			add = i + wa;
			rem = i - wa - 1;
			if (add < na) { sum += sumBuf[(size_t)add * nr + j]; cnt += cntBuf[(size_t)add * nr + j]; }
			if (rem >= 0) { sum -= sumBuf[(size_t)rem * nr + j]; cnt -= cntBuf[(size_t)rem * nr + j]; }
			out[(size_t)i * nr + j] = (cnt > 0) ? (float)(sum / cnt) : in[(size_t)i * nr + j];
		}
	}
}

void computeSmoothRadiusOffsets(smrParams *p)
{
	size_t n = (size_t)p->nr * (size_t)p->na;
	size_t q;
	int32_t r, maxRadius;
	unsigned char *valid, *locked, *okAll;

	maxRadius = (p->maxRadiusR > p->maxRadiusA) ? p->maxRadiusR : p->maxRadiusA;
	valid = (unsigned char *)malloc(n);
	locked = (unsigned char *)malloc(n);
	okAll = (unsigned char *)malloc(n * (size_t)maxRadius);
	if (valid == NULL || locked == NULL || okAll == NULL)
	{
		fprintf(stderr, "computeSmoothRadiusOffsets: out of memory for %ld pixels x %d radii\n",
				(long)n, maxRadius);
		exit(1);
	}
	for (q = 0; q < n; q++)
	{
		valid[q] = (p->dr[q] > SMR_NODATA_THRESH && p->da[q] > SMR_NODATA_THRESH &&
					isfinite(p->dr[q]) && isfinite(p->da[q])) ? 1 : 0;
		locked[q] = valid[q] ? 0 : 1;
		p->radius[q] = 0;
	}

#pragma omp parallel for schedule(dynamic, 1) num_threads(p->ompThreads)
	for (r = 1; r <= maxRadius; r++)
	{
		float *sDr, *sDa, *tmp;
		double *sumBuf;
		int32_t *cntBuf, wr, wa, k;
		size_t m;
		unsigned char *ok = okAll + (size_t)(r - 1) * n;
		sDr = (float *)malloc(n * sizeof(float));
		sDa = (float *)malloc(n * sizeof(float));
		tmp = (float *)malloc(n * sizeof(float));
		sumBuf = (double *)malloc(n * sizeof(double));
		cntBuf = (int32_t *)malloc(n * sizeof(int32_t));
		if (sDr == NULL || sDa == NULL || tmp == NULL || sumBuf == NULL || cntBuf == NULL)
		{
			fprintf(stderr, "computeSmoothRadiusOffsets: out of memory in sweep r=%d "
							"(reduce -ompThreads)\n", r);
			exit(1);
		}
		/* Each cap clamps its own axis, so a sweep past one cap keeps widening the other. */
		wr = (r < p->maxRadiusR) ? r : p->maxRadiusR;
		wa = (r < p->maxRadiusA) ? r : p->maxRadiusA;
		memcpy(sDr, p->dr, n * sizeof(float));
		memcpy(sDa, p->da, n * sizeof(float));
		for (k = 0; k < p->nIter; k++)
		{
			maskedBoxMean(sDr, valid, tmp, p->nr, p->na, wr, wa, sumBuf, cntBuf);
			memcpy(sDr, tmp, n * sizeof(float));
			maskedBoxMean(sDa, valid, tmp, p->nr, p->na, wr, wa, sumBuf, cntBuf);
			memcpy(sDa, tmp, n * sizeof(float));
		}
		for (m = 0; m < n; m++)
		{
			ok[m] = (fabsf(sDr[m] - p->dr[m]) <= p->tolDr[m] &&
					 fabsf(sDa[m] - p->da[m]) <= p->tolDa[m]) ? 1 : 0;
		}
		free(sDr); free(sDa); free(tmp); free(sumBuf); free(cntBuf);
	}

	/* Fold in radius order: the first violation locks the pixel at r-1. */
	for (r = 1; r <= maxRadius; r++)
	{
		unsigned char *ok = okAll + (size_t)(r - 1) * n;
		for (q = 0; q < n; q++)
		{
			if (locked[q])
			{
				continue;
			}
			if (ok[q])
			{
				p->radius[q] = (unsigned char)r;
			}
			else
			{
				locked[q] = 1;
			}
		}
	}
	free(valid);
	free(locked);
	free(okAll);
}
