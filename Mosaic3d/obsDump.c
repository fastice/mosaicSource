#include "stdio.h"
#include "string.h"
#include "math.h"
#include <stdlib.h>
#include "mosaicSource/common/common.h"
/* MAXOBSDUMPPTS is defined in mosaic3d.h alongside the other mosaic3d limits */
#include "mosaic3d.h"

/*
  Per-observation dump at a small list of points (mosaic3d -obsDump <file>).

  Purpose: let an external tool rebuild mosaic3d's own estimate from its own
  inputs, and then rebuild it again with synthetic observables derived from a GPS
  time series over EACH observation's own window.  Comparing the two cancels
  temporal sampling exactly -- see Documents/gpsForwardModelPlan.md.

  Everything needed for that (a_i, w_i, d_i, t_i, dt_i) exists inside the solver's
  accumulation loop and nowhere else, which is why this is dumped from the C
  rather than reconstructed in Python: reimplementing the weight chain
  (sig2Base, demError, the .sr band) is the one thing guaranteed to diverge
  silently and look like signal.

  Input file: one point per line, "lat lon [name]".  Blank lines and lines
  beginning with # or ; are ignored.
  Output: whitespace-separated rows, one per contributing observation.

  ROW CONVENTION.  The row emitted here is the internal accumulation-buffer
  index, which runs bottom-up (row 0 = southernmost, y increasing with row).
  The GeoTIFF and VRT products are written north-up, so indexing a product
  with this row reads the wrong pixel -- silently, and the error grows with
  distance from the grid centre.  Use tifRow = nRows - 1 - row; col matches.
  Verified 2026-09-02 against the geotransform at ten Greenland GPS sites,
  where the unflipped read was wrong by up to 397 m/yr on fast ice.
*/

static int32_t nPts = 0;
static int32_t *ptRow = NULL, *ptCol = NULL;
static char **ptName = NULL;
static unsigned char *rowWanted = NULL; /* per output row: any point on it? */
static int32_t nRowWanted = 0;
static FILE *dumpFp = NULL;

int32_t obsDumpActive(void) { return (dumpFp != NULL); }

/*
  Read the point list and map each to an output-grid cell.  Called once from
  mosaic3d main AFTER the output grid is defined.
*/
void obsDumpInit(char *pointsFile, char *outName, outputImageStructure *outputImage,
				 double stdLat)
{
	FILE *fp;
	char line[1024], name[256];
	double lat, lon, x, y;
	int32_t i, r, c, n = 0;
	extern double Rotation;

	if (pointsFile == NULL)
		return;
	fp = openInputFile(pointsFile);
	if (fp == NULL)
		error("obsDumpInit: cannot open %s", pointsFile);
	ptRow = (int32_t *)malloc(sizeof(int32_t) * MAXOBSDUMPPTS);
	ptCol = (int32_t *)malloc(sizeof(int32_t) * MAXOBSDUMPPTS);
	ptName = (char **)malloc(sizeof(char *) * MAXOBSDUMPPTS);
	while (fgets(line, sizeof(line), fp) != NULL && n < MAXOBSDUMPPTS)
	{
		if (line[0] == '#' || line[0] == ';' || line[0] == '\n')
			continue;
		name[0] = '\0';
		if (sscanf(line, "%lf %lf %255s", &lat, &lon, name) < 2)
			continue;
		if (name[0] == '\0')
			sprintf(name, "pt%i", n);
		/* lat/lon -> output grid cell.  lltoxy1 returns KILOMETRES, but originX/deltaX are
		   METRES -- mosaicHopper3D.c reads the grid as (originX + jj*deltaX) * MTOKM to get
		   km.  So convert to metres here; treating both as km silently places every point
		   outside the sector. */
		/* stdLat is passed in: outputImage->slat is flagged "not fully implemented"
		   in geocode.h, so the caller supplies the DEM's value instead. */
		lltoxy1(lat, lon, &x, &y, Rotation, stdLat);
		x /= MTOKM;
		y /= MTOKM;
		c = (int32_t)((x - outputImage->originX) / outputImage->deltaX + 0.5);
		r = (int32_t)((y - outputImage->originY) / outputImage->deltaY + 0.5);
		if (r < 0 || c < 0 || r >= outputImage->ySize || c >= outputImage->xSize)
		{
			fprintf(stderr, "obsDumpInit: %s (%.5f %.5f) is outside this sector -- skipped\n",
					name, lat, lon);
			continue;
		}
		ptRow[n] = r;
		ptCol[n] = c;
		ptName[n] = strdup(name);
		fprintf(stderr, "obsDump: %s lat %.5f lon %.5f -> row %i col %i\n", name, lat, lon, r, c);
		n++;
	}
	fclose(fp);
	nPts = n;
	if (nPts == 0)
	{
		fprintf(stderr, "obsDump: no points fall in this sector -- dump disabled\n");
		return;
	}
	/* Cheap per-row gate so the inner loop costs one array lookup */
	nRowWanted = outputImage->ySize;
	rowWanted = (unsigned char *)calloc((size_t)nRowWanted, sizeof(unsigned char));
	for (i = 0; i < nPts; i++)
		rowWanted[ptRow[i]] = 1;
	dumpFp = fopen(outName, "w");
	if (dumpFp == NULL)
		error("obsDumpInit: cannot open %s for writing", outName);
	fprintf(dumpFp, "# mosaic3d -obsDump : one row per contributing observation\n");
	fprintf(dumpFp, "# OBS   point row col obsType frame julDay nDays ax ay az d sigma w\n");
	fprintf(dumpFp, "# PIXEL point row col sx sy\n");
	fprintf(dumpFp, "#   d and sigma are in m/yr (phase already scaled by "
					"365.25/(twok*nDays*sin(psi)))\n");
	fprintf(dumpFp, "#   the 3D solve is v = (sum w a a^T)^-1 (sum w a d)\n");
	fprintf(dumpFp, "#   the surface-parallel (2D) solve projects that with "
					"C = [[1,0],[0,1],[sx,sy]]: N2 = C^T N3 C, b2 = C^T b3.\n");
	fprintf(dumpFp, "#   sx,sy come from the PIXEL row (0,0 on shelf); az is 0 for azimuth rows.\n");
	fprintf(dumpFp, "#   row is the INTERNAL bottom-up index (row 0 = southernmost).  The\n");
	fprintf(dumpFp, "#   GeoTIFF/VRT products are written north-up, so to index them use\n");
	fprintf(dumpFp, "#       tifRow = nRows - 1 - row      (col is the same in both).\n");
}

/*
  Record the per-pixel surface slope used by the surface-parallel projection.
  Without this the dump cannot reproduce a forced-2D solve, because sx/sy come
  from the DEM inside the solver and recomputing them in Python is exactly the
  silent divergence this dump exists to avoid.
*/
void obsDumpPixel(int32_t iRow, int32_t jj, double sx, double sy)
{
	int32_t i;
	if (dumpFp == NULL || iRow < 0 || iRow >= nRowWanted || rowWanted[iRow] == 0)
		return;
	for (i = 0; i < nPts; i++)
	{
		if (ptRow[i] != iRow || ptCol[i] != jj)
			continue;
#pragma omp critical(obsDump)
		{
			fprintf(dumpFp, "PIXEL %s %i %i %.8e %.8e\n", ptName[i], iRow, jj, sx, sy);
			fflush(dumpFp);
		}
		return;
	}
}

/*
  Record one observation.  `iRow`/`jj` are the output-grid cell the solver is
  accumulating into.  Cheap no-op for the overwhelming majority of pixels.
*/
void obsDumpRecord(int32_t iRow, int32_t jj, const char *obsType, const char *frame,
				   double julDay, double nDays, double ax, double ay, double az,
				   double d, double sigma, double w)
{
	int32_t i;
	if (dumpFp == NULL || iRow < 0 || iRow >= nRowWanted || rowWanted[iRow] == 0)
		return;
	for (i = 0; i < nPts; i++)
	{
		if (ptRow[i] != iRow || ptCol[i] != jj)
			continue;
/* Rarely reached, so serialising here costs nothing; the solver's pixel loop
   is OpenMP-partitioned by row and would otherwise interleave writes. */
#pragma omp critical(obsDump)
		{
			fprintf(dumpFp, "OBS %s %i %i %s %s %.6f %.4f %.8e %.8e %.8e %.8e %.8e %.8e\n",
					ptName[i], iRow, jj, obsType, frame != NULL ? frame : "none",
					julDay, nDays, ax, ay, az, d, sigma, w);
			fflush(dumpFp);
		}
		return;
	}
}

void obsDumpClose(void)
{
	if (dumpFp != NULL)
	{
		fclose(dumpFp);
		dumpFp = NULL;
	}
}
