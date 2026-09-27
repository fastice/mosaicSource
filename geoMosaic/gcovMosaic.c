#include "stdio.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <glob.h>
#include <unistd.h>
#include <omp.h>
#include <hdf5.h>
#include "mosaicSource/common/common.h"
#include "geomosaic.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include "gdal.h"
#include "ogr_srs_api.h"

/*
  Support for already geocoded NISAR GCOV (L2) products in geomosaic.

  GCOV products are read directly from the HDF5 file through GDAL's HDF5 driver. GDAL does not
  attach a geotransform or CRS to the GCOV subdatasets, so the grid is reconstructed from the
  xCoordinates/yCoordinates datasets (pixel centres) and the projection epsg_code attribute.

  For each GCOV the covariance term (default HHHH, i.e. RTC gamma0) is block averaged (not
  decimated) to roughly the output resolution, then bilinearly sampled at each output grid point.
  The result is placed in imageTmp/scaleTmp/psiBufTmp/gBufTmp exactly as the range/Doppler loop
  in makeGeoMosaic does, so the normal feathering and geoMosaicScaling accumulation are reused.
*/

#define GCOVINVALID -1.0
#define NEDGE 101

typedef struct
{
	GDALDatasetH hFile;	   /* file-level dataset, for attributes */
	GDALDatasetH hGamma;   /* covariance term, e.g. HHHH */
	GDALDatasetH hFactor;  /* rtcGammaToSigmaFactor */
	GDALDatasetH hMask;	   /* valid sample subswath mask */
	/* Direct HDF5 handles, opened ONLY for float16 bands (slim products). GDAL reports such a
	   band as Float32 and has libhdf5 convert it, which is a generic software path: measured
	   1.273 s against 0.179 s for reading the raw 16-bit values and expanding them here, on
	   the same 16 Mpx window. H5I_INVALID_HID means "not float16, use the GDAL band". */
	hid_t h5Gamma, h5GammaFile;
	hid_t h5Factor, h5FactorFile;
	int32_t nx, ny;
	double x0, y0, dx, dy; /* pixel centre of (0,0) and spacing, in the GCOV CRS */
	int32_t epsg;
} gcovFile;

static void addGCOVFile(gcovInputs *gcov, char *file, float weight);
static void dropSupersededGCOVs(gcovInputs *gcov);
static int32_t openGCOV(gcovInputs *gcov, char *file, gcovFile *g, int32_t needFactor);
static void closeGCOV(gcovFile *g);
static int32_t readCoordinates(char *file, char *dataset, int32_t n, double **coords);
static int32_t getEPSGAttribute(GDALDatasetH hFile, char *key);
static int32_t parseGCOVName(char *file, int32_t *year, int32_t *month, int32_t *day, int32_t *passType);
static int32_t canonicalBandwidth(int32_t mhz);
static int32_t gcovNameBandwidth(char *file, char freq);
static void dropGCOVByBandwidth(gcovInputs *gcov);
static OGRCoordinateTransformationH makeTransform(int32_t fromEPSG, int32_t toEPSG);
static int32_t outputEPSG();
static int32_t transformBox(OGRCoordinateTransformationH ct, double xa, double xb, double ya, double yb,
							double *minX, double *maxX, double *minY, double *maxY);
static float **reduceGCOV(gcovFile *g, int32_t c0, int32_t r0, int32_t nbc, int32_t nbr, int32_t kx, int32_t ky,
						  int32_t useMask, float ***sigma, int32_t maskOnly);
static float **allocFloat2D(int32_t nr, int32_t nc);
static int32_t alignOffset(int32_t start, double anchor, int32_t k);
static void gcovAnchor(gcovFile *g, outputImageStructure *outputImage, double *ax, double *ay);
static void freeFloat2D(float **a);
static int32_t readIncidenceCube(char *file, gcovInputs *gcov, gcovFile *g, float ****cube, double **hgt,
								 int32_t *nh, int32_t *cnx, int32_t *cny, double *cx0, double *cy0, double *cdx, double *cdy);
static float interpIncidence(float ***cube, double *hgt, int32_t nh, int32_t cnx, int32_t cny,
							 double cx0, double cy0, double cdx, double cdy, double X, double Y, double h);

/* ------------------------------------------------------------------------------------------
   float16 fast path.

   Slim GCOV products store the covariance term and the RTC factor as float16. GDAL reports
   such a band as Float32 and asks libhdf5 to convert, which lands in a generic software
   conversion: on a 16 Mpx window, 1.273 s against 0.277 s for a real float32 band. Reading
   the raw 16-bit values and expanding them here takes 0.179 s -- 7x faster than the GDAL
   path, and faster than float32, since it moves 39% fewer bytes off disk.

   A 65536-entry lookup table is used instead of F16C intrinsics so that no architecture
   specific compiler flag is needed and the code stays portable to arm64. The table is built
   once; every GCOV read is serial within a process (the loop in makeGeoMosaic is a plain
   for, and the OpenMP regions here start after the read), so no locking is required.
   ------------------------------------------------------------------------------------------ */
static float f16Table[65536];
static int32_t f16TableReady = FALSE;

static void buildF16Table()
{
	uint32_t h, sign, exp, mant, bits;
	if (f16TableReady == TRUE)
	{
		return;
	}
	for (h = 0; h < 65536; h++)
	{
		sign = (h & 0x8000u) << 16;
		exp = (h >> 10) & 0x1fu;
		mant = h & 0x3ffu;
		if (exp == 0)
		{
			if (mant == 0)
			{
				bits = sign; /* +-0 */
			}
			else
			{ /* subnormal: renormalise into a float32 normal */
				exp = 127 - 15 + 1;
				while ((mant & 0x400u) == 0)
				{
					mant <<= 1;
					exp--;
				}
				mant &= 0x3ffu;
				bits = sign | (exp << 23) | (mant << 13);
			}
		}
		else if (exp == 31)
		{
			bits = sign | 0x7f800000u | (mant << 13); /* inf / NaN */
		}
		else
		{
			bits = sign | ((exp - 15 + 127) << 23) | (mant << 13);
		}
		memcpy(&f16Table[h], &bits, sizeof(float));
	}
	f16TableReady = TRUE;
}

/*
  Open dataset `path` in `file` directly with HDF5, but only if it is 16-bit float. Returns
  TRUE and sets *fileId/*dsetId in that case; otherwise leaves them invalid and the caller
  keeps using GDAL.
*/
static int32_t openIfFloat16(char *file, char *path, hid_t *fileId, hid_t *dsetId)
{
	hid_t f, d, t;
	int32_t isF16;
	*fileId = *dsetId = H5I_INVALID_HID;
	H5Eset_auto2(H5E_DEFAULT, NULL, NULL); /* probing: a miss is not an error */
	if ((f = H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
	{
		return FALSE;
	}
	if ((d = H5Dopen2(f, path, H5P_DEFAULT)) < 0)
	{
		H5Fclose(f);
		return FALSE;
	}
	t = H5Dget_type(d);
	isF16 = (H5Tget_class(t) == H5T_FLOAT && H5Tget_size(t) == 2) ? TRUE : FALSE;
	H5Tclose(t);
	if (isF16 == FALSE)
	{
		H5Dclose(d);
		H5Fclose(f);
		return FALSE;
	}
	buildF16Table();
	*fileId = f;
	*dsetId = d;
	return TRUE;
}

/*
  Read a window of a float16 dataset into float32. The read uses the FILE's own datatype as
  the memory type, so HDF5 copies the bits rather than converting them -- that conversion is
  the whole cost being avoided.
*/
static int32_t readF16Window(hid_t dset, int32_t x, int32_t y, int32_t w, int32_t h,
							 float *out, unsigned short *tmp)
{
	hid_t fs, ms, ft;
	hsize_t off[2], cnt[2];
	size_t i, n = (size_t)w * (size_t)h;
	int32_t ok;
	off[0] = (hsize_t)y;
	off[1] = (hsize_t)x;
	cnt[0] = (hsize_t)h;
	cnt[1] = (hsize_t)w;
	fs = H5Dget_space(dset);
	if (H5Sselect_hyperslab(fs, H5S_SELECT_SET, off, NULL, cnt, NULL) < 0)
	{
		H5Sclose(fs);
		return FALSE;
	}
	ms = H5Screate_simple(2, cnt, NULL);
	ft = H5Dget_type(dset);
	ok = (H5Dread(dset, ft, ms, fs, H5P_DEFAULT, tmp) >= 0) ? TRUE : FALSE;
	H5Tclose(ft);
	H5Sclose(ms);
	H5Sclose(fs);
	if (ok == TRUE)
	{
		for (i = 0; i < n; i++)
		{
			out[i] = f16Table[tmp[i]];
		}
	}
	return ok;
}

/*
  Read one window of a band as float32, taking the float16 path when the band has one.
*/
static int32_t readGCOVWindow(hid_t h5, GDALRasterBandH band, int32_t x, int32_t y,
							  int32_t w, int32_t h, float *out, unsigned short *tmp)
{
	if (h5 != H5I_INVALID_HID)
	{
		return readF16Window(h5, x, y, w, h, out, tmp);
	}
	return (GDALRasterIO(band, GF_Read, x, y, w, h, out, w, h, GDT_Float32, 0, 0) == CE_None)
			   ? TRUE : FALSE;
}

/*
  Read the yaml file specifying the GCOV inputs. Hand parsed (no libyaml), same approach as
  getBaseline.c. Format:

	polarization: HH        # HH -> HHHH, HV -> HVHV, or a full term such as HHHH
	frequency: A           # A (default) or B; B is the ionosphere band
	bandwidth: 5, 20, 40, 80  # optional MHz list for the SELECTED frequency;
	                          # omitted = all. 80 and 77 mean the same mode.
	useMask: true           # drop samples flagged invalid (0) or fill (255) in the GCOV mask
	glob: /path/to/*.h5     # optional, may repeat; matches are sorted
	files:                  # optional list, "- path [weight]"
	  - /path/a.h5
	  - /path/b.h5 0.5
*/
void readGCOVYaml(char *yamlFile, gcovInputs *gcov)
{
	FILE *fp;
	char line[2048], key[256], value[2048], *c, *p;
	float weight;
	int32_t inFiles, i;
	glob_t globResult;

	strcpy(gcov->polarization, "HHHH");
	gcov->nBandwidths = 0; /* no filter: every bandwidth accepted */
	strcpy(gcov->frequency, "A");
	gcov->useMask = TRUE;
	gcov->factorDir[0] = '\0';
	gcov->nFiles = 0;
	gcov->files = NULL;
	gcov->weights = NULL;
	fp = openInputFile(yamlFile);
	inFiles = FALSE;
	while (fgets(line, sizeof(line), fp) != NULL)
	{
		/* Strip comments and trailing newline */
		if ((c = strchr(line, '#')) != NULL)
		{
			*c = '\0';
		}
		line[strcspn(line, "\r\n")] = '\0';
		for (p = line; *p == ' ' || *p == '\t'; p++)
		{
		}
		if (*p == '\0')
		{
			continue;
		}
		/* List item under files: */
		if (*p == '-')
		{
			if (inFiles == FALSE)
			{
				error("readGCOVYaml: list item outside of files: in %s\n%s", yamlFile, line);
			}
			weight = 1.0;
			if (sscanf(p + 1, "%s %f", value, &weight) < 1)
			{
				error("readGCOVYaml: missing file name in %s\n%s", yamlFile, line);
			}
			addGCOVFile(gcov, value, weight);
			continue;
		}
		inFiles = FALSE;
		if ((c = strchr(p, ':')) == NULL)
		{
			error("readGCOVYaml: cannot parse line in %s\n%s", yamlFile, line);
		}
		*c = '\0';
		if (sscanf(p, "%255s", key) != 1)
		{
			error("readGCOVYaml: missing key in %s", yamlFile);
		}
		value[0] = '\0';
		sscanf(c + 1, "%2047s", value);
		if (strcmp(key, "files") == 0)
		{
			inFiles = TRUE;
		}
		else if (strcmp(key, "bandwidth") == 0)
		{
			/* A list on one line: "bandwidth: 5, 20, 40, 80". Parsed from the rest of the
			   line rather than from value, which holds only the first token. */
			char *q = c + 1;
			gcov->nBandwidths = 0;
			while (*q != '\0')
			{
				if (isdigit((unsigned char)*q) == 0)
				{
					q++;
					continue;
				}
				if (gcov->nBandwidths >= MAXGCOVBANDWIDTHS)
				{
					error("readGCOVYaml: more than %i bandwidths in %s", MAXGCOVBANDWIDTHS, yamlFile);
				}
				gcov->bandwidths[gcov->nBandwidths] = canonicalBandwidth(atoi(q));
				gcov->nBandwidths++;
				while (isdigit((unsigned char)*q) != 0)
				{
					q++;
				}
			}
			if (gcov->nBandwidths == 0)
			{
				error("readGCOVYaml: bandwidth: given but no values in %s", yamlFile);
			}
		}
		else if (strcmp(key, "polarization") == 0)
		{
			for (i = 0; i < (int32_t)strlen(value); i++)
			{
				value[i] = toupper(value[i]);
			}
			if (strlen(value) == 2)
			{
				sprintf(gcov->polarization, "%s%s", value, value);
			}
			else if (strlen(value) == 4)
			{
				strcpy(gcov->polarization, value);
			}
			else
			{
				error("readGCOVYaml: invalid polarization %s", value);
			}
		}
		else if (strcmp(key, "frequency") == 0)
		{
			if (strlen(value) != 1)
			{
				error("readGCOVYaml: invalid frequency %s", value);
			}
			gcov->frequency[0] = toupper(value[0]);
			gcov->frequency[1] = '\0';
		}
		else if (strcmp(key, "useMask") == 0)
		{
			gcov->useMask = (strcasecmp(value, "true") == 0 || strcmp(value, "1") == 0) ? TRUE : FALSE;
		}
		else if (strcmp(key, "factorFrom") == 0)
		{
			if (strlen(value) >= sizeof(gcov->factorDir))
			{
				error("readGCOVYaml: factorFrom path too long: %s", value);
			}
			strcpy(gcov->factorDir, value);
		}
		else if (strcmp(key, "glob") == 0)
		{
			if (glob(value, 0, NULL, &globResult) == 0)
			{
				for (i = 0; i < (int32_t)globResult.gl_pathc; i++)
				{
					addGCOVFile(gcov, globResult.gl_pathv[i], 1.0);
				}
			}
			else
			{
				fprintf(stderr, "readGCOVYaml: warning, no files match glob %s\n", value);
			}
			globfree(&globResult);
		}
		else
		{
			error("readGCOVYaml: unknown key %s in %s", key, yamlFile);
		}
	}
	fclose(fp);
	dropGCOVByBandwidth(gcov);
	dropSupersededGCOVs(gcov);
	fprintf(stderr, "GCOV inputs: %i files, frequency%s/%s, useMask %i\n",
			gcov->nFiles, gcov->frequency, gcov->polarization, gcov->useMask);
	fprintf(stderr, "GCOV bandwidth filter: ");
	if (gcov->nBandwidths == 0)
	{
		fprintf(stderr, "none (all bandwidths)\n");
	}
	else
	{
		for (i = 0; i < gcov->nBandwidths; i++)
		{
			fprintf(stderr, "%i MHz%s", gcov->bandwidths[i],
					(i < gcov->nBandwidths - 1) ? ", " : "\n");
		}
	}
	for (i = 0; i < gcov->nFiles; i++)
	{
		fprintf(stderr, "  GCOV %i: %s %f\n", i + 1, gcov->files[i], gcov->weights[i]);
	}
}

/*
  The 77 MHz mode is called "80 MHz" in most mission documents and either spelling turns up,
  in file names and in hand-written yaml alike. Fold them together so the comparison cannot
  depend on which name the writer happened to use.
*/
static int32_t canonicalBandwidth(int32_t mhz)
{
	return (mhz == 80) ? 77 : mhz;
}

/*
  Bandwidth in MHz of frequency `freq` ('A' or 'B') taken from the granule NAME: field 8
  (0-based) is a four-character pair "AABB" giving each band's bandwidth, e.g. 4005 = A 40 MHz
  + B 5 MHz, 7700 = A 77 MHz + B absent. Returns -1 if the name does not parse and 0 if that
  frequency is absent.

  Deliberately read from the name rather than from
  metadata/sourceData/swaths/frequency<X>/acquiredRangeBandwidth (which is authoritative, and
  agrees -- checked on both modes): the point of the filter is to reject a granule without
  opening it, which for a remote input means without a network round trip.
*/
static int32_t gcovNameBandwidth(char *file, char freq)
{
	char name[1024], *base, *tok, *fields[20];
	int32_t n, i, mhz;
	base = strrchr(file, '/');
	strncpy(name, (base == NULL) ? file : base + 1, sizeof(name) - 1);
	name[sizeof(name) - 1] = '\0';
	n = 0;
	for (tok = strtok(name, "_"); tok != NULL && n < 20; tok = strtok(NULL, "_"))
	{
		fields[n++] = tok;
	}
	if (n < 9 || strlen(fields[8]) != 4)
	{
		return -1;
	}
	for (i = 0; i < 4; i++)
	{
		if (isdigit((unsigned char)fields[8][i]) == 0)
		{
			return -1;
		}
	}
	i = (freq == 'B') ? 2 : 0;
	mhz = (fields[8][i] - '0') * 10 + (fields[8][i + 1] - '0');
	return canonicalBandwidth(mhz);
}

/*
  Apply the optional bandwidth filter. A no-op unless the yaml gave a bandwidth: list, so
  mixed bandwidths are the default.
*/
static void dropGCOVByBandwidth(gcovInputs *gcov)
{
	int32_t i, k, n, keep, bw, nBefore, nBad;
	if (gcov->nBandwidths == 0)
	{
		return;
	}
	nBefore = gcov->nFiles;
	nBad = 0;
	n = 0;
	for (i = 0; i < gcov->nFiles; i++)
	{
		bw = gcovNameBandwidth(gcov->files[i], gcov->frequency[0]);
		keep = FALSE;
		if (bw < 0)
		{
			fprintf(stderr, "GCOV: cannot read bandwidth from name %s -- dropped by the "
							"bandwidth filter\n", gcov->files[i]);
			nBad++;
		}
		else
		{
			for (k = 0; k < gcov->nBandwidths; k++)
			{
				if (gcov->bandwidths[k] == bw)
				{
					keep = TRUE;
				}
			}
		}
		if (keep == TRUE)
		{
			gcov->files[n] = gcov->files[i];
			gcov->weights[n] = gcov->weights[i];
			n++;
		}
		else
		{
			if (bw >= 0)
			{
				fprintf(stderr, "GCOV: dropping %s (frequency%s bandwidth %i MHz)\n",
						gcov->files[i], gcov->frequency, bw);
			}
			free(gcov->files[i]);
		}
	}
	gcov->nFiles = n;
	/* Every granule failing the filter is far more likely to be a wrong bandwidth list or a
	   changed naming convention than a real empty tile, and an empty GCOV list otherwise
	   produces an all-nodata mosaic that exits 0. Same reasoning as the remote-open failure. */
	if (nBefore > 0 && gcov->nFiles == 0)
	{
		error("readGCOVYaml: the bandwidth filter rejected all %i granules (%i had an\n"
			  "  unparseable name). Check the bandwidth: list against frequency%s.",
			  nBefore, nBad, gcov->frequency);
	}
}

static void addGCOVFile(gcovInputs *gcov, char *file, float weight)
{
	char control[2048];
	/* a download still in progress (aria2 control file alongside) may open but be incomplete */
	snprintf(control, sizeof(control), "%s.aria2", file);
	if (access(control, F_OK) == 0)
	{
		fprintf(stderr, "GCOV: skipping %s (download in progress: %s exists)\n", file, control);
		return;
	}
	gcov->files = (char **)realloc(gcov->files, sizeof(char *) * (gcov->nFiles + 1));
	gcov->weights = (float *)realloc(gcov->weights, sizeof(float) * (gcov->nFiles + 1));
	gcov->files[gcov->nFiles] = strdup(file);
	gcov->weights[gcov->nFiles] = weight;
	gcov->nFiles++;
}

/*
  ASF can deliver the same granule more than once with an incremented product counter (the final
  _NNN field, e.g. _001 and a reprocessed _002). Keep only the highest counter so a frame is not
  mosaicked twice. Names are compared on the basename minus the counter.
*/
static void dropSupersededGCOVs(gcovInputs *gcov)
{
	char stemI[1024], stemJ[1024], *base, *u;
	int32_t i, j, n, keep;
	n = 0;
	for (i = 0; i < gcov->nFiles; i++)
	{
		base = strrchr(gcov->files[i], '/');
		strncpy(stemI, (base == NULL) ? gcov->files[i] : base + 1, sizeof(stemI) - 1);
		stemI[sizeof(stemI) - 1] = '\0';
		if ((u = strrchr(stemI, '_')) != NULL)
		{
			*u = '\0';
		}
		keep = TRUE;
		for (j = 0; j < gcov->nFiles; j++)
		{
			if (j == i)
			{
				continue;
			}
			base = strrchr(gcov->files[j], '/');
			strncpy(stemJ, (base == NULL) ? gcov->files[j] : base + 1, sizeof(stemJ) - 1);
			stemJ[sizeof(stemJ) - 1] = '\0';
			if ((u = strrchr(stemJ, '_')) != NULL)
			{
				*u = '\0';
			}
			/* same granule: drop i if j has a higher counter, or the same name listed earlier */
			if (strcmp(stemI, stemJ) == 0)
			{
				base = strrchr(gcov->files[i], '_');
				u = strrchr(gcov->files[j], '_');
				if (strcmp(u, base) > 0 || (strcmp(u, base) == 0 && j < i))
				{
					keep = FALSE;
				}
			}
		}
		if (keep == TRUE)
		{
			gcov->files[n] = gcov->files[i];
			gcov->weights[n] = gcov->weights[i];
			n++;
		}
		else
		{
			fprintf(stderr, "GCOV: dropping superseded/duplicate %s\n", gcov->files[i]);
		}
	}
	gcov->nFiles = n;
}

/*
  Expand minX..maxY (km, output projection) to include the footprints of the GCOV rasters.
  Used to autosize the output grid.
*/
void gcovBounds(gcovInputs *gcov, double *minX, double *maxX, double *minY, double *maxY)
{
	gcovFile g;
	OGRCoordinateTransformationH ct;
	double x1, x2, y1, y2;
	int32_t i;
	if (gcov == NULL)
	{
		return;
	}
	for (i = 0; i < gcov->nFiles; i++)
	{
		if (openGCOV(gcov, gcov->files[i], &g, FALSE) == FALSE)
		{
			continue;
		}
		ct = makeTransform(g.epsg, outputEPSG());
		if (transformBox(ct, g.x0 - 0.5 * g.dx, g.x0 + (g.nx - 0.5) * g.dx,
						 g.y0 - 0.5 * g.dy, g.y0 + (g.ny - 0.5) * g.dy, &x1, &x2, &y1, &y2) == TRUE)
		{
			*minX = min(*minX, x1 * MTOKM);
			*maxX = max(*maxX, x2 * MTOKM);
			*minY = min(*minY, y1 * MTOKM);
			*maxY = max(*maxY, y2 * MTOKM);
		}
		OCTDestroyCoordinateTransformation(ct);
		closeGCOV(&g);
	}
}

/*
  Sample GCOV file number iFile onto the output grid. Fills imageTmp, scaleTmp, psiBufTmp and gBufTmp
  over [iMin,iMax) x [jMin,jMax), and gcovImage (weight, date, passType) for geoMosaicScaling.
  Returns FALSE if the file is skipped (missing, out of date range, no overlap, zero weight).
*/
/*
  Pass 1 of the range selection over the GCOV inputs: fold each frame's ellipsoidal incidence
  angle into the coarse extremum buffer. Runs the same gcovToOutputGrid setup as pass 2 (region,
  block alignment, skips) so the two passes cannot drift apart, but reads only the mask and the
  incidence cube rather than HHHH and rtcGammaToSigmaFactor.
*/
void gcovIncidenceCoarse(gcovInputs *gcov, outputImageStructure *outputImage, void *dem, incBuffer *incBuf)
{
	inputImageStructure gcovImage;
	int32_t i, imageDate, iMin, iMax, jMin, jMax;
	if (gcov == NULL)
	{
		return;
	}
	for (i = 0; i < gcov->nFiles; i++)
	{
		gcovToOutputGrid(gcov, i, outputImage, dem, NULL, NULL, NULL, NULL,
						 &gcovImage, &imageDate, &iMin, &iMax, &jMin, &jMax, incBuf, NULL, 1);
	}
}

int32_t gcovToOutputGrid(gcovInputs *gcov, int32_t iFile, outputImageStructure *outputImage, void *dem,
						 float **imageTmp, float **scaleTmp, float **psiBufTmp, float **gBufTmp,
						 inputImageStructure *gcovImage, int32_t *imageDate,
						 int32_t *iMin, int32_t *iMax, int32_t *jMin, int32_t *jMax,
						 incBuffer *incBuf, unsigned char **selTmp, int32_t pass)
{
	extern int32_t rangeSelect;
	int32_t padC, padR;
	extern int32_t S1Cal;
	extern int32_t calOutput;
	extern int32_t gcovOnly;
	int32_t doy[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 333};
	gcovFile g;
	OGRCoordinateTransformationH ct, *ctThread;
	float **gamma, **sigma, ***cube;
	double *hgt, cx0, cy0, cdx, cdy;
	double x1, x2, y1, y2, xc, yc, X[3], Y[3], Z[3];
	double spanC, spanR, jd, tRead, tSample, anchorX, anchorY;
	int32_t year, month, day, passType, needSigma, doPsi, needCube;
	int32_t c0, c1, r0, r1, kx, ky, nbc, nbr, nh, cnx, cny, nThreads, t;
	char *file;

	file = gcov->files[iFile];
	fprintf(stderr, "GCOV %s\n", file);
	*iMin = 0;
	*iMax = 0;
	*jMin = 0;
	*jMax = 0;
	if (gcov->weights[iFile] < 1e-20)
	{
		fprintf(stderr, "Skip (0 weight)\n");
		return FALSE;
	}
	/* Date and pass direction from the NISAR file name */
	if (parseGCOVName(file, &year, &month, &day, &passType) == FALSE)
	{
		error("gcovToOutputGrid: cannot parse date/direction from NISAR name %s", file);
	}
	jd = juldayDouble(month, day, year);
	if (jd < outputImage->jd1 || jd > outputImage->jd2)
	{
		fprintf(stderr, "Skip (outside date range)\n");
		return FALSE;
	}
	memset(gcovImage, 0, sizeof(inputImageStructure));
	gcovImage->file = file;
	gcovImage->weight = gcov->weights[iFile];
	gcovImage->passType = passType;
	gcovImage->year = year;
	gcovImage->month = month;
	gcovImage->day = day;
	gcovImage->noData = 0.0;
	*imageDate = year * 365. + doy[month - 1] + day;
	/* sigma0 needed for calibrated output or if requested for uncalibrated */
	needSigma = ((S1Cal & TRUE) == TRUE || calOutput == CALOUTPUT_SIGMA0) ? TRUE : FALSE;
	/* gamma0 is what the GCOV stores, so when it is the only output wanted AND GCOVs are the
	   only inputs, sigma0 is dead work: value carries gamma0 directly and the gBuf offset is
	   0, giving the identical gamma0 at output (image holds dB, gamma[][] is added to it).
	   That skips the whole rtcGammaToSigmaFactor band -- an extra read the size of the data
	   band, which dominates when the granule is remote.

	   Conditional on gcovOnly because under -S1Cal the range/Doppler images put SIGMA0 in the
	   same accumulator; mixing gamma0 from GCOVs into it would average two different
	   quantities. */
	if (calOutput == CALOUTPUT_GAMMA0 && gcovOnly == TRUE)
	{
		needSigma = FALSE;
	}
	doPsi = ((S1Cal & TRUE) == TRUE && (S1Cal & PSISAVE) > 0) ? TRUE : FALSE;
	/* the range selection needs the incidence cube whether or not psi is being written out */
	needCube = (doPsi == TRUE || rangeSelect != RANGESELECT_NONE) ? TRUE : FALSE;
	if (pass == 1)
	{
		/* footprint and incidence only: no value, no sigma, nothing written to the tmp buffers */
		needSigma = FALSE;
		doPsi = FALSE;
		needCube = TRUE;
	}
	if (openGCOV(gcov, file, &g, needSigma) == FALSE)
	{
		fprintf(stderr, "Skip (could not open)\n");
		return FALSE;
	}
	/*
	  Output region covered by the GCOV raster
	*/
	ct = makeTransform(g.epsg, outputEPSG());
	if (transformBox(ct, g.x0 - 0.5 * g.dx, g.x0 + (g.nx - 0.5) * g.dx,
					 g.y0 - 0.5 * g.dy, g.y0 + (g.ny - 0.5) * g.dy, &x1, &x2, &y1, &y2) == FALSE)
	{
		error("gcovToOutputGrid: could not transform bounds for %s", file);
	}
	OCTDestroyCoordinateTransformation(ct);
	*jMin = max(0, (int32_t)floor((x1 - outputImage->originX) / outputImage->deltaX));
	*jMax = min(outputImage->xSize, (int32_t)ceil((x2 - outputImage->originX) / outputImage->deltaX) + 1);
	*iMin = max(0, (int32_t)floor((y1 - outputImage->originY) / outputImage->deltaY));
	*iMax = min(outputImage->ySize, (int32_t)ceil((y2 - outputImage->originY) / outputImage->deltaY) + 1);
	if (*jMax <= *jMin || *iMax <= *iMin)
	{
		fprintf(stderr, "Skip (out of bounds)\n");
		closeGCOV(&g);
		*iMin = *iMax = *jMin = *jMax = 0;
		return FALSE;
	}
	/*
	  Source window and averaging factors
	*/
	ct = makeTransform(outputEPSG(), g.epsg);
	transformBox(ct, outputImage->originX + *jMin * outputImage->deltaX, outputImage->originX + (*jMax - 1) * outputImage->deltaX,
				 outputImage->originY + *iMin * outputImage->deltaY, outputImage->originY + (*iMax - 1) * outputImage->deltaY,
				 &x1, &x2, &y1, &y2);
	/* GCOV pixels spanned by one output pixel, at the centre of the region */
	/* an output pixel centre near the middle of the region (also the block-alignment anchor) */
	xc = outputImage->originX + ((*jMin + *jMax) / 2) * outputImage->deltaX;
	yc = outputImage->originY + ((*iMin + *iMax) / 2) * outputImage->deltaY;
	X[0] = xc;
	Y[0] = yc;
	X[1] = xc + outputImage->deltaX;
	Y[1] = yc;
	X[2] = xc;
	Y[2] = yc + outputImage->deltaY;
	Z[0] = Z[1] = Z[2] = 0.0;
	OCTTransform(ct, 3, X, Y, Z);
	OCTDestroyCoordinateTransformation(ct);
	spanC = sqrt(pow((X[1] - X[0]) / g.dx, 2) + pow((X[2] - X[0]) / g.dx, 2));
	spanR = sqrt(pow((Y[1] - Y[0]) / g.dy, 2) + pow((Y[2] - Y[0]) / g.dy, 2));
	kx = max(1, (int32_t)floor(spanC));
	ky = max(1, (int32_t)floor(spanR));
	/* window in source pixels, padded, clipped, and trimmed to a multiple of kx, ky */
	padC = 2 * kx + 2;
	padR = 2 * ky + 2;
	if (pass == 1)
	{
		/* pass 1 samples cell CENTRES, which can sit up to half a cell outside the region; the
		   window has to cover them or a cell's value would depend on how the mosaic is tiled */
		padC += (incBuf->stride / 2 + 1) * kx;
		padR += (incBuf->stride / 2 + 1) * ky;
	}
	c0 = (int32_t)floor(min((x1 - g.x0) / g.dx, (x2 - g.x0) / g.dx)) - padC;
	c1 = (int32_t)ceil(max((x1 - g.x0) / g.dx, (x2 - g.x0) / g.dx)) + padC;
	r0 = (int32_t)floor(min((y1 - g.y0) / g.dy, (y2 - g.y0) / g.dy)) - padR;
	r1 = (int32_t)ceil(max((y1 - g.y0) / g.dy, (y2 - g.y0) / g.dy)) + padR;
	c0 = max(0, c0);
	r0 = max(0, r0);
	c1 = min(g.nx, c1);
	r1 = min(g.ny, r1);
	/*
	  Align the averaging blocks so a block centre falls on an output pixel centre. The anchor
	  is the output-lattice point nearest the GCOV origin, which depends only on the GCOV and
	  the output lattice (not on this tile's extent), so every tile uses the same phase and
	  results do not depend on how a mosaic is tiled -- including when the output spacing is
	  not a whole number of GCOV pixels (e.g. 25 m from 10 m). When it is (same CRS), each
	  output pixel is exactly the box average of the GCOV pixels inside it.
	*/
	gcovAnchor(&g, outputImage, &anchorX, &anchorY);
	c0 += alignOffset(c0, (anchorX - g.x0) / g.dx - 0.5 * (kx - 1), kx);
	r0 += alignOffset(r0, (anchorY - g.y0) / g.dy - 0.5 * (ky - 1), ky);
	nbc = (c1 - c0) / kx;
	nbr = (r1 - r0) / ky;
	fprintf(stderr, "GCOV window cols %i-%i rows %i-%i, averaging %i x %i -> %i x %i\n",
			c0, c1, r0, r1, kx, ky, nbc, nbr);
	if (nbc < 2 || nbr < 2)
	{
		fprintf(stderr, "Skip (no overlap)\n");
		closeGCOV(&g);
		*iMin = *iMax = *jMin = *jMax = 0;
		return FALSE;
	}
	tRead = omp_get_wtime();
	gamma = reduceGCOV(&g, c0, r0, nbc, nbr, kx, ky, gcov->useMask, needSigma ? &sigma : NULL,
					   (pass == 1) ? TRUE : FALSE);
	tRead = omp_get_wtime() - tRead;
	if (needCube == TRUE)
	{
		if (readIncidenceCube(file, gcov, &g, &cube, &hgt, &nh, &cnx, &cny, &cx0, &cy0, &cdx, &cdy) == FALSE)
		{
			if (rangeSelect != RANGESELECT_NONE)
			{
				/* silently mosaicking it unfiltered would quietly corrupt the selection */
				error("no incidence angle cube in %s, required by nearRange/farRange", file);
			}
			fprintf(stderr, "Warning: no incidence angle cube for %s, psi set to 0\n", file);
			doPsi = FALSE;
			needCube = FALSE;
		}
	}
	closeGCOV(&g);
	/*
	  One coordinate transform per thread (OGR transforms are not thread safe)
	*/
	nThreads = omp_get_max_threads();
	ctThread = (OGRCoordinateTransformationH *)malloc(sizeof(OGRCoordinateTransformationH) * nThreads);
	for (t = 0; t < nThreads; t++)
	{
		ctThread[t] = makeTransform(outputEPSG(), g.epsg);
	}
	if (pass == 1)
	{
		int32_t cMinX, cMaxX, cMinY, cMaxY;
		incCellSpan(incBuf, TRUE, *jMin, *jMax - 1, &cMinX, &cMaxX);
		incCellSpan(incBuf, FALSE, *iMin, *iMax - 1, &cMinY, &cMaxY);
		nThreads = omp_get_max_threads();
		ctThread = (OGRCoordinateTransformationH *)malloc(sizeof(OGRCoordinateTransformationH) * nThreads);
		for (t = 0; t < nThreads; t++)
		{
			ctThread[t] = makeTransform(outputEPSG(), g.epsg);
		}
#pragma omp parallel
		{
			OGRCoordinateTransformationH myCt = ctThread[omp_get_thread_num()];
			int32_t ci, cj, i1, j1;
			double xw, yw, zw, col, row, h;
			/* one thread per CELL ROW, so incUpdate touches a row no other thread touches */
#pragma omp for schedule(static)
			for (ci = cMinY; ci <= cMaxY; ci++)
			{
				i1 = incCellCentre(incBuf, FALSE, ci);
				for (cj = cMinX; cj <= cMaxX; cj++)
				{
					j1 = incCellCentre(incBuf, TRUE, cj);
					xw = outputImage->originX + j1 * outputImage->deltaX;
					yw = outputImage->originY + i1 * outputImage->deltaY;
					zw = 0.0;
					OCTTransform(myCt, 1, &xw, &yw, &zw);
					col = ((xw - g.x0) / g.dx - c0 - 0.5 * (kx - 1)) / kx;
					row = ((yw - g.y0) / g.dy - r0 - 0.5 * (ky - 1)) / ky;
					if (bilinearInterp(gamma, col, row, nbc, nbr, 0.0, GCOVINVALID) <= 0)
					{
						continue;
					}
					h = getXYHeightXY((outputImage->originX + j1 * outputImage->deltaX) * MTOKM,
									  (outputImage->originY + i1 * outputImage->deltaY) * MTOKM, 0.0,
									  (xyDEM *)dem, 0.0, ELLIPSOIDAL);
					incUpdate(incBuf, i1, j1, interpIncidence(cube, hgt, nh, cnx, cny, cx0, cy0, cdx, cdy,
															  xw, yw, h));
				}
			}
		}
		for (t = 0; t < nThreads; t++)
		{
			OCTDestroyCoordinateTransformation(ctThread[t]);
		}
		free(ctThread);
		freeFloat2D(gamma);
		for (t = 0; t < nh; t++)
		{
			freeFloat2D(cube[t]);
		}
		free(cube);
		free(hgt);
		return TRUE;
	}
	/*
	  Invalidate a 1-pixel ring around the region. computeScale (feathering) scans the whole
	  buffer, so stale values left by earlier inputs would otherwise hide a data edge that
	  coincides with the region boundary (e.g. frame ends along track).
	*/
	{
		int32_t i1, j1;
		for (i1 = max(*iMin - 1, 0); i1 < min(*iMax + 1, outputImage->ySize); i1++)
		{
			for (j1 = max(*jMin - 1, 0); j1 < min(*jMax + 1, outputImage->xSize); j1++)
			{
				if (i1 < *iMin || i1 >= *iMax || j1 < *jMin || j1 >= *jMax)
				{
					imageTmp[i1][j1] = -LARGEINT;
				}
			}
		}
	}
	tSample = omp_get_wtime();
#pragma omp parallel
	{
		int32_t i1, j1, n;
		double *xr, *yr, *zr, col, row, h;
		float gam, sig, value, psi;
		OGRCoordinateTransformationH myCt = ctThread[omp_get_thread_num()];
		n = *jMax - *jMin;
		xr = (double *)malloc(sizeof(double) * n);
		yr = (double *)malloc(sizeof(double) * n);
		zr = (double *)malloc(sizeof(double) * n);
#pragma omp for schedule(dynamic, 8)
		for (i1 = *iMin; i1 < *iMax; i1++)
		{
			for (j1 = *jMin; j1 < *jMax; j1++)
			{
				xr[j1 - *jMin] = outputImage->originX + j1 * outputImage->deltaX;
				yr[j1 - *jMin] = outputImage->originY + i1 * outputImage->deltaY;
				zr[j1 - *jMin] = 0.0;
			}
			OCTTransform(myCt, n, xr, yr, zr);
			for (j1 = *jMin; j1 < *jMax; j1++)
			{
				/* fractional index in reduced buffer: block b is centred on source pixel c0 + b*kx + (kx-1)/2 */
				col = ((xr[j1 - *jMin] - g.x0) / g.dx - c0 - 0.5 * (kx - 1)) / kx;
				row = ((yr[j1 - *jMin] - g.y0) / g.dy - r0 - 0.5 * (ky - 1)) / ky;
				gam = bilinearInterp(gamma, col, row, nbc, nbr, 0.0, GCOVINVALID);
				value = -LARGEINT;
				if (gam > 0)
				{
					sig = (needSigma == TRUE) ? bilinearInterp(sigma, col, row, nbc, nbr, 0.0, GCOVINVALID) : gam;
					if (sig > 0)
					{
						if ((S1Cal & TRUE) == TRUE)
						{
							/* sigma0 is mosaicked; gamma0 = sigma0 + gBuf (dB) at output. Same clamp as makeGeoMosaic */
							value = sig;
							gBufTmp[i1][j1] = round(10.0 * log10(gam / sig) * 100.) / 100.;
							gBufTmp[i1][j1] = min(max(gBufTmp[i1][j1], -29.9), 35.0);
						}
						else
						{
							value = (calOutput == CALOUTPUT_SIGMA0) ? sig : gam;
						}
						if (needCube == TRUE)
						{
							/* DEM is in the output projection (km); lat and Re are unused for ELLIPSOIDAL */
							h = getXYHeightXY((outputImage->originX + j1 * outputImage->deltaX) * MTOKM,
											  (outputImage->originY + i1 * outputImage->deltaY) * MTOKM, 0.0,
											  (xyDEM *)dem, 0.0, ELLIPSOIDAL);
							psi = interpIncidence(cube, hgt, nh, cnx, cny, cx0, cy0, cdx, cdy,
												  xr[j1 - *jMin], yr[j1 - *jMin], h);
							if (doPsi == TRUE)
							{
								psiBufTmp[i1][j1] = psi;
							}
						}
					}
				}
				if (rangeSelect != RANGESELECT_NONE)
				{
					selTmp[i1][j1] = (unsigned char)incWeight(incBuf, i1, j1, psi);
				}
				imageTmp[i1][j1] = value;
				if (value > 0)
				{
					scaleTmp[i1][j1] = 1;
				}
			}
		}
		free(xr);
		free(yr);
		free(zr);
	}
	for (t = 0; t < nThreads; t++)
	{
		OCTDestroyCoordinateTransformation(ctThread[t]);
	}
	free(ctThread);
	fprintf(stderr, "GCOV read+average %.1f s, resample %.1f s\n", tRead, omp_get_wtime() - tSample);
	freeFloat2D(gamma);
	if (needSigma == TRUE)
	{
		freeFloat2D(sigma);
	}
	if (needCube == TRUE)
	{
		for (t = 0; t < nh; t++)
		{
			freeFloat2D(cube[t]);
		}
		free(cube);
		free(hgt);
	}
	return TRUE;
}

/*
  Block average kx x ky source pixels into an nbr x nbc buffer. Samples that are NaN, <= 0, or
  flagged by the mask (0 = invalid, 255 = fill) are excluded. If sigma != NULL, gamma*factor is
  averaged over the same samples (those with a valid positive factor). Empty blocks are GCOVINVALID.
*/
static float **reduceGCOV(gcovFile *g, int32_t c0, int32_t r0, int32_t nbc, int32_t nbr, int32_t kx, int32_t ky,
						  int32_t useMask, float ***sigma, int32_t maskOnly)
{
	float **gamma, *gBuf, *fBuf;
	unsigned char *mBuf;
	unsigned short *tBuf = NULL; /* raw float16 staging, only when a band is float16 */
	double *gSum, *sSum;
	int32_t *count, width, br, br0, bc, r, c, k, nbStrip, nbThis;
	float gv, fv;
	GDALRasterBandH hG, hF, hM;

	width = nbc * kx;
	/*
	  Read strips of ~1024 rows. The GCOV layers are gzip compressed in 512x512 chunks, so reading
	  only ky rows at a time decompresses every chunk row ~512/ky times (measured ~20x slower).
	*/
	nbStrip = max(1, 1024 / ky);
	gamma = allocFloat2D(nbr, nbc);
	if (sigma != NULL)
	{
		*sigma = allocFloat2D(nbr, nbc);
	}
	gBuf = (float *)malloc(sizeof(float) * (size_t)width * ky * nbStrip);
	fBuf = (float *)malloc(sizeof(float) * (size_t)width * ky * nbStrip);
	mBuf = (unsigned char *)malloc(sizeof(unsigned char) * (size_t)width * ky * nbStrip);
	gSum = (double *)malloc(sizeof(double) * nbc);
	sSum = (double *)malloc(sizeof(double) * nbc);
	count = (int32_t *)malloc(sizeof(int32_t) * nbc);
	if (g->h5Gamma != H5I_INVALID_HID || g->h5Factor != H5I_INVALID_HID)
	{
		tBuf = (unsigned short *)malloc(sizeof(unsigned short) * (size_t)width * ky * nbStrip);
		if (tBuf == NULL)
		{
			error("reduceGCOV: malloc failed");
		}
	}
	if (gamma == NULL || gBuf == NULL || fBuf == NULL || mBuf == NULL)
	{
		error("reduceGCOV: malloc failed");
	}
	/*
	  maskOnly (pass 1 of the range selection): only the footprint is needed, so read just the
	  mask - 1.3 MB compressed against 1.5 GB for HHHH - and emit 1 where a block has any valid
	  sample, using the same per-sample test as pass 2. mask == 1 is a strict subset of "HHHH is
	  valid" (measured: zero false positives), and under-claiming is the safe direction - it only
	  makes the filter more permissive. With no mask there is nothing cheap to test, so fall back
	  to reading gamma.
	*/
	if (maskOnly == TRUE && g->hMask == NULL)
	{
		maskOnly = FALSE;
	}
	hG = GDALGetRasterBand(g->hGamma, 1);
	hF = (sigma != NULL) ? GDALGetRasterBand(g->hFactor, 1) : NULL;
	hM = (g->hMask != NULL && (useMask == TRUE || maskOnly == TRUE)) ? GDALGetRasterBand(g->hMask, 1) : NULL;
	for (br0 = 0; br0 < nbr; br0 += nbStrip)
	{
		nbThis = min(nbStrip, nbr - br0);
		if (maskOnly == FALSE &&
			readGCOVWindow(g->h5Gamma, hG, c0, r0 + br0 * ky, width, ky * nbThis, gBuf, tBuf) == FALSE)
		{
			error("reduceGCOV: read error");
		}
		if (hF != NULL && readGCOVWindow(g->h5Factor, hF, c0, r0 + br0 * ky, width, ky * nbThis,
										 fBuf, tBuf) == FALSE)
		{
			error("reduceGCOV: factor read error");
		}
		if (hM != NULL && GDALRasterIO(hM, GF_Read, c0, r0 + br0 * ky, width, ky * nbThis, mBuf, width, ky * nbThis,
									   GDT_Byte, 0, 0) != CE_None)
		{
			error("reduceGCOV: mask read error");
		}
		for (br = br0; br < br0 + nbThis; br++)
		{
			for (bc = 0; bc < nbc; bc++)
			{
				gSum[bc] = 0.0;
				sSum[bc] = 0.0;
				count[bc] = 0;
			}
			for (r = (br - br0) * ky; r < (br - br0 + 1) * ky; r++)
			{
				for (c = 0; c < width; c++)
				{
					k = r * width + c;
					if (maskOnly == TRUE)
					{
						/* exactly the per-sample test pass 2 applies below, so the two passes agree on
						   the footprint. Over-claiming is the harmful direction: it pulls the extremum to
						   a frame with no data there, pushing real data onto the fallback. */
						if (mBuf[k] == 0 || mBuf[k] == 255)
						{
							continue;
						}
						gSum[c / kx] += 1.0;
						count[c / kx]++;
						continue;
					}
					gv = gBuf[k];
					if (!(gv > 0) || isinf(gv))
					{
						continue;
					}
					if (hM != NULL && (mBuf[k] == 0 || mBuf[k] == 255))
					{
						continue;
					}
					if (hF != NULL)
					{
						fv = fBuf[k];
						if (!(fv > 0) || isinf(fv))
						{
							continue;
						}
						sSum[c / kx] += gv * fv;
					}
					gSum[c / kx] += gv;
					count[c / kx]++;
				}
			}
			for (bc = 0; bc < nbc; bc++)
			{
				gamma[br][bc] = (count[bc] > 0) ? gSum[bc] / count[bc] : GCOVINVALID;
				if (sigma != NULL)
				{
					(*sigma)[br][bc] = (count[bc] > 0) ? sSum[bc] / count[bc] : GCOVINVALID;
				}
			}
		}
	}
	free(gBuf);
	free(fBuf);
	free(mBuf);
	free(tBuf);
	free(gSum);
	free(sSum);
	free(count);
	return gamma;
}

/*
  Read the metadata incidenceAngle cube (heights x rows x cols) and its coordinates.
*/
static int32_t readIncidenceCube(char *file, gcovInputs *gcov, gcovFile *g, float ****cube, double **hgt,
								 int32_t *nh, int32_t *cnx, int32_t *cny, double *cx0, double *cy0, double *cdx, double *cdy)
{
	char path[2048];
	GDALDatasetH hCube;
	double *xc, *yc;
	int32_t k;
	if (getEPSGAttribute(g->hFile, "science_LSAR_GCOV_metadata_radarGrid_projection_epsg_code") != g->epsg)
	{
		return FALSE;
	}
	sprintf(path, "HDF5:\"%s\"://science/LSAR/GCOV/metadata/radarGrid/incidenceAngle", file);
	hCube = GDALOpen(path, GA_ReadOnly);
	if (hCube == NULL)
	{
		return FALSE;
	}
	*cnx = GDALGetRasterXSize(hCube);
	*cny = GDALGetRasterYSize(hCube);
	*nh = GDALGetRasterCount(hCube);
	if (readCoordinates(file, "metadata/radarGrid/heightAboveEllipsoid", *nh, hgt) == FALSE ||
		readCoordinates(file, "metadata/radarGrid/xCoordinates", *cnx, &xc) == FALSE ||
		readCoordinates(file, "metadata/radarGrid/yCoordinates", *cny, &yc) == FALSE)
	{
		GDALClose(hCube);
		return FALSE;
	}
	*cx0 = xc[0];
	*cy0 = yc[0];
	*cdx = (xc[*cnx - 1] - xc[0]) / (*cnx - 1);
	*cdy = (yc[*cny - 1] - yc[0]) / (*cny - 1);
	free(xc);
	free(yc);
	*cube = (float ***)malloc(sizeof(float **) * (*nh));
	for (k = 0; k < *nh; k++)
	{
		(*cube)[k] = allocFloat2D(*cny, *cnx);
		if (GDALRasterIO(GDALGetRasterBand(hCube, k + 1), GF_Read, 0, 0, *cnx, *cny, (*cube)[k][0],
						 *cnx, *cny, GDT_Float32, 0, 0) != CE_None)
		{
			error("readIncidenceCube: read error %s", path);
		}
	}
	GDALClose(hCube);
	return TRUE;
}

/*
  Trilinear interpolation of the incidence angle cube at map coords X, Y and height h.
*/
static float interpIncidence(float ***cube, double *hgt, int32_t nh, int32_t cnx, int32_t cny,
							 double cx0, double cy0, double cdx, double cdy, double X, double Y, double h)
{
	double col, row, t, u, w, v[2];
	int32_t i, j, k, m;
	col = (X - cx0) / cdx;
	row = (Y - cy0) / cdy;
	j = min(max((int32_t)col, 0), cnx - 2);
	i = min(max((int32_t)row, 0), cny - 2);
	t = min(max(col - j, 0.0), 1.0);
	u = min(max(row - i, 0.0), 1.0);
	for (k = 0; k < nh - 2 && h > hgt[k + 1]; k++)
	{
	}
	w = min(max((h - hgt[k]) / (hgt[k + 1] - hgt[k]), 0.0), 1.0);
	for (m = 0; m < 2; m++)
	{
		v[m] = (1 - t) * (1 - u) * cube[k + m][i][j] + t * (1 - u) * cube[k + m][i][j + 1] +
			   t * u * cube[k + m][i + 1][j + 1] + (1 - t) * u * cube[k + m][i + 1][j];
	}
	return (float)((1 - w) * v[0] + w * v[1]);
}

/*
  Open the GCOV datasets and reconstruct the grid geometry.
*/
static int32_t openGCOV(gcovInputs *gcov, char *file, gcovFile *g, int32_t needFactor)
{
	char path[2048], key[256], grid[256], h5path[2048];
	double *xc, *yc;
	memset(g, 0, sizeof(gcovFile));
	/* memset gives 0, which is a valid hid_t; the "no handle" sentinel is H5I_INVALID_HID */
	g->h5Gamma = g->h5GammaFile = g->h5Factor = g->h5FactorFile = H5I_INVALID_HID;
	g->hFile = GDALOpen(file, GA_ReadOnly);
	if (g->hFile == NULL)
	{
		/* A local granule that will not open is one bad product among many, and skipping it
		   so the rest of the tile still builds is deliberate.  A /vsi path is different: the
		   usual cause is a transient network or an expired presigned URL, there is no
		   redundancy to fall back on, and skipping every input yields an all-nodata mosaic
		   that exits 0 and looks finished.  Fail loudly instead. */
		if (strncmp(file, "/vsi", 4) == 0)
		{
			error("openGCOV: could not open remote input %s\n"
				  "  (transient network error, or an expired presigned URL).", file);
		}
		fprintf(stderr, "Warning: could not open %s\n", file);
		return FALSE;
	}
	sprintf(grid, "science/LSAR/GCOV/grids/frequency%s", gcov->frequency);
	sprintf(path, "HDF5:\"%s\"://%s/%s", file, grid, gcov->polarization);
	g->hGamma = GDALOpen(path, GA_ReadOnly);
	if (g->hGamma == NULL)
	{
		fprintf(stderr, "Warning: no %s in %s\n", path, file);
		closeGCOV(g);
		return FALSE;
	}
	g->nx = GDALGetRasterXSize(g->hGamma);
	g->ny = GDALGetRasterYSize(g->hGamma);
	sprintf(h5path, "/science/LSAR/GCOV/grids/frequency%s/%s", gcov->frequency, gcov->polarization);
	openIfFloat16(file, h5path, &(g->h5GammaFile), &(g->h5Gamma));
	if (needFactor == TRUE)
	{
		char factorFile[4096], *base;
		/* With factorFrom set, the factor lives outside the granule (slim products share one
		   per track/frame/grid across cycles). The downloader leaves a symlink named exactly
		   like the granule, so the basename is all we need. */
		strcpy(factorFile, file);
		if (gcov->factorDir[0] != '\0')
		{
			char shared[4096];
			base = strrchr(file, '/');
			base = (base == NULL) ? file : base + 1;
			if (snprintf(shared, sizeof(shared), "%s/%s", gcov->factorDir, base) >= (int32_t)sizeof(shared))
			{
				error("openGCOV: factor path too long for %s", base);
			}
			/* Fall back to the granule when there is no external factor for it, so one
			   directory can hold slim products (factor shared, symlink present) and
			   archive products (factor in the file) at the same time. Still loud if
			   neither has one: the GDALOpen below fails. */
			if (access(shared, R_OK) == 0)
			{
				strcpy(factorFile, shared);
			}
		}
		sprintf(path, "HDF5:\"%s\"://%s/rtcGammaToSigmaFactor", factorFile, grid);
		if ((g->hFactor = GDALOpen(path, GA_ReadOnly)) == NULL)
		{
			error("openGCOV: could not open %s", path);
		}
		/* A factor on a different grid than the gamma band reads at the same pixel window and
		   silently misregisters the gamma->sigma conversion. Partial frames make this real:
		   a track/frame that arrives partial in one cycle and full in the next has a
		   different grid, so a shared factor must be rejected rather than reused. */
		if (GDALGetRasterXSize(g->hFactor) != g->nx || GDALGetRasterYSize(g->hFactor) != g->ny)
		{
			error("openGCOV: factor grid %i x %i does not match gamma grid %i x %i in %s",
				  GDALGetRasterXSize(g->hFactor), GDALGetRasterYSize(g->hFactor), g->nx, g->ny, factorFile);
		}
		sprintf(h5path, "/science/LSAR/GCOV/grids/frequency%s/rtcGammaToSigmaFactor", gcov->frequency);
		openIfFloat16(factorFile, h5path, &(g->h5FactorFile), &(g->h5Factor));
	}
	if (gcov->useMask == TRUE)
	{
		sprintf(path, "HDF5:\"%s\"://%s/mask", file, grid);
		if ((g->hMask = GDALOpen(path, GA_ReadOnly)) == NULL)
		{
			error("openGCOV: could not open %s", path);
		}
	}
	sprintf(key, "science_LSAR_GCOV_grids_frequency%s_projection_epsg_code", gcov->frequency);
	g->epsg = getEPSGAttribute(g->hFile, key);
	if (g->epsg <= 0)
	{
		error("openGCOV: no %s in %s", key, file);
	}
	sprintf(path, "grids/frequency%s/xCoordinates", gcov->frequency);
	sprintf(key, "grids/frequency%s/yCoordinates", gcov->frequency);
	if (readCoordinates(file, path, g->nx, &xc) == FALSE || readCoordinates(file, key, g->ny, &yc) == FALSE)
	{
		error("openGCOV: could not read coordinates from %s", file);
	}
	/* Coordinates are pixel centres */
	g->x0 = xc[0];
	g->y0 = yc[0];
	g->dx = (xc[g->nx - 1] - xc[0]) / (g->nx - 1);
	g->dy = (yc[g->ny - 1] - yc[0]) / (g->ny - 1);
	free(xc);
	free(yc);
	fprintf(stderr, "GCOV grid %i x %i EPSG:%i origin %f %f spacing %f %f\n",
			g->nx, g->ny, g->epsg, g->x0, g->y0, g->dx, g->dy);
	return TRUE;
}

static void closeGCOV(gcovFile *g)
{
	if (g->hGamma != NULL)
	{
		GDALClose(g->hGamma);
	}
	if (g->hFactor != NULL)
	{
		GDALClose(g->hFactor);
	}
	if (g->hMask != NULL)
	{
		GDALClose(g->hMask);
	}
	if (g->hFile != NULL)
	{
		GDALClose(g->hFile);
	}
	g->hGamma = g->hFactor = g->hMask = g->hFile = NULL;
	/* the float16 handles, when this product had any */
	if (g->h5Gamma != H5I_INVALID_HID)
	{
		H5Dclose(g->h5Gamma);
	}
	if (g->h5GammaFile != H5I_INVALID_HID)
	{
		H5Fclose(g->h5GammaFile);
	}
	if (g->h5Factor != H5I_INVALID_HID)
	{
		H5Dclose(g->h5Factor);
	}
	if (g->h5FactorFile != H5I_INVALID_HID)
	{
		H5Fclose(g->h5FactorFile);
	}
	g->h5Gamma = g->h5GammaFile = g->h5Factor = g->h5FactorFile = H5I_INVALID_HID;
}

/*
  Read a 1-D coordinate dataset (exposed by GDAL as an n x 1 raster) under science/LSAR/GCOV/.
*/
static int32_t readCoordinates(char *file, char *dataset, int32_t n, double **coords)
{
	char path[2048];
	GDALDatasetH hDS;
	sprintf(path, "HDF5:\"%s\"://science/LSAR/GCOV/%s", file, dataset);
	hDS = GDALOpen(path, GA_ReadOnly);
	if (hDS == NULL)
	{
		return FALSE;
	}
	if (GDALGetRasterXSize(hDS) * GDALGetRasterYSize(hDS) != n || n < 2)
	{
		fprintf(stderr, "readCoordinates: size mismatch for %s\n", path);
		GDALClose(hDS);
		return FALSE;
	}
	*coords = (double *)malloc(sizeof(double) * n);
	if (GDALRasterIO(GDALGetRasterBand(hDS, 1), GF_Read, 0, 0, GDALGetRasterXSize(hDS), GDALGetRasterYSize(hDS),
					 *coords, GDALGetRasterXSize(hDS), GDALGetRasterYSize(hDS), GDT_Float64, 0, 0) != CE_None)
	{
		GDALClose(hDS);
		return FALSE;
	}
	GDALClose(hDS);
	return TRUE;
}

static int32_t getEPSGAttribute(GDALDatasetH hFile, char *key)
{
	const char *value;
	value = GDALGetMetadataItem(hFile, key, NULL);
	if (value == NULL)
	{
		return -1;
	}
	return atoi(value);
}

/*
  NISAR_L2_PR_GCOV_<cycle>_<track>_<A|D>_<frame>_<mode>_<pol>_<A>_<YYYYMMDDTHHMMSS>_...
*/
static int32_t parseGCOVName(char *file, int32_t *year, int32_t *month, int32_t *day, int32_t *passType)
{
	char name[1024], *base, *tok, *fields[20];
	int32_t n;
	base = strrchr(file, '/');
	strncpy(name, (base == NULL) ? file : base + 1, sizeof(name) - 1);
	name[sizeof(name) - 1] = '\0';
	n = 0;
	for (tok = strtok(name, "_"); tok != NULL && n < 20; tok = strtok(NULL, "_"))
	{
		fields[n++] = tok;
	}
	if (n < 12 || strcmp(fields[3], "GCOV") != 0)
	{
		return FALSE;
	}
	if (sscanf(fields[11], "%4d%2d%2d", year, month, day) != 3)
	{
		return FALSE;
	}
	if (fields[6][0] == 'A')
	{
		*passType = ASCENDING;
	}
	else if (fields[6][0] == 'D')
	{
		*passType = DESCENDING;
	}
	else
	{
		return FALSE;
	}
	return TRUE;
}

static int32_t outputEPSG()
{
	extern double Rotation;
	extern double SLat;
	extern int32_t HemiSphere;
	const grimpProj *p = grimpDefaultProj();
	/* Prefer the resolved output projection's own code. Deriving it from Rotation/SLat only
	   ever produces a polar stereographic EPSG, so a UTM or geographic output grid would be
	   transformed as if it were polar stereographic -- silently, with plausible-looking
	   coordinates. A geographic grid has no rot/stdLat to derive from at all. */
	if (p != NULL && p->epsg != 0)
	{
		return p->epsg;
	}
	return atoi(getEPSGFromProjectionParams(Rotation, SLat, HemiSphere));
}

static OGRCoordinateTransformationH makeTransform(int32_t fromEPSG, int32_t toEPSG)
{
	OGRSpatialReferenceH from, to;
	OGRCoordinateTransformationH ct;
	from = OSRNewSpatialReference(NULL);
	to = OSRNewSpatialReference(NULL);
	if (OSRImportFromEPSG(from, fromEPSG) != OGRERR_NONE || OSRImportFromEPSG(to, toEPSG) != OGRERR_NONE)
	{
		error("makeTransform: invalid EPSG %i or %i", fromEPSG, toEPSG);
	}
	/* x=lon, y=lat for geographic CRSs */
	OSRSetAxisMappingStrategy(from, OAMS_TRADITIONAL_GIS_ORDER);
	OSRSetAxisMappingStrategy(to, OAMS_TRADITIONAL_GIS_ORDER);
	ct = OCTNewCoordinateTransformation(from, to);
	if (ct == NULL)
	{
		error("makeTransform: could not create transform %i -> %i", fromEPSG, toEPSG);
	}
	OSRDestroySpatialReference(from);
	OSRDestroySpatialReference(to);
	return ct;
}

/*
  Transform the edges of the box [xa,xb] x [ya,yb] and return the bounding box of the result.
*/
static int32_t transformBox(OGRCoordinateTransformationH ct, double xa, double xb, double ya, double yb,
							double *minX, double *maxX, double *minY, double *maxY)
{
	double x[4 * NEDGE], y[4 * NEDGE], z[4 * NEDGE], f;
	int32_t ok[4 * NEDGE], i, n, nGood;
	n = 0;
	for (i = 0; i < NEDGE; i++)
	{
		f = (double)i / (NEDGE - 1);
		x[n] = xa + f * (xb - xa);
		y[n++] = ya;
		x[n] = xa + f * (xb - xa);
		y[n++] = yb;
		x[n] = xa;
		y[n++] = ya + f * (yb - ya);
		x[n] = xb;
		y[n++] = ya + f * (yb - ya);
	}
	for (i = 0; i < n; i++)
	{
		z[i] = 0.0;
	}
	OCTTransformEx(ct, n, x, y, z, ok);
	*minX = *minY = 1.e30;
	*maxX = *maxY = -1.e30;
	nGood = 0;
	for (i = 0; i < n; i++)
	{
		if (ok[i])
		{
			*minX = min(*minX, x[i]);
			*maxX = max(*maxX, x[i]);
			*minY = min(*minY, y[i]);
			*maxY = max(*maxY, y[i]);
			nGood++;
		}
	}
	return (nGood > 0) ? TRUE : FALSE;
}

/*
  Smallest non-negative shift that moves start onto the grid anchor + n*k (anchor rounded to
  the nearest source pixel).
*/
/*
  Output-lattice pixel centre nearest the GCOV origin, returned in GCOV coordinates.
*/
static void gcovAnchor(gcovFile *g, outputImageStructure *outputImage, double *ax, double *ay)
{
	OGRCoordinateTransformationH ct;
	double x, y, z = 0.0;
	x = g->x0;
	y = g->y0;
	ct = makeTransform(g->epsg, outputEPSG());
	OCTTransform(ct, 1, &x, &y, &z);
	OCTDestroyCoordinateTransformation(ct);
	x = outputImage->originX + round((x - outputImage->originX) / outputImage->deltaX) * outputImage->deltaX;
	y = outputImage->originY + round((y - outputImage->originY) / outputImage->deltaY) * outputImage->deltaY;
	ct = makeTransform(outputEPSG(), g->epsg);
	z = 0.0;
	OCTTransform(ct, 1, &x, &y, &z);
	OCTDestroyCoordinateTransformation(ct);
	*ax = x;
	*ay = y;
}

static int32_t alignOffset(int32_t start, double anchor, int32_t k)
{
	int32_t phase = (int32_t)lround(anchor - floor(anchor / k) * k) % k;
	return ((phase - start) % k + k) % k;
}

static float **allocFloat2D(int32_t nr, int32_t nc)
{
	float **a, *buf;
	int32_t i;
	buf = (float *)malloc(sizeof(float) * (size_t)nr * nc);
	a = (float **)malloc(sizeof(float *) * nr);
	if (buf == NULL || a == NULL)
	{
		error("allocFloat2D: malloc failed %i x %i", nr, nc);
	}
	for (i = 0; i < nr; i++)
	{
		a[i] = &(buf[(size_t)i * nc]);
	}
	return a;
}

static void freeFloat2D(float **a)
{
	free(a[0]);
	free(a);
}
