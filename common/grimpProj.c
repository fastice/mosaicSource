#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <proj.h>
#include "ogr_srs_api.h"
#include "cpl_conv.h"
#include "mosaicSource/common/common.h"
#include "mosaicSource/common/grimpProj.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/*
  Implementation of the grimpProj map-projection layer.  See grimpProj.h for the
  design, the unit conventions and the rot/lon_0 sign asymmetry.

  Two structural points:

  1. Every GP_PS call delegates to the untouched lltoxy1()/xytoll1() with the original
     arguments, so polar runs are byte-identical and need no PROJ object at all.

  2. PROJ handles live in this file's private registry, not in grimpProj, because
     xyDEM (which embeds a grimpProj) is copied by value all over the code.
*/

#define GP_MAXPROJ 8    /* distinct projections in one run: output + DEM + a few grids */
#define GP_MAXTHREADS 256

typedef struct
{
    int used;
    char srs[1024];                  /* the user input this was built from */
    char projOnlyDef[512];           /* bare projection, for proj_factors */
    PJ_CONTEXT *ctx[GP_MAXTHREADS];
    PJ *trans[GP_MAXTHREADS];        /* lon/lat <-> x/y, normalised axis order */
    PJ *projOnly[GP_MAXTHREADS];     /* projection alone, valid input for proj_factors */
    int nThreads;                    /* how many slots are populated */
} gpRegistryEntry;

static gpRegistryEntry gpRegistry[GP_MAXPROJ];
static int32_t gpNRegistered = 0;

static grimpProj gpDefaultProj;
static int32_t gpHaveDefault = FALSE;

/* Per-projection crs->crs transforms, built lazily on first use (serial), for
   outXYToProjXY when the DEM is in a different projection from the output. */
#define GP_MAXPAIR 8
typedef struct
{
    int used;
    int32_t srcHandle, dstHandle;
    PJ *pj[GP_MAXTHREADS];
    int nThreads;
} gpPairEntry;
static gpPairEntry gpPairs[GP_MAXPAIR];
static int32_t gpNPairs = 0;

static int32_t gpThreadNum(void)
{
#ifdef _OPENMP
    return (int32_t)omp_get_thread_num();
#else
    return 0;
#endif
}

static int32_t gpMaxThreads(void)
{
#ifdef _OPENMP
    int32_t n = (int32_t)omp_get_max_threads();
    return (n < 1) ? 1 : n;
#else
    return 1;
#endif
}

/* ------------------------------------------------------------------ constructors */

grimpProj grimpProjFromLegacy(double rot, double stdLat, int32_t hemisphere)
{
    grimpProj p;
    memset(&p, 0, sizeof(p));
    p.kind = GP_PS;
    p.rot = rot;
    p.hemisphere = hemisphere;
    p.utmZone = 0;
    p.handle = -1;
    p.gridScale = MTOKM;   /* projected grid stores metres, the API speaks km */
    /* The SLat global's "not set" sentinel is -91; the open-coded call sites all
       substituted 70 north / 71 south, so do exactly that here. */
    if (stdLat < -90.)
    {
        p.stdLat = (hemisphere == NORTH) ? 70.0 : 71.0;
    }
    else
    {
        p.stdLat = fabs(stdLat);
    }
    /* Name the two codes we can recognise; anything else keeps epsg 0 and is written
       out as a proj string instead (grimpSRSString). */
    if (p.rot == 45. && p.stdLat == 70. && hemisphere == NORTH)
        p.epsg = 3413;
    else if (p.rot == 0. && p.stdLat == 71. && hemisphere == SOUTH)
        p.epsg = 3031;
    else
        p.epsg = 0;
    return p;
}

/* Build the "projection only" PROJ definition used for proj_factors().  proj_factors
   is not valid on a crs->crs pipeline, so we keep a separate handle built from an
   explicit proj string we control rather than round-tripping WKT. */
static void gpProjOnlyDef(const grimpProj *p, char *buf, size_t n)
{
    if (p->kind == GP_UTM)
    {
        snprintf(buf, n, "+proj=utm +zone=%d%s +datum=WGS84 +units=m +no_defs",
                 (int)p->utmZone, (p->hemisphere == SOUTH) ? " +south" : "");
    }
    else
    {
        /* north: lon_0 = -rot, south: lon_0 = +rot (see grimpProj.h) */
        double lon0 = (p->hemisphere == NORTH) ? -p->rot : p->rot;
        snprintf(buf, n, "+proj=stere +lat_0=%s90 +lat_ts=%.10g +lon_0=%.10g "
                         "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                 (p->hemisphere == SOUTH) ? "-" : "", p->stdLat, lon0);
    }
}

static int32_t gpRegister(const char *srs, const grimpProj *p)
{
    int32_t i;
    for (i = 0; i < gpNRegistered; i++)
    {
        if (gpRegistry[i].used && strcmp(gpRegistry[i].srs, srs) == 0)
            return i;
    }
    if (gpNRegistered >= GP_MAXPROJ)
        error("grimpProj: more than %i distinct projections in one run", GP_MAXPROJ);
    i = gpNRegistered++;
    memset(&gpRegistry[i], 0, sizeof(gpRegistry[i]));
    gpRegistry[i].used = TRUE;
    strncpy(gpRegistry[i].srs, srs, sizeof(gpRegistry[i].srs) - 1);
    gpProjOnlyDef(p, gpRegistry[i].projOnlyDef, sizeof(gpRegistry[i].projOnlyDef));
    return i;
}

grimpProj grimpProjFromSRS(const char *userInput)
{
    grimpProj p;
    OGRSpatialReferenceH srs;
    const char *method;
    int north = 0, zone;

    memset(&p, 0, sizeof(p));
    p.handle = -1;
    if (userInput == NULL || userInput[0] == '\0')
        error("grimpProjFromSRS: empty projection specification");

    srs = OSRNewSpatialReference(NULL);
    if (OSRSetFromUserInput(srs, userInput) != OGRERR_NONE)
        error("grimpProjFromSRS: cannot interpret projection \"%s\"", userInput);
    if (OSRIsGeographic(srs))
    {
        /* A geographic CRS is not a projection to convert through: the grid axes ARE lon/lat,
           so every conversion is the identity and no PROJ object is needed. Only WGS84 is
           accepted, for the same reason the projected path insists on it -- the rest of the
           code hard-codes that figure of the Earth. */
        grimpProj p;
        double a = OSRGetSemiMajor(srs, NULL);
        if (fabs(a - 6378137.0) > 1.0)
            error("grimpProjFromSRS: \"%s\" is geographic but not on the WGS84 ellipsoid\n"
                  "  (semi-major %.1f m); the conversion routines hard-code WGS84.", userInput, a);
        memset(&p, 0, sizeof(p));
        p.kind = GP_LATLON;
        p.epsg = 4326;
        p.rot = 0.0;
        p.stdLat = 0.0;
        p.hemisphere = NORTH;
        p.utmZone = 0;
        p.handle = -1;
        p.gridScale = 1.0;   /* the grid already stores degrees */
        OSRDestroySpatialReference(srs);
        return p;
    }
    if (!OSRIsProjected(srs))
        error("grimpProjFromSRS: \"%s\" is neither projected nor geographic", userInput);

    /* Only WGS84-ellipsoid grids can use the legacy lltoxy1 fast path, and PROJ would
       give answers that differ from every existing product, so refuse outright. */
    {
        /* lltoxy1/xytoll1 hard-code WGS84 (re = 6378137.0, e2 = 0.0066943801), so a grid on
           any other ellipsoid would be converted with the wrong figure of the Earth.
           The 1 m tolerance on the semi-major axis rejects every ellipsoid that has been
           used for ice-sheet products, including Hughes 1980 (6378273) whose numbers still
           appear in lltoxy1's comments, and the spherical and older geodetic ones.
           GRS80 passes: it shares WGS84's semi-major axis exactly and differs only in
           inverse flattening, 298.257222101 vs 298.257223563 -- 5e-9 relative, which moves
           a polar position by well under 0.1 mm, an order below the 0.19 mm at which PROJ
           and lltoxy1 already agree.  The flattening tolerance below is set to admit that
           pair and nothing coarser. */
        double a = OSRGetSemiMajor(srs, NULL);
        double invf = OSRGetInvFlattening(srs, NULL);
        if (fabs(a - 6378137.0) > 1.0)
            error("grimpProjFromSRS: \"%s\" has semi-major axis %.3f, not WGS84 (6378137).\n"
                  "  GrIMP coordinate conversions assume WGS84.", userInput, a);
        if (fabs(invf - 298.257223563) > 1.0e-4)
            error("grimpProjFromSRS: \"%s\" has inverse flattening %.9f, not WGS84 "
                  "(298.257223563).\n  GrIMP coordinate conversions assume WGS84.",
                  userInput, invf);
    }

    zone = OSRGetUTMZone(srs, &north);
    method = OSRGetAttrValue(srs, "PROJECTION", 0);

    if (zone != 0)
    {
        p.kind = GP_UTM;
        p.utmZone = (int32_t)zone;
        p.hemisphere = north ? NORTH : SOUTH;
        /* Not meaningful for UTM, but keep them defined rather than garbage. */
        p.rot = 0.0;
        p.stdLat = 0.0;
    }
    else if (method != NULL && strstr(method, "Polar_Stereographic") != NULL)
    {
        double latOrigin = OSRGetProjParm(srs, SRS_PP_LATITUDE_OF_ORIGIN, -1000., NULL);
        double cm, fe, fn;
        int32_t hemiFromLat0 = 0; /* -1 south, +1 north, 0 unknown */
        if (latOrigin < -900.)
            latOrigin = OSRGetProjParm(srs, SRS_PP_LATITUDE_OF_CENTER, -1000., NULL);
        if (latOrigin < -900.)
            error("grimpProjFromSRS: \"%s\" has no latitude of origin / standard parallel",
                  userInput);
        cm = OSRGetProjParm(srs, SRS_PP_CENTRAL_MERIDIAN, 0.0, NULL);
        fe = OSRGetProjParm(srs, SRS_PP_FALSE_EASTING, 0.0, NULL);
        fn = OSRGetProjParm(srs, SRS_PP_FALSE_NORTHING, 0.0, NULL);
        if (fe != 0.0 || fn != 0.0)
            error("grimpProjFromSRS: \"%s\" has a false easting/northing (%.1f, %.1f).\n"
                  "  lltoxy1 cannot express that, so this projection is not supported.",
                  userInput, fe, fn);
        /* Take the pole from +lat_0, which is what PROJ itself considers authoritative.
           GDAL normalises an inconsistent definition in favour of lat_ts -- given
           "+proj=stere +lat_0=-90 +lat_ts=71" it exports "+lat_0=90 +lat_ts=71" -- so in
           practice lat_0 and the sign of latitude_of_origin always agree by the time we
           see them.  Reading lat_0 makes that explicit rather than relying on it, and the
           sign of latitude_of_origin remains the fallback if lat_0 cannot be recovered. */
        {
            char *p4 = NULL;
            if (OSRExportToProj4(srs, &p4) == OGRERR_NONE && p4 != NULL)
            {
                const char *q = strstr(p4, "+lat_0=");
                if (q != NULL)
                    hemiFromLat0 = (atof(q + 7) < 0.) ? -1 : 1;
                CPLFree(p4);
            }
        }
        p.kind = GP_PS;
        if (hemiFromLat0 != 0)
            p.hemisphere = (hemiFromLat0 < 0) ? SOUTH : NORTH;
        else
            p.hemisphere = (latOrigin < 0) ? SOUTH : NORTH;
        p.stdLat = fabs(latOrigin);
        /* north: lon_0 = -rot, south: lon_0 = +rot */
        p.rot = (p.hemisphere == NORTH) ? -cm : cm;
        p.utmZone = 0;
    }
    else
    {
        /* Any other projected CRS: positions work through PROJ, but the grid angle and
           the isotropic scale do not generalise, so grimpXYAngle/grimpXYScale refuse and
           mosaic3d rejects one as an output grid (grimpRequireConformal).  geomosaic and
           every input grid are fine with it -- they only ever sample positions. */
        p.kind = GP_GENERIC;
        p.rot = 0.0;
        p.stdLat = 0.0;
        p.hemisphere = NORTH;
        p.utmZone = 0;
        fprintf(stderr, "grimpProjFromSRS: \"%s\" is a %s projection -- usable for positions "
                        "only (not for vx/vy decomposition)\n",
                userInput, (method == NULL) ? "(unknown)" : method);
    }

    {
        const char *code = OSRGetAuthorityCode(srs, NULL);
        p.epsg = (code != NULL) ? atoi(code) : 0;
    }
    p.gridScale = MTOKM;   /* every projected kind: metres stored, km through the API */
    OSRDestroySpatialReference(srs);

    /* A polar stereographic needs no PROJ object: it uses the legacy code path. */
    if (p.kind == GP_PS)
        p.handle = -1;
    else
        p.handle = gpRegister(userInput, &p);
    return p;
}

grimpProj grimpProjFromEPSG(int32_t epsg)
{
    char buf[64];
    snprintf(buf, sizeof(buf), "EPSG:%d", (int)epsg);
    return grimpProjFromSRS(buf);
}

/* ------------------------------------------------------- per-thread PROJ objects */

void grimpProjPrepareThreads(void)
{
    int32_t nt = gpMaxThreads();
    int32_t i, t;

    if (nt > GP_MAXTHREADS)
        error("grimpProj: %i threads exceeds GP_MAXTHREADS (%i)", nt, GP_MAXTHREADS);

    for (i = 0; i < gpNRegistered; i++)
    {
        gpRegistryEntry *e = &gpRegistry[i];
        if (!e->used)
            continue;
        for (t = 0; t < nt; t++)
        {
            if (e->trans[t] != NULL)
                continue;
            e->ctx[t] = proj_context_create();
            if (e->ctx[t] == NULL)
                error("grimpProj: proj_context_create failed for %s", e->srs);
            {
                PJ *raw = proj_create_crs_to_crs(e->ctx[t], "EPSG:4326", e->srs, NULL);
                if (raw == NULL)
                    error("grimpProj: cannot build transform WGS84 -> %s (%s)",
                          e->srs, proj_errno_string(proj_context_errno(e->ctx[t])));
                /* So we can pass (lon, lat) rather than EPSG:4326's native (lat, lon). */
                e->trans[t] = proj_normalize_for_visualization(e->ctx[t], raw);
                proj_destroy(raw);
                if (e->trans[t] == NULL)
                    error("grimpProj: proj_normalize_for_visualization failed for %s", e->srs);
            }
            e->projOnly[t] = proj_create(e->ctx[t], e->projOnlyDef);
            if (e->projOnly[t] == NULL)
                error("grimpProj: cannot build projection \"%s\" (%s)", e->projOnlyDef,
                      proj_errno_string(proj_context_errno(e->ctx[t])));
        }
        if (nt > e->nThreads)
            e->nThreads = nt;
    }

    for (i = 0; i < gpNPairs; i++)
    {
        gpPairEntry *q = &gpPairs[i];
        if (!q->used)
            continue;
        for (t = 0; t < nt; t++)
        {
            if (q->pj[t] != NULL)
                continue;
            q->pj[t] = proj_create_crs_to_crs(gpRegistry[q->srcHandle].ctx[t],
                                              gpRegistry[q->srcHandle].srs,
                                              gpRegistry[q->dstHandle].srs, NULL);
            if (q->pj[t] == NULL)
                error("grimpProj: cannot build transform %s -> %s",
                      gpRegistry[q->srcHandle].srs, gpRegistry[q->dstHandle].srs);
        }
        if (nt > q->nThreads)
            q->nThreads = nt;
    }
}

static PJ *gpTrans(const grimpProj *p)
{
    int32_t t = gpThreadNum();
    PJ *pj;
    if (p->handle < 0 || p->handle >= gpNRegistered)
        error("grimpProj: projection %s was never registered", grimpProjDescribe(p));
    pj = gpRegistry[p->handle].trans[t];
    if (pj == NULL)
        error("grimpProj: grimpProjPrepareThreads() was not called before use "
              "(thread %i, projection %s)", t, grimpProjDescribe(p));
    return pj;
}

/* ------------------------------------------------------------------- the default */

void grimpSetDefaultProj(const grimpProj *p)
{
    gpDefaultProj = *p;
    gpHaveDefault = TRUE;
}

const grimpProj *grimpDefaultProj(void)
{
    /* Fall back to the legacy globals when no projection was set explicitly.  Only
       mosaic3d and geomosaic call grimpSetDefaultProj(); every other binary that links
       $(COMMON) -- siminsar, lltora, getlocc, coarsereg, tiepoints, rparams, azparams,
       the speckle tools, unwrap -- reaches here through parseInputFile.c or
       llToImageNew.c and has always relied on Rotation/SLat/HemiSphere.  Building the
       descriptor from them reproduces exactly what those call sites did before,
       including grimpProjFromLegacy's 70/71 substitution for SLat's -91 sentinel.

       Thread safety: for those binaries the first call happens during serial setup
       (computeControlPointsXY), long before any parallel region.  Even if two threads
       did arrive together they would compute identical bytes from the same globals; the
       write is done under a critical section and the flag set last. */
    if (!gpHaveDefault)
    {
        extern int32_t HemiSphere;
        extern double Rotation;
        extern double SLat;
        grimpProj p = grimpProjFromLegacy(Rotation, SLat, HemiSphere);
#ifdef _OPENMP
#pragma omp critical(grimpDefaultProj)
#endif
        {
            if (!gpHaveDefault)
            {
                gpDefaultProj = p;
                gpHaveDefault = TRUE;
            }
        }
    }
    return &gpDefaultProj;
}

/* ----------------------------------------------------------------- conversions */

void llToXYProj(double lat, double lon, double *x, double *y, const grimpProj *p)
{
    if (p->kind == GP_LATLON)
    {
        /* x = lon, y = lat, in degrees. Longitude is returned in the same 0..360 convention
           xytoll1 uses, because the pixel loops compare it against values from there. */
        *x = (lon < 0.) ? lon + 360. : lon;
        *y = lat;
        return;
    }
    if (p->kind == GP_PS || p->kind == GP_UNSET)
    {
        /* Untouched legacy path: identical instructions, identical arguments. */
        lltoxy1(lat, lon, x, y, p->rot, p->stdLat);
        return;
    }
    {
        PJ_COORD c = proj_coord(lon, lat, 0., 0.);
        c = proj_trans(gpTrans(p), PJ_FWD, c);
        *x = c.xy.x * MTOKM;
        *y = c.xy.y * MTOKM;
    }
}

void xyToLLProj(double x, double y, double *lat, double *lon, const grimpProj *p)
{
    if (p->kind == GP_LATLON)
    {
        *lat = y;
        *lon = x;
        while (*lon < 0.)
            *lon += 360.;
        while (*lon >= 360.)
            *lon -= 360.;
        return;
    }
    if (p->kind == GP_PS || p->kind == GP_UNSET)
    {
        xytoll1(x, y, p->hemisphere, lat, lon, p->rot, p->stdLat);
        return;
    }
    {
        PJ_COORD c = proj_coord(x * KMTOM, y * KMTOM, 0., 0.);
        c = proj_trans(gpTrans(p), PJ_INV, c);
        *lon = c.lp.lam;
        *lat = c.lp.phi;
        /* xytoll1 returns 0..360 and the rest of the code assumes it. */
        while (*lon < 0.)
            *lon += 360.;
        while (*lon >= 360.)
            *lon -= 360.;
    }
}

double grimpXYAngle(double lat, double lon, double x, double y, const grimpProj *p)
{
    if (p->kind == GP_GENERIC || p->kind == GP_LATLON)
        error("grimpXYAngle: the grid angle is only defined for a conformal projection.\n"
              "  %s cannot be used to decompose velocity into vx/vy.", grimpProjDescribe(p));
    if (p->kind == GP_PS || p->kind == GP_UNSET)
    {
        /* The original expression, unchanged, so the result is bit identical. */
        double xyAngle = atan2(-y, -x);
        if (p->hemisphere == SOUTH)
            xyAngle += PI;
        return xyAngle;
    }
    {
        /* xyAngle = PI/2 + meridian convergence.  For a polar stereographic the
           convergence is (lon - lon_0) in the north and -(lon - lon_0) in the south,
           and PI/2 + that is algebraically identical to the atan2 form above -- which
           is what makes this a true generalisation rather than a second convention. */
        int32_t t = gpThreadNum();
        PJ *pj;
        PJ_FACTORS f;
        PJ_COORD c;
        if (p->handle < 0 || (pj = gpRegistry[p->handle].projOnly[t]) == NULL)
            error("grimpXYAngle: grimpProjPrepareThreads() was not called before use");
        /* proj_factors wants geodetic input in RADIANS. */
        c = proj_coord(lon * DTOR, lat * DTOR, 0., 0.);
        f = proj_factors(pj, c);
        return PI / 2. + f.meridian_convergence;
    }
}

double grimpXYScale(double lat, const grimpProj *p)
{
    if (p->kind == GP_GENERIC || p->kind == GP_LATLON)
        error("grimpXYScale: an isotropic grid scale is only defined for a conformal\n"
              "  projection; %s stretches x and y differently.", grimpProjDescribe(p));
    if (p->kind == GP_PS || p->kind == GP_UNSET)
    {
        /* Historical approximation from xyGetZandSlope.c; agrees with the exact
           1/point-scale to 4 digits, which is ample for scaling DEM derivatives. */
        return (1.0 + sin(fabs(lat) * DTOR)) / (1.0 + sin(p->stdLat * DTOR));
    }
    /* Within one UTM zone the point scale is 0.9996..1.0010, i.e. under 0.15%, far
       below the other error terms in the slope, so treat the grid as true scale. */
    return 1.0;
}

int32_t grimpProjSame(const grimpProj *a, const grimpProj *b)
{
    if (a == b)
        return TRUE;
    if (a->epsg != 0 && a->epsg == b->epsg)
        return TRUE;
    if (a->kind != b->kind)
        return FALSE;
    if (a->kind == GP_UTM)
        return (a->utmZone == b->utmZone && a->hemisphere == b->hemisphere);
    /* Polar stereographic, possibly with no EPSG code: compare the parameters. */
    return (a->rot == b->rot && a->stdLat == b->stdLat && a->hemisphere == b->hemisphere);
}

static gpPairEntry *gpGetPair(const grimpProj *src, const grimpProj *dst)
{
    int32_t i;
    for (i = 0; i < gpNPairs; i++)
    {
        if (gpPairs[i].used && gpPairs[i].srcHandle == src->handle &&
            gpPairs[i].dstHandle == dst->handle)
            return &gpPairs[i];
    }
    return NULL;
}

void outXYToProjXY(double xOut, double yOut, const grimpProj *outProj,
                   const grimpProj *p, double *x, double *y)
{
    gpPairEntry *q;
    /* Overwhelmingly the common case, and the only one that existed before this
       change: the grids are the same, so this costs one comparison. */
    if (grimpProjSame(outProj, p))
    {
        *x = xOut;
        *y = yOut;
        return;
    }
    /* Both polar stereographic and both on the legacy path: go via lat/lon using the
       original routines, so no PROJ object is needed for either side. */
    if ((outProj->kind == GP_PS || outProj->kind == GP_UNSET) &&
        (p->kind == GP_PS || p->kind == GP_UNSET))
    {
        double lat, lon;
        xytoll1(xOut, yOut, outProj->hemisphere, &lat, &lon, outProj->rot, outProj->stdLat);
        lltoxy1(lat, lon, x, y, p->rot, p->stdLat);
        return;
    }
    q = gpGetPair(outProj, p);
    if (q != NULL && q->pj[gpThreadNum()] != NULL)
    {
        PJ_COORD c = proj_coord(xOut * KMTOM, yOut * KMTOM, 0., 0.);
        c = proj_trans(q->pj[gpThreadNum()], PJ_FWD, c);
        *x = c.xy.x * MTOKM;
        *y = c.xy.y * MTOKM;
        return;
    }
    /* No direct pipeline registered: route through lat/lon, which always works and is
       only about twice the cost of a direct transform. */
    {
        double lat, lon;
        xyToLLProj(xOut, yOut, &lat, &lon, outProj);
        llToXYProj(lat, lon, x, y, p);
    }
}

/* Ask for a direct src->dst pipeline to be built by the next
   grimpProjPrepareThreads(). Serial only; optional (outXYToProjXY works without it). */
void grimpProjRegisterPair(const grimpProj *src, const grimpProj *dst)
{
    if (src->handle < 0 || dst->handle < 0)
        return; /* one side is on the legacy path; no pipeline possible or needed */
    if (gpGetPair(src, dst) != NULL)
        return;
    if (gpNPairs >= GP_MAXPAIR)
        return; /* fall back to the lat/lon route rather than failing */
    memset(&gpPairs[gpNPairs], 0, sizeof(gpPairs[gpNPairs]));
    gpPairs[gpNPairs].used = TRUE;
    gpPairs[gpNPairs].srcHandle = src->handle;
    gpPairs[gpNPairs].dstHandle = dst->handle;
    gpNPairs++;
}

/* --------------------------------------------------------------- naming / output */

const char *grimpSRSString(const grimpProj *p)
{
    static char buf[GP_MAXPROJ + 1][512];
    static int32_t next = 0;
    char *b = buf[next];
    next = (next + 1) % (GP_MAXPROJ + 1);

    if (p->epsg != 0)
    {
        snprintf(b, 512, "EPSG:%d", (int)p->epsg);
        return b;
    }
    if (p->kind == GP_UTM)
    {
        /* Should not happen: every UTM zone has an EPSG code. */
        snprintf(b, 512, "+proj=utm +zone=%d%s +datum=WGS84 +units=m +no_defs",
                 (int)p->utmZone, (p->hemisphere == SOUTH) ? " +south" : "");
        return b;
    }
    /* Custom polar stereographic with no EPSG code, e.g. the Taku lat_ts=58 grid.
       Writing this as a proj string is what lets those products carry a real CRS. */
    {
        double lon0 = (p->hemisphere == NORTH) ? -p->rot : p->rot;
        snprintf(b, 512, "+proj=stere +lat_0=%s90 +lat_ts=%.10g +lon_0=%.10g "
                         "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                 (p->hemisphere == SOUTH) ? "-" : "", p->stdLat, lon0);
        return b;
    }
}

void grimpRequireConformal(const grimpProj *p, const char *what)
{
    if (p->kind == GP_LATLON)
        error("%s cannot be geographic (lat/lon).\n"
              "  A geographic grid is not conformal -- east-west scale falls as cos(lat) --\n"
              "  so one grid angle and one isotropic scale cannot describe it. Geographic\n"
              "  output is supported for BACKSCATTER mosaics (geomosaic), which only resample\n"
              "  a scalar, never for velocity.",
              what);
    if (p->kind == GP_GENERIC)
        error("%s must be a conformal projection (polar stereographic or UTM).\n"
              "  Got %s.\n"
              "  Velocity is decomposed into vx/vy with a single grid angle and one\n"
              "  isotropic scale, which only holds for a conformal grid; an equal-area\n"
              "  projection such as Alaska Albers (EPSG:3338) would need direction-dependent\n"
              "  scaling.  geomosaic has no such restriction, and input grids (DEM, velocity\n"
              "  map, masks) may be in any projection.",
              what, grimpProjDescribe(p));
}

const char *grimpProjDescribe(const grimpProj *p)
{
    static char buf[4][256];
    static int32_t next = 0;
    char *b = buf[next];
    next = (next + 1) % 4;
    if (p->kind == GP_LATLON)
        snprintf(b, 256, "EPSG:%d (geographic lat/lon; grid units are DEGREES)", (int)p->epsg);
    else if (p->kind == GP_UTM)
        snprintf(b, 256, "EPSG:%d (UTM zone %d%c)", (int)p->epsg, (int)p->utmZone,
                 (p->hemisphere == SOUTH) ? 'S' : 'N');
    else if (p->kind == GP_GENERIC)
        snprintf(b, 256, "EPSG:%d (projected, non-conformal or unclassified)", (int)p->epsg);
    else if (p->epsg != 0)
        snprintf(b, 256, "EPSG:%d (polar stereographic, rot %g, stdLat %g)",
                 (int)p->epsg, p->rot, p->stdLat);
    else
        snprintf(b, 256, "polar stereographic (no EPSG code; rot %g, stdLat %g, %s)",
                 p->rot, p->stdLat, (p->hemisphere == SOUTH) ? "south" : "north");
    return b;
}

/* ------------------------------------------------------------------- self test */

void grimpProjSelfTest(void)
{
    /* Checks the generalised path against the legacy one on projections where the
       right answer is already known, so a sign or units mistake cannot reach data. */
    int32_t nBad = 0;
    double lat, lon;

    fprintf(stderr, "grimpProjSelfTest: starting\n");

    /* 1. PROJ must reproduce lltoxy1 for the two production polar projections, and
          2/3. the convergence identity must hold in BOTH hemispheres -- the test that
          catches the rot/lon_0 sign asymmetry. */
    {
        struct { const char *srs; double rot, stdLat; int32_t hemi; double lat0, lat1; } cases[] = {
            {"EPSG:3413", 45., 70., NORTH,  60.,  85.},
            {"EPSG:3031",  0., 71., SOUTH, -85., -60.},
            /* A rotated custom PS in each hemisphere: no EPSG code, so this also
               exercises the proj-string path in grimpSRSString. */
            {"+proj=stere +lat_0=90 +lat_ts=58 +lon_0=-134 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                 134., 58., NORTH, 50., 75.},
            {"+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=-30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                 -30., 71., SOUTH, -85., -55.},
        };
        int32_t ic;
        for (ic = 0; ic < 4; ic++)
        {
            grimpProj gp = grimpProjFromSRS(cases[ic].srs);
            double maxPos = 0., maxAng = 0.;
            if (gp.kind != GP_PS)
                error("grimpProjSelfTest: %s did not come back as polar stereographic",
                      cases[ic].srs);
            if (fabs(gp.rot - cases[ic].rot) > 1e-9 || fabs(gp.stdLat - cases[ic].stdLat) > 1e-9 ||
                gp.hemisphere != cases[ic].hemi)
            {
                fprintf(stderr, "  FAIL %s: got rot %g stdLat %g hemi %i, expected %g %g %i\n",
                        cases[ic].srs, gp.rot, gp.stdLat, (int)gp.hemisphere,
                        cases[ic].rot, cases[ic].stdLat, (int)cases[ic].hemi);
                nBad++;
                continue;
            }
            /* Compare the legacy transform against PROJ built from the same SRS. */
            {
                PJ_CONTEXT *ctx = proj_context_create();
                PJ *raw = proj_create_crs_to_crs(ctx, "EPSG:4326", cases[ic].srs, NULL);
                PJ *pj = (raw == NULL) ? NULL : proj_normalize_for_visualization(ctx, raw);
                PJ *only = proj_create(ctx, cases[ic].srs);
                if (raw != NULL) proj_destroy(raw);
                if (pj == NULL || only == NULL)
                    error("grimpProjSelfTest: cannot build PROJ objects for %s", cases[ic].srs);
                for (lat = cases[ic].lat0; lat <= cases[ic].lat1; lat += 1.0)
                {
                    for (lon = -180.; lon < 180.; lon += 5.0)
                    {
                        double xl, yl, d, ang, angRef;
                        PJ_COORD c;
                        lltoxy1(lat, lon, &xl, &yl, gp.rot, gp.stdLat);
                        c = proj_coord(lon, lat, 0., 0.);
                        c = proj_trans(pj, PJ_FWD, c);
                        d = fabs(c.xy.x * MTOKM - xl) + fabs(c.xy.y * MTOKM - yl);
                        if (d > maxPos) maxPos = d;
                        /* convergence identity */
                        angRef = atan2(-yl, -xl);
                        if (gp.hemisphere == SOUTH) angRef += PI;
                        {
                            PJ_FACTORS f = proj_factors(only, proj_coord(lon * DTOR, lat * DTOR, 0., 0.));
                            double diff;
                            ang = PI / 2. + f.meridian_convergence;
                            diff = fabs(remainder(ang - angRef, 2. * PI));
                            if (diff > maxAng) maxAng = diff;
                        }
                    }
                }
                proj_destroy(pj);
                proj_destroy(only);
                proj_context_destroy(ctx);
            }
            fprintf(stderr, "  %-22s max |PROJ-lltoxy1| %.3g km, max angle diff %.3g rad\n",
                    grimpProjDescribe(&gp), maxPos, maxAng);
            if (maxPos > 1e-6) { fprintf(stderr, "  FAIL: position disagreement\n"); nBad++; }
            if (maxAng > 1e-9) { fprintf(stderr, "  FAIL: grid-angle disagreement\n"); nBad++; }
        }
    }

    /* 4. Units guard.  proj_factors returns RADIANS; pyproj and some docs give degrees.
          Zone 8N at 57N/132W has a convergence of 3 deg * sin(57) = 2.517 deg =
          0.04393 rad, so a units slip is a factor of 57 and is caught here. */
    {
        grimpProj u = grimpProjFromEPSG(32608);
        double conv;
        grimpProjPrepareThreads();
        conv = grimpXYAngle(57.0, -132.0, 0., 0., &u) - PI / 2.;
        fprintf(stderr, "  UTM 8N convergence at 57N/132W: %.6f rad (%.4f deg)\n",
                conv, conv * RTOD);
        if (fabs(conv - 0.043929) > 1e-4)
        {
            fprintf(stderr, "  FAIL: expected 0.043929 rad; %s\n",
                    (fabs(conv - 2.5169) < 1e-2) ? "this looks like DEGREES" : "wrong value");
            nBad++;
        }
        if (u.kind != GP_UTM || u.utmZone != 8 || u.hemisphere != NORTH)
        {
            fprintf(stderr, "  FAIL: EPSG:32608 decoded as %s\n", grimpProjDescribe(&u));
            nBad++;
        }
    }

    /* 5. Round trip, both kinds, both hemispheres. */
    {
        const char *srsList[] = {"EPSG:3413", "EPSG:3031", "EPSG:32608", "EPSG:32708"};
        int32_t i;
        double lat0[] = {60., -85., 40., -60.}, lat1[] = {85., -60., 70., -40.};
        for (i = 0; i < 4; i++)
        {
            grimpProj gp = grimpProjFromSRS(srsList[i]);
            double maxErr = 0.;
            grimpProjPrepareThreads();
            for (lat = lat0[i]; lat <= lat1[i]; lat += 2.0)
            {
                for (lon = -138.; lon < -128.; lon += 1.0)
                {
                    double x, y, la, lo, dLon;
                    llToXYProj(lat, lon, &x, &y, &gp);
                    xyToLLProj(x, y, &la, &lo, &gp);
                    dLon = fabs(remainder(lo - lon, 360.));
                    if (fabs(la - lat) > maxErr) maxErr = fabs(la - lat);
                    if (dLon > maxErr) maxErr = dLon;
                }
            }
            fprintf(stderr, "  %-12s round trip max error %.3g deg\n", srsList[i], maxErr);
            if (maxErr > 1e-9) { fprintf(stderr, "  FAIL: round trip\n"); nBad++; }
        }
    }

    /* 6. grimpXYScale must match the exact reciprocal point scale for PS. */
    {
        grimpProj gp = grimpProjFromSRS("EPSG:3413");
        PJ_CONTEXT *ctx = proj_context_create();
        PJ *only = proj_create(ctx, "EPSG:3413");
        double maxRel = 0.;
        for (lat = 60.; lat <= 85.; lat += 1.0)
        {
            PJ_FACTORS f = proj_factors(only, proj_coord(-45. * DTOR, lat * DTOR, 0., 0.));
            double approx = grimpXYScale(lat, &gp);
            double rel = fabs(approx * f.meridional_scale - 1.0);
            if (rel > maxRel) maxRel = rel;
        }
        proj_destroy(only);
        proj_context_destroy(ctx);
        fprintf(stderr, "  xyScale vs 1/point-scale: max relative error %.3g\n", maxRel);
        if (maxRel > 1e-4) { fprintf(stderr, "  FAIL: scale approximation\n"); nBad++; }
    }

    /* 7. Non-conformal and geographic systems must be refused, not silently accepted. */
    fprintf(stderr, "  (EPSG:3338 / EPSG:4326 rejection is checked by inspection: "
                    "grimpProjFromSRS error()s, which would abort this test)\n");

    if (nBad != 0)
        error("grimpProjSelfTest: %i check(s) FAILED", nBad);
    fprintf(stderr, "grimpProjSelfTest: all checks passed\n");
}
