#ifndef GRIMPPROJ_H
#define GRIMPPROJ_H
#include <stdint.h>

/*
  Map-projection descriptor for the GrIMP mosaickers.

  Background
  ----------
  Historically every projected coordinate in this code went through lltoxy1()/xytoll1()
  (common/lltoxy1.c, common/xytoll1.c): a hand-coded Snyder polar stereographic
  parameterised by a rotation `dlam` and a standard parallel `slat`, working in
  KILOMETRES on a hard-coded WGS84 ellipsoid.  That covers any polar stereographic but
  nothing else, and the output side could only name EPSG:3413 and EPSG:3031.

  grimpProj generalises this to polar stereographic *and* UTM -- both conformal, so the
  existing "rotate by the grid angle, scale isotropically" machinery stays valid.
  Equal-area projections (e.g. Alaska Albers, EPSG:3338) are anisotropic and are
  deliberately NOT supported.

  Backwards compatibility
  -----------------------
  For kind == GP_PS every routine here delegates to the original lltoxy1()/xytoll1()
  with the original arguments, so existing polar runs execute exactly the same
  instructions and produce byte-identical output.  PROJ is reached only for GP_UTM.
  That includes custom polar stereographics with no EPSG code, such as the Taku
  (lat_ts=58, lon_0=-134) grid already in production.

  Units and conventions, unchanged from lltoxy1/xytoll1
  ----------------------------------------------------
  - x, y are in KILOMETRES.
  - longitude is returned in 0..360.
  - `stdLat` is ALWAYS POSITIVE; the hemisphere is carried separately.
  - `rot` is lltoxy1's `dlam`, which maps to the projection's central meridian
    ASYMMETRICALLY, because lltoxy1 negates the longitude in the southern hemisphere:
        north:  lon_0 = -rot
        south:  lon_0 = +rot
    Getting this wrong only shows up for a *rotated* southern polar stereographic
    (EPSG:3031 has lon_0 = 0), which is why it went unnoticed for years.

  Threading
  ---------
  PROJ transformation objects are not thread safe, so one is needed per OpenMP thread.
  The handles deliberately do NOT live in this struct: xyDEM (which embeds a grimpProj)
  is passed BY VALUE throughout the code -- computeXYangle(), interpXYDEM(), and
  computeB()'s `xyDEM dem = *xydem;` all copy it -- so a pointer cached here would be
  cached into a temporary copy and recreated, and leaked, on every call.  Instead the
  struct holds an integer `handle` into a private registry inside grimpProj.c.

  Call grimpProjPrepareThreads() once, SERIALLY, before any parallel region that will
  convert coordinates.  This mirrors the existing rule for svAzOffset()/svInterpBnBp()
  (see the root CLAUDE.md): lazily creating this state from several threads races.
*/

typedef enum
{
    GP_UNSET = 0, /* not yet resolved; treated as polar stereographic */
    GP_PS = 1,    /* polar stereographic, incl. custom (possibly no EPSG code) */
    GP_UTM = 2    /* universal transverse Mercator, one zone */
} grimpProjKind;

typedef struct
{
    int32_t epsg;          /* EPSG code, or 0 when the grid has none */
    grimpProjKind kind;
    double rot;            /* lltoxy1 dlam; see the sign note above */
    double stdLat;         /* always positive */
    int32_t hemisphere;    /* NORTH or SOUTH (geocode.h) */
    int32_t utmZone;       /* GP_UTM only, 1..60; 0 otherwise */
    int32_t handle;        /* index into grimpProj.c's registry, -1 if no PROJ needed */
} grimpProj;

/* ---- constructors (SERIAL only) ------------------------------------------------- */

/* From the legacy triple.  Always GP_PS, never needs PROJ.  A stdLat of -91 (the
   "unset" sentinel carried by the SLat global) is replaced by the historical
   hemisphere default, 70 north / 71 south, exactly as the open-coded call sites did. */
grimpProj grimpProjFromLegacy(double rot, double stdLat, int32_t hemisphere);

/* From anything OSRSetFromUserInput accepts: "EPSG:32608", a proj string, or WKT.
   Errors out unless the result is polar stereographic or UTM. */
grimpProj grimpProjFromSRS(const char *userInput);

/* Convenience wrapper: grimpProjFromSRS("EPSG:<epsg>"). */
grimpProj grimpProjFromEPSG(int32_t epsg);

/* Create the per-thread PROJ objects for every projection registered so far.
   Idempotent, and a cheap no-op when nothing new has been registered. */
void grimpProjPrepareThreads(void);

/* ---- the process-wide "current output projection" ------------------------------- */
/* For the handful of common/ routines that have no outputImage in scope
   (parseInputFile.c, llToImageNew.c's checkLL/computeFootprintPolygon, getIrregData.c).
   Set once at startup, read-only afterwards. */
void grimpSetDefaultProj(const grimpProj *p);
const grimpProj *grimpDefaultProj(void);

/* ---- conversions (thread safe after grimpProjPrepareThreads) -------------------- */

/* lat/lon (degrees) -> x/y (KM) */
void llToXYProj(double lat, double lon, double *x, double *y, const grimpProj *p);

/* x/y (KM) -> lat/lon (degrees, lon in 0..360) */
void xyToLLProj(double x, double y, double *lat, double *lon, const grimpProj *p);

/* The grid angle the mosaickers difference against the SAR heading:
      xyAngle = PI/2 + meridian convergence.
   For GP_PS this is the original atan2(-y,-x) (+PI in the south) evaluated on the
   supplied x,y, bit for bit.  For GP_UTM it comes from PROJ's meridian convergence.
   x,y must be in this projection, in km, and correspond to lat,lon. */
double grimpXYAngle(double lat, double lon, double x, double y, const grimpProj *p);

/* Ratio of grid distance to true ground distance, used to normalise DEM slopes.
   GP_PS keeps the historical approximation (1+sin|lat|)/(1+sin stdLat); GP_UTM
   returns 1.0, which is within 0.15% over a single zone. */
double grimpXYScale(double lat, const grimpProj *p);

/* Convert an x/y (km) in the OUTPUT projection into the x/y (km) of projection `p`.
   A no-op, with no PROJ call at all, when the two describe the same grid -- which is
   every case that existed before this change. */
void outXYToProjXY(double xOut, double yOut, const grimpProj *outProj,
                   const grimpProj *p, double *x, double *y);

/* Optional: ask for a direct src->dst PROJ pipeline to be built by the next
   grimpProjPrepareThreads(), rather than routing through lat/lon.  Serial only.
   outXYToProjXY works correctly whether or not this was called. */
void grimpProjRegisterPair(const grimpProj *src, const grimpProj *dst);

/* Something OSRSetFromUserInput understands, for the GeoTIFF/VRT writers:
   "EPSG:nnnnn" when known, otherwise a proj string for a custom polar stereographic.
   The returned string is owned by grimpProj.c; do not free it. */
const char *grimpSRSString(const grimpProj *p);

/* True when the two descriptors name the same grid. */
int32_t grimpProjSame(const grimpProj *a, const grimpProj *b);

/* Short human-readable description for log lines, e.g. "EPSG:32608 (UTM zone 8N)". */
const char *grimpProjDescribe(const grimpProj *p);

/* Consistency checks. error()s on failure. */
void grimpProjSelfTest(void);

#endif /* GRIMPPROJ_H */
