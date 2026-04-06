# common — Shared Library Reference

All non-static functions in `mosaicSource/common/`. Grouped by category.

---

## Recurring Variables

These symbols appear throughout the library with consistent meaning:

| Symbol | Type | Units | Meaning |
|--------|------|-------|---------|
| `Re` | double | km | Local Earth radius (spherical approximation at pixel latitude) |
| `ReH` | double | km | Satellite orbital radius = $R_e + H$ (varies along azimuth) |
| `ReZ` | double | km | Surface radius = $R_e + h$ (where $h$ is the terrain height) |
| `R`, `Range` | double | m | Slant range from satellite to target |
| `range` | double | pixels | Fractional range pixel coordinate in image |
| `azimuth` | double | pixels | Fractional azimuth pixel coordinate in image |
| `theta` | double | rad | Look angle at the satellite (angle between nadir and look vector) |
| `thetaD` | double | rad | Deviation of look angle from scene centre look angle: $\theta - \theta_c$ |
| `thetaC` | double | rad | Centre look angle of the image |
| `psi` | double | rad | Incidence angle at the surface (angle between look vector and local vertical) |
| `aPsi`, `dPsi` | double | rad | Incidence angles for ascending / descending images |
| `lat`, `lon` | double | degrees | Geodetic latitude (+N) and longitude (0–360 E) |
| `x`, `y` | double | km | Polar-stereographic map coordinates |
| `h`, `hWGS84` | double | m | Surface elevation on the WGS84 ellipsoid |
| `zSp` | double | m | Surface elevation referenced to a local sphere of radius `Re` |
| `Bn`, `Bp` | double | m | Normal and parallel interferometric baseline components |
| `dBn`, `dBp` | double | m | Along-track linear derivatives of `Bn`, `Bp` |
| `dBnQ`, `dBpQ` | double | m | Along-track quadratic derivatives of `Bn`, `Bp` |
| `normAzimuth` | double | — | Normalised azimuth position $x = (\text{az} - N/2) / N \in [-0.5, 0.5]$ |
| `nDays` | float | days | Temporal baseline between image acquisitions |
| `lambda` | double | m | SAR carrier wavelength |
| `twok` | double | rad/m | $4\pi/\lambda$ (twice the wavenumber) |
| `HemiSphere` | int | — | `NORTH` or `SOUTH` — selects polar-stereographic hemisphere |
| `Rotation` | double | degrees | Grid rotation angle (0° longitude offset for PS projection) |
| `fl` | float | m | Feathering length for image edge blending |
| `LARGEINT` | macro | — | Large sentinel value (~2×10⁹) used as no-data |
| `MTOKM` | macro | — | Conversion factor metres → km (0.001) |
| `DTOR` | macro | — | Degrees to radians ($\pi/180$) |
| `RTOD` | macro | — | Radians to degrees ($180/\pi$) |

---

## Coordinate Transformations

### `lltoxy1.c` / `lltoxy.c`

**`lltoxy1`**(alat, alon, \*x, \*y, dlam, slat) → void  
**`lltoxy`**(alat, alon, \*x, \*y, dlam) → void

Convert geodetic latitude/longitude to polar stereographic (x, y) coordinates in km.
`lltoxy` fixes the standard parallel at 70°; `lltoxy1` accepts it as `slat`.
Uses the WGS84 ellipsoid (a = 6378.137 km, e² = 0.006694380).

The conformal latitude factor:

$$
t(\varphi) = \tan\!\left(\frac{\pi}{4} - \frac{\varphi}{2}\right)
             \left(\frac{1 - e\sin\varphi}{1 + e\sin\varphi}\right)^{e/2}
$$

Scale factor at the standard parallel $\varphi_s$:

$$
c_m = \frac{\cos\varphi_s}{\sqrt{1 - e^2\sin^2\varphi_s}}, \qquad
\rho = \frac{a\, c_m\, t(\varphi)}{t(\varphi_s)}
$$

Map coordinates (with optional grid rotation $\Delta\lambda$):

$$
x = \rho\sin(\lambda + \Delta\lambda), \qquad
y = -\rho\cos(\lambda + \Delta\lambda)
$$

Signs are flipped for the southern hemisphere. Reference: Snyder (1982), USGS Bulletin 1532.

---

### `xytoll1.c` / `xytoll.c`

**`xytoll1`**(x, y, hemi, \*alat, \*alon, dlam, slat) → void  
**`xytoll`**(x, y, hemi, \*alat, \*alon, dlam) → void

Inverse of `lltoxy1`/`lltoxy`. Converts polar stereographic (x, y) in km to geodetic lat/lon.
`xytoll` fixes the standard parallel at 70°; `xytoll1` accepts it as `slat`.

$$
\rho = \sqrt{x^2 + y^2}, \qquad
t = \frac{\rho\, t(\varphi_s)}{a\, c_m}
$$

Conformal latitude:

$$
\chi = \frac{\pi}{2} - 2\arctan t
$$

Geodetic latitude by series expansion (Snyder Eq. 3-5):

$$
\varphi = \chi
  + \left(\frac{e^2}{2} + \frac{5e^4}{24} + \frac{e^6}{12}\right)\sin 2\chi
  + \left(\frac{7e^4}{48} + \frac{29e^6}{240}\right)\sin 4\chi
  + \frac{7e^6}{120}\sin 6\chi
$$

Longitude: $\lambda = \arctan(x / {-y}) - \Delta\lambda$.

---

### `llToImageNew.c`

**`initllToImageNew`**(inputImage) → void  
Precomputes satellite state vectors, Earth radius arrays, and look-angle geometry tables needed
for fast lat/lon → range/azimuth conversion. Must be called once before `llToImageNew`.  
*Calls:* `polintVec`, `psiRReZReH`, `thetaRReZReH`, `earthRadiusWGS84`

**`llToImageNew`**(lat, lon, h, \*range, \*azimuth, inputImage) → void  
Converts geodetic (lat, lon, height) to SAR image (range, azimuth) pixel coordinates using
polynomial-interpolated satellite state vectors. Iterates on azimuth time until the Doppler
zero-Doppler condition is satisfied, then solves for slant range:  
*Calls:* `llToECEF`, `polintVec`

$$
R = \sqrt{(x_s - x_t)^2 + (y_s - y_t)^2 + (z_s - z_t)^2}
$$

where $(x_s, y_s, z_s)$ is the interpolated satellite ECEF position and $(x_t, y_t, z_t)$ is the
target ECEF position at height $h$.

**`llToECEF`**(lat, lon, h, \*x, \*y, \*z) → void  
Converts geodetic coordinates to Earth-Centred Earth-Fixed (ECEF) Cartesian coordinates
using the WGS84 ellipsoid:

$$
N = \frac{a}{\sqrt{1 - e^2\sin^2\varphi}}, \qquad
\begin{cases}
x = (N + h)\cos\varphi\cos\lambda \\
y = (N + h)\cos\varphi\sin\lambda \\
z = \left(\dfrac{b^2}{a^2} N + h\right)\sin\varphi
\end{cases}
$$

---

### `groundRangeToLLNew.c`

**`groundRangeToLLNew`**(groundRange, azimuth, \*lat, \*lon, inputImage, recycle) → double  
Converts SAR ground range and azimuth pixel coordinates to geodetic lat/lon by iterating
on an ellipsoidal Earth model. Returns the computed slant range. The `recycle` flag reuses
the previously interpolated state vectors for efficiency when called repeatedly at the same
azimuth time.  
*Calls:* `polintVec`, `smlocateZD`, `psiRReZReH`, `getReH`

---

### `rangeAzimuthToLL.c`

**`rangeAzimuthToLL`**(rg, range, iFloat, rhoSp, ReH, Re, \*lat, \*lon, \*hWGS, inputImage, xyDem, tol, step) → int32_t  
Iterative inversion of slant range + azimuth to lat/lon/height on the WGS84 ellipsoid,
using a geoid/DEM model. Starts from a first guess and steps in lat/lon until the residual
in slant range and along-track position is below `tol`. Returns iteration count.

---

### `smlocateZD.c`

**`smlocateZD`**(rxs, rys, rzs, rvsx, rvsy, rvsz, rsl, \*lat, \*lon, lookDir, trgalt) → void  
*Calls:* `earthRadiusWGS84`  
Locates a target on the WGS84 ellipsoid given a spacecraft ECEF position and velocity and
a slant range `rsl` (km). Iteratively solves for the look angle $\theta$ using:

$$
R = \sqrt{R_E^2 + R_{EH}^2 - 2 R_E R_{EH} \cos\rho}
$$

where $R_E = R_e + h_\text{tgt}$, $R_{EH} = |\mathbf{r}_s|$, and $\rho$ is the Earth central
angle. The look direction (LEFT/RIGHT) selects the correct solution branch.

---

## Earth Geometry

### `earthRadiusFunctions.c`

**`earthRadius`**(lat, rp, re) → double  
Radial distance from Earth's centre to the ellipsoid surface (km) at geodetic latitude
$\varphi$ (radians), for semi-major axis $a$ and semi-minor axis $b$:

$$
N = \frac{a^2}{\sqrt{(a\cos\varphi)^2 + (b\sin\varphi)^2}}, \qquad
R = \sqrt{(N\cos\varphi)^2 + \!\left(\tfrac{b^2}{a^2}N\sin\varphi\right)^2}
$$

**`earthRadiusWGS84`**(lat) → double  
Calls `earthRadius` with WGS84 constants ($a$ = 6378.137 km, $b$ = 6356.752 km).

**`earthRadiusCurvatureWGS84`**(lat) → double  
Prime-vertical radius of curvature $N$ at latitude $\varphi$:

$$
N = \frac{a^2}{\sqrt{(a\cos\varphi)^2 + (b\sin\varphi)^2}}
$$

Used for converting geodetic to ECEF coordinates.

---

### `initRoutines.c` (geometry helpers — declared in `common.h`)

**`psiRReZReH`**(R, ReZ, ReH) → double  
Incidence angle $\psi$ at the surface from slant range $R$, surface radius $R_{eZ} = R_e + h$,
and satellite radius $R_{eH} = R_e + H$, using the spherical law of sines:

$$
\theta = \arccos\!\left(\frac{R^2 + R_{eH}^2 - R_{eZ}^2}{2\,R\,R_{eH}}\right), \qquad
\psi = \arcsin\!\left(\frac{R_{eH}}{R_{eZ}}\sin\theta\right)
$$

**`thetaRReZReH`**(R, ReZ, ReH) → double  
Look angle $\theta$ at the satellite (intermediate result used by `psiRReZReH`).

**`rhoRReZReH`**(R, ReZ, ReH) → double  
Earth central angle $\rho$ between sub-satellite point and target:

$$
\rho = \arccos\!\left(\frac{R^2 - R_{eH}^2 - R_{eZ}^2}{-2\,R_{eZ}\,R_{eH}}\right)
$$

**`getReH`**(cP, inputImage, azimuth) → double  
Returns the azimuth-varying satellite radius $R_{eH}(\text{az})$ by nearest-neighbour lookup
in a precomputed table. Falls back to a fixed value if the table is absent.

**`slantRange`**(rho, ReZ, ReH) → double  
Computes slant range from Earth central angle:

$$
R = \sqrt{R_{eH}^2 + R_{eZ}^2 - 2\,R_{eZ}\,R_{eH}\cos\rho}
$$

---

### `computeHeading.c`

**`computeHeading`**(lat, lon, z, inputImage, cP) → double  
Computes the SAR flight-track heading angle (radians, clockwise from north) at a given
surface point by evaluating the satellite velocity vector at the corresponding azimuth time
and projecting onto the local north/east plane.

---

### `computeXYangle.c`

**`computeXYangle`**(lat, lon, \*xyAngle, xydem) → void  
**`computeXYangleNoDem`**(lat, lon, \*xyAngle, stdLat) → void  
Compute the angle (radians) between the polar-stereographic x-axis and geographic north
at a given lat/lon. Used to rotate between map and SAR coordinate frames. The DEM version
reads the standard latitude and grid rotation from the `xyDEM` structure; the `NoDem`
version takes `stdLat` directly.

---

## DEM Access

### `readXYDEM.c`

**`readXYDEM`**(xyFile, xydem) → void  
Reads a complete XY polar-stereographic DEM (binary float, big-endian) and its `.geodat`
metadata file into an `xyDEM` structure.

**`readXYDEMcrop`**(xyFile, xyDEM, xmin, xmax, ymin, ymax) → void  
Reads only the subset of the DEM that overlaps the bounding box [xmin, xmax] × [ymin, ymax]
(km). Reduces memory for large DEMs when only a regional subset is needed.

**`readXYDEMGeoInfo`**(xyFile, xydem, resetProjection) → void  
Reads geometric metadata only (origin, pixel spacing, size, projection) without loading
the elevation data.

**`readXYGeoInfoGDAL`**(xyFile, xyImage, type) → void  
**`readXYProjInfoGDAL`**(xyFile, obj, type) → void  
Read geometric and projection metadata from GDAL-compatible files (GeoTIFF, VRT, NetCDF).
Sets up polar stereographic parameters from the file's coordinate reference system.

**`readXYVel`**(xyvel, velFile) → void  
**`readXYCropVel`**(xyvel, velFile, xmin, xmax, ymin, ymax) → void  
Read vx/vy velocity map pairs (full or cropped) into an `xyVEL` structure.

---

### `getXYHeight.c`

**`getXYHeight`**(lat, lon, xydem, Re, heightFlag) → double  
Bilinearly interpolates the XY DEM at (lat, lon) and returns either the ellipsoidal height
(`ELLIPSOIDAL`) or converts to spherical height above a sphere of radius `Re`
(`SPHERICAL`): $h_\text{sph} = h_\text{WGS84} + R_e(\varphi) - R_e$.  
*Calls:* `lltoxy1`, `interpXYDEM`, `earthRadiusWGS84`

---

### `getHeight.c`

**`getHeight`**(lat, lon, dem, Re, heightFlag) → double  
Same as `getXYHeight` but for a lat/lon-indexed `demStructure` rather than an XY-projected
DEM.

---

### `interpXYDEM.c`

**`interpXYDEM`**(x, y, xydem) → double  
Bilinear interpolation of elevation directly in XY (km) coordinates, without lat/lon
conversion. Faster than `getXYHeight` when x/y are already available.

---

### `interpVCorrect.c`

**`interpVCorrect`**(x, y, vCorrect) → double  
Bilinear interpolation of a vertical velocity correction field (e.g., submergence/emergence
rate) at XY coordinates. Used to adjust tidal or dynamic topography on floating ice.

---

### `interpTideDiff.c`

**`interpTideDiff`**(x, y, xydem) → double  
Bilinear interpolation of a tide-difference map at XY coordinates. The map stores the
difference in tidal displacement between the two SAR acquisition epochs, used to correct
phase or offset measurements over floating ice.

---

### `xyGetZandSlope.c`

**`xyGetZandSlope`**(lat, lon, x, y, \*zSp, \*zWGS84, \*da, \*dr, cP, vhParam, image) → void  
Retrieves DEM height at (lat, lon) and computes the local surface slope projected onto
the SAR range ($dr$) and azimuth ($da$) directions. The slope terms are used in the
speckle-tracking velocity inversion to account for terrain-induced apparent motion.

---

### `getShelfMask.c` / `readShelf.c`

**`getShelfMask`**(shelfMask, x, y) → unsigned char  
Returns the mask value (`GROUNDED` or `SHELF`) for XY coordinates, by nearest-neighbour
lookup in the shelf mask grid. Used to gate tidal corrections to floating ice only.

**`readGeodatFile`**(geodatFile, \*x0, \*y0, \*deltaX, \*deltaY, \*sx, \*sy) → void  
Parses a `.geodat` file and returns grid origin (km), pixel spacing (km), and raster size.

**`readShelf`**(outputImage, shelfMaskFile) → void  
Reads a binary shelf mask and its `.geodat` into the output image structure.

---

## InSAR Phase and Baseline

### `computePhiZ.c`

**`computePhiZ`**(phiZ, azimuth, vhParam, phaseImage, thetaD, Range, ReH, ReHfixed, Re, thetaCfixed, \*phaseError) → void  
*Calls:* `thetaRReZReH`  
Computes the topographic phase contribution $\phi_Z$ and associated phase error for a
pixel at slant range $R$ and look angle $\theta_D$, given the along-track baseline
polynomial model. The baseline components are evaluated at normalised azimuth position
$x = (\text{az} - N/2)/N$:

$$
B_n(x) = B_{n0} + \delta B_n\, x + \delta B_{nQ}\, x^2, \qquad
B_p(x) = B_{p0} + \delta B_p\, x + \delta B_{pQ}\, x^2
$$

Phase (after flat-Earth removal):

$$
\phi_Z = \frac{4\pi}{\lambda}
\Bigl[\sqrt{R^2 - 2R(B_n\sin\theta_D + B_p\cos\theta_D) + B^2} - R\Bigr]
- \frac{4\pi}{\lambda}\Bigl[-B_n\sin\theta_\text{flat} - B_p\cos\theta_\text{flat} + \tfrac{B^2}{2R}\Bigr]
$$

Phase error (combines $\pi/4$ baseline noise and 15 m DEM uncertainty):

$$
\sigma_\phi = \sqrt{\left(\frac{\pi}{4}\right)^2 + \left(\frac{4\pi}{\lambda}\frac{|B_n| \cdot 30}{R\sin\theta_\text{flat}}\right)^2}
$$

---

### `getBaseline.c`

**`getBaseline`**(baselineFile, params, noPhase) → void  
Parses a CW-format baseline file and populates a `vhParams` structure with:
$B_n$, $B_p$, $\delta B_n$, $\delta B_p$, constant range bias, $\delta B_{nQ}$, $\delta B_{pQ}$,
and the 6×6 parameter covariance matrix. If `noPhase` is set, only baseline parameters
relevant to offset mosaicking are loaded.

---

### `svBase.c`

**`svBnBp`**(myTime, theta, dt1t2, sv1, sv2, \*bn, \*bp, lookDir) → double  
*Calls:* `svBaseTCN`, `polintVec`  
Computes interferometric baseline components $(B_n, B_p)$ in the normal/parallel (cross-track
and height) frame from the state vectors of two acquisitions:

$$
\mathbf{b} = \mathbf{r}_2(t) - \mathbf{r}_1(t)
$$

Decomposed into TCN (Track/Cross-track/Normal) components via `svBaseTCN`, then rotated
to $(B_n, B_p)$ using the look angle $\theta$:

$$
B_n = -B_C\cos\theta + B_N\sin\theta, \qquad B_p = B_C\sin\theta + B_N\cos\theta
$$

**`svBaseTCN`**(myTime, dt1t2, sv1, sv2, bTCN[3]) → void  
*Calls:* `polintVec`, `norm`, `cross`  
Computes the baseline vector in Track/Cross-track/Normal coordinates. The TCN frame
is defined as:
- $\hat{N} = -\mathbf{r}/|\mathbf{r}|$ (radial, outward)
- $\hat{C} = (\hat{N} \times \hat{V})/|\hat{N} \times \hat{V}|$ (cross-track)
- $\hat{T} = \hat{C} \times \hat{N}$ (along-track)

**`svInitBnBp`**(inputImage, offsets) → void  
Precomputes $B_n$ and $B_p$ as a function of azimuth index across the full image extent
and stores in `offsets->bnS`, `offsets->bpS`.

**`svInterpBnBp`**(inputImage, offsets, azimuth, \*bnS, \*bpS) → void  
Nearest-neighbour lookup into the precomputed $B_n(az)$, $B_p(az)$ tables.

**`svInitAzParams`**(inputImage, offsets) → void  
Fits a polynomial (constant + range + azimuth + range×azimuth + range²) to the
along-track azimuth offset between two images, derived from state vectors, over a
50-point grid. Used to compute and subtract the geometric zero-offset in azimuth
speckle tracking.

**`svAzOffset`**(inputImage, offsets, range, azimuth) → double  
Evaluates the polynomial azimuth offset model at a given range/azimuth pixel coordinate.

**`svOffsets`**(image1, image2, offsets, \*cnstR, \*cnstA) → void  
Computes bulk range and azimuth registration offsets between two images using
their control-point state vectors.

**`dmatrixRecycle`**(nrl, nrh, ncl, nch, mR, mBuf) → double\*\*  
Allocates a 2D double matrix view over pre-allocated row-pointer (`mR`) and data (`mBuf`)
buffers. Used to avoid repeated `malloc`/`free` in tight loops.

---

## Offset Interpolation

### `interpOffsets.c`

**`interpAzOffset`**(range, azimuth, offsets, inputImage, Range, theta, azSLPixSize) → float  
*Calls:* `bilinearInterp`, `svAzOffset`  
Interpolates the azimuth speckle-tracking offset map at (range, azimuth) and subtracts
the geometric zero-offset:

$$
\Delta a_\text{corr} = \Delta a_\text{raw} \cdot \delta s - \left[c_1 + R\sin\theta\,\frac{dB_c}{ds} - R\cos\theta\,\frac{dB_h}{ds} + \ell\cdot x\right] - \text{SV offset}
$$

where $\delta s$ is the single-look azimuth pixel size, $c_1$ is the constant offset,
and $\ell\cdot x$ is the optional linear along-track trend correction.

**`interpAzSigma`**(range, azimuth, offsets, inputImage, Range, theta, azSLPixSize) → float  
Interpolates the azimuth offset sigma map and combines with the streak noise floor:

$$
\sigma_a = \sqrt{\sigma_\text{map}^2 + \sigma_\text{streaks}^2} \cdot \delta s
$$

Sigma is capped at `MAXSIG` = 0.2 pixels to prevent runaway errors.

**`interpRangeOffset`**(range, azimuth, offsets, inputImage, Range, thetaD, rSLPixSize, theta, \*demError) → float  
*Calls:* `bilinearInterp`, `svInterpBnBp`  
Interpolates the range speckle-tracking offset map and subtracts the geometric (baseline)
zero-offset using the quadratic baseline model (Joughin et al., *J. Glaciol.*, 1996, Eq. 7):

$$
\Delta r_0 = \sqrt{R^2 - 2R(B_n\sin\theta_D + B_p\cos\theta_D) + B^2} - R + c
$$

Also returns the DEM-induced range error estimate:

$$
\sigma_\text{DEM} = \frac{|B_n| \cdot 15}{R\sin\theta}
$$

(15 m nominal DEM uncertainty.)

**`interpRangeSigma`**(range, azimuth, offsets, inputImage, Range, thetaD, rSLPixSize) → float  
Interpolates range offset sigma and applies a minimum floor (`offsets->sigmaRange`).
Sigma is capped at `MAXSIG` = 0.2 pixels.

---

## Vector Mathematics

### `vectorFunc.c`

**`dot`**(x1, y1, z1, x2, y2, z2) → double  
Dot product: $\mathbf{a} \cdot \mathbf{b} = a_x b_x + a_y b_y + a_z b_z$

**`norm`**(x1, y1, z1) → double  
Euclidean norm: $|\mathbf{a}| = \sqrt{a_x^2 + a_y^2 + a_z^2}$

**`cross`**(a1, a2, a3, b1, b2, b3, \*c1, \*c2, \*c3) → void  
Cross product: $\mathbf{c} = \mathbf{a} \times \mathbf{b}$

---

### `polintVec.c`

**`polintVec`**(xa[], y1[]…y6[], x, yr1…yr6) → void  
Neville's algorithm for polynomial interpolation of six dependent variables simultaneously.
Given tabulated values at nodes $x_a$, evaluates all six outputs at $x$ in one pass.
Used for efficient state-vector interpolation (position x, y, z and velocity vx, vy, vz).

---

### `bilinearInterp.c`

**`bilinearInterp`**(fimage, range, azimuth, nr, na, minvalue, noData) → float  
Bilinear interpolation on a 2D float array at fractional pixel coordinates (range, azimuth).
Returns `noData` if the integer-pixel neighbours are out-of-bounds or below `minvalue`:

$$
f(r, a) = (1-\delta r)(1-\delta a)\,f_{ij}
         + \delta r(1-\delta a)\,f_{i+1,j}
         + (1-\delta r)\delta a\,f_{i,j+1}
         + \delta r\,\delta a\,f_{i+1,j+1}
$$

where $\delta r = r - \lfloor r \rfloor$, $\delta a = a - \lfloor a \rfloor$.

---

## Coordinate Rotation

### `rotateFlowDirectionToRA.c`

**`rotateFlowDirectionToRA`**(dxn, dyn, \*dan, \*drn, xyAngle, hAngle) → void  
Rotates a displacement/velocity vector from the polar-stereographic (x, y) frame into the
SAR (range, azimuth) frame. The rotation angle is $\alpha = h_\text{angle} - \text{xyAngle}$:

$$
\begin{pmatrix} d_r \\ d_a \end{pmatrix}
= \begin{pmatrix} \cos\alpha & -\sin\alpha \\ \sin\alpha & \cos\alpha \end{pmatrix}
\begin{pmatrix} d_x \\ d_y \end{pmatrix}
$$

---

### `rotateFlowDirectionToXY.c`

**`rotateFlowDirectionToXY`**(drn, dan, \*dxn, \*dyn, xyAngle, hAngle) → void  
Inverse of `rotateFlowDirectionToRA`. Rotates from SAR (range, azimuth) back to
polar-stereographic (x, y).

---

## Mosaicking and Scaling

### `computeScale.c`

**`computeScale`**(inImage, scale, azimuthSize, rangeSize, fl, weight, minVal) → void  
*Calls:* `fillRadialKernel`  
Computes per-pixel feathering weights for blending overlapping SAR images.
For each valid pixel, the weight decays from `weight` to 0 over a distance `fl` (m)
from the image boundary. Interior pixels far from any edge retain the full weight.

**`fillRadialKernel`**(rDist, fl, weight) → void  
Pre-fills a radial distance lookup kernel used internally by `computeScale` to accelerate
the feathering computation.

---

### `scalingFunctions.c`

**`redoNormalization`**(myWeight, outputImage, iMin, iMax, jMin, jMax, vX, vY, vZ, eX, eY, sX, sY, sZ, fScale, vxTmp, vyTmp, vzTmp, sxTmp, syTmp, statsFlag) → void  
Accumulates a new velocity/error estimate into the running weighted sum, applying both the
per-pixel feather weight `fScale` and the per-image weight `myWeight`:

$$
v_x^\text{acc} \mathrel{+}= \hat{v}_x \cdot f \cdot w, \qquad
s_x^\text{acc} \mathrel{+}= \sigma_x^{-2} \cdot f \cdot w, \qquad
e_x^\text{acc} \mathrel{+}= \sigma_x^{-2} \cdot f^2 \cdot w^2
$$

In `statsFlag` mode, accumulates $v_x^2$ for variance computation instead.

**`undoNormalization`**(outputImage, vX, vY, vZ, eX, eY, sX, sY, sZ, fScale, statsFlag) → void  
Reverses `endScale` prior to ingesting a new mosaicking round, multiplying stored
normalised velocities back by their accumulated scale weights.

**`endScale`**(outputImage, vX, vY, vZ, eX, eY, sX, sY, sZ, statsFlag) → void  
Finalises the weighted average after all inputs have been accumulated:

$$
\bar{v}_x = \frac{v_x^\text{acc}}{s_x^\text{acc}}, \qquad
\sigma_x^2 = \frac{e_x^\text{acc}}{(s_x^\text{acc})^2}
$$

In `statsFlag` mode, computes the sample variance:
$\sigma_x^2 = \langle v_x^2 \rangle - \bar{v}_x^2$.

---

## Irregular Data

### `parseIrregFile.c`

**`parseIrregFile`**(irregFile, \*\*irregData) → void  
Parses a text file that lists paths to individual irregular scattered-data files and
builds a linked list of `irregularData` nodes.

---

### `getIrregData.c`

**`getIrregData`**(irregData) → void  
Reads each irregular data file listed in `irregData`, converts point (lat, lon) positions
to polar-stereographic XY coordinates, and performs Delaunay triangulation on the point
set. The triangulation is stored for later use by `addIrregData`.

---

### `addIrregData.c`

**`addIrregData`**(irregDat, outputImage, fl) → void  
*Calls:* `getIrregData`, `lltoxy1`  
Blends triangulated irregular velocity observations into the output grid. For each output
pixel, iterates over triangles to find containment (cross-product sign test), then
performs **linear barycentric interpolation** within the enclosing triangle:

$$
v_x = a_x\, x + b_x\, y + c_x, \qquad v_y = a_y\, x + b_y\, y + c_y
$$

Triangles wider than 15 km or larger than 75 km² are rejected. A fixed error of
$\sigma = 200$ m/yr is assigned to all irregular data to prevent them from overriding
higher-quality SAR estimates in the weighted combination.

---

## Input File Parsing

### `parseInputFile.c`

**`parseInputFile`**(inputFile, inputImage) → void  
Parses a SAR image `.geodat` parameter file. Reads near range, PRF, azimuth pixel size,
wavelength, look direction, heading, number of looks, and the full state-vector table.
Also supports GeoJSON input via `parseGeojson`.

**`parseGeojson`**(inputFile, inputImage) → void  
Parses a GeoJSON file containing SAR image metadata (stored as feature properties) and
state vectors (stored as arrays of position/velocity fields). Populates the same
`inputImageStructure` as `parseInputFile`.

**`dumpParsedInput`**(inputImage, fp) → void  
Writes a human-readable summary of parsed image parameters to `fp`, useful for validation.

---

### `getMVhInputFile.c`

**`getMVhInputFile`**(inputFile, phaseFiles, geodatFiles, baselineFiles, offsetFiles,
azParamsFiles, rOffsetFiles, rParamsFiles, outputImage, nDays, weights, crossFlags,
\*nFiles, offsetFlag, rOffsetFlag, threeDOffFlag) → void  
Parses the multi-product input list file for `mosaic3d`. Reads the output grid geometry
from line 1, file count from line 2, then allocates and fills arrays for all file types
(phase, geodat, baseline, azimuth offsets, azimuth params, range offsets, range params).
Handles the `rOffsetFlag` (`.vrt` extension triggers GDAL read path) and `crossFlag`
per image.

---

### `parseSLCvrt.c`

**`parseSLCVrt`**(vrtFile, sarD, sv, \*byteOrder) → void  
Parses SAR image metadata and state vectors embedded in a VRT XML file. Populates
a `SARData` and `stateV` structure.

**`parseSLCVrtNew`**(vrtFile, sarD, sv, \*byteOrder, \*hDS, \*hBand) → void  
Like `parseSLCVrt` but also returns open GDAL dataset and raster band handles for
subsequent pixel reads.

---

## Tie Points

### `readTiePoints.c`

**`readTiePoints`**(fp, tiePoints, noDEM) → void  
Reads tie-point records from an open file. Each record contains lat, lon, elevation,
and velocity (vx, vy, vz) components, plus optional per-point weights (e.g., rock vs.
ice). If `noDEM` is set, elevation is not required.

---

### `computeTiePoints.c`

**`computeTiePoints`**(inputImage, tiePoints, dem, noDEM, geodatFile, shelfMask, suppressOutput) → void  
*Calls:* `llToImageNew`, `getXYHeight`, `interpTideDiff`, `getShelfMask`  
Projects all tie points from geodetic coordinates to SAR image range/azimuth using
`llToImageNew`. Applies DEM height corrections and optionally tidal corrections on
floating ice (via `shelfMask`). Discards tie points that fall outside the image boundaries.

---

## Output

### `outputGeocodedImage.c`

**`outputGeocodedImage`**(outputImage, outputFile) → void  
Writes the geocoded backscatter (or other) image to disk as a flat binary big-endian float
file, followed by a `.geodat` text file containing grid geometry (origin in km, pixel
spacing in km, raster size).

**`outputGeocodedImageTiff`**(outputImage, outputFile, driverType, epsg, summaryMetaData, noDataValue, dataType) → void  
Writes the output image as a GeoTIFF or COG using GDAL. Sets the coordinate reference
system via EPSG code, attaches summary metadata as TIFF tags, and applies the no-data
value and data type.

---

## Offset Corrections (Ionospheric)

### `readOffsetCorrectionImpl.c` / `readOffsetCorrection.c`

**`readOffsetCorrection`**(correctionFile, offsets, bufferMode) → void  
Reads an ionospheric range-offset correction field from a VRT file into a shared memory
pool. The `bufferMode` selects between the image-1 and image-2 correction buffers when
two corrections are needed (e.g., for crossing-orbit pairs). NaN values are preserved.

**`mallocCorrectionBuffer`**(bufferMode) → void  
Lazily allocates the correction buffer pool on first use.

---

### `interpolateOffsetCorrection.c`

**`interpolateOffsetCorrection`**(corr, range, azimuth, minValue, noData) → float  
Bilinear interpolation of the ionospheric correction at SLC range/azimuth coordinates.
Returns `noData` if the correction is missing or below `minValue`.

---

## Date and Time

### `julianDay.c`

**`julday`**(mm, id, iyyy) → int32_t  
Converts a calendar date (month, day, year) to an integer Julian day number.

**`juldayDouble`**(mm, id, iyyy) → double  
Same as `julday` but returns a double with a $-0.5$ offset (Julian date epoch convention).

**`julian_to_gregorian`**(jd, \*year, \*month, \*day) → void  
Converts an integer Julian day number to Gregorian calendar date.

**`jd_to_date_and_time`**(jd, \*year, \*month, \*day, \*hour, \*minute, \*second) → void  
Converts a fractional Julian day to Gregorian date and time of day.

---

## Miscellaneous Utilities

### `initMatrix.c`

**`initFloatMatrix`**(x, nr, nc, initValue) → void  
**`initDoubleMatrix`**(x, nr, nc, initValue) → void  
Fill a pre-allocated 2D float (or double) array of size `nr` × `nc` with `initValue`.

---

### `getDataStringSpecial.c`

**`getDataStringSpecial`**(fp, lineCount, line, \*eod, special, \*specialFound) → int32_t  
Line-by-line reader for structured ASCII data files. Skips blank lines and comment lines
(`;` prefix). Recognises an end-of-data marker (`&`) and an optional user-defined
`special` character. Returns the count of data lines read.

---

### `getRegion.c`

**`getRegion`**(image, \*iMin, \*iMax, \*jMin, \*jMax, outputImage) → void  
*Calls:* `llToImageNew`, `xytoll1`  
Determines the output grid pixel bounding box that intersects the footprint of `image`.
Projects the four SAR image corners (near/far range × start/end azimuth) to the output
grid and clips to valid bounds. Used to restrict inner loops to the relevant output region.

---

### `getAzimuthBoundsForXYBox.c`

**`getAzimuthBoundsForXYBox`**(imin, imax, jmin, jmax, currentImage, outputImage, \*azimuthMin, \*azimuthMax) → float  
*Calls:* `xytoll1`, `llToImageNew`  
Finds the range of SAR azimuth pixel values that correspond to an XY output grid
sub-region. Samples corner and edge points of the box, projects to SAR image coordinates,
and returns the min/max azimuth. Used to limit pre-loading of offset data to the
relevant azimuth extent.

---

### `readOffsets.c`

**`readOffsetDataAndParams`**(offsets, azimuthMin, azimuthMax) → void  
Reads offset data arrays (range offset, azimuth offset, and their sigma maps) and the
associated parameter files (`.dat` metadata), loading only the azimuth rows between
`azimuthMin` and `azimuthMax` to reduce memory. Handles both flat binary and GDAL
(VRT/GeoTIFF) input formats.

---

### `readOldPar.c`

**`readOldPar`**(parFile, sarD, stateV) → void  
Parses a legacy SAR parameter file format and populates `SARData` and state-vector
structures. Primarily used to support older ERS/RADARSAT parameter files predating
the `.geodat` format.

---

### `geojsonCode.c`

**`getGeojsonDataSet`**(geojsonFile) → OGRDataSourceH  
Creates an OGR GeoJSON data source for writing vector feature data.

**`createGeometry`**(lat, lon) → OGRGeometryH  
Builds a polygon OGR geometry from arrays of lat/lon boundary points.

**`svTag`**(i, svType) → const char\*  
Generates a field-name string (e.g., `"sv_x_003"`) for encoding state-vector
position or velocity components as GeoJSON properties.

**`createFeatureDef`**(nState) → OGRFeatureDefnH  
Creates an OGR feature definition schema for a SAR image metadata record, including
fields for orbit parameters, timing, and `nState` state-vector position/velocity entries.

---

### `buffers.c`

Declares external buffer variables shared across the mosaicking pipeline. No public functions.

---

## Key Data Structures (Summary)

| Structure | Contents |
|-----------|----------|
| `inputImageStructure` | SAR image geometry: near range, PRF, pixel sizes, wavelength, look direction, state vectors, antenna pattern |
| `outputImageStructure` | Output grid geometry: origin (km), pixel size (m), raster size, projection |
| `xyDEM` | XY-projected DEM: origin, pixel spacing, size, standard latitude, rotation, elevation array |
| `vhParams` | Interferometric parameters: $B_n$, $B_p$, $\delta B_n$, $\delta B_p$, $\delta B_{nQ}$, $\delta B_{pQ}$, constant bias, covariance matrix |
| `Offsets` | Offset image geometry and baseline parameters: $B_n$, $B_p$, $c$, $dB_c/ds$, $dB_h/ds$, range/azimuth offset arrays and sigmas |
| `tiePointsStructure` | Tiepoint arrays: lat, lon, height, vx, vy, vz, per-point weights |
| `irregularData` | Triangulated scattered velocity data: XY positions, velocities, Delaunay triangulation |
| `stateV` | Satellite state vectors: times, ECEF position (x, y, z) and velocity (vx, vy, vz) arrays |
