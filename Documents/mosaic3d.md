# mosaic3d — 3D Velocity Mosaic Builder

## Purpose

Builds a mosaicked surface velocity field (vx, vy, vz) by combining measurements from
multiple input data types — InSAR phase, speckle-tracked range/azimuth offsets,
crossing-orbit range offsets, and Landsat feature tracking. Each data type contributes
an independent estimate at each output pixel; the final result is the error-weighted
combination of all available inputs, producing an optimal velocity mosaic with
associated uncertainty estimates.

---

## Usage

```
mosaic3d [options] inputFile demFile outFileBase
```

### Required Arguments

| Argument      | Description |
|---------------|-------------|
| `inputFile`   | ASCII input list file describing all SAR/Landsat inputs (geodats, phase files, baseline files, offset files, etc.) |
| `demFile`     | DEM file in XY (polar stereographic) format — sets projection and provides elevations |
| `outFileBase` | Base name for output files |

### Options

> **Note:** Some flags are experimental or may be obsolete. Options marked ⚠️ below
> are known to be experimental, rarely used, or may not be fully supported in current
> builds.

| Option                          | Description |
|---------------------------------|-------------|
| `-fl <length>`                  | Feathering length (m) at image edges for blending overlapping inputs |
| `-north` / `-south`             | Force northern / southern hemisphere projection |
| `-shelfMask <file>`             | Ice-shelf and grounding-zone mask for tidal corrections |
| `-tideFile <file>`              | Tidal model file for shelf corrections |
| `-verticalCorrection <file>`    | Vertical correction field (e.g. submergence/emergence) |
| `-irregFile <file>`             | Irregularly spaced supplemental velocity data |
| `-landSat <file>`               | Landsat feature-tracking input list |
| `-extraTieFile <file>`          | ⚠️ Extra tiepoints for output tiepoint file (used with `-makeTies`) |
| `-tieThresh <val>`              | ⚠️ Threshold for tiepoint selection (used with `-makeTies`) |
| `-refVel <file>`                | ⚠️ Reference velocity map for clipping outliers (experimental) |
| `-date1 <date>` / `-date2`      | Date range for output (YYYY-MM-DD) |
| `-timeThresh <days>`            | Maximum time separation for crossing-orbit offset pairs (`-3dOff`) |
| `-timeThreshPhase <days>`       | Maximum time separation for crossing-orbit InSAR pairs |
| `-no3d`                         | Skip the 3D InSAR phase solution |
| `-noVh`                         | Skip the InSAR phase + azimuth offset (Vh) solution |
| `-3dOff`                        | Include the crossing-orbit range offset 3D solution |
| `-makeTies`                     | ⚠️ Output a tiepoint file rather than a velocity mosaic (experimental) |
| `-statsFlag`                    | ⚠️ Use unit weights for computing spatial statistics (experimental) |
| `-GTiff`                        | Write output as GeoTIFF |
| `-COG`                          | Write output as Cloud-Optimised GeoTIFF |
| `-quiet`                        | Suppress verbose output |

### Output Files

**Flat binary output** (default):

| File                        | Contents |
|-----------------------------|----------|
| `outFileBase.vx`            | East velocity component (m/yr) — 32-bit MSB float |
| `outFileBase.vy`            | North velocity component (m/yr) — 32-bit MSB float |
| `outFileBase.vz`            | Vertical velocity component (m/yr) — 32-bit MSB float |
| `outFileBase.ex`            | East velocity error (m/yr) — 32-bit MSB float |
| `outFileBase.ey`            | North velocity error (m/yr) — 32-bit MSB float |
| `outFileBase.vx.geodat`     | Grid geometry descriptor (one per component, see below) |
| `outFileBase.meta`          | Input file log and metadata |

Each binary file is accompanied by a `.geodat` text file with the following format:

```
# 2
;
;  Image size (pixels) nx ny
;
nx  ny
;
;  Pixel size (m) deltaX deltaY
;
deltaX  deltaY
;
;  Origin, lower left corner (km) Xo  Yo
;
Xo  Yo
&
```

**GeoTIFF / COG output** (`-GTiff` or `-COG`):

| File                        | Contents |
|-----------------------------|----------|
| `outFileBase.vx.tif`        | East velocity component (m/yr) |
| `outFileBase.vy.tif`        | North velocity component (m/yr) |
| `outFileBase.vz.tif`        | Vertical velocity component (m/yr) |
| `outFileBase.ex.tif`        | East velocity error (m/yr) |
| `outFileBase.ey.tif`        | North velocity error (m/yr) |

---

## Input File Format

The input file is an ASCII file parsed by `getMVhInputFile`. Comments begin with `;`.

### Header (2 lines)

**Line 1 — Output grid geometry:**
```
x0  y0  xSize  ySize  deltaX  deltaY
```
| Field    | Description |
|----------|-------------|
| `x0 y0`  | Origin of output grid in km (polar stereographic) |
| `xSize ySize` | Grid extent in km |
| `deltaX deltaY` | Pixel spacing in km |

**Line 2 — Number of input images:**
```
nFiles
```

### Per-Image Lines

One line per image. Fields depend on which data types are active:

**Phase only** (InSAR, no offsets):
```
phaseFile  geodatFile  baselineFile  nDays  [weight]
```

**Phase + azimuth offsets** (`-noVh` mode or Vh mode):
```
phaseFile  geodatFile  baselineFile  nDays  [weight]  offsetFile  azParamsFile
```

**Full speckle tracking** (range + azimuth offsets, `-3dOff` or `rOffsetFlag`):
```
phaseFile  geodatFile  baselineFile  nDays  [weight]  offsetFile  azParamsFile  rOffsetFile  rParamsFile  [crossFlag]
```

| Field          | Description |
|----------------|-------------|
| `phaseFile`    | Geocoded interferometric phase or `nophase` to skip phase for this image |
| `geodatFile`   | SAR image geometry parameter file |
| `baselineFile` | CW state-vector baseline file (Bn, Bp, dBn, dBp, const, dBnQ, dBpQ) |
| `nDays`        | Temporal baseline in days |
| `weight`       | Image weight (default 1.0); used for time-overlap weighting |
| `offsetFile`   | Azimuth offset file from speckle tracking |
| `azParamsFile` | Azimuth parameter file from `azparams` |
| `rOffsetFile`  | Range offset file from speckle tracking |
| `rParamsFile`  | Range parameter file from `rparams` |
| `crossFlag`    | 1 = eligible for crossing-orbit pairs, 0 = same-orbit only (default 1) |

Use `none` or `None` for any optional file field to indicate it is not available.

> **Offset file formats:** `offsetFile` and `rOffsetFile` can be flat binary MSB
> (big-endian) float files, which is the default format. If the filename ends in
> `.vrt`, the file is read via GDAL using the accompanying VRT descriptor, allowing
> any GDAL-supported format (GeoTIFF, NetCDF, etc.) to be used as input.

---

## Algorithm

### Overview — Processing Pipeline

mosaic3d runs up to five successive mosaicking steps. Each step accumulates its result
into a shared output grid using error-weighted averaging:

| Step | Routine               | Data type used | Flag required |
|------|-----------------------|----------------|---------------|
| 0    | `makeLandSatMosaic`   | Landsat optical feature tracking | `-landSat` |
| 1    | `make3DMosaic`        | Crossing-orbit InSAR phase pairs | default (disable with `-no3d`) |
| 2    | `make3DOffsets`       | Crossing-orbit speckle-tracked range offset pairs | `-3dOff` |
| 3    | `makeVhMosaic`        | Single-pass InSAR phase + azimuth offsets | default (disable with `-noVh`) |
| 4    | `speckleTrackMosaic`  | Single-pass range + azimuth speckle-tracked offsets | `rOffsetFlag` in input file |

After all steps, irregularly gridded supplemental data (`-irregFile`) are blended in.

---

### Step 0 — Landsat Feature Tracking (`makeLandSatMosaic`)

Optical feature-tracking displacements from Landsat image pairs are processed first so
that subsequent SAR steps can build on them. For each Landsat image in the input list:

1. The feature-tracking match result grid (x, y displacements in pixels) is read and
   bilinearly interpolated to each output pixel within the image's bounding box.
2. A polynomial correction (up to linear in x and y) determined from tiepoints is
   subtracted to remove residual registration offsets:

$$
\Delta x_\text{corr} = \Delta x - \left(p_{x0} + p_{x1}(x - x_0) + p_{x2}(y - y_0)\right)
$$

3. Displacements are scaled to velocity (m/yr) using the image pair temporal baseline
   $\Delta t$ and a latitude-dependent map-projection scale factor $\lambda_\text{lat}$:

$$
v_x = \frac{365.25}{\Delta t}\,\Delta x_\text{corr} \cdot \delta x_\text{pix} \cdot \lambda_\text{lat}
$$

4. The per-pixel error variance combines the global tiepoint residual (with the
   mean per-image tracking sigma removed to avoid double-counting), local per-pixel
   tracking sigma, and a quantisation floor $\sigma_Q = 0.05$ px:

$$
\sigma_x^2 = \max\!\left(\sigma_\text{tie}^2 - \bar{\sigma}_\text{img}^2,\, 0\right) + \sigma_\text{local}^2 + \sigma_Q^2
$$

   with a cap at $\sigma_\text{max} = 0.2$ px. The weight $w = 1/\sigma_v^2$ is used
   in the error-weighted accumulation.

---

### Step 1 — Crossing InSAR Phase Pairs (`make3DMosaic`)

All image pairs with sufficient heading difference (typically ascending × descending)
are looped over. For each output pixel the program:

1. Geocodes the pixel to lat/lon using the DEM, then projects to range/azimuth in both
   images.
2. Interpolates the phase value from each image.
3. Subtracts the topographic phase contribution $\phi_Z$ (flat-earth removed baseline
   phase) computed from the polynomial baseline model:

$$
B_n(x) = B_{n0} + \delta B_n\, x + \delta B_{nQ}\, x^2, \quad
B_p(x) = B_{p0} + \delta B_p\, x + \delta B_{pQ}\, x^2
$$

$$
\phi_Z = \frac{4\pi}{\lambda}\left(\sqrt{R^2 - 2R(B_n\sin\theta_D + B_p\cos\theta_D) + B^2} - R\right) - \phi_\text{flat}
$$

4. Applies tidal and submergence/emergence corrections on floating ice.
5. Constructs the 2×2 geometric conversion matrix $\mathbf{A}$ from the two look
   directions, and the surface-slope correction matrix $\mathbf{B}$ from the DEM.
6. Solves for $(v_x, v_y)$ and derives $v_z = v_x \partial z/\partial x + v_y \partial z/\partial y$.
7. Propagates baseline covariance to a per-pixel phase error $\sigma_\phi$.

#### Matrix A — Geometric Conversion (`computeA`)

Define the following angles:

| Symbol | Meaning |
|--------|---------|
| $H_A$, $H_D$ | Satellite heading angles (radians from north, CW) for ascending and descending images |
| $\alpha = H_A - H_D$ | Heading difference between the two images |
| $\phi = \text{atan2}(-y, -x)$ | Azimuth angle of the output pixel in polar-stereographic coordinates |
| $\beta = \phi - H_A$ | Pixel azimuth angle relative to the ascending heading |

The A matrix maps scaled phase measurements to horizontal velocity components:

$$
\mathbf{A} = \frac{1}{\sin^2\!\alpha}
\begin{pmatrix}
\cos\beta - \cos\alpha\cos(\alpha+\beta) & \cos(\alpha+\beta) - \cos\alpha\cos\beta \\
\sin\beta - \cos\alpha\sin(\alpha+\beta) & \sin(\alpha+\beta) - \cos\alpha\sin\beta
\end{pmatrix}
$$

A minimum heading difference of $|\alpha| \geq 0.8$ rad ($\approx 46°$) is required for a
well-conditioned solution; pixels where $|\alpha| < 0.8$ are skipped.

#### Matrix B — Surface-Slope Correction (`computeB`)

Surface slopes $\partial z/\partial x$ and $\partial z/\partial y$ are computed from the DEM
by centred finite differences over a spacing of at least 90 m, and capped at $\pm 0.1$
($\approx 5.7°$). The B matrix accounts for the vertical velocity component
$v_z = v_x\,\partial z/\partial x + v_y\,\partial z/\partial y$ contributing to the
line-of-sight phase through the $\cos\psi / \sin\psi$ projection:

$$
\mathbf{B} =
\begin{pmatrix}
\dfrac{\partial z/\partial x}{\tan\psi_A} & \dfrac{\partial z/\partial y}{\tan\psi_A} \\[8pt]
\dfrac{\partial z/\partial x}{\tan\psi_D} & \dfrac{\partial z/\partial y}{\tan\psi_D}
\end{pmatrix}
$$

where $\psi_A$, $\psi_D$ are the local incidence angles for the ascending and descending images.

#### Full Inversion (`computeVxy`)

Each phase is scaled to velocity units (m/yr):

$$
p_i = \frac{365.25}{\frac{4\pi}{\lambda_i}\,N_{\text{days},i}\,\sin\psi_i}\,\phi_i
$$

The forward model including the slope correction is:

$$
\begin{pmatrix} p_A \\ p_D \end{pmatrix} =
\left(\mathbf{A}^{-1} - \mathbf{B}\right)
\begin{pmatrix} v_x \\ v_y \end{pmatrix}
$$

which rearranges to the solution:

$$
\begin{pmatrix} v_x \\ v_y \end{pmatrix} =
\left(\mathbf{I} - \mathbf{A}\mathbf{B}\right)^{-1} \mathbf{A}
\begin{pmatrix} p_A \\ p_D \end{pmatrix}
$$

If $\det(\mathbf{I} - \mathbf{AB}) < 0.25$ (poorly conditioned due to extreme slopes), no
solution is assigned.

Error variances are propagated as the diagonal of the output covariance matrix:

$$
\sigma_{v_x}^2 = D_{00}^2\,\sigma_{p_A}^2 + D_{01}^2\,\sigma_{p_D}^2, \qquad
\sigma_{v_y}^2 = D_{10}^2\,\sigma_{p_A}^2 + D_{11}^2\,\sigma_{p_D}^2
$$

where $\mathbf{D} = (\mathbf{I} - \mathbf{AB})^{-1}\mathbf{A}$.

#### Phase Error (`computePhiZM3d`)

The per-pixel phase error combines baseline parameter uncertainty (propagated through the
6×6 covariance matrix $\mathbf{C}$) and tiepoint noise $\sigma_\text{tp}$:

$$
\sigma_\phi = \sqrt{\mathbf{v}^T \mathbf{C}\, \mathbf{v} + \min(\pi,\,\sigma_\text{tp})^2}
$$

where $\mathbf{v} = \frac{4\pi}{\lambda}(-\sin\theta_D,\,-\cos\theta_D,\,-x\sin\theta_D,\,-x\cos\theta_D,\,-x^2\sin\theta_D,\,-x^2\cos\theta_D)$
is the Jacobian of the topographic phase with respect to the six baseline parameters
$(B_n,\,B_p,\,\delta B_n,\,\delta B_p,\,\delta B_{nQ},\,\delta B_{pQ})$.

---

### Step 2 — Crossing-Orbit Range Offsets (`make3DOffsets`)

For image pairs flagged as crossing orbits (`crossFlag`) whose acquisition dates are
within `-timeThresh` days of each other, the program solves for $(v_x, v_y)$ using the
speckle-tracked range offsets from two different look directions — the same geometry as
Step 1 but without phase. For each pixel in the intersection region:

1. Geocodes the pixel and projects to range/azimuth in both images.
2. Bilinearly interpolates the range offset $\Delta r$ (m) from each offset field.
   Ionospheric range corrections are subtracted if provided.
3. Applies tidal and submergence corrections on floating ice (scaled by $\cos\psi$).
4. Scales offsets to horizontal velocity units:

$$
p_A = \frac{365.25\, \Delta r_A}{\Delta t_A \sin\psi_A}, \qquad
p_D = \frac{365.25\, \Delta r_D}{\Delta t_D \sin\psi_D}
$$

5. Constructs $\mathbf{A}$ and $\mathbf{B}$ matrices (same as Step 1) and calls
   `computeVxy` to solve for $(v_x, v_y)$.
6. Per-pixel errors combine the interpolated range offset sigma, DEM-induced range
   error, and baseline parameter uncertainty:

$$
\sigma_R = \sqrt{\sigma_\text{off}^2 + \sigma_\text{dem}^2 + \sigma_\text{base}^2}
$$

---

### Step 3 — InSAR Phase + Azimuth Offsets (`makeVhMosaic`)

For each single-pass image, the InSAR phase (range component) is combined with the
azimuth speckle-tracking offset (along-track component). The range velocity is derived
from the topography-corrected phase, and the azimuth offset provides the along-track
velocity. The two components are rotated from the radar frame into the map
(polar-stereographic x, y) frame using the local satellite heading angle.

---

### Step 4 — Pure Speckle Tracking (`speckleTrackMosaic`)

For images with both range and azimuth offset fields (but no phase), the program
computes velocity from the two speckle-tracked offset components. The range-direction
velocity accounts for surface slope and the azimuth-range coupling:

$$
v_r = \frac{d_r \cdot \frac{365.25}{\Delta t} / \sin\psi + v_a \cot\psi \cdot \partial z/\partial a}{1 - \cot\psi \cdot \partial z/\partial r}
$$

where $d_r$ is the range offset (m), $d_a$ the azimuth offset (m), $\psi$ the incidence
angle, and $\partial z/\partial r$, $\partial z/\partial a$ the surface slopes in range
and azimuth. On floating ice the slope correction is suppressed. Ionospheric range
corrections are subtracted if provided.

---

### Step 5 — Irregularly Gridded Supplemental Data (`addIrregData`)

After all raster-based steps, optionally spaced velocity observations (e.g. GPS,
stake measurements) can be blended in. The irregular data file (`-irregFile`) contains
a list of per-dataset files; each dataset holds point observations at arbitrary $(x, y)$
locations with associated $(v_x, v_y)$ values.

1. A Delaunay triangulation is pre-computed for each dataset (`getIrregData`).
2. For each output pixel, the bounding box of each triangle is computed and the pixel's
   $(x, y)$ position is tested for membership using the cross-product sign test.
3. If the pixel lies inside a valid triangle (max edge length ≤ 15 km, max area ≤ 75 km²),
   velocity is estimated by **linear barycentric interpolation** over the triangle plane:

$$
v_x = a_x\, x + b_x\, y + c_x, \qquad v_y = a_y\, x + b_y\, y + c_y
$$

   where the plane coefficients are solved from the three triangle vertex values.

4. A fixed error of $\sigma = 200$ m/yr ($w = 1/\sigma^2$) is assigned to all irregular
   data points. This large uncertainty ensures the irregular data fills gaps but does
   not override higher-quality SAR or Landsat estimates in the weighted combination.

---

### Error-Weighted Combination

After each step, the new result is accumulated into the output mosaic via
`redoNormalization`. For each pixel, the running weighted sum is updated:

$$
\bar{v}_x = \frac{\sum_i w_i^x \cdot \hat{v}_x^i}{\sum_i w_i^x}, \qquad
w_i^x = \frac{1}{\sigma_{x,i}^2}
$$

where $\sigma_{x,i}$ is the per-pixel velocity error from step $i$. The final
`endScale` pass converts accumulated weighted sums to normalised velocities and
error estimates. Edge blending between overlapping images uses a distance-weighted
feather zone of length `fl`.

---

## Dependencies

- `setup3D` — parses all input geodat, phase, baseline, and offset files
- `make3DMosaic` — crossing InSAR phase solution
- `make3DOffsets` — crossing-orbit range offset solution
- `makeVhMosaic` — single-pass InSAR phase + azimuth offset solution
- `speckleTrackMosaic` — pure speckle-tracking solution
- `makeLandSatMosaic` — Landsat feature-tracking solution
- `addIrregData` / `parseIrregFile` — irregularly spaced supplemental data
- `computeA`, `computeB`, `computeVxy` — 3D inversion geometry
- `computePhiZM3d` — topographic phase and baseline error propagation
- `readXYDEM` — DEM I/O and projection setup
- `svdfit` / `svdvar` — SVD least squares (Numerical Recipes)
- GDAL — GeoTIFF / COG output
