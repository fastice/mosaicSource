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

| Option                                 | Description |
|----------------------------------------|-------------|
| `-fl <length>`                         | Feathering length (m) at image edges for blending overlapping inputs |
| `-north`                               | Force northern hemisphere projection |
| `-shelfMask <file>`                    | Ice-shelf and grounding-zone mask for tidal corrections |
| `-noTide`                              | Skip tidal correction even on shelf pixels |
| `-tideFile <file>`                     | Tidal model file for shelf corrections |
| `-verticalCorrection <file>`           | Vertical correction field (e.g. submergence/emergence velocity) in m/yr |
| `-verticalCorrectionSuffix <suffix>`   | Suffix appended to baseline/rParams filenames only (not phase) for an alternate vertical-correction scenario |
| `-irreg <file>`                        | File listing irregularly spaced supplemental velocity datasets |
| `-landSat <file>`                      | File containing list of Landsat feature-tracking offset inputs |
| `-rOffsets`                            | Use range+azimuth speckle-tracked offsets for both velocity components where needed |
| `-offsets`                             | *(obsolete — azimuth offsets are always used; flag is silently ignored)* |
| `-no3d`                                | Skip the crossing-orbit InSAR phase (3D) solution |
| `-noVh`                                | Skip the single-pass InSAR phase + azimuth offset (Vh) solution |
| `-3dOff`                               | Enable crossing-orbit range offset 3D solution |
| `-noSepAscDesc`                        | Use heading difference rather than ascending/descending classification for crossing pairs |
| `-timeThresh <days>`                   | Max time separation for crossing-orbit range-offset pairs (`-3dOff`); default 12 days |
| `-timePhaseThresh <days>`              | Max time separation for crossing-orbit InSAR pairs; default 548 days (1.5 years) |
| `-date1 <MM-DD-YYYY>`                  | Start of output date range |
| `-date2 <MM-DD-YYYY>`                  | End of output date range |
| `-timeOverlap`                         | Include images that partially overlap the date range (weighted by fractional overlap); default is fully-contained-only |
| `-SVConst`                             | Apply constant-only state-vector correction on top of orbital solution |
| `-SVAlongTrack`                        | Apply along-track quadratic state-vector correction |
| `-refVel <file>`                       | Reference velocity map used to clip large residuals |
| `-initMap`                             | Interpolate `-refVel` map as the starting point for the mosaic |
| `-clipThresh <val>`                    | With `-refVel`: clip differences exceeding `val` m/yr for slow regions (< 100 m/yr) |
| `-extraTies <file>`                    | ⚠️ Extra tiepoints to include (used with `-makeTies`) |
| `-tieThresh <val>`                     | ⚠️ Velocity threshold for tiepoint selection (default 100 m/yr) |
| `-makeTies`                            | ⚠️ Output a tiepoint file rather than a velocity mosaic |
| `-sigmaAThresh <val>`                  | Skip azimuth offsets for an image if the residual from the azimuth parameter fit exceeds `val` metres (default 1000 m — effectively no filtering) |
| `-stats`                               | ⚠️ Compute unweighted mean vx/vy and std dev of ex/ey; output point count in vz channel (speckle-track inputs only) |
| `-vzFlag <val>`                        | Select vertical-channel output: 0 = vz (default), 1 = horizontal 1/sin ψ, 2 = 1/cos ψ, 3 = LOS scaled to m/yr, 4 = incidence angle |
| `-writeBlank`                          | Force output to be written even when no valid data exist |
| `-GTiff`                               | Write output as GeoTIFF |
| `-COG`                                 | Write output as Cloud-Optimised GeoTIFF |
| `-ompThreads <N>`                      | Number of OpenMP threads for parallel pixel processing (default: 4; overridden by `OMP_NUM_THREADS` environment variable) |
| `-center`                              | *(obsolete — silently ignored)* |
| `-useSquint`                           | Apply per-image squint(r,a) heading correction before building the crossing-orbit solving matrix; phase (`make3DMosaic`) only, default off — see "Squint (Residual Doppler) Correction" below |

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

4. Applies tidal (`-tideFile`, floating ice only) and submergence/emergence
   (`-verticalCorrection`) corrections to each phase:

$$
\phi_A \mathrel{+}= v_z^{\text{SMB}}\,\cos\psi_A\,\frac{4\pi}{\lambda_A}\,\frac{N_{\text{days},A}}{365.25},
\qquad
\phi_D \mathrel{+}= v_z^{\text{SMB}}\,\cos\psi_D\,\frac{4\pi}{\lambda_D}\,\frac{N_{\text{days},D}}{365.25}
$$

   where $v_z^{\text{SMB}}$ is the vertical rate (m/yr, **positive = up**) read directly,
   unmodified, from the `-verticalCorrection` grid (`interpVCorrect`, a pure bilinear
   interpolation with no sign flip) or the tide-height-rate grid (`-tideFile`), and
   $\psi_A,\psi_D$ are the local incidence angles. **Do not confuse with the solved output**
   $v_z$ **in step 6 below** — that's flow-driven vertical motion derived from slope×horizontal
   velocity; $v_z^{\text{SMB}}$ here is an independent, externally supplied vertical rate being
   *removed* from the phase before solving for horizontal motion. See "Vertical-Motion
   Correction Sign Convention" below for why `+=` (not `-=`) is correct given the up-positive
   convention, and "Look-Direction Sign Convention" for why no left/right-looking adjustment
   is needed.
5. Constructs the 2×2 geometric conversion matrix $\mathbf{A}$ from the two look
   directions (optionally squint-corrected, see "Squint (Residual Doppler) Correction"
   below — off by default), and the surface-slope correction matrix $\mathbf{B}$ from the DEM.
6. Solves for $(v_x, v_y)$ and derives $v_z = v_x \partial z/\partial x + v_y \partial z/\partial y$.
7. Propagates baseline covariance to a per-pixel phase error $\sigma_\phi$.

#### Heading Angle Convention (`computeHeading`)

$H_A$, $H_D$ (used below) are the **cross-track** heading at the pixel, not the along-track
flight direction. For a small latitude offset $\delta\text{lat}$, project $(\text{lat}\mp
\delta\text{lat},\,\text{lon})$ into the image's range/azimuth coordinates to get an azimuth
displacement $da$ and ground-range displacement $dgr$, then:

$$
H = \begin{cases}
\text{atan2}(da,\ dgr) & \text{right-looking} \\
\text{atan2}(da,\ -dgr) & \text{left-looking}
\end{cases}
$$

This sign flip (`common/computeHeading.c`) is the **only** place look direction enters the
velocity-inversion geometry. It is folded into $H_A$/$H_D$ before $\alpha$/$\beta$ — and hence
$\mathbf{A}$ — are computed, so no further look-direction handling is needed downstream of this
point.

#### Squint (Residual Doppler) Correction (`-useSquint`, off by default)

`computeHeading` returns the heading of the idealized, **zero-Doppler (broadside)** cross-track
direction — real NISAR acquisitions carry a small residual squint (~1.5°–1.7°) that this
geometry doesn't model. With `-useSquint`, each image's own measured squint, a 6-parameter
polynomial in range $r$ and azimuth $a$ fit upstream (`nisarhdf`/`SetupNISAR`, see their
CLAUDE.md) and threaded into the geodat, is evaluated and added directly to that image's own
heading **before** $\alpha$/$\beta$ — and hence $\mathbf{A}$ — are computed:

$$
H_A \to H_A + \text{squint}_A(r_A, a_A), \qquad H_D \to H_D + \text{squint}_D(r_D, a_D)
$$

This is exact (not a post-hoc rotation of the output $(v_x,v_y)$) because $\mathbf{A}$ is a pure
function of $\alpha,\beta$ with no other squint dependence — correcting the headings and
changing nothing else reproduces the matrix that would have been built from the true geometry.
**Phase only** (`make3DMosaic`'s `computeA` call) — `make3DOffsets`'s crossing-orbit range-offset
solution never applies this, flag or no flag, since offsets are self-consistent regardless of
squint by construction (the zero-Doppler condition forces true LOS ⊥ true velocity at the
assigned time, independent of squint). Ships off by default: the underlying analysis found no
NISAR data processed today needs it; verified (numerically and on a real overlapping
ascending/descending pair) to shift the recovered direction by the predicted ~1.5°–1.7° with the
correct sign — see `mosaicSource/CLAUDE.md`'s squint section for the full derivation and
verification.

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

#### Vertical-Motion Correction Sign Convention

$v_z^{\text{SMB}} > 0$ means the surface is moving **up** — toward the satellite, which
**decreases** slant range. Over time interval $\Delta t$ (years), the range change contributed
by vertical motion alone is therefore:

$$
\Delta R_{\text{vert}} = -\,v_z^{\text{SMB}}\,\cos\psi\,\Delta t
$$

(negative because moving toward the satellite shortens the range). This contaminates the raw
measurement:

$$
\Delta R_{\text{measured}} = \Delta R_{\text{horizontal}} + \Delta R_{\text{vert}}
= \Delta R_{\text{horizontal}} - v_z^{\text{SMB}}\,\cos\psi\,\Delta t
$$

so isolating the horizontal-motion-only range change requires **adding back**
$v_z^{\text{SMB}}\cos\psi\,\Delta t$:

$$
\Delta R_{\text{horizontal}} = \Delta R_{\text{measured}} + v_z^{\text{SMB}}\,\cos\psi\,\Delta t
$$

This is exactly the code's `+=` (the `X -= -Y` idiom seen in the source is algebraically `X += Y`)
— confirmed correct, not a sign bug, given $v_z^{\text{SMB}}$ is read unmodified from the grid
(positive = up). The phase form (Step 1) carries the identical sign: this codebase's own
topographic-phase formula (`computePhiZM3d`, $\phi_Z$ above) increases $\phi$ with increasing
path length, i.e. $\phi = +\frac{4\pi}{\lambda}\Delta R$, so the same $\Delta R \to
\Delta R + v_z^{\text{SMB}}\cos\psi\,\Delta t$ correction carries straight through with the
$4\pi/\lambda$ factor multiplied in, matching the `aPhase`/`dPhase` formula above. Step 2
(`make3DOffsets`, range offsets in metres) needs no phase-conversion factor and applies the
identical $\Delta r \mathrel{+}= v_z^{\text{SMB}}\cos\psi\,\Delta t$ form directly. The existing
`-tideFile` correction (tide height, also positive = up, by the standard oceanographic
convention) uses this same `+=` form, which `dzdtSubmergence`/`tideCorrection` share verbatim
in the source (`Mosaic3d/make3DMosaic.c`, `Mosaic3d/make3DOffsets.c`).

#### Look-Direction Sign Convention

The tide/submergence-emergence correction (above) needs **no** adjustment for left- vs
right-looking sensors. Verified directly against the geometry code:

- $\psi$ (incidence angle) comes from `psiRReZReH`/`thetaRReZReH` (`common/initRoutines.c`) —
  pure range/Earth-radius/satellite-height geometry (`acos`/`asin` of always-positive ratios),
  with no `lookDir` term. Incidence angle is the same physical quantity regardless of which side
  of the track the radar looks, so $\cos\psi$ is always positive and look-direction-independent.
- $4\pi/\lambda$ is a function of wavelength only — no look-direction term.
- The one and only look-direction-dependent sign in this whole inversion is the
  $\text{atan2}(da,\pm dgr)$ flip inside `computeHeading` (see above), which is already baked
  into $H_A$/$H_D$ — and hence $\mathbf{A}$ — before the correction terms are ever applied to
  $\phi_A$/$\phi_D$.

In other words, look direction changes how a given LOS phase gets decomposed into
$(v_x, v_y)$ (via $\mathbf{A}$), not the sign or magnitude of the vertical-motion correction
applied to that LOS phase beforehand.

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

#### Error Analysis: Crossing-Geometry Sensitivity

How errors in $p_A$, $p_D$ propagate into $(v_x, v_y)$ depends strongly on the crossing
geometry ($\alpha = H_A - H_D$) and on pixel position, and the dependence is qualitatively
**different** for two distinct error sources. Define $H_{\text{mean}} = (H_A+H_D)/2$ and
$\gamma = \phi - H_{\text{mean}}$ (pixel azimuth relative to the *mean* track heading, rather
than to $H_A$ alone as $\beta$ is). Both results below use the flat-terrain approximation
($\mathbf{B}=0$, i.e. $\mathbf{D}=\mathbf{A}$); a companion tool,
`insarScripts/bin/plotVerticalSensitivity.py`, visualizes both for real crossing-pair
geometries.

**1. Common-mode errors** (the same physical quantity contributing to both $p_A$ and $p_D$ —
e.g. uncompensated vertical motion, see "Vertical-Motion Correction Sign Convention" above,
where $p_A,p_D$ pick up $\delta\cdot\cot\psi_A,\ \delta\cdot\cot\psi_D$ for a common vertical
rate $\delta$). For equal contributions $\delta$ to both measurements, $\mathbf{A}\cdot
(\delta,\delta)^T$ reduces (via sum-to-product identities) to a clean closed form:

$$
\Delta v_x = \frac{\delta\,\cos\gamma}{\cos(\alpha/2)}, \qquad
\Delta v_y = \frac{\delta\,\sin\gamma}{\cos(\alpha/2)}
\qquad\Rightarrow\qquad
\frac{\Delta v_y}{\Delta v_x} = \tan\gamma
$$

The split between $v_x$ and $v_y$ depends **only on $\gamma$** (pixel position relative to
mean heading) — $\alpha$ drops out of the ratio entirely. $\alpha$ instead sets the *overall
amplitude* via $1/\cos(\alpha/2)$, which diverges as $\alpha\to180°$ (near-antiparallel
headings — the typical same-platform ascending/descending case) and is modest
($1/\cos(45°)=\sqrt2$) at $\alpha=90°$ (orthogonal crossing). So near-antiparallel geometry
amplifies a common-mode bias the most, regardless of where the pixel sits; pixel position only
determines *which* of $v_x$/$v_y$ absorbs more of it.

**2. Independent (uncorrelated) errors** in $p_A$ and $p_D$ separately (e.g. random
range/phase measurement noise, $\sigma_{p_A},\sigma_{p_D}$ uncorrelated) follow the
$\sigma_{v_x},\sigma_{v_y}$ formula above — the **row norms** of $\mathbf{A}$, a different
combination than case 1's row *sums*. At $\alpha=90°$ exactly, $\sin\alpha=1,\cos\alpha=0$ and
$\mathbf{A}$ collapses to $\begin{pmatrix}\cos\beta&-\sin\beta\\\sin\beta&\cos\beta\end{pmatrix}$
— an orthonormal **rotation matrix**, so $\sigma_{v_x}=\sigma_{v_y}$ (for
$\sigma_{p_A}=\sigma_{p_D}$) at *every* pixel position: perfectly isotropic error propagation.
Departing from $\alpha=90°$ introduces $\gamma$-dependent anisotropy — e.g. at $\alpha=170°$
(near-antiparallel), the $\sigma_{v_x}/\sigma_{v_y}$ ratio ranges from $\approx0.12$ to
$\approx8$ depending on $\gamma$, vs. exactly $1$ at every $\gamma$ for $\alpha=90°$.

**Practical implication:** crossing geometry near $\alpha=90°$ is doubly favorable — it
minimizes amplification of common-mode bias *and* gives uniform, isotropic sensitivity to
independent measurement noise, regardless of pixel position. Same-platform
ascending/descending pairs ($\alpha$ near $180°$) are the worst case on both counts, though
whether $v_x$ or $v_y$ is hit harder by a given common-mode bias depends entirely on $\gamma$,
not on $\alpha$ itself.

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
3. Applies tidal and submergence/emergence corrections (same $v_z^{\text{SMB}}\cos\psi$ form as
   Step 1, but without the $4\pi/\lambda$ phase-scaling term, since these are range offsets in
   metres):

$$
\Delta r_A \mathrel{+}= v_z^{\text{SMB}}\,\cos\psi_A\,\frac{N_{\text{days},A}}{365.25}, \qquad
\Delta r_D \mathrel{+}= v_z^{\text{SMB}}\,\cos\psi_D\,\frac{N_{\text{days},D}}{365.25}
$$

   See "Vertical-Motion Correction Sign Convention" and "Look-Direction Sign Convention" under
   Step 1 — both apply identically here ($v_z^{\text{SMB}} > 0$ = up; no left/right-looking
   adjustment needed).
4. Scales offsets to horizontal velocity units:

$$
p_A = \frac{365.25\, \Delta r_A}{\Delta t_A \sin\psi_A}, \qquad
p_D = \frac{365.25\, \Delta r_D}{\Delta t_D \sin\psi_D}
$$

5. Constructs $\mathbf{A}$ and $\mathbf{B}$ matrices (same as Step 1) and calls
   `computeVxy` to solve for $(v_x, v_y)$. The crossing-geometry error sensitivity
   ("Error Analysis: Crossing-Geometry Sensitivity" under Step 1) applies identically here —
   same $\mathbf{A}$, same $\gamma$/$\alpha$ dependence — $\sigma_R$ below plays the role of
   $\sigma_{p_A}/\sigma_{p_D}$ there.
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

## Performance

### Multi-threading (OpenMP)

All five pixel-loop routines are parallelised with OpenMP:

| Routine | Per-thread state |
|---|---|
| `makeLandSatMosaic` | No shared mutable image state in pixel path — no per-thread copies needed |
| `make3DMosaic` | Per-thread copies of ascending and descending `inputImageStructure` |
| `make3DOffsets` | Per-thread copies of ascending and descending `inputImageStructure`; per-thread `Aset` flag |
| `makeVhMosaic` | Per-thread copy of `inputImageStructure` |
| `speckleTrackMosaic` | Per-thread copy of `inputImageStructure` |

The outer pixel-row loop (`i`) uses `schedule(dynamic, 8)` — rows are handed out in chunks of 8 to threads as they become free, which handles the non-uniform work distribution (pixels outside the image footprint exit cheaply; interior pixels do geocoding, interpolation, and SVD evaluation).

Per-thread image copies are required because `llToImageNew` writes a warm-start cache field (`lastTime`) into the image struct, and `interpTideError` writes a tide correction field. Making each thread work on its own copy of the struct eliminates these write races.

Before the parallel region for routines that use SVD-based offset interpolation, the lazy-init routines (`svAzOffset`, `svInterpBnBp`) are called once in the serial section to ensure global SVD workspace buffers are allocated before any thread enters the loop.

Thread count is controlled by `-ompThreads N` (default 4). Note that GDAL uses its own internal thread pool for I/O decompression, which can push CPU usage above 100% even at `-ompThreads 1`.

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
