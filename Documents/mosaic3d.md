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
| `-useSquint`                           | Apply per-image squint(r,a) heading correction before building the crossing-orbit solving matrix; phase path only, default off — see "Squint (Residual Doppler) Correction" |
| **Gate flags** — see the **Gates** section for equations and which solver honours each ||
| `-jointMaxSigma <X>`                   | **The sigma cap** (m/yr). Under the hopper it gates the one system holding phase, range *and* azimuth rows; default 35, `0` disables |
| `-jointMaxSigmaPhase <X>`              | Legacy alias for `-jointMaxSigma`, retained so existing templates keep working |
| `-jointMaxSigmaRange <X>`              | The legacy crossing-offsets round's own cap. **Not read in a hopper run**; only meaningful under `-legacyCode`; default 100 |
| `-speckleTrackJoint`, `-jointMaxSigmaSpeckle <X>` | *(retired 2026-09 — accepted and ignored; the hopper supersedes the joint speckle solver. Use `-noPhaseRows` for a speckle-only solve)* |
| `-gateNEff`                            | Use Kish's effective sample size $n_\text{eff}=(\sum w)^2/\sum w^2$ in place of the raw row count in the cap. **On by default**; hopper solvers only |
| `-gateAbsolute`                        | Drop the $\sqrt{n}$ normalisation, making the cap a plain limit on the pixel's own formal sigma. Default off; hopper solvers only |
| `-gateSpeedFrac <F>`                   | Speed-aware cap: effective cap becomes $\max(X,\ F\lVert v\rVert)$, so fast ice is not rejected for having proportionally small error. Default 0.0 = inert; hopper solvers only |
| `-maxChi2 <X>`                         | Reject solved pixels with reduced $\chi^2 > X$. Default −1 = off. **Behaves as a fast-ice filter — keep off or very loose** |
| `-hopper3DMaxSigma <S>`                | 2D/3D toggle for `mosaicHopper3D`; `S < 0` (default) forces the 2D surface-parallel projection everywhere. Not a rejection gate — no pixel is lost to it |
| `-noErrorGate`                         | *(inert — parsed but never read; see "Gates that are not gates")* |
| `-noMask`                              | Ignore any embedded VRT dataset mask band on offset inputs (e.g. `autocleanNISAR.py`'s `range/azimuth.offsets.good` masks); default off, so a mask is honored when present — masked pixels are read as no-data |

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

Since 2026-09 `mosaic3d` solves velocity in **one per-pixel weighted least-squares system** (the
"hopper"). Landsat feature tracking and irregular supplemental data remain separate steps, and the
legacy four-round pipeline is still available via `-legacyCode` (Appendix B).

| Step | Routine | Data | Selected by |
|------|---------|------|-------------|
| 0 | `makeLandSatMosaic` | Landsat optical feature tracking | `-landSat` |
| 1 | `mosaicHopper` / `mosaicHopper3D` | InSAR phase + range offsets + azimuth offsets, one solve | **default** / `-hopper3D` |
| 2 | `addIrregData` | irregularly gridded point observations | `-irregFile` |

Legacy alternative (`-legacyCode`): four successive rounds — crossing phase, crossing range
offsets, phase + azimuth offsets, and pure speckle tracking — each accumulating into the shared
grid by error-weighted averaging. See Appendix B.

#### What a "round" means

The word appears throughout this document and in the source, and it is **legacy vocabulary**. It
is worth being explicit, because the two architectures differ in what they combine:

- A **round** is one complete pass over the data that produces its *own independent velocity
  estimate* on the output grid, then adds it into shared accumulator buffers as
  $\sum_r v_r/\sigma_r^2$ and $\sum_r 1/\sigma_r^2$. After the last round those are divided out
  (`endScale`, `common/scalingFunctions.c`), so a pixel covered by several rounds receives the
  **inverse-variance average of several separately-solved answers**. Each round gates its own
  answer before contributing it.
- The **hopper collapses all four SAR rounds into one.** Every SAR observation — phase, range
  offset, azimuth offset, from every image — is a single row in one per-pixel normal-equation
  system, solved once. There is nothing to average among them, because nothing was solved
  separately.

So the legacy pipeline averages *solutions*; the hopper solves *measurements* jointly. That is the
substance of the 2026-09 change, and it is why the over-count correction ($\rho$) is needed for
one and meaningless for the other.

**But a hopper run is not necessarily single-round overall.** Landsat (`makeLandSatMosaic`) and
irregular supplemental data (`addIrregData`) are still separate passes, and they use the identical
`undoNormalization` → `redoNormalization` → `endScale` accumulation, so their answers are
inverse-variance averaged against the hopper's. A run with `-landSat` has two rounds; with
`-irregFile` as well, three.

What matters for the gates is narrower: **the sigma caps apply only to the SAR-solve round.**
Landsat and irregular data are not gated by them at all. So in a default SAR-only run there is one
gated round and one cap in play, and that stays true however many Landsat or irregular passes are
added.

**`-legacyCode` selects the pipeline SHAPE, not the estimators.** It turns off the hopper and
restores the four-round structure, but each round still uses the *joint* (normal-equation) solver
unless its own legacy-pair flag is also given:

| round | default estimator under `-legacyCode` | legacy pair estimator |
|---|---|---|
| crossing phase | `make3DMosaicJoint` | `make3DMosaic` (`-legacyPairPhase`) |
| crossing range offsets | `make3DOffsetsJoint` | `make3DOffsets` (`-legacyPairRange`) |
| speckle tracking | `speckleTrackMosaic` | — (the joint speckle solver was removed in 2026-09) |

So reproducing a genuinely pre-2026 product needs `-legacyCode -legacyPairPhase -legacyPairRange`,
not `-legacyCode` alone. Appendix B documents the *pair* estimators specifically, which is why the
matrices $\mathbf{A}$ and $(\mathbf{I}-\mathbf{AB})^{-1}$ live there.

Two dispatch traps worth knowing:

- The crossing-range-offsets round runs **only** when `-legacyCode` is set. `-3dOff` on its own
  does not enable it — under the hopper the range offsets are already consumed as rows, and
  running the round as well would enter the same measurements twice.
- `-no3d` short-circuits the hopper (it returns immediately), producing a full set of empty
  results that look like a successful run rather than a failure.

---

### Step 1 — The Hopper (`mosaicHopper`, `mosaicHopper3D`)

Every observation that overlaps a pixel contributes **one row** to a single normal-equation
system. Nothing is paired, and each measurement enters exactly once.

$$ \mathbf N \;=\; \sum_i w_i\,\mathbf a_i\mathbf a_i^{\mathsf T}, \qquad
   \mathbf b \;=\; \sum_i w_i\,\mathbf a_i d_i, \qquad
   \hat{\mathbf v} \;=\; \mathbf N^{-1}\mathbf b $$

with $\mathrm{Cov}(\hat{\mathbf v}) = \mathbf N^{-1}$. Because no pairs are formed there is no
combinatorial over-count to correct, and the $\rho$ parameter of the legacy scheme
(Appendix B, "The crossing-pair over-count correction") is not needed.

#### The three row types

| Observable | $\mathbf a_i$ | $d_i$ | $\sigma_i$ |
|---|---|---|---|
| InSAR phase | LOS direction $(u_x, u_y, u_z)$ | $\phi\cdot 365.25/(2k\,\Delta t)$ | $\sqrt{\texttt{sig2Base} + \min(6\pi,\sigma_{vh})^2}$ scaled |
| Range offset | same LOS direction | $\Delta r\cdot 365.25/\Delta t$ | $\sqrt{\sigma_{\rm off}^2+\sigma_{\rm dem}^2+\texttt{sig2Base}+\sigma_{\rm acc}^2}$ |
| Azimuth offset | $(-\sin\gamma,\;\cos\gamma,\;0)$ | $\Delta a\cdot 365.25/\Delta t$ | $\sqrt{\sigma_{\rm az}^2+\texttt{sig2Off}+\sigma_{\rm acc}^2}$ |

**Note the $\sin\psi$ convention.** In the legacy solvers the phase scaling carries
$1/\sin\psi$ in the *data*; in the hopper that factor lives in the *row* instead, so
$d_i$ omits it. The two are algebraically identical, and mixing them produces a clean
multiplicative scale error. Weights are $w_i = 1/\sigma_i^2$ times the frame's
temporal-overlap weight.

Azimuth rows have $u_z \equiv 0$ — azimuth offsets sense horizontal motion only, which is why
they contribute nothing to the vertical and why the 3D gate (below) does not count them.

#### Row selection

All three types are used by default. `-noPhaseRows`, `-noRangeRows` and `-noAzimuthRows` drop a
type; they are the reduction-test switches and always override a derived value.

Legacy round-selection flags are **translated** rather than ignored, so an existing pair template
runs unchanged and measures what it always measured:

```
phase rows   ON  unless (-no3d AND -noVh)      # crossing round, or vh round
range rows   ON  if (-3dOff OR -rOffsets)      # crossing-offsets round, or speckle round
azimuth rows ON  if (NOT -noVh OR -rOffsets)   # vh round, or speckle round
```

A deprecation line names the modern equivalent. The translation fires only when a legacy flag was
actually given; a bare run gets all three row types.

**Worked examples.** The most common misreading is that a legacy flag still selects a *solver*.
Under the hopper it does not — it selects **rows**, and the solver is always the hopper:

| flags given | rows kept | what actually runs |
|---|---|---|
| *(none)* | phase, range, azimuth | `mosaicHopper` |
| `-rOffsets` | phase, range, azimuth | `mosaicHopper` — **identical to a bare run** |
| `-noVh -rOffsets` | phase, range, azimuth | `mosaicHopper` |
| `-noVh -3dOff` | phase, range | `mosaicHopper` (no azimuth) |
| `-no3d -noVh -rOffsets` | range, azimuth | `mosaicHopper` (no phase) |
| `-hopper3D -hopper3DMaxSigma -1 -noVh -rOffsets` | phase, range, azimuth | `mosaicHopper3D`, forced 2D |
| `-legacyCode -rOffsets` | n/a | `make3DMosaicJoint` **and** `speckleTrackMosaic` — two rounds |
| `-legacyCode -legacyPairPhase -legacyPairRange -3dOff -rOffsets` | n/a | the four legacy pair rounds |

Two consequences worth stating plainly, because both surprise people:

- **`-rOffsets` on its own changes nothing about which rows are used.** All three types are on in
  a bare run, and `-rOffsets` turns all three on. It is not "speckle tracking only" any more —
  under the hopper there is no separate speckle round to select. What it still does is make
  `setup3D` parse the offset files (`setup3D.c:382-400`), which is why it remains **required**
  when you want range or azimuth rows at all.
- **`-offsets` is obsolete and silently ignored** (`mosaic3d.c:1924` prints
  `ignoring obsolet offsets flag - always enabled`). A template carrying it is not doing what its
  name suggests; the range rows in such a run come from `-3dOff` or `-rOffsets`, not from
  `-offsets`.

To get the speckle-tracking round as a *separate solve* you need `-legacyCode`. Without it, step 1 and step 3 are both
skipped entirely — they are gated on `legacyCode == TRUE` — because the hopper has already
consumed the same range and azimuth offsets as rows, and running them again would enter the same
measurements twice.

#### Surface-parallel constraint, and the 3D option

`mosaicHopper` accumulates the $2\times2$ system directly in `float`. `mosaicHopper3D`
accumulates the full $3\times3$ in `double` and recovers the 2D system by projection:

$$ C = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ s_x & s_y \end{pmatrix}, \qquad
   \mathbf N_2 = C^{\mathsf T}\mathbf N_3 C, \qquad
   \mathbf b_2 = C^{\mathsf T}\mathbf b_3 $$

This is an identity, not an approximation, so one accumulator serves both solutions. $s_x, s_y$
come from the DEM at solve time and are zeroed on shelf pixels.

`-hopper3DMaxSigma` selects between them: **$<0$ forces the 2D projection everywhere (the
default)**, $=0$ disables the gate, $>0$ applies an $n$-normalised conditioning gate counting only
rows with vertical sensitivity. Unconstrained 3D is off by default because it was measured to cost
horizontal accuracy in every sector tried (+15.6 %, +70.7 %, +144.4 %) — dropping a mostly-true
constraint can only add variance.

Full derivations: `mosaickingDocuments/hopperDerivation/` (theory.md, implementation.md).

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

## Shared Geometry and Conventions

These apply to **every** solver — the hopper and the legacy rounds alike. They describe how a
lat/lon position becomes a radar geometry and how the sign conventions are fixed; nothing here
depends on which estimator consumes the result.

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
velocity-inversion geometry. It is folded into $H_A$/$H_D$ before anything downstream consumes
them — the per-image sensitivity row $\mathbf{a}_i$ in the joint and hopper solvers, or
$\alpha$/$\beta$ and hence $\mathbf{A}$ on the legacy path — so no further look-direction
handling is needed after this point.

#### Squint (Residual Doppler) Correction (`-useSquint`, off by default)

`computeHeading` returns the heading of the idealized, **zero-Doppler (broadside)** cross-track
direction — real NISAR acquisitions carry a small residual squint (~1.5°–1.7°) that this
geometry doesn't model. With `-useSquint`, each image's own measured squint, a 6-parameter
polynomial in range $r$ and azimuth $a$ fit upstream (`nisarhdf`/`SetupNISAR`, see their
CLAUDE.md) and threaded into the geodat, is evaluated and added directly to that image's own
heading **before** that heading is used to build any sensitivity row or matrix:

$$
H_A \to H_A + \text{squint}_A(r_A, a_A), \qquad H_D \to H_D + \text{squint}_D(r_D, a_D)
$$

This is exact (not a post-hoc rotation of the output $(v_x,v_y)$) because every downstream
geometry term is a pure function of the headings with no other squint dependence — correcting the
headings and changing nothing else reproduces exactly what would have been built from the true
geometry. `evaluateSquint()` (`common/initRoutines.c`) is the single choke point all paths use.
**Phase only** (`make3DMosaic`'s `computeA` call) — `make3DOffsets`'s crossing-orbit range-offset
solution never applies this, flag or no flag, since offsets are self-consistent regardless of
squint by construction (the zero-Doppler condition forces true LOS ⊥ true velocity at the
assigned time, independent of squint). Ships off by default: the underlying analysis found no
NISAR data processed today needs it; verified (numerically and on a real overlapping
ascending/descending pair) to shift the recovered direction by the predicted ~1.5°–1.7° with the
correct sign — see `mosaicSource/CLAUDE.md`'s squint section for the full derivation and
verification.

#### Surface-Slope Coupling (`computeB`) — used by **every** solver

Ice flowing over a sloping surface has a vertical component
$v_z = v_x\,\partial z/\partial x + v_y\,\partial z/\partial y$, which projects into the
line of sight and must be removed to recover the horizontal velocity. `computeB`
(`common/initRoutines.c`) supplies that coupling for all solvers — legacy, joint and hopper
alike — though they consume it in two different shapes.

Slopes come from the DEM by centred finite differences over a spacing of at least 90 m and are
clamped by `limitSlope` to $\pm 0.25$ ($\approx 14°$):

$$
\frac{\partial z}{\partial x} = \mathrm{clamp}\!\left(\frac{z(x{+}\tfrac{dx}{2}) - z(x{-}\tfrac{dx}{2})}{dx},\ \pm 0.25\right)
$$

If any of the four DEM samples is invalid, **both** slopes are set to zero rather than the pixel
being rejected. The clamp was raised from 0.1 in 2026-07: `detC` (below) is a purely geometric
quantity, so letting the real slope through lets the conditioning test reject genuinely singular
geometry instead of the clamp hiding it.

The returned matrix is

$$
\mathbf{B} =
\begin{pmatrix}
\dfrac{\partial z/\partial x}{\tan\psi_A} & \dfrac{\partial z/\partial y}{\tan\psi_A} \\[8pt]
\dfrac{\partial z/\partial x}{\tan\psi_D} & \dfrac{\partial z/\partial y}{\tan\psi_D}
\end{pmatrix}
$$

**Two consumption patterns.** The distinction matters when reading the source:

| caller | call | how it is used |
|---|---|---|
| legacy pair (`make3DMosaic`, `make3DOffsets`) | `computeB(..., aPsi, dPsi, ...)` | the full 2×2, one row per image of the pair, inside $(\mathbf{I}-\mathbf{AB})^{-1}\mathbf{A}$ — see Appendix B |
| joint and hopper solvers | `computeB(..., psi, psi, ...)` | **one image's** $\psi$ in both slots, so the two rows are identical; only row 0 is read, as the per-image row correction $\mathbf{a}_i = (\cos\gamma_i,\ \sin\gamma_i) - \mathbf{B}_i$ |

So there is no separate "3D slope model": the same $\cot\psi\,(\partial z/\partial x,\ \partial z/\partial y)$
coupling appears in all of them, subtracted from each image's own sensitivity row rather than
assembled into a pair matrix.

**Ice-shelf carve-out.** Where the shelf mask marks a pixel as `SHELF`, $\mathbf{B}$ is zeroed
(the slopes themselves are left alone for the $v_z$ calculation). True shelf-interior slope is
near zero, so a large DEM slope there is almost always a stale rift or calving front — a DEM
error, which unclamping would make worse rather than better.


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
  into $H_A$/$H_D$ — and hence into every sensitivity row or matrix built from them — before the
  correction terms are ever applied to $\phi_A$/$\phi_D$.

In other words, look direction changes how a given LOS phase gets decomposed into
$(v_x, v_y)$, not the sign or magnitude of the vertical-motion correction applied to that LOS
phase beforehand.

### Step 2 — Irregularly Gridded Supplemental Data (`addIrregData`)

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

### Interpreting the formal errors (`ex`, `ey`)

`mosaic3d` propagates per-observation sigmas through the inversion and reports the result as
`ex`/`ey`. Under `-legacyCode` a crossing-pair over-count correction is also applied — that
correction is specific to the pair solvers and is documented in Appendix B.

**In one line: the reported error is globally correct, spatially wrong, and not Gaussian.** Its
overall level is calibrated against an independent reference; it does not identify *which* pixels
are bad; and it does not define a confidence interval at the usual multipliers.

#### What the reported error does and does not mean

**Globally correct.** Across a 2.5× change in pair count and three solution types, $k$ stays within
a few percent of 1 for phase and the combined product. Aggregated over a region, `ex`/`ey` predicts
the observed scatter.

**Spatially wrong.** The budget contains per-observation terms (matching noise, tie-point fit
residual, baseline covariance, DEM error) but **no term for gross failures** — unwrapping errors,
correlation mismatches, bad frames. Those are spatially localised and carry a disproportionate share
of the variance (for crossing offsets in $v_y$, the worst 1% of pixels hold ~48%). No
per-observation model can identify *which* pixels those are, so `ex`/`ey` is nearly flat — p10 to
p90 spans a factor of 2.5 — where the true error varies by more than an order of magnitude.

**Not Gaussian.** Below, $z=d/e$ is rescaled so $\mathrm{RMS}(z)=1$ exactly — i.e. assuming the
scale has been tuned so $k=1$ — which isolates the distribution's *shape* from the question of
overall scale.

**Table — Coverage and confidence multipliers.** Left: fraction of pixels whose actual error falls
within 1, 2 and 3 times the reported sigma. Right: the multiplier of the reported sigma needed to
enclose the stated fraction. $T=10000$, 409,476 px of stable ground, $v_x$ and $v_y$ pooled.
Measured on legacy pair products at $\rho$ = 0.6, but the *shape* is a property of the
per-observation error distribution and carries over to the hopper.

| product | <1σ | <2σ | <3σ | 68.3% | 90% | 95% | 99% | 99.9% |
|---|---|---|---|---|---|---|---|---|
| phase | 86.4% | 96.3% | 98.3% | 0.55 | 1.19 | 1.71 | 3.89 | 9.66 |
| offsets | 88.0% | 95.9% | 98.0% | 0.62 | 1.09 | 1.75 | 4.11 | 8.60 |
| both | 83.6% | 95.9% | 98.3% | 0.63 | 1.30 | 1.82 | 3.75 | 8.46 |
| *Gaussian* | *68.3%* | *95.4%* | *99.7%* | *1.00* | *1.65* | *1.96* | *2.58* | *3.29* |

The core is **tighter** than Gaussian — 84–88% of pixels inside 1σ against 68% — while the tail is
much heavier: 99% coverage needs 3.8–4.1σ against 2.58, and 99.9% needs 8.5–9.7σ against 3.29.
Roughly 6× as many pixels sit beyond 3σ as normality predicts.

This is **orthogonal to any overall scaling.** The legacy $\rho$ (Appendix B) governs how the error
scales with *pair count*; the tail is a property of the per-observation error distribution and would
be present with a single pair — or, under the hopper, with a single row.

Practical guidance:

- Use `ex`/`ey` for **relative weighting** — inverse-variance combination depends only on the ratio
  between contributions, which is the part that is reliable.
- Treat a single pixel's value as an estimate of the **local noise level**, not of that pixel's
  actual error.
- A 2σ envelope is roughly honest (95.9–96.3%). Anything quoted at **99% or beyond is wrong by a
  factor of several** unless the multipliers above are used.
- For **RMS-based requirement verification**, quote the measured difference against an independent
  reference together with the outlier fraction, rather than substituting `ex`/`ey`.
- A large `ex`/`ey` is informative; a small one is **not a guarantee**, since the dominant tail
  failure modes are invisible to the budget.

---

## Gates

A "gate" here is any test that can **reject a pixel** (or a whole pair) after the geometry and
the measurements are in hand. They are the main reason two runs over the same data return
different coverage, so they are collected here rather than scattered through the solver
descriptions.

Gates fall into three groups by **where** they act:

| level | acts on | survives to output? |
|---|---|---|
| **image / frame** | one input product, before any of its pixels are used | that product contributes nothing anywhere; others still do |
| **pair / scene** | a candidate pair of images, before any pixel is visited | pair is skipped entirely; other pairs still contribute |
| **pixel, per solver round** | one output pixel within one round (phase, range, speckle…) | pixel gets no contribution from that round; other rounds may still fill it |
| **final product** | the assembled mosaic | pixel is blank in the delivered product |

**Nothing here is a quality flag on the data itself.** Every gate below tests *geometry* or
*internal consistency*, not whether the measurement agrees with any external truth.

There are currently **no final-product gates inside `mosaic3d`** — the third row is listed because
it is where one would go, and because `-noErrorGate`'s help text wrongly implies one exists (see
"Gates that are not gates"). Masking and gap-filling of the delivered product happen downstream in
`mosaicworkflow`.

### Image-level gates

These reject a whole input product rather than a pixel, and they are **the only sigma-type gates
the legacy path has.** They are cheap, they fire before any geometry is computed, and they apply
across solvers unless noted.

**A. Azimuth-residual threshold** (`-sigmaAThresh X`, metres; default **1000**, i.e. effectively
off). If an image's azparams tie-point fit residual exceeds $X$, its azimuth information is
dropped:

$$
\text{drop if } \sigma_A^{\text{residual}} > X
$$

What "dropped" means depends on the solver, and the difference matters:

| solver | effect |
|---|---|
| `speckleTrackMosaic` (legacy step 4), `mosaicTrue3D` | **the whole image is skipped** — its range offsets are discarded too |
| `mosaicHopper`, `mosaicHopper3D` | only the azimuth row is dropped |

So the same flag costs strictly more coverage on the legacy path than under the hopper.

**B. Azimuth threshold in velocity units** (`-sigmaAThreshVel X`, m/yr; default **150**, `<0`
disables). The same residual rescaled by the pair separation, $\sigma_A^{\text{residual}}\cdot
365.25/n_{\text{days}}$, which is the quantity that actually matters for a velocity. **Hopper
solvers only** — the legacy path has no equivalent, so a 12-day and a 90-day pair with the same
metre-level residual are treated identically there.

**C. No-solution sentinels.** `rparams`/`azparams` write $\sigma < 0$ when their fit found no
solution. Consumers skip that product: `make3DOffsets` nulls `rFile` on
`sigmaRresidual < 0`, `speckleTrackMosaic` skips on `sigmaAresidual < 0` (a separate test,
because a negative value is never greater than a positive threshold). Not tunable, and not
really a gate so much as a validity check.

**D. Partner weight floor.** `make3DOffsets` skips a descending partner with
`weight < 0.05` — a hard-coded constant, legacy crossing-offsets path only.

### Pair-level gates (legacy path only)

Both are inside the legacy pair solvers; the joint and hopper solvers loop over single images and
never form a pair, so neither test exists for them.

**1. Scene heading separation** (`computeSceneAlpha`, `common/initRoutines.c`). Before any pixel
of a candidate crossing pair is visited, the two scene-centre cross-track headings are compared:

$$
\alpha = H_A - H_D, \qquad \text{skip the pair if } \lvert\alpha\rvert < 30°
\ \text{ or }\ \lvert\alpha\rvert > 360° - 30°
$$

`MINCROSSINGHEADINGSEP` (`common/common.h`) is currently 30°, reduced from 40° to recover
coastal coverage. Two near-parallel look directions cannot separate $v_x$ from $v_y$ at any noise
level, so the pair is rejected wholesale rather than per pixel.

**2. Per-pixel heading separation** (`computeA`). The same test, re-applied at each pixel with
that pixel's own local headings, at a tighter threshold:

$$
\text{reject the pixel if } \lvert\alpha\rvert < 0.8\ \text{rad}\ (\approx 46°)
$$

The two thresholds differ on purpose: the scene test is a cheap reject of hopeless pairs, the
pixel test is the one that actually protects the inversion.

### Pixel-level gates

**3. Conditioning / minimum eigenvalue.** Every non-legacy solver builds a per-pixel normal
system $\mathbf{N} = \sum_i w_i\,\mathbf{a}_i\mathbf{a}_i^{\mathsf T}$ and rejects the pixel
outright if it is not positive definite:

$$
\det \mathbf{N} > 0 \quad\text{and}\quad \lambda_{\min}(\mathbf{N}) > 0
$$

For the 2×2 case $\lambda_{\min} = \tfrac{1}{2}(N_{xx}+N_{yy}) - \sqrt{\tfrac14 (N_{xx}-N_{yy})^2 + N_{xy}^2}$.
The worst-direction formal sigma follows directly:

$$
\sigma_{\text{worst}} = \frac{1}{\sqrt{\lambda_{\min}}}
$$

**4. The sigma cap** (`-jointMaxSigma`). The main coverage control. A pixel is kept only if

$$
\sigma_{\text{worst}} \cdot g \ \le\ X
$$

where $X$ is the cap in m/yr ($X = 0$ disables the gate) and $g$ is a normalisation factor set by
the modifiers below.

**`-jointMaxSigma` is the flag to use.** Under the hopper there is one system holding phase, range
*and* azimuth rows, and one cap gating it — so a single name is the honest description. The
per-round variants below exist because the *legacy* pipeline solved each observable in its own
round and gated each separately (see *What a "round" means* in the Overview); they are **not** a
per-observable choice.

| variable | flag | read by | governs |
|---|---|---|---|
| `jointMaxSigma` | `-jointMaxSigma`<br>`-jointMaxSigmaPhase` *(legacy alias)* | **both hoppers**, `make3DMosaicJoint` | the hopper's single combined system; under `-legacyCode`, the crossing-phase round |
| `jointMaxSigmaRange` | `-jointMaxSigma` *(sets both)*<br>`-jointMaxSigmaRange` *(per-round)* | `make3DOffsetsJoint` | the legacy crossing-offsets round only |

Consequences, in order of how often they bite:

- **In a default (hopper) run, only `jointMaxSigma` is ever read.** `jointMaxSigmaRange` is not
  consulted at all — neither hopper source file mentions it. Setting `-jointMaxSigmaRange` in a
  hopper template does nothing.
- **Two caps apply in the same run only under `-legacyCode`** (without the `-legacyPair*` flags),
  where the rounds run in sequence and each gates *its own contribution* before it is
  error-weighted into the shared grid. Even then they never both act on the same measurement: a
  pixel can be accepted by the phase round and rejected by the offsets round, in which case it
  survives carrying phase information only.
- **Under `-legacyCode -legacyPairPhase -legacyPairRange` neither applies.** The legacy pair
  estimators use $|\alpha|$ and $\det\mathbf{C}$ instead.

The defaults differ (35 and 100) because they were sized against different observables, not
because one is stricter in spirit: an offset error is proportionally smaller on fast ice, so the
same absolute number means something different.

**History (2026-09).** The variable behind `-jointMaxSigma` was called `jointMaxSigmaPhase`, and
the log header printed `jointMaxSigPhase`, `jointMaxSigRange` and `jointMaxSigSpeck`
unconditionally — which made a per-round split look like a per-observable one, in runs where only
the first was read. The variable is now `jointMaxSigma`, `-jointMaxSigmaPhase` is retained as an
alias so existing templates keep working, and `jointMaxSigmaSpeckle` went with the joint speckle
solver. **Logs written before this change carry the old key names.**

> **Reading a log:** the header prints `; jointMaxSigma :` — the cap that acted. A
> `; jointMaxSigRange :` line appears only under `-legacyCode`, where it can genuinely differ.
> To confirm which solver ran, read the `; SOLVER:` line and that solver's own cap line, e.g.
> `; mosaicHopper3D jointMaxSigma : 50.000000`.

**Modifiers to $g$.** These are **two independent choices, not three alternatives** — a point the
flag names obscure. One picks *which* $n$; the other decides whether $n$ is used at all:

```
nGate = gateNEff ? (sum w)^2 / sum(w^2)   /* Kish effective count -- ON by default */
                 : nObs                   /* raw row count       -- -noGateNEff    */

g     = gateAbsolute ? 1.0                /* -gateAbsolute: nGate is computed and DISCARDED */
                     : sqrt(nGate)        /* default                                        */
```

So the four combinations collapse to three distinct behaviours, and `-gateAbsolute` makes
`-gateNEff` moot:

| flags | $g$ | |
|---|---|---|
| *(compiled default)* | $\sqrt{n_{\text{eff}}}$ | `gateNEff` is **on**, so the default is the effective count, not the raw one |
| `-noGateNEff` | $\sqrt{n_{\text{obs}}}$ | raw row count |
| `-gateAbsolute` | $1$ | plain cap on the pixel's own formal sigma; the `nEff`/`nObs` choice has no effect |

The n-normalised forms read $\sigma_{\text{worst}}\sqrt{n}$ as the *per-measurement* sigma, so the
effective cut tightens as $X/\sqrt{n}$ — it asks "is coverage thin?". The absolute form asks "is
this measurement noisy?". Because $\sigma_{\text{worst}}$ carries geometric dilution as well as
noise, the n-normalised form rejects geometry-limited pixels however good their data is.

`-gateNEff` exists because the $\sqrt{n}$ reading is only valid when rows carry comparable weight,
and they do not: azimuth rows sit ~4 orders of magnitude below phase rows (measured
$7\times10^{-6}$ against $4\times10^{-2}$), contribute nothing to $\lambda_{\min}$, yet pad $n$ and
inflate the gate ~1.7×.

**Note for reading production logs:** the Greenland and Antarctic recipes pass `-gateAbsolute`, so
in those runs the logged `gateNEff : 1` is inert.

The n-normalised default asks *"is coverage thin?"*; the absolute form asks *"is this measurement
noisy?"*. Because $\sigma_{\text{worst}}$ carries geometric dilution as well as noise, the
n-normalised form rejects geometry-limited pixels however good their data is.

**5. The speed-aware cap** (`-gateSpeedFrac F`, default 0.0 = inert). A fixed cap in m/yr is a
tightening constraint as speed rises: 50 m/yr is 0.7% of a 7 km/yr velocity but 50% of a
100 m/yr one. With $F > 0$ the cap becomes

$$
X_{\text{eff}} = \max\!\left(X,\ F\,\lVert \mathbf{v} \rVert\right)
$$

so a 7000 m/yr pixel at $F = 0.03$ is allowed 210 m/yr of error instead of 50. This is why the
test is evaluated **after** the solve — the speed is not known before it. With $F = 0$ the
ordering is immaterial and the result is identical to the pre-solve form.

**6. Reduced chi-square** (`-maxChi2 X`, default −1 = off). A blunder screen applied after the
solve, asking whether the measurements at a pixel agree *with each other*:

$$
\chi^2_\nu = \frac{S_{dd} - \mathbf{v}\cdot\mathbf{b}}{n_{\text{obs}} - n_{\text{par}}},
\qquad \text{reject if } \chi^2_\nu > X
$$

$S_{dd} = \sum_i w_i d_i^2$, so this is free — no extra pass. **Keep it loose or off.** Measured
on real mosaics it behaves as a fast-ice filter rather than a blunder screen: fast ice genuinely
disagrees more between observations, so the screen removes the glaciers preferentially, and the
fast pixels it keeps have *higher* formal error than the ones it discarded. Typical measured
retention at $X = 100$: >99% below 100 m/yr, but 3–16% between 1 and 2 km/yr.

**7. The 3D/2D toggle** (`-hopper3DMaxSigma`, default −1). Not a rejection gate — **no pixel is
lost to it.** It decides whether `mosaicHopper3D` solves all three components at a pixel or
projects onto the surface-parallel plane:

$$
\text{3D if } \lambda_{\min}(\mathbf{N}_3) > 0,\ n_{\text{obs}} \ge 3,\ n_{\text{range}} \ge 2,\ \text{and } \frac{g}{\sqrt{\lambda_{\min}(\mathbf{N}_3)}} \le S
$$

with $S$ = `-hopper3DMaxSigma`; $S < 0$ forces 2D everywhere (the default), $S = 0$ means "3D
whenever it is solvable". $n_{\text{range}}$ counts **phase + range** rows only: azimuth rows
have $u_z \equiv 0$ and say nothing about the vertical, so crediting them would let measurements
vouch for a solve they cannot inform. The `.mode` diagnostic band records which branch ran (3 or
2) per pixel.

**8. Reference-velocity clip** (`-refVel <file> -clipThresh X`, off unless both are given). An
outlier screen against an independent velocity map rather than against the data's own statistics.
At each solved pixel the reference is interpolated and the pixel rejected when

$$
\lVert \mathbf{v} - \mathbf{v}_{\text{ref}} \rVert > X
\quad\text{and}\quad
\bigl(\lVert \mathbf{v}_{\text{ref}} \rVert < 100 \ \text{ or }\ \lVert \mathbf{v} \rVert < 100\bigr)
$$

The speed condition confines it to slow ice, so fast outlets are never clipped against a reference
that may be from a different epoch. A pixel with no reference value (outside the grid, or no-data)
is **kept**. Applied after the solve and after `-maxChi2`, since it needs the final velocity;
counted separately in the logs as `rejClip`.

Available in the hopper solvers since 2026-09 (`clipVelChecked`, `speckleTrackMosaic.c`) as well as
the legacy speckle round, which has always had it (`clipVel`). The two differ in one respect: the
legacy `clipVel()` ignores `refVelInterp()`'s return value, so where there is no reference value it
compares against uninitialised stack. That is a real defect, left in place so the legacy path stays
bit-for-bit reproducible; the hopper version checks the return.

**9. Legacy conditioning** (`computeVxy`). The legacy pair path's own version of gate 3:

$$
\mathbf{C} = \mathbf{I} - \mathbf{A}\mathbf{B}, \qquad \text{reject if } \det\mathbf{C} < 0.5
$$

Raised from 0.25 in 2026-07 alongside the slope-clamp change; the two were tuned together.

### Which solver honours which gate

This table is the one worth checking before interpreting a run — **the modifiers are honoured by
the hopper solvers only.** Passing `-gateAbsolute` to a joint or legacy run is silently inert.

**The legacy pair solvers have no per-pixel sigma gate at all.** Nothing in `make3DMosaic` or
`make3DOffsets` rejects a pixel because its computed error came out large; their per-pixel
rejections are purely geometric ($\lvert\alpha\rvert$ and $\det\mathbf{C}$). Consequently
`-jointMaxSigma`, `-jointMaxSigmaRange`, `-maxChi2` and every gate modifier have **zero
effect** on a full `-legacyCode -legacyPairPhase -legacyPairRange` run. What legacy does have is
the image-level gates above — it screens its *inputs* on their fit residuals, where the hopper
screens its *output* on the solution's formal sigma.

| gate | hopper / hopper3D | joint (`make3DMosaicJoint`, `make3DOffsetsJoint`) | legacy pair |
|---|:--:|:--:|:--:|
| `-sigmaAThresh` (image) | ✓ (az row) | ✓ (az row) | ✓ (**whole image**) |
| `-sigmaAThreshVel` (image) | ✓ | — | — |
| $\sigma<0$ sentinels (image) | ✓ | ✓ | ✓ |
| partner `weight < 0.05` | — | — | ✓ (offsets) |
| scene heading separation (30°) | — | — | ✓ |
| per-pixel $\lvert\alpha\rvert \ge 0.8$ rad | — | — | ✓ |
| $\det\mathbf{C} \ge 0.5$ | — | — | ✓ |
| positive-definite $\mathbf{N}$ | ✓ | ✓ | — |
| sigma cap $X$ | ✓ (`jointMaxSigma`) | ✓ (per round) | — |
| `-gateNEff` | ✓ | — | — |
| `-gateAbsolute` | ✓ | — | — |
| `-gateSpeedFrac` | ✓ | — | — |
| `-maxChi2` | ✓ | — | — |
| `-hopper3DMaxSigma` | ✓ (3D only) | — | — |
| `-refVel`/`-clipThresh` clip | ✓ *(since 2026-09)* | — | ✓ (speckle round) |

### Gates that are not gates

- **`-noErrorGate` is inert.** The flag is parsed and sets a global that **nothing reads**
  (`grep noErrorGate` finds only the declaration, the parse arm and the usage string). Its usage
  text — "default is to remove them" — describes behaviour that does not exist: `toSigma()`
  deliberately does *not* gate velocity on the error, because the downstream workflow gap-fills
  `.vx`/`.vy` but not `.ex`/`.ey`, so a missing error is the flag marking an interpolated pixel.
  The flag and its help text should be removed or implemented; documented here so nobody
  concludes from the usage string that a gate exists.
- **`-timeThresh` / `-timePhaseThresh`** are pair-formation windows on the legacy path, not
  gates. Under the hopper they are **ignored** — the solver forms no pairs, and temporal extent
  is set by `-date1`/`-date2`. The hopper prints a note to stderr when a non-default value is
  passed. Before 2026-09-07 they were mis-applied there as a per-image window against the mosaic
  centre date, which silently discarded most phase rows for small values.
- **The slope clamp** ($\pm 0.25$) bounds an input rather than rejecting a pixel; a pixel with
  extreme slope is still solved, and it is `detC`/$\lambda_{\min}$ that reject it if the
  resulting geometry is singular.


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

---

## Appendix B — The legacy pair solvers (`-legacyCode`)

The four rounds below were the default until 2026-09, and remain available via `-legacyCode`.
They are retained in full because the equations are still correct for that path, and because
the hopper's reduction tests are defined against them.

They share the geometry and sign conventions of the main text; only the estimator differs.

### The crossing-pair over-count correction

Applies **only** to the two pair solvers below: `inflatePairOverCount()`
(`common/scalingFunctions.c`) is called from exactly two places, `make3DMosaic.c` and
`make3DOffsets.c`. The joint and hopper solvers accumulate each measurement once and form no
pairs, so there is nothing to over-count and none of `-rhoPhase`, `-rhoOffsets`,
`-pairCountLegacy` or `-noPairOverCount` has any effect on them.

#### Quantities

Per output pixel, within one mosaicking round:

| symbol | meaning |
|---|---|
| $n_A$ | distinct **ascending** images contributing (outer-loop images, deduped via `aContrib`) |
| $n_D$ | distinct **descending** images contributing, derived as $P/n_A$ |
| $P$ | contributing **pairs** — the quantity actually accumulated |
| $\rho$ | share of the per-pixel variance **common** to every pair, so never averaged away |
| $f$ | factor by which this round's error accumulator is inflated |

#### The over-count correction

Every crossing pair is accumulated as an independent observation, but a pixel seen by $n_A$
ascending and $n_D$ descending images yields only about $n_A+n_D$ independent measurements, not
$n_A n_D$. After each round, that round's own contribution to the error accumulator is inflated by

$$f = \rho\,P + (1-\rho)\,\frac{n_A+n_D}{2}, \qquad f \ge 1$$

$\rho=0$ credits the full averaging that $(n_A+n_D)/2$ implies — correct when the error is
independent per observation. $\rho=1$ asserts the error is entirely common to every image at the
pixel, so nothing averages and the whole pair count is spurious. Both default to **0.6**
(`-rhoPhase`, `-rhoOffsets`).

Nominally $\rho$ is the share of variance common to *every* image, but in practice it does a larger
job: it stands in for correlation between pairs that **share an image**. A frame's baseline
residual attaches to that frame, so all $d$ pairs the frame joins inherit it — the shared-error
floor derived in `crossingOrbitRedundancy.md` §2.2. The exact treatment is $f=\sum_i d_i^2/(2P)$
over contributing images (which reduces to $(n_A+n_D)/2$ for a complete bipartite graph, i.e. the
$\rho=0$ form above), but that needs a per-pixel record of *which* images contributed. $\rho$ is
the fixed-cost stand-in for it.

**$\rho$ corrects the reported error, not the weighting.** The accumulated weight still grows like
$P$, so the crossing-orbit method still out-votes the other methods in the error-weighted blend by
a bookkeeping artefact — `crossingOrbitRedundancy.md` §2.2 "damage 2". That document sets out two
practical routes to fixing the underlying problem rather than rescaling its symptom: a block-local
minimum edge cover over the pair list (§4), and a per-pixel normal-equation accumulator that makes
double-counting impossible by construction (§5). Neither is implemented; $\rho$ is what ships today.

Only the error accumulator is scaled — the velocity is $v_{X}^{\text{image}}/\text{scale}_X$, so
$\rho$ cannot move it at any value (verified $\max|dv_x| = 0$). `-noPairOverCount` disables the
correction entirely; `-pairCountLegacy` restores the pre-2026-08-26 formula.

#### How $\rho$ = 0.6 was chosen

Raising the crossing time threshold $T$ admits more **pairs** while the contributing **images** stay
fixed — image count moves 0.5% from $T=12$ to $T=37$ while the pair count moves 2.5×. Those extra
pairs carry no new information, so the measured accuracy is flat across $T$, and a correctly scaled
formal error must be flat too. That makes the threshold response a measurement of $\rho$ that is
independent of the absolute error level. See the appendix, *Calibrating the over-count parameter*.

![Formal error at the shipped rho against the measured difference](rhoThresholdShipped.png)

**Figure — the shipped $\rho$ = 0.6, measured.** One column per solution type. Top: $\mathrm{RMS}(e)$,
the population prediction from the per-pixel formal errors, against $\mathrm{std}(d)$, the measured
scatter of the difference from a Sentinel-1 reference over stable ground (S1 speed $\le$ 50 m/yr,
406,845 px common to all twelve runs, 1600 m grid). Bottom: $k=\mathrm{std}(d)/\mathrm{RMS}(e)$;
1.0 is calibrated. Every point is a real run at $\rho$ = 0.6 — nothing is interpolated. $T=6$ is
shaded because ~10% of images have no partner within 6 days and drop out entirely, breaking the
same-images premise.

Phase is calibrated to within 1–3% and flat; the combined product to within 1–6% and flat. Crossing
offsets alone are flat but over-corrected by ~20%; no single $\rho$ lifts them without breaking the
other two, and offsets-only is not a shipped NISAR product.


### Calibrating the over-count parameter

How $\rho$ = 0.6 was arrived at. Full write-up, including the runs that were discarded and why, in
`Release/velocity/errCal/report/errCalResults.md` in the Greenland project directory. Related design
work on fixing the underlying redundancy rather than rescaling it: `crossingOrbitRedundancy.md`.

**The experiment.** Raising the crossing time threshold $T$ admits more pairs while the contributing
images stay fixed — from $T=12$ to $T=37$ the image count moves 0.5% while the pair count moves
2.5×. Those pairs carry no new information, so the measured error is flat across $T$; a correctly
scaled formal error must be flat too. The threshold response therefore measures $\rho$ *independently
of the absolute error level*, which matters because the absolute level is also affected by terms the
budget is missing entirely.

**Statistic.** $k=\mathrm{std}(d)/\mathrm{RMS}(e)$, where $d$ is the difference from a Sentinel-1
reference over stable ground and $\mathrm{RMS}(e)=\sqrt{\overline{e^2}}$ is the population
prediction from the per-pixel formal errors. $\mathrm{RMS}(e)$, not $\mathrm{mean}(e)$ or
$\mathrm{median}(e)$: variances average, standard deviations do not. Using a median here biases $k$
high by 25–40%, and an early version of this analysis did exactly that.

**Result 1 — the endpoints are excluded.** At $\rho=0$ the reported error falls ~2× across the
threshold range while the real error does not; at $\rho=1$ it overshoots ($k$ = 0.66–0.88). The
answer is interior and constrained from both sides.

**Result 2 — bracketing.** Comparing the formal-error ratio $\mathrm{RMS}(e)_{T=10000}/
\mathrm{RMS}(e)_{T=12}$ against the measured $\mathrm{std}(d)$ ratio gives a crossing at
$\rho$ = 0.76 (phase), 0.85 (both), 0.87 (offsets) — bracketed by real runs, not extrapolated.
$\mathrm{RMS}(e)^2$ is linear in $\rho$ to better than 0.5%, which makes the interpolation exact in
practice.

**Result 3 — an independent check on a different pairing graph.** A sandbox build that consumes each
image at most once per pixel forms a *matching*, which has no over-counting by construction. Its
over-count factor is not $(n_A+n_D)/2$ — that form is specific to the complete bipartite graph —
but $f=\rho m + (1-\rho)$ over the $m$ disjoint pairs. Calibrated that way it gives $k=1$ at
$\rho$ = 0.584 (phase) and 0.501 (offsets). Two different criteria on two different graphs both land
near 0.6.

Note this places 0.6 at the **top** of the supported 0.50–0.87 range rather than at its centre, and
the two criteria do not agree exactly — $\rho$ is a single scalar standing in for correlation
structure that is not really one number. A matching is also a worse *estimator* (it strands the
surplus when the two sides are unbalanced) and is a diagnostic only, not a candidate design; see
`crossingOrbitRedundancy.md` §4.3 for why an edge cover rather than a matching is the right
selection rule.

**What $\rho$ cannot fix.** The heavy tail. It is a property of the per-observation error
distribution — localised blunders the budget contains no term for — and would be present with a
single pair. Raising $\rho$ far enough would make the reported sigma match the observed scatter, but
only by using a redundancy parameter to absorb a blunder population instead of modelling it.


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
   directions (optionally squint-corrected, see "Squint (Residual Doppler) Correction" in
   *Shared Geometry and Conventions* — off by default), and the surface-slope correction matrix
   $\mathbf{B}$ from the DEM.
6. Solves for $(v_x, v_y)$ and derives $v_z = v_x \partial z/\partial x + v_y \partial z/\partial y$.
7. Propagates baseline covariance to a per-pixel phase error $\sigma_\phi$.

#### Matrix A — Geometric Conversion (`computeA`) — legacy pair path only

Used **only** by `make3DMosaic` and `make3DOffsets`. The joint and hopper
solvers never call it: a per-image sensitivity row replaces the pair matrix
(see *Algorithm*), and the identity $\mathbf{A} = \mathbf{N}^{-1}$ is what makes
the two agree exactly at $n = 2$.

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

If $\det(\mathbf{I} - \mathbf{AB}) < 0.5$ (poorly conditioned due to extreme slopes), no
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
