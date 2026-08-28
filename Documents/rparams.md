# rparams — Range Parameter Estimator

## Purpose

Calibrates range offset measurements derived from speckle tracking by estimating
interferometric baseline parameters (Bn, Bp, and their along-track derivatives) using
tiepoints with known locations and surface velocities. Output is written to stdout in
a format suitable for downstream mosaicking tools.

---

## Usage

```
rparams [options] geodatFile tiepointsFile offsetFile baselineFile
```

### Required Arguments (positional, in order)

| Argument        | Description |
|-----------------|-------------|
| `geodatFile`    | SAR image geometry parameter file |
| `tiepointsFile` | Tiepoint file with lat, lon, elevation, and velocity (vx, vy, vz) |
| `offsetFile`    | Binary range offset file (a companion `.dat` metadata file must also exist) |
| `baselineFile`  | CW state-vector baseline file |

### Options

| Option               | Description |
|----------------------|-------------|
| `-nDays <days>`      | Temporal baseline in days (default: 24) |
| `-shelfMask <file>`  | Ice-shelf mask file for tidal displacement corrections |
| `-constOnly`         | Estimate only the constant range offset term |
| `-bnbpOnly`          | Estimate Bn, Bp, and dBp (hold dBn fixed) |
| `-bnbpdBpOnly`       | Estimate Bn, Bp, and dBp |
| `-bpdBpOnly`         | Estimate Bp and dBp only |
| `-quadB`             | Estimate quadratic along-track baseline terms (dBnQ, dBpQ) |
| `-deltaBQ`           | Estimate quadratic correction to state-vector baseline |
| `-deltaBC`           | Estimate constant correction to Bp component of baseline |
| `-quiet`             | Suppress tiepoint echo to stdout |
| `-outputFile <path>` | Write the solution to `<path>` instead of stdout. Mutually exclusive with `-runFile` (which already names an output per run). |
| `-debug`             | Write every tie point used in the fit, plus its residual, to a GeoPackage (see "Debug output" below) |
| `-noMask`            | Ignore any embedded VRT dataset mask band on the offset file (e.g. `autocleanNISAR.py`'s `range.offsets.good` mask); default off, so a mask is honored when present |

> **Note:** `-bnbpOnly`, `-bpdBpOnly`, and `-bnbpdBpOnly` are mutually exclusive.

### Output

Results are written to **stdout** (or to `-outputFile`, if given) in a commented format:

- Tiepoint locations and extracted range offsets (with `;` prefix)
- Fit residual sigma (`sigma*sqrt(X2/n)`)
- 6×6 parameter covariance matrix
- Final estimated baseline parameters: `Bn  Bp  dBn  dBp  const  dBnQ  dBpQ`

### Debug output (`-debug`)

Writes a `.gpkg` point layer named `residuals` containing every tie point used in the
fit: `id`, `lat`, `lon`, `x_km`/`y_km` (polar-stereo), `range`, `azimuth`, `z`,
`weight`, and `range_residual_m` (the fit residual in meters). Filename defaults to
`rparams.<mode>.gpkg` (`<mode>` reflects the active `deltaB`/`bnbpOnly`/etc. flag); if
`-outputFile <path>` is given, the debug file is `<path>.residuals.gpkg` instead. In
`-runFile` mode the debug filename is always derived from each run's own `outfile`
(`<outfile>.residuals.gpkg`), and for the ION_AUTO ionosphere-comparison mode only the
winning attempt's residuals are written. See `mosaicSource/CLAUDE.md` "Debug residual
GeoPackage output" for the shared writer implementation (also used by `tiepoints` and
`azparams`).

---

## Algorithm

### 1. Initialization

The geodat file is parsed to obtain SAR image geometry (near range, PRF, azimuth size,
pixel sizes, state vectors). The hemisphere is determined from tiepoint latitudes to
select the appropriate polar stereographic map projection.

### 2. Reading Range Offsets (`getROffsets`)

Range offsets are read from the binary offset file and bilinearly interpolated to each
tiepoint's range/azimuth image coordinates. Offsets are scaled by the range pixel size
and, if present, an ionospheric range-offset correction field is subtracted. Valid
tiepoint offsets are stored as pseudo-phase values in metres.

### 3. Velocity Correction (`addVelCorrections`)

Surface motion contributes a range displacement between the two acquisition epochs.
For each tiepoint, the horizontal velocity vector (vx, vy) is rotated from the
polar-stereographic grid into the radar line-of-sight / along-track frame using the
local satellite heading angle. The range-direction component (vyra) combined with the
vertical velocity (vz) gives a displacement:

$$
\Delta r_\text{motion} = \Delta t \left( v_\text{ra} \sin\psi - v_z \cos\psi \right)
$$

where $\psi$ is the local incidence angle and $\Delta t = N_\text{days} / 365.25$ years.
This displacement is subtracted from the observed range offset at each tiepoint,
isolating the purely geometric (baseline) signal.

### 4. State-Vector Baseline Initialisation (optional)

If two geodat files are provided in the offset `.dat` file, the program uses the
satellite state vectors from both acquisitions to compute initial estimates of the
normal (Bn) and parallel (Bp) baseline components and their linear along-track
derivatives (dBn, dBp) at the image start and end times.

### 5. Baseline Parameter Estimation (`computeRParams`)

The range offset at a tiepoint is related to the interferometric baseline through the
linearised geometric model (Joughin et al., *J. Glaciol.*, 1996, Eq. 7):

$$
\Delta r \approx -B_n \sin(\theta - \theta_c) - B_p \cos(\theta - \theta_c) + \frac{B^2}{2r} + c
$$

where $\theta$ is the local look angle, $\theta_c$ is the central look angle, $r$ is the slant
range, and $c$ is a constant range bias. Known non-linear terms are precomputed
from the current baseline estimate and subtracted from the observations, keeping the
inversion linear.

The along-track baseline variation is modelled as a polynomial in normalised azimuth
position $x \in [0, 1]$:

$$
B_n(x) = B_{n0} + \delta B_n \, x + \delta B_{nQ} \, x^2
$$
$$
B_p(x) = B_{p0} + \delta B_p \, x + \delta B_{pQ} \, x^2
$$

The system is solved by **Singular Value Decomposition (SVD)** least squares. Up to
six parameters are solved simultaneously depending on the selected mode:

| Mode (flag)        | Parameters solved           | Min tiepoints |
|--------------------|-----------------------------|---------------|
| Default            | Bn, dBn, dBp, const         | 4 |
| `-constOnly`       | const                        | 2 |
| `-bnbpdBpOnly`     | Bn, dBp, const              | 3 |
| `-bpdBpOnly`       | dBp, const                  | 2 |
| `-quadB`           | Bn, dBn, dBp, const, dBnQ, dBpQ | 6 |
| `-deltaBQ`         | Bn, dBn, dBp, dBnQ, dBpQ, Bp | 6 |
| `-deltaBC`         | Bp correction only          | 2 |

The fit iterates three times (k = 0..2). On iterations after the first, tiepoint
weights are scaled inversely by the residual standard deviation from the previous
iteration, providing a mild robustness weighting. Point weights read from the tiepoint
file (e.g., rock vs. ice) further modulate the per-point sigma, but are normalised so
the mean sigma tracks the observed residual level.

The final covariance matrix is computed via `svdvar` and reported as a 6×6 matrix
(mapped from the reduced parameter set to the canonical [Bn, dBn, dBp, dBnQ, dBpQ,
const/Bp] ordering).

---

## Dependencies

- `parseInputFile` / `initllToImageNew` — geodat parsing and image geometry
- `readTiePoints` / `computeTiePoints` — tiepoint I/O and image-coordinate projection
- `getBaselineFile` — state-vector baseline file reader
- `svdfit` / `svdvar` — SVD least-squares (Numerical Recipes)
- `bilinearInterp` — offset interpolation at tiepoint locations
- GDAL — registered at startup for any GDAL-format inputs
