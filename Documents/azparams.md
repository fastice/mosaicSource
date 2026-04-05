# azparams — Azimuth Parameter Estimator

## Purpose

Calibrates azimuth offset measurements derived from speckle tracking by estimating
along-track baseline parameters (constant offset, along-track baseline rate, and
optionally a linear along-track trend) using tiepoints with known locations and
surface velocities. Output is written to stdout in a format suitable for downstream
mosaicking tools.

---

## Usage

```
azparams [options] -nDays nDays geodatFile tiepointsFile offsetFile baselineFile
```

### Required Arguments (positional, in order)

| Argument        | Description |
|-----------------|-------------|
| `geodatFile`    | SAR image geometry parameter file |
| `tiepointsFile` | Tiepoint file with lat, lon, elevation, and velocity (vx, vy, vz) |
| `offsetFile`    | Binary azimuth offset file (a companion `.dat` metadata file must also exist) |
| `baselineFile`  | CW state-vector baseline file |

### Options

| Option          | Description |
|-----------------|-------------|
| `-nDays nDays`  | Temporal baseline in days **(required)** |
| `-constOnly`    | Estimate only the constant offset term (do not solve for baseline-dependent terms) |
| `-linear`       | Add a linear along-track trend to either the `constOnly` or baseline parameter solution |
| `-useSV`        | Estimate a correction after removing a state-vector-determined azimuth offset |
| `-quiet`        | Suppress tiepoint echo to stdout |

> **Note:** `-nDays` is required; the program will error if it is omitted.

### Output

Results are written to **stdout** in a commented format:

- Tiepoint locations and extracted azimuth offsets (with `;` prefix)
- Number of tiepoints used / given
- Fit residual sigma (`sigma*sqrt(X2/n)`)
- 4×4 parameter covariance matrix
- Final estimated parameters on one line: `const  dBc/ds  dBh/ds  linConst`

---

## Algorithm

### 1. Initialization

The geodat file is parsed to obtain SAR image geometry (near range, PRF, azimuth pixel
size, state vectors). The hemisphere is determined from tiepoint latitudes to select
the appropriate polar stereographic map projection.

### 2. Reading Azimuth Offsets (`getOffsets`)

Azimuth offsets are read from the binary offset file and bilinearly interpolated to
each tiepoint's range/azimuth image coordinates. Valid offsets are scaled by the
azimuth pixel size to convert from pixels to metres and stored as pseudo-phase values.

### 3. Along-Track Baseline Rates

The azimuth offset at a tiepoint depends on the cross-track ($B_c$) and height ($B_h$)
components of the interferometric baseline and their along-track rates of change. These
rates are obtained in one of two ways:

- **From the baseline file** (default): $dB_c/ds$ and $dB_h/ds$ are read directly from
  the state-vector baseline file and normalised by PRF and single-look azimuth pixel
  size.
- **From state vectors** (`-useSV`): rates are computed from the TCN baseline evaluated
  at the image start and end times using both acquisitions' state vectors:

$$
\frac{dB_c}{ds} = \frac{B_c(t_2) - B_c(t_1)}{N_\text{az} \cdot \delta s}, \qquad
\frac{dB_h}{ds} = \frac{B_h(t_2) - B_h(t_1)}{N_\text{az} \cdot \delta s}
$$

where $N_\text{az}$ is the number of single-look azimuth pixels and $\delta s$ is the
single-look azimuth pixel size.

### 4. Offset Corrections (`addOffsetCorrections`)

Previously removed baseline contributions are added back to the tiepoint observations
before the inversion, ensuring the fit is performed against the full observed offset.

### 5. Azimuth Parameter Estimation (`computeAzParams`)

The azimuth offset at a tiepoint is related to the along-track baseline geometry by
(analogous to Joughin et al., *J. Glaciol.*, 1996):

$$
\Delta a = c + r \sin\theta \cdot \frac{dB_c}{ds} - r \cos\theta \cdot \frac{dB_h}{ds} + \ell \cdot x
$$

where $r$ is slant range, $\theta$ is the local look angle, $c$ is a constant azimuth
bias, $\ell$ is an optional linear along-track trend coefficient, and $x$ is the
normalised along-track position. Since $dB_h/ds$ is never solved for directly, its
contribution is precomputed and subtracted from the observations before the inversion.

The system is solved by **Singular Value Decomposition (SVD)** least squares. The
number of parameters depends on the selected mode:

| Mode                        | Parameters solved             | Min tiepoints |
|-----------------------------|-------------------------------|---------------|
| Default                     | $c$, $dB_c/ds$               | 2 |
| `-constOnly`                | $c$ only                     | 1 |
| `-linear`                   | $c$, $dB_c/ds$, $\ell$       | 3 |
| `-constOnly -linear`        | $c$, $\ell$                  | 2 |

The fit iterates three times. On each iteration, tiepoint sigmas are updated using the
residual standard deviation from the previous iteration, providing mild robustness
weighting. Per-point weights (e.g., rock vs. ice) further modulate the sigmas but are
normalised to keep the mean sigma consistent with the observed residual level (minimum
sigma floor of $0.1 \cdot \sigma_P$).

The final 4×4 covariance matrix is computed via `svdvar` and reported in the canonical
parameter ordering: $[c,\ dB_c/ds,\ dB_h/ds,\ \ell]$.

---

## Dependencies

- `parseInputFile` / `initllToImageNew` — geodat parsing and image geometry
- `readTiePoints` / `computeTiePoints` — tiepoint I/O and image-coordinate projection
- `getBaselineFile` / `getBaselineRates` — state-vector baseline file reader
- `svdfit` / `svdvar` — SVD least-squares (Numerical Recipes)
- `bilinearInterp` — offset interpolation at tiepoint locations
- `svInitAzParams` / `svAzOffset` — state-vector azimuth offset model (`-useSV` mode)
- GDAL — registered at startup for any GDAL-format inputs
