# Azimuth Parameter Estimation: Model and Equations

## Overview

`computeAzParams` estimates azimuth offset calibration parameters by fitting a forward model to tie-point azimuth offset observations using SVD least squares. The goal is to calibrate the along-track component of SAR image offsets so that they can be converted to ice velocity.

---

## Geometry and Notation

| Symbol | Description |
|--------|-------------|
| `r0` | Slant range to tie point: `r0 = R_near + i * dr` |
| `θ` | Look angle at tie point: `θ = θ(r0, Re + z, ReH)` |
| `θC` | Center look angle |
| `Re` | Earth radius |
| `ReH` | Earth radius + satellite altitude |
| `z` | Elevation of tie point |
| `x` | Normalized along-track position: `x = azimuth_sample / (N_az * n_looks * slpA)` |
| `Bc`, `Bh` | Along-track (cross) and height baseline components |
| `dBc/ds`, `dBh/ds` | Rates of change of Bc and Bh along track (per single-look sample) |

The look angle is computed from slant range geometry:

```
θ = arccos( (ReH² + r0² - (Re+z)²) / (2 · ReH · r0) )
```

---

## Observation Equation

Each tie point provides a measured azimuth offset `φ_i` (in slant-range pixels, scaled to metres). The full model is:

```
φ_i = C0 + r0 · sin(θ) · (dBc/ds) - r0 · cos(θ) · (dBh/ds) + C1 · x_i
```

where:
- `C0` — constant azimuth offset (metres)
- `dBc/ds` — along-track rate of the cross-track baseline (rad/sample or equivalent)
- `dBh/ds` — along-track rate of the height baseline (fixed, not solved for)
- `C1` — linear along-track variation of the constant offset

The `dBh/ds` term is **never solved for** — it is always taken from orbit state vectors or an external file, and removed from the data before fitting:

```
y_i = φ_i - ( -r0 · cos(θ) · dBh/ds )
    = φ_i + r0 · cos(θ) · dBh/ds
```

If solving as a correction to an SV (state vector) predicted offset, the SV-predicted contribution is also subtracted:

```
y_i -= δ_SV(r_i, a_i)
```

---

## Baseline Rates

### From orbit state vectors (`initWithSV = TRUE`)

The cross and height baseline components are evaluated at the start and end of the image using two-epoch state vectors:

```
dBc/ds = ( Bc(t_end) - Bc(t_start) ) / ( N_az · n_looks · slpA )
dBh/ds = ( Bh(t_end) - Bh(t_start) ) / ( N_az · n_looks · slpA )
```

where `slpA` is the single-look azimuth pixel size in seconds (1/PRF).

### From external baseline file

```
dBc/ds = x2 / (PRF · slpA)
dBh/ds = x3 / (PRF · slpA)
```

where `x2`, `x3` are columns from the baseline parameter file.

If `deltaB != DELTABNONE`, both rates are forced to zero regardless of source.

---

## Fit Modes

Four fit modes are available depending on flags `constOnlyFlag` and `linFlag`:

### 1. Constant only (`constOnlyFlag=TRUE`, `linFlag=FALSE`) — 1 parameter

Solves for `C0` only. `dBc/ds` is held fixed from the orbit/file.

**Design matrix row:** `[1]`

**Solution:**
```
result = [ C0,  dBc/ds_fixed,  dBh/ds_fixed,  0 ]
```

### 2. Constant + linear along-track (`constOnlyFlag=TRUE`, `linFlag=TRUE`) — 2 parameters

Solves for `C0` and `C1`. `dBc/ds` is held fixed.

**Design matrix row:** `[1,  x_i]`

**Solution:**
```
result = [ C0,  dBc/ds_fixed,  dBh/ds_fixed,  C1 ]
```

### 3. Constant + dBc/ds (`constOnlyFlag=FALSE`, `linFlag=FALSE`) — 2 parameters

Solves for `C0` and `dBc/ds` simultaneously.

To condition the solution, `C0` is normalised by `AZCONST = 800000`:

**Design matrix row:** `[AZCONST,  r0 · sin(θ)]`

SVD solves for scaled parameters `a[1]`, `a[2]`:
```
C0       = a[1] · AZCONST
dBc/ds   = a[2]
```

**Solution:**
```
result = [ C0,  dBc/ds,  dBh/ds_fixed,  0 ]
```

### 4. Full: constant + dBc/ds + linear along-track (`constOnlyFlag=FALSE`, `linFlag=TRUE`) — 3 parameters

**Design matrix row:** `[AZCONST,  r0 · sin(θ),  x_i · AZCONST]`

SVD solves for `a[1]`, `a[2]`, `a[3]`:
```
C0       = a[1] · AZCONST
dBc/ds   = a[2]
C1       = a[3] · AZCONST
```

**Solution:**
```
result = [ C0,  dBc/ds,  dBh/ds_fixed,  C1 ]
```

---

## Iterative Weighting

The fit is run 3 times. After each iteration, the RMS residual `σ` is computed:

```
σ = sqrt( var(y_i - ŷ_i) )
```

Individual point weights are applied via:

```
sig_i = max( σ · w_i · N / Σw_j,  0.1 · σ )
```

where `w_i` is the tie-point weight and `N / Σw_j` renormalises so the mean sigma matches the data. The floor of `0.1 · σ` prevents the covariance matrix from being dominated by a few very low-weight points.

---

## Output

The four output parameters written to stdout are:

```
C0   dBc/ds   dBh/ds   C1
```

followed by a 4×4 covariance matrix (only the solved-for entries are non-zero; fixed parameters have zero covariance). The covariance is scaled by `azconst[i] · azconst[j]` to undo the conditioning normalisation.

---

## Left-Looking Geometry

The model equations require one correction for left-looking acquisitions, in how `dBc/ds` is computed.

### TCN coordinate frame

Baseline components are computed in the TCN (along-Track, Cross-track, Nadir) frame:
- **C** = N × V, cross-track — points **right** (e.g. east for a descending orbit) regardless of look direction
- **N** = −R/|R|, nadir (toward Earth)

For right-looking acquisitions the target is in the +C direction, so `bTCN[1]` (the C component) enters the model with the correct sign. For left-looking, the target is on the **−C** side, so the cross-track baseline contribution reverses sign.

This is exactly what `svBnBp` does when computing the perpendicular/normal baseline components:

```c
if (lookDir == LEFT)
    bTCN[1] = -bTCN[1];   // flip C component before projecting onto look direction
*bp = bTCN[2]*cos(theta) + bTCN[1]*sin(theta);
```

### Model validity for left-looking

The model formula is **correct as-is** for both look directions using raw TCN[1] components — no sign flip is needed in `computeBaselineRates`. Empirically: removing the measurement sign flips (see below) and using raw `dBc/ds` recovers the correct C0; adding a sign flip to `dBc/ds` for left-looking shifts C0 by `+2·r0·sin(θ)·dBc/ds`, which for typical imagery is hundreds of thousands of metres — clearly wrong.

The reason `svBnBp` flips bTCN[1] for left-looking is specific to how the perpendicular baseline is projected onto the slant-range look direction for interferometric phase. The azimuth offset model uses bTCN[1] in a different way (as a predictor variable in a least-squares fit), where the correct sign is naturally captured through the correlation structure of the fit itself.

### Measurement and application sign flips (incorrect workaround)

The original code contained heuristic sign flips (with a `????` comment) in two places:

- **`getOffsets.c`**: `if(inputImage.lookDir==LEFT) tiePoints->phase[i] *= -1;`
- **`interpOffsets.c`**: `if (inputImage->lookDir == LEFT) result *= -1.0;`

These are **unnecessary and incorrect**. Both are **commented out**. Azimuth offset measurements have the same sign convention (positive = forward along track) for both look directions, and the model equations need no modification.
