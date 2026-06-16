# geomosaic — SAR Backscatter Geocoder and Mosaicker

## Purpose

Geocodes and mosaics one or more SAR amplitude (power) images onto a polar
stereographic output grid using a DEM for terrain correction. Supports both
uncalibrated (cosmetic) and fully calibrated outputs, including Sentinel-1 σ₀ and
radiometrically terrain-corrected (RTC) γ₀.

### Layover Behaviour

geomosaic's RTC pipeline follows a **map-to-SLC** architecture: each output pixel is
independently geocoded back to radar coordinates and reads intensity from the
corresponding SLC range bin. In layover zones, where multiple terrain surfaces fold onto
the same range bin, several output map pixels each receive the full energy from that bin
and independently apply their own $A_\beta / A_\gamma$ terrain correction. The result
is that layover pixels appear several dB brighter than in products built by the
forward (SLC-to-map) approach — such as the NISAR GCOV product — which distributes bin
energy across all contributing map pixels naturally.

Because layover pixels represent geometrically ambiguous scattering with no recoverable
surface information regardless of how they are processed, this brightness excess is an
acceptable artefact of the architecture rather than a correctable error. Comparisons
with NISAR GCOV show good agreement over flat to moderately sloped terrain; differences
are confined to slopes steep enough to produce true geometric layover.

Where clean suppression of layover-affected pixels is required, the `-maskLayover` flag
(used in combination with any sub-pixel RTC mode) classifies an output pixel as layover
when more than half of its sub-pixel samples show a reversed Jacobian sign, and sets
those pixels to the $-30$ dB no-data level — consistent with the treatment of shadowed
pixels.

---

## Usage

```
geomosaic [options] inputFile demFile outFile
```

### Required Arguments

| Argument    | Description |
|-------------|-------------|
| `inputFile` | ASCII list of input images, geodat files, weights, and optional antenna pattern files |
| `demFile`   | DEM in XY polar stereographic format — provides elevations and projection |
| `outFile`   | Output file base name |

### Options

> **Note:** Some flags are experimental or may be obsolete. Options marked ⚠️ are
> rarely used or not fully supported in current builds.

| Option                        | Description |
|-------------------------------|-------------|
| `-fl <length>`                | Feathering length (m) at image edges |
| `-date1 MM-DD-YYYY`           | Start date for mosaic |
| `-date2 MM-DD-YYYY`           | End date for mosaic |
| `-descending`                 | Give descending-pass images priority in overlap regions |
| `-ascending`                  | Give ascending-pass images priority in overlap regions |
| `-nearestDate YYYY:MM:DD`     | ⚠️ Use only the image nearest to this date at each pixel |
| `-hybridZ <zthresh>`          | ⚠️ Apply `-nearestDate` only below elevation `zthresh` m |
| `-smoothL <n>`                | Multilook input images n×n before geocoding |
| `-smoothOut <n>`              | ⚠️ Oversample output (1, 2, or 3) — test only, use 1 in practice |
| `-removePad <n>`              | Zero-out first and last `n` columns of each input image |
| `-noData <val>`               | No-data value for input images |
| `-noPower`                    | Input is not power data (suppresses sin³ correction) |
| `-rsatFineCal`                | Apply RADARSAT fine-beam antenna pattern calibration |
| `-S1Cal`                      | Apply Sentinel-1 calibration to produce σ₀ and γ₀ (dB) |
| `-S1Psi`                      | Also output incidence angle map (`outFile.inc`) with `-S1Cal` |
| `-S1GammaCorr`                | Also output raw γ₀ terrain correction field (`outFile.gamcor`) |
| `-GTiff`                      | Write output as GeoTIFF |
| `-COG`                        | Write output as Cloud-Optimised GeoTIFF |
| `-subPixelRTC`                | Sub-pixel RTC: sample each output pixel on a uniform grid matched to the input resolution (see [Sub-pixel RTC](#sub-pixel-radiometric-terrain-correction)) |
| `-linearSubPixelRTC`          | Implies `-subPixelRTC`; replaces per-sub-pixel Newton solves with a first-order Jacobian approximation (~2× wall-time speedup; see [Linear Sub-pixel RTC](#linear-sub-pixel-rtc-linearsubpixelrtc)) |
| `-jacobianSubPixelRTC`        | Implies `-subPixelRTC`; weights each sub-pixel by $A_\beta \cdot |J|$ where $J$ is the map→radar Jacobian determinant — corrects foreshortening bias when reading SLC directly; suppresses γ₀ in pure layover (see [Jacobian Sub-pixel RTC](#jacobian-weighted-sub-pixel-rtc-jacobiansubpixelrtc)) |
| `-maskLayover`                | Used with any sub-pixel RTC mode: suppresses output pixels classified as layover (more than half of sub-pixels have reversed Jacobian sign) to $-30$ dB no-data, instead of passing through their elevated but unreliable brightness values |
| `-ompThreads N`               | Set the OpenMP thread count at runtime (default: 10, or `OMP_NUM_THREADS` if set in the shell) |
| `-byteScale`                  | ⚠️ Write output as scaled 8-bit byte image |
| `-BSlowerBound <val>`         | ⚠️ Lower bound for byte scaling [default: 0.53] |
| `-BSupperBound <val>`         | ⚠️ Upper bound for byte scaling [default: 2.4] |
| `-BSscale <val>`              | ⚠️ Scale factor for byte scaling [default: 0.001154] |
| `-BSexponent <val>`           | ⚠️ Exponent for byte scaling [default: 0.2] |

### Output Files

**Default (flat binary MSB float + `.geodat`)**:

| File               | Contents |
|--------------------|----------|
| `outFile`          | Geocoded mosaicked backscatter (linear power or dB) |
| `outFile.geodat`   | Grid geometry descriptor |
| `outFile.sigma0`   | σ₀ in dB (`-S1Cal`) |
| `outFile.gamma0`   | RTC γ₀ in dB (`-S1Cal`) |
| `outFile.inc`      | Incidence angle in degrees (`-S1Cal -S1Psi`) |
| `outFile.gamcor`   | Raw γ₀ terrain correction field in dB (`-S1Cal -S1GammaCorr`) |

**GeoTIFF/COG** (`-GTiff` / `-COG`): same names with `.tif` extension appended.

---

## Input File Format

ASCII file parsed by `processInputFileGeo`. Comments begin with `;`.

**Line 1 — Output grid geometry:**
```
x0  y0  xSize  ySize  deltaX  deltaY
```
Same format as `mosaic3d` (origin in km, size in km, pixel spacing in km).

**Line 2 — Number of input images:**
```
nFiles
```

**Per-image lines:**
```
imageFile  geodatFile  weight  [antPatFile]
```

| Field        | Description |
|--------------|-------------|
| `imageFile`  | SAR power image file (flat binary MSB float) |
| `geodatFile` | SAR image geometry parameter file |
| `weight`     | Image weight for feathered averaging (typically 1.0) |
| `antPatFile` | Optional antenna pattern file; use `poly` for RADARSAT fine-beam polynomial, `alos` for ALOS |

---

## Algorithm

### 1. Geocoding

For each output pixel at polar-stereographic coordinates $(x, y)$:

1. Convert $(x, y)$ to lat/lon using the map projection.
2. Look up the surface elevation $h$ from the DEM.
3. Project (lat, lon, $h$) to SAR image range/azimuth coordinates via the geodat
   state vectors.
4. Bilinearly interpolate the input power image at the computed range/azimuth.

### 2. Mosaicking

Overlapping images are combined by feather-weighted averaging. The feather weight at
each pixel decays with distance from the image edge over a length `-fl`. The final
mosaic value is the weighted sum divided by the accumulated weight. Alternatively,
`-ascending` / `-descending` priority or `-nearestDate` selection can be used instead
of averaging.

---

## Radiometric Calibration

### Uncalibrated Mode (default)

No antenna pattern file, no `-S1Cal`, no `-rsatFineCal`. An empirical incidence-angle
correction is applied to improve cosmetic appearance across the swath:

$$
\text{DN}_\text{out} = \text{DN}_\text{in} \cdot \sin^3(\psi)
$$

where $\psi$ is the local incidence angle. The first factor of $\sin\psi$ approximates
the $\sigma^0$ geometric normalisation; the remaining $\sin^2\psi$ provides an
empirical across-track shading correction. **This is not a true radiometric calibration.**

### RADARSAT Fine-Beam (`-rsatFineCal`)

The raw power is divided by a 9th-degree polynomial antenna pattern $G(\theta)$
fitted to the ASF-provided gain table (referenced to look angle $\theta \approx 33.4°$):

$$
\sigma^0_\text{linear} = \frac{\text{DN}_\text{in}}{G(\theta)}
$$

A final linear calibration is then applied:

$$
\sigma^0_\text{cal} = \frac{\sigma^0_\text{linear}}{27.3} - 0.0058
$$

followed by conversion to dB:

$$
\sigma^0\,[\text{dB}] = 10 \log_{10}\!\left(\sigma^0_\text{cal}\right)
$$

### Sentinel-1 Calibrated (`-S1Cal`)

#### σ₀

The Sentinel-1 calibration constant $\beta_0$ is read from a file named `betaNought`
in the same directory as the input image. The calibrated sigma-nought is:

$$
\sigma^0 = \frac{|\text{DN}|^2 \cdot \sin\psi}{\beta_0^2}
$$

where $\psi$ is the local incidence angle. Converted to dB:

$$
\sigma^0\,[\text{dB}] = 10 \log_{10}\!\left(\sigma^0\right)
$$

Values below $-29.9$ dB are clipped to a no-data value of $-30$ dB.

#### γ₀ — Radiometric Terrain Correction (RTC)

Gamma-nought accounts for the actual terrain geometry rather than the ideal flat-Earth
incidence angle. For each output pixel, the 3D ECEF coordinates of its four corners
are projected onto two reference planes:

- **Beta plane** ($A_\beta$): the plane containing the satellite velocity vector and
  the look vector (the range-Doppler plane). This represents the area as seen by the
  radar in the slant-range azimuth geometry.
- **Gamma plane** ($A_\gamma$): the plane perpendicular to the look vector. This
  represents the area contributing to the backscatter coefficient.

The γ₀ terrain correction stored per pixel is (in dB):

$$
\Delta\gamma\,[\text{dB}] = 10\log_{10}\!\left(\frac{A_\beta / A_\gamma}{\sin\psi}\right)
$$

The final RTC gamma-nought is:

$$
\gamma^0\,[\text{dB}] = \sigma^0\,[\text{dB}] + \Delta\gamma\,[\text{dB}]
$$

In linear units this is equivalent to:

$$
\gamma^0 = \sigma^0 \cdot \frac{A_\beta}{A_\gamma \cdot \sin\psi}
           = \frac{|\text{DN}|^2}{\beta_0^2} \cdot \frac{A_\beta}{A_\gamma}
$$

This matches the standard RTC formula of Small (2011, *IEEE TGRS*):

$$
\gamma^0_T = \beta^0 \cdot \frac{A_\beta}{A_\gamma}
$$

where $\beta^0 = |\text{DN}|^2 / \beta_0^2$ is the beta-nought intensity.
**The implementation is correct.**

For shadowed pixels (surface normal facing away from the radar, detected by a
negative dot product of the surface normal with the look vector), the correction is
set to $-30$ dB (no-data).

The area ratio $A_\beta / A_\gamma$ is computed by projecting the ECEF corner
positions of a neighbourhood of output pixels onto the two planes and summing the
3D triangle areas. For bright pixels near flat terrain, surrounding pixels are
aggregated to account for multiple ground samples contributing to the same range bin.

---

## Sub-pixel Radiometric Terrain Correction

### Motivation

The standard RTC pipeline (`-S1Cal` without `-subPixelRTC`) samples the input image
once at the centre of each output pixel and computes the beta/gamma area ratio
$A_\beta / A_\gamma$ for that single point. This is accurate when the output and input
pixel sizes are similar, but when the output grid is coarser than the input resolution
the output pixel footprint spans several input pixels. A layover or shadow zone that
covers part of the output footprint is then either counted in full or missed entirely,
depending on whether the centre sample falls inside it. The result is systematic
over- or under-correction of $\gamma^0$ at resolution transitions.

Sub-pixel RTC solves this by sampling every output pixel on a uniform grid, accumulating
both power and area consistently from the same set of sub-pixel positions.

---

### Standard RTC (single-sample)

For output pixel centred at $(x, y)$ in polar-stereographic km coordinates:

1. $(x, y) \xrightarrow{\texttt{xytoll1}} (\phi, \lambda)$
2. $(\phi, \lambda) \xrightarrow{\text{DEM}} h$
3. $(\phi, \lambda, h_\text{WGS}) \xrightarrow{\texttt{llToImageNew}} (r_0, az_0)$  — one Newton solve
4. $p_0 = \texttt{bilinearInterp}(r_0, az_0)$
5. $\psi_0, p_0 \leftarrow \texttt{applyCorrections}(r_0, az_0, h)$
6. $(A_\beta, A_\gamma) \leftarrow \texttt{AbAg}(x, y, az_0)$  — 4 DEM lookups + orbital normal

$$
\gamma^0 = \frac{p_0}{\beta_0^2} \cdot \frac{A_\beta}{A_\gamma}
$$

```
  ┌─────────────────────┐
  │                     │
  │          ×          │   × = single sample at pixel centre
  │                     │
  └─────────────────────┘
```

---

### Sub-pixel RTC (`-subPixelRTC`)

#### Sub-pixel grid

The output pixel footprint (size $\Delta X_\text{out} \times \Delta Y_\text{out}$) is
subdivided into an $n_r \times n_a$ grid where the sub-pixel spacing matches the input
slant-range and azimuth resolution:

$$
n_r = \max\!\left(1,\; \text{round}\!\left(\frac{\Delta X_\text{out}}{\Delta r_\text{in}}\right)\right),
\qquad
n_a = \max\!\left(1,\; \text{round}\!\left(\frac{\Delta Y_\text{out}}{\Delta az_\text{in}}\right)\right)
$$

Sub-pixel spacings in km:

$$
\delta x = \frac{\Delta X_\text{out} \cdot \kappa}{n_r}, \qquad
\delta y = \frac{\Delta Y_\text{out} \cdot \kappa}{n_a}
$$

where $\kappa = 10^{-3}$ converts metres to km. Sub-pixel centres:

$$
x_k = x + \left(k - \frac{n_r - 1}{2}\right)\delta x, \quad k = 0, \ldots, n_r - 1
$$

$$
y_l = y + \left(l - \frac{n_a - 1}{2}\right)\delta y, \quad l = 0, \ldots, n_a - 1
$$

```
  ┌─────────────────────────────┐
  │  ×    ×    ×    ×    ×    × │
  │                             │
  │  ×    ×    ×    ×    ×    × │    nr = 6, na = 4
  │                             │    × = sub-pixel sample
  │  ×    ×    ×    ×    ×    × │
  │                             │
  │  ×    ×    ×    ×    ×    × │
  └─────────────────────────────┘
  δx ←──→                       δy ↕
```

#### Per-sub-pixel computation

For each sub-pixel $(k, l)$:

1. $(x_{kl}, y_{kl}) \xrightarrow{\texttt{xytoll1}} (\phi_{kl}, \lambda_{kl})$
2. $(\phi_{kl}, \lambda_{kl}) \xrightarrow{\text{DEM}} h_{kl}$
3. $(\phi_{kl}, \lambda_{kl}, h_{kl}) \xrightarrow{\texttt{llToImageNew}} (r_{kl}, az_{kl})$  — full Newton solve
4. $p_{kl} = \texttt{bilinearInterp}(r_{kl}, az_{kl})$
5. $\psi_{kl}, p_{kl} \leftarrow \texttt{applyCorrections}(r_{kl}, az_{kl}, h_{kl})$
6. $(A_{\beta,kl}, A_{\gamma,kl}) \leftarrow \texttt{AbAg}(x_{kl}, y_{kl}, az_{kl}, \texttt{subOutput})$

where `subOutput` is a copy of the output image structure with $\Delta X$ and $\Delta Y$
scaled by $1/n_r$ and $1/n_a$ respectively, so that `AbAg` computes the area for the
sub-pixel footprint rather than the full output pixel.

Sub-pixels for which $p_{kl} \leq 0$ or that fall in shadow ($A_\beta < -0.001$) are
excluded.

#### Accumulation

Power is accumulated with $A_\beta$-weighting, which is consistent with the
beta-nought reference used by `applyCorrections`:

$$
\bar{p} = \frac{\displaystyle\sum_{k,l} p_{kl}\, A_{\beta,kl}}
               {\displaystyle\sum_{k,l} A_{\beta,kl}}
$$

$$
A_\beta = \sum_{k,l} A_{\beta,kl}, \qquad
A_\gamma = \sum_{k,l} A_{\gamma,kl}, \qquad
\bar\psi = \frac{1}{N_\text{valid}}\sum_{k,l} \psi_{kl}
$$

The final $\gamma^0$ is assembled from these accumulated quantities exactly as in the
single-sample case:

$$
\gamma^0 = \frac{\bar{p}}{\beta_0^2} \cdot \frac{A_\beta}{A_\gamma}
$$

#### Cost

Each output pixel requires $n_r \times n_a$ Newton solves and $4 \times n_r \times n_a$
DEM lookups (4 ECEF corner lookups per `AbAg` call). For a 20 m output grid over
NISAR L-band SLC input with ~3 m slant-range and ~6 m azimuth spacing,
$n_r \approx 7$, $n_a \approx 3$, giving ~21 Newton solves and ~84 DEM lookups per
output pixel. This is accurate but expensive; see
[Linear Sub-pixel RTC](#linear-sub-pixel-rtc-linearsubpixelrtc) for the fast path.

---

### Linear Sub-pixel RTC (`-linearSubPixelRTC`)

#### Motivation

The mapping $(x, y) \mapsto (r, az)$ is the composition of a smooth map projection
inversion, a DEM lookup, and an orbital Newton solve. Over the 20 m extent of a
single output pixel the mapping is very nearly affine; a first-order Taylor expansion
is indistinguishable from the full solve at sub-pixel precision.

#### Jacobian computation

Three anchor points are solved with the full Newton pipeline — the pixel centre and
two offset points one full output pixel away in each map direction:

$$
P_0 = (x,\; y), \qquad
P_x = (x + \Delta X \cdot \kappa,\; y), \qquad
P_y = (x,\; y + \Delta Y \cdot \kappa)
$$

```
        y + ΔY·κ
           ×  P_y
           │
    ΔY·κ  │
           │
  P_0 ×───────────────× P_x
  (x, y)      ΔX·κ    (x + ΔX·κ, y)
```

Each anchor is resolved through the full chain
$\texttt{xytoll1} \to \texttt{getXYHeight} \to \texttt{sphericalToWGSElev} \to \texttt{llToImageNew}$,
giving $(r_0, az_0)$, $(r_x, az_x)$, $(r_y, az_y)$.

The first-order Jacobian is estimated by finite difference:

$$
\frac{\partial r}{\partial x} \approx \frac{r_x - r_0}{\Delta X \cdot \kappa}, \qquad
\frac{\partial r}{\partial y} \approx \frac{r_y - r_0}{\Delta Y \cdot \kappa}
$$

$$
\frac{\partial az}{\partial x} \approx \frac{az_x - az_0}{\Delta X \cdot \kappa}, \qquad
\frac{\partial az}{\partial y} \approx \frac{az_y - az_0}{\Delta Y \cdot \kappa}
$$

#### Sub-pixel approximation

For sub-pixel $(k, l)$ at displacement $(\delta x_k, \delta y_l)$ from the pixel centre:

$$
r_{kl} \approx r_0 + \frac{\partial r}{\partial x}\,\delta x_k + \frac{\partial r}{\partial y}\,\delta y_l
$$

$$
az_{kl} \approx az_0 + \frac{\partial az}{\partial x}\,\delta x_k + \frac{\partial az}{\partial y}\,\delta y_l
$$

No Newton solve is needed for any sub-pixel. The height $h_0$ from the centre anchor is
reused for all sub-pixels; height variation across a single output pixel causes
negligible error in the depression angle correction.

#### Normal-vector recycling

The look-direction normals $(\hat{n}_l)$ and range-plane normals $(\hat{n}_r)$ computed
inside `AbAg` via `satNorm` / `rangePlane` depend on the azimuth time. This
varies by at most a few tenths of a millisecond across a single output pixel —
far below the scale at which the orbital geometry changes. The normals are therefore
computed once on the first valid sub-pixel (`recycle = FALSE`) and reused for all
subsequent sub-pixels (`recycle = TRUE`), eliminating $n_r n_a - 1$ orbital
interpolations per output pixel.

#### Cost comparison

| Path | Newton solves | `satNorm` calls | DEM lookups (AbAg) |
|---|---|---|---|
| Standard RTC | 1 | 1 | 4 |
| `subPixelRTC` | $n_r n_a$ | $n_r n_a$ | $4\,n_r n_a$ |
| `linearSubPixelRTC` | **3** | **1** | $4\,n_r n_a$ |

The dominant remaining cost with `-linearSubPixelRTC` is the $4 n_r n_a$ DEM corner
lookups inside `AbAg`. For a 20 m output grid with $n_r \approx 7$, $n_a \approx 3$
this is ~84 lookups per output pixel, compared to 3 Newton solves — the DEM accesses
are the primary bottleneck.

#### Accuracy

The linear approximation introduces a truncation error of order $O(\delta x^2)$ in
range and azimuth. For a 20 m output pixel with $n_r = 7$, $n_a = 3$
($\delta x \approx 3$ m, $\delta y \approx 7$ m) and a typical along-track range-rate
of ~$6 \times 10^{-3}$ pixels/m, the range error at the far corner of the pixel is:

$$
\varepsilon_r \sim \frac{1}{2}\,\frac{\partial^2 r}{\partial x^2}\,(\delta x)^2 \ll 0.01\text{ pixels}
$$

This is well below the bilinear interpolation floor and negligible in the final $\gamma^0$.

Observed wall-time reduction: **~53%** (26.5 min → 12.5 min) on a representative
NISAR mosaic with a 20 m output grid over Greenland.

---

### Jacobian-Weighted Sub-pixel RTC (`-jacobianSubPixelRTC`)

#### Motivation

`-subPixelRTC` and `-linearSubPixelRTC` accumulate power and area with uniform
$A_\beta$-weighting: every sub-pixel that passes the shadow and data-quality checks
contributes equally to $\sum A_\beta$. This is correct when the output grid is coarser
than the input (pre-multilooked) product because each sub-pixel then represents a
distinct ground patch.

When geomosaic reads an **SLC directly** (no pre-multilooking), the sub-pixel
accumulation *is* the multilooking step. In foreshortened terrain multiple map-space
sub-pixels $(x_k, y_l)$ project to the same or adjacent SLC range bins. Without
Jacobian weighting those duplicate samples are counted equally, artificially inflating
$\sum A_\beta$ and biasing the area ratio $A_\beta / A_\gamma$.

The map→radar Jacobian determinant

$$
J_{kl} = \frac{\partial r}{\partial x}\frac{\partial az}{\partial y}
        - \frac{\partial r}{\partial y}\frac{\partial az}{\partial x}
$$

measures how much unique radar-space area each map-space sub-pixel represents:

- **Normal terrain:** $|J| \approx 1$ (one map sub-pixel → one radar pixel).
- **Foreshortening:** many map sub-pixels share the same SLC bin; each has small $|J|$;
  their sum $\sum |J| \approx 1$ unique radar-pixel worth of area.
- **Layover:** $J$ changes sign relative to normal terrain (the range ordering of the
  terrain reverses). Heavy layover gives $J$ with the *opposite* sign to normal terrain
  for the majority of sub-pixels.

Weighting by $A_\beta \cdot |J|$ makes the accumulation equivalent to GCOV's
area-weighted SLC multilooking in radar geometry.

#### Algorithm

**Pass 1 — build the range/azimuth grid.**

Full Newton solves at every sub-pixel position store the radar coordinates in
three arrays:

$$
\text{rGrid}[k,l] = r_{kl}, \quad
\text{azGrid}[k,l] = az_{kl}, \quad
\text{hGrid}[k,l] = h_{kl}
$$

**Reference sign.**

Before the accumulation loop, the Jacobian is evaluated at the grid centre
$(k_c, l_c) = (\lfloor n_r/2 \rfloor,\, \lfloor n_a/2 \rfloor)$ using central
differences on the pre-filled grids. Its sign is stored as `refSign` ($\pm 1$). This
step is orbit-independent: `refSign` is $-1$ for descending right-looking passes
(where $\partial az/\partial y < 0$ in GrIMP's northward-$y$ convention) and $+1$ for
ascending passes. Using the centre pixel avoids hard-coding either convention.

**Pass 2 — accumulate with $A_\beta \cdot |J|$ weighting.**

For each sub-pixel $(k, l)$, the Jacobian is estimated by central finite differences
of the pre-filled grids (forward/backward differences at the grid edges):

$$
\frac{\partial r}{\partial x}\bigg|_{kl}
  \approx \frac{\text{rGrid}[k{+}1,l] - \text{rGrid}[k{-}1,l]}{(2)\,\delta x}
$$

$$
\frac{\partial az}{\partial y}\bigg|_{kl}
  \approx \frac{\text{azGrid}[k,l{+}1] - \text{azGrid}[k,l{-}1]}{(2)\,\delta y}
$$

(and similarly for the off-diagonal terms). The composite weight is

$$
w_{kl} = A_{\beta,kl} \cdot |J_{kl}|
$$

and the accumulation becomes

$$
\bar{p} = \frac{\displaystyle\sum_{k,l} p_{kl}\, w_{kl}}
               {\displaystyle\sum_{k,l} w_{kl}}, \qquad
A_\beta = \sum_{k,l} w_{kl}, \qquad
A_\gamma = \sum_{k,l} A_{\gamma,kl}\,|J_{kl}|
$$

$$
\bar\psi = \frac{\displaystyle\sum_{k,l} \psi_{kl}\,|J_{kl}|}
                {\displaystyle\sum_{k,l} |J_{kl}|}
$$

The normal-vector recycling optimisation from `-linearSubPixelRTC` is also applied:
`AbAg` is called with `recycle = FALSE` only on the first valid sub-pixel and
`recycle = TRUE` for all subsequent ones.

#### Layover treatment

A sub-pixel is in **layover** if $J_{kl}$ has the *opposite* sign to `refSign`.
If more than half of the valid sub-pixels are in layover the pixel is classified as a
pure-layover zone: the function returns `TRUE` (shadow-equivalent), suppressing the
$\gamma^0$ terrain correction for that output pixel. The $\sigma^0$ path is unaffected
(the shadow flag in the caller only controls the $\gamma^0$ correction buffer).

For **mixed pixels** (foreshortening edge — some sub-pixels in layover, some not)
the $|J|$ weighting handles the transition continuously and the function returns
`FALSE` with the accumulated values.

```
  Map space                     Radar space
  ─────────────────────        ────────────────────
  × × × × ×  ← slope tip       ╔══╗
  × × × × ×  ← mid-slope  →    ║  ║  ← all collapsed into ~1 bin
  × × × × ×  ← slope base      ╚══╝
             each |J| ≈ 0.2     Σ|J| ≈ 1.0
```

#### Cost

Pass 1 performs the same $n_r \times n_a$ Newton solves as `-subPixelRTC`. Pass 2 adds
only array reads and arithmetic for the finite-difference Jacobian (negligible), plus
one additional malloc/free for the three $n_r \times n_a$ grids. The `AbAg` recycling
reduces the full orbital-geometry computation to a single call per output pixel,
comparable to `-linearSubPixelRTC`.

| Path | Newton solves | Full `AbAg` calls | Recycled `AbAg` |
|---|---|---|---|
| `subPixelRTC` | $n_r n_a$ | $n_r n_a$ | 0 |
| `linearSubPixelRTC` | 3 | 1 | $n_r n_a - 1$ |
| `jacobianSubPixelRTC` | $n_r n_a$ | 1 | $n_r n_a - 1$ |

**When to use:** `-jacobianSubPixelRTC` when reading SLC input directly and accurate
foreshortening/layover treatment matters. `-linearSubPixelRTC` when speed is the
priority and the input is already multilooked to approximately the output resolution.

---

## Parallelism

The inner geocoding loop over output pixels is parallelised with OpenMP
(`#pragma omp parallel for schedule(dynamic, 8)`). Each thread works on a contiguous
block of rows; `llToImageNew`'s warm-start state is kept thread-private by giving each
thread its own copy of the `inputImageStructure`.

The thread count is resolved in this priority order:

1. `-ompThreads N` command-line flag (highest priority)
2. `OMP_NUM_THREADS` environment variable
3. Built-in default of **10 threads**

---

## Dependencies

- `parseInputFile` / `initllToImageNew` — geodat parsing and SAR image geometry
- `readXYDEM` — DEM I/O and projection
- `llToImageNew` — lat/lon to SAR range/azimuth projection
- `getXYHeight` — DEM elevation lookup
- `computeScale` — feathering weight computation
- `outputGeocodedImage` / `outputGeocodedImageTiff` — binary and GeoTIFF output
- `subPixelGammaRTC` (`subpixelRTC.c`) — full sub-pixel RTC accumulation loop
- `subPixelGammaRTCLinear` (`linearSubpixelRTC.c`) — Jacobian-approximated sub-pixel RTC
- `subPixelGammaRTCJacobian` (`jacobianSubpixelRTC.c`) — $|J|$-weighted sub-pixel RTC with layover detection
- `AbAg` (`makeGeoMosaic.c`) — beta/gamma projected area computation
- GDAL — GeoTIFF/COG output
