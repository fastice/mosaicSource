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
| `-fl <length>`                | Feathering length in **output pixels** (not metres), tapering each input's weight linearly to its data edge; e.g. `-fl 100` at 100 m = 10 km. Allowed with `-S1Cal` since 2026-09-23 (was refused), not with `-ascending/-descending` |
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
| `-nearRange` / `-farRange`    | Keep only inputs whose ellipsoidal incidence angle is within `-angleTolerance` of the per-pixel minimum (near) or maximum (far). Mutually exclusive; not allowed with `-min`/`-max`, `-ascending`/`-descending` or `-nearestDate` (see [Range selection](#near-range--far-range-selection-nearrange--farrange)) |
| `-angleTolerance <deg>`       | Tolerance for the above [1.0] |
| `-angleStride <n>`            | Output pixels per cell in the selection pre-pass [10] |
| `-gcov <file.yaml>`           | Also mosaic already-geocoded NISAR GCOV HDF5 products listed in the YAML file (see [NISAR GCOV Inputs](#nisar-gcov-inputs-gcov)). With GCOVs, `inputFile` may list 0 range/Doppler images |
| `-calOutput sigma0\|gamma0\|both` | With `-S1Cal`: write only `.sigma0`, only `.gamma0`, or both [both, the previous behaviour]. Without `-S1Cal` it only selects the quantity GCOVs contribute [gamma0] |
| `-int16`                     | With `-S1Cal` and `-GTiff`/`-COG`: write `.sigma0`/`.gamma0` as Int16 round(dB×100), with scale 0.01 and nodata −3000 in the file. Lossless and about 2× smaller as a COG (a predictor is used) |
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

`nFiles` may be 0 when `-gcov` supplies the inputs; the grid line is still required.

---

## Near-range / far-range selection (`-nearRange` / `-farRange`)

Normally every input covering a pixel is averaged, so the mosaic mixes near- and far-range looks
and its viewing geometry follows the track layout rather than the ice. These flags keep, at each
output pixel, only the inputs whose incidence angle is within `-angleTolerance` of the **minimum**
(`-nearRange`) or **maximum** (`-farRange`) over all inputs there. Works for range/Doppler and
GCOV inputs, together or separately. Without either flag the output is byte-identical to before.

### What the tolerance actually does

Across a NISAR swath the incidence angle spans ~20 deg and is strictly monotone, so **cross-track
overlaps differ by 20-26 deg** (median 9.1 deg over 1358 measured overlapping pairs) while
**adjacent frames of the same track differ by exactly 0**. A 1 deg tolerance therefore rejects the
neighbouring track outright; its real job is to keep **repeat passes of the same track and
along-track frame neighbours**, which then average and feather together as usual.

The consequence is that cross-track seams become deliberate hard boundaries. Feathering cannot
soften them, because the selection removes the overlap that feathering needs — on either side of
the seam a different track is the only contributor, and the weighted average normalises its taper
away. That step is physical (different geometry), not an artifact, and widening the tolerance
moves the seam rather than removing it.

**But the step is small, and shrinks with overlap — as does the whole point of the mode.** With
swath width `W`, track spacing `D` and incidence spanning `Δ` across a swath, each track is used
only over the strip where it is nearest (or farthest), of width `D`. So the output incidence
sawtooths over `Δ·D/W`, and the seam step is that *same* quantity — not `Δ`. Heavy overlap shrinks
both together; at 50% overlap you get half of each.

Measured on the PIG box below, against a ~20 deg single-frame swath span:

| | incidence span (p1-p99) | seam step p99.9 | pixels with a step > 1 deg |
|---|---|---|---|
| plain | 9.58 deg | 0.01 deg | 0.009% |
| `-nearRange` | 6.17 deg | 3.42 deg | 0.230% |
| `-farRange` | 2.77 deg | 5.91 deg | 0.331% |

so `Δ·D/W` ≈ 6.2 deg implies `D/W` ≈ 0.31, i.e. ~69% overlap at this latitude. Two notes:

- `-farRange` gives the *narrower* span, because `dψ/dx` flattens with range: a strip of width `D`
  at far range covers fewer degrees than the same strip at near range. The trade is the far-range
  end of the SNR and resolution range.
- `D/W` is set by the orbit and latitude, not by how much data accumulates. Tracks converge toward
  the pole, so the mode tightens up exactly over the Antarctic interior; extra cycles add repeats
  at the same geometry rather than narrowing the span further.

### Why two passes

The obvious single pass — track the running extremum and test each input against it as it
arrives — is order-dependent. A pixel seen by track 4 at 46 deg and track 5 at 35 deg, with files
in name order, has extremum 46 when track 4 is tested, so track 4 is included and averaged in;
reverse the order and it is rejected. So pass 1 establishes the extremum over every input before
any accumulation, and pass 2 mosaics and filters.

### Pass 1 is cheap

Pass 1 reads no image data. Incidence varies ~0.05 deg/km and is monotone, so it is evaluated on a
coarse grid (`-angleStride`, default 10 output pixels) — 0.03 deg of quantisation at 100 m against
a 1 deg tolerance. Measured: stride 5 vs 10 moves 0.6% of pixels, all at selection boundaries.

- Range/Doppler: geometry only (`xytoll1` -> `getXYHeight` -> `llToImageNew` -> `psiRReZReH`), about
  1/stride^2 of a geocoding pass. `llToImageNew`'s warm start (`lastTime`) is saved and restored so
  pass 2 solves exactly as it would alone — verified byte-identical.
- GCOV: the `metadata/radarGrid/incidenceAngle` cube (11 MB, 0.4 s per frame, against 15 s for
  HHHH) plus the `mask` layer for the footprint. The footprint test is essential: the cube is valid
  over the **whole rectangle**, including the ~39% of it that is no-data, and since psi is monotone
  the extremum over the rectangle sits at a corner outside the swath.

### It cannot leave a hole

Pass 1 sees geometry, not data, so it can pick an input that turns out to have no valid pixel
there — GCOV NaN corners, `bilinearInterp` rejecting a pixel whose 4 neighbours are not all valid,
`smoothImage`, `removePad`, interior dropouts. Rather than trying to predict all of that, pass 2
accumulates a **second, unfiltered mosaic**, and any pixel the filter empties falls back to it,
carrying the matching incidence and gamma correction across. Coverage is therefore identical to the
plain mosaic by construction, which is checked in testing.

The fallback rate is logged per run (`rangeSelect: N of M valid pixels ... fell back`). It should
be a fraction of a percent; a large or structured count means the pre-pass footprint is wrong. This
is also why pass 1 applies **exactly** pass 2's mask test: an earlier version accepted `mask != 255`
instead of `mask == 1`, which over-claims coverage and drove the far-range fallback to 11.6%
against 0.3% once matched.

### Incidence angle used

The **ellipsoidal** angle, not the local (terrain) one: `psiRReZReH(aRange, Re + h, ReH)` for
range/Doppler, and the GCOV cube, both of which are free of any terrain-slope term and so cannot
toggle on slopes. Forcing `h = 0` was tried and rejected — the true sensitivity at fixed ground
position is only 0.075 deg per 2 km of elevation, but zeroing `h` while the range still comes from
geocoding at the DEM height decouples the two and injects 0.15-0.42 deg of elevation-dependent
error instead.

The angle is computed explicitly rather than reusing the `psiE` in the pixel loop, which
`applyCorrections` leaves at 0 on the antenna-pattern path. `-linearSubPixelRTC` without `-S1Cal`
is refused, because that branch leaves `range`/`azimuth`/`h` unset.

### Refused combinations

`-min`/`-max` compare the *feathered* value, so pixels tapered at a selection boundary would win
MIN and draw lines along every boundary. `-ascending`/`-descending` and `-nearestDate` key off the
scale buffer (a passType or a date), which the fallback test would misread. All are rejected at
startup rather than silently doing something odd.

### Measured behaviour (PIG, cycle 30, 150 km box at 200 m)

| | valid px | incidence p5/50/95 | sigma0 median | fell back |
|---|---|---|---|---|
| plain | 562500 | 35.63 / 39.67 / 43.38 | -16.50 dB | - |
| `-nearRange` | 562500 | 34.36 / 36.47 / 40.64 | -15.96 dB | 0.47% |
| `-farRange` | 562500 | 44.37 / 45.64 / 46.85 | -17.09 dB | 0.30% |

Coverage identical; the gamma correction tracks the selection (gBuf median 0.95 / 1.14 / 1.56 dB).
Tiled (2x2) and untiled agree exactly, as do `-ompThreads 1` and `8`.

---

## NISAR GCOV Inputs (`-gcov`)

NISAR L2 GCOV products are already geocoded, RTC-corrected covariance terms (γ₀ for the
diagonal terms). `geomosaic -gcov file.yaml` adds them to the mosaic in a stage that runs after
the range/Doppler loop (`geoMosaic/gcovMosaic.c`, called from `makeGeoMosaic`). GCOVs can
supplement range/Doppler images or replace them entirely. Without `-gcov`, output is
byte-identical to earlier builds.

### YAML file

```yaml
polarization: HH      # HH -> HHHH, HV -> HVHV, or a full covariance term (HHHH) [HH]
frequency: A          # A or B [A]
useMask: true         # drop samples flagged 0 (invalid) or 255 (fill) in the GCOV mask [true]
glob: /Volumes/insar4/ian/Data/NISAR/Antarctica-GCOV/*.h5   # optional, may repeat; sorted
files:                # optional, "- path [weight]", weight default 1.0
  - /path/a.h5
  - /path/b.h5 0.5
```

If the same granule appears more than once with different product counters (the final `_NNN`
field, e.g. an `_001` and a reprocessed `_002`), only the highest counter is kept, and the
dropped files are logged. The date and pass direction come from the NISAR file name, so
`-date1/-date2`, `-nearestDate`, and `-ascending/-descending` behave as they do for
range/Doppler inputs.

### Reading

The HDF5 file is read directly through GDAL's HDF5 driver (no extra library). GDAL does not
attach a geotransform or CRS to GCOV subdatasets, so the grid is rebuilt from
`grids/frequencyX/xCoordinates`/`yCoordinates` (pixel centres) and the
`..._projection_epsg_code` attribute. Output grid points are converted into the GCOV CRS with
OGR (EPSG 3031/3413 → the GCOV EPSG; one transform per thread). Antarctic GCOVs are EPSG:3031
at 10 m; Greenland ones seen so far are EPSG:4326.

For each GCOV, only the window covering the output region is read. It is **block-averaged** in
linear power, `k×k` source pixels per block with `k = floor(output spacing / GCOV spacing)` per
axis, which averages rather than decimates. NaN, non-positive, and masked samples are excluded.
The averaged grid is then bilinearly sampled at each output pixel. The blocks are **aligned
to the output grid**: a block centre falls on an output pixel centre. Results therefore do not
depend on how a mosaic is tiled (tiled and untiled were identical, 0 of 10⁶ pixels different).
When the output spacing is an exact multiple of the GCOV spacing in the same CRS, each output
pixel is exactly the box average of the GCOV pixels inside it. Before this alignment, 35% of
pixels differed between tiled and untiled runs, by up to 2.5 dB. The phase anchor is the
output-lattice pixel nearest the GCOV origin (`gcovAnchor`), not the tile centre. With a tile-
centre anchor, a spacing that is not a whole number of GCOV pixels (25 m from 10 m, ratio 2.5)
still gave different phases in different tiles: 99% of pixels differed between 4×5 and 12×12
tilings, median 0.23 dB. Now 0 of 4×10⁶ pixels differ at 25 m.

**Feathering cost:** `computeScale` touches a (2·fl+1)² window for every data-edge pixel. At
25 m with a 10 km feather (fl = 400), feathering took 475 s of an 841 s tile (29 frames); with
`-fl 0` the tile took 366 s. The cost is proportional to the number of edge pixels times fl²,
and tiling does not reduce it.

**Fast feathering (2026-09-23), geomosaic only:** `computeScaleFast` (`makeGeoMosaic.c`) gets the
same weights from an exact squared Euclidean distance transform (Felzenszwalb–Huttenlocher),
over the output region plus fl. Its cost is proportional to the pixel count, and the kernel
formula is evaluated as in `fillRadialKernel`, so output is **byte-identical**. That was verified
for: a 25 m GCOV tile at fl 400 (841 s → 397 s, of which 275 s is reading); range/Doppler inputs,
uncalibrated and `-S1Cal`, at fl 20; and a mixed-weight input.

Weight quirk: `computeScale` caches its taper (`rDistSave`) built with the FIRST image's input
weight. A later image with a different weight has its interior at its own weight but is tapered
with the first image's weight. `computeScaleFast` handles the equal-weight case and falls back
to `computeScale` otherwise, so that case is identical by construction. `common/computeScale.c`,
used by every mosaic3d solver, is unchanged. Reads are made in strips of
about 1024 rows, because the layers are gzip-compressed in 512×512 chunks. Reading only `k` rows
at a time decompresses each chunk row about 512/k times, and was measured to be 14× slower
(493 s vs 34 s for one 32k×32k frame at 200 m).

### Radiometry

| Mode | Contribution to the mosaic |
|------|----------------------------|
| `-S1Cal` | σ₀ = γ₀·`rtcGammaToSigmaFactor` (averaged over the same samples) is mosaicked. The γ₀ correction `gBufTmp = 10·log10(γ₀/σ₀)` uses the same rounding and [-29.9, 35] dB clamp as the range/Doppler path, so `.sigma0` and `.gamma0` come out consistent. With `-S1Psi`, `.inc` is trilinearly interpolated from the `metadata/radarGrid/incidenceAngle` cube at the DEM height |
| no `-S1Cal` | GCOV γ₀ in linear power (σ₀ with `-calOutput sigma0`). Range/Doppler images in this mode are `power·sin³ψ`, a different quantity, so mixing the two is cosmetic only |

GCOV pixels are accumulated with the same `geoMosaicScaling` call as range/Doppler pixels,
applying weights, `-fl` feathering, min/max, orbit priority, and nearest date in the same way.
Each input counts once, so a pixel covered by one S1 image and two overlapping GCOVs is
(R + G₁ + G₂)/3 in linear power.

The `.gamma0` limitation noted elsewhere applies here too: in average mode, `gBuf` is
"last image wins" rather than averaged.

### Verification (2026-09-23, PIG, NISAR cycle 30 ascending)

- **Range/Doppler regression.** Checked against a build from git HEAD on an S1 PIG tile.
  - Uncalibrated `.unc`: byte-identical single-threaded. The multi-threaded HEAD binary itself
    varies in about 50 pixels at 1e-5 relative between runs.
  - `-S1Cal` `.sigma0/.gamma0/.inc/.gamcor`: byte-identical.
  - `-calOutput both` gives the same output as no flag.
  - `-calOutput sigma0` and `-calOutput gamma0` each write only the selected file, and it is
    identical to the default.
- **Radiometry.** For one GCOV at 200 m, against an independent Python block average of the raw
  HDF5, the median difference was -0.002 dB (σ₀) and -0.006 dB (γ₀).
- **Geolocation.** A ±1-pixel shift search is minimised at zero shift (0.075 dB RMS). The four
  one-pixel shifts give a symmetric 0.19 dB, so there is no offset.
- **Blending.** A GCOV blended with one S1 image reproduces the linear mean to within the
  0.01 dB output rounding. Pixels covered by only one of the two are reproduced exactly.
- **Threads.** `-ompThreads 1` and `-ompThreads 8` give byte-identical output.

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
- GDAL HDF5 driver + OGR/PROJ — NISAR GCOV reading and CRS transforms (`-gcov`)
- `outputGeocodedImage` / `outputGeocodedImageTiff` — binary and GeoTIFF output
- `subPixelGammaRTC` (`subpixelRTC.c`) — full sub-pixel RTC accumulation loop
- `subPixelGammaRTCLinear` (`linearSubpixelRTC.c`) — Jacobian-approximated sub-pixel RTC
- `subPixelGammaRTCJacobian` (`jacobianSubpixelRTC.c`) — $|J|$-weighted sub-pixel RTC with layover detection
- `AbAg` (`makeGeoMosaic.c`) — beta/gamma projected area computation
- GDAL — GeoTIFF/COG output
