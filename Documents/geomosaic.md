# geomosaic — SAR Backscatter Geocoder and Mosaicker

## Purpose

Geocodes and mosaics one or more SAR amplitude (power) images onto a polar
stereographic output grid using a DEM for terrain correction. Supports both
uncalibrated (cosmetic) and fully calibrated outputs, including Sentinel-1 σ₀ and
radiometrically terrain-corrected (RTC) γ₀.

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

## Dependencies

- `parseInputFile` / `initllToImageNew` — geodat parsing and SAR image geometry
- `readXYDEM` — DEM I/O and projection
- `llToImageNew` — lat/lon to SAR range/azimuth projection
- `getXYHeight` — DEM elevation lookup
- `computeScale` — feathering weight computation
- `outputGeocodedImage` / `outputGeocodedImageTiff` — binary and GeoTIFF output
- GDAL — GeoTIFF/COG output
