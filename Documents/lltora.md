# lltora — Lat/Lon to SAR Range/Azimuth Converter

## Purpose

Converts geographic point locations (latitude, longitude, elevation) to SAR image
pixel coordinates (single-look range and azimuth). Useful for projecting tiepoints,
GPS positions, or DEM-derived XY grids into SAR image space for use by `rparams`,
`azparams`, or manual inspection.

---

## Usage

```
lltora geoInputFile llFile
lltora geoInputFile llFile outFile
lltora geoInputFile DEM llFile outFile
```

### Required Arguments

| Argument       | Description |
|----------------|-------------|
| `geoInputFile` | SAR image geodat parameter file |
| `llFile`       | Input point file (see format below) |
| `DEM`          | XY polar-stereographic DEM file (required for 2-column binary mode only) |
| `outFile`      | Binary output file (omit for ASCII stdout mode) |

---

## Input / Output Modes

### Mode 1 — ASCII tiepoints → stdout

```
lltora geoInputFile llFile
```

`llFile` is a standard tiepoint file with lat, lon, elevation, and optional velocities
(read by `readTiePoints`). Output is written to stdout:

```
# 5
    range  azimuth  lat  lon  elevation
    ...
&
```

Range and azimuth are in **single-look pixels** (sub-pixel float), corrected for the
multilook offset: $r_\text{SL} = r_\text{ML} \cdot n_{lr} + (n_{lr}-1)/2$.

---

### Mode 2 — Binary lat/lon/z → binary range/azimuth/z

```
lltora geoInputFile llFile outFile
```

`llFile` is a big-endian binary file:

```
[nPts nCols]  (uint32, 2 values, nCols = 3)
[lat lon z]   (float64, nPts × 3 rows)
```

Output `outFile` is big-endian binary:

```
[nPts 3]      (uint32)
[range azimuth z]  (float64, nPts × 3 rows)
```

Points that fall outside the image boundaries are written as `(-9999, -9999, z)`.

---

### Mode 3 — Binary XY grid + DEM → binary range/azimuth/z

```
lltora geoInputFile DEM llFile outFile
```

`llFile` contains only 2 columns (lat, lon) with no elevation. Heights are looked up
from the XY polar-stereographic DEM by converting lat/lon to XY via `lltoxy1` and
interpolating with `interpXYDEM`. The hemisphere, standard latitude, and grid rotation
are read from the DEM metadata.

---

## Algorithm

1. Parse the geodat file with `parseInputFile`.
2. Read tiepoints or binary lat/lon with `readTiePoints` or `readLLinput`.
3. Determine hemisphere from the first tiepoint latitude; set `HemiSphere` and `Rotation`.
4. Call `computeTiePoints`, which calls `llToImageNew` for each point:
   - Converts (lat, lon, z) to ECEF.
   - Iterates on azimuth time to satisfy the Doppler zero-Doppler condition.
   - Returns fractional multilook (range, azimuth) pixel coordinates.
5. Output single-look coordinates (multiply by number of looks, add half-look offset).

*Calls:* `parseInputFile`, `readTiePoints` / `readLLinput`, `computeTiePoints`,
`llToImageNew`, `lltoxy1`, `interpXYDEM`, `readXYDEM`

---

## Supporting Function

### `readLLinput`(fp, tiePoints, DEM)

Reads a binary lat/lon (or lat/lon/z) tiepoint file into a `tiePointsStructure`.

- Header: two `uint32` values giving `[nPts, nCols]` where `nCols` is 2 or 3.
- If `nCols == 2` (no elevation): reads the XY DEM, converts lat/lon to XY via
  `lltoxy1`, and looks up height with `interpXYDEM`. Points with negative DEM height
  are marked invalid.
- Hemisphere is inferred from the first valid latitude; `HemiSphere` and `Rotation`
  globals are set accordingly.

*Calls:* `readXYDEM`, `lltoxy1`, `interpXYDEM`

---

## Dependencies

| Function | Source | Purpose |
|----------|--------|---------|
| `parseInputFile` | `common/parseInputFile.c` | Read geodat and populate image structure |
| `readTiePoints` | `common/readTiePoints.c` | ASCII tiepoint file reader |
| `computeTiePoints` | `common/computeTiePoints.c` | Project lat/lon to SAR range/azimuth |
| `llToImageNew` | `common/llToImageNew.c` | Geodetic → SAR pixel coordinates |
| `readXYDEM` | `common/readXYDEM.c` | Load XY polar-stereographic DEM |
| `lltoxy1` | `common/lltoxy1.c` | Lat/lon → polar-stereographic XY |
| `interpXYDEM` | `common/interpXYDEM.c` | Bilinear DEM height interpolation |
