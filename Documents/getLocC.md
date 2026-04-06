# getlocc — SAR Geodat / GeoJSON Generator

## Purpose

Converts a legacy CW-processor SAR parameter file (`.par`) into the GrIMP geodat
format (`.geodat`) and a companion GeoJSON file (`.geojson`). The geodat file provides
the image geometry (range extents, look angles, timing, state vectors, corner
coordinates) used by all downstream mosaicking and geocoding programs.

---

## Usage

```
getlocc [options] nlr nla noffset squintAngle parFile outFile passType lookDir
```

### Required Arguments (positional, last 8)

| Argument      | Description |
|---------------|-------------|
| `nlr`         | Number of range looks |
| `nla`         | Number of azimuth looks |
| `noffset`     | Integer azimuth offset of first record (lines) |
| `squintAngle` | Squint angle (degrees), or squint time (s) if `-squintTime` |
| `parFile`     | Input CW processor parameter file |
| `outFile`     | Output geodat file (`.geojson` written alongside) |
| `passType`    | `0` = descending, `1` = ascending |
| `lookDir`     | `+1.0` = right-looking, `-1.0` = left-looking |

### Options

| Option          | Description |
|-----------------|-------------|
| `-lambda <val>` | Override wavelength (m); default = ERS-1 wavelength |
| `-squintTime`   | Interpret `squintAngle` as squint time (s) rather than angle (deg) |
| `-ersFix`       | Apply 54 m near-range correction for known ERS timing error |
| `-survey`       | Divide PRF by 8 (survey-mode data) |

---

## Output Files

| File            | Contents |
|-----------------|----------|
| `outFile`       | ASCII geodat: image geometry, corner coordinates, state vectors |
| `outFile.geojson` | GeoJSON polygon feature with all parameters as properties and state vectors as array fields |

---

## Algorithm and Supporting Functions

### `correctTime`(sarD, \*squint, noffset, \*tskew, \*toffset, squintTime)

Corrects the image start time for squint and azimuth line offset. Two modes:

- **Angle mode** (`squintTime = FALSE`): converts squint angle to a time skew:

$$
t_\text{skew} = \frac{R_c \sin(\psi_\text{squint})}{\delta s \cdot f_\text{PRF}}
$$

where $R_c$ is centre range, $\delta s$ is azimuth pixel size, and $f_\text{PRF}$ is the
pulse repetition frequency.

- **Time mode** (`squintTime = TRUE`): uses the squint value directly as $t_\text{skew}$
  and back-computes the equivalent squint angle.

The line offset contribution is $t_\text{offset} = N_\text{offset} / f_\text{PRF}$.
Both corrections are added to `sarD->sec`, with carry propagated into minutes and hours.  
*Calls:* (none from common — arithmetic only)

---

### `centerLL`(sarD, sv, nla, \*lat, \*lon, deltaT)

Finds the geodetic lat/lon of the image centre by interpolating state vectors at the
centre azimuth time and iterating on the Earth radius:

1. Centre azimuth time: $t_c = t_\text{start} + \frac{(N_{az}/2)\, n_{la}}{f_\text{PRF}}$  
   (Wrapped by +86400 s if state vectors straddle midnight.)
2. Satellite ECEF position $(x_s, y_s, z_s)$ interpolated via `polintVec`.
3. Satellite radius: $R_{eH} = |{\mathbf{r}_s}|$ (corrected for ellipsoidal shape).
4. Iterates 3 times: updates $R_e = \text{earthRadius}(\varphi)$, recomputes central
   angles $\rho_\text{near}$, $\rho_\text{far}$, $\rho_\text{mid}$ from:

$$
\rho = \arccos\!\left(\frac{R^2 - R_{eH}^2 - R_e^2}{-2\, R_e\, R_{eH}}\right)
$$

   and calls `smlocateZD` to get lat/lon at the centre slant range.

*Calls:* `polintVec`, `earthRadius`, `rhoRReZReH`, `smlocateZD`

---

### `glatlon`(sarD, sv, nlr, nla, \*\*lat1, \*\*lon1, ma, mr, deltaT)

Computes geodetic lat/lon at an $ma \times mr$ grid of pixel positions across the image
(used for writing corner coordinates to the geodat). For each grid point:

1. Pixel azimuth time: $t_{ij} = t_\text{start} + \frac{a_{ij}\, n_{la}}{f_\text{PRF}} + \Delta t$  
   (Wrapped by +86400 s if state vectors straddle midnight.)
2. Slant range: $R_{ij} = R_\text{near} + r_{ij}\, \delta r_\text{SL}\, n_{lr}$
3. State vectors interpolated via `polint` (scalar Neville's algorithm, 5-point).
4. Lat/lon computed via `smlocateZD`.

*Calls:* `polint`, `smlocateZD`

---

### GeoJSON Output

The `.geojson` output encodes the image footprint as a polygon and stores all geodat
parameters as GeoJSON feature properties, including:

- Image name, date, nominal time, PRF, wavelength
- Near/centre/far range (m), Earth radii, spacecraft altitude
- Number of range and azimuth looks, SLC pixel sizes
- Corrected start time, squint offset, pass type, look direction
- All state vectors as `sv_Pos_NNN` and `sv_Vel_NNN` array fields

The raw GeoJSON output (one line from GDAL/OGR) is reformatted with indentation
by `formatGeoJSON`.

*Calls:* `readOldPar`, `correctTime`, `centerLL`, `glatlon`, `earthRadius`,
`getGeojsonDataSet`, `createFeatureDef`, `createGeometry`, `svTag`

---

## Dependencies

| Function | Source | Purpose |
|----------|--------|---------|
| `readOldPar` | `common/readOldPar.c` | Parse legacy CW `.par` file |
| `correctTime` | `getLocC/correctTime.c` | Squint and offset timing correction |
| `centerLL` | `getLocC/centerLL.c` | Image centre lat/lon |
| `glatlon` | `getLocC/glatlon.c` | Corner/grid lat/lon from state vectors |
| `polintVec` / `polint` | `common/polintVec.c` | State-vector interpolation |
| `smlocateZD` | `common/smlocateZD.c` | Slant-range → lat/lon on ellipsoid |
| `earthRadius` | `common/earthRadiusFunctions.c` | Ellipsoidal Earth radius |
| `rhoRReZReH` | `common/initRoutines.c` | Earth central angle from ranges |
| `getGeojsonDataSet` / `createFeatureDef` / `createGeometry` / `svTag` | `common/geojsonCode.c` | GeoJSON OGR output |
| GDAL / OGR | — | Vector feature I/O |
