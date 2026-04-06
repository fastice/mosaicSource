# siminsar — InSAR Phase Simulator

## Purpose

Simulates a synthetic interferogram from a DEM and optionally a surface velocity field.
For each output pixel the program iterates from ground range to find the geodetic lat/lon
on the actual DEM surface, then computes the interferometric phase from the topographic
path-length difference and (optionally) the line-of-sight surface displacement.
Also supports output of DEM heights, shelf/ice masks, and lat/lon grids.

---

## Usage

```
siminsar [options] demFile displacementFile sceneFile outputImage
```

### Required Arguments

| Argument          | Description |
|-------------------|-------------|
| `demFile`         | XY polar-stereographic DEM (binary float + `.geodat`) |
| `displacementFile`| Velocity field (`.vx`/`.vy` binary pair, or mask/VRT) |
| `sceneFile`       | SAR image geodat file defining the output image geometry |
| `outputImage`     | Output file base name |

### Options

| Option | Description |
|--------|-------------|
| `-bn <val>` | Normal baseline component $B_n$ (m) |
| `-bp <val>` | Parallel baseline component $B_p$ (m) |
| `-dBn <val>` | Linear change in $B_n$ over scene (m) |
| `-dBp <val>` | Linear change in $B_p$ over scene (m) |
| `-bnStart <val>` / `-bnEnd <val>` | $B_n$ at start / end of scene (alternative to `-dBn`) |
| `-bpStart <val>` / `-bpEnd <val>` | $B_p$ at start / end of scene (alternative to `-dBp`) |
| `-bParamsFile <file>` | Read $B_n$, $B_p$, $\delta B_n$, $\delta B_p$ from a baseline parameter file |
| `-dT <days>` | Temporal baseline for velocity displacement (default: 12 days) |
| `-flat` | Output flat-Earth-removed (flattened) phase |
| `-height` | Output DEM height values instead of phase |
| `-mask` | Output shelf/ice mask using `displacementFile` as the mask |
| `-velocity` | Use velocity field for displacement component |
| `-rPix <val>` | Override range single-look pixel size (m) |
| `-aPix <val>` | Override azimuth single-look pixel size (m) |
| `-saveLL` | Save lat/lon arrays on the geodat-defined grid (`.lat`, `.lon`) |
| `-toLL <file.dat>` | Save lat/lon on an offset-defined sub-grid |
| `-LSB` | Write output as little-endian (default: MSB big-endian) |

---

## Output Files

| File | Contents |
|------|----------|
| `outputImage` | Simulated phase (float32, radians) or height (float32, m) or mask (byte) |
| `outputImage.vrt` | VRT descriptor for the output image |
| `outputImage.simdat` | ASCII header: image size, baseline start/end, DEM and displacement filenames |
| `outputImage.lat` / `.lon` | Lat/lon grids (float64) if `-saveLL` or `-toLL` |
| `outputImage.ll.vrt` | VRT for lat/lon pair |

---

## Algorithm

### `simInSARimage`(scene, dem, xyVel)

Main simulation loop. Iterates over azimuth lines then range pixels.

*Calls:* `initllToImageNew`, `getReH`, `thetaRReZReH`, `groundRangeToLLNew`,
`earthRadius`, `rhoRReZReH`, `getXYHeight`, `slantRange`, `psiRReZReH`,
`lltoxy1`, `interpXYVel`, `computeXYangle`, `computeHeading`, `getShelfMask`

#### 1. Ground-Range Iteration

For each output pixel, the program iterates from a starting ground range $\rho$
on a spherical Earth to find the lat/lon on the actual DEM surface:

1. `groundRangeToLLNew` converts spherical ground range to lat/lon and returns slant range $R_0$.
2. The WGS84 Earth radius at that lat/lon is computed: $R_{eWGS} = \text{earthRadius}(\varphi)$.
3. The Earth central angle on the ellipsoidal earth: $\rho_{WGS} = \text{rhoRReZReH}(R_0, R_{eWGS}, R_{eH})$.
4. DEM height $h_{WGS}$ is looked up at (lat, lon).
5. Height-corrected slant range:

$$
R = \sqrt{R_{eH}^2 + (R_{eWGS} + h_{WGS})^2 - 2\,R_{eH}(R_{eWGS} + h_{WGS})\cos\rho_{WGS}}
$$

6. Ground range is incremented by $\frac{1}{2}(R_\text{target} - R)$ until $R \geq R_\text{target} - 0.1$ m.

#### 2. Precise Lat/Lon at DEM Height

After convergence, `withHeight` is called to refine lat/lon at the actual DEM surface
elevation using `smlocateZD` with the target height passed directly.  
*Calls:* `polintVec`, `smlocateZD`

#### 3. Topographic Phase

The look angle at DEM height $h_{sp}$ (spherical):

$$
\theta = \text{thetaRReZReH}(R, R_e + h_{sp}, R_{eH}), \qquad \theta_D = \theta - \theta_c
$$

Topographic path-length difference (same formula as baseline phase in InSAR):

$$
\Delta = \sqrt{R^2 - 2R(B_n\sin\theta_D + B_p\cos\theta_D) + B_n^2 + B_p^2} - R
$$

If `-flat` is set, the flat-Earth contribution is subtracted:

$$
\Delta \mathrel{-}= -B_n\sin\theta_D - B_p\cos\theta_D + \frac{B_n^2 + B_p^2}{2R}
$$

Phase: $\phi_\text{topo} = \Delta \cdot \frac{4\pi}{\lambda}$

#### 4. Velocity Displacement (optional, `-velocity`)

If a velocity field is provided, the range-direction surface displacement is computed:

1. Convert (lat, lon) to XY and interpolate $(v_x, v_y)$ with `interpXYVel`.
2. Compute the polar-stereographic x-axis angle $\psi_{xy}$ and satellite heading $H$.
3. Project velocity into the ground-range direction:

$$
v_r = v_x\cos(H - \psi_{xy}) - v_y\sin(H - \psi_{xy})
$$

4. DEM slope in the range direction $\partial z/\partial r$.
5. Line-of-sight displacement:

$$
\delta r = v_r\sin\psi - v_r\,\frac{\partial z}{\partial r}\cos\psi
$$

where $\psi$ is the local incidence angle. Phase contribution:

$$
\phi_\text{vel} = \delta r \cdot \frac{\Delta t}{365.25} \cdot \frac{4\pi}{\lambda}
$$

Total output phase: $\phi = \phi_\text{topo} + \phi_\text{vel}$

#### 5. Baseline Along-Track Variation

The baseline varies linearly with azimuth line index:

$$
B_n(i) = B_{n,\text{start}} + i \cdot \frac{B_{n,\text{end}} - B_{n,\text{start}}}{N_{az} - 1}
$$

(and equivalently for $B_p$). This allows simulation of realistic along-track baseline
drift from orbit geometry.

---

## Supporting Functions

### `parseSceneFile`(sceneFile, scene) → void
Parses the SAR geodat file (`sceneFile`) into `scene->I` via `parseInputFile`.
Sets the baseline step sizes from start/end values and image size.
Allocates output arrays (`scene->image`, `scene->latImage`, `scene->lonImage`).
If `-toLL` mode, reads the offset grid geometry from a `.dat` or `.vrt` parameter file
(via `parseOffsetParamFile`) and sizes the output to match.  
*Calls:* `parseInputFile`, GDAL (`GDALOpen`, `readDataSetMetaData`)

---

### `interpXYVel`(x, y, xyVel, \*vx, \*vy) → void
Bilinearly interpolates $v_x$ and $v_y$ from an `xyVEL` velocity map at XY coordinates.
Returns `-LARGEINT` components if outside the valid grid.  
*Calls:* `bilinearInterp`

---

### `outputSimulatedImage`(scene, outputFile, demFile, displacementFile) → void
Writes all output files:
- Calls `outputLL` if lat/lon output is requested (`-saveLL` or `-toLL`), writing
  `.lat` / `.lon` binary files and a `.ll.vrt` descriptor.
- Calls `outputSimImage` to write the phase / height / mask image and its `.vrt`.
- Writes a `.simdat` ASCII header with image size and baseline information.

*Calls:* `writeSingleVRT`, `appendSuffix`, `fwriteOptionalBS`

---

### `parseBnBpParamsFile`(bpParamsFile, \*bn, \*bp, \*dbn, \*dbp) → int
Reads $B_n$, $B_p$, $\delta B_n$, $\delta B_p$ from a 5-column baseline parameter file
(same format as `computeBaseline` output: `Bn Bp dBn x dbp`). Returns 0 on success,
-1 on error.

---

## Dependencies

| Function | Source | Purpose |
|----------|--------|---------|
| `parseInputFile` | `common/parseInputFile.c` | Read SAR geodat into scene image structure |
| `initllToImageNew` | `common/llToImageNew.c` | Initialise geocoding geometry |
| `groundRangeToLLNew` | `common/groundRangeToLLNew.c` | Ground range → lat/lon iteration |
| `getXYHeight` | `common/getXYHeight.c` | DEM height at lat/lon |
| `thetaRReZReH` / `psiRReZReH` | `common/initRoutines.c` | Look / incidence angles |
| `rhoRReZReH` / `slantRange` | `common/initRoutines.c` | Earth central angle, slant range |
| `earthRadius` | `common/earthRadiusFunctions.c` | WGS84 Earth radius |
| `getReH` | `common/initRoutines.c` | Azimuth-varying satellite radius |
| `smlocateZD` | `common/smlocateZD.c` | Satellite position → surface lat/lon |
| `polintVec` | `common/polintVec.c` | State-vector interpolation |
| `lltoxy1` | `common/lltoxy1.c` | Lat/lon → polar-stereographic XY |
| `interpXYDEM` | `common/interpXYDEM.c` | DEM interpolation at XY |
| `computeHeading` | `common/computeHeading.c` | Satellite heading at pixel |
| `computeXYangle` | `common/computeXYangle.c` | Map x-axis angle at pixel |
| `getShelfMask` | `common/getShelfMask.c` | Shelf/grounded mask lookup |
| `writeSingleVRT` | `gdalIO/gdalIO/gdalIO.c` | Write VRT descriptor for output |
| `readXYDEM` / `readXYVel` | `common/readXYDEM.c` | Read DEM and velocity maps |
