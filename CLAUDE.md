# CLAUDE.md — mosaicSource

Supplements the root `/Users/ian/progs/GIT64/CLAUDE.md`. This directory holds the velocity-mosaic
programs and the `common/` library they all link against.

## Programs

Each program has its own subdirectory and its own `make <target>` rule in `mosaicSource/Makefile`
(run from `mosaicSource/`, e.g. `make mosaic3d ROOTDIR=... PROGDIR=... BINDIR=...`).

| Binary | Source dir | Purpose | Doc |
|---|---|---|---|
| `mosaic3d` | `Mosaic3d/` | Main velocity mosaic driver (speckle, phase, crossing-orbit, Landsat) | `Documents/mosaic3d.md` |
| `geomosaic` | `geoMosaic/` | Calibrated/uncalibrated SAR image mosaics | `Documents/geomosaic.md` |
| `siminsar` | `simInSAR/` | Simulates offset/phase products | `Documents/simInSAR.md` |
| `rparams` | `rParams/` | Baseline params for range offsets | `Documents/rparams.md` |
| `azparams` | `azParams/` | Calibration params for azimuth offsets | `Documents/azparams.md` |
| `tiepoints` | `tiePoints/` | Baseline params from unwrapped interferogram + GCPs | — |
| `lltora` | `LLtoRA/` | lat/lon/z → range/azimuth for a SAR image | `Documents/lltora.md` |
| `getlocc` | `getLocC/` | Builds a `geodat` metadata file for a multi-look product | `Documents/getLocC.md` |
| `coarsereg` | `coarseReg/` | Coarse range/azimuth alignment of two SLCs from metadata | — |
| `computebaseline` | `computeBaseline/` | Standalone baseline computation | `Documents/computeBaseline.md` |
| `offsetvrt` | `offsetVRT/` | VRT offset utility | — |
| `testgeo` | `testGeo/` | Geocoding test harness (not a production tool) | — |

`make all` builds the `TARGETS` set: `mosaic3d siminsar rparams azparams coarsereg tiepoints lltora
getlocc geomosaic` (does not include `computebaseline`, `offsetvrt`, `testgeo` — build those
individually).

## Shared code

- **`common/`** — ~45 source files used by every program above (geocoding, interpolation,
  baselines, DEM/tie-point/coordinate handling). See `Documents/common.md`. Changes here ripple
  into every binary's `make` target — rebuild affected programs, not just the one you're working on.
- **`landsatMosaic/`** — Landsat velocity integration, linked into `mosaic3d` and `geomosaic` only.

## RTC / sub-pixel radiometric terrain correction (`geoMosaic/`)

`subpixelRTC.c`, `linearSubpixelRTC.c`, `jacobianSubpixelRTC.c` implement three flavors of
sub-pixel RTC, all building cleanly. `-maskLayover` (global `int32_t maskLayover` in
`geomosaic.c`) suppresses layover pixels (output MINS1DB = -30 dB) in all three:
- `jacobianSubpixelRTC.c`: without flag, suppress when `nLayover > nValid/2`; with flag, suppress
  when `nLayover > 0` (any layover sub-pixel).
- `linearSubpixelRTC.c` / `subpixelRTC.c`: with flag, compute pixel-level Jacobian sign and
  compare against `inputImage->passType` (see root CLAUDE.md sign-convention note).

**Open investigation:** all three RTC flavors remain several dB brighter than NISAR GCOV in
steep/layover terrain. Root cause: map→SLC architecture — multiple output map pixels
independently read the same SLC bin in layover, each applying its own Ab/Ag correction, so
energy is double-counted. NISAR GCOV uses SLC→map accumulation (distributes energy once). A
two-pass multiplicity fix would be theoretically correct but ~2x cost and requires restructuring
the sub-pixel functions; not implemented — may not be worth it since layover pixels have no
recoverable signal anyway. Flat-terrain comparisons (no layover) are good; residual foreshortening
differences are likely DEM-related (datum/geoid/resampling vs Copernicus GLO-30/90), not algorithm.
Next step when revisiting: run with `-maskLayover` and re-compare against NISAR GCOV.

## ISCE/NISAR flat-earth baseline (`tiePoints/`, `common/`)

See root CLAUDE.md "ISCE/NISAR flat-earth baseline path" for the overview. Implementation:
- `tiePoints/tiepoints.c` — `-yaml` flag
- `tiePoints/computeBaseline.c` — YAML output branch (`applyFlatEarth: true` + Bn/Bp/dBn/dBp/sigma/nTiepoints/nDays)
- `common/getBaseline.c` — `.yaml` extension detection + key:value parsing
- `common/common.h` — `vhParams.applyFlatEarth` (int32_t)
- `common/computePhiFlatEarth.c` — flat-earth phase function used by `makeVhMosaic.c`
- `Mosaic3d/make3DMosaic.c` — `computePhiFlatEarthM3d`, pixel loop branches on the flag

Sample CLI: `tiepoints -bpOnly -nDays $nDays -motion -center -noDEM -dBp -yaml ../geodat30x6.in
$tiefile $phasefile baselines.orig > baseline.30x6.yaml`, then pass the original ISCE phase +
`baseline.30x6.yaml` directly to `mosaic3d` (no `changeflat` step).

## Vertical correction for submergence/emergence (`-verticalCorrection`, `Mosaic3d/`)

`mosaic3d -verticalCorrection vcFile` supplies a vertical-velocity grid (m/yr, `xyDEM` format read
by `readXYDEM`) — typically the ice-equivalent submergence/emergence rate derived from SMB, i.e.
the rate at which the ice surface moves vertically relative to horizontal flow. This vertical
motion contributes a LOS component to InSAR phase and speckle-tracked range offsets that must be
removed before solving for horizontal velocity.

- CLI flags: `-verticalCorrection vcFile` (`mosaic3d.c` ~1212) and `-verticalCorrectionSuffix
  suffix` (~1207) — the suffix check must come first since `verticalCorrection` is a substring of
  `verticalCorrectionSuffix`.
- `vcFile` → `outputImage->verticalCorrection` (`xyDEM *`, `common/geocode.h`), read once in
  `mosaic3d.c` main via `readXYDEM`.
- `common/interpVCorrect.c` — `interpVCorrect(x, y, vCorrect)` bilinearly interpolates
  `vCorrect->z` at (x,y) in km; returns `0.0` if the interpolated value is `<= MINVCORRECT` (`-100`,
  the nodata sentinel for this grid, `common/common.h`).

### Sign convention

Applied in all four velocity pixel loops with the same `X -= -dzdtSubmergence * cos(psi) * scale`
idiom (≡ `X += dzdtSubmergence * cos(psi) * scale`) used for the SHELF `tideCorrection` in the same
loops — `dzdtSubmergence` follows `tideCorrection`'s sign convention. `psi`/`aPsi`/`dPsi` is the
local incidence angle from nadir, so `cos(psi)` projects the vertical rate onto the LOS
(complementary to `sin(psi)`, used to convert LOS phase/offset to horizontal velocity).

- `Mosaic3d/make3DMosaic.c` (~308-310): `aPhase += dzdtSubmergence*cos(aPsi)*twokA*nDays/365.25`,
  and the same for `dPhase`/`dPsi`/`twokD` — applied to both ascending and descending phases
  before `computeVxy`.
- `Mosaic3d/make3DOffsets.c` (~321-323): same form applied to `aDelta`/`dDelta` (no `twok`, since
  these are range offsets in metres, not phase).
- `Mosaic3d/speckleTrackMosaic.c` (~257-258): `dr += dzdtSubmergence*cos(psi)*nDays/365.25`,
  applied to the range offset before computing `vr`.
- `Mosaic3d/makeVhMosaic.c` (~229-230): `delta += dzdtSubmergence*cos(psi)` — no `nDays/365.25`
  here because `delta` is already annualized (m/yr) via `scalePhase`.

`Mosaic3d/writeTieFile.c` (~105) uses the same grid (passed in as `vzCorrect`) additively when
writing tie points: `vz = vZimage[i][j] + interpVCorrect(x, y, vzCorrect)`.

### `-verticalCorrectionSuffix`

Companion to the per-pixel grid above, not independent of it: `rParams`/baseline fitting (and the
tiepoints behind it) is itself sensitive to the assumed submergence/emergence rate, so each
scenario needs its own baseline solution. The standard run uses the default tiepoints/baseline
fit; an alternate scenario (different `-verticalCorrection` grid) is run with tiepoints chosen for
that scenario, producing a second, suffix-named set of phase/baseline/rParams files alongside the
originals. `-verticalCorrectionSuffix` selects that alternate set so it stays paired with the
matching `-verticalCorrection vcFile`, without duplicating the whole input list.

Mechanically: appends `.<suffix>` to the **phase and baseline** filenames in
`common/getMVhInputFile.c` (`appendSuffix`), and to the **rParams** filename via
`RgOffsetsParamName` in `common/readOffsets.c` (inserted before the existing `.deltabp`/`.quad`
deltaB suffix). Propagated via `dumParams->offsets.verticalCorrectionSuffix =
outputImage->verticalCorrectionSuffix` in `Mosaic3d/setup3D.c`.

## OpenMP in `simInSAR/` (siminsar / simoffsets)

`-ompThreads N` (default 4) parallelizes the outer azimuth `iLoop` in `simInSARimage.c` with
`schedule(dynamic, 8)`. Requires per-thread `inputImageStructure` copies (`localImgs[nthreads]`,
same pattern as `make3DMosaic.c`) because `groundRangeToLLNew`/`withHeight` cache conversion state
in `cpAll`. `nIterTotal` uses `reduction(+: nIterTotal)`; progress prints are guarded to
`myThread == 0`. Remember the two-Makefile `-fopenmp` pattern (root CLAUDE.md). Not yet validated
against serial output (`-ompThreads 1` vs `N` diff) on a real scene.

## Other contents

- **`Documents/`** — per-program design notes; check here before `Documents/` in the GIT64 root.
- **`BUGFIXES_*.md`** — dated records of fixes applied from past automated code reviews
  (mosaic3d/makeVhMosaic, azParams, geoMosaic, rParams, simInSAR — all 2026-05). Historical record,
  not living docs; safe to ignore unless debugging a regression in one of those files.
- **`cloneAll`** — legacy csh script to `git clone` sibling repos (clib, cRecipes, fft,
  speckleSource, etc.) when this was a multi-repo checkout. Not needed in the GIT64 monorepo.
- **`baseline.26x16.yaml`** — sample/test config for baseline computation.
