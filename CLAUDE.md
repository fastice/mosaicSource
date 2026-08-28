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

## Ionosphere correction on the phase path (`tiepoints -ionosphere`, `Mosaic3d/`)

The phase-side counterpart of `rparams`' ION_AUTO (see `rParams/rparams.c:294-350`), added for
Sentinel-1: ISCE's `topsApp` split-spectrum step produces an ionospheric phase screen
(`merged/topophase.ion`), which `Sentinel1Phase/setupisceuw.py` already exports into the GrIMP
product directory as `ionosphere.tif`. Nothing consumed it before. NISAR is unaffected — without
`-ionosphere` the behaviour is exactly as it was.

**`tiepoints` flags** (names deliberately match `rparams`; defined in `tiePoints/tiePoints.h`):
`-ionosphere <file>` (absent = old behaviour), `-noIonosphere` (ION_NONE), `-forceIonosphere`
(ION_FORCE), `-ionSigmaMargin <frac>` (default `IONSIGMAMARGIN` = 0.05). Default with a file
given is ION_AUTO. Note `readArgs`'s `argc` guard had to go 20 → 28, and the flag loop uses
`strstr`, so the four new branches are ordered longest-first.

**Sign: SUBTRACT.** The file holds the ionospheric phase itself, in the same sign convention as
the phase image it accompanies, so `phaseCorrected = phase - ion`. This is ISCE's own convention —
`runMergeBursts.py` applies its screen as `a*exp(-1.0*J*b)` with `b = topophase.ion`. It is the
**opposite** of the range-offset ionosphere correction (`rParams/getROffsets.c`,
`Mosaic3d/make3DOffsets.c`), which is a pre-negated *correction* that gets added — see root
CLAUDE.md. Both new apply sites carry a comment saying so.

**Python-side gotcha (not fixable in C):** `setupisceuw.py` writes the GrIMP phase as
`mySign * unw` where `mySign = sign(date2 - date1)`, but writes `ionosphere.tif` with **no**
`mySign`. For pairs with `mySign == -1` the exported screen is sign-flipped relative to its own
phase. The C code assumes the two share a convention and cannot detect a violation.

**Flow:**
- `tiePoints/getPhases.c` — the image read is factored into `readRadarImage()` (GDAL for
  `.vrt`/`.tif`, `freadBS` otherwise, hard error on a dimension mismatch); new `getIonosphere()`
  samples the screen with the *same* `interpolatePhase(r[i], a[i], ...)` bilinear call into
  `tiePoints->ionPhase` (`common/common.h`, alongside `phaseSquint`).
- Sampling once up front is equivalent to sampling after the corrections, because every downstream
  correction (`addBaselineCorrections.c:261`, `addMotionCorrections.c:106`) is a purely additive
  perturbation of `phase[]` independent of its value.
- `tiePoints/computeBaseline.c` — `computeBaseline()` runs `fitBaseline()` a second time on
  `phase - ionPhase` and keeps the winner in `noSquintFit`, so everything downstream is unchanged.
  Rule: ION_FORCE, or the uncorrected fit found no solution → use ion; otherwise
  `sigmaIon < sigmaNoIon * (1 - margin)`. `-debug` writes only the winner's residuals.
- Output — YAML gains four flat top-level keys before `noSquint:` (`ionosphereCorrectionFile`,
  which is `nil` when unused; `usingIon`; and, when both fits ran, `sigmaWithIonCorrection` /
  `sigmaWithoutIonCorrection`), emitted on the `sigma: -1` sentinel path too. Legacy text output
  gains `;* ionosphereCorrectionFile <path>` after the `&`, only when the correction won.
- `common/getBaseline.c` — `setIonospherePhaseFile()` parses both formats into
  `vhParams.ionospherePhaseFile`, resolving a relative name against the **baseline file's own
  directory** (absolute used as-is) and hard-erroring if it does not exist.
- `Mosaic3d/setup3D.c` copies it onto `inputImageStructure.ionospherePhaseFile`
  (`common/geocode.h`, next to `image`), clears it on the `sigma<0` → `nophase` downgrade, and
  `allocateOffsetBuffers()`/`setupADImageBuffers()` allocate the `AIonBuffer`/`DIonBuffer` pools
  (`common/buffers.c`) **only** when some baseline actually named a correction — so runs without
  one have exactly their old memory footprint. `ionospherePhase == NULL` is what tells the pixel
  loops there is no correction.
- `common/readOffsets.c:getIonospherePhaseImage()` reads it, honoring the same
  `yMin/nAzimuthLooks .. yMax/nAzimuthLooks` window as `getMosaicInputImage()`. It opens with
  `GDALOpen` **directly** rather than via `checkForVrt()`, so a plain `.tif` works with no `.vrt`
  wrapper — `checkForVrt()` would have silently fallen through to the raw big-endian path. A NULL
  `buf` means "keep the pool `setupADImageBuffers` assigned", which is what `makeVhMosaic.c` needs
  since it never calls `setBuffer()`.
- Applied via `interpIonPhaseImage()` (`common/interpPhaseImage.c`) in the **radians** domain,
  right after the `phase - phiZ` topography removal and before the velocity scaling:
  `make3DMosaic.c` (both `aPhase` and `dPhase`) and `makeVhMosaic.c`. No mosaic3d CLI flag —
  it unconditionally honours what the baseline file recorded, exactly as the offsets side does.
  `make3DOffsets.c`/`speckleTrackMosaic.c` are untouched.

**Verified** (2026-08-20, real NISAR frame `shadow/track-1/4252_0000`): a no-flag rerun differs
from the committed `baseline.26x16.yaml` only by the two new lines, with numeric drift inside
`tiepoints`' own run-to-run spread (two consecutive runs of the *same* binary differ by ~8e-5 in
sigma — the nondeterminism documented for `rparams`/`azparams` applies here too). A zero screen
gives bit-identical sigmas and correctly loses under AUTO / wins under `-forceIonosphere`. A
synthetic azimuth-varying screen added to the phase is recovered exactly: sigma 7.036 (uncorrected)
→ 3.271 and `Bp` 0.02968 vs the uncorrupted frame's 3.275 / 0.02959; `-ionSigmaMargin 0.99`
correctly rejects it. `mosaic3d` round-trips both baseline formats, resolves the relative path,
allocates the ion pools only when needed, and shifts velocities by ~0.94 m/yr per radian of screen.

**Not yet exercised:** the `make3DMosaic.c` (crossing asc/desc pair) apply site — the verification
above only drove `makeVhMosaic.c`, since no crossing phase pair with a correction was available.
The code is symmetric with the `makeVhMosaic.c` path that was verified.

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
that scenario, producing a second, suffix-named set of baseline/rParams files alongside the
originals. The phase file itself is unaffected by the vertical-correction scenario and is shared
unchanged between them. `-verticalCorrectionSuffix` selects the alternate baseline/rParams set so
it stays paired with the matching `-verticalCorrection vcFile`, without duplicating the whole
input list.

Mechanically: appends `.<suffix>` to the **baseline** filename only (not phase) in
`common/getMVhInputFile.c` (`appendBaselineSuffix`), and to the **rParams** filename via
`RgOffsetsParamName` in `common/readOffsets.c` (inserted before the existing `.deltabp`/`.quad`
deltaB suffix). Propagated via `dumParams->offsets.verticalCorrectionSuffix =
outputImage->verticalCorrectionSuffix` in `Mosaic3d/setup3D.c`.

## Squint (residual Doppler) sensitivity — analyzed and implemented, off by default

Real NISAR acquisitions carry a small residual squint (~1.5°–1.7°, measured directly — see
below) that `mosaic3d`'s geometry code does not model by default: `common/llToImageNew.c`/
`computeHeading.c` assume exactly-zero-Doppler (broadside) LOS. **Conclusion of the original
analysis: no fix is warranted for the current NISAR processing chain** — the correction below
exists and is verified, but ships **off by default** (`-useSquint`), since real NISAR data doesn't
need it today. Full derivation, both error-propagation figures, and the reasoning live in the
`nisarErrors` Python package (`~/PycharmProjects/packages/nisarErrors/Documents/plotSquintError.md`)
— this section is a pointer + the minimum needed to reconstruct or extend that analysis without
re-deriving it.

**What was measured:** squint = (heading of `losUnitVector` − heading of
`alongTrackUnitVector`, wrapped) − 90°, read directly from a real RUNW product's
`science/LSAR/RUNW/metadata/geolocationGrid` HDF5 cube (`losUnitVectorX/Y`,
`alongTrackUnitVectorX/Y` — NISAR ATBD, JPL D-95677 Rev A, §3.4–3.8). Not computable from orbit
state vectors alone — NISAR ATBD §3.8 splits Doppler centroid into a geometric component (needs
spacecraft *attitude*, not just orbit) and a measured component (empirically calibrated from raw
data) — so reproducing it in C would require either reading this HDF5 group from Python and
passing a derived value through, or accepting it as a documented, bounded limitation (the
recommendation here).

**How it propagates (the math, if extending the analysis):** `mosaic3d` builds its solving
matrix from assumed zero-squint headings, $\mathbf{A}_0=\text{computeA}(H_A(0),H_D(0))$
(`common/initRoutines.c` `computeA`). The true matrix is
$\mathbf{A}_{\text{true}}=\text{computeA}(H_A(0)+\text{squint}_A,\,H_D(0)+\text{squint}_D)$, and
mosaic3d effectively computes $\vec v_{\text{computed}}=\mathbf{A}_0\mathbf{A}_{\text{true}}^{-1}\vec v_{\text{true}}=\mathbf{M}\vec v_{\text{true}}$.
At the real measured squint, $\mathbf{M}$ is dominated by its antisymmetric part (close to a
pure small-angle rotation) — net effect is a near-constant ~1.5°–1.7° direction bias in computed
flow direction (regardless of true flow orientation) plus a sub-0.3% speed error.

**Why offsets (range/azimuth, `make3DOffsets.c`) don't need a fix:** the zero-Doppler condition
$\vec V_{sat}(\eta_{0,T})\cdot(\vec T-\vec R_{sat}(\eta_{0,T}))=0$ forces the true LOS to be
*exactly* perpendicular to the true velocity at the assigned time, by construction — so
range-sensitivity ⊥ azimuth-sensitivity regardless of squint. The only residual risk is whether
`computeHeading.c`'s own state-vector-derived heading matches the true acquisition geometry
closely enough — bounded above by the same analysis (sub-percent).

**Why phase (`computePhiZ.c`/`computePhiFlatEarth.c`) has a separate exposure, but it's also
negligible for NISAR:** these convert the inter-orbit baseline ($B_n$,$B_p$) to a range
correction via `thetaRReZReH(Range, Re, ReH)` (`common/initRoutines.c:918`) — a pure spherical
law-of-cosines on three scalar distances, with no heading/along-track term at all (a hard-coded
zero-squint assumption, independent of `computeHeading`). Two findings limit this in practice:
(1) `common/svBase.c`'s `svBaseTCN()`/`svBnBp()` (lines 271–385) already iteratively solve for
the second image's matching time specifically to drive the baseline's along-track component to
~0 (`dt = dot(b,T)/C1`, converged to `1e-11`) before rotating only the cross-track/normal
components by `thetaRReZReH`'s angle — structurally guards against the missing along-track term
mattering much; (2) for NISAR's flat-earth path (`applyFlatEarth=1`, see "ISCE/NISAR flat-earth
baseline" above), $B_n$/$B_p$ are already just the small residual orbit-error ramp left after the
official processor removed the bulk geometric phase — not the full physical baseline — so any
squint-induced projection error acts on an already-small input. **This argument does not extend
to legacy/non-NISAR full-physical-baseline ISCE data through `computePhiZ.c` without
`applyFlatEarth`** — flagged as a known, unquantified gap, not a current concern.

**The exact fix, if one is ever needed (not a post-hoc approximation):** `computeA(α,β)` is the
matrix inverse of $N(\alpha,\beta)=\begin{pmatrix}\cos\beta&\sin\beta\\\cos(\alpha+\beta)&\sin(\alpha+\beta)\end{pmatrix}$,
whose rows are the two images' ground-projected look directions; squint rotates each image's own
row by that image's own squint angle. A tempting shortcut — rotate the *output* $(v_x,v_y)$ by
$-(\text{squint}_A+\text{squint}_D)/2$ after `computeVxy` — works well only when
$\text{squint}_A\approx\text{squint}_D$ (true for the measured Greenland values, ~14× error
reduction) and degrades badly once they diverge (a hypothetical 3°/−1° case barely improves). The
**correct, exact fix instead applies each measured squint to its own image's heading before**
$\alpha,\beta$ **are computed**: $H_A\to H_A+\text{squint}_A$, $H_D\to H_D+\text{squint}_D$, then
proceed through `computeHeading.c`'s existing $\alpha=H_A-H_D$, $\beta=\phi-H_A$ unchanged. This
gives $\mathbf{A}_{\text{true}}$ directly with zero residual (to the precision of the measured
squint), with no degradation as $\text{squint}_A,\text{squint}_D$ diverge — see
`nisarErrors/Documents/plotSquintError.md` §"Exact structure of the error..." for the proof.

**Design (implemented — extraction, merge, and C-code consumption all done):**

*Functional form — a 6-parameter polynomial in range $r$ and azimuth $a$, not a single scalar or
a 1D fit in either variable alone:*

$$
\text{squint}(r,a) = c_0 + c_1 r + c_2 a + c_3 a^2 + c_4\,ra + c_5 r^2
$$

Range and azimuth need different treatment because they differ in kind, not just degree.
**Range** is bounded by one swath width (only ~2× wider for 40 MHz vs. narrower bandwidths) and
already confirmed close to linear within a frame (linear-fit residual std ~0.0085°, measured
directly from real RUNW `geolocationGrid` data) — $c_1 r$ (+ $c_5 r^2$ as cheap insurance for
wider swaths than tested) captures this, and it's the *dominant* within-frame effect (std ~0.086°
vs. ~0.009° for azimuth in the one frame measured so far). **Azimuth** varies little within a
single sub-frame, but a virtual frame can span 20+ merged sub-frames — enough along-track
distance that slow Doppler/heading drift, or sharper behavior approaching a track's high-latitude
turning point, could plausibly introduce real curvature that Greenland data structurally cannot
rule out (these tracks never approach that regime). $c_2 a + c_3 a^2$ covers linear drift plus
that possible curvature; $c_4\,ra$ allows the range-slope itself to drift along the track. If
$c_3,c_4,c_5$ fit out near zero on real data, that's confirmation the simpler form would have
sufficed — no harm in carrying them.

*Pipeline path (all three stages done):*
1. **`RUNWtoGrimp`/`nisarhdf`** (`~/PycharmProjects/packages/nisargrimpworkflow`/`nisarhdf` — see
   their CLAUDE.md files): `nisarBaseHDF.getSquintAnglePolynomial()` samples squint from the
   RUNW's `geolocationGrid` cube via `losUnitVectorCube()`/`alongTrackUnitVectorCube()` at a grid
   of $(r,a)$ points spanning each sub-frame, and fits the local 6 coefficients. Reference image
   only (not secondary — it shares the reference's HDF5/cube, which can't describe the
   secondary's own acquisition geometry).
2. **`SetupNISAR.mergeSquintAnglePolynomial()`**: since the 2D fit can't be merged by averaging
   per-sub-frame coefficients (each is only valid over its own local azimuth window), it resamples
   every sub-frame's own fitted surface onto the merged range span × that sub-frame's own azimuth
   span, pools the samples, and refits the single 6-parameter polynomial over the combined domain
   — the same "redistribute onto one common reference" idea `mergeStateVectors()` already uses.
3. **C-code consumption** (`mosaic3d`, **phase path only — `make3DMosaic.c`**, not
   `make3DOffsets.c`, since offsets are self-consistent regardless of squint by construction — see
   above): `common/parseInputFile.c:parseGeojson()` reads three new flat geodat fields
   (`squintCoefficients`, `squintRefRange`, `squintRefAzimuthTime` — flat because GDAL's OGR
   GeoJSON driver can't read the nested `squintAnglePolynomial` dict the Python side also keeps)
   into four new `inputImageStructure` fields (`common/geocode.h`). `computeA()`
   (`common/initRoutines.c`) takes a new `applySquint` parameter; when true, it calls a new static
   `evaluateSquint()` helper (right next to `computeA`) for each image immediately after the
   existing `computeHeading()` calls, **before** $\alpha,\beta$ are computed — exactly the
   insertion point derived above, confirmed correct by `computeHeading.c`'s own header comment
   ("the angle of the cross-track... direction from north"), i.e. exactly the idealized
   broadside-LOS heading that squint is defined as a correction to. `make3DMosaic.c`'s call site
   passes the new global `useSquint` flag (`-useSquint`, default off, mirrors `-maskLayover`'s
   pattern); `make3DOffsets.c`'s call site passes `FALSE` explicitly, with a comment — never
   applies the correction, flag or no flag, not an oversight.

**Verified** (2026-06-28/29): a standalone numeric test (real North Greenland ascending/descending
geodats, squint_A≈1.658°, squint_D≈1.496°) found the rotation implied purely by the `A` matrices
(`M = A_off @ inv(A_on)`, decoupled from any specific test vector) is +1.587°, matching the
predicted $(\text{squint}_A+\text{squint}_D)/2≈1.577°$ to within 0.6%, with `M` confirmed dominated
by its antisymmetric/rotation part (det≈0.9994) exactly as predicted above. A real end-to-end
`mosaic3d` run on an actual overlapping crossing pair (track-58 ascending × track-97 descending)
gave a consistent-sign rotation (median -1.17° over all crossing-affected pixels, tightening to
median -1.26°, std 0.33° when restricted to faster/less-noise-dominated pixels) and a speed
difference consistent with the sub-0.3% prediction (median +0.03% for the fast-pixel subset).
**`tiepoints`'s two separate squint exposures — one analyzed and closed, one fixed.**
`tiepoints`/`computeBaseline.c`'s baseline-decomposition (the `Bn`/`Bp` fit via `thetaRReZReH`)
was analyzed and found to need **no fix**: `common/svBase.c`'s `svBaseTCN()` already drives the
baseline's along-track TCN component ($B_T$) to ≈0 via an iterative solve that only depends on
the inter-satellite position difference and satellite 1's own velocity direction — never the
antenna pointing direction — so $B_T\approx0$ regardless of squint. Decomposing the true
(squint-tilted) LOS in TCN as tilted by angle $\delta\approx$squint away from the broadside C-N
plane toward T, $\hat L_{\text{true}}=(\sin\delta,\cos\delta\sin\theta,\cos\delta\cos\theta)$, the
baseline-LOS dot product becomes $\vec B\cdot\hat L_{\text{true}} = B_T\sin\delta + B_p\cos\delta$
— the term squint would leak into ($B_T\sin\delta$) is already suppressed, leaving only a
**second-order** ($\propto\delta^2/2$) residual on $B_p$, unlike the first-order heading effect
that mattered for `mosaic3d`. Checked against a real flat-earth baseline (`track-58/3271_0000`,
$B_p$=3.7cm): residual range error ≈15 microns, phase error ≈0.0008 rad — about 1/8000th of a
fringe. This sharpens the legacy non-flat-earth gap mentioned below: the same $B_p\delta^2/2$
formula gives ~2 rad for a 100m baseline, real but out of scope for the current NISAR/flat-earth
project.

`tiePoints/addMotionCorrections.c` (enabled by `-motion`, used by every real tie-point run) is a
**separate, first-order exposure that does need a fix, and got one**: it projects each tie point's
*known* ice-velocity vector onto the LOS via `hAngle = computeHeading(...)` — the same idealized,
broadside heading `mosaic3d` had to correct — to build `vyra`, then converts that directly to a
phase correction subtracted from the tie point's measured phase before the baseline fit ever runs.
Fixed the same way as `mosaic3d`: `hAngle += evaluateSquint(&inputImage, range1, azimuth1) *
DTOR` right after `computeHeading()` returns (`evaluateSquint()` was made non-`static` in
`initRoutines.c` and prototyped in `common.h` so both binaries can call it). **Verified** on a real
geodat (track-58, squint≈1.658°) with a synthetic 224 m/yr tie-point velocity: the exact derivative
prediction `d(vyra)/d(rotAngle) * squint_rad` = -3.349 vs. measured `vyra` change = -3.292 (1.7%
agreement, correct sign).

**Dual-solution design (2026-06-30): `tiepoints` computes both solutions; `mosaic3d -useSquint`
picks one.** Originally `tiepoints`' own `-useSquint` flag gated which single corrected phase was
computed, separately from `mosaic3d`'s own `-useSquint`. This was redesigned so the two flags
aren't independent: `tiepoints` (`tiepoints.c` main, around the motion-correction block) now
always runs `addMotionCorrections()` with `applySquint=FALSE` into `tiePoints.phase` (unsquinted,
default), and additionally — whenever `inputImage.hasSquintPolynomial` is true and `-motion`/`-vr`
is given — runs it a second time with `applySquint=TRUE` into a new parallel array
`tiePoints.phaseSquint`, setting `tiePoints.hasSquintSolution=TRUE`
(`addMotionCorrections()`'s signature changed to take an explicit `applySquint` flag and write
into a caller-supplied `phaseOut` array rather than mutating `tiePoints->phase` in place, so it
can run twice without interference; `tiePoints`/`common.h`). `tiepoints`' own `-useSquint` CLI
flag is now a deprecated no-op (prints a note) — both solutions are always computed when squint
data is available, so old callers that still pass it keep working unmodified.

`tiePoints/computeBaseline.c`'s `-yaml` output (`computeBaseline()`) fits both phase arrays
(`fitBaseline()`, factored out of the old single-solution function) and writes them as two
labeled, nested blocks under top-level `hasSquintSolution: true/false`:
```yaml
applyFlatEarth: true
hasSquintSolution: true
nDays: 12.000000
noSquint:
  nTiepoints: 412
  sigma: 0.083000
  Bn: ...
  ...
  C:
    - [...]
squint:
  nTiepoints: 412
  sigma: 0.081000
  Bn: ...
  ...
```
The sigma<0 "no solution" sentinel (see root CLAUDE.md) is emitted per-block the same way. Legacy
(non-`-yaml`) text output is unchanged — unsquinted solution only, since that format and its
consumers (`rparams`, `computebaseline` standalone) have no notion of squint.

`common/getBaseline.c` (`getBaseline()`, now takes a 4th `useSquint` argument) parses this by
stripping leading whitespace before its existing key-matching `sscanf` calls — so the indented
child keys under `noSquint:`/`squint:` parse with identical code to old flat (un-nested,
single-solution) files, which lack section headers entirely and are routed into the `noSquint`
block by default for backward compatibility. It selects the `squint` block when `useSquint` is
true and `hasSquintSolution` was true in the file, otherwise `noSquint`; requesting `-useSquint`
against a baseline file with no squint block is an error (`error(...)`), not a silent fallback.
The sole call site is `Mosaic3d/setup3D.c` (`extern int32_t useSquint;`, mosaic3d's existing
global, same one that gates the phase-path heading correction in `make3DMosaic.c`) — so
`mosaic3d -useSquint` now does two things from one flag: selects the squint-corrected baseline
solution from the YAML *and* applies the squint heading correction in `computeA()`.

**Plumbed end-to-end (Python side still on the old single-flag design — not yet updated for the
dual-solution `tiepoints` behavior above; deferred, see below)**: `project.yaml`'s
`applySquintCorrection` (default `false`) → `setupNISARTracks.py --useSquint`/`--noUseSquint`
(CLI overrides the project default when given, resolved once at this top level) →
`refreshties.py --useSquint` → `makeframetie.py --useSquint` → `tie_script`/`tieScript.py
--useSquint`, which appends `-useSquint` to `info['extra_flags']` (reaches `mosaic3d`) and to
`info['default_p_flags']`/`info['default_yaml_p_flags']` (reaches the `tiepoints -motion` calls)
— three places, since `tieScript.py` is the only point in this chain invoking both binaries.
Deliberately **not** threaded into `maketies -run` (baseline planning, a separate concern from
where the binaries actually run). See `nisargrimpworkflow`/`mosaicworkflow`/`insarworkflow`
CLAUDE.md files for the package-level pointers. **Not yet updated**: since `tiepoints -useSquint`
is now a no-op, these Python callers can stop passing it to `tiepoints` (only `mosaic3d` needs it
now) — left as-is for this pass since it still works unmodified (no-op flag), to be cleaned up in
a follow-up.

**Open gaps, not yet resolved:** whether squint varies non-monotonically over a much wider
heading/latitude range than the three Greenland crossing pairs analyzed span (e.g. a long
Antarctic track); the legacy non-flat-earth `computePhiZ.c` exposure, now quantified above as
potentially ~2 rad for a 100m baseline but still unaddressed since out of scope for NISAR.

A related, structurally-identical-geometry question — how common-mode vertical-motion errors and
independent range/phase noise propagate through `computeA`/`computeVxy` — is documented in
`Documents/mosaic3d.md` §"Error Analysis: Crossing-Geometry Sensitivity" and visualized in
`nisarErrors`' `plotVerticalSensitivity.py`/`Documents/plotVerticalSensitivity.md`.

## OpenMP in `simInSAR/` (siminsar / simoffsets)

`-ompThreads N` (default 4) parallelizes the outer azimuth `iLoop` in `simInSARimage.c` with
`schedule(dynamic, 8)`. Requires per-thread `inputImageStructure` copies (`localImgs[nthreads]`,
same pattern as `make3DMosaic.c`) because `groundRangeToLLNew`/`withHeight` cache conversion state
in `cpAll`. `nIterTotal` uses `reduction(+: nIterTotal)`; progress prints are guarded to
`myThread == 0`. Remember the two-Makefile `-fopenmp` pattern (root CLAUDE.md). Not yet validated
against serial output (`-ompThreads 1` vs `N` diff) on a real scene.

## Debug residual GeoPackage output (`-debug`/`-outputFile`, tiepoints/rparams/azparams)

`tiepoints`, `rparams`, and `azparams` each fit a least-squares baseline/calibration
solution against SAR tie points but historically only reported an aggregate sigma. A
`-debug` flag on all three now writes every tie point used in the fit, plus its
per-point residual, to a GeoPackage — viewable directly in QGIS/geopandas with both
geographic (lat/lon, polar-stereo x/y) and image (range/azimuth) coordinates, so
points can be overlaid on the original SAR frame.

- **Shared writer**: `common/writeTieResidualsGpkg.{c,h}` (new files) — plain OGR/GDAL
  vector writer modeled on the existing `common/geojsonCode.c` pattern, targeting the
  `GPKG` driver instead of `GeoJSONSeq`. `libgdal` (already linked into every binary via
  `$(GDAL)`) provides both GPKG read/write and the SRS machinery; no new dependency.
  Layer schema: `id`, `lat`, `lon`, `x_km`, `y_km`, `range`, `azimuth`, `z`, `weight`,
  plus a program-specific residual field (`phase_residual_rad` / `range_residual_m` /
  `azimuth_residual_m`). CRS comes from `getEPSGFromProjectionParams()`
  (`gdalIO/gdalIO/tiffWriteCode.c`) — pass `tiePoints->stdLat`, **not** the global
  `SLat`, since none of these three programs ever populate `SLat` (it stays at its
  `-91.0` sentinel; only `tiePoints->stdLat`, set via
  `setTiePointsMapProjectionForHemisphere()`, is real).
- **Per-point residuals**: none of the three programs previously retained them —
  `svdfit()` (`cRecipes/svdfit.c`) computes them internally just to accumulate `chisq`
  and discards them. Each `compute{Baseline,RParams,AzParams}()` now recomputes
  residuals after the winning `svdfit()` call by re-invoking the same coefficient
  function and dotting with the fitted parameter vector — cheap (O(npts·nParams)), no
  `svdfit()` changes needed. Requires a small `origIndex[]` array mapping the 1-based
  SVD fit index back to the tie point's original array index; **stored 0-based**
  (`origIndex[i1 - 1] = i`) to match how the writer reads it (`origIndex[k]` for
  `k = 0..npts-1`) — this indexing mismatch was a real bug once (off-by-one, `origIndex[0]`
  uninitialized, occasionally segfaulting far outside the tie-point arrays); watch for
  it if this code is touched again.
- **`tiepoints`** writes two layers (`residuals_noSquint`/`residuals_squint`) when
  `tiePoints->hasSquintSolution` is set (see "Squint" above), one dataset open across
  both `writeTieResidualsLayer()` calls.
- **`rparams`**'s ION_AUTO dual-attempt mode (tries with/without ionosphere correction,
  picks the lower-sigma one) computes residuals for *both* attempts but only writes the
  winning one's GPKG — the write happens at the same point the winning attempt's stdout
  is already selected, not inside `computeRParams()` itself.
- **`-outputFile <path>`**: new on all three, redirects the program's stdout solution
  to a file via the same `dup2`/`STDOUT_FILENO` idiom `rparams.c` already used for its
  ION_AUTO winner-copy — no internal `fprintf(stdout, ...)` call needed to change. Also
  names the debug GPKG (`<path>.residuals.gpkg`); without it, the debug filename falls
  back to `<progname>.<mode>.gpkg` in the cwd. `rparams -runFile` is the one exception:
  it always derives the debug filename from each run's own `outfile`
  (`<outfile>.residuals.gpkg`), since one invocation can process many tie-point
  files/modes — `-outputFile` and `-runFile` together are a CLI error.
- **Makefile**: `writeTieResidualsGpkg.c` added to `common/Makefile`'s `SRCS` and to
  the shared `COMMON=` object list in `mosaicSource/Makefile` (same treatment as
  `geojsonCode.o` — linked into every binary, only actually called from the three that
  use it).

## Wrong-side / footprint-clipping fixes in `llToImageNew.c`

`common/llToImageNew.c`'s `llToImageNew()` (lat/lon/h → range/azimuth) is called by
every program that geocodes points against a SAR frame (tiepoints, rparams, azparams,
mosaic3d, siminsar, coarsereg, getlocc, lltora, computebaseline). Two real bugs here
were found and fixed while investigating spurious/missing tie points on a long,
near-polar Antarctic frame (track spanning ~56,000 azimuth lines / ~160° of longitude
after sub-frame merging — see "Squint" above for other context on such virtual
frames):

- **Left/right (wrong-side) ambiguity**: the zero-Doppler Newton solve
  (`dot(target-satellite, velocity) = 0`) and the range computation (`range =
  |target-satellite|`) are both symmetric under reflection of the target through the
  orbital plane — nothing in the solve itself can tell a real point from its mirror
  image on the wrong side of the ground track. Fixed with a post-solve cross-track
  sign check: `ph = cross(velocity, satellite_unit_position)` is the cross-track
  direction, and for a correctly-illuminated target `dot(target-satellite, ph)` always
  carries the same sign as `inputImage->lookDir` (`LEFT`=-1.0/`RIGHT`=1.0,
  `geocode.h`) — this is the identical sign convention already used to steer the
  forward geolocation in `smlocateZD.c` (`elook = lookDir * acos(...)`) and to flip the
  cross-track baseline component in `svBase.c`'s `svBnBp()`; reused rather than
  re-derived.
- **`checkLL()`'s bounding box was too loose**: this pre-filter (rejects candidate
  points far from the frame before running the expensive Newton iteration) was an
  axis-aligned rectangle — min/max of the 4 corner points, ±50km pad. For a long,
  curving, near-polar frame that box clips along a straight constant-x/y line wherever
  the true (curved/rotated) footprint diverges from the box edge — cutting off real
  swath data at one end while, thanks to the padding, still letting through off-swath
  points elsewhere. Replaced with a real point-in-polygon test against the 4 actual
  corner points (`inputImage->latControlPoints[1..4]`/`lonControlPoints[1..4]`), with
  the same 50km pad now applied as a distance-to-boundary tolerance instead of a box
  expansion. The corners have to be angle-sorted around their centroid before the
  polygon test — `parseInputFile.c`'s `parseControlPointsGeoJson()` index remap
  (`index={0,0,3,1,2}`) was chosen to make the *old* axis-aligned min/max
  order-independent, not to trace the corners in perimeter order.
  - **A 4-corner straight-sided polygon is still only an approximation** of a curved
    orbital footprint — verified adequate on the two Antarctic frames tested (track-64,
    ~56,000-line merged frames), but not stress-tested against ascending passes,
    right-looking data, Northern-hemisphere/Greenland frames, or a track crossing much
    closer to the pole than ~-76°. If clipping/spurious-point symptoms reappear on a
    different geometry, this is the first place to look.
  - The left/right check and the polygon fix are **not fully redundant in theory**
    (a point close to a long polygon edge could in principle be on the wrong side and
    still pass the pad tolerance) even though empirically, on the frame tested,
    disabling the left/right check made no difference once the polygon fix was in
    place. Both are kept since the left/right check is now cheap (see caching below).
- **Performance**: the corner-to-x/y conversion + angle-sort was originally done
  *inside* `checkLL()` — i.e. recomputed on every candidate point (up to 500k times for
  tie points, far more for mosaic3d's per-pixel geocoding), even though the frame's
  footprint polygon is constant. Moved to a one-time computation
  (`computeFootprintPolygon()`) called from `initllToImageNew()` (already invoked
  exactly once per `inputImage` by every caller, before any per-point geocoding calls),
  cached in new `inputImage->footprintX[4]`/`footprintY[4]` fields (`geocode.h`).
  `checkLL()` now just does the query point's own `lltoxy1()` conversion and reads the
  cached array. Verified as a pure performance change (identical point counts/extents
  before and after).

## Single-image fast path, per-row footprint tightening, and DEM cropping (`Mosaic3d/`)

For a genuinely single-track `-rOffsets` run (no phase/Landsat/irregular data merged in, `no3d`/
`noVh` both set), three related changes cut wall-clock time and peak memory roughly in half on a
real cross-Antarctica strip (track-64: 9995×17555 grid, 200m spacing) — 88.6s/12.9GiB →
~44s/7.2GiB, verified pixel-for-pixel against the pre-change baseline.

- **Per-row footprint tightening** (`common/getRegion.c:getRowBounds()`, called only from
  `Mosaic3d/speckleTrackMosaic.c`): `getRegion()`'s single axis-aligned bbox can be far larger
  than a long diagonal track's true footprint (~L²/2 vs true L·W area). `getRowBounds()` reuses
  the same clipped swath polygon (factored into a shared `clipSwathToOutputRect()` helper) to
  compute a tight per-row column range instead. Padding must be at least as generous as
  `checkLL()`'s (`llToImageNew.c`) 50km distance-to-boundary tolerance, dilated in *both* x and
  the row direction (a 2D buffer, not just an x-pad on each row's own crossing) — confirmed
  necessary by two real regressions during verification: an x-only 50km pad still dropped 280k
  valid pixels near the swath's along-track end caps (shallow-angle edges relative to rows) before
  the 2D dilation fixed it to an exact match; an under-sized (~10km) pad dropped 2.5M pixels before
  that. Deliberately scoped to `speckleTrackMosaic.c` only — the crossing-pair loops
  (`make3DMosaic.c`/`make3DOffsets.c`) already bound to the AABB *overlap* of two swaths via
  `getIntersect()`, which is already small/dense for real crossing geometries.
- **Single-image fast path** (`outputImageStructure.singleImageFastPath`, `geocode.h`; set in
  `mosaic3d.c` main before `mallocOutputImage()`, auto-detected, no CLI flag): with exactly one
  contributing image and nothing else touching the grid, the multi-image weighted-accumulation
  buffers (`image`/`image2`/`image3`/`scale`/`scale2`/`scale3`) and per-image weighting scratch
  (`sxTmp`/`syTmp`/`fScale`) are provably redundant — `redoNormalization()`/`endScale()`'s math
  (`common/scalingFunctions.c`) telescopes to `vX_final = vx` (the raw per-pixel value, before the
  `*scX` weighting that only matters for averaging against a second image) and
  `errorX_final = ex`. `speckleTrackMosaic()` writes these final values directly in the per-pixel
  loop instead — 5 of 14 buffers, for the whole run, not just a transient peak. Precondition
  (`mosaic3d.c`, right before `mallocOutputImage()`): `rOffsetFlag && (nAsc+nDesc)==1 &&
  !threeDOffFlag && noVhFlag && no3d && !landSatFile && !irregFile && !makeTies && !statsFlag &&
  !initMapFlag && !timeOverlapFlag` — `no3d` matters because `make3DMosaic()` is called
  unconditionally whenever `(nAsc+nDesc)>0` (not gated by `noVhFlag`) and only skips touching the
  (now-NULL) accumulation buffers via its own early `if (no3d) return;`; found as a real crash risk
  during implementation, not exercised by the verification run itself since `-no3d` was already in
  the test command.
- **DEM cropping** (`mosaic3d.c` main, DEM read before `make3DMosaic`/`speckleTrackMosaic`):
  `readXYDEM()` ("Read a full XY DEM") always passes the `0,0,0,0` "no crop" sentinel to
  `readXYDEMcrop()` — for a circumpolar reference DEM (e.g. the full ~20741×20741 Antarctic
  Copernicus 270m DEM) this loaded ~10x more than a single-track output grid needs, regardless of
  mode (not specific to the fast path above). `mosaic3d.c` now calls `readXYDEMcrop()` directly
  with the output grid's own extent (+10km pad for the interpolation stencil and any DEM/output
  projection-parameter mismatch). `readXYDEMcrop()` itself needed no changes — it already fully
  supported this; `readXYDEM()` just never exercised it. Verified: crop window tightly wraps the
  output grid (7478×13078 vs the full 20741×20741 for track-64), ~2.5GiB additional memory saved.

**Open, unexplained (tiny, deferred) residual**: after the DEM crop fix, 8 of 7,539,859 valid
pixels (track-64, `va` band) differ from the pre-crop baseline by up to 0.0074 (vs ~48 typical
magnitude, <0.02% relative) — everywhere else is an exact match. Ruled out: GDAL read mechanics (a
direct windowed-vs-full-array read of the exact same DEM region is bit-identical, tested directly);
run-to-run nondeterminism (two runs of the identical post-fix binary are bit-identical); float
narrowing in the coordinate math (`xyDEM.x0/y0/deltaX/deltaY`, `interpXYDEM()`,
`computeXYslope()` are `double` throughout — the only float narrowing found,
`common/bilinearInterp.c`'s interpolation weights `t`/`u` and its `float` return value, is
identical code/precision for both the crop and full-DEM read, so it can't explain a
crop-vs-full-specific difference). Leading (unconfirmed) theory: `readXYDEMcrop()`'s cropped
origin `x0Crop = xMinDEM + xOffset*dx` (`getCropBounds()`, `common/readXYDEM.c`) is a real
floating-point add+multiply that rounds relative to the full-DEM case's exact `x0 = xMinDEM`
(`xOffset=0`, no rounding), and at query points where the resulting `xi = (x-x0)/deltaX` lands
very close to an integer, this can flip which DEM pixel pair `interpXYDEM()`'s `floor()` selects —
but a back-of-envelope estimate of double-precision ULP-level rounding suggests this should be far
rarer than 8-in-7.5M, so the mechanism isn't fully confirmed. **If this ever shows up as a larger,
more consequential discrepancy** (not just this ~0.02%-relative curiosity): the fix is to make
`interpXYDEM()` compute `xi` from the DEM's native origin and only subtract the crop's integer
pixel offset afterward (exact, since subtracting an integer from a double doesn't round), rather
than computing it from a pre-shifted (rounded) cropped origin — requires threading the native
origin/offset through `xyDEM` (`geocode.h`), `getCropBounds()`, and
`interpXYDEM()`/`computeXYslope()`. Not done now because `interpXYDEM()` is a shared utility used
by every DEM-height-lookup call across `mosaicSource/` (tiepoints, rparams, azparams, siminsar,
etc.), so the fix's blast radius is much larger than the ~0.02%-relative, 1-in-a-million-pixel
issue it would close.

## Embedded VRT mask band on offset inputs (`-noMask`, `common/readOffsets.c`) — autocleanNISAR

Offset VRTs can carry a GDAL **dataset mask band** (embedded by `autocleanNISAR.py`'s
`range.offsets.good.tif`/`azimuth.offsets.good.tif` products) that flags per-pixel offset
validity. `readGDALOffsets()` (`common/readOffsets.c`) honors it by default:

- After the raster read and the NaN→`-LARGEINT` conversion, if a mask band is present
  (`GDALGetMaskFlags(hBand)` lacks `GMF_ALL_VALID`), it reads the mask via `GDALGetMaskBand()`
  and sets every pixel with `mask==0` to `-LARGEINT`. Downstream code gates validity purely on
  `dr[i][j]/da[i][j] > -LARGEINT`, so masked-out offsets simply drop out — the same path as
  nodata.
- Applied to the **primary value bands only** (`AZIMUTHBUFF`, `RANGEBUFF`,
  `RANGEUSEAZIMUTHBUFF`), not the sigma/error buffers — masking the value band already removes
  the pixel, so the sigma bands don't need it. The partial-read window (`iAzMin..iAzMax`) is
  respected: the mask read uses the same row range.
- `readGDALOffsets()` returns a `maskApplied` flag purely for the `"+mask"` log annotation
  (see the caller at `readOffsets.c:~1214`); it has no effect on control flow.

**`-noMask`** turns this off (`noMask=TRUE` → mask band ignored, all pixels read as-is). The flag
lives in `mosaic3d`, `rparams`, and `azparams` — each program's own CLI parser sets the single
shared global `int32_t noMask` (defined once in `common/getRegion.c`, `extern`-declared in
`readOffsets.c` and the three programs), since all three link the same `readOffsets.c`. Default is
off (mask honored when present). `mosaic3d` logs the flag state (`; noMask Flag : %i`).

## Crossing-orbit-offset slope handling (`make3DOffsets.c`/`make3DMosaic.c`, `common/initRoutines.c`)

**Symptom (2026-07-16):** on a real Antarctic crossing pair (track-24 × track-32), individual
single-track `vr` speeds of ~70-100 m/yr combined into `-3dOff` `vx`/`vy` solutions of
~1000-1200 m/yr — confirmed via `maskTest` (small test area, `/Volumes/insar4/ian/Antarctica2026/maskTest`).
Root cause and fix, in order of investigation:

1. **Not the `speckleTrackMosaic.c`/`makeVhMosaic.c` slope-instability bug found earlier the same
   day** (see git history around this date for the `zSp`-vs-`zWGS84` elevation-validity fix and the
   `vr = (dr*scaleDr + va*cotanpsi*dzda)/(1.0-cotanpsi*dzdr)` unguarded-denominator issue) — those
   two functions aren't even called in a `-3dOff` run (no `-rOffsets`), and `make3DOffsets.c`/
   `make3DMosaic.c` already validated elevation on `zWGS84` correctly.
2. **Real mechanism**: `computeB()` (`common/initRoutines.c`) builds `B[i][j] = dzdx or dzdy /
   tan(psi_i)`. At the blown-up pixels, both `dzdx` and `dzdy` were pinned at `±limitSlope`'s clamp
   with opposite sign (steep, diagonally-trending real terrain — confirmed genuine, not a geocoding
   artifact). This makes each row of `B` sum to zero, which forces each row of `C = I - A·B`
   (`computeVxy()`) to sum to exactly 1, collapsing `detC` to the simple identity `detC = C[0][0] -
   C[1][0]` — structurally fragile for this crossing pair's `A` (from track headings) at the
   original slope clamp of 0.1, landing right at `detC≈0.28`, just above the old `detC<0.25`
   rejection guard.
3. **Why phase-only didn't show it**: `make3DMosaic.c` (phase) and `make3DOffsets.c` (offsets) share
   the *identical* `computeA()`/`computeB()`/`computeVxy()` call — same `A`, `B`, `detC` at the same
   pixel. The ill-conditioning amplifies whatever noise/inconsistency exists between the two
   look-direction measurements, not the signal alone; phase's cm-level precision (vs. offset
   tracking's decimeter-to-meter level, further degraded by the same steep/layover terrain that
   saturates the slope clamp) simply gives the shared amplification much less to work with.
   Quantified below (bias-from-zeroing-slope vs. noise-from-unbounded-slope, swept over slope
   magnitude, using this crossing pair's actual `A` and incidence angles) — bias scales linearly
   with both slope and true velocity, noise amplification is velocity-independent and explodes
   near this geometry's own singularity (`detC=0` at slope≈0.141):

   ![Zero-slope bias vs. unbounded-slope noise amplification](Documents/slopeAmplificationTradeoff.png)
4. **Fix applied** (`common/initRoutines.c`):
   - `limitSlope()` clamp in `computeB()` raised from 0.1 to **0.25** — since `detC` is a purely
     geometric quantity (satellite headings + real terrain slope, no dependence on measurement
     noise), the *right* fix for genuinely steep terrain is to let the real DEM slope through
     rather than truncate it, so the existing `detC` guard can correctly reject pixels whose real
     geometry is singular instead of the old clamp masking that fact. Verified on `maskTest`:
     removes all 10 originally-flagged blown-up pixels; StdDev of `vx` dropped 60.6→39.8, but a new
     (less extreme) worst-case tail appears at different pixels each time the clamp is raised further
     (tried 0.25 and 0.5) — diminishing returns, not a value that makes the failure mode disappear,
     since fragility depends on slope *and* each pixel's own local crossing angle.
   - `detC` rejection threshold in `computeVxy()` raised from 0.25 to **0.5** — combined with the
     slope-clamp change, cuts `maskTest`'s worst-case `vx` from 1082→446 m/yr and max `ex` from
     152→82, at the cost of ~0.6 percentage points of coverage (14.15%→13.50%) — all genuinely
     fragile pixels, not good data.
   - **Ice-shelf carve-out added** to both `make3DOffsets.c` and `make3DMosaic.c`, matching a
     precedent that already existed in `speckleTrackMosaic.c`/`makeVhMosaic.c` (dated 10/13/17,
     "Turn off slope correction for shelf to avoid shelf front or rift artifacts") but had never
     been extended to the crossing-pair path: `if (sMask == SHELF) { B[0][0]=B[0][1]=B[1][0]=B[1][1]=0.0; }`
     right after `computeB()`, leaving `dzdx`/`dzdy` themselves untouched for the `vz` calculation.
     Rationale: true ice-shelf-interior slope is expected to be near-zero, so a large DEM slope
     reading there is almost always a stale/misplaced rift or calving front (DEM vs. current ice
     extent mismatch) rather than real terrain — unclamping (step above) is the right call for real
     grounded steep terrain, but would actively make things worse for this DEM-error case, hence the
     separate mask-driven override.
5. **Open / not yet done**: the ice-shelf carve-out is implemented per the established precedent but
   **not verified against real shelf/rift data** — `maskTest` has no `-shelfMask` flag at all, so
   `sMask` stays `GROUNDED` unconditionally there and the new branch is never exercised. Re-verify
   against a scene with `-shelfMask` and real shelf/calving-front coverage before trusting this in
   production. Also not done: `common/initRoutines.c` is shared across every `mosaicSource/` binary,
   not just `mosaic3d` — only `mosaic3d` has been rebuilt/tested with these changes so far; run
   `csh Makeall` before relying on this elsewhere.

## Crossing-pair error over-counting (`-rhoPhase` / `-rhoOffsets`; `-noPairOverCount` disables)

`make3DMosaic.c`/`make3DOffsets.c` accumulate every crossing pair as an independent observation,
but a pixel covered by n_A ascending and n_D descending images yields only n_A + n_D independent
measurements, not n_A x n_D.

How much the naive `1/sum(w)` understates the variance depends on how CORRELATED the per-pair
errors are, which is what **`rho`** parameterises (`-rhoPhase X` / `-rhoOffsets X`):

```
f(rho) = rho*nPairs + (1 - rho)*(n_A + n_D)/2
```

- `rho = 0` -- purely per-image error, averaging down as the naive count assumes minus the
  `(n_A+n_D)/2` shortfall. Right when random measurement noise dominates.
- `rho = 1` -- error entirely common to every image at the pixel (DEM error, slope error,
  unmodelled vertical rate); such a term does not average at all, so the shortfall is the full
  pair count.
- `rho = 0.5` -- for n_A = n_D = n this equals `(n + n^2)/2`, i.e. **exactly** the legacy
  `-pairCountLegacy` formula, so that historical choice was implicitly a rho = 0.5 claim.

**Both default to 0** (defined in `common/getRegion.c`, passed at `make3DMosaic.c` and
`make3DOffsets.c` respectively), which reproduces the pre-2026-08-27 behaviour numerically.

The value comes from the error budget, not from fitting an observed scatter. `rho` is the share of
per-pixel variance common to EVERY image; the physical candidates -- DEM height and slope, the SMB
grid, tides, geolocation -- are small (the SMB term measures ~0.16 m/yr against a total phase sigma
of ~1.24, i.e. rho ~ 0.02), while the dominant terms (matching noise, per-frame baseline residual,
ionosphere) attach to individual acquisitions and DO average down -- exactly what
`(n_A + n_D)/2` already encodes. Independent confirmation: at rho = 0 the reported error matches
the robust scale of the S1 difference to a few percent (k(MAD) 0.97/1.08 phase, 1.01/1.16
combined), i.e. a budget built from per-observation terms correctly predicts the non-blunder
population it actually models.

**A mixed default (`rhoPhase = 0.5`, `rhoOffsets = 0`) was tried and rejected**, for two reasons:

1. It breaks the combination. See the limitation note below -- the correction inflates the error
   accumulator per round while the weights are untouched, so the coded combined variance
   `(S_p f_p + S_o f_o)/(S_p+S_o)^2` equals the correct `1/(S_p/f_p + S_o/f_o)` only when
   `f_p == f_o`. Mixed rho drives `f_p/f_o` to ~16, past the `> 2 + S_p/S_o` threshold, and
   **2.49% of pixels in a real `both` mosaic came back with a formal error larger than one of the
   inputs that produced it** -- impossible for a genuine inverse-variance combination. Keeping the
   two equal holds `f_p/f_o` at ~1.9, below the `> 2` bar, so the violation is unreachable.
2. Raising rho to ~0.5 does make the reported sigma match the observed *standard* deviation
   (k(std) ~ 1.0), but only by absorbing a blunder population the budget does not model at any rho
   -- unwrapping errors, mismatches, bad frames, which are localised (the worst 1% of crossing-offset
   `vy` pixels hold ~48% of the variance). That is a redundancy parameter used to hide an outlier
   problem. See `Documents/mosaic3d.md` "Interpreting the formal errors" for what `ex`/`ey`
   consequently do and do not mean.

Note the spatial correlation of the error field (rho ~ 0.27 at 170 km) does **not** argue for
rho > 0: it says the field is smooth, not that it fails to average over images. A per-image
ionospheric screen is spatially smooth *and* averages down with more acquisitions.

`-pairCountLegacy` is retained bit-for-bit, ignores `rho`, and is the *asymmetric* relative of
f(rho = 0.5) -- the two coincide only when n_A = n_D. Keep it for reproducing July-vintage
products. Logged per run as `; rhoPhase :` / `; rhoOffsets :`.

**Beware the CLI parser:** `mosaic3d.c` dispatches on `strstr`, and `-rhoOffsets` contains
`offsets`, an existing flag tested 5th in the chain. It escapes only because of the capital O, so
the rho tests are deliberately placed *ahead* of `rOffsets`/`offsets`. Adding two value-taking
flags also required the `argc` cap 44 -> 52 and two more `%s` in `usage()`'s format string.

**Limitation: this corrects the reported error, NOT the weighting.** `inflatePairOverCount()`
writes only to `errorX`/`errorY`, and the velocity is `vXimage/scaleX`, so it cannot move the
velocity for any solution type -- verified by an identical-configuration rerun, which differs from
a legacy-count run by the same ~0.27% of pixels (OpenMP threshold flips), i.e. indistinguishably.
For a COMBINED product the correct treatment would weight each round by `1/f`, giving
`(V_p/f_p + V_o/f_o)/(S_p/f_p + S_o/f_o)`, which differs whenever `f_phase != f_offsets` -- and
with `timePhaseThresh` 10000 against `timeThresh` 37 they differ ~1.7x. Not implemented: with
these formulas it would shift weight *toward* the noisier observable (offsets MAD 2.16 vs phase
1.17). Arguably the separation is a feature -- the relative per-observation sigma drives the
weighting, the absolute scale drives the reported error, and they can be calibrated independently.

`inflatePairOverCount()` (`common/scalingFunctions.c`) applies it. Both callers carry per-pixel
counts through the pair loop -- `nAtmp` (distinct outer-loop images, deduped via the `aContrib`
byte map folded in at the end of each outer iteration) and `nDtmp` (contributing **pairs**) --
snapshot `errorX`/`errorY` into `errorX0`/`errorY0` right after `undoNormalization()`, and call
the helper before `endScale()`. `-noPairOverCount` skips all of it (buffers not even allocated).
Every mosaic log records the state as `; pairOverCount : 0|1`.

**Two invariants, both of which an earlier always-on version violated:**

1. **Inflate only this round's own contribution, and before `endScale()`.** `undoNormalization()`
   is the exact inverse of `endScale()` for `errorX` (`*= scaleX^2` vs `/= scaleX^2`), so a
   multiplier applied *after* `endScale()` to the whole buffer is folded back into the accumulator
   and inflated again by the next round's factor. The buffers are cumulative across rounds
   (Landsat -> make3DMosaic -> make3DOffsets -> makeVhMosaic -> speckleTrackMosaic), so "the whole
   buffer" is never just this round's data. **This was the bug** behind the reported symptom:
   adding crossing offsets to a phase mosaic made errors grow 2.4x while the GNSS velocity
   accuracy simultaneously improved.
2. **`nDtmp` is a pair count, not n_D.** The helper derives `n_D = nPairs / nOuter`, exact for the
   standard setup (`consolodateLists()` orders all-ascending-then-all-descending; `sepAscDesc`
   defaults TRUE, so outer-loop contributors are exactly the ascending images). Using
   `nAtmp + nDtmp` gives n_A + n_A*n_D -- for a Greenland sector with n~28 per direction that is
   ~423 rather than ~28, i.e. sigma too large by ~3.9x.

Velocities are unaffected either way -- only the error accumulator is scaled (verified
max|dvx| = 0).

### Calibration evidence (2026-08-25) -- and a dating trap

Measured, sector 003.003 of the Greenland Release phase-only mosaic, median `ex`:

| variant | median ex |
|---|---|
| no correction at all | **0.136** |
| Jul 21 validated product (old n_A + n_A*n_D count) | **2.803** |
| (n_A+n_D)/2, this implementation | ~0.73 (inferred; ratio implies n ~ 28.6) |

At the five NIU interior-Greenland GNSS sites
(`Notebooks/calValReport/preliminaryAssessment`, Figure 4; two vintages in `figures/` Jul 22 and
`figures2026-08-25/`), the **vy component is the only clean read on random error** -- its
residuals are small and mixed-sign (+0.9, -0.7, -0.3, +0.3, +1.0), MAE **0.65 m/yr**. Against
that: the old count gave a formal error of **3.48** (5.4x too large), and no correction at all
implies roughly **0.17** (~4x too small). The derived (n_A+n_D)/2 falls between and matches the
observed scatter for n ~ 56 at those heavily-overlapped interior sites. vx and vv residuals are
all one sign (common-mode bias ~ +3 / -3 m/yr) and so cannot calibrate a random-error bar.

**Dating trap -- do not repeat this.** These products carry no record of which error calibration
produced them (the `; pairOverCount :` log line was added precisely to fix that). An attempt to
date the correction's introduction from output statistics alone concluded, wrongly, that the July
products predated it -- reasoning from all-squint Jul 7 (2.859) vs Aug 25 (12.677) being a 4.37x
jump in identical configuration. The real explanation is that the phase-side correction was
present in both, and what changed between Jul 21 and Aug 25 was on the offsets side / the
compounding. The only reliable test is to **run the uncorrected code on current data and compare**
(0.136 vs 2.803 settles it in one run); an inference from ratios between products of different
vintages does not.

**Still open:** the correction's magnitude is only as good as the input sigmas it scales.
`phaseError = sqrt(sig2Base + sigma^2)` (`computePhiZM3d`) is built from per-image systematics
(baseline covariance, tie-point residual). Calibrating those against the GNSS scatter directly --
rather than tuning the pair-count factor to compensate -- is the cleaner fix. The theoretically
exact treatment (average ascending and descending measurements per pixel, then solve once) would
change velocities too and needs the inner loop restructured; not done.

## Range-offset error budget: the missing accuracy term (`-noRSigmaResidual`, ON by default)

The range-offset formal error was ~6x too small, and — the part that actually corrupted
products — *smaller than the phase formal error*, though phase is far the more precise
observable. That inverted the inverse-variance weighting in `computeVxy`, so crossing offsets
dominated the combined solution and dragged it.

### The mechanism

The two error budgets were not analogous:

- Phase (`make3DMosaic.c:639,703`):
  `phaseError = sqrt(sig2Base + min(6*PI, vhParam->sigma)^2)` — `vhParam->sigma` is the tie-point
  fit residual, sampled across the whole frame, so it sees long-wavelength error.
- Offsets (`make3DOffsets.c`, `speckleTrackMosaic.c`):
  `sigmaR = sqrt(interpRangeSigma^2 + demError^2 + sig2Base)`. `interpRangeSigma`
  (`common/interpOffsets.c:169-190`) reads only the `.sr` band, which Cullst computes as a LOCAL
  neighbourhood scatter **after removing a local plane** (`Cullst/cullStats.c:78-90`) and divides
  by `sqrt(nLooksEff)` (`cullSmooth.c:94-95`). It measures matching noise and is **structurally
  blind to long-wavelength error** — ionosphere, orbit ramps — at any amplitude.

`offsets.sigmaRresidual` — the rparams tie-point fit residual in metres, the exact analogue of
`vhParam->sigma`, written by `rParams/computeRParams.c:458` — existed, was parsed
(`readOffsets.c:239`) and printed (`:280`), and was used **only** as the `< 0` no-solution
sentinel (`make3DOffsets.c:174,243`). It never reached the error line.

Compounding it, the ionospheric correction applied at `make3DOffsets.c:326,337` carries **no
uncertainty at all**: `offsetCorrection` (`common.h`) has no sigma field. `sigmaRresidual` covers
that empirically, since rparams fits on ion-corrected offsets (`rParams/getROffsets.c:91` applies
the correction, `:96` converts to metres, then `computeRParams` forms `sigP` from that array) and
`checkForIonosphereCorrection()` (`readOffsets.c:941-980`) hard-errors unless mosaic3d applies the
same correction file rparams fitted against. Two caveats: `meanP` is subtracted
(`computeRParams.c:295`), so a pure constant range bias is excluded by construction; and the
baseline polynomial absorbs smooth ionospheric ramps, so `sigP` carries only the ionosphere's
*departure* from that polynomial (the absorbed part is corrected, and its uncertainty is in
`sig2Base`).

### The fix

`sigmaRresidual` added in quadrature at `make3DOffsets.c` (ascending and descending) and
`speckleTrackMosaic.c`, via `rangeAccuracyVar()` (`common/interpOffsets.c`). Both terms are metres
of slant range at that point — `interpRangeSigma` has already been multiplied by `rSLPixSize` — so
they combine directly. The `< 0` sentinel is filtered upstream, so it can never enter the sum.
Default ON; `-noRSigmaResidual` restores the old budget. Logged per run as `; rSigmaResidual :`.

**Why it must be in the WEIGHT, not only the reported error.** A tempting variant is to weight on
`.sr` alone and add the residual only to the reported `ex`/`ey`, which would preserve `.sr`'s
spatial discrimination. That is wrong, and dangerously so: `.sr` cannot see ionosphere by
construction, so an ionosphere-contaminated frame still reports a *small* `.sr` and would be given
*high* weight. The worst frames would dominate the average. The residual is the only signal the
code has that a frame is contaminated, so it has to drive the weighting.

**The division of labour is sensor-dependent, and that is the point.** The term adapts on its own,
which a fixed constant (`-rSigmaConst`, diagnostic only) cannot:

| sensor | residual | behaviour |
|---|---|---|
| NISAR | ~0.22 m, ~5x `.sr` | ionosphere dominates; budget becomes essentially per-frame, and noisy frames are down-weighted wholesale — correct, the data really is frame-limited |
| TSX | X-band, negligible ionosphere | `.sr` captures essentially all the noise and should mirror the residual; worst case the two double-count, i.e. sigma over-estimated by ~sqrt(2) |
| S1 | intermediate | a mix of the two regimes |

So weighting is per-frame where the frame-level error dominates and per-pixel where it does not,
sliding between them automatically. A sqrt(2) over-estimate in the X-band limit is an acceptable
price.

**Expected side effect, not a defect.** `.sr` varies spatially within a frame (measured p90/p10
~7x in sigma, ~48x in weight). Adding a term much larger than it flattens that: 99% of the spatial
weighting is lost for NISAR, 76% for an S1-like 0.04 m residual. Where the frame-level error
genuinely dominates, flat weighting IS the right answer.

**Open, unexplained.** Against the S1 reference over stable ground the offsets-ONLY product got
worse with the term (difference MAD 3.98 -> 7.58 m/yr) even as its calibration improved
(k 6.4 -> 3.5), while the COMBINED product improved as expected (vy MAD 2.26 -> 1.48, k 4.35 ->
1.75). If the residual tracks contamination as intended, the offsets-only velocity should have
improved too. Worth understanding before this is trusted quantitatively — candidates are that the
residual is a single scalar summarising a spatially varying ionospheric field, or that the
degradation is concentrated where down-weighted frames were the only coverage.

Note `sigP` is raw `sqrt(varP - meanP^2)` with **no** `sqrt(chi2/n)` scaling despite the legacy
label at `computeRParams.c:493` (see the comment at `:453`), and it is per-frame, not spatial.

### Evidence (2026-08-26, Sentinel-1 comparison over stable ground, S1 speed < 50 m/yr)

Built with `nisarErrors`' `makeS1Reference` / `compareToS1` (see that package's Documents).
`k` = robust scale of the S1 difference / median formal error; 1.0 is calibrated.

| product | diff MAD | formal median | k | \|z\|>3 |
|---|---|---|---|---|
| phase | 1.18 | 0.863 | 1.37 | 9.6% |
| offsets | 3.27 | 0.568 | 5.75 | 55.5% |
| both | 1.52 | 0.530 | 2.88 | 33.3% |

Three independent lines of evidence identified the mechanism:

1. **kDecile.** `k` within deciles of the formal error declined monotonically for offsets
   (6.29 -> 3.82) while phase's was flat near 1. A multiplicative scale error gives flat `k`; a
   **missing term added in quadrature** gives exactly this decline. So rescaling `.sr` would have
   been the wrong shape of fix.
2. **Robust variogram.** The offsets difference grows 0.89 -> 3.16 m/yr between 1.6 and 51 km
   without saturating, and the total MAD (3.82) exceeds even the 51 km value. The error is
   long-correlation, which `.sr` cannot represent.
3. **Short-lag comparison.** The offsets formal error (0.568) matches the 1.6 km difference (0.89)
   to within 35% — it is a *correct* estimate of local matching noise and nothing else. Phase
   shows the mirror image: formal 0.863 against a 1.6 km difference of 0.21, i.e. 4x larger,
   because its budget already contains the frame-scale residual.

Averaging invariance confirms it: a natively-coarse 1600 m mosaic and an 8x8-averaged 200 m one
agree to within 5% at every lag. White noise would have been suppressed up to 8x by the averaging.

**Cost of the mis-weighting**, on matched pixels (`--commonMask`): the combined product was 29%
worse than phase alone in vx and 55% worse in vy, while reporting the lowest formal error of the
three. Correct weighting of true errors 1.18 and 3.27 gives 1.11 — slightly better than phase.

**Predicted effect**: measured over the 60 sampled frames, `sigmaRresidual` (median 0.2275 m) is
5.6x the `.sr` term (median 0.0416 m; note `rSLPixSize` is ~1.56 m for these products, not the
~7 m an ML pixel might suggest), so the quadrature multiplier is ~5.6x — taking offsets `k` from
5.81 to ~1.03 and putting the offsets formal error ~4x *above* phase's, so phase takes ~93% of
the weight where both exist. All 595 frames in the master input name an `rBaseline.deltabp.yaml`
and none is missing; 591 have a usable sigma and the 4 with the `< 0` sentinel were already
skipped.

**Still open**: the offsets carry a heavy outlier tail independent of this (std/MAD ~2.5-3, and
the tail worsens under point sampling), which is a localised-bad-areas problem rather than a
global scaling one and is not addressed here.

## Other contents

- **`Documents/`** — per-program design notes; check here before `Documents/` in the GIT64 root.
- **`BUGFIXES_*.md`** — dated records of fixes applied from past automated code reviews
  (mosaic3d/makeVhMosaic, azParams, geoMosaic, rParams, simInSAR — all 2026-05). Historical record,
  not living docs; safe to ignore unless debugging a regression in one of those files.
- **`cloneAll`** — legacy csh script to `git clone` sibling repos (clib, cRecipes, fft,
  speckleSource, etc.) when this was a multi-repo checkout. Not needed in the GIT64 monorepo.
- **`baseline.26x16.yaml`** — sample/test config for baseline computation.
