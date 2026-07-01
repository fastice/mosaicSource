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

## Other contents

- **`Documents/`** — per-program design notes; check here before `Documents/` in the GIT64 root.
- **`BUGFIXES_*.md`** — dated records of fixes applied from past automated code reviews
  (mosaic3d/makeVhMosaic, azParams, geoMosaic, rParams, simInSAR — all 2026-05). Historical record,
  not living docs; safe to ignore unless debugging a regression in one of those files.
- **`cloneAll`** — legacy csh script to `git clone` sibling repos (clib, cRecipes, fft,
  speckleSource, etc.) when this was a multi-repo checkout. Not needed in the GIT64 monorepo.
- **`baseline.26x16.yaml`** — sample/test config for baseline computation.
