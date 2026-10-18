# Per-sensor row switches for the hopper — design notes (not implemented)

Status: ideas only, written 2026-10-13 at the user's request ("not now, but in the near future").
**Update 2026-10-14:** the azimuth part is implemented as the flag form
(`-noAzimuthRowsSensors`, sensor inferred from the geodat label -- option 2 of section 1), driven by
`excludeRows: {azimuth: [...]}` in project.yaml/template; see mosaic3d.md "Row selection". The
`-rowConfig` yaml and the other per-sensor parameters below remain ideas.
Motivation: the experimental 2026 annual merges S1/TSX and NISAR into one hopper solve. The row
switches (`-noPhaseRows`, `-noRangeRows`, `-noAzimuthRows`) and the azimuth residual cut
(`-sigmaAThreshVel`) are global, but the right choice differs by sensor:

| sensor | azimuth rows | evidence |
|---|---|---|
| S1 Greenland | ON | collapse a spurious fast tail; with the 2026-10-11 az.est sigma fix they help NIU 2023-24 (vy 3.92 -> 1.89) |
| NISAR Greenland | OFF in production (`-noAzimuthRows`) | ~0.1 % of the weight; LEFT-look sign history |
| NISAR / S1 Antarctica fast ice | OFF | azimuth inflated speeds (median 3,783 m/yr above 5 km/yr, Thwaites deck) |

## 1. Knowing each image's sensor

Nothing in `inputImageStructure` records the sensor today. Options, best first:

1. **Tag it in the master.** `makeMaster` already knows each line's source (`masterInput.yaml`
   `sources:` -> S1-insar9, TSX, phase, ...). Add a `sensor:` key per source and write it as an
   OPTIONAL trailing token on the input line (e.g. `sensor=NISAR`). `parseInputFile` reads it when
   present; old masters (no token) fall through to (2). This is explicit and survives merged
   masters like `experimentalAnnual/mergeMasters.py`.
2. **Infer from the geodat** as a fallback: NISAR geojson `ImageName` starts `NISAR_`; S1 `.in`
   `; Image name: s1[abcd]-...`; TSX names likewise. Store an enum in `inputImageStructure`
   (`SENSOR_S1`, `SENSOR_NISAR`, `SENSOR_TSX`, `SENSOR_UNKNOWN`).

Unknown sensor -> global switches only (exactly today's behaviour).

## 2. Interface

A small yaml passed as `-rowConfig rows.yaml`, keyed by sensor, so it can grow without new flags:

```yaml
default: {phase: true, range: true, azimuth: true}
NISAR:   {azimuth: false}
S1:      {sigmaAThreshVel: 150}
TSX:     {}
# later, same file: per-sensor sigma scales (see 4)
```

Rules: a global `-no*Rows` still wins (backwards compatible, keeps every existing template
working); otherwise sensor entry > `default` > compiled-in. Log the resolved table and the per-sensor
row counts actually used (`; rows NISAR: phase n, range n, azimuth 0`).

A flag form (`-noAzimuthRowsFor NISAR`) is simpler for one switch but does not scale to the
per-sensor parameters in 4.

## 3. Where it goes in the code

The hopper already decides rows PER IMAGE, so this is a narrow change:
`mosaicHopper3D.c` ~356-375 (and `mosaicHopper.c` ~310-325) set `hasPhase`, `useRangeRow`,
`useAzRow` per frame. AND each with the image's resolved sensor switch, and skip the raster read
when a row type is off for that frame (as `-noAzimuthRows` already does globally — gated on the
hopper running, since `makeVhMosaic` under `-legacyCode` would consume an unread buffer).
The legacy rounds need not support it.

## 4. Natural extensions (same config)

- **Per-sensor sigma scale** — the "tie-accuracy scale" idea: a multiplicative calibration on each
  row type's sigma (e.g. S1 phase residual includes the offsets residual; NISAR does not). Could
  default from the baseline file (`tieSigmaPhase`) and be overridden here.
- **Per-sensor `sigmaAThreshVel`** and ionosphere options (`useAzIonosphere` is NISAR-only today).
- **Per-sensor date windows**, replacing the ad hoc master filtering used for the experimental
  annual (NISAR only from 2026-06-17).

## 5. Checks when implemented

- No `-rowConfig`: byte-identical output to the current binary.
- A config that matches the global flags reproduces them exactly (e.g. `default: {azimuth:false}`
  == `-noAzimuthRows`).
- obsDump row census per sensor confirms the switches took effect.
