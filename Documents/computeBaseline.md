# computeBaseline — Interferometric Baseline Estimator

## Purpose

Computes the interferometric baseline between two SAR acquisitions directly from their
satellite state vectors, without requiring tiepoints or phase data. Outputs a baseline
parameter file suitable for use by `rparams`, `azparams`, and `mosaic3d`.

---

## Usage

```
computeBaseline geodat1 geodat2 baselineFile
```

### Required Arguments

| Argument        | Description |
|-----------------|-------------|
| `geodat1`       | Geodat file for the first (reference) SAR image |
| `geodat2`       | Geodat file for the second SAR image |
| `baselineFile`  | Output baseline parameter file |

---

## Algorithm

Both geodat files are parsed and their state vectors initialised via `parseInputFile` and
`initllToImageNew`. The Earth radius and satellite radius at image centre are computed:

$$
R_e = \text{earthRadius}(\varphi_c), \qquad
R_{eH} = \text{getReH}(\text{cP},\, \text{image},\, N_{az}/2)
$$

The centre look angle $\theta_c$ is derived from the centre slant range $R_c$:

$$
\theta_c = \text{thetaRReZReH}(R_c,\, R_e,\, R_{eH})
$$

Baselines are evaluated at the start and end azimuth times of image 1 using
`svBnBp` (which calls `svBaseTCN` to form the TCN baseline vector and then projects
it onto the normal/parallel frame at look angle $\theta_c$):

$$
B_n^{(1)},\, B_p^{(1)} \quad\text{at azimuth start time } t_1
$$
$$
B_n^{(2)},\, B_p^{(2)} \quad\text{at azimuth end time } t_2
$$

The output parameters are:

$$
B_n = \tfrac{1}{2}(B_n^{(1)} + B_n^{(2)}), \qquad
B_p = \tfrac{1}{2}(B_p^{(1)} + B_p^{(2)})
$$
$$
\delta B_n = B_n^{(2)} - B_n^{(1)}, \qquad
\delta B_p = B_p^{(2)} - B_p^{(1)}
$$

The along-track quadratic terms $\delta B_{nQ}$, $\delta B_{pQ}$ and the constant range
bias are set to zero (written as `0.0`).

---

## Output File Format

A single ASCII line:

```
Bn  Bp  dBn  0.0  dbp
```

This is a subset of the full 7-parameter baseline format used by `rparams` output
(`Bn Bp dBn dBp const dBnQ dBpQ`). The quadratic terms and constant bias are absent
(effectively zero) since they cannot be estimated without tiepoints.

---

## Dependencies

| Function | Source | Purpose |
|----------|--------|---------|
| `parseInputFile` | `common/parseInputFile.c` | Read geodat and populate `inputImageStructure` |
| `initllToImageNew` | `common/llToImageNew.c` | Initialise state-vector geometry |
| `getReH` | `common/initRoutines.c` | Azimuth-varying satellite radius |
| `thetaRReZReH` | `common/initRoutines.c` | Centre look angle from ranges |
| `svBnBp` | `common/svBase.c` | Normal/parallel baseline from state vectors |
| `svBaseTCN` | `common/svBase.c` | Baseline in Track/Cross-track/Normal frame |
| `polintVec` | `common/polintVec.c` | Polynomial state-vector interpolation |
