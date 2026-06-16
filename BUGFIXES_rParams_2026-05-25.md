# Bug fixes applied to rParams 2026-05-25

Applied findings from automated code review of mosaicSource/rParams/.
All changes compile clean with no warnings.

---

## computeRParams.c

### 1. computeLinearBaseline writes dBp into dBpQorig, then immediately zeroes it
**Problem:** `computeLinearBaseline` was passed `&(tiePoints->dBpQorig)` as its fourth
output (`*dbp`), so the linear-Bp slope was written to the quadratic field. The very next
line zeroed `dBpQorig`, discarding the computed value. Meanwhile `tiePoints->dBporig` was
never assigned in this code path (`initWithSV == TRUE`), so the read at line 169
(`dBp = tiePoints->dBporig`) consumed garbage — corrupting the initial baseline estimate
for the entire run when no baseline file is present.

**Fix:** Changed the fourth argument from `&(tiePoints->dBpQorig)` to `&(tiePoints->dBporig)`.
Also fixed the subsequent `fprintf` to print `dBporig` (the actual dBp) instead of
`dBpQorig` (which is now always 0 at that point).

**To revert:** Change `&(tiePoints->dBporig)` back to `&(tiePoints->dBpQorig)` and
restore `dBpQorig` in the fprintf.

---

### 2. bpdBpFlag output: constant term written to dBpQ column instead of const column
**Problem:** With `bpdBpFlag == TRUE`, the SVD fit has two free parameters: `a[1]` (dBp
slope) and `a[2]` (constant offset). The baseline file column order is
`Bn Bp dBn dBp const dBnQ dBpQ`. The fprintf placed `a[2]` in the seventh column (dBpQ)
with the fifth (const) hardcoded to `0.0`:
```c
fprintf(stdout, "...  %f %f %f %f\n&\n", ..., a[1], 0.0, 0.0, a[2]);
//                                              dBp  const dBnQ  dBpQ ← a[2] in wrong column
```
Any downstream tool reading the baseline file would pick up zero for the constant and
garbage for dBpQ.

**Fix:** Swapped `a[2]` and `0.0` to put the constant in column 5 and zero in columns 6–7:
```c
fprintf(stdout, "...  %f %f %f %f\n&\n", ..., a[1], a[2], 0.0, 0.0);
```

**To revert:** Swap back to `a[1], 0.0, 0.0, a[2]`.

---

### 3. nParams declared double, used as integer count
**Problem:** `double nParams` was passed to `dvector(1, nParams)` and `dmatrix(1, npts, 1, nParams)`
(which take `int32_t`), compared with integer `ma`, and used in `if (i1 < nParams)`. The
implicit `double → int32_t` truncation was numerically harmless for the values used (1–6)
but generated type-mismatch warnings and obscured intent. The `fprintf` on line 105 also
printed it with `%f`.

**Fix:** Changed declaration to `int32_t nParams`; changed `%f` to `%i` in the fprintf.

**To revert:** Change `int32_t nParams` back to `double nParams`; restore `%f`.

---

### 4. Insufficient-points guard uses uninitialized i1 and appears after division by npts
**Problem:** `i1` is only assigned inside the valid-point loop body. If `npts == 0`,
`i1` is never assigned, and `if (i1 < nParams)` reads garbage — undefined behaviour.
Additionally, the guard appeared after the divisions `varP / (double)npts` and
`meanP / (double)npts`, so a zero-npts case would divide by zero before the guard
could fire.

**Fix:** Moved the guard to before the divisions, and replaced `i1` with `npts`
(the pre-computed valid-point count, always valid):
```c
if (npts < nParams) { fprintf(...); fewPoints(); }
varP  = varP  / (double)npts;
meanP = meanP / (double)npts;
sigP  = sqrt(max(0.0, varP - meanP * meanP));
```
Also guarded `sqrt` against floating-point cancellation producing a tiny negative argument.

**To revert:** Move the guard back after the sqrt; change `npts` back to `i1`;
remove `max(0.0, ...)`.

---

## getBaselineFile.c

### 5. BpvC / dBpv uninitialized when dBpFlag == FALSE
**Problem:** `BpvC` and `dBpv` were only assigned inside
`if (tiePoints->dBpFlag == TRUE) { ... }`. Immediately after, `tiePoints->BpCorig = BpvC`
and `tiePoints->dBporig = dBpv` always executed unconditionally. If `dBpFlag == FALSE`,
both struct fields are written with uninitialized stack values. Also, `Bpv = Bp1 + Bp2`
was assigned immediately after the `if` block but never read — dead code.

**Fix:** Moved `BpvC`, `dBpv`, and `dBpQv` assignments unconditionally outside the `if`
block; removed the dead `Bpv` line.

**To revert:** Restore the `if (tiePoints->dBpFlag == TRUE)` guard around the three Bp
assignments, and restore `Bpv = Bp1 + Bp2` after the block.

---

## addOffsetCorrections.c

### 6. inputImage->lookDir uses pointer syntax for a pass-by-value parameter
**Problem:** The function signature is `addOffsetCorrections(inputImageStructure inputImage, ...)`
— `inputImage` is passed by value. Line 63 used arrow notation:
```c
if (inputImage->lookDir == LEFT)
```
This is a type error (dereferencing a struct value as a pointer). The function is not
currently called anywhere (`rparams.c` does not call it), so there is no runtime effect,
but it would prevent compilation of any future caller.

**Fix:** Changed to `if (inputImage.lookDir == LEFT)`.

**To revert:** Change `.lookDir` back to `->lookDir`.

---

## rparams.c

### 7. mkstemp return value unchecked — crash if temp file creation fails
**Problem:** `mkstemp(tmp1)` and `mkstemp(tmp2)` each return `-1` on failure. The
return values were not checked. If either call fails, `dup2(-1, STDOUT_FILENO)` silently
fails, output from both runs is interleaved on real stdout, and the subsequent `fopen`
on the unmodified template string (with literal `X`es) returns `NULL`. Reading from a
NULL FILE pointer is undefined behaviour — typically a crash.

**Fix:** Added:
```c
if (fd1 < 0 || fd2 < 0) error("rparams: mkstemp failed\n");
```
immediately after both `mkstemp` calls.

**To revert:** Remove the `if (fd1 < 0 || fd2 < 0)` check.
