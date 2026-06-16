# Bug fixes applied to azParams 2026-05-25

Applied findings from automated code review of mosaicSource/azParams/.
All changes compile clean with no warnings.

---

## azparams.c

### 1. constOnlyFlag / linFlag read from local (uninitialized) variables
**Problem:** `readArgs` parses `-constOnly` and `-linear` into local variables
`constOnlyFlag` and `linFlag`, then stores them in `tiePoints->constOnlyFlag` and
`tiePoints->linFlag`. The diagnostic `fprintf` calls at the end of `main` printed
the local variables — which are out of scope and were already assigned to the struct
at that point, but the real issue was that the original code tested the local variables
directly:
```c
if (constOnlyFlag == TRUE)
    fprintf(stderr, "\n(****Constant only fit*****\n");
if (linFlag == TRUE)
    fprintf(stderr, "\n(****Including linear term fit*****\n");
```
These two locals are declared inside `readArgs` and not accessible in `main`. The
original code used the same names in `main` before the `readArgs` call, so the
diagnostics always printed from uninitialized stack values.

**Fix:** Changed both checks to read from the struct:
```c
if (tiePoints.constOnlyFlag == TRUE)
if (tiePoints.linFlag == TRUE)
```

**To revert:** Change `tiePoints.constOnlyFlag` back to `constOnlyFlag` and
`tiePoints.linFlag` back to `linFlag` in the two `fprintf` calls in `main`.

---

## computeAzparams.c

### 2. azimuth uninitialized before first getReH() call
**Problem:** The sequence in `computeAzParams` was:
```c
double azimuth, ReH, ...;
...
ReH = getReH(cP, inputImage, azimuth);   // azimuth is uninitialized garbage here
...
azimuth = inputImage->azimuthSize / 2.0;
ReH = getReH(cP, inputImage, azimuth);   // correct
thetaC = thetaRReZReH(cP->RCenter, (Re + 0), (ReH));
```
The first call to `getReH` used an uninitialized `azimuth`, so `ReH` was garbage.
The result of that call was then overwritten by the second call, so `thetaC` itself
was computed correctly — but the dead first call still reads undefined memory and
would exhibit undefined behaviour if `azimuth` happened to be out of range.

**Fix:** Added `azimuth = inputImage->azimuthSize / 2.0;` immediately before the
first `getReH` call, and kept the second assignment (which was already correct):
```c
azimuth = inputImage->azimuthSize / 2.0;
ReH = getReH(cP, inputImage, azimuth);
...
azimuth = inputImage->azimuthSize / 2.0;
ReH = getReH(cP, inputImage, azimuth);
thetaC = thetaRReZReH(cP->RCenter, (Re + 0), (ReH));
```

**To revert:** Remove the first `azimuth = inputImage->azimuthSize / 2.0;` assignment
(the one before the first `getReH` call).

---

### 3. lineCount uninitialized in getBaselineRates
**Problem:** `getBaselineRates` declared `int32_t lineCount, eod;` without
initialising `lineCount`. The very first call is
`lineCount = getDataString(fp, lineCount, ...)`, which passes the uninitialized
value as the input line counter. `getDataString` uses that counter only for
error messages, so the practical impact is a garbage line number in any error
message printed for a malformed baseline file. Not a correctness issue at runtime
but undefined behaviour per the C standard.

**Fix:** Changed declaration to `int32_t lineCount = 0, eod;`.

**To revert:** Remove the `= 0` initialiser from `lineCount`.

---

### 4. Insufficient-points guard uses uninitialized i1; appears after division by npts
**Problem:** Mirroring the same bug found in `computeRParams.c`. `i1` is the
per-iteration inner-loop index and is only ever assigned inside the valid-point
loop body. If `npts == 0`, `i1` is never assigned and `if (i1 < ma)` reads garbage
— undefined behaviour. Additionally the guard appeared **after** the divisions
`varP / (double)npts` and `meanP / (double)npts`, so a zero-npts case would divide
by zero before the guard could fire:
```c
varP  = varP  / npts;
meanP = meanP / npts;
sigP  = sqrt(varP - meanP * meanP);
if (i1 < ma) error(...);      // too late; also i1 is uninitialized
```

**Fix:** Moved the guard to before the divisions, replaced `i1` with `npts`
(the pre-computed valid-point count, always valid), and guarded `sqrt` against
floating-point cancellation producing a tiny negative argument:
```c
if (npts < ma)
    error("azparams: Insufficient Number (%i) of Valid tie points \n", npts);
varP  = varP  / npts;
meanP = meanP / npts;
sigP  = sqrt(max(0.0, varP - meanP * meanP));
```

**To revert:** Move the guard back after the `sqrt`; change `npts` back to `i1`;
remove `max(0.0, ...)`.

---

### 5. Division by weightSum without zero-check
**Problem:** After accumulating `weightSum` over all valid tie points, the code
used it immediately as `npts / weightSum` inside the `sig[i1]` computation:
```c
sig[i1] = max(sigP * tiePoints->weight[i] * npts / weightSum, 0.1 * sigP);
```
If all tie-point weights are zero (a configuration error, but one that can arise
from a corrupt tiepoint file), `weightSum` is zero and the division produces `Inf`
or `NaN`, silently corrupting all sigma values and the SVD solution.

**Fix:** Added an explicit guard after the accumulation loop:
```c
if (weightSum == 0.0)
    error("azparams: all tiepoint weights are zero\n");
```

**To revert:** Remove the `if (weightSum == 0.0)` check.
