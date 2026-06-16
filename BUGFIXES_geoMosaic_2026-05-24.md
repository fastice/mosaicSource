# Bug fixes applied to geoMosaic 2026-05-24

Applied findings from automated code review of mosaicSource/geoMosaic/.
All changes compile clean with no warnings.

---

## makeGeoMosaic.c

### 1. Uninitialised `i` used as array index in getGeoMosaicImage()
**Problem:** `int32_t k1, k2, i;` — `i` declared but never initialised.
`smoothImage(&(inputImage[i]), smoothL)` used it as a subscript into the `inputImage`
array, but `inputImage` is a pointer to a single struct (the parameter). With an
uninitialised `i` (typically non-zero stack garbage) this either crashes or silently
calls `smoothImage` on the wrong struct.

**Fix:** Changed to `smoothImage(inputImage, smoothL)` — direct pointer, no subscript.

**To revert:** Change back to `smoothImage(&(inputImage[i]), smoothL)`.

---

### 2. ASCENDING orbit-priority: first pixel never written
**Problem:** `scale` is initialised to `−1` for all pixels when `orbitPriority >= 0`
(line ~792). The ASCENDING write-condition was:
```c
if (scale[i1][j1] > 0.1 || scale[i1][j1] == DESCENDING)
```
With ASCENDING=1, DESCENDING=0: an uninitialised pixel has `scale == −1`, which satisfies
neither `> 0.1` nor `== 0`. So no pixel is ever written with `-ascending`. By contrast the
DESCENDING branch correctly uses `< −0.1` to catch the no-data case.

**Fix:** Changed `> 0.1` to `< -0.1`:
```c
if (scale[i1][j1] < -0.1 || scale[i1][j1] == DESCENDING)
```

**To revert:** Change `< -0.1` back to `> 0.1`.

---

### 3. `sum` printed uninitialised in smoothImage()
**Problem:** `float sum, sumF;` — `sum` is declared without initialisation.
`fprintf(stderr, "hw = %i %f\n", hw, sum)` at line ~973 prints it before it is first
assigned (which happens inside the smoothing loop at line ~986). Printed value is garbage.

**Fix:** Changed declaration to `float sum = 0.0, sumF;`

**To revert:** Remove the `= 0.0` initialiser.

---

### 4. `filt` base pointer lost; memory leaked every call to smoothImage()
**Problem:**
```c
filt = (float *)malloc((size_t)((hw * 2 + 1) * sizeof(float)));
filt = &(filt[hw]);   // original malloc pointer overwritten — leaked
```
The pointer returned by `malloc` is immediately discarded. `smoothImage` has no `free`
call, so every invocation (once per input image) leaks `(2*hw+1) * 4` bytes.

**Fix:** Kept the original pointer in `filtBase`:
```c
float *filtBase = (float *)malloc((size_t)((hw * 2 + 1) * sizeof(float)));
filt = &(filtBase[hw]);
```
Added `free(filtBase);` at the end of the function (before the closing `}`).

**To revert:** Remove `filtBase`; restore `filt = (float *)malloc(...); filt = &(filt[hw]);`;
remove the `free(filtBase)` line.

---

### 5. Division by zero: sin(psiE * DTOR) can be 0 in gamma correction
**Problem:**
```c
gBufTmp[i1][j1] = 10.0 * log10((AbCum / AgCum) / sin(psiE * DTOR));
```
`psiE` is the depression angle in degrees. At exact nadir (`psiE == 0`) `sin(0) == 0`,
producing ±infinity which then propagates into `gBuf` and corrupts the output.

**Fix:** Guard the division:
```c
double sinPsi = sin(psiE * DTOR);
gBufTmp[i1][j1] = (sinPsi > 1e-6) ? 10.0 * log10((AbCum / AgCum) / sinPsi) : MINS1DB;
```

**To revert:** Remove the `sinPsi` variable and restore the bare `/ sin(psiE * DTOR)` expression.

---

## geomosaic.c

### 6. Wrong NULL check: `buf1` tested instead of `buf1s` after scale buffer malloc
**Problem:** After allocating the scale buffer:
```c
buf1s = (float *)malloc(...);
if (buf1 != NULL)          // ← tests buf1 (already verified), not buf1s
    fprintf(stderr, "Malloced ... scale buffer ...");
else
    error("Malloc failed for scale buffer ...");
```
A failed `buf1s` allocation is not detected. Execution continues, and the row-pointer
setup for `outputImage->scale` dereferences a NULL `buf1s`, causing a segfault.
The `error()` format string also used `%i` for a `size_t` argument (UB on 64-bit).

**Fix:** Changed `buf1 != NULL` to `buf1s != NULL`; fixed the `error()` format to `%lu`
with a `(size_t)` cast.

**To revert:** Change `buf1s != NULL` back to `buf1 != NULL`; restore `%i` format.

---

### 7. Pointer compared instead of dereferenced: hybridZ/nearestDate guard never fires
**Problem:** `hybridZ` and `nearestDate` are `int32_t *` parameters. Line 713:
```c
if (hybridZ > 0 && nearestDate < 0)
    error("hybrid Z requires a nearest date ");
```
compares the pointer addresses, not the pointed-to values. On a 64-bit system, stack
pointers are always positive, so `hybridZ > 0` is always true and `nearestDate < 0` is
always false. The validation guard never fires, and `hybridZ` without a `nearestDate`
produces silently wrong output.

**Fix:** Added dereference:
```c
if (*hybridZ > 0 && *nearestDate < 0)
```

**To revert:** Remove the `*` dereferences: `hybridZ > 0 && nearestDate < 0`.

---

## processInputFileGeo.c

### 8. `lineCount` uninitialised before first `getDataString` call
**Problem:** `int lineCount, eod;` — `lineCount` is passed uninitialised to
`getDataString` at line 29. `getDataString` uses it as an incrementing line counter
(same convention as all other callers in the codebase, which all initialise `lineCount = 0`).
Error messages would print a garbage line number; on unlucky stack values the
parse could behave incorrectly.

**Fix:** Changed to `int lineCount = 0, eod;`

**To revert:** Remove the `= 0` initialiser.

---

### 9. Non-"poly"/non-"alos" antenna pattern filenames silently discarded
**Problem:**
```c
if (strlen(antPat) > 0 && ((strstr(antPat, "poly") != NULL) || (strstr(antPat, "alos") != NULL)))
    (*antPatFiles)[i] = strdup(antPat);
else
    (*antPatFiles)[i] = NULL;
```
Any antenna pattern filename that doesn't contain the substrings `"poly"` or `"alos"`
is silently nulled out — no warning, no error. The `parseAntPat` function (which this
array feeds) already handles the two special-case tokens internally; the outer guard
here causes real files to be dropped.

**Fix:** Removed the `poly`/`alos` substring check; all non-empty antenna pattern
filenames are now stored:
```c
if (strlen(antPat) > 0)
    (*antPatFiles)[i] = strdup(antPat);
else
    (*antPatFiles)[i] = NULL;
```

**To revert:** Restore the `&& ((strstr(antPat, "poly") != NULL) || (strstr(antPat, "alos") != NULL))`
condition.

---

## readComplexAsPower.c

### 10. Buffer initialisation uses SLC dimensions, overflowing the multilooked buffer
**Problem:**
```c
nPix  = (int64_t)xSize * ySize;   // xSize, ySize from GDAL = full SLC dimensions
power = (float *)inputImage->image[0];
for (i = 0; i < nPix; i++)
    power[i] = (float)-LARGEINT;
```
`xSize` and `ySize` come from `GDALGetRasterBandXSize/YSize`, which return the SLC
(full-resolution) pixel counts. `inputImage->image[0]` is allocated as
`azimuthSize * rangeSize` floats — the multilooked dimensions. When `nAzimuthLooks > 1`,
`xSize * ySize = rangeSize * azimuthSize * nAzimuthLooks`, causing the initialisation
loop to write `nAzimuthLooks` times beyond the end of the heap buffer.

**Fix:** Changed `nPix` to use the allocated (multilooked) dimensions:
```c
nPix  = (int64_t)inputImage->azimuthSize * inputImage->rangeSize;
```

**To revert:** Change back to `nPix = (int64_t)xSize * ySize;`
