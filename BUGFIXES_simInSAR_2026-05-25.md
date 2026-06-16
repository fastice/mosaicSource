# Bug fixes applied to simInSAR 2026-05-25

Applied findings from automated code review of mosaicSource/simInSAR/.
All changes compile clean with no warnings.

---

## siminsar.c

### 1. bnStartFlag / bpStartFlag swapped — wrong baseline component assigned
**Problem:** The two flag guards at the end of `readArgs` were transposed:
```c
if (bpStartFlag == TRUE) {      // WRONG: guarded by bp flag
    scene->bnStart = bnStart;   // but assigning bn values
    scene->bnEnd   = bnEnd;
}
if (bnStartFlag == TRUE) {      // WRONG: guarded by bn flag
    scene->bpStart = bpStart;   // but assigning bp values
    scene->bpEnd   = bpEnd;
}
```
When the user supplies `-bnStart`/`-bnEnd`, the values landed in `bpStart`/`bpEnd` and
vice versa, corrupting the baseline ramp for both components and producing wrong
interferogram phase throughout the scene.

**Fix:** Swapped the flag names so each flag guards its own component:
```c
if (bnStartFlag == TRUE) { scene->bnStart = bnStart; scene->bnEnd = bnEnd; }
if (bpStartFlag == TRUE) { scene->bpStart = bpStart; scene->bpEnd = bpEnd; }
```

**To revert:** Swap `bnStartFlag`/`bpStartFlag` back.

---

### 2. parseBnBpParamsFile return value ignored — silent failure on bad file
**Problem:** `parseBnBpParamsFile(bpParamsFile, &bn, &bp, &dBn, &dBp)` returns -1 on
failure but the return value was discarded. The simulation would continue silently with
uninitialised baseline parameters if the file was missing or malformed.

**Fix:**
```c
if (parseBnBpParamsFile(bpParamsFile, &bn, &bp, &dBn, &dBp) != 0)
    error("readArgs: failed to parse bParamsFile %s\n", bpParamsFile);
```

**To revert:** Remove the `if (...)` wrapper and restore the bare call.

---

## getSlantRangeDEM.c

### 3. freadBS nitems/size inverted — byteswap disabled, all floats read wrong
**Problem:**
```c
freadBS(dem->demData[i], scene.I.rangeSize * sizeof(float), 1, fp, FLOAT32FLAG);
```
`freadBS` signature is `freadBS(ptr, nitems, element_size, fp, flags)`. Passing
`rangeSize * sizeof(float)` as `nitems` and `1` as `element_size` tells the byteswap
loop that each element is 1 byte, so no swapping occurs. The file is read correctly
byte-for-byte, but each 4-byte float is left in its original big-endian byte order,
producing garbage float values on x86_64.

**Fix:** `freadBS(dem->demData[i], scene.I.rangeSize, sizeof(float), fp, FLOAT32FLAG);`

**To revert:** Change back to `scene.I.rangeSize * sizeof(float), 1`.

---

### 4. Rows allocated independently, violating contiguous-block convention
**Problem:** The codebase requires `float **` 2D arrays to be a single contiguous
allocation (via `mallocImage`) with a separate row-pointer array. The original code
allocated rows with individual `malloc` calls:
```c
dem->demData = (float **)malloc(scene.I.azimuthSize * sizeof(float *));
for (i = 0; i < scene.I.azimuthSize; i++)
    dem->demData[i] = (float *)malloc(scene.I.rangeSize * sizeof(float));
```
Any code that treats the array as a flat buffer (e.g., `memcpy(dem->demData[0], ...,
azimuthSize * rangeSize * sizeof(float))`) would access only row 0's allocation.

**Fix:** Replaced with `dem->demData = mallocImage(scene.I.azimuthSize, scene.I.rangeSize);`

**To revert:** Restore the `malloc(azimuthSize * sizeof(float *))` + inner `malloc` loop.

---

## parseSceneFile.c

### 5. getVRTOffsetMeta() — missing return statement and no NULL check after GDALOpen
**Problem (a):** Function declared as `static GDALRasterBandH getVRTOffsetMeta(...)` but
has no `return` statement. The caller does not use the return value, so no runtime crash,
but the missing `return` is undefined behaviour in C.

**Fix (a):** Changed declaration and definition to `static void getVRTOffsetMeta(...)`.

**Problem (b):** `GDALOpen` at line 18 can return `NULL` if the file is not found or
unreadable. The return value was passed directly to `GDALGetRasterXSize(hDS)` without a
NULL check, causing a crash on a bad path.

**Fix (b):** Added:
```c
if (hDS == NULL)
    error("getVRTOffsetMeta: could not open %s\n", datFile);
```

**To revert:** Restore `static GDALRasterBandH` return type; remove the NULL check.

---

### 6. lineCount uninitialised before first getDataString call
**Problem:** `int32_t lineCount, eod, nRead;` in `parseOffsetParamFile` — `lineCount`
is passed uninitialised to `getDataString`. Error messages would print a garbage line
number.

**Fix:** Changed to `int32_t lineCount = 0, eod, nRead;`

**To revert:** Remove the `= 0` initialiser.

---

### 7. Division by zero when azimuthSize == 1
**Problem:**
```c
scene->bnStep = (scene->bnEnd - scene->bnStart) / (double)(scene->I.azimuthSize - 1.0);
scene->bpStep = (scene->bpEnd - scene->bpStart) / (double)(scene->I.azimuthSize - 1.0);
```
When `azimuthSize == 1`, this divides by zero, producing ±Inf which propagates into
`localBn`/`localBp` and corrupts every simulated pixel.

**Fix:**
```c
double azDenom = (scene->I.azimuthSize > 1) ? (double)(scene->I.azimuthSize - 1) : 1.0;
scene->bnStep = (scene->bnEnd - scene->bnStart) / azDenom;
scene->bpStep = (scene->bpEnd - scene->bpStart) / azDenom;
```

**To revert:** Remove `azDenom` and restore the direct `/ (double)(azimuthSize - 1.0)` division.

---

## outputSimulatedImage.c

### 8. Log message always prints .lat filename for both .lat and .lon iterations
**Problem:** Inside the `for (k=0; k<2; k++)` loop, the log message was:
```c
fprintf(stderr, "writing %s\n", buf1);
```
`buf1` holds only the `.lat` filename. For `k==1` (`.lon`), `file` points to `buf2`
but the message still printed `buf1`, making the log misleading.

**Fix:** Changed to `fprintf(stderr, "writing %s\n", file);` — `file` is correctly
set to `buf1` or `buf2` on each iteration by `appendSuffix`.

**To revert:** Change `file` back to `buf1`.

---

## simInSARimage.c

### 9. rgSave assigned but never read — dead code
**Problem:** `double rgSave` was declared, listed in `#pragma omp parallel private(...)`,
and assigned `rgSave = rg` when `jLoop == 0`, but never read anywhere. Remnant of a
removed feature.

**Fix:** Removed `rgSave` from the declaration, from the `private(...)` clause, and
deleted the `rgSave = rg` assignment.

**To revert:** Re-add `double iFloat, rgSave;`, add `rgSave` to the `private` clause,
and restore `if (jLoop == 0) rgSave = rg;`.
