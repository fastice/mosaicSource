# mosaicSource

This repository contains many of the C-language programs used to produce Greenland Ice Mapping Project (GrIMP) velocity mosaics. To compile, code from several other repositories is needed as described below.

## Code Summary

The code is divided into several subdirectories, which are:

**Mosaic3d:** This directory contains much of the main mosaicking code for **mosaic3d** the main program for producing velocity mosaics. Depending on the inputs, it produces velocity estimates for up to 4 configurations:
- Range and azimuth speckle tracked offsets for each individual pair;
- Range from phase and azimuth from offsets for each individual pair;
- Range offsets from a region with crossing ascending and descending orbit coverage using a surface parallel flow assumption; 
- Phase data (interferograms) for a region with crossing ascending and descending orbit coverage using a surface parallel flow assumption;
- Landsat 8 or other optical feature tracked offsets.

The code will build up a mosaic using anywhere from 1 to many 1000s products. At any given point, the final velocity estimate will consist of the 
error-weighted average from whatever combination of the five methods can be applied. See [Joughin 2002](https://www.cambridge.org/core/journals/annals-of-glaciology/article/icesheet-velocity-mapping-a-combined-interferometric-and-speckletracking-approach/816A13E46FCFC2570B435E25EACBA367) for further detail.

**geoMosaic** This directory contains much the mosaicking code for produce calibrated and uncalibrated SAR image mosaics.

**coarseReg** This directory contains a program, **coarsereg**, which is used to coarsely align two SLC images with single range/azimuth offset computed from the metadata.

**getLocC** This directory contains a utility program, **getlocc**, creates a *geodat* file, which has the essential meta data for a multi-look SAR product, which is used by many of the programs listed here.

**tiePoints** This program (**tiepoints**) estimates the baseline parameters for an unwrapped interferogram for a given set of control points.

**rParams** Similar to **tiepoints**, except **rparams** that it calculates the baseline for a range offset rather than phase product.

**azParams** This program (**azParams**) in this directory computes the parameters to calibrate the azimuth offsets.

**simInSAR** This directory contains a program, **siminsar** that can simulate offset and phase products. 

**LLtoRA** This directory contains a program for converting lat/lon/z coordinates to range/azimuths for a SAR Image.

**common** This directory contains several functions that are common to more than one of the above programs.

**landsatMosaic** This directory contains the Landsat mosaicking modules for **mosaic3d**

# Other Repositories Needed to Compile

This mosaicking code requires code from the following repositories:

[cRecipes](https://github.com/fastice/cRecipes): A handful of Numerical Recipes from C routines with modifications as needed.

[clib](https://github.com/fastice/clib/tree/master): Some utilities for opening files used by the mosaicking code.

[geotiff](https://github.com/fastice/geotiff): Some geotiff library include files for cases where they are not part of the system installation.

[landsatSource64](https://github.com/fastice/landsatSource64/tree/master): Code for offset tracking in pairs of Landsat image with the same path row. Some of the code is required by the mosaicker.

[triangle](https://github.com/fastice/triangle/blob/master/README.md): Ancient triangle code from the 1990s need to compile (maintains some rarely used functionality).


