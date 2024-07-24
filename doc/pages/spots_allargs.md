# Spot Detection Module - Full Argument List
This page lists all arguments recognized by `Main_RNASpots.m`. Some additional arguments may be found directly in the script as well, but those have been deprecated and are nonfunctional.

[Back](../dochome.md)

## Required

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-input` | *Path* - Path to input image file.  | The input image or image stack. This parameter may be substituted with `-tif` for TIFF images, or `-matimg` and `-matvar` for MATLAB save files. |

## Input Options

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-tif` | *Path* - Path to TIFF input image file.  | The input image or image stack TIFF file. This argument is equivalent to `-input`, but is included for convenience and compatibility with previous versions. |
| `-matimg` | *Path* - Path to MATLAB input image file.  | A MATLAB save file containing the image/image stack to process. This argument is equivalent to `-input`, but is included for convenience and compatibility with previous versions. |
| `-matvar` | *String* - Name of variable containing image in MATLAB input file.  | The name of the MATLAB variable in the input MATLAB save file containing the image/image stack data to process. Defaults to "imgdat". |
| `-cellseg` | *Path* - Path to file containing cell segmentation mask.  | A file containing a cell segmentation mask. Accepts CellSeg module MAT outputs as well as single channel TIFFs, and numeric tsvs and csvs. |
| `-ctrlimg` | *Path* - Path to control image file.  | A control image or image stack file. Only TIFF or MAT files are recognized at this time. |
| `-ctrltif` | *Path* - Path to control TIFF file.  | A control image or image stack file in TIFF format. This argument is equivalent to `-ctrlimg`, but is included for convenience and compatibility with previous versions. |
| `-chtotal` | *Integer* - Channel Count  | The total number of channels in the input image stack. (Default: 1) |
| `-chsamp` | *Integer* - Channel Index (1-based) | The channel in the input image stack to process as sample. (Default: 1) |
| `-chtrans` | *Integer* - Channel Index (1-based) | The TRANS or passthrough light channel in the input image stack. Only used for background extraction along with cell segmentation mask. |
| `-chctrtotal` | *Integer* - Channel Count  | The total number of channels in the control image stack. (Default: 1) |
| `-chctrsamp` | *Integer* - Channel Index (1-based) | The channel in the control image stack to process as sample control. (Default: 1) |

## Output Options
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-outdir` | *Path* - Directory path | Path to directory to place output files. File names will be automatically generated from image name. (Default: Directory of input image.) |
| `-outstem` | *Path Stem* - Output path stem | Filename stem for output files including output directory. (Default: (Directory of input image)/(input image name)_(detection strategy) |
| `-csvout` | *Path* - Output file path | Path to output callset as a csv file. |
| `-runparamout` | *Path* - Output file path | Path dump run parameters to a plain text file. |
| `-ovrw` | - | **Flag** - Overwrite any existing output. |
| `-csvzero` | - | **Flag** - Output csv table uses 0-based coordinates instead of MATLAB default 1-based coordinates. |
| `-csvfull` | - | **Flag** - Output csv table contains all calls at all thresholds. |
| `-csvrange` | - | **Flag** - Output csv table contains only calls at or above the bottom of the auto-threshold range. |
| `-csvthonly` | - | **Flag** - Output csv table contains only calls at or above the automatically called threshold. |

## Metadata
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-imgname` | *String* - Name | Name to assign to image and output files. Useful for images with long unwieldy file names. (Default: Input file name.) |
| `-voxelsize` | *Int pair or triplet* - Voxel size formatted "(x,y,z)" | Voxel dimensions in nanometers. This is only recorded as metadata and is not used or required for tool function. Functions the same as `-pixelsize`. |
| `-pixelsize` | *Int pair or triplet* - Pixel size formatted (x,y) | Pixel dimensions in nanometers. This is only recorded as metadata and is not used or required for tool function. Functions the same as `-voxelsize`. |
| `-expspotsize` | *Int pair or triplet* - 2D or 3D size formatted (x,y) or (x,y,z) | Expected size of box dimensions enveloping a single spot in nanometers. This is only recorded as metadata and is not used or required for tool function. |
| `-probetype` | *String* - Probe name | Name of fluorescent or other probe used to illuminate sample in target channel. |
| `-target` | *String* - Target name | Name of target biomolecule. |
| `-targettype` | *String* - Target type descriptor | The type of biomolecule target is (eg. "lncRNA" or "Protein") |
| `-species` | *String* - Species name | The species the sample was sourced from. |
| `-celltype` | *String* - Cell type name/descriptor | The cell type of the sample. |

## Detection Tuning
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-thmin` | *Integer* - Threshold value | Minimum filtered intensity threshold value to scan, inclusive (Default: 10) |
| `-thmax` | *Integer* - Threshold value | Maximum filtered intensity threshold value to scan, inclusive (Default: 500) |
| `-ztrim` | *Integer* - Z slice count | Number of Z slices to trim out from top and bottom of stack (Default: 0) |
| `-zmin` | *Integer* - Z slice index (1-based) | Bottom Z slice to include in processing. |
| `-zmax` | *Integer* - Z slice index (1-based) | Top Z slice to include in processing. |
| `-gaussrad` | *Integer* - Radius in pixels | Radius of Gaussian filter (Default: 7) |
| `-nodpc` | - | **Flag** - Skip dead pixel cleaning step. |
| `-maxzproj` | - | **Flag** - Use maximum Z projection for maxima detection instead of the full 3D stack. |
| `-autominth` | - | **Flag** - Determine threshold scan minimum automatically from image properties. |
| `-automaxth` | - | **Flag** - Determine threshold scan maximum automatically from image properties. |

## Thresholding Tuning

### Presets
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-sensitivity` | *Integer* - Level (0-5) | Set thresholding parameters to preset reflecting desired level of sensitivity (Default: 0) |
| `-precision` | *Integer* - Level (0-5) | Set thresholding parameters to preset reflecting desired level of precision (Default: 0) |
| `-sensitive` | - | **Flag** - Set thresholding parameters to a preset with higher sensitivity. |
| `-precise` | - | **Flag** - Set thresholding parameters to a preset with higher precision. |

### Fine Tune
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-mfmin` | *Float* - Multiplier | Minimum MAD factor to use for threshold range determination (Default: -1.0) |
| `-mfmax` | *Float* - Multiplier | Maximum MAD factor to use for threshold range determination (Default: 1.0) |
| `-wszmin` | *Integer* - Window size | Minimum Fano factor window size to scan for threshold range determination (Default: 3) |
| `-wszmax` | *Integer* - Window size | Maximum Fano factor window size to scan for threshold range determination (Default: 21) |
| `-wszincr` | *Integer* - Window size | Increment between window sizes to try (Default: 3) |
| `-thmw` | *Float* - Weight | Weight factor for inclusion of median/MAD derived threshold candidates in final call (Default: 0.25) |
| `-thfw` | *Float* - Weight | Weight factor for inclusion of two-piece fit knot threshold candidates in final call (Default: 0.0) |
| `-thiw` | *Float* - Weight | Weight factor for inclusion of two-piece fit right intercept threshold candidates in final call (Default: 0.75) |
| `-fitlog` | *Boolean* - True/False | Whether to use log10 transformations of spot count derived curves for threshold determination (Default: true) |
| `-stdfac` | *Float* - Standard deviations | Brute-force shift threshold selection by this many standard deviations of the threshold candidate pool (Default: 0.0) |
| `-usespc` | - | **Flag** - Include spot count curve directly in threshold seletion |
| `-usedfc` | - | **Flag** - Include absolute value of first deriv. approximate transformation directly in threshold seletion. |
| `-fitwavg` | - | **Flag** - Use a different weighting algorithm to determine two-piece fits |

## Verbosity
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-verbose` | - | **Flag** - Turn on default log verbosity. |
| `-quiet` | - | **Flag** - Minimize log output messages. |
| `-debug` | - | **Flag** - Run with basic debug logging/output. |
| `-debugv` | - | **Flag** - Run with extra verbose debug logging/output. |

## Performance
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-threads` | *Integer* - Thread number | Number of threads to request from MATLAB. As many threads may run in parallel as there are available cores. Parallelizing detection speeds it up considerably. (Default: 1) |
