# Cell Segmentation Module - Full Argument List
This page lists all arguments recognized by `Main_CSCellpose.m`. Some additional arguments may be found directly in the script as well, but those are either experimental or have been deprecated and are nonfunctional.

CellDissect (TrueSpot's default cell segmentation algorithm) can optionally be used for nuclear segmentation with Cellpose handling cell/cytoplasmic segmentation. If CellDissect nuclear segmentation is used, the nuclear parameters are the same as for the standard CellSeg interface.

[Back](../dochome.md)

## Required

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-input` | *Path* - Path to input image file.  | The input image or image stack containing light/TRANS channel. |

## Input Options

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-innuc` | *Path* - Path to image file.  | The input image stack containing the nuclear marker channel, if different from the main input stack. |
| `-chtotal` | *Integer* - Channel Count  | The total number of channels in the input image stack. (Default: 1) |
| `-chtotnuc` | *Integer* - Channel Count  | The total number of channels in the `innuc` image stack, if provided. (Default: 1) |
| `-chlight` | *Integer* - Channel Index (1-based) | The TRANS or passthrough light channel in the input image stack. Only used for background extraction along with cell segmentation mask. |
| `-chnuc` | *Integer* - Channel Index (1-based) | The nuclear stain (eg. DAPI) channel in the input image stack, or the `innuc` image stack if one is provided. (Default: 1) |
| `-nuczmin` | *Integer* - Slice Index (1-based) | The minimum slice in the nuclear marker channel to include in nucleus detection. (Default: 1) |
| `-nuczmax` | *Integer* - Slice Index (1-based) | The maximum slice in the nuclear marker channel to include in nucleus detection. (Default: Z) |
| `-lightzmin` | *Integer* - Slice Index (1-based) | The minimum slice in the TRANS channel to include in cell boundary detection. (Default: 1) |
| `-lightzmax` | *Integer* - Slice Index (1-based) | The maximum slice in the TRANS channel to include in cell boundary detection. (Default: Z) |
| `-skipnuc` | - | **Flag** - Skip nuclear segmentation |

## Output Options
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-outpath` | *Path* - Directory path | Path to directory to place output files. File names will be automatically generated from image name. (Default: Directory of input image.) |
| `-ocellmask` | *Path* - Path to image (tif) file output | Path to write rendered TIF of cell mask. (Default: None) |
| `-onucmask` | *Path* - Path to image (tif) file output | Path to write rendered TIF of nuclear mask. (Default: None) |
| `-osettings` | *Path* - Path to text file output | Path to write text file with run settings information. (Default: None.) |
| `-ovrw` | - | **Flag** - Overwrite any existing output. |
| `-dumpsummary` | - | **Flag** - Dump text file containing input parameter summary. |

## Metadata
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-imgname` | *String* - Name | Name to assign to image and output files. Useful for images with long unwieldy file names. (Default: Input file name.) |

## Segmentation Tuning

### Presets
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-template` | *String* - Template name | Name of preset template to use. |
| `-savetmpl` | *String* - Template name | Save settings for this run as a preset template with the specified name. |

### Fine Tune - Common

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-voxelsize` | *Int pair or triplet* - Voxel size formatted "(x,y,z)" | Voxel dimensions in nanometers. This is used to calculate the z to xy ratio for 3D nuclear segmentation. |
| `-ensemble` | - | **Flag** - Use ensemble models for both cell/cytoplasmic and nuclear segmentation |
| `-norm` | - | **Flag** - Normalize image for both cell/cytoplasmic and nuclear segmentation |

### Fine Tune - Nucleus

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-navgdia` | *Integer* - Size in pixels | Average expected nuclear diameter, in pixels (Default: NaN) |
| `-nmodel` | *String* - Model name | Name of Cellpose model to use (Default: 'nuclei') |
| `-ncth` | *Number* - Threshold value | Nuclear [cell threshold value](https://www.mathworks.com/help/medical-imaging/ref/cellpose.segmentcells2d.html) (see Cellpose documentation) (Default: 0) |
| `-nfth` | *Number* - Threshold value | Nuclear [flow threshold value](https://www.mathworks.com/help/medical-imaging/ref/cellpose.segmentcells2d.html) (see Cellpose documentation) (Default: 0.4) |
| `-nszmin` | *Integer* - Size in pixels | Minimum expected area of nucleus, in pixels (Default: 15) |
| `-nszmax` | *Integer* - Size in pixels | Maximum expected area of nucleus, in pixels (Default: NaN) |
| `-nensemble` | - | **Flag** - Use ensemble model for nuclear segmentation |
| `-nnorm` | - | **Flag** - Normalize image for nuclear segmentation |

### Fine Tune - Nucleus (CellDissect)

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-cdnuc` | - | **Flag** - Use CellDissect (3D) for nuclear segmentation instead of Cellpose |
| `-xtrim` | *Integer* - Pixel count | Number of pixels to trim off each edge in the x direction (Default: 4) |
| `-ytrim` | *Integer* - Pixel count | Number of pixels to trim off each edge in the y direction (Default: 4) |
| `-nzrange` | *Integer* - Number of Z slices | Number of z slices around plane with highest standard deviation to use for nucleus segmentation (Default: 3) |
| `-nthsmpl` | *Integer* - Test count | Number of DAPI/nuclear stain thresholds to sample (Default: 200) |
| `-nszmin` | *Integer* - Size in pixels | Minimum expected area of nucleus, in pixels (Default: 40) |
| `-nszmax` | *Integer* - Size in pixels | Maximum expected area of nucleus, in pixels (Default: 200) |
| `-ncutoff` | *Float* - Proportion | Proportion of tested thresholds a pixel must be present in nuclear mask of to not be filtered out. (Default: 0.05) |
| `-ndxy` | *Float* - Pixel count | Radius of pixels to look around a putative nucleus center. (Default: round(sqrt(nszmax/pi))) |

### Fine Tune - Cell/Cytoplasm

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-cavgdia` | *Integer* - Size in pixels | Average expected cell diameter, in pixels (Default: NaN) |
| `-cmodel` | *String* - Model name | Name of Cellpose model to use (Default: 'cyto') |
| `-ccth` | *Number* - Threshold value | Cytoplasmic [cell threshold value](https://www.mathworks.com/help/medical-imaging/ref/cellpose.segmentcells2d.html) (see Cellpose documentation) (Default: 0) |
| `-cfth` | *Number* - Threshold value | Cytoplasmic [flow threshold value](https://www.mathworks.com/help/medical-imaging/ref/cellpose.segmentcells2d.html) (see Cellpose documentation) (Default: 0.4) |
| `-cszmin` | *Integer* - Size in pixels | Minimum expected area of cell, in pixels (Default: 15) |
| `-cszmax` | *Integer* - Size in pixels | Maximum expected area of cell, in pixels (Default: NaN) |
| `-censemble` | - | **Flag** - Use ensemble model for cell/cytoplasmic segmentation |
| `-cnorm` | - | **Flag** - Normalize image for cell/cytoplasmic segmentation |
