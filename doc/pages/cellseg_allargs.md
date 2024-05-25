# Cell Segmentation Module - Full Argument List
This page lists all arguments recognized by `Main_CellSegConsole.m`. Some additional arguments may be found directly in the script as well, but those have been deprecated and are nonfunctional.

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

NOTE: The GUI version appears to have 13 hardcoded as the equivalent to `nuczmin` (see `A1_segment_predefined_variables_streamlined_generalized.m`).

## Output Options
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-outpath` | *Path* - Directory path | Path to directory to place output files. File names will be automatically generated from image name. (Default: Directory of input image.) |
| `-ocellmask` | *Path* - Path to image (tif) file output | Path to write rendered TIF of cell mask. (Default: None) |
| `-onucmask` | *Path* - Path to image (tif) file output | Path to write rendered TIF of nuclear mask. (Default: None) |
| `-osettings` | *Path* - Path to text file output | Path to write text file with run settings information. (Default: None.) |
| `-ovrw` | - | **Flag** - Overwrite any existing output. |
| `-dumpsummary` | - | **Flag** - Dump text file containing input parameter summary. |
| `-savebk` | - | **Flag** - Use previous version file naming conventions for output. |

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

### Fine Tune
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-cszmin` | *Integer* - Size in pixels | Minimum expected area of cell, in pixels (Default: 600) |
| `-cszmax` | *Integer* - Size in pixels | Maximum expected area of cell, in pixels (Default: 1200) |
| `-fplstrat` | *String* - Strategy name | Strategy for selecting focus plane for cell segmentation. (Default: 'specify') |
| `-fzmin` | *Integer* - 1-based Z slice index| Minimum z plane of focus region. Ignored if fplstrat is not 'specify'. (Default: 1) |
| `-fzmax` | *Integer* - 1-based Z slice index| Maximum z plane of focus region. Ignored if fplstrat is not 'specify'. (Default: Z) |
| `-foffmin` | *Integer* - Number of Z slices | Offset from focus plane to start bottom of focus region. Ignored if fplstrat is not 'midplane' or 'midplane2'. (Default: 3) |
| `-foffmax` | *Integer* - Number of Z slices | Offset from focus plane to start top of focus region. Ignored if fplstrat is not 'midplane' or 'midplane2'. (Default: 7) |
| `-xtrim` | *Integer* - Pixel count | Number of pixels to trim off each edge in the x direction (Default: 4) |
| `-ytrim` | *Integer* - Pixel count | Number of pixels to trim off each edge in the y direction (Default: 4) |
| `-nzrange` | *Integer* - Number of Z slices | Number of z slices around plane with highest standard deviation to use for nucleus segmentation (Default: 3) |
| `-nthsmpl` | *Integer* - Test count | Number of DAPI/nuclear stain thresholds to sample (Default: 200) |
| `-nszmin` | *Integer* - Size in pixels | Minimum expected area of nucleus, in pixels (Default: 40) |
| `-nszmax` | *Integer* - Size in pixels | Maximum expected area of nucleus, in pixels (Default: 200) |
| `-ncutoff` | *Float* - Proportion | Proportion of tested thresholds a pixel must be present in nuclear mask of to not be filtered out. (Default: 0.05) |
| `-ndxy` | *Float* - Pixel count | Radius of pixels to look around a putative nucleus center. (Default: round(sqrt(nszmax/pi))) |

### fplstrat Valid Values
* `specify`
* `max_cells`
* `midplane`
* `midplane2`
* `first5`
* `last5`

## Built-In Presets
* `jurkat_20x`
* `jurkat_100x`
* `mESC_20x`
* `mESC_100x`
* `sacCer_100x`
* `sPombe_20x`
* `sPombe_100x`



