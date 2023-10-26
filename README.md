# TrueSpot
A MATLAB pipeline for automatically processing TIFF image stacks. Functions include cell segmentation, RNA/marker spot detection (with automatic thresholding), and RNA/marker quantification.

## Citations
Cell Segmentation:
Local Maxima Detection:
Spot & Cloud Quantification:

## Dependencies
[MATLAB 2018b (or higher)](https://www.mathworks.com/products/get-matlab.html?s_tid=gn_getml)

[MATLAB Image Processing Toolkit](https://www.mathworks.com/solutions/image-video-processing.html)

[Java JRE 8 (or higher)](https://www.java.com/en/download/manual.jsp)

## Compatibility
Scripts should run on any system with compatible versions of MATLAB and Java.

## Usage
This section covers basic usage commands for each module. Additional usage options are detailed in the repo wiki. We recommend using the wrapper bash scripts for Linux command line usage just to keep things clean, but direct usage for the MATLAB scripts will be outlined here as well.

You do not need to build anything - MATLAB is an interpreter/JIT compiling virtual environment. All you need to use these scripts is MATLAB.

### MATLAB Commmand Line

### Cell Segmentation
MATLAB GUI Interface Script: `Main_CellSegGUI.m`
MATLAB Command Line Interface Script: `Main_CellSegConsole.m`
Bash Wrapper Script (Command Line): `TrueSpot_CellSeg.sh`

### Local Maxima Detection
MATLAB Interface Script: `Main_RNASpots.m`
Bash Wrapper Script: `TrueSpot_RNASpots.sh`

The only required argument is `-input`. See [documentation](https://github.com/neuertlab/RNA-FISH-Auto/blob/main/doc/pages/spots_allargs.md) for full argument list.

**Common Arguments - Input** (Options are the same for both interfaces)
| Argument | Parameter | Description |
| ----- | ----- | ----- |
| `-input` | *Path* - Path to input image file. | (Required) The input image or image stack. |
| `-cellseg` | *Path* - Path to MATLAB file containing cell segmentation mask.  | A MATLAB save file containing a cell segmentation mask output by the cell segementation module. File name will usually start with "Lab_". |
| `-ctrlimg` | *Path* - Path to control image file.  | A control image or image stack file. Only TIFF or MAT files are recognized at this time. |

**Common Arguments - Output** (Options are the same for both interfaces)
| Argument | Parameter | Description |
| ----- | ----- | ----- |
| `-outstem` | *Path Stem* - Output path stem | Filename stem for output files including output directory. (Default: (Directory of input image)/(input image name)_(detection strategy) |
| `-csvout` | *Path* - Output file path | Path to output callset as a csv file. |
| `-runparamout` | *Path* - Output file path | Path dump run parameters to a plain text file. |
| `-csvzero` | - | **Flag** - Output csv table uses 0-based coordinates instead of MATLAB default 1-based coordinates. |

**Common Arguments - Channels**
| Argument | Parameter(s) | Description |
| ----- | ----- | ----- |
| `-chtotal` | *Integer* - Channel Count  | The total number of channels in the input image stack. (Default: 1) |
| `-chsamp` | *Integer* - Channel Index (1-based) | The channel in the input image stack to process as sample. (Default: 1) |
| `-chtrans` | *Integer* - Channel Index (1-based) | The TRANS or passthrough light channel in the input image stack. Only used for background extraction along with cell segmentation mask. |
| `-chctrtotal` | *Integer* - Channel Count  | The total number of channels in the control image stack. (Default: 1) |
| `-chctrsamp` | *Integer* - Channel Index (1-based) | The channel in the control image stack to process as sample control. (Default: 1) |

**Common Arguments - Tuning**
| Argument | Parameter(s) | Description |
| ----- | ----- | ----- |
| `-thmin` | *Integer* - Threshold value | Minimum filtered intensity threshold value to scan, inclusive (Default: 10) |
| `-thmax` | *Integer* - Threshold value | Maximum filtered intensity threshold value to scan, inclusive (Default: 500) |
| `-autominth` | - | **Flag** - Determine threshold scan minimum automatically from image properties. |
| `-automaxth` | - | **Flag** - Determine threshold scan maximum automatically from image properties. |
| `-ztrim` | *Integer* - Z slice count | Number of Z slices to trim out from top and bottom of stack (Default: 0) |
| `-sensitivity` | *Integer* - Level (0-2) | Set thresholding parameters to preset reflecting desired level of sensitivity (Default: 0) |
| `-precision` | *Integer* - Level (0-2) | Set thresholding parameters to preset reflecting desired level of precision (Default: 0) |

**Common Arguments - Misc.**
| Argument | Parameter(s) | Description |
| ----- | ----- | ----- |
| `-imgname` | *String* - Name | Name to assign to image and output files. Useful for images with long unwieldy file names. (Default: Input file name.) |
| `-verbose` | - | **Flag** - Turn on default log verbosity. |
| `-quiet` | - | **Flag** - Minimize log output messages. |

### Spot & Cloud Quantification
MATLAB Interface Script: `Main_RNAQuant.m`
Bash Wrapper Script: `TrueSpot_RNAQuant.sh`

## Contact
