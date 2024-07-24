# Quantification Module - Full Argument List
This page lists all arguments recognized by `Main_RNAQuant.m`. Some additional arguments may be found directly in the script as well, but those have been deprecated or are for debug purposes and it is recommended not to use them to avoid confusion.

No single options are required, but some form of input (ie. image path or spots run path) must be provided.

[Back](../dochome.md)

## Input/Output

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-runinfo` | *Path* - Path to spots run info file  | MAT file (generally ending in "_rnaspotsrun.mat") produced by the Spot Detection module. Contains links and run parameters. |
| `-tif` | *Path* - Path to TIF image file  | File containing image/stack to be analyzed. (Default: Path stored in provided run info)|
| `-outdir` | *Path* - Output directory path  | Directory to place results files in (Default: Directory of primary input file)|
| `-cellsegpath` | *Path* - Path to cell mask  | File containing a 2D or 3D cell segmentation mask. Accepts CellSeg module MAT outputs, single channel TIF, and unheaded csv or tsv. |
| `-nucsegpath` | *Path* - Path to nucleus mask | File containing a 2D or 3D nuclear segmentation mask, if different from cell mask. Accepts CellSeg module MAT outputs, single channel TIF, and headerless csv or tsv. (Default: cellsegpath) |
| `-coordtable` | *Path* - Path to maxima coordinate table  | File containing xyz coordinates of maxima calls made by a spot detector. Accepts call tables output by spot detect module or headerless text tables (csv or tsv) with columns ordered x,y,z. (Default: Path referenced by spotsrun, if provided) |
| `-ch` | *Int* - Channel Index (1-based)  | The channel in the input image stack to process as sample. (Default: 1, or value stored in spotsrun) |
| `-chcount` | *Int* - Channel Count | The total number of channels in the input image stack. (Default: 1, or value stored in spotsrun) |

## Threshold

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-mthresh` | *Int* - Threshold value  | Threshold value to use to obtain starting calls (Default: Value stored in spotsrun) |
| `-rethresh` | - | **Flag** - Rerun automatic thresholding using default parameters or parameters stored in spotsrun. |

## Thresholding Tuning

### Presets
| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-sensitivity` | *Integer* - Level (0-5) | Set thresholding parameters to preset reflecting desired level of sensitivity (Default: 0) |
| `-precision` | *Integer* - Level (0-5) | Set thresholding parameters to preset reflecting desired level of precision (Default: 0) |
| `-sensitive` | - | **Flag** - Set thresholding parameters to a preset with higher sensitivity. |
| `-precise` | - | **Flag** - Set thresholding parameters to a preset with higher precision. |

## Tuning

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-zadj` | *Float* - Adjustment factor  | Ratio of distance between z-planes to xy pixel size. (Default: 1.0, or z/x ratio from runinfo voxel dimensions if stored) |
| `-smobjsz` | *Int* - Size in pixels  | Maximum size of "small objects" during noise cleanup. (Default: 3) |
| `-gaussrad` | *Int* - Size in pixels | Gaussian filter radius. (Default: 7, or value stored in spotsrun) |
| `-radxy` | *Int* - Size in pixels  | Pixel radius in xy to look around a maximum call during gaussian fitting. (Default: 4) |
| `-radz` | *Int* - Size in pixels  | Radius in z to look around a maximum call during gaussian fitting. (Default: 2) |
| `-norefilter` | - | **Flag** - Skip refiltering to rederive maxima list. If this flag is set, a coordinate table path is expected either directly or via the spotsrun. |
| `-nmhi` | - | **Flag** - Use the nuclear mask derived from the high threshold (only applicable to output from CellDissect/TrueSpot's cellseg module). |
| `-nmlo` | - | **Flag** - Use the nuclear mask derived from the low threshold (only applicable to output from CellDissect/TrueSpot's cellseg module). |
| `-nm2d` | - | **Flag** - Use the 2D nuclear mask (only applicable to output from CellDissect/TrueSpot's cellseg module). |
| `-nocells` | - | **Flag** - Evaluate input image as a whole instead of cell-by-cell. If this flag is set, cell segmentation paths are ignored. |
| `-noclouds` | - | **Flag** - Skip cloud detection and only do gaussian fitting. |

## Performance

| Name | Parameter | Description |
| ----- | ----- | ----- |
| `-workers` | *Int* - Thread number  | Number of threads to request from MATLAB. As many threads may run in parallel as there are available cores. When there is more than one worker thread, fitting/quantification can be done on multiple cells at once. (Default: 1) |



