# RNA-FISH-Auto

A MATLAB pipeline for automatically processing TIFF image stacks. Functions include cell segmentation, RNA/marker spot detection (with automatic thresholding), and RNA/marker quantification.

## Dependencies
[MATLAB 2018b](https://www.mathworks.com/products/get-matlab.html?s_tid=gn_getml) (Newer versions may work, but older versions probably will not)
[MATLAB Image Processing Toolkit](https://www.mathworks.com/solutions/image-video-processing.html)
[Java JRE 8 (or higher)](https://www.java.com/en/download/manual.jsp)

# Compatibility
Scripts should run on any system with compatible versions of MATLAB and Java.

**Caveat**
The dead pixel detection routine calls the native function `imhistc`, which is normally wrapped by `imhist` and hidden from the MATLAB user interface. To get around this, we copied `imhistc.m` and `imhistc.mexw64` from a Windows 8 installation to the workspace folder. These files have been removed from the public repo since they are part of MATLAB itself. At this time, the dead pixel detection function should attempt to search for `imhistc` and if it doesn't find it, call `imhist` instead. I do not know how well this works yet.
If you want to use `imhistc` directly as we do, you can try to copy the files from your installation to the workspace directory. On Windows, they can be found at:
`%MATLAB_INSTALLATION_DIR%\`
``

## Usage
