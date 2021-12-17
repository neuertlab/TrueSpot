%Main script for automatic background mask generation
%Blythe Hospelhorn
%Version 1.0.0
%Updated Feb 19, 2021

%Update history
%   2021.02.19 | 1.0.0
%       - Initial documentation

%%
%   Automatically detect an intensity threshold for background from the
%   light channel of an image stack and cell segmentation mask and output a
%   new 2D boolean mask indicating what parts of the image are and are not
%   background.
%   The background can then be isolated to use as a control for spot
%   detection.
%
% ARGS
%   img_path (string) - Path on local file system to image stack to extract
%       background mask from. Must be tiff file.
%   seg_path (string) - Path on local file system to cell segmentation data
%       from previous pipeline (usually "Lab_(image label).mat")
%   out_stem (string) - Prefix for output files. Should include directory
%       path as well as desired file name. Output file will be out_stem.mat
%   totalCh (int) - Total number of channels in input image stack.
%   lightCh (int) - One-based index of white light channel in image stack.
%   outputpng (bool) - True to output png images of the masked image and
%       histogram in addition to other outputs (use for QC).
%
% OUTPUTS
%   - Background mask saved in a .mat file specified by out_stem argument.
%
% NOTES
%
function Main_BackgroundMask(img_path, seg_path, out_stem, totalCh, lightCh, outputpng)

    addpath('./core');

    %Force to numbers
    totalCh = Force2Num(totalCh);
    lightCh = Force2Num(lightCh);
    outputpng = Force2Bool(outputpng);

    [tif, ~] = LoadTif(img_path, totalCh, [lightCh], true);
    Bkg_Mask_Core(tif{lightCh,1}, seg_path, out_stem, outputpng);

end