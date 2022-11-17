%
%%
% Update AUG 2022:
%   This method now takes an RNASpotsRun for most parameters.
%   There are just too many parameters for listing them all to be
%   practical.

% To generate an RNASpotsRun with default params all you gotta do is:
%   rna_spot_run = RNASpotsRun.initDefault();

%Debug levels: 
%   0 - None
%   1 - Regular verbosity
%   2 - Regular verbosity + output plots
%   3 - High verbosity + output plots

%%
%The preloaded_images param expects a ImageChannelSet
%(./core/ImageChannelSet.m) object. 
%To init one such, you can use the following syntax:
%   preloaded_images = ImageChannelSet;
%Then just set your image channels to its variables!

%%

function rna_spot_run = Adapter_RNASpots(rna_spot_run, guimode, preloaded_images, debug_lvl, limitSaveSize, thread_request)

addpath('./core');
RNA_Fisher_State.setGUIMode(guimode);
bPreloaded = (nargin > 2) & ~isempty(preloaded_images);

if nargin < 4
    debug_lvl = 1;
end

if nargin < 5
    limitSaveSize = true;
end

if nargin < 6
    thread_request = 1;
end

%Check required arguments
if (isempty(rna_spot_run.tif_path)) || (~ischar(rna_spot_run.tif_path))
    if ~bPreloaded
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Input image path is required."), true);
        return;
    else
        %Check preload
        if isempty(preloaded_images)
            RNA_Fisher_State.outputMessageLineStatic(sprintf("Provided preloaded image struct is empty!"), true);
            return;
        end
    end
end
if (~isempty(rna_spot_run.out_dir)) && (~ischar(rna_spot_run.out_dir))
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Output directory argument is invalid."), true);
    return;
end
if (~isempty(rna_spot_run.cellseg_path)) && (~ischar(rna_spot_run.cellseg_path))
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Cellseg data path argument is invalid."), true);
    return;
end
if rna_spot_run.total_ch < 1
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Input image must have at least one channel."), true);
    return;
end
if rna_spot_run.light_ch > rna_spot_run.total_ch
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Light (TRANS) channel index is invalid."), true);
    return;
end
if rna_spot_run.rna_ch > rna_spot_run.total_ch
    RNA_Fisher_State.outputMessageLineStatic(sprintf("RNA/Signal channel index is invalid."), true);
    return;
end

%Set defaults
rna_spot_run.ztrim_auto = 5;
% if ttune_winsize < 1
%     spotsrun.ttune_winsize = 10;
% end
%if ttune_wscorethresh < 0
%    spotsrun.ttune_wscorethresh = 0.9;
%end
if rna_spot_run.t_min < 1
    rna_spot_run.t_min = 1;
end
if rna_spot_run.t_max < 1
    rna_spot_run.t_max = 300;
end
if isempty(rna_spot_run.out_dir)
    %Defaults to input directory
    [rna_spot_run.out_dir, ~, ~] = fileparts(rna_spot_run.tif_path);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Output path not provided. Set to %s", rna_spot_run.out_dir), true);
end
if isempty(rna_spot_run.img_name)
    %Defaults to input file name
    [~, rna_spot_run.img_name, ~] = fileparts(rna_spot_run.tif_path);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Image name not provided. Set to %s", rna_spot_run.img_name), true);
end
if ~isempty(rna_spot_run.ctrl_path)
    if ~ischar(rna_spot_run.ctrl_path)
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Control path is not valid. Setting to empty."), true);
        rna_spot_run.ctrl_path = [];
    end
end
if rna_spot_run.ctrl_ch > rna_spot_run.ctrl_chcount
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Control RNA/signal channel index is invalid. Control will not be used."), true);
    rna_spot_run.ctrl_path = [];
end

%Call core
if bPreloaded
    RNA_Pipeline_Core(rna_spot_run, debug_lvl, preloaded_images, limitSaveSize, thread_request);
else
    RNA_Pipeline_Core(rna_spot_run, debug_lvl, [], limitSaveSize, thread_request);
end

RNA_Fisher_State.clearStaticState();

end