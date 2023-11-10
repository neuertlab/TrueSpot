%
%%
%The preloaded_images param expects a ImageChannelSet
%(./core/ImageChannelSet.m) object. 
%To init one such, you can use the following syntax:
%   preloaded_images = ImageChannelSet;
%Then just set your image channels to its variables!

%%

function rna_spot_run = Adapter_RNASpotsv2(img_name, tif_path, rna_ch, light_ch, total_ch,...
    out_dir, t_min, t_max, ztrim, cellseg_path, ctrl_path, ctrl_ch, ctrl_chcount, ttune_winsize, ttune_wscorethresh,...
    overwrite_output, guimode, preloaded_images)

addpath('./core');
RNA_Fisher_State.setGUIMode(guimode);
bPreloaded = (nargin > 17);

%Build a Run object for easy save/load
spotsrun = RNASpotsRun;
rna_spot_run = spotsrun;

%Save string args to run obj
spotsrun.img_name = img_name;
spotsrun.tif_path = tif_path;
spotsrun.out_dir = out_dir;
spotsrun.cellseg_path = cellseg_path;
spotsrun.ctrl_path = ctrl_path;

%Check arguments
%Force to numbers
spotsrun.rna_ch = Force2Num(rna_ch);
spotsrun.light_ch = Force2Num(light_ch);
spotsrun.total_ch = Force2Num(total_ch);
spotsrun.ctrl_ch = Force2Num(ctrl_ch);
spotsrun.ctrl_chcount = Force2Num(ctrl_chcount);
spotsrun.t_min = Force2Num(t_min);
spotsrun.t_max = Force2Num(t_max);
spotsrun.ztrim = Force2Num(ztrim);
spotsrun.ttune_winsize = Force2Num(ttune_winsize);
spotsrun.ttune_wscorethresh = Force2Num(ttune_wscorethresh);
spotsrun.overwrite_output = Force2Bool(overwrite_output);

%Check required arguments
if (isempty(tif_path)) || (~ischar(tif_path))
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
if (~isempty(out_dir)) && (~ischar(out_dir))
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Output directory argument is invalid."), true);
    return;
end
if (~isempty(cellseg_path)) && (~ischar(cellseg_path))
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Cellseg data path argument is invalid."), true);
    return;
end
if total_ch < 1
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Input image must have at least one channel."), true);
    return;
end
if light_ch > total_ch
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Light (TRANS) channel index is invalid."), true);
    return;
end
if rna_ch > total_ch
    RNA_Fisher_State.outputMessageLineStatic(sprintf("RNA/Signal channel index is invalid."), true);
    return;
end

%Set defaults
spotsrun.ztrim_auto = 5;
if ttune_winsize < 1
    spotsrun.ttune_winsize = 10;
end
if ttune_wscorethresh < 0
    spotsrun.ttune_wscorethresh = 0.9;
end
if t_min < 1
    spotsrun.t_min = 1;
end
if t_max < 1
    spotsrun.t_max = 300;
end
if isempty(spotsrun.out_dir)
    %Defaults to input directory
    [spotsrun.out_dir, ~, ~] = fileparts(spotsrun.tif_path);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Output path not provided. Set to %s", spotsrun.out_dir), true);
end
if isempty(spotsrun.img_name)
    %Defaults to input file name
    [~, spotsrun.img_name, ~] = fileparts(spotsrun.tif_path);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Image name not provided. Set to %s", spotsrun.img_name), true);
end
if ~isempty(spotsrun.ctrl_path)
    if ~ischar(spotsrun.ctrl_path)
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Control path is not valid. Setting to empty."), true);
        spotsrun.ctrl_path = [];
    end
end
if spotsrun.ctrl_ch > spotsrun.ctrl_chcount
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Control RNA/signal channel index is invalid. Control will not be used."), true);
    spotsrun.ctrl_path = [];
end

%Call core
verbosity = 1;
if bPreloaded
    RNA_Pipeline_Core(spotsrun, verbosity, preloaded_images);
else
    RNA_Pipeline_Core(spotsrun, verbosity);
end

RNA_Fisher_State.clearStaticState();

end