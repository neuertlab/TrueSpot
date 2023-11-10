%%
function Main_CellSegConsole(varargin)
addpath('./core');
addpath('./thirdparty');
addpath('./celldissect');

% ========================== Process args ==========================

arg_debug = true; %CONSTANT used for debugging arg parser.
cellseg_options = genOptionsStruct();

lastkey = [];
for i = 1:nargin
    argval = varargin{i};
    if ischar(argval) & startsWith(argval, "-")
        %Key
        if size(argval,2) >= 2
            lastkey = argval(2:end);
        else
            lastkey = [];
        end
        
        %Account for boolean keys...
        if strcmp(lastkey, "ovrw")
            cellseg_options.overwrite_output = true;
            if arg_debug; fprintf("Overwrite Output: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "dumpsummary")
            cellseg_options.dump_summary = true;
            if arg_debug; fprintf("Dump Parameter Summary: On\n"); end
            lastkey = [];
        end
        
    else
        if isempty(lastkey)
            fprintf("Value without key: %s - Skipping...\n", argval);
            continue;
        end
        
        %Value
        if strcmp(lastkey, "input")
            cellseg_options.input_path = argval;
            if arg_debug; fprintf("Input Path Set: %s\n", cellseg_options.input_path); end
        elseif strcmp(lastkey, "innuc")
            cellseg_options.input_nuc = argval;
            if arg_debug; fprintf("Input Path (Nuclear Dye Channel) Set: %s\n", cellseg_options.input_nuc); end
        elseif strcmp(lastkey, "imgname")
            cellseg_options.imgname = argval;
            if arg_debug; fprintf("Image Name Set: %s\n", cellseg_options.imgname); end
        elseif strcmp(lastkey, "outdir")
            cellseg_options.output_dir = argval;
            if arg_debug; fprintf("Output Directory Set: %s\n", cellseg_options.output_dir); end
        elseif strcmp(lastkey, "ocellmask")
            cellseg_options.outpath_cell_mask = argval;
            if arg_debug; fprintf("Cell Mask Output Path Set: %s\n", cellseg_options.outpath_cell_mask); end
        elseif strcmp(lastkey, "onucmask")
            cellseg_options.outpath_nuc_mask = argval;
            if arg_debug; fprintf("Nucleus Mask Output Path Set: %s\n", cellseg_options.outpath_nuc_mask); end
        elseif strcmp(lastkey, "onucth")
            cellseg_options.outpath_nuc_th = argval;
            if arg_debug; fprintf("Nucleus TH Output Path Set: %s\n", cellseg_options.outpath_nuc_th); end
        elseif strcmp(lastkey, "osettings")
            cellseg_options.outpath_settings = argval;
            cellseg_options.dump_summary = true;
            if arg_debug; fprintf("Settings Output Path Set: %s\n", cellseg_options.outpath_settings); end
        elseif strcmp(lastkey, "chtotal")
            cellseg_options.total_ch = Force2Num(argval);
            if arg_debug; fprintf("Primary Input Channel Count Set: %s\n", cellseg_options.total_ch); end
        elseif strcmp(lastkey, "chtotnuc")
            cellseg_options.total_ch_nuc = Force2Num(argval);
            if arg_debug; fprintf("Nuclear Input Channel Count Set: %s\n", cellseg_options.total_ch_nuc); end
        elseif strcmp(lastkey, "chlight")
            cellseg_options.ch_light = Force2Num(argval);
            if arg_debug; fprintf("Light Channel Index Set: %s\n", cellseg_options.ch_light); end
        elseif strcmp(lastkey, "chnuc")
            cellseg_options.ch_nuc = Force2Num(argval);
            if arg_debug; fprintf("Nuc Channel Index Set: %s\n", cellseg_options.ch_nuc); end
        else
            fprintf("Key not recognized: %s - Skipping...\n", lastkey);
        end
    end
end

%--- Check args (Fill in defaults based on inputs)
if isempty(cellseg_options.input_path)
    fprintf("ERROR: At least one input image (TIF format) is required! Exiting...\n");
    return;
end

if isempty(cellseg_options.imgname)
    %Just grab the file name.
    [~, cellseg_options.imgname, ~] = fileparts(cellseg_options.input_path);
    while contains(cellseg_options.imgname, ".")
        [~, cellseg_options.imgname, ~] = fileparts(cellseg_options.imgname);
    end
end

if isempty(cellseg_options.output_dir)
    %Use input
    [cellseg_options.output_dir, ~, ~] = fileparts(cellseg_options.input_path);
    fprintf("Output directory path was not provided. Using input dir: %s...\n", cellseg_options.output_dir);
end

if isempty(cellseg_options.outpath_cell_mask)
    cellseg_options.outpath_cell_mask = [cellseg_options.output_dir filesep 'Lab_' cellseg_options.imgname '.mat'];
    fprintf("Output cell mask path was not provided. Using: %s...\n", cellseg_options.outpath_cell_mask);
end

if isempty(cellseg_options.outpath_nuc_mask)
    cellseg_options.outpath_nuc_mask = [cellseg_options.output_dir filesep 'nuclei_' cellseg_options.imgname '.mat'];
    fprintf("Output nuc mask path was not provided. Using: %s...\n", cellseg_options.outpath_nuc_mask);
end

if isempty(cellseg_options.outpath_nuc_th)
    cellseg_options.outpath_nuc_th = [cellseg_options.output_dir filesep 'nuclei_TH_' cellseg_options.imgname '.mat'];
    fprintf("Output nuc TH data path was not provided. Using: %s...\n", cellseg_options.outpath_nuc_th);
end

if and(cellseg_options.dump_summary, isempty(cellseg_options.outpath_settings))
    cellseg_options.outpath_settings = [cellseg_options.output_dir filesep 'CSegSettings_' cellseg_options.imgname '.mat'];
    fprintf("Output settings summary path was not provided. Using: %s...\n", cellseg_options.outpath_settings);
end

end

% ========================== Helper Functions ==========================

function cellseg_options = genOptionsStruct()
    cellseg_options = struct('input_path', []);
    cellseg_options.input_nuc = [];
    cellseg_options.output_dir = [];
    cellseg_options.outpath_cell_mask = [];
    cellseg_options.outpath_nuc_mask = [];
    cellseg_options.outpath_nuc_th = [];
    cellseg_options.outpath_settings = [];

    cellseg_options.imgname = [];

    cellseg_options.total_ch = 1;
    cellseg_options.total_ch_nuc = 0;
    cellseg_options.ch_nuc = 0;
    cellseg_options.ch_light = 1;

    cellseg_options.overwrite_output = false;
    cellseg_options.dump_summary = false;
end