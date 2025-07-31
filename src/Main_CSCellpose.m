%
%%

function Main_CSCellpose(varargin)

addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2025.07.31.00';
VERSION_STRING = 'v1.3.1';

% ========================== Process args ==========================

arg_debug = true; %CONSTANT used for debugging arg parser.
cellpose_options = struct();
cellpose_options.overwrite_output = false;
cellpose_options.dump_summary = false;
cellpose_options.input_path = [];
cellpose_options.input_nuc = [];
cellpose_options.imgname = [];
cellpose_options.output_path = [];
cellpose_options.outpath_cell_mask = [];
cellpose_options.outpath_nuc_mask = [];
cellpose_options.outpath_settings = [];
cellpose_options.total_ch = 1;
cellpose_options.total_ch_nuc = 0;
cellpose_options.ch_light = 0;
cellpose_options.ch_nuc = 0;
cellpose_options.use_template = [];
cellpose_options.save_template_as = [];

cellpose_settings = CellPoseTS.genCellposeParamStruct();

override_checklist = struct();
OVERRIDE_OPS = {'cnorm' 'nnorm' 'censemble' 'nensemble' 'voxelsize' 'cmodel' 'nmodel' ...
    'cavgdia' 'navgdia' 'ccth' 'ncth' 'cfth' 'nfth' 'cszmin' 'nszmin' ...
    'czmin' 'czmax' 'nzmin' 'nzmax'};
ovrcount = size(OVERRIDE_OPS, 2);
for i = 1:ovrcount
    override_checklist.(OVERRIDE_OPS{i}) = false;
end
clear i

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
            cellpose_options.overwrite_output = true;
            if arg_debug; fprintf("Overwrite Output: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "dumpsummary")
            cellpose_options.dump_summary = true;
            if arg_debug; fprintf("Dump Parameter Summary: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "cnorm")
            cellpose_settings.cyto_params.normalize_bool = true;
            if arg_debug; fprintf("Normalize Image for Cell/Cyto Analysis: On\n"); end
            override_checklist.(lastkey) = true;
            lastkey = [];
        elseif strcmp(lastkey, "nnorm")
            cellpose_settings.nuc_params.normalize_bool = true;
            if arg_debug; fprintf("Normalize Image for Nuc Analysis: On\n"); end
            override_checklist.(lastkey) = true;
            lastkey = [];
        elseif strcmp(lastkey, "norm")
            cellpose_settings.cyto_params.normalize_bool = true;
            cellpose_settings.nuc_params.normalize_bool = true;
            if arg_debug; fprintf("Normalize Image for Cell & Nuc Analysis: On\n"); end
            override_checklist.cnorm = true;
            override_checklist.nnorm = true;
            lastkey = [];
        elseif strcmp(lastkey, "censemble")
            cellpose_settings.cyto_params.ensemble_bool = true;
            if arg_debug; fprintf("Ensemble Model for Cell/Cyto Analysis: On\n"); end
            override_checklist.(lastkey) = true;
            lastkey = [];
        elseif strcmp(lastkey, "nensemble")
            cellpose_settings.nuc_params.ensemble_bool = true;
            if arg_debug; fprintf("Ensemble Model for Nuc Analysis: On\n"); end
            override_checklist.(lastkey) = true;
            lastkey = [];
        elseif strcmp(lastkey, "ensemble")
            cellpose_settings.cyto_params.ensemble_bool = true;
            cellpose_settings.nuc_params.ensemble_bool = true;
            if arg_debug; fprintf("Ensemble Model for Cell & Nuc Analysis: On\n"); end
            override_checklist.censemble = true;
            override_checklist.nensemble = true;
            lastkey = [];
        end
        
    else
        if isempty(lastkey)
            fprintf("Value without key: %s - Skipping...\n", argval);
            continue;
        end
        
        %Value
        if strcmp(lastkey, "input")
            cellpose_options.input_path = argval;
            if arg_debug; fprintf("Input Path Set: %s\n", cellpose_options.input_path); end
        elseif strcmp(lastkey, "innuc")
            cellpose_options.input_nuc = argval;
            if arg_debug; fprintf("Input Path (Nuclear Dye Channel) Set: %s\n", cellpose_options.input_nuc); end
        elseif strcmp(lastkey, "imgname")
            cellpose_options.imgname = argval;
            if arg_debug; fprintf("Image Name Set: %s\n", cellpose_options.imgname); end
        elseif strcmp(lastkey, "outpath")
            cellpose_options.output_path = argval;
            if arg_debug; fprintf("Output Path Set: %s\n", cellpose_options.output_path); end
        elseif strcmp(lastkey, "ocellmask")
            cellpose_options.outpath_cell_mask = argval;
            if arg_debug; fprintf("Cell Mask TIF Output Path Set: %s\n", cellpose_options.outpath_cell_mask); end
        elseif strcmp(lastkey, "onucmask")
            cellpose_options.outpath_nuc_mask = argval;
            if arg_debug; fprintf("Nucleus Mask TIF Output Path Set: %s\n", cellpose_options.outpath_nuc_mask); end
        elseif strcmp(lastkey, "osettings")
            cellpose_options.outpath_settings = argval;
            cellpose_options.dump_summary = true;
            if arg_debug; fprintf("Settings Output Path Set: %s\n", cellpose_options.outpath_settings); end
        elseif strcmp(lastkey, "chtotal")
            cellpose_options.total_ch = Force2Num(argval);
            if arg_debug; fprintf("Primary Input Channel Count Set: %d\n", cellpose_options.total_ch); end
        elseif strcmp(lastkey, "chtotnuc")
            cellpose_options.total_ch_nuc = Force2Num(argval);
            if arg_debug; fprintf("Nuclear Input Channel Count Set: %d\n", cellpose_options.total_ch_nuc); end
        elseif strcmp(lastkey, "chlight")
            cellpose_options.ch_light = Force2Num(argval);
            if arg_debug; fprintf("Light Channel Index Set: %d\n", cellpose_options.ch_light); end
        elseif strcmp(lastkey, "chnuc")
            cellpose_options.ch_nuc = Force2Num(argval);
            if arg_debug; fprintf("Nuc Channel Index Set: %d\n", cellpose_options.ch_nuc); end
        elseif strcmp(lastkey, "voxelsize")
            idims_voxel = parseDimsTo(argval, cellpose_options.idims_voxel);
            if idims_voxel.z > 0
                cellpose_settings.z2xy = idims_voxel.z / idims_voxel.x;
                if arg_debug; fprintf("Z to XY Anisotropy Ratio Set: %.3f\n", cellpose_settings.z2xy); end
            end
            override_checklist.(lastkey) = true;
            clear idims_voxel
        elseif strcmp(lastkey, "cmodel")
            cellpose_settings.cyto_params.model_name = argval;
            if arg_debug; fprintf("Cell/Cyto Model Name Set: %s\n", cellpose_settings.cyto_params.model_name); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "nmodel")
            cellpose_settings.nuc_params.model_name = argval;
            if arg_debug; fprintf("Nuc Model Name Set: %s\n", cellpose_settings.nuc_params.model_name); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "cavgdia")
            cellpose_settings.cyto_params.avg_dia = Force2Num(argval);
            if arg_debug; fprintf("Cell Average Diameter Set (Pixels): %d\n", cellpose_settings.cyto_params.avg_dia); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "navgdia")
            cellpose_settings.nuc_params.avg_dia = Force2Num(argval);
            if arg_debug; fprintf("Nuc Average Diameter Set (Pixels): %d\n", cellpose_settings.nuc_params.avg_dia); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "ccth")
            cellpose_settings.cyto_params.cell_threshold = Force2Num(argval);
            if arg_debug; fprintf("Cell/Cyto Cell Threshold Set: %d\n", cellpose_settings.cyto_params.cell_threshold); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "ncth")
            cellpose_settings.nuc_params.cell_threshold = Force2Num(argval);
            if arg_debug; fprintf("Nuc Cell Threshold Set: %d\n", cellpose_settings.nuc_params.cell_threshold); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "cfth")
            cellpose_settings.cyto_params.flow_threshold = Force2Num(argval);
            if arg_debug; fprintf("Cell/Cyto Flow Threshold Set: %.3f\n", cellpose_settings.cyto_params.flow_threshold); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "nfth")
            cellpose_settings.nuc_params.flow_threshold = Force2Num(argval);
            if arg_debug; fprintf("Nuc Flow Threshold Set: %.3f\n", cellpose_settings.nuc_params.flow_threshold); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "cszmin")
            cellpose_settings.cyto_params.min_size = Force2Num(argval);
            if arg_debug; fprintf("Min Cell Size Set: %d\n", cellpose_settings.cyto_params.min_size); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "nszmin")
            cellpose_settings.nuc_params.min_size = Force2Num(argval);
            if arg_debug; fprintf("Min Nuc Size Set: %d\n", cellpose_settings.nuc_params.min_size); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "czmin")
            cellpose_settings.czmin = Force2Num(argval);
            if arg_debug; fprintf("Cell Channel Z Min Set: %d\n", cellpose_settings.czmin); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "czmax")
            cellpose_settings.czmax = Force2Num(argval);
            if arg_debug; fprintf("Cell Channel Z Max Set: %d\n", cellpose_settings.czmax); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "nzmin")
            cellpose_settings.nzmin = Force2Num(argval);
            if arg_debug; fprintf("Nuc Channel Z Min Set: %d\n", cellpose_settings.nzmin); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "nzmax")
            cellpose_settings.nzmax = Force2Num(argval);
            if arg_debug; fprintf("Nuc Channel Z Max Set: %d\n", cellpose_settings.nzmax); end
            override_checklist.(lastkey) = true;
        elseif strcmp(lastkey, "template")
            cellpose_options.use_template = argval;
            if arg_debug; fprintf("Template Set: %s\n", cellpose_options.use_template); end
        elseif strcmp(lastkey, "savetmpl")
            cellpose_options.save_template_as = argval;
            if arg_debug; fprintf("Template Save Name Set: %s\n", cellpose_options.save_template_as); end
        else
            fprintf("Key not recognized: %s - Skipping...\n", lastkey);
        end
    end
end

%--- Print start message

fprintf("-------------- TRUESPOT Cellpose Wrapper --------------\n");
fprintf("TrueSpot Version: %s\n", VERSION_STRING);
fprintf("CSCellpose Version: %s\n", BUILD_STRING);
fprintf("Timestamp: %s\n", datetime);

%--- Check args (Fill in defaults based on inputs)
if isempty(cellpose_options.input_path)
    fprintf("Input file is required! Exiting... \n");
    return;
end

if isempty(cellpose_options.imgname)
    [~, fn, ~] = fileparts(cellpose_options.input_path);
    cellpose_options.imgname = fn;
    fprintf("Image name not specified. Set to %s\n", cellpose_options.imgname);
end

if isempty(cellpose_options.outpath)
    [indir, ~, ~] = fileparts(cellpose_options.input_path);
    cellpose_options.outpath = [indir filesep 'TSCP_' cellpose_options.imgname '.mat'];
    fprintf("Output path not specified. Set to %s\n", cellpose_options.outpath);
end

if ~cellpose_options.overwrite_output & isfile(cellpose_options.outpath)
    fprintf("Output file %s already exists! Exiting...\n", cellpose_options.outpath);
    return;
end

%--- Load template (if applicable)
if ~isempty(cellpose_options.use_template)
    cellpose_settings = loadTemplateInto(cellpose_options.use_template, cellpose_settings, override_checklist);
end

%--- Save template (if applicable)
if ~isempty(cellpose_options.save_template_as)
    saveTemplate(cellpose_options.save_template_as, cellpose_settings);
end

%--- Run Cellpose
nuc_channel = [];
nuc_lbl = [];
if cellpose_options.ch_nuc > 0
    fprintf("[%s] Attempting nuclear segmentation...\n", datetime);
    fprintf("[%s] Loading nuclear marker channel...\n", datetime);
    npath = cellpose_options.input_nuc;
    ntifch = cellpose_options.total_ch_nuc;
    if isempty(npath)
        npath = cellpose_options.input_path; 
        ntifch = cellpose_options.total_ch; 
    end

    [channels, ~] = LoadTif(npath, ntifch, cellpose_options.ch_nuc, 1);
    nuc_channel = channels{cellpose_options.ch_nuc, 1};
    clear channels

    nuc_lbl = CellPoseTS.runNuclearSegmentation(nuc_channel, cellpose_settings);
end

cell_lbl = [];
if cellpose_options.ch_light > 0
    fprintf("[%s] Attempting cell segmentation...\n", datetime);
    fprintf("[%s] Loading TRANS/cyto marker channel...\n", datetime);

    [channels, ~] = LoadTif(cellpose_options.input_path, cellpose_options.total_ch, cellpose_options.ch_light, 1);
    cyto_channel = channels{cellpose_options.ch_light, 1};
    clear channels

    cell_lbl = CellPoseTS.runCytoSegmentation(cyto_channel, nuc_channel, cellpose_settings);
end
fprintf("[%s] Segmentation attempt completed. Saving...\n", datetime);

%--- Save Results

runMeta = struct();
runMeta.modifiedDate = datetime;
runMeta.tsCellSegBuild = buildString;
runMeta.tsCellSegVersion = versionString;

runMeta.srcImage = cellpose_options.input_path;
runMeta.srcImageChTotal = cellpose_options.total_ch;
runMeta.srcImageChTrans = cellpose_options.ch_light;
if ~isempty(cellpose_options.input_nuc)
    runMeta.srcNucImage = cellpose_options.input_nuc;
    runMeta.srcNucImageChTotal = cellpose_options.total_ch_nuc;
end
runMeta.srcImageChNuc = cellpose_options.ch_nuc;

if isfolder(cellpose_options.output_path)
    %It's a directory. Generate a file name, then.
    outpath = [cellpose_options.output_path filesep 'TSCP_' cellpose_options.imgname '.mat'];
else
    [outdir, ~, ~] = fileparts(cellpose_options.output_path);
    if ~isfolder(outdir)
        mkdir(outdir);
    end
    outpath = cellpose_options.output_path;
end

cellSeg = struct('cell_mask', cell_lbl);
cellSeg.cell_info = CellSeg.getCellInfo(cell_lbl, nuc_lbl);
nucleiSeg = struct();
nucleiSeg.results = CellSeg.genNucSegResultsStruct();
nucleiSeg.results.nuc_label = max(double(nuc_lbl), [], 3, 'omitnan');
nucleiSeg.results.nuclei = uint16(nucleiSeg.results.nuc_label);
nucleiSeg.results.lbl_mid = (nuc_lbl ~= 0);
save(outpath, 'cellSeg', 'nucleiSeg', 'runMeta', 'cellpose_settings', '-v7.3');

clear cellSeg nucleiSeg



%--- Output Extra Stuff
if cellpose_options.dump_summary
    printSummary(cellpose_options, cellpose_settings);
end

%TODO Masks as tif or png


end %--- END MAIN

function dumpMaskToImageFile(outpath, maskdata)
%TODO
if endsWith(outpath, '.png')
    %TODO
elseif endsWith(outpath, '.tif')
    %TODO
end

end

function printSummary(options, cellpose_settings)
    if isempty(options.outpath_settings)
        options.outpath_settings = [options.output_path filesep 'cellpose_settings.txt'];
    end

    if ~options.overwrite_output & isfile(options.outpath_settings)
        return;
    end

    fileHandle = fopen(options.outpath_settings, 'w');
    fprintf(fileHandle, 'input_path=%s\n', options.input_path);
    fprintf(fileHandle, 'input_nuc=%s\n', options.input_nuc);
    fprintf(fileHandle, 'output_path=%s\n', options.output_path);
    fprintf(fileHandle, 'outpath_cell_mask=%s\n', options.outpath_cell_mask);
    fprintf(fileHandle, 'outpath_nuc_mask=%s\n', options.outpath_nuc_mask);
    fprintf(fileHandle, 'outpath_settings=%s\n', options.outpath_settings);
    fprintf(fileHandle, 'imgname=%s\n', options.imgname);
    fprintf(fileHandle, 'total_ch=%d\n', options.total_ch);
    fprintf(fileHandle, 'total_ch_nuc=%d\n', options.total_ch_nuc);
    fprintf(fileHandle, 'ch_nuc=%d\n', options.ch_nuc);
    fprintf(fileHandle, 'ch_light=%d\n', options.ch_light);
    fprintf(fileHandle, 'overwrite_output=%d\n', options.overwrite_output);
    fprintf(fileHandle, 'dump_summary=%d\n', options.dump_summary);
    fprintf(fileHandle, 'save_template_as=%s\n', options.save_template_as);
    fprintf(fileHandle, 'use_template=%s\n', options.use_template);

    fprintf(fileHandle, 'cellpose_settings.z2xy=%.3f\n', cellpose_settings.z2xy);
    fprintf(fileHandle, 'cellpose_settings.nzmin=%d\n', cellpose_settings.nzmin);
    fprintf(fileHandle, 'cellpose_settings.nzmax=%d\n', cellpose_settings.nzmax);
    fprintf(fileHandle, 'cellpose_settings.czmin=%d\n', cellpose_settings.czmin);
    fprintf(fileHandle, 'cellpose_settings.czmax=%d\n', cellpose_settings.czmax);

    fprintf(fileHandle, 'cellpose_settings.cyto_params.model_name=%s\n', cellpose_settings.cyto_params.model_name);
    fprintf(fileHandle, 'cellpose_settings.cyto_params.avg_dia=%d\n', cellpose_settings.cyto_params.avg_dia);
    fprintf(fileHandle, 'cellpose_settings.cyto_params.cell_threshold=%d\n', cellpose_settings.cyto_params.cell_threshold);
    fprintf(fileHandle, 'cellpose_settings.cyto_params.flow_threshold=%.3f\n', cellpose_settings.cyto_params.flow_threshold);
    fprintf(fileHandle, 'cellpose_settings.cyto_params.min_size=%d\n', cellpose_settings.cyto_params.min_size);
    fprintf(fileHandle, 'cellpose_settings.cyto_params.normalize_bool=%d\n', cellpose_settings.cyto_params.normalize_bool);
    fprintf(fileHandle, 'cellpose_settings.cyto_params.ensemble_bool=%d\n', cellpose_settings.cyto_params.ensemble_bool);
    fprintf(fileHandle, 'cellpose_settings.cyto_params.do3D=%d\n', cellpose_settings.cyto_params.do3D);

    fprintf(fileHandle, 'cellpose_settings.nuc_params.model_name=%s\n', cellpose_settings.nuc_params.model_name);
    fprintf(fileHandle, 'cellpose_settings.nuc_params.avg_dia=%d\n', cellpose_settings.nuc_params.avg_dia);
    fprintf(fileHandle, 'cellpose_settings.nuc_params.cell_threshold=%d\n', cellpose_settings.nuc_params.cell_threshold);
    fprintf(fileHandle, 'cellpose_settings.nuc_params.flow_threshold=%.3f\n', cellpose_settings.nuc_params.flow_threshold);
    fprintf(fileHandle, 'cellpose_settings.nuc_params.min_size=%d\n', cellpose_settings.nuc_params.min_size);
    fprintf(fileHandle, 'cellpose_settings.nuc_params.normalize_bool=%d\n', cellpose_settings.nuc_params.normalize_bool);
    fprintf(fileHandle, 'cellpose_settings.nuc_params.ensemble_bool=%d\n', cellpose_settings.nuc_params.ensemble_bool);
    fprintf(fileHandle, 'cellpose_settings.nuc_params.do3D=%d\n', cellpose_settings.nuc_params.do3D);
    fclose(fileHandle);
end

function idims = parseDimsTo(dims_arg, trg_struct)
    idims = trg_struct;
	if ischar(dims_arg)
        %Read as string
        %'(x,y,z)'
        repl_str = replace(dims_arg, {'(', ')'}, '');
        split_str = split(repl_str, ',');
        ndims = size(split_str,1);
        idims.x = Force2Num(split_str{1,1});
        if ndims > 1; idims.y = Force2Num(split_str{2,1}); end
        if ndims > 2; idims.z = Force2Num(split_str{3,1}); end
	else
        %Try to read as vector
        if isvector(dims_arg)
            ndims = size(dims_arg,2);
            idims.x = dims_arg(1);
            if ndims > 1; idims.y = dims_arg(2); end
            if ndims > 2; idims.z = dims_arg(3); end
        end
	end
end

function cellpose_settings = loadTemplateInto(template_name, cellpose_settings, override_checklist)
    templateDir = './cellsegTemplates/cellpose';
    if ~isfolder(templateDir)
        fprintf("WARNING: Template %s requested, but template directory does not exist!\n", template_name);
        return;
    end
    savePath = [templateDir '/' template_name '.mat'];
    if ~isfile(savePath)
        fprintf("WARNING: Template %s requested, but template file %s does not exist!\n", template_name, savePath);
    end
    load(savePath, 'params');
    t_cyto = params.cyto_params;
    t_nuc = params.nuc_params;

    %Copy into cellpose_settings if NOT overridden
    if ~override_checklist.cnorm; cellpose_settings.cyto_params.normalize_bool = t_cyto.normalize_bool; end
    if ~override_checklist.nnorm; cellpose_settings.nuc_params.normalize_bool = t_nuc.normalize_bool; end
    if ~override_checklist.censemble; cellpose_settings.cyto_params.ensemble_bool = t_cyto.ensemble_bool; end
    if ~override_checklist.nensemble; cellpose_settings.nuc_params.ensemble_bool = t_nuc.ensemble_bool; end
    if ~override_checklist.voxelsize; cellpose_settings.z2xy = params.z2xy; end
    if ~override_checklist.cmodel; cellpose_settings.cyto_params.model_name = t_cyto.model_name; end
    if ~override_checklist.nmodel; cellpose_settings.nuc_params.model_name = t_nuc.model_name; end
    if ~override_checklist.cavgdia; cellpose_settings.cyto_params.avg_dia = t_cyto.avg_dia; end
    if ~override_checklist.navgdia; cellpose_settings.nuc_params.avg_dia = t_nuc.avg_dia; end
    if ~override_checklist.ccth; cellpose_settings.cyto_params.cell_threshold = t_cyto.cell_threshold; end
    if ~override_checklist.ncth; cellpose_settings.nuc_params.cell_threshold = t_nuc.cell_threshold; end
    if ~override_checklist.cfth; cellpose_settings.cyto_params.flow_threshold = t_cyto.flow_threshold; end
    if ~override_checklist.nfth; cellpose_settings.nuc_params.flow_threshold = t_nuc.flow_threshold; end
    if ~override_checklist.cszmin; cellpose_settings.cyto_params.min_size = t_cyto.min_size; end
    if ~override_checklist.nszmin; cellpose_settings.nuc_params.min_size = t_nuc.min_size; end
    if ~override_checklist.czmin; cellpose_settings.czmin = params.czmin; end
    if ~override_checklist.czmax; cellpose_settings.czmax = params.czmax; end
    if ~override_checklist.nzmin; cellpose_settings.nzmin = params.nzmin; end
    if ~override_checklist.nzmax; cellpose_settings.nzmax = params.nzmax; end
end

function saveTemplate(template_name, cellpose_settings)
    templateDir = './cellsegTemplates/cellpose';
    if ~isfolder(templateDir)
        mkdir(templateDir);
    end
    savePath = [templateDir '/' template_name '.mat'];
    params = cellpose_settings; clear cellpose_settings; %Change var name
    ModifiedTime = datetime;
    save(savePath, 'params', 'ModifiedTime');
end