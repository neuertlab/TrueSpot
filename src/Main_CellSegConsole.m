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
        elseif strcmp(lastkey, "savebk")
            cellseg_options.save_fmt_BK = true;
            if arg_debug; fprintf("Use CellDissect Output Organization: On\n"); end
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
        elseif strcmp(lastkey, "outpath")
            cellseg_options.output_path = argval;
            if arg_debug; fprintf("Output Path Set: %s\n", cellseg_options.output_path); end
        elseif strcmp(lastkey, "ocellmask")
            cellseg_options.outpath_cell_mask = argval;
            if arg_debug; fprintf("Cell Mask TIF Output Path Set: %s\n", cellseg_options.outpath_cell_mask); end
        elseif strcmp(lastkey, "onucmask")
            cellseg_options.outpath_nuc_mask = argval;
            if arg_debug; fprintf("Nucleus Mask TIF Output Path Set: %s\n", cellseg_options.outpath_nuc_mask); end
        elseif strcmp(lastkey, "osettings")
            cellseg_options.outpath_settings = argval;
            cellseg_options.dump_summary = true;
            if arg_debug; fprintf("Settings Output Path Set: %s\n", cellseg_options.outpath_settings); end
        elseif strcmp(lastkey, "chtotal")
            cellseg_options.total_ch = Force2Num(argval);
            if arg_debug; fprintf("Primary Input Channel Count Set: %d\n", cellseg_options.total_ch); end
        elseif strcmp(lastkey, "chtotnuc")
            cellseg_options.total_ch_nuc = Force2Num(argval);
            if arg_debug; fprintf("Nuclear Input Channel Count Set: %d\n", cellseg_options.total_ch_nuc); end
        elseif strcmp(lastkey, "chlight")
            cellseg_options.ch_light = Force2Num(argval);
            if arg_debug; fprintf("Light Channel Index Set: %d\n", cellseg_options.ch_light); end
        elseif strcmp(lastkey, "chnuc")
            cellseg_options.ch_nuc = Force2Num(argval);
            if arg_debug; fprintf("Nuc Channel Index Set: %d\n", cellseg_options.ch_nuc); end
        elseif strcmp(lastkey, "cszmin")
            cellseg_options.cell_params.min_cell_size = Force2Num(argval);
            if arg_debug; fprintf("Min Cell Size Set: %d\n", cellseg_options.cell_params.min_cell_size); end
        elseif strcmp(lastkey, "cszmax")
            cellseg_options.cell_params.max_cell_size = Force2Num(argval);
            if arg_debug; fprintf("Max Cell Size Set: %d\n", cellseg_options.cell_params.max_cell_size); end
        elseif strcmp(lastkey, "fplstrat")
            cellseg_options.cell_params.focus_plane_strat = argval;
            if arg_debug; fprintf("Focus Plane Strategy Set: %s\n", cellseg_options.cell_params.focus_plane_strat); end
        elseif strcmp(lastkey, "fzmin")
            cellseg_options.cell_params.min_plane = Force2Num(argval);
            if arg_debug; fprintf("Focus Z Min Set: %d\n", cellseg_options.cell_params.min_plane); end
        elseif strcmp(lastkey, "fzmax")
            cellseg_options.cell_params.max_plane = Force2Num(argval);
            if arg_debug; fprintf("Focus Z Max Set: %d\n", cellseg_options.cell_params.max_plane); end
        elseif strcmp(lastkey, "foffmin")
            cellseg_options.cell_params.focus_offset_min = Force2Num(argval);
            if arg_debug; fprintf("Focus Offset Min Set: %d\n", cellseg_options.cell_params.focus_offset_min); end
        elseif strcmp(lastkey, "foffmax")
            cellseg_options.cell_params.focus_offset_max = Force2Num(argval);
            if arg_debug; fprintf("Focus Offset Max Set: %d\n", cellseg_options.cell_params.focus_offset_max); end
        elseif strcmp(lastkey, "xtrim")
            cellseg_options.cell_params.x_trim = Force2Num(argval);
            if arg_debug; fprintf("X Trim Set: %d\n", cellseg_options.cell_params.x_trim); end
        elseif strcmp(lastkey, "ytrim")
            cellseg_options.cell_params.y_trim = Force2Num(argval);
            if arg_debug; fprintf("Y Trim Set: %d\n", cellseg_options.cell_params.y_trim); end
        elseif strcmp(lastkey, "nzrange")
            cellseg_options.nuc_params.range = Force2Num(argval);
            if arg_debug; fprintf("NucSeg Z Range Set: %d\n", cellseg_options.nuc_params.range); end
        elseif strcmp(lastkey, "nthsmpl")
            cellseg_options.nuc_params.threshold_sampling = Force2Num(argval);
            if arg_debug; fprintf("NucSeg Threshold Sampler Set: %d\n", cellseg_options.nuc_params.threshold_sampling); end
        elseif strcmp(lastkey, "nszmin")
            cellseg_options.nuc_params.min_nucleus_size = Force2Num(argval);
            if arg_debug; fprintf("Min Nucleus Size Set: %d\n", cellseg_options.nuc_params.min_nucleus_size); end
        elseif strcmp(lastkey, "nszmax")
            cellseg_options.nuc_params.max_nucleus_size = Force2Num(argval);
            if arg_debug; fprintf("Max Nucleus Size Set: %d\n", cellseg_options.nuc_params.max_nucleus_size); end
        elseif strcmp(lastkey, "ncutoff")
            cellseg_options.nuc_params.cutoff = Force2Num(argval);
            if arg_debug; fprintf("NucSeg Cutoff Factor: %f\n", cellseg_options.nuc_params.cutoff); end
        elseif strcmp(lastkey, "ndxy")
            cellseg_options.nuc_params.dxy = Force2Num(argval);
            if arg_debug; fprintf("NucSeg Radius Factor: %f\n", cellseg_options.nuc_params.dxy); end
        elseif strcmp(lastkey, "template")
            cellseg_options.use_template = argval;
            if arg_debug; fprintf("Template Set: %s\n", cellseg_options.use_template); end
        elseif strcmp(lastkey, "savetmpl")
            cellseg_options.save_template_as = argval;
            if arg_debug; fprintf("Template Save Name Set: %s\n", cellseg_options.save_template_as); end
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

if isempty(cellseg_options.output_path)
    %Use input
    [cellseg_options.output_path, ~, ~] = fileparts(cellseg_options.input_path);
    fprintf("Output directory path was not provided. Using input dir: %s...\n", cellseg_options.output_path);
end

if and(cellseg_options.dump_summary, isempty(cellseg_options.outpath_settings))
    cellseg_options.outpath_settings = [cellseg_options.output_dir filesep 'CSegSettings_' cellseg_options.imgname '.mat'];
    fprintf("Output settings summary path was not provided. Using: %s...\n", cellseg_options.outpath_settings);
end

if ~isempty(cellseg_options.use_template)
    %Try to find and load template
    [cellseg_options, okay] = loadParamTemplate(cellseg_options, cellseg_options.use_template);
    if ~okay
        fprintf('Template with name "%s" could not be found! Exiting...\n', cellseg_options.use_template);
        return;
    end
end

okay = runCellseg(cellseg_options);
if ~okay
    fprintf('ERROR: Cell segmentation run failed! Exiting...\n');
    return;
end

if ~isempty(cellseg_options.save_template_as)
    [cellseg_options, okay] = saveParamTemplate(cellseg_options, cellseg_options.save_template_as);
    if ~okay
        fprintf('WARNING: Save of template "%s" failed!\n', cellseg_options.save_template_as);
    end
end

end

% ========================== Helper Functions ==========================

function okay = runCellseg(options)
    okay = false;

    %Load image channels
    nuc_ch_data = [];
    light_ch_data = [];

    if ~isempty(options.input_nuc)
        %Nuc in its own file
        fprintf('> Loading nuclear marker channel...\n');
        [channels, ~] = LoadTif(options.input_nuc, options.total_ch_nuc, [options.ch_nuc], 1);
        if isempty(channels); return; end
        nuc_ch_data = channels{options.ch_nuc, 1};
        clear channels
    else
        %Load both at once
        fprintf('> Loading nuclear marker and light channels...\n');
        [channels, ~] = LoadTif(options.input_path,...
            options.total_ch, [options.ch_nuc options.ch_light], 1);
        if isempty(channels); return; end
        nuc_ch_data = channels{options.ch_nuc, 1};
        light_ch_data = channels{options.ch_light, 1};
        clear channels
    end

    %Attempt nuclear segmentation
    nucseg = CellSeg.AutosegmentNuclei(nuc_ch_data, options.nuc_params);
    if isempty(nucseg.nuc_label)
        fprintf('Nuclei segmentation failed!\n');
        return;
    end
    clear nuc_ch_data
    
    %Load light channel, if not alread loaded.
    if isempty(light_ch_data)
        fprintf('> Loading light channel...\n');
        [channels, ~] = LoadTif(options.input_path, options.total_ch, [options.ch_light], 1);
        if isempty(channels); return; end
        light_ch_data = channels{options.ch_light, 1};
        clear channels
    end

    %Attempt cell segmentation
    [cell_mask, cell_info, trans_plane, cellseg_info] = ...
        CellSeg.AutosegmentCells(light_ch_data, nucseg.nuc_label, options.cell_params);
    clear light_ch_data
    if isempty(cell_mask)
        fprintf('Cell segmentation failed!\n');
        return;
    end

    %Save
    if ~isempty(options.outpath_cell_mask)
        %Save cell mask TIF
        tiffops = struct('overwrite', true);
        saveastiff(cell_mask, options.outpath_cell_mask, tiffops);
        clear tiffops
    end
    if ~isempty(options.outpath_nuc_mask)
        %Save nuc mask TIF
        tiffops = struct('overwrite', true);
        saveastiff(nucseg.nuc_label, options.outpath_nuc_mask, tiffops);
        clear tiffops
    end

    if options.save_fmt_BK
        if ~isfolder(options.output_path)
            mkdir(options.output_path);
        end

        lab_path = [options.output_path filesep 'Lab_' options.imgname '.mat'];
        nuc_path = [options.output_path filesep 'nuclei_' options.imgname '.mat'];
        nucth_path = [options.output_path filesep 'nuclei_TH_' options.imgname '.mat'];

        fprintf('Saving cell mask...\n');
        if options.overwrite_output | ~isfile(lab_path)
            CellInfo = cell_info;
            cells = cell_mask;
            segment_mode = cellseg_info.focus_plane_strat;

            save(lab_path, 'CellInfo', 'cells', 'segment_mode', 'trans_plane', '-v7.3');

            clear CellInfo Cells segment_mode
        end

        fprintf('Saving nucleus mask...\n');
        if options.overwrite_output | ~isfile(nuc_path)
            Label_hi = nucseg.lbl_hi;
            Label_low = nucseg.lbl_lo;
            Label_mid = nucseg.lbl_mid;
            Maj_axis = nucseg.nuc_axis_major;
            Min_axis = nucseg.nuc_axis_minor;
            Nuc_int = nucseg.nuc_int;
            Nuc_vol = nucseg.nuc_vol;
            nuclei = nucseg.nuclei;

            save(nuc_path, 'Label_hi', 'Label_low', 'Label_mid',...
                'Maj_axis', 'Min_axis', 'Nuc_int', 'Nuc_vol', 'nuclei', '-v7.3');

            clear Label_hi Label_low Label_mid Maj_axis Min_axis Nuc_int Nuc_vol
        end

        fprintf('Saving nucleus threshold data...\n');
        if options.overwrite_output | ~isfile(nucth_path)
            DAPI_ims_added = nucseg.test_sum;
            dapi_label = nucseg.nuc_label;
            dapi_label_low_1 = nucseg.nuc_label_lo;
            dapi_threshold = nucseg.nuc_threshold;

            save(nucth_path, 'DAPI_ims_added', 'dapi_label', 'dapi_label_low_1',...
                'dapi_threshold', '-v7.3');

            clear DAPI_ims_added dapi_label dapi_label_low_1 dapi_threshold
        end
    else
        if isfolder(options.output_path)
            %It's a directory. Generate a file name, then.
            outpath = [options.output_path filesep 'CellSeg_' options.imgname '.mat'];
        else
            [outdir, ~, ~] = fileparts(options.output_path);
            if ~isfolder(outdir)
                mkdir(outdir);
            end
            outpath = options.output_path;
        end

        cellSeg = cellseg_info;
        cellSeg.cell_mask = cell_mask;
        cellSeg.cell_info = cell_info;
        cellSeg.trans_plane = trans_plane;
        nucleiSeg = nucseg;
        save(outpath, 'cellSeg', 'nucleiSeg', '-v7.3');

        clear cellSeg nucleiSeg
    end

    %Dump Settings
    if options.dump_summary
        printSummary(options);
    end

    okay = true;
end

function printSummary(options)
    if isempty(options.outpath_settings)
        return;
    end

    if ~options.overwrite_output & isfile(options.outpath_settings)
        return;
    end

    fileHandle = fopen(options.outpath_settings, 'w');
    fprintf('input_path=%s\n', options.input_path);
    fprintf('input_nuc=%s\n', options.input_nuc);
    fprintf('output_path=%s\n', options.output_path);
    fprintf('outpath_cell_mask=%s\n', options.outpath_cell_mask);
    fprintf('outpath_nuc_mask=%s\n', options.outpath_nuc_mask);
    fprintf('outpath_settings=%s\n', options.outpath_settings);
    fprintf('imgname=%s\n', options.imgname);
    fprintf('total_ch=%d\n', options.total_ch);
    fprintf('total_ch_nuc=%d\n', options.total_ch_nuc);
    fprintf('ch_nuc=%d\n', options.ch_nuc);
    fprintf('ch_light=%d\n', options.ch_light);
    fprintf('overwrite_output=%d\n', options.overwrite_output);
    fprintf('dump_summary=%d\n', options.dump_summary);
    fprintf('save_fmt_BK=%d\n', options.save_fmt_BK);
    fprintf('save_template_as=%s\n', options.save_template_as);
    fprintf('use_template=%s\n', options.use_template);
    fprintf('cell_params.focus_plane_strat=%s\n', options.cell_params.focus_plane_strat);
    fprintf('cell_params.min_cell_size=%d\n', options.cell_params.min_cell_size);
    fprintf('cell_params.max_cell_size=%d\n', options.cell_params.max_cell_size);
    fprintf('cell_params.min_plane=%d\n', options.cell_params.min_plane);
    fprintf('cell_params.max_plane=%d\n', options.cell_params.max_plane);
    fprintf('cell_params.focus_offset_min=%d\n', options.cell_params.focus_offset_min);
    fprintf('cell_params.focus_offset_max=%d\n', options.cell_params.focus_offset_max);
    fprintf('cell_params.x_trim=%d\n', options.cell_params.x_trim);
    fprintf('cell_params.y_trim=%d\n', options.cell_params.y_trim);
    fprintf('nuc_params.range=%d\n', options.nuc_params.range);
    fprintf('nuc_params.threshold_sampling=%d\n', options.nuc_params.threshold_sampling);
    fprintf('nuc_params.min_nucleus_size=%d\n', options.nuc_params.min_nucleus_size);
    fprintf('nuc_params.max_nucleus_size=%d\n', options.nuc_params.max_nucleus_size);
    fprintf('nuc_params.cutoff=%f\n', options.nuc_params.cutoff);
    fprintf('nuc_params.dxy=%d\n', options.nuc_params.dxy);
    fclose(fileHandle);
end

function [options, okay] = loadParamTemplate(options, templateId)
    okay = false;
    templateDir = './cellsegTemplates';
    if ~isfolder(templateDir)
        return;
    end

    templatePath = [templateDir filesep templateId '.mat'];
    if ~isfile(templatePath)
        return;
    end

    load(templatePath, 'cell_params', 'nuc_params');
    if isempty(cell_params); return; end
    if isempty(nuc_params); return; end

    options.cell_params = cell_params;
    options.nuc_params = nuc_params;

    okay = true;
end

function [options, okay] = saveParamTemplate(options, templateId)
    okay = false;
    templateDir = './cellsegTemplates';
    if ~isfolder(templateDir)
        mkdir(templateDir);
    end

    templatePath = [templateDir filesep templateId '.mat'];
    if ~options.overwrite_output & isfile(templatePath)
        return;
    end

    cell_params = options.cell_params;
    nuc_params = options.nuc_params;
    save(templatePath, 'cell_params', 'nuc_params', '-v7.3');

    okay = true;
end

function cellseg_options = genOptionsStruct()
    cellseg_options = struct('input_path', []);
    cellseg_options.input_nuc = [];
    cellseg_options.output_path = []; %If savebk, this is a dir. Otherwise, mat file.
    cellseg_options.outpath_cell_mask = []; %As TIF
    cellseg_options.outpath_nuc_mask = []; %As TIF
    cellseg_options.outpath_settings = [];

    cellseg_options.imgname = [];

    cellseg_options.total_ch = 1;
    cellseg_options.total_ch_nuc = 0;
    cellseg_options.ch_nuc = 0;
    cellseg_options.ch_light = 1;

    cellseg_options.overwrite_output = false;
    cellseg_options.dump_summary = false;
    cellseg_options.save_fmt_BK = false;

    cellseg_options.save_template_as = [];
    cellseg_options.use_template = [];

    cellseg_options.cell_params = CellSeg.genCellSegParameterStruct();
    cellseg_options.nuc_params = CellSeg.genNucSegStruct();
end