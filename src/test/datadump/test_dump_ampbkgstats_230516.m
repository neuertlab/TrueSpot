%
%%  !! UPDATE TO YOUR BASE DIR
%BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
BaseDir = 'D:\usr\bghos\labdat\imgproc';

%ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

%ImgDir = 'C:\Users\hospelb\labdata\imgproc';
ImgDir = 'D:\usr\bghos\labdat\imgproc';

TblOutDir = [BaseDir filesep 'tables'];
ResultsDir = [BaseDir filesep 'data' filesep 'results'];

addpath('./core');
addpath('./test');

% ========================== Constants ==========================

TablePath_Main = [BaseDir filesep 'test_images.csv'];
OutTablePath = [TblOutDir filesep 'stats_est_bkgamp.tsv'];

ImageTableCols = {'IMGNAME', 'AMPMIN', 'AMPMAX', 'AMPAVG', 'AMPSTD',...
    'BKGMIN', 'BKGMAX', 'BKGAVG', 'BKGSTD'};

ImageTableColCount = size(ImageTableCols,2);

% ========================== Prep ==========================

image_table = testutil_opentable(TablePath_Main);

if ~isfolder(TblOutDir)
    mkdir(TblOutDir);
end

OutTableFile = fopen(OutTablePath, 'w');

for i = 1:ImageTableColCount
    if i ~= 1; fprintf(OutTableFile, '\t'); end
    fprintf(OutTableFile, ImageTableCols{i});
end
fprintf(OutTableFile, '\n');

% ========================== Do Things ==========================

entry_count = size(image_table,1);
for r = 1:entry_count

    myname = getTableValue(image_table, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);
    if startsWith(myname, 'sim')
        fprintf('Skipping sim image...\n');
        continue;
    elseif startsWith(myname, 'rsfish_sim')
        fprintf('Skipping sim image...\n');
        continue;
    end

    %Get HB stem and bkg file path
    hb_stem_raw = getTableValue(image_table, r, 'OUTSTEM');
    hb_stem = [BaseDir replace(hb_stem_raw, '/', filesep)];
    [hb_dir, ~, ~] = fileparts(hb_stem);
    bkg_path = [hb_dir filesep 'bkgmask' filesep myname '_bkg.mat'];

    %Get res file path
    set_group_dir = getSetOutputDirName(myname);
    ResFilePath = [ResultsDir filesep set_group_dir filesep myname '_summary.mat'];

    %Get amplitude values
    a_min = NaN;
    a_max = NaN;
    a_avg = NaN;
    a_std = NaN;
    has_amp_data = false;
    loaded_channel = [];
    if isfile(ResFilePath)
        load(ResFilePath, 'analysis');
        if isfield(analysis, 'truthset_BH')
            loaded_channel = loadTIFChannel(image_table, r);
            if ~isempty(loaded_channel)
                has_amp_data = true;
                ref_coords = analysis.truthset_BH.exprefset;
                coords_1D = sub2ind(size(loaded_channel),...
                    ref_coords(:,1).', ref_coords(:,2).', ref_coords(:,3).');

                intensities = double(loaded_channel(coords_1D));
                a_min = min(intensities, [], 'all', 'omitnan');
                a_max = max(intensities, [], 'all', 'omitnan');
                a_avg = mean(intensities, 'all', 'omitnan');
                a_std = std(intensities, 0, 'all', 'omitnan');
                clear coords_1D ref_coords intensities
            end
        end
        clear analysis;
    end

    %Get background values
    b_min = NaN;
    b_max = NaN;
    b_avg = NaN;
    b_std = NaN;
    has_bkg_data = false;
    if isfile(bkg_path)
        if isempty(loaded_channel)
            loaded_channel = loadTIFChannel(image_table, r);
        end
        if ~isempty(loaded_channel)
            load(bkg_path, 'background_mask');
            Z = size(loaded_channel,3);
            Y = size(loaded_channel,1);
            X = size(loaded_channel,2);
            mask3 = false(Y,X,Z);
            for z = 1:Z
                mask3(:,:,z) = background_mask;
            end

            all_bkg_vox = loaded_channel(find(mask3));
            b_min = min(all_bkg_vox, [], 'all', 'omitnan');
            b_max = max(all_bkg_vox, [], 'all', 'omitnan');
            b_avg = mean(all_bkg_vox, 'all', 'omitnan');
            b_std = std(all_bkg_vox, 0, 'all', 'omitnan');
            clear background_mask mask3 all_bkg_vox X Y Z
            has_bkg_data = true;
        end
    end

    clear loaded_channel;

    %Write to table
    if has_amp_data | has_bkg_data
        fprintf(OutTableFile, '%s\t', myname);
        fprintf(OutTableFile, '%f\t%f\t%f\t%f\t', a_min, a_max, a_avg, a_std);
        fprintf(OutTableFile, '%f\t%f\t%f\t%f\n', b_min, b_max, b_avg, b_std);
    end
end

fclose(OutTableFile);

% ========================== Helper Functions ==========================

function channel = loadTIFChannel(mytable, row_index)
    tif_path_raw = getTableValue(mytable, row_index, 'IMAGEPATH');
    tif_path = [ImgDir replace(tif_path_raw, '/', filesep)];

    ch_tot = getTableValue(mytable, row_index, 'CH_TOTAL');
    ch_samp = getTableValue(mytable, row_index, 'CHANNEL');

    allch = LoadTif(tif_path, ch_tot, [ch_samp], 0);
    channel = allch{ch_samp, 1};

    clear allch;
end

function dirname = getSetOutputDirName(imgname)
    inparts = split(imgname, '_');
    groupname = inparts{1,1};
    if strcmp(groupname, 'sctc')
        dirname = [groupname filesep inparts{2,1}];
    elseif strcmp(groupname, 'simvarmass')
        if contains(imgname, 'TMRL') | contains(imgname, 'CY5L')
            dirname = 'simytc';
        else
            dirname = groupname;
        end
    else
        if startsWith(groupname, 'ROI')
            dirname = 'munsky_lab';
        else
            dirname = groupname;
        end
    end
end

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end