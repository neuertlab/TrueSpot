%
%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%BaseDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
addpath('./test');
addpath('./test/datadump');

% ========================== Constants ==========================

START_INDEX = 1;
END_INDEX = 1000;

ResultsDir = [BaseDir filesep 'data' filesep 'results'];

% ========================== Load csv Table ==========================

%InputTablePath = [BaseDir filesep 'test_images_simytc.csv'];
InputTablePath = [BaseDir filesep 'test_images_simvarmass.csv'];
%InputTablePath = [BaseDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

% ========================== Iterate through table entries ==========================

entry_count = size(image_table,1);

if START_INDEX < 1; START_INDEX = 1; end
if END_INDEX > entry_count; END_INDEX = entry_count; end

for r = START_INDEX:END_INDEX

    myname = getTableValue(image_table, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);

    hb_stem_base = getTableValue(image_table, r, 'OUTSTEM');
    hb_stem = [BaseDir replace(hb_stem_base, '/', filesep)];
    hbq_stem = replace(hb_stem, '_all_3d', '');

    GroupSummaryDir = [ResultsDir filesep getSetOutputDirName(myname)];
    SummaryFilePath = [GroupSummaryDir filesep myname '_summary.mat'];

    if ~isfile(SummaryFilePath)
        fprintf('Summary data not found for %s! Skipping...\n', myname);
        continue;
    end

    %Check for quant file
    QuantReportPath = [hbq_stem '_quantData.mat'];
    if ~isfile(QuantReportPath)
        QuantReportPath = [hbq_stem '_quantData_varThresh.mat'];
        if ~isfile(QuantReportPath)
            fprintf('Quant data not found for %s! Skipping...\n', myname);
            continue;
        end
    end

    load(SummaryFilePath, 'analysis');
    load(QuantReportPath, 'quant_results');
    cell_results = quant_results.cell_rna_data;
    clear quant_results;

    image_size = [analysis.image_dims.y analysis.image_dims.x analysis.image_dims.z];

    if isfield(analysis, 'results_hb')
        analysis.results_hb.callset = ...
            RNACoords.addFitDataFromQuant(analysis.results_hb.callset, cell_results, analysis.results_hb.threshold, image_size);

        %Get ref set
        ref_coords = [];
        if isfield(analysis, 'simkey')
            if isstruct(analysis.simkey)
                refspotcount = size(analysis.simkey,2);
                ref_coords = NaN(refspotcount,3);
                ref_coords(:,1) = [analysis.simkey.x];
                ref_coords(:,2) = [analysis.simkey.y];
                ref_coords(:,3) = [analysis.simkey.z];
            else
                ref_coords = analysis.simkey;
            end
        elseif isfield(analysis, 'exprefset')
            ref_coords = analysis.exprefset;
        end

        if ~isempty(ref_coords)
            analysis.results_hb.callset =...
                RNACoords.updateRefDistancesToUseFits(analysis.results_hb.callset,...
                ref_coords, analysis.results_hb.ref_call_map);
        end

    end

    save(SummaryFilePath, 'analysis');
    clear analysis cell_results ref_coords;
end

% ========================== Helper functions ==========================

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