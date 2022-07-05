%

%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Paths ==========================

InputDir = [BaseDir '\data\bigfish\mESC_4d\Tsix-AF594'];
RefStem = [BaseDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d'];

BFStem = [InputDir filesep 'BIGFISH_Tsix-AF594'];

% ========================== Params ==========================

tmin = 10;
tmax = 400;
trange = [tmin:1:tmax];
T = size(trange,2);

zmin = 14;
bf_thresh = 191;

madf = -0.25;
winsz = 15;

% ========================== Read BIGFISH Data ==========================
addpath('./core');
spot_counts = zeros(T,2);
spot_coords = cell(T,1);

%Read cellseg masks
cell_n = uint16(csvread([InputDir filesep 'cell_n.csv']));
cell_t = uint16(csvread([InputDir filesep 'cell_t.csv']));
nuc_n = uint16(csvread([InputDir filesep 'nuc_n.csv']));
nuc_t = uint16(csvread([InputDir filesep 'nuc_t.csv']));

save([BFStem '_bfcellseg'], 'cell_n', 'cell_t', 'nuc_n', 'nuc_t');

%Read spot detection results
%Remember to add 1 to all BIGFISH coords (I think it is 0 - based)
i = 1;
for t = tmin:tmax
    filetblname = [InputDir filesep 'spots_' sprintf('%04d', t) '.csv'];
    if isfile(filetblname)
        raw_coord_table = uint16(csvread(filetblname));
        spot_counts(i,1) = t;
        spot_counts(i,2) = size(raw_coord_table,1);
        scount = spot_counts(i,2);
        
        tcoords = uint16(zeros(scount,3));
        tcoords(:,1) = raw_coord_table(:,3) + 1;
        tcoords(:,2) = raw_coord_table(:,2) + 1;
        tcoords(:,3) = raw_coord_table(:,1) + 1 + zmin;
        
        spot_coords{i} = tcoords;
    else
        break;
    end
    i = i+1;
end

% ========================== Save BIGFISH Data for MATLAB ==========================

coord_table = spot_coords;
spot_table = spot_counts;
save([BFStem '_coordTable'], 'coord_table');
save([BFStem '_spotTable'], 'spot_table');

% ========================== See what threshold finder calls ==========================

[threshold, win_out, score_thresh, scanst] = RNA_Threshold_Common.estimateThreshold(spot_counts, [], winsz, 0.5, madf);
fprintf("Threshold w/ MADFactor = %.3f, WinSize = %d: %d\n", madf, winsz, threshold);

% ========================== Import Ref Set & Comparison Data ==========================

fprintf("Loading reference SpotSelector data...\n");
src_selector = RNA_Threshold_SpotSelector.openSelector(RefStem);
[~, bf_selector] = src_selector.makeCopy();

fprintf("Loading reference spot table...\n");
load([RefStem '_spotTable'], 'spot_table');
spot_counts_ref = spot_table;

fprintf("Importing new data into anno obj...\n");
bf_selector.save_stem = BFStem;
bf_selector = bf_selector.loadNewSpotset(spot_counts, spot_coords);
bf_selector.toggle_del_unsnapped = false;
bf_selector = bf_selector.refSnapToAutoSpots();
bf_selector = bf_selector.updateFTable();
bf_selector.saveMe();

% ========================== Export Plots ==========================

plotdir = [InputDir filesep 'plots'];
mkdir(plotdir);

%Spots circled image in BF set at BF and homebrew thresholds
bf_selector.toggle_singleSlice = false; %Set to max proj
bf_selector.threshold_idx = bf_thresh - tmin + 1;
bf_selector = bf_selector.drawImages();
saveas(bf_selector.fh_filter, [plotdir filesep 'circlespots_bfthresh.png']);
bf_selector.threshold_idx = threshold - tmin + 1;
bf_selector = bf_selector.drawImages();
saveas(bf_selector.fh_filter, [plotdir filesep 'circlespots_hbthresh.png']);
close(bf_selector.fh_filter);
close(bf_selector.fh_raw);
if bf_selector.slice_drawn ~= 0
    close(bf_selector.fh_origslice);
end
if bf_selector.mode_3d
    close(bf_selector.fh_3d);
end

%Spots vs. thresh in homebrew spot set and BF spot set
legend_names = cell(1, 2);
color1 = [0.849 0.633 0.778]; %#F2BBE0
color2 = [0.759 0.665 0.802]; %#DBC3E6
color3 = [0.416 0.427 0.678]; %#6A6DAD
color4 = [0.416 0.588 0.678]; %#6A96AD
grey = [0.608 0.621 0.628]; %#A8ABAD
dim = [.45, .4, .5, .5];
figh = figure(989);
clf;
ax = axes;        
plot(spot_counts(:,1),log10(spot_counts(:,2)),'LineWidth',2,'Color',color1);
hold on;
plot(spot_counts_ref(:,1),log10(spot_counts_ref(:,2)),'-.','LineWidth',2,'Color',color2);
line([bf_thresh bf_thresh], get(ax,'YLim'),'Color',grey,'LineStyle','--');
legend_names{1,1} = 'BIG-FISH';
legend_names{1,2} = 'Homebrew';
legend(legend_names);
xlabel('Threshold');
ylabel('log10(# Spots Detected)');
saveas(figh, [plotdir filesep 'logspots.png']);
close(figh);
            
%Fscore curves of both spot sets
src_selector = src_selector.updateFTable();

%Window score plot
%Cellseg masks visualized

