%

%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Paths ==========================

% InputDir = [BaseDir '\data\bigfish\mESC_4d\Tsix-AF594'];
% RefStem = [BaseDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_Tsix-AF594'];
% zmin = 14; zmax = 68;
% bf_thresh = 191;
% %bf_thresh = 36;

% InputDir = [BaseDir '\data\bigfish\mESC_4d\Xist-CY5'];
% RefStem = [BaseDir '\data\preprocess\feb2018\Xist_CY5\Xist\Xist-CY5_IMG1_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_Xist-CY5'];
% zmin = 14; zmax = 68;
% bf_thresh = 184;

% InputDir = [BaseDir '\data\bigfish\yeast_rna\E1R2I4\STL1'];
% RefStem = [BaseDir '\data\preprocess\YeastFISH\E1R2\STL1\E1R2-STL1-TMR_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_STL1-TMR'];
% zmin = 4; zmax = 23;
% bf_thresh = 154;

InputDir = [BaseDir '\data\bigfish\yeast_rna\E1R2I4\CTT1'];
RefStem = [BaseDir '\data\preprocess\YeastFISH\E1R2\CTT1\E1R2-CTT1-CY5_all_3d'];
BFStem = [InputDir filesep 'BIGFISH_CTT1-CY5'];
zmin = 1; zmax = 20;
% bf_thresh = 124;

% InputDir = [BaseDir '\data\bigfish\yeast_rna\E2R1I2\STL1'];
% RefStem = [BaseDir '\data\preprocess\YeastFISH\E2R1\STL1\E2R1-STL1-TMR_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_STL1-TMR'];
% zmin = 3;
% bf_thresh = 231;

% InputDir = [BaseDir '\data\bigfish\yeast_rna\E2R1I2\CTT1'];
% RefStem = [BaseDir '\data\preprocess\YeastFISH\E2R1\CTT1\E2R1-CTT1-CY5_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_CTT1-CY5'];
% zmin = 2;
% bf_thresh = 805;

% InputDir = [BaseDir '\data\bigfish\yeast_rna\E2R2I3\STL1'];
% RefStem = [BaseDir '\data\preprocess\YeastFISH\E2R2\Img3\STL1\E2R2I3-STL1-TMR_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_STL1-TMR'];
% zmin = 5;
% bf_thresh = 231;

% InputDir = [BaseDir '\data\bigfish\yeast_rna\E2R2I3\CTT1'];
% RefStem = [BaseDir '\data\preprocess\YeastFISH\E2R2\Img3\CTT1\E2R2I3-CTT1-CY5_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_CTT1-CY5'];
% zmin = 1;
% bf_thresh = 281;

% InputDir = [BaseDir '\data\bigfish\yeast_rna\E2R2I5\STL1'];
% RefStem = [BaseDir '\data\preprocess\YeastFISH\E2R2\Img5\STL1\E2R2I5-STL1-TMR_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_STL1-TMR'];
% zmin = 5;
% bf_thresh = 208;

% InputDir = [BaseDir '\data\bigfish\yeast_rna\E2R2I5\CTT1'];
% RefStem = [BaseDir '\data\preprocess\YeastFISH\E2R2\Img5\CTT1\E2R2I5-CTT1-CY5_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_CTT1-CY5'];
% zmin = 2;
% bf_thresh = 340;

% InputDir = [BaseDir '\data\bigfish\histones\D0I4\Xist'];
% RefStem = [BaseDir '\data\preprocess\histones\D0_I4\Ch2\Histone_D0_img4_ch2_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_Xist-CY5'];
% zmin = 17;
% bf_thresh = 250;

% InputDir = [BaseDir '\data\bigfish\histones\D0I4\Tsix'];
% RefStem = [BaseDir '\data\preprocess\histones\D0_I4\Ch3\Histone_D0_img4_ch3_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_Tsix-TMR'];
% zmin = 17;
% bf_thresh = 246;

% InputDir = [BaseDir '\data\bigfish\histones\D0I4\H3K36me3'];
% RefStem = [BaseDir '\data\preprocess\histones\D0_I4\Ch4\Histone_D0_img4_ch4_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_H3K36me3-AF488'];
% zmin = 17;
% bf_thresh = 356;

% InputDir = [BaseDir '\data\bigfish\histones\D0I6\Xist'];
% RefStem = [BaseDir '\data\preprocess\histones\D0_I6\Ch2\Histone_D0_img6_ch2_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_Xist-CY5'];
% zmin = 17;
% bf_thresh = 226;

% InputDir = [BaseDir '\data\bigfish\histones\D0I6\Tsix'];
% RefStem = [BaseDir '\data\preprocess\histones\D0_I6\Ch3\Histone_D0_img6_ch3_all_3d'];
% BFStem = [InputDir filesep 'BIGFISH_Tsix-TMR'];
% zmin = 14; zmax = 80;
% bf_thresh = 208;

% ========================== Params ==========================

tmin = 10;
tmax = 1000;
trange = [tmin:1:tmax];
T = size(trange,2);

% ========================== Read BIGFISH Data ==========================
addpath('./core');
spot_counts = zeros(T,2);
spot_coords = cell(T,1);

%Read cellseg masks
fprintf("Reading BIG-FISH cell seg data...\n");
%cell_n = uint16(csvread([InputDir filesep 'cell_n.csv']));
cell_t = uint16(csvread([InputDir filesep 'cell_t.csv']));
%nuc_n = uint16(csvread([InputDir filesep 'nuc_n.csv']));
nuc_t = uint16(csvread([InputDir filesep 'nuc_t.csv']));

%save([BFStem '_bfcellseg'], 'cell_n', 'cell_t', 'nuc_n', 'nuc_t');
save([BFStem '_bfcellseg_tonly'], 'cell_t', 'nuc_t');

%Read spot detection results
%Remember to add 1 to all BIGFISH coords (I think it is 0 - based)
fprintf("Reading BIG-FISH spot detection data...\n");
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

%[thresh_hbbf, win_out, score_thresh, scanst] = RNA_Threshold_Common.estimateThreshold(spot_counts, [], winsz, 0.5, madf);
%fprintf("Threshold w/ MADFactor = %.3f, WinSize = %d: %d\n", madf, winsz, thresh_hbbf);
thresh_hbbf = RNAThreshold.runDefaultParameters(spot_counts);
save([BFStem '_threshres.mat'], 'thresh_hbbf');

% ========================== Import Ref Set & Comparison Data ==========================
%TODO - Is the selector properly masking out spots in zslices trimmed by bf
%   in the reference set?

fprintf("Loading reference SpotSelector data...\n");
src_selector = RNA_Threshold_SpotSelector.openSelector(RefStem, true);
[~, bf_selector] = src_selector.makeCopy();

fprintf("Loading reference spot table...\n");
load([RefStem '_spotTable'], 'spot_table');
spot_counts_ref = spot_table;

fprintf("Importing new data into anno obj...\n");
bf_selector.save_stem = BFStem;
bf_selector = bf_selector.loadNewSpotset(spot_counts, spot_coords);
bf_selector.toggle_del_unsnapped = false;
bf_selector = bf_selector.refSnapToAutoSpots();
bf_selector.z_min = zmin;
bf_selector.z_max = zmax;
bf_selector = bf_selector.updateFTable();
bf_selector.saveMe();

%[thresh_hb, win_out_hb, score_thresh_hb, scanst_hb] = RNA_Threshold_Common.estimateThreshold(spot_counts_ref, [], winsz, 0.5, madf);
%fprintf("Threshold on Ref w/ MADFactor = %.3f, WinSize = %d: %d\n", madf, winsz, thresh_hb);
spotsrun_path = [RefStem '_rnaspotsrun.mat'];
if isfile(spotsrun_path)
    %This will apply the same z trim as has been being used for run
    spotsrun = RNASpotsRun.loadFrom(RefStem);
    thresh_hb = RNAThreshold.runSavedParameters(spotsrun, 1);
else
    thresh_hb = RNAThreshold.runDefaultParameters(spot_counts_ref);
end

% ========================== Export Plots ==========================

plotdir = [InputDir filesep 'plots'];
mkdir(plotdir);

%Spots circled image in BF set at BF and homebrew thresholds
bf_selector = bf_selector.generateColorTables();
bf_selector.toggle_singleSlice = false; %Set to max proj
bf_selector.threshold_idx = bf_thresh - tmin + 1;
bf_selector = bf_selector.drawImages();
saveas(bf_selector.fh_filter, [plotdir filesep 'circlespots_bfthresh.png']);
bf_selector.threshold_idx = thresh_hbbf.threshold - tmin + 1;
bf_selector = bf_selector.drawImages();
saveas(bf_selector.fh_filter, [plotdir filesep 'circlespots_hbbfthresh.png']);
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
color1 = [0.071 0.094 0.529]; %#121887 - Blue
color2 = [0.529 0.051 0.051]; %#870d0d - Red
color3 = [0.369 0.035 0.427]; %#5e096d - Purple
color1_light = [0.627 0.647 0.969]; %#a0a5f7
color2_light = [0.969 0.588 0.588]; %#f79696
color3_light = [0.910 0.588 0.969]; %#e896f7
dim = [.45, .4, .5, .5];
figh = figure(989);
clf;
ax = axes;        
plot(spot_counts(:,1),log10(spot_counts(:,2)),'LineWidth',2,'Color',color1);
hold on;
plot(spot_counts_ref(:,1),log10(spot_counts_ref(:,2)),'-.','LineWidth',2,'Color',color2);
line([bf_thresh bf_thresh], get(ax,'YLim'),'Color',color1_light,'LineStyle','--','LineWidth',2);
line([thresh_hb.threshold thresh_hb.threshold], get(ax,'YLim'),'Color',color2_light,'LineStyle','--','LineWidth',2);
line([thresh_hbbf.threshold thresh_hbbf.threshold], get(ax,'YLim'),'Color',color3_light,'LineStyle','--','LineWidth',2);
legend_names{1,1} = 'BIG-FISH';
legend_names{1,2} = 'Homebrew';
legend(legend_names);
xlabel('Threshold');
ylabel('log10(# Spots Detected)');
saveas(figh, [plotdir filesep 'logspots.png']);
close(figh);
            
%Fscore curves of both spot sets
src_selector = src_selector.updateFTable();

figh = figure(988);
clf;
ax = axes;
plot(spot_counts(:,1),bf_selector.f_scores(:,1),'LineWidth',2,'Color',color1);
hold on;
plot(spot_counts_ref(:,1),src_selector.f_scores(:,1),'-.','LineWidth',2,'Color',color2);
ylim([0.0 1.0]);
line([bf_thresh bf_thresh], get(ax,'YLim'),'Color',color1_light,'LineStyle','--','LineWidth',2);
line([thresh_hb.threshold thresh_hb.threshold], get(ax,'YLim'),'Color',color2_light,'LineStyle','--','LineWidth',2);
line([thresh_hbbf.threshold thresh_hbbf.threshold], get(ax,'YLim'),'Color',color3_light,'LineStyle','--','LineWidth',2);
legend(legend_names);
xlabel('Threshold');
ylabel('FScore');
saveas(figh, [plotdir filesep 'fscoreplot.png']);
close(figh);

%Window score plot
%figh = RNA_Threshold_Common.drawWindowscorePlot(spot_counts(:,1), win_out, score_thresh, thresh_hbbf);
%saveas(figh, [plotdir filesep 'winscore.png']);
%close(figh);

%Cellseg masks visualized
%TODO: Gen an overlay over light channel. Also import Ben's output and gen
%overlay for that. Maybe write a function that does that in core
figh = figure(987);
clf;
% imshow(cell_n,[]);
% saveas(figh, [plotdir filesep 'cell_n.png']);
% clf;
imshow(cell_t,[]);
saveas(figh, [plotdir filesep 'cell_t.png']);
clf;
% imshow(nuc_n,[]);
% saveas(figh, [plotdir filesep 'nuc_n.png']);
% clf;
imshow(nuc_t,[]);
saveas(figh, [plotdir filesep 'nuc_t.png']);
close(figh);

%TODO
%Plot of spot counts overlapping between callsets.
%Also output the overlaps & non-overlap counts at thresholds
%   Could even gen an image with non-overlapped circled for each picked
%   threshold.
