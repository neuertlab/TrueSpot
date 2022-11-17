%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';

addpath('./core');
% ========================== Constants ==========================

DETECT_THREADS = 1;

% ========================== Load csv Table ==========================

AllFigDir = [ImgProcDir filesep 'figures' filesep 'curves'];

InputTablePath = [ImgDir filesep 'test_images.csv'];
imgtbl = readtable(InputTablePath,'Delimiter',',','ReadVariableNames',true,'Format','%s%s%d%d%d%d%s%s%s%s%s%d%s%s%s%d%d%d%d%d%d%d%d%s');

%For one image
SingleImgName = [];
%SingleImgName = 'mESC4d_Tsix-AF594';

% ========================== Iterate through table entries ==========================
entry_count = size(imgtbl,1);

for i = 1:entry_count
    if ~isempty(SingleImgName)
        myname = getTableValue(imgtbl, i, 'IMGNAME');
        if ~strcmp(myname, SingleImgName); continue; end
    end
    
    mystem = replace(getTableValue(imgtbl, i, 'OUTSTEM'), '/', filesep);
    mystem = [ImgDir mystem];
    spotsrun = RNASpotsRun.loadFrom(mystem);
    if isempty(spotsrun)
        fprintf("Image %d of %d - Run could not be found for %s! Skipping...\n", i, entry_count, mystem);
        continue;
    end
    fprintf("Image %d of %d - Thresholding %s...\n", i, entry_count, mystem);
    
    spotsrun.out_stem = mystem;
	if ~isempty(spotsrun.ctrl_stem)
        spotsrun.ctrl_stem = [erase(mystem, 'all_3d') 'Control_all_3d'];
	end
	[spotsrun.out_dir, ~, ~] = fileparts(mystem);
    plotdir = [spotsrun.out_dir filesep 'plots'];

    spotsrun = RNAThreshold.applyPreset(spotsrun, imgtbl{i,'THRESH_SETTING'});
    spotsrun.saveMe();
    
    spotsrun = RNA_Pipeline_Core(spotsrun, 2, [], true, DETECT_THREADS);
    
    %Generate new threshold plots...
    [spotsrun, spotcounts, ~] = spotsrun.loadZTrimmedTables_Sample();
    spotcounts = double(spotcounts);
    spotcounts(:,2)  = log10(spotcounts(:,2));
    
    figh = RNAThreshold.plotThreshRanges(spotsrun,spotcounts,'log10(#Spots)',[],1);
    saveas(figh, [plotdir filesep 'threshold_range_spots.png']);
    saveas(figh, [AllFigDir filesep 'spotcount' filesep spotsrun.img_name '_spotcount.png']);
    close(figh);
    
    if RNA_Threshold_SpotSelector.refsetExists(mystem)
        T = size(spotcounts,1);
        fscores_vals = RNA_Threshold_SpotSelector.loadFScores(mystem);
        fscores = NaN(T,2);
        fscores(:,1) = spotcounts(:,1);
        fscores(:,2) = fscores_vals(:,1);
        figh = RNAThreshold.plotThreshRanges(spotsrun,fscores,'F-Score',[0 1],1);
        saveas(figh, [plotdir filesep 'threshold_range_fscore.png']);
        saveas(figh, [AllFigDir filesep 'fscore' filesep spotsrun.img_name '_fscore.png']);
        close(figh);
    end
    
    %Do BIG-FISH compare, if there is BIG-FISH output...
    bfstem = [ImgDir replace(getTableValue(imgtbl, i, 'BIGFISH_OUTSTEM'), '/', filesep)];
    [bf_dir, bf_imgname, ~] = fileparts(bfstem);
    [success_bool, fighandles] = BigfishCompare.doBigfishCompare(mystem, bf_dir, bf_imgname, true, false);
    if success_bool
        fprintf("Big-FISH run processed. Now saving figures...\n");
        if ~isempty(fighandles.spotsfig)
            saveas(fighandles.spotsfig, [AllFigDir filesep 'bf_spotcount' filesep spotsrun.img_name '_bfspots.png']);
            close(fighandles.spotsfig);
        end
        if ~isempty(fighandles.fscorefig)
            saveas(fighandles.fscorefig, [AllFigDir filesep 'bf_fscore' filesep spotsrun.img_name '_bffscore.png']);
            close(fighandles.fscorefig);
        end
        if ~isempty(fighandles.bfhbimg)
            close(fighandles.bfhbimg);
        end
        
        if ~isempty(fighandles.bfhbbfimg)
            close(fighandles.bfhbbfimg);
        end
    end
    
end

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    val = val{1,1};
end