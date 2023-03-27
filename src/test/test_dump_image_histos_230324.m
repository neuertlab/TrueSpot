%
%%  !! UPDATE TO YOUR BASE DIR
DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

ImgDir = 'C:\Users\hospelb\labdata\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';

% ========================== Constants ==========================

START_INDEX = 39;
END_INDEX = 68;

% ========================== Load csv Table ==========================

AllFigDir = [ImgProcDir filesep 'figures' filesep 'histos'];
OutTablePath = [AllFigDir filesep 'image_stats.csv'];

InputTablePath = [DataDir filesep 'test_images.csv'];
imgtbl = testutil_opentable(InputTablePath);

% ========================== Iterate through entries ==========================
entry_count = size(imgtbl,1);

if START_INDEX < 1; START_INDEX = 1; end
if END_INDEX > entry_count; END_INDEX = entry_count; end

OutTableFile = fopen(OutTablePath, 'w');
fprintf(OutTableFile, 'IMGNAME,MIN_RAW,MAX_RAW,MEAN_RAW,STD_RAW,MED_RAW,MODE_RAW,RAW_75,RAW_80,RAW_85,RAW_90,RAW_95,RAW_99,RAW_99OF99,');
fprintf(OutTableFile, 'MIN_F,MAX_F,MEAN_F,STD_F,MED_F,MODE_F,FILT_75,FILT_80,FILT_85,FILT_90,FILT_95,FILT_99,FILT_99OF99\n');
for r = START_INDEX:END_INDEX
    myname = getTableValue(imgtbl, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);
    fprintf(OutTableFile, '%s', myname);
    
    tifpath_raw = getTableValue(imgtbl, r, 'IMAGEPATH');
    tifpath = [ImgDir replace(tifpath_raw, '/', filesep)];
    if endsWith(tifpath, '.mat')
        [tifdir, tifname, ~] = fileparts(tifpath);
        tifpath = [tifdir filesep 'tif' filesep tifname '.tif'];
    end
    
    if ~isfile(tifpath)
        fprintf('Tif file %s not found. Skipping...', tifpath);
        fprintf(OutTableFile, '\n');
        continue;
    end
    
    ch_total = getTableValue(imgtbl, r, 'CH_TOTAL');
    ch_sample = getTableValue(imgtbl, r, 'CHANNEL');
    
    [channels, ~] = LoadTif(tifpath, ch_total, [ch_sample], 1);
    sample_image = channels{ch_sample,1};
    clear channels;
    
    %Apply filter
    gaussian_rad = 7;
    
    dead_pix_path = 'deadpix.mat';
    RNA_Threshold_Common.saveDeadPixels(sample_image, dead_pix_path, true);
    
    IMG3D = uint16(sample_image);
    clear sample_image;
    IMG3D = RNA_Threshold_Common.cleanDeadPixels(IMG3D, dead_pix_path, true);
    delete(dead_pix_path);
    
    IMG_filtered = RNA_Threshold_Common.applyGaussianFilter(IMG3D, gaussian_rad, 2);
    IMG_filtered = RNA_Threshold_Common.applyEdgeDetectFilter(IMG_filtered);
    IMG_filtered = RNA_Threshold_Common.blackoutBorders(IMG_filtered, gaussian_rad+1, 0);
    IMG_filtered = uint16(IMG_filtered);
    
    %Histograms
    imin = min(IMG3D, [], 'all');
    imax = max(IMG3D, [], 'all');
    [hbins, ~] = imhist(IMG3D, imax+1);
    hbins_x = [0:1:imax];
    hbins_y = log10(hbins);
    iperctl = prctile(IMG3D,[50 75 80 85 90 95 99],'all');
    top1 = find(IMG3D >= iperctl(7));
    iperctl2 = prctile(IMG3D(top1),[50 75 80 85 90 95 99],'all');
    
    figh1 = figure(1);
    plot(hbins_x, hbins_y);
    xlabel('Intensity');
    ylabel('log10(Count)');
    title([myname ' (Raw)']);
    xline(iperctl(7), ':', '99th Percentile', 'LineWidth',1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
    %xline(iperctl2(1), ':', 'Top 50% of 99th Percentile', 'LineWidth',1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
    xline(iperctl2(7), ':', 'Top 1% of 99th Percentile', 'LineWidth',1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
    saveas(figh1, [AllFigDir filesep myname '_raw.png']);
    close(figh1);
    
    fprintf(OutTableFile, ',%d,%d', imin, imax);
    fprintf(OutTableFile, ',%.2f', nanmean(double(IMG3D), 'all'));
    fprintf(OutTableFile, ',%.2f', nanstd(double(IMG3D), 0, 'all'));
    fprintf(OutTableFile, ',%d', nanmedian(IMG3D, 'all'));
    fprintf(OutTableFile, ',%d', mode(IMG3D, 'all'));
    for i = 2:7
        fprintf(OutTableFile, ',%d', iperctl(i));
    end
    fprintf(OutTableFile, ',%d', iperctl2(7));
    
    %Repeat with filtered 
    imax = max(IMG_filtered, [], 'all');
    [hbins, ~] = imhist(IMG_filtered, imax+1);
    hbins_x = [0:1:imax];
    hbins_y = log10(hbins);
    iperctl = prctile(IMG_filtered,[50 75 80 85 90 95 99],'all');
    top1 = find(IMG_filtered >= iperctl(7));
    iperctl2 = prctile(IMG_filtered(top1),[50 75 80 85 90 95 99],'all');
    
    figh2 = figure(2);
    plot(hbins_x, hbins_y);
    xlabel('Intensity');
    ylabel('log10(Count)');
    title([myname ' (Filtered)']);
    xline(iperctl(7), ':', '99th Percentile', 'LineWidth',1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
    %xline(iperctl2(1), ':', 'Top 50% of 99th Percentile', 'LineWidth',1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
    xline(iperctl2(7), ':', 'Top 1% of 99th Percentile', 'LineWidth',1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
    saveas(figh2, [AllFigDir filesep myname '_filtered.png']);
    close(figh2);
    
    fprintf(OutTableFile, ',%d,%d', min(IMG_filtered, [], 'all'), imax);
    fprintf(OutTableFile, ',%.2f', nanmean(double(IMG_filtered), 'all'));
    fprintf(OutTableFile, ',%.2f', nanstd(double(IMG_filtered), 0, 'all'));
    fprintf(OutTableFile, ',%d', nanmedian(IMG_filtered, 'all'));
    fprintf(OutTableFile, ',%d', mode(IMG_filtered, 'all'));
    for i = 2:7
        fprintf(OutTableFile, ',%d', iperctl(i));
    end
    fprintf(OutTableFile, ',%d', iperctl2(7));
    
    fprintf(OutTableFile, '\n');
end
fclose(OutTableFile);

% ========================== Helper functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end
