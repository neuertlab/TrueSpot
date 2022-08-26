%
%%

ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Image Channels ==========================
%De-comment one at a time

%----- mESC Set 1
%save_stem_rna = [ImgDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2018\Xist_CY5\Xist\Xist-CY5_IMG1_all_3d'];

%----- mESC Set 2
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\1Day\Tsix\mESC_1d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\1Day\Xist\mESC_1d_Xist_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\2Day\Tsix\mESC_2d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\2Day\Xist\mESC_2d_Xist_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\3Day\Tsix\mESC_3d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\3Day\Xist\mESC_3d_Xist_all_3d'];

%----- yeast RNA
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E1R2\STL1\E1R2-STL1-TMR_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E1R2\CTT1\E1R2-CTT1-CY5_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R1\STL1\E2R1-STL1-TMR_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R1\CTT1\E2R1-CTT1-CY5_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\STL1\E2R2I3-STL1-TMR_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\CTT1\E2R2I3-CTT1-CY5_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\STL1\E2R2I5-STL1-TMR_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\CTT1\E2R2I5-CTT1-CY5_all_3d'];

%----- yeast proteins
%save_stem_rna = [ImgDir '\data\preprocess\msb2\2M1m_img2\Msb2_02M_1m_img2_GFP_all_3d'];
save_stem_rna = [ImgDir '\data\preprocess\msb2\2M5m_img3\Msb2_02M_5m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\4M1m_img3\Msb2_04M_1m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\4M5m_img2\Msb2_04M_5m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\NoSalt_img1\Msb2_NoSalt_img1_GFP_all_3d'];

%save_stem_rna = [ImgDir '\data\preprocess\opy2\2M1m_img3\Opy2_02M_1m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\2M5m_img2\Opy2_02M_5m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\4M1m_img1\Opy2_04M_1m_img1_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\4M5m_img3\Opy2_04M_5m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\NoSalt_img2\Opy2_NoSalt_img2_GFP_all_3d'];

%----- Xist/Tsix + histones
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I4\Ch2\Histone_D0_img4_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I4\Ch3\Histone_D0_img4_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I4\Ch4\Histone_D0_img4_ch4_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I6\Ch2\Histone_D0_img6_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I6\Ch3\Histone_D0_img6_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I6\Ch4\Histone_D0_img6_ch4_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I3\Ch2\Histone_D2_img3_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I3\Ch3\Histone_D2_img3_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I3\Ch4\Histone_D2_img3_ch4_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I8\Ch2\Histone_D2_img8_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I8\Ch3\Histone_D2_img8_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I8\Ch4\Histone_D2_img8_ch4_all_3d'];

% ========================== Load Run Info ==========================

addpath('./core');
spotsrun = RNASpotsRun.loadFrom(save_stem_rna);

% ========================== Gen Graphs? ==========================

%Load spot count tables
if spotsrun.ztrim > spotsrun.ztrim_auto
    load([spotsrun.out_stem '_ztrim' num2str(spotsrun.ztrim)], 'trimmed_coords');
    [~, spots_sample] = spotsrun.loadSpotsTable();
    T = size(trimmed_coords,1);
    for i = 1:T
        spots_sample(i,2) = size(trimmed_coords{i},1);
    end
    
    load([spotsrun.ctrl_stem '_ztrim' num2str(spotsrun.ztrim)], 'trimmed_coords');
    [~, spots_control] = spotsrun.loadControlSpotsTable();
    T = size(trimmed_coords,1);
    for i = 1:T
        spots_control(i,2) = size(trimmed_coords{i},1);
    end
else
    [~, spots_sample] = spotsrun.loadSpotsTable();
    [~, spots_control] = spotsrun.loadControlSpotsTable();
end

%Try to load fscore table
fscores = [];
fscore_path = [save_stem_rna '_fscores.csv'];
Tmax = 0;
if isfile(fscore_path)
    fscores = csvread(fscore_path);
    Tmax = size(fscores,1);
end

%Determine window sizes to try
win_sizes = [5, 10, 15, 20, 25];
wincount = size(win_sizes,2);
winscore_plots = cell(wincount,1);
picked_threshes = NaN(wincount,2); %Col 1 is intensity thresh, col 2 is winscore thresh

%MAD factors to try
mad_factors = [-2.0, -1.5, -1.0, -0.75,-0.5, -0.25, 0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1.0, 1.5, 2.0];
mfcount = size(mad_factors,2);
th_mtx = NaN(mfcount, wincount);
fs_mtx = NaN(mfcount, wincount);

%T range
range_test = true;
range_starts = [1, 5, 10, 25, 50, 75, 100, 150, 250, 300, 500, 750];
range_sizes = [10, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500, 750, 1000];
count1 = size(range_starts,2);
count2 = size(range_sizes,2);
%th_mtxs = cell(count1, count2);
%fs_mtxs = cell(count1, count2);
th_mtxs = NaN(count1, count2, mfcount, wincount);
fs_mtxs = NaN(count1, count2, mfcount, wincount);
T = size(spots_sample,1);
rst_sub_min = 1;
rst_sub_max = 4;
rsz_sub_min = 7;
rsz_sub_max = count2;
rst_count_sub = rst_sub_max-rst_sub_min + 1;
rsz_count_sub = rsz_sub_max-rsz_sub_min + 1;

%Generate winscore plots.
if range_test
    rt_dir = [ImgDir '\tables\rangetest\'];
    fsmax = NaN(count2, count1);
    fsmean = NaN(count2, count1);
    fsstd = NaN(count2, count1);
    thstd = NaN(count2, count1);
    %fsmax_I = cell(count2, count1);
    for m = 1:count1
        rs = range_starts(1,m);
        for n = 1:count2
            rz = range_sizes(1,n);
            ed = rs + rz;
            if ed > T
                ed = T;
            end
            rz = ed - rs;
            spots_sample_sub = spots_sample(rs:ed,:);
            if ~isempty(spots_control)
                spots_control_sub = spots_control(rs:ed,:);
            else
                spots_control_sub = [];
            end
            fprintf("DEBUG -- rs = %d, ed = %d, m = %d, n = %d\n", rs, ed, m, n);
            
            for i = 1:wincount
                for j = 1:mfcount
                    th_mtx(j,i) = NaN;
                    fs_mtx(j,i) = NaN;
                    test_winsize = win_sizes(1,i);
                    if rz >= test_winsize
                        [thresh, win_scores, score_thresh, ~] = RNA_Threshold_Common.estimateThreshold(spots_sample_sub, spots_control_sub, test_winsize, 0.5, mad_factors(1,j));
                        if mad_factors(1,j) == 0.0
                            picked_threshes(i,1) = thresh;
                            picked_threshes(i,2) = score_thresh;
                            winscore_plots{i,1} = win_scores;
                        end
                        %th_mtx(j,i) = thresh;
                        if ~isempty(fscores) && thresh > 0 && thresh <= Tmax
                            th_mtx(j,i) = thresh;
                            fs_mtx(j,i) = fscores(thresh,1);
                        end
                    end
                end
            end
            %th_mtxs{m,n} = th_mtx;
            %fs_mtxs{m,n} = fs_mtx;
            th_mtxs(m,n,:,:) = th_mtx(:,:);
            fs_mtxs(m,n,:,:) = fs_mtx(:,:);
            fsmax(n,m) = max(fs_mtx,[],'all','omitnan');
            fsmean(n,m) = nanmean(fs_mtx,'all');
            fsstd(n,m) = nanstd(fs_mtx,0, 'all');
            thstd(n,m) = nanstd(th_mtx,0, 'all');
            tbl_out_path = [rt_dir num2str(n) '_' num2str(m) '_th.csv'];
            csvwrite(tbl_out_path, th_mtx);
            tbl_out_path = [rt_dir num2str(n) '_' num2str(m) '_fs.csv'];
            csvwrite(tbl_out_path, fs_mtx);
        end
    end
    
    %Generate plot.
    %X = Win size, Y = MAD factor, Z = thresh, C = fscore
    %Tiled: Y = range start, X = range size
    spi = 1;
    figure(16);
    [X,Y] = meshgrid(win_sizes,mad_factors);

    somevec = ones(1,256);
    somevec(:,:) = 1.0/255.0;
    sclr_r = [255:-1:0].*somevec;
    sclr_g = [0:1:255] .* somevec;
    sclr_b = zeros(1,256);
    sclr_r(1,256) = 0.0;
    sclr_g(1,1) = 0.0;
    clrmap = zeros(256,3);
    clrmap(:,1) = sclr_r(1,:);
    clrmap(:,2) = sclr_g(1,:);
    clrmap(:,3) = sclr_b(1,:);
    
    Z = NaN(mfcount,wincount);
    %C_idxs = zeros(mfcount,wincount);
    %C = NaN(mfcount,wincount,3);
    C = NaN(mfcount,wincount);
    for m = rst_sub_min:rst_sub_max
        for n = rsz_sub_min:rsz_sub_max
            subplot(rst_count_sub, rsz_count_sub, spi);
            spi = spi + 1;
            k = 1;
            Z(:,:) = th_mtxs(m,n,:,:);
%             C_idxs(:,:) = round(fs_mtxs(m,n,:,:).*(255.0)) + 1;
%             C_idxs(isnan(C_idxs)) = 1;
%             for i = 1:wincount
%                 for j = 1:mfcount
%                     C(j,i,1) = sclr_r(1,C_idxs(j,i));
%                     C(j,i,2) = sclr_g(1,C_idxs(j,i));
%                     C(j,i,3) = sclr_b(1,C_idxs(j,i));
%                 end
%             end
            rst = range_starts(1,m);
            red = rst + range_sizes(1,n);
            if red > T
                red = T;
            end
            C(:,:) = fs_mtxs(m,n,:,:);
            surf(X,Y,Z,C,'FaceColor','interp','FaceLighting','gouraud');
            title(['Range ' num2str(rst) ' - ' num2str(red)]);
            %xlabel('Window Size');
            %ylabel('MAD Factor');
            %zlabel('Threshold');
            colormap(clrmap);
            ax = gca;
            ax.CLim = [0,1];
        end
    end
    
    figure(17);
    m = rst_sub_min; n = rsz_sub_max;
    Z(:,:) = th_mtxs(m,n,:,:);
%     C_idxs(:,:) = round(fs_mtxs(m,n,:,:).*(255.0)) + 1;
%     C_idxs(isnan(C_idxs)) = 1;
%     for i = 1:wincount
%         for j = 1:mfcount
%             C(j,i,1) = sclr_r(1,C_idxs(j,i));
%             C(j,i,2) = sclr_g(1,C_idxs(j,i));
%             C(j,i,3) = sclr_b(1,C_idxs(j,i));
%         end
%     end
    rst = range_starts(1,m);
    red = rst + range_sizes(1,n);
    if red > T
        red = T;
    end
    C(:,:) = fs_mtxs(m,n,:,:);
    surf(X,Y,Z,C,'FaceColor','interp','FaceLighting','gouraud');
    title(['Range ' num2str(rst) ' - ' num2str(red)]);
    xlabel('Window Size');
    ylabel('MAD Factor');
    zlabel('Threshold');
    colormap(clrmap);
    ax = gca;
    ax.CLim = [0,1];
    c = colorbar;
    c.Label.String = 'fscore';
    
    
else
    for i = 1:wincount
        for j = 1:mfcount
            [thresh, win_scores, score_thresh, ~] = RNA_Threshold_Common.estimateThreshold(spots_sample, spots_control, win_sizes(1,i), 0.5, mad_factors(1,j));
            if mad_factors(1,j) == 0.0
                picked_threshes(i,1) = thresh;
                picked_threshes(i,2) = score_thresh;
                winscore_plots{i,1} = win_scores;
            end
            th_mtx(j,i) = thresh;
            if ~isempty(fscores) && thresh > 0 && thresh <= Tmax
                fs_mtx(j,i) = fscores(thresh,1);
            end
        end
    end
end

%Take a look at the intensity distribution of the filtered image...
load([save_stem_rna '_prefilteredIMG'], 'img_filter');
max_intensity = max(img_filter, [], 'all', 'omitnan');
fprintf("Maximum intensity value in filtered image: %d\n", max_intensity);

%figure(1);
%histogram(img_filter);
