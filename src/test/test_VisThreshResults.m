%
%%

addpath('./core');
ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Paths ==========================

%----- mESC Set 1
%save_stem_rna = [ImgDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d'];
save_stem_rna = [ImgDir '\data\preprocess\feb2018\Xist_CY5\Xist\Xist-CY5_IMG1_all_3d'];

%----- mESC Set 2
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\1Day\Tsix\mESC_1d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\1Day\Xist\mESC_1d_Xist_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\2Day\Tsix\mESC_2d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\2Day\Xist\mESC_2d_Xist_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\3Day\Tsix\mESC_3d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\3Day\Xist\mESC_3d_Xist_all_3d'];

%----- mESC Set 3

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
%save_stem_rna = [ImgDir '\data\preprocess\msb2\2M5m_img3\Msb2_02M_5m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\4M1m_img3\Msb2_04M_1m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\4M5m_img2\Msb2_04M_5m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\NoSalt_img1\Msb2_NoSalt_img1_GFP_all_3d'];

%save_stem_rna = [ImgDir '\data\preprocess\opy2\2M1m_img3\Opy2_02M_1m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\2M5m_img2\Opy2_02M_5m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\4M1m_img1\Opy2_04M_1m_img1_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\4M5m_img3\Opy2_04M_5m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\NoSalt_img2\Opy2_NoSalt_img2_GFP_all_3d'];

%----- Xist/Tsix + histones (Nov 2020)
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

%----- Xist/Tsix + histones (Feb 2020)

% ========================== Load Spotsrun & Threshold Results ==========================

tinfo_path = [save_stem_rna '_threshinfo.mat'];
spotsrun = RNASpotsRun.loadFrom(save_stem_rna);

[datadir, ~, ~] = fileparts(save_stem_rna);
outdir = [datadir filesep 'plots'];

if isfile(tinfo_path)
    load([save_stem_rna '_threshinfo.mat'], 'threshold_results');
end

% ========================== Run Updated Thresholder & Merge ==========================

[spotsrun, sample_spot_table] = spotsrun.loadSpotsTable();
if ~isfile(tinfo_path) | ~isfield(threshold_results, 'struct_ver')
    %Re-run.
	parameter_info = RNA_Threshold_Common.genEmptyThresholdParamStruct();
	parameter_info.verbosity = 1;
	parameter_info.test_data = true;
	parameter_info.test_diff = true;
    parameter_info.sample_spot_table = sample_spot_table;
    [spotsrun, parameter_info.control_spot_table] = spotsrun.loadControlSpotsTable();
    threshold_results = RNA_Threshold_Common.estimateThreshold(parameter_info);
    save([save_stem_rna '_threshinfo.mat'], 'threshold_results');
end

% ========================== Load Spot Selector ==========================

selector = RNA_Threshold_SpotSelector.openSelector(save_stem_rna, true);
selector = selector.updateFTable();

% ========================== Setup Render Parameters ==========================

color0 = [0.290, 0.290, 0.290];
color1A = [0.929, 0.529, 0.529];
color1B = [0.431, 0.020, 0.020];
color2A = [0.529, 0.929, 0.529];
color2B = [0.020, 0.431, 0.020];

color3A1 = [0.529, 0.929, 0.929];
color3B1 = [0.020, 0.431, 0.431];
color3A2 = [0.529, 0.729, 0.929];
color3B2 = [0.020, 0.226, 0.431];
color3A3 = [0.529, 0.529, 0.929];
color3B3 = [0.020, 0.020, 0.431];
color3A4 = [0.729, 0.529, 0.929];
color3B4 = [0.226, 0.020, 0.431];
color3A5 = [0.929, 0.529, 0.929];
color3B5 = [0.431, 0.020, 0.431];

color3A = [color3A1 
           color3A2
           color3A3
           color3A4
           color3A5];
color3B = [color3B1 
           color3B2
           color3B3
           color3B4
           color3B5];

% ========================== Render Plots ==========================

figh1 = figure(505);
clf;
ax = axes;
plot(sample_spot_table(:,1),selector.f_scores(:,1),'LineWidth',2,'Color',color0);
hold on;
ylim([0.0 1.0]);

%Overall
linex = threshold_results.threshold;
if linex > 0
    line([linex linex], get(ax,'YLim'),'Color',color0,'LineStyle','--','LineWidth',2);
end
xlabel('Threshold');
ylabel('FScore');
saveas(figh1, [outdir filesep 'thtest_00_overall_fscore.png']);

figh2 = figure(1);
clf;
ax = axes;
plot(sample_spot_table(:,1),log10(sample_spot_table(:,2)),'LineWidth',2,'Color',color0);
hold on;
if linex > 0
    line([linex linex], get(ax,'YLim'),'Color',color0,'LineStyle','--','LineWidth',2);
end
xlabel('Threshold');
ylabel('log(# Spots)');
saveas(figh2, [outdir filesep 'thtest_00_overall_logplot.png']);

%Data
if ~isempty(threshold_results.test_data)
    figure(figh1);
    clf;
    ax = axes;
    plot(sample_spot_table(:,1),selector.f_scores(:,1),'LineWidth',2,'Color',color0);
    hold on;
    ylim([0.0 1.0]);
    
    linex = threshold_results.test_data.spline_knot_x;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color1B,'LineStyle','--','LineWidth',1);
    end
    linex = threshold_results.test_data.medth_min;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color1A,'LineStyle',':','LineWidth',1);
    end
    linex = threshold_results.test_data.medth_max;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color1A,'LineStyle',':','LineWidth',1);
    end
    
    xlabel('Threshold');
    ylabel('FScore');
    saveas(figh1, [outdir filesep 'thtest_01_data_fscore.png']);
    
    figure(figh2);
    clf;
    ax = axes;
    plot(sample_spot_table(:,1),log10(sample_spot_table(:,2)),'LineWidth',2,'Color',color0);
    hold on;
    
    linex = threshold_results.test_data.spline_knot_x;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color1B,'LineStyle','--','LineWidth',1);
    end
    linex = threshold_results.test_data.medth_min;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color1A,'LineStyle',':','LineWidth',1);
    end
    linex = threshold_results.test_data.medth_max;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color1A,'LineStyle',':','LineWidth',1);
    end
    
    xlabel('Threshold');
    ylabel('log(# Spots)');
    saveas(figh2, [outdir filesep 'thtest_01_data_logplot.png']);
end

%Diff
if ~isempty(threshold_results.test_diff)
    figure(figh1);
    clf;
    ax = axes;
    plot(sample_spot_table(:,1),selector.f_scores(:,1),'LineWidth',2,'Color',color0);
    hold on;
    ylim([0.0 1.0]);
    
    linex = threshold_results.test_diff.spline_knot_x;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color2B,'LineStyle','--','LineWidth',1);
    end
    linex = threshold_results.test_diff.medth_min;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color2A,'LineStyle',':','LineWidth',1);
    end
    linex = threshold_results.test_diff.medth_max;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color2A,'LineStyle',':','LineWidth',1);
    end
    xlabel('Threshold');
    ylabel('FScore');
    saveas(figh1, [outdir filesep 'thtest_02_diff_fscore.png']);
    
    txform = diff(sample_spot_table(:,2));
    txform = smooth(txform);
    txform = abs(txform);
    dsz = size(txform,1);
    figure(figh2);
    clf;
    ax = axes;
    plot(sample_spot_table(1:dsz,1),log10(txform(:,1)),'LineWidth',2,'Color',color0);
    hold on;
    
    linex = threshold_results.test_diff.spline_knot_x;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color2B,'LineStyle','--','LineWidth',1);
    end
    linex = threshold_results.test_diff.medth_min;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color2A,'LineStyle',':','LineWidth',1);
    end
    linex = threshold_results.test_diff.medth_max;
    if linex > 0
        line([linex linex], get(ax,'YLim'),'Color',color2A,'LineStyle',':','LineWidth',1);
    end
    xlabel('Threshold');
    ylabel('log(|diff(Spots)|)');
    saveas(figh2, [outdir filesep 'thtest_02_diff_logplot.png']);
end

%Winscores
if ~isempty(threshold_results.test_winsc)
    wincount = size(threshold_results.test_winsc,2);
    if wincount > 5; wincount = 5; end
    for i = 1:wincount
        tinfo = threshold_results.test_winsc(1,i);
        if ~isempty(tinfo)
            figure(figh1);
            clf;
            ax = axes;
            plot(sample_spot_table(:,1),selector.f_scores(:,1),'LineWidth',2,'Color',color0);
            hold on;
            ylim([0.0 1.0]);
            
            winsize = threshold_results.window_sizes(1,i);
            winsize_str = sprintf('%02d',winsize);
            linex = tinfo.spline_knot_x;
            if linex > 0
                line([linex linex], get(ax,'YLim'),'Color',color3B(i,:),'LineStyle','--','LineWidth',1);
            end
            linex = tinfo.medth_min;
            if linex > 0
                line([linex linex], get(ax,'YLim'),'Color',color3A(i,:),'LineStyle',':','LineWidth',1);
            end
            linex = tinfo.medth_max;
            if linex > 0
                line([linex linex], get(ax,'YLim'),'Color',color3A(i,:),'LineStyle',':','LineWidth',1);
            end
            
            xlabel('Threshold');
            ylabel('FScore');
            saveas(figh1, [outdir filesep 'thtest_03_winsz' winsize_str '_fscore.png']);
            
            dsz = size(sample_spot_table,1)-1;
            figure(figh2);
            clf;
            ax = axes;
            plot(sample_spot_table(1:dsz,1),log10(threshold_results.window_scores(:,i)),'LineWidth',2,'Color',color0);
            hold on;
            
            linex = tinfo.spline_knot_x;
            if linex > 0
                line([linex linex], get(ax,'YLim'),'Color',color3B(i,:),'LineStyle','--','LineWidth',1);
            end
            linex = tinfo.medth_min;
            if linex > 0
                line([linex linex], get(ax,'YLim'),'Color',color3A(i,:),'LineStyle',':','LineWidth',1);
            end
            linex = tinfo.medth_max;
            if linex > 0
                line([linex linex], get(ax,'YLim'),'Color',color3A(i,:),'LineStyle',':','LineWidth',1);
            end
            
            xlabel('Threshold');
            ylabel('log(Window Score)');
            saveas(figh2, [outdir filesep 'thtest_03_winsz' winsize_str '_logplot.png']);
        end
    end
end

close(figh1);
close(figh2);

