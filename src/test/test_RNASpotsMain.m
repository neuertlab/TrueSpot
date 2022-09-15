%
%%
addpath('./core');

%ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
ImgDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Image Channels ==========================

%----- mESC Set 1
save_stem_rna = [ImgDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2018\Xist_CY5\Xist\Xist-CY5_IMG1_all_3d'];

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

% ========================== Load Existing RNASpotsRun for some params ==========================

path_spotsrun = [save_stem_rna '_rnaspotsrun.mat'];
if isfile(path_spotsrun)
    spotsrun = RNASpotsRun.loadFrom(path_spotsrun);
end

% ========================== Prepare Args ==========================

tif_path = spotsrun.tif_path;
img_name = spotsrun.img_name;
outdir = spotsrun.out_dir;
cellseg_path = spotsrun.cellseg_path;
ch_sample = num2str(spotsrun.rna_ch);
ch_light = num2str(spotsrun.light_ch);
ch_total = num2str(spotsrun.total_ch);
thmin = num2str(spotsrun.t_min);
thmax = num2str(spotsrun.t_max);
ztrim = num2str(spotsrun.ztrim);
probetype = spotsrun.type_probe;
target_str = spotsrun.type_target;
ttype = spotsrun.type_targetmol;
species = spotsrun.type_species;
celltype = spotsrun.type_cell;

if ~isempty(spotsrun.ctrl_path)
    ctrl_tif = spotsrun.ctrl_path;
    ch_ctrl = num2str(spotsrun.ctrl_ch);
    ch_ctrltot = num2str(spotsrun.ctrl_chcount);
end

tuning_str = '-specific';

% ========================== Call Main ==========================

if ~isempty(spotsrun.ctrl_path)
    if isempty(tuning_str)
        Main_RNASpots('-imgname', img_name, '-tif', tif_path, '-outdir', outdir, '-cellseg', cellseg_path, ...
    '-chsamp', ch_sample, '-chtrans', ch_light, '-chtotal', ch_total, '-ctrltif', ctrl_tif, '-chctrsamp', ch_ctrl, ...
        '-chctrtotal', ch_ctrltot, '-thmin', thmin, '-thmax', thmax, '-ztrim', ztrim, '-probetype', probetype, ...
        '-target', target_str, '-targettype', ttype, '-species', species, '-celltype', celltype, '-debug');
    else
        Main_RNASpots('-imgname', img_name, '-tif', tif_path, '-outdir', outdir, '-cellseg', cellseg_path, ...
        '-chsamp', ch_sample, '-chtrans', ch_light, '-chtotal', ch_total, '-ctrltif', ctrl_tif, '-chctrsamp', ch_ctrl, ...
         '-chctrtotal', ch_ctrltot, '-thmin', thmin, '-thmax', thmax, '-ztrim', ztrim, '-probetype', probetype, ...
         '-target', target_str, '-targettype', ttype, '-species', species, '-celltype', celltype, '-debug', tuning_str);
    end
else
    if isempty(tuning_str)
        Main_RNASpots('-imgname', img_name, '-tif', tif_path, '-outdir', outdir, '-cellseg', cellseg_path, ...
        '-chsamp', ch_sample, '-chtrans', ch_light, '-chtotal', ch_total, ...
        '-thmin', thmin, '-thmax', thmax, '-ztrim', ztrim, '-probetype', probetype, ...
        '-target', target_str, '-targettype', ttype, '-species', species, '-celltype', celltype, '-debug');
    else
        Main_RNASpots('-imgname', img_name, '-tif', tif_path, '-outdir', outdir, '-cellseg', cellseg_path, ...
        '-chsamp', ch_sample, '-chtrans', ch_light, '-chtotal', ch_total, ...
        '-thmin', thmin, '-thmax', thmax, '-ztrim', ztrim, '-probetype', probetype, ...
        '-target', target_str, '-targettype', ttype, '-species', species, '-celltype', celltype, '-debug', tuning_str);
    end
end

