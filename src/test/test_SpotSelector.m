%
%TODO: We have a new problem - some spots appear to be in the data but not
%rendering in the selector randomly? At least only on some z slices.
%But spots on ALL z slices should always render unless specified
%otherwise...
%This doesn't appear to be due to the new masking, so why.

%%

%ImgDir = 'D:\usr\bghos\labdat\imgproc'; %chromat
%ImgDir = 'C:\Users\Blythe\labdata\imgproc'; %aelec
ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc'; %workstation

save_stem_rna = [ImgDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d'];
%nucl_seg_path = [ImgDir '\data\cell_seg\mESC_4d\nuclei_20180202_4d_mESC_Tsix-AF594_img_1.mat'];

%save_stem_rna = [ImgDir '\data\preprocess\feb2018\Xist_CY5\Xist'];
%nucl_seg_path = [ImgDir '\data\cell_seg\mESC_4d\nuclei_20180202_4d_mESC_Xist-CY5_img_1.mat'];

%save_stem_rna = [ImgDir '\data\preprocess\feb2019\1Day\Tsix\mESC_1d_Tsix_all_3d'];

%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E1R2\STL1\E1R2-STL1-TMR_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\STL1\E2R2I3-STL1-TMR_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\CTT1\E2R2I3-CTT1-CY5_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\STL1\E2R2I5-STL1-TMR_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\CTT1\E2R2I5-CTT1-CY5_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E1R2\Ch1\all_3d\E1R2-CH1_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E1R2\Ch2\all_3d\E1R2-CH2_all_3d'];
%nucl_seg_path = [ImgDir '\data\cell_seg\yeast\nuclei_Exp1_rep2_10min_im4.mat'];

%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R1\Ch1\all_3d\E2R1-CH1_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R1\Ch2\all_3d\E2R1-CH2_all_3d'];
%nucl_seg_path = [ImgDir '\data\cell_seg\yeast\nuclei_Exp2_rep1_10min_im2.mat'];

%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\Ch1\all_3d\E2R2-IM3-CH1_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\Ch2\all_3d\E2R2-IM3-CH2_all_3d'];
%nucl_seg_path = [ImgDir '\data\cell_seg\yeast\nuclei_Exp2_rep2_10min_im3.mat'];

%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\Ch1\all_3d\E2R2-IM5-CH1_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\Ch2\all_3d\E2R2-IM5-CH2_all_3d'];
%nucl_seg_path = [ImgDir '\data\cell_seg\yeast\nuclei_Exp2_rep2_10min_im5.mat'];

%save_stem_rna = [ImgDir '\data\preprocess\msb2\2M1m_img2\all_3d\Msb2_02M_1m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\2M5m_img3\all_3d\Msb2_02M_5m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\4M1m_img3\all_3d\Msb2_04M_1m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\4M5m_img2\all_3d\Msb2_04M_5m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\backup\2M1m_img2\Msb2_02M_1m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\2M1m_img2\Msb2_02M_1m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\2M5m_img3\Msb2_02M_5m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\4M1m_img3\Msb2_04M_1m_img3_GFP_all_3d'];

%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I4\Ch2\all_3d\Histone_D0_img4_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I4\Ch3\Histone_D0_img4_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I4\Ch4\Histone_D0_img4_ch4_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I6\Ch2\Histone_D0_img6_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I6\Ch3\Histone_D0_img6_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I6\Ch4\Histone_D0_img6_ch4_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I3\Ch3\Histone_D2_img3_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I3\Ch4\Histone_D2_img3_ch4_all_3d'];

addpath('./core');
spotsrun = RNASpotsRun.loadFrom(save_stem_rna);
%th_tbl = spotsrun.t_min:1:spotsrun.t_max;
%th_tbl = transpose(th_tbl);

%selector = RNA_Threshold_SpotSelector;
%selector = selector.initializeNew(save_stem_rna, spotsrun.intensity_threshold - spotsrun.t_min + 1);
%selector.ztrim = spotsrun.ztrim;
selector = RNA_Threshold_SpotSelector.openSelector(save_stem_rna);
%selector = RNA_Threshold_SpotSelector.openSelector(save_stem_rna, true);

%selector.ztrim = spotsrun.ztrim;
% selector.selmcoords = zeros(4,1);
% selector.selmcoords(1,1) = uint16(207);
% selector.selmcoords(2,1) = uint16(719);
% selector.selmcoords(3,1) = uint16(233);
% selector.selmcoords(4,1) = uint16(745);

%selector = selector.launchGUI(); %Select f+/f- spots from auto detect results
selector = selector.launchRefSelectGUI(); %Agnostic selection

%Nucll Mask
%figure(234);
%imshow(nuclei);

%Counts
%[selector, tpos, fpos, fneg, maskedout] = selector.takeCounts([]); %Unmasked

%Counts - trim z
%ztrimtop = 7;
%ztrimbot = 7;

%Counts - 1/4 sampling

%Counts - 1/8 sampling
