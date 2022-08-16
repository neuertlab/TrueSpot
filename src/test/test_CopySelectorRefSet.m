%
%%

ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

%FromStem = [ImgDir '\data\preprocess\YeastFISH\E1R2\Ch1\E1R2-CH1_all_3d'];
%FromStem = [ImgDir '\data\preprocess\YeastFISH\E1R2\Ch2\E1R2-CH2_all_3d'];
%FromStem = [ImgDir '\data\preprocess\YeastFISH\E2R1\Ch1\E2R1-CH1_all_3d'];
%FromStem = [ImgDir '\data\preprocess\YeastFISH\E2R1\Ch2\E2R1-CH2_all_3d'];
%FromStem = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\Ch1\E2R2-IM3-CH1_all_3d'];
%FromStem = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\Ch2\E2R2-IM3-CH2_all_3d'];
%FromStem = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\Ch1\E2R2-IM5-CH1_all_3d'];
FromStem = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\Ch2\E2R2-IM5-CH2_all_3d'];
%FromStem = [ImgDir '\data\preprocess\msb2\backup\2M1m_img2\Msb2_02M_1m_img2_GFP_all_3d'];
%FromStem = [ImgDir '\data\preprocess\msb2\backup\2M5m_img3\Msb2_02M_5m_img3_GFP_all_3d'];

%ToStem = [ImgDir '\data\preprocess\YeastFISH\E1R2\STL1\E1R2-STL1-TMR_all_3d'];
%ToStem = [ImgDir '\data\preprocess\YeastFISH\E1R2\CTT1\E1R2-CTT1-CY5_all_3d'];
%ToStem = [ImgDir '\data\preprocess\YeastFISH\E2R1\STL1\E2R1-STL1-TMR_all_3d'];
%ToStem = [ImgDir '\data\preprocess\YeastFISH\E2R1\CTT1\E2R1-CTT1-CY5_all_3d'];
%ToStem = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\STL1\E2R2I3-STL1-TMR_all_3d'];
%ToStem = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\CTT1\E2R2I3-CTT1-CY5_all_3d'];
%ToStem = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\STL1\E2R2I5-STL1-TMR_all_3d'];
ToStem = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\CTT1\E2R2I5-CTT1-CY5_all_3d'];
%ToStem = [ImgDir '\data\preprocess\msb2\2M1m_img2\Msb2_02M_1m_img2_GFP_all_3d'];
%ToStem = [ImgDir '\data\preprocess\msb2\2M5m_img3\Msb2_02M_5m_img3_GFP_all_3d'];

addpath('./core');
selector_from = RNA_Threshold_SpotSelector.openSelector(FromStem, true);
selector_to = RNA_Threshold_SpotSelector.openSelector(ToStem, true);
selector_to.ref_coords = selector_from.ref_coords;

%selector_to.selmcoords = zeros(4,1);
%selector_to.selmcoords(1,1) = uint16(30);
%selector_to.selmcoords(2,1) = uint16(1054);
%selector_to.selmcoords(3,1) = uint16(450);
%selector_to.selmcoords(4,1) = uint16(1474);
selector_to.selmcoords = uint16(selector_from.selmcoords);

selector_to = selector_to.saveMe();