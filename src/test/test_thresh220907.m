%
%%
addpath('./core');

%ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
ImgDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Image Channels ==========================

%----- mESC Set 1
img_path = [ImgDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1'];
%img_path = [ImgDir '\data\preprocess\feb2018\Xist_CY5\Xist\Xist-CY5_IMG1'];

% ========================== Run ==========================

test_number = 10;

outstem = [img_path '_all_3d'];
ctrlstem = [img_path 'Control_all_3d'];

path_spotsrun = [outstem '_rnaspotsrun.mat'];
spotsrun = RNASpotsRun.loadFrom(path_spotsrun);
spotsrun.out_stem = outstem;
if ~isempty(spotsrun.ctrl_stem)
    spotsrun.ctrl_stem = ctrlstem;
end

spotsrun.ttune_winsz_min = 3;
spotsrun.ttune_winsz_max = 21;
spotsrun.ttune_winsz_incr = 3;

spotsrun.ttune_madf_min = -1.0;
spotsrun.ttune_madf_max = 1.0;
spotsrun.ttune_spline_itr = 3;

spotsrun.ttune_use_rawcurve = false;
spotsrun.ttune_use_diffcurve = false;
spotsrun.ttune_fit_to_log = false;
spotsrun.ttune_reweight_fit = true;
spotsrun.ttune_fit_strat = 1;

for i = test_number:-1:1
    thres_vec(i) = RNAThreshold.runSavedParameters(spotsrun, 1);
end