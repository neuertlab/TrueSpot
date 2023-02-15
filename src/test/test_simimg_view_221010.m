%
%%  !! UPDATE TO YOUR BASE DIR
%ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Image Channels ==========================

%----- Test
data_file = [ImgDir '\img\sim\yeast_proteinGFP_100x_1_blur.mat'];

% ========================== Run ==========================

addpath('./core');
load(data_file, 'imgdat');

% Z = size(imgdat, 3);
% Y = size(imgdat, 1);
% X = size(imgdat, 2);
% 
% xy_rad = 7;
% z_rad = 2;
% 
% %Extra Gaussian?
% useimg = ones(Y, X, Z);
% g_imgs = zeros(Y, X, Z, 5);
% g_imgs(:,:,:,1) = imgdat(:,:,:);
% for g = 4:-1:1
%     g_imgs(:,:,:,g+1) = RNA_Threshold_Common.applyGaussianFilter(imgdat, xy_rad, g, z_rad, false);
% end
% ptiles = prctile(imgdat, [20 40 60 80], 'all');
% useimg(find(imgdat < ptiles(2,1) & imgdat >= ptiles(1,1))) = 2;
% useimg(find(imgdat < ptiles(3,1) & imgdat >= ptiles(2,1))) = 3;
% useimg(find(imgdat < ptiles(4,1) & imgdat >= ptiles(3,1))) = 4;
% useimg(find(imgdat >= ptiles(4,1))) = 5;
% 
% imgdat = zeros(Y,X,Z);
% for g = 1:5
%     imgdat = imgdat + immultiply(g_imgs(:,:,:,g), (useimg == g));
% end

% xy_rad = 7;
% z_rad = 2;
% g_amt = 4;
% 
% ptiles_o = prctile(imgdat, [10 20 30 40 50 60 70 80 90], 'all');
% min_o = min(imgdat, [], 'all');
% max_o = max(imgdat, [], 'all');
% ao = double(min_o);
% do = double(max_o - min_o);
% 
% imgdat = RNA_Threshold_Common.applyGaussianFilter(imgdat, xy_rad, g_amt, z_rad, false);
% imgdat = RNAUtils.medianifyBorder(imgdat, [xy_rad xy_rad z_rad]);
% 
% min_f = min(imgdat, [], 'all');
% max_f = max(imgdat, [], 'all');
% af = double(min_f);
% df = double(max_f - min_f);
% props = (double(imgdat) - af) ./ df;
% imgdat = round((props .* do) + ao);
% 
% ptiles_f = prctile(imgdat, [10 20 30 40 50 60 70 80 90], 'all');

fig_handle = MatImages.viewMaxProjection(imgdat, false, 505);
