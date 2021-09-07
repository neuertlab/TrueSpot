%

%%

base_dir = 'D:\usr\bghos\labdat\imgproc';
%save_dir = [base_dir '\data\preprocess\feb2018\Xist_CY5'];
%save_stem = [save_dir '\Xist'];

save_dir = [base_dir '\data\preprocess\feb2018\Tsix_AF594\Tsix\all_3d'];
save_stem = [save_dir '\Tsix-AF594_IMG1_all_3d'];

imgstructs_suffix = '_imgviewstructs.mat';
istruct_path = [save_stem imgstructs_suffix];

load(istruct_path, 'my_images');
img1 = my_images(1);

figure(1);
fig = imshow(img1.image,[img1.Lmin img1.Lmax]);
%imsave();
saveas(fig,[save_stem '_FilteredImage.png']);
