%%

% Paths
rootdir = 'C:\Users\Blythe\labdata\imgproc';

%tif_path = [rootdir '\img\mESC_4d\20180202_4d_mESC_Tsix-AF594_img_1_MMStack.ome.tif'];
%save_stem = [rootdir '\data\preprocess\feb2018\Tsix_AF594\Tsix\all_3d\Tsix-AF594_IMG1_all_3d'];
%chID = 3;
%chCount = 5;

%tif_path = [rootdir '\img\mESC_4d\20180202_4d_mESC_No-probe_img_2_MMStack.ome.tif'];
%save_stem = [rootdir '\data\preprocess\feb2018\NoProbe\AF594\all_3d\NoProbe-AF594_all_3d'];
%chID = 3;
%chCount = 5;

%tif_path = [rootdir '\img\mESC_4d\20180202_4d_mESC_No-probe_img_2_MMStack.ome.tif'];
%save_stem = [rootdir '\data\preprocess\feb2018\NoProbe\CY5\all_3d\NoProbe-CY5_all_3d'];
%chID = 2;
%chCount = 5;

tif_path = [rootdir '\img\mESC_4d\20180202_4d_mESC_SCR-CY5_img_1_MMStack.ome.tif'];
save_stem = [rootdir '\data\preprocess\feb2018\Xist_CY5\SCR\all_3d\SCR-CY5_all_3d'];
chID = 2;
chCount = 5;

%Load image
[stack, img_read] = tiffread2(tif_path);
Z = img_read/chCount;
Y = size(stack(1,1).data,1);
X = size(stack(1,1).data,2);

img_channel = NaN(Y,X,Z);
idx = chID;
for z = 1:Z
    img_channel(:,:,z) = stack(1,idx).data;
    idx = idx + chCount;
end

%Detect
addpath('./core');
RNA_Threshold_Common.saveDeadPixels(img_channel);

%Build list of dead pixels
load('recurring pixels 3 out of 6', 'recurring_pixels');
rp_count = size(recurring_pixels,1);
tmp_sz = 0;
tmp = zeros(rp_count, 3);
X
Y
Z
dead_mask = true(Y,X);
for i = 1:rp_count
    p = recurring_pixels(i);
    p = p-1;
    
    %Convert to x,y coords
    x = floor(p/Y);
    y = floor(p - (x*Y));
    
    %Skip if on the border.
    if (x == 0) | (y == 0) | (x >= (X-1)) | (y >= (Y-1))
        continue;
    end
    
    %fprintf("Dead pixel found: %d (%d,%d)\n", p, x, y);
    tmp_sz = tmp_sz + 1;
    tmp(tmp_sz, 1) = (x+1);
    tmp(tmp_sz, 2) = (y+1);
    tmp(tmp_sz, 3) = (p+1);
    dead_mask(y+1,x+1) = false;
end
dead_pix = tmp(1:tmp_sz,:);
dp_count = size(dead_pix,1);
%dead_pix

%Clean from coord and spot count tables
load([save_stem '_spotTable2d.mat'], 'spot_table_2D');
load([save_stem '_coordTable2d.mat'], 'coord_table_2D');
T = size(spot_table_2D, 1);
coord_table_new = cell(T,1);
rcount = 0;
for t = 1:T
    ctbl = coord_table_2D{t};
    S = size(ctbl,1);
    tmp = zeros(S,2);
    ctr = 0;
    rcount = 0;
    for s = 1:S
        x = ctbl(s,1);
        y = ctbl(s,2);
        if dead_mask(y,x)
           %Keep
           ctr = ctr+1;
           tmp(ctr,1) = x;
           tmp(ctr,2) = y;
        else
            rcount = rcount+1;
        end
    end
    coord_table_new{t,1} = tmp(1:ctr,:);
    spot_table_2D(t,2) = ctr;
    if rcount > 0
        fprintf("2D spots removed from threshold %d: %d\n", t, rcount);
    end
end
coord_table_2D = coord_table_new;
save([save_stem '_spotTable2d.mat'], 'spot_table_2D');
save([save_stem '_coordTable2d.mat'], 'coord_table_2D');

load([save_stem '_spotTable.mat'], 'spot_table');
load([save_stem '_coordTable.mat'], 'coord_table');
for t = 1:T
    ctbl = coord_table{t};
    S = size(ctbl,1);
    tmp = zeros(S,3);
    ctr = 0;
    rcount = 0;
    for s = 1:S
        x = ctbl(s,1);
        y = ctbl(s,2);
        if dead_mask(y,x)
           %Keep
           ctr = ctr+1;
           tmp(ctr,1) = x;
           tmp(ctr,2) = y;
           tmp(ctr,3) = ctbl(s,3);
        else
            rcount = rcount+1;
        end
    end
    coord_table{t} = tmp(1:ctr,:);
    spot_table(t,2) = ctr;
    if rcount > 0
        fprintf("3D spots removed from threshold %d: %d\n", t, rcount);
    end
end
save([save_stem '_spotTable.mat'], 'spot_table');
save([save_stem '_coordTable.mat'], 'coord_table');

%Clean from spot selector save
selector = RNA_Threshold_SpotSelector.openSelectorSetPaths(save_stem);
img = selector.loaded_ch;
%   Images
for dp = 1:dp_count
    x = dead_pix(dp,1);
    y = dead_pix(dp,2);
    for z = 1:Z
        sum = 0;
        sum = sum + img(y-1,x,z);
        sum = sum + img(y-1,x-1,z);
        sum = sum + img(y-1,x+1,z);
        sum = sum + img(y,x-1,z);
        sum = sum + img(y,x+1,z);
        sum = sum + img(y+1,x-1,z);
        sum = sum + img(y+1,x,z);
        sum = sum + img(y+1,x+1,z);
        avg = sum/8;
        img(y,x,z) = avg;
    end
end
selector.loaded_ch = img;
istructs = selector.img_structs;
img = istructs(1).image;
for dp = 1:dp_count
    x = dead_pix(dp,1);
    y = dead_pix(dp,2);
    sum = 0;
    sum = sum + img(y-1,x);
    sum = sum + img(y-1,x-1);
    sum = sum + img(y-1,x+1);
    sum = sum + img(y,x-1);
    sum = sum + img(y,x+1);
    sum = sum + img(y+1,x-1);
    sum = sum + img(y+1,x);
    sum = sum + img(y+1,x+1);
    avg = sum/8;
    img(y,x) = avg;
end
istructs(1).image = img;
img = istructs(2).image;
for dp = 1:dp_count
    x = dead_pix(dp,1);
    y = dead_pix(dp,2);
    sum = 0;
    sum = sum + img(y-1,x);
    sum = sum + img(y-1,x-1);
    sum = sum + img(y-1,x+1);
    sum = sum + img(y,x-1);
    sum = sum + img(y,x+1);
    sum = sum + img(y+1,x-1);
    sum = sum + img(y+1,x);
    sum = sum + img(y+1,x+1);
    avg = sum/8;
    img(y,x) = avg;
end
istructs(2).image = img;
selector.img_structs = istructs;

%   Spots
for t = 1:T
    %Pos table
    intbl = selector.positives{t};
    if ~isempty(intbl)
        dim1 = size(intbl,1);
        dim2 = size(intbl,2);
        tmp = zeros(dim1,dim2);
        ctr = 0;
        rcount = 0;
        for s = 1:dim1
            x = intbl(s,1);
            y = intbl(s,2);
            if dead_mask(y,x)
                %Keep
                ctr = ctr+1;
                tmp(ctr,:) = intbl(s,:);
            else
                rcount = rcount+1;
            end
        end
        selector.positives{t} = tmp(1:ctr,:);
        if rcount > 0
            fprintf("Positives removed from threshold %d: %d\n", t, rcount);
        end
    end
    
    %Neg table
    intbl = selector.false_negs{t};
    if ~isempty(intbl)
        dim1 = size(intbl,1);
        dim2 = size(intbl,2);
        tmp = zeros(dim1,dim2);
        ctr = 0;
        rcount = 0;
        for s = 1:dim1
            x = intbl(s,1);
            y = intbl(s,2);
            if dead_mask(y,x)
                %Keep
                ctr = ctr+1;
                tmp(ctr,:) = intbl(s,:);
            else
                rcount = rcount+1;
            end
        end
        selector.false_negs{t} = tmp(1:ctr,:);
        if rcount > 0
            fprintf("Negatives removed from threshold %d: %d\n", t, rcount);
        end
    end
end

intbl = selector.ref_coords;
if ~isempty(intbl)
	dim1 = size(intbl,1);
	dim2 = size(intbl,2);
	tmp = zeros(dim1,dim2);
	ctr = 0;
	rcount = 0;
	for s = 1:dim1
        x = intbl(s,1);
        y = intbl(s,2);
        if dead_mask(y,x)
            %Keep
            ctr = ctr+1;
            tmp(ctr,:) = intbl(s,:);
        else
            rcount = rcount+1;
        end
	end
    selector.ref_coords = tmp(1:ctr,:);
	if rcount > 0
        fprintf("Spots removed from ref set: %d\n", rcount);
	end
end
selector.saveMe();
