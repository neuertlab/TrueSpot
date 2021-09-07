%%
%

base_dir = 'C:\Users\Blythe\labdata\imgproc';

data_stem = [base_dir '\data\preprocess\feb2018\Tsix_AF594\Tsix\all_3d\Tsix-AF594_IMG1_all_3d'];
bkg_stem = [data_stem '_bkgmasked'];

load([bkg_stem '_coordTable'], 'coord_table');
bkg_coord_table = coord_table;

load([data_stem '_coordTable'], 'coord_table');
zmin = 1;
zmax = 69;

T = size(coord_table,1);
spot_counts = zeros(T,2);
bkg_spot_counts = zeros(T,2);
spot_counts(:,1) = 1:1:T;
bkg_spot_counts(:,1) = 1:1:T;

for t = 1:T
    ttbl = coord_table{t};
    btbl = bkg_coord_table{t};
    S = size(ttbl,1);
    count = 0;
    for s = 1:S
        z = ttbl(s,3);
        if z >= zmin
            if z <= zmax
                count = count+1;
            end
        end
    end
    spot_counts(t,2) = count;
    
    S = size(btbl,1);
    bcount = 0;
    for s = 1:S
        z = btbl(s,3);
        if z >= zmin
            if z <= zmax
                bcount = bcount+1;
            end
        end
    end
    
    bkg_spot_counts(t,2) = bcount;
end

%deriv_1 = smooth(diff(spot_counts));
%deriv_2 = smooth(diff(deriv_1));

%Try the method...
addpath('..');  
addpath('./core');  
[thresh, ~] = RNA_Threshold_Common.estimateThreshold(spot_counts, bkg_spot_counts, 5, 0.0, 0.9);
fprintf("Selected threshold: %d\n", thresh);

%Find values of interest

%Draw some plots