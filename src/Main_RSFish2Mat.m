%%
%%

function Main_RSFish2Mat(input_dir, output_stem)
writerver = 23021301;
writer_ver_str = 'v 23.02.13.1';
fprintf('Main_RSFish2Mat (%s)\n', writer_ver_str);

addpath('./core');

[outputdir, ~, ~] = fileparts(output_stem);
mkdir(outputdir);

%Determine how many thresholds were scanned
thmax = 0;
i = 1;
tblpath = [input_dir filesep 'thitr_' num2str(i) '.csv'];
while ~isfile(tblpath)
    i = i + 1;
    tblpath = [input_dir filesep 'thitr_' num2str(i) '.csv'];
end
thmin = i;

while isfile(tblpath)
    thmax = i;
    finfo = dir(tblpath);
    if finfo.bytes < 1
        break;
    end
    i = i + 1;
    tblpath = [input_dir filesep 'thitr_' num2str(i) '.csv'];
end

if thmin < 1 | thmax < 1
    fprintf('No files found!\n');
    return;
end

%Allocate spot and coord tables
coord_table = cell(thmax, 1);
spot_table = NaN(thmax,2);
fit_table = cell(thmax, 1);

%Read in
for i = thmin:thmax
    tblpath = [input_dir filesep 'thitr_' num2str(i) '.csv'];
    
    finfo = dir(tblpath);
    if finfo.bytes < 1
        coord_table{i,1} = [];
        break;
    end
    
    c_table = readtable(tblpath,'Delimiter',',','ReadVariableNames',true,'Format',...
    '%f%f%f%f%f%f');
    scount = size(c_table,1);
    spot_table(i,1) = i;
    spot_table(i,2) = scount;
    
    cth_table = NaN(scount, 3);
    fth_table = NaN(scount, 4);
    for r = 1:scount
        %Add to coord table (round)
        cth_table(r,1) = max(1, round(c_table{r,1}));
        cth_table(r,2) = max(1, round(c_table{r,2}));
        cth_table(r,3) = max(1, round(c_table{r,3}));
        
        fth_table(r,1) = c_table{r,1};
        fth_table(r,2) = c_table{r,2};
        fth_table(r,3) = c_table{r,3};
        fth_table(r,4) = c_table{r,6}; %Intensity
    end
    
    fit_table{i,1} = fth_table;
    coord_table{i,1} = uint16(cth_table);
end

%Save
save([output_stem '_coordTable.mat'], 'coord_table', 'writerver', 'writer_ver_str');
save([output_stem '_spotTable.mat'], 'spot_table', 'writerver', 'writer_ver_str');
save([output_stem '_fitTable.mat'], 'fit_table', 'writerver', 'writer_ver_str');

%Delete original csvs
%Comment this out until we know the above works correctly
i = thmin;
tblpath = [input_dir filesep 'thitr_' num2str(i) '.csv'];
while isfile(tblpath)
    delete(tblpath);
    i = i + 1;
    tblpath = [input_dir filesep 'thitr_' num2str(i) '.csv'];
end

end