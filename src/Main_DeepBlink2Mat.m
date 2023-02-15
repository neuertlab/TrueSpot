%
%%

%Ver 23020900 -- Updated so that table is shifted to 1 based indices.
%   Tables written with prior versions should be shifted up 1 before used.

function Main_DeepBlink2Mat(input_file, output_stem)
writerver = 23020900;
writer_ver_str = 'v 23.02.09.0';
fprintf('Main_DeepBlink2Mat (%s)\n', writer_ver_str);

addpath('./core');

[outputdir, ~, ~] = fileparts(output_stem);
mkdir(outputdir);

%See if file exists at input path. If not, just look for the first csv in
%the directory
if ~isfile(input_file)
    fprintf('Requested input file "%s" does not exist.\n', input_file);
    fprintf('Searching for csvs in parent directory...\n');
    [parent_dir,~,~] = fileparts(input_file);
    dir_contents = dir(parent_dir);
    
    content_count = size(dir_contents,1);
    for i = 1:content_count
        if endsWith(dir_contents(i,1).name, '.csv')
            input_file = [parent_dir filesep dir_contents(i,1).name];
            fprintf('Using "%s"...\n', input_file);
            break;
        end
    end
end

import_table = readtable(input_file,'Delimiter',',','ReadVariableNames',true,'Format',...
    '%f%f%f%f');
import_mtx = table2array(import_table);

prob_cutoffs = [0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99];
prob_count = size(prob_cutoffs,2);
coord_table = cell(prob_count,1);

for i = 1:prob_count
    cutoff = prob_cutoffs(i);
    find_res = find(import_mtx(:,3) >= cutoff);
    if ~isempty(find_res)
        found_count = size(find_res,1);
        these_coords = NaN(found_count,3);
        these_coords(:,1:2) = round(import_mtx(find_res,1:2));
        these_coords(:,3) = import_mtx(find_res,4);
        these_coords = these_coords + 1; %Deepblink output is 0 based
        
        coord_table{i,1} = uint16(these_coords);
    else
    end
end

save([output_stem '_coordTable.mat'], 'coord_table', 'writer_ver_str', 'writerver');

end
