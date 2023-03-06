%
%%

%Ver 23020900 -- Updated so that table is shifted to 1 based indices.
%   Tables written with prior versions should be shifted up 1 before used.
%Ver 23022100 -- Exits if it doesn't find a csv
%Ver 23022300 -- Now runs from prob 0.01 to 1.00
%Ver 23022400 -- x and y were switched :)
%Ver 23022700 -- So I DID have x and y switched by accident, but it looks
%   like DeepBlink has been reading the x and y dims the other way round??

function Main_DeepBlink2Mat(input_file, output_stem)
writerver = 23022700;
writer_ver_str = 'v 23.02.27.0';
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

if ~isfile(input_file)
    fprintf('DeepBlink2Mat - No importable files found! Exiting...\n');
    return;
end

import_table = readtable(input_file,'Delimiter',',','ReadVariableNames',true,'Format',...
    '%f%f%f%f');
import_mtx = table2array(import_table);

prob_cutoffs = [0.01:0.01:1.00];
prob_count = size(prob_cutoffs,2);
coord_table = cell(prob_count,1);
spot_table = NaN(prob_count,2);

for i = 1:prob_count
    cutoff = prob_cutoffs(i);
    find_res = find(import_mtx(:,3) >= cutoff);
    if ~isempty(find_res)
        found_count = size(find_res,1);
        these_coords = NaN(found_count,3);
        these_coords(:,1) = round(import_mtx(find_res,1));
        these_coords(:,2) = round(import_mtx(find_res,2));
        these_coords(:,3) = import_mtx(find_res,4);
        these_coords = these_coords + 1; %Deepblink output is 0 based
        
        spot_table(i,1) = prob_cutoffs(i);
        spot_table(i,2) = found_count;
        coord_table{i,1} = uint16(these_coords);
    else
    end
end

ceilinged_bool = false;
save([output_stem '_coordTable.mat'], 'coord_table', 'writer_ver_str', 'writerver','ceilinged_bool');
save([output_stem '_spotTable.mat'], 'spot_table', 'writer_ver_str', 'writerver');

end
