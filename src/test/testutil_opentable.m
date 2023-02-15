%
%%

function image_table = testutil_opentable(tblfile_path)
image_table = readtable(tblfile_path,'Delimiter',',','ReadVariableNames',true,'Format',...
    '%s%s%d%d%d%d%s%s%s%s%s%s%s%d%s%s%s%d%d%d%d%d%d%d%d%d%d%d%d%s');
end