%
%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%BaseDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
addpath('./test');

% ========================== Constants ==========================

START_INDEX = 1;
END_INDEX = 38;

% ========================== Load csv Table ==========================

DataFilePath = [BaseDir filesep 'data' filesep 'DataComposite.mat']; %Use this to get th idxs

InputTablePath = [BaseDir filesep 'test_images.csv'];
imgtbl = testutil_opentable(InputTablePath);

% ========================== Iterate through table entries ==========================
entry_count = size(imgtbl,1);

if START_INDEX < 1; START_INDEX = 1; end
if END_INDEX > entry_count; END_INDEX = entry_count; end

for r = START_INDEX:END_INDEX
    myname = getTableValue(imgtbl, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);
    
    %TODO
    
end

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end


