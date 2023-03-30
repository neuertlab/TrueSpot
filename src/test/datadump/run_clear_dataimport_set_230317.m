%
%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%BaseDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
addpath('./test');

% ========================== Constants ==========================

START_INDEX = 1017;
END_INDEX = 1017;

TOOL_CODE = 'all';

% ========================== Other paths ==========================

DataFilePath = [BaseDir filesep 'data' filesep 'DataComposite.mat'];

% ========================== Load csv Table ==========================

AllFigDir = [ImgProcDir filesep 'figures' filesep 'curves'];

InputTablePath = [BaseDir filesep 'test_images.csv'];
imgtbl = testutil_opentable(InputTablePath);

% ========================== Iterate through table entries ==========================
entry_count = size(imgtbl,1);

if ~isfile(DataFilePath)
    image_analyses(entry_count) = struct('imgname', '', 'analysis', []);
    %image_analyses(entry_count) = ImageResults;
    for r = 1:entry_count
        %Initialize
        image_analyses(r).imgname = getTableValue(imgtbl, r, 'IMGNAME');
        image_analyses(r).analysis = ImageResults.initializeNew();
        image_analyses(r).analysis.image_name = image_analyses(r).imgname;
    end
else
    load(DataFilePath, 'image_analyses');
    %If entry_count is larger, add new entries to end.
    save_count = size(image_analyses,2);
    if entry_count > save_count
        %image_analyses(entry_count) = ImageResults; %Expand.
        image_analyses(entry_count) = struct('imgname', '', 'analysis', []);
        for r = save_count+1:entry_count
            %Initialize
            image_analyses(r).imgname = getTableValue(imgtbl, r, 'IMGNAME');
            image_analyses(r).analysis = ImageResults.initializeNew();
            image_analyses(r).analysis.image_name = image_analyses(r).imgname;
        end
    end
end

if START_INDEX < 1; START_INDEX = 1; end
if END_INDEX > entry_count; END_INDEX = entry_count; end

for r = START_INDEX:END_INDEX
    myname = getTableValue(imgtbl, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);

    image_analyses(r).analysis = image_analyses(r).analysis.clearResultsForTool(TOOL_CODE);
end
save(DataFilePath, 'image_analyses');

% ========================== Helper functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end
