%
%%
%Dump AUC, selected F-Scores, and peak F-Scores

%!! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

addpath('./core');

% ========================== Other paths ==========================

DataFilePath = [BaseDir filesep 'data' filesep 'DataComposite.mat'];
ImportListPath = [BaseDir filesep 'dumplist.txt'];

OutputTablePath = [BaseDir filesep 'statsdump.csv'];

% ========================== Prepare files ==========================

dumpFiles = readtable(ImportListPath, 'Delimiter', ',', 'ReadVariableNames',false);
load(DataFilePath, 'image_analyses');

OutputTableFile = fopen(OutputTablePath, 'w');
fprintf(OutputTableFile,'ImageName,SPOTS_TS,SPOTS_HB,SPOTS_BF,AUC_HB,AUC_BF,AUC_RS,AUC_DB,Fscore_HB,Fscore_BF\n');

% ========================== Get data ==========================

%This find is SUPER inefficient but idc
output_count = size(dumpFiles,1);
image_count = size(image_analyses,2);
for i = 1:output_count
    %Find its analysis data
    iname = getTableValue(dumpFiles, i, 1);
    afound = false;
    for j = 1:image_count
        if strcmp(iname,image_analyses(j).imgname)
            afound = true;
            myAnalysis = image_analyses(j).analysis;
        end
    end
    
    fprintf(OutputTableFile,'%s', iname);
    if ~afound
        fprintf('Analysis data for requested image "%s" not found!\n', iname);
        for k = 1:9
            fprintf(OutputTableFile,',NaN');
        end
        fprintf(OutputTableFile,'\n');
        continue;
    end
    
    fprintf('Working on %s...\n', iname);
    
    %-------- Spot counts
    %-> Truthset
    if ~isempty(myAnalysis.truthset)
        ts1 = myAnalysis.truthset{1};
        if ~isempty(ts1)
            fprintf(OutputTableFile,',%d', size(ts1,1));
        else
            fprintf(OutputTableFile,',NaN');
        end
    else
        fprintf(OutputTableFile,',NaN');
    end
    
    %-> HB
    if ~isempty(myAnalysis.callset_homebrew)
        callset = myAnalysis.callset_homebrew;
        if ~isempty(callset)
            fprintf(OutputTableFile,',%d', size(callset,1));
        else
            fprintf(OutputTableFile,',NaN');
        end
    else
        fprintf(OutputTableFile,',NaN');
    end
    
    %-> BF
    if ~isempty(myAnalysis.callset_bigfish)
        callset = myAnalysis.callset_bigfish;
        if ~isempty(callset)
            fprintf(OutputTableFile,',%d', size(callset,1));
        else
            fprintf(OutputTableFile,',NaN');
        end
    else
        fprintf(OutputTableFile,',NaN');
    end
    
    %-------- AUCs
    if ~isempty(myAnalysis.res_homebrew)
        restbl = myAnalysis.res_homebrew{1};
        if ~isempty(restbl)
            x = table2array(restbl(:,'sensitivity'));
            y = table2array(restbl(:,'precision'));
            auc_value = RNAUtils.calculateAUC(x, y);
            fprintf(OutputTableFile,',%f',auc_value);
        else
            fprintf(OutputTableFile,',NaN');
        end
    else
        fprintf(OutputTableFile,',NaN');
    end
    
    if ~isempty(myAnalysis.res_bigfish)
        restbl = myAnalysis.res_bigfish{1};
        if ~isempty(restbl)
            x = table2array(restbl(:,'sensitivity'));
            y = table2array(restbl(:,'precision'));
            auc_value = RNAUtils.calculateAUC(x, y);
            fprintf(OutputTableFile,',%f',auc_value);
        else
            fprintf(OutputTableFile,',NaN');
        end
    else
        fprintf(OutputTableFile,',NaN');
    end
    
    if ~isempty(myAnalysis.res_rsfish)
        restbl = myAnalysis.res_rsfish{1};
        if ~isempty(restbl)
            x = table2array(restbl(:,'sensitivity'));
            y = table2array(restbl(:,'precision'));
            auc_value = RNAUtils.calculateAUC(x, y);
            fprintf(OutputTableFile,',%f',auc_value);
        else
            fprintf(OutputTableFile,',NaN');
        end
    else
        fprintf(OutputTableFile,',NaN');
    end
    
    if ~isempty(myAnalysis.res_deepblink)
        restbl = myAnalysis.res_deepblink{1};
        if ~isempty(restbl)
            x = table2array(restbl(:,'sensitivity'));
            y = table2array(restbl(:,'precision'));
            auc_value = RNAUtils.calculateAUC(x, y);
            fprintf(OutputTableFile,',%f',auc_value);
        else
            fprintf(OutputTableFile,',NaN');
        end
    else
        fprintf(OutputTableFile,',NaN');
    end
    
    
    %-------- F-Scores
    if ~isempty(myAnalysis.res_homebrew)
        restbl = myAnalysis.res_homebrew{1};
        if ~isempty(restbl)
            fscore = getTableValue(restbl,myAnalysis.threshold_index_hb,'fScore');
            fprintf(OutputTableFile,',%f',fscore);
        else
            fprintf(OutputTableFile,',NaN');
        end
    else
        fprintf(OutputTableFile,',NaN');
    end
    
    if ~isempty(myAnalysis.res_bigfish)
        restbl = myAnalysis.res_bigfish{1};
        if ~isempty(restbl)
            fscore = getTableValue(restbl,myAnalysis.threshold_index_bf,'fScore');
            fprintf(OutputTableFile,',%f',fscore);
        else
            fprintf(OutputTableFile,',NaN');
        end
    else
        fprintf(OutputTableFile,',NaN');
    end
    
    %-------- Formatting
    fprintf(OutputTableFile,'\n');
end

% ========================== Close ==========================

fclose(OutputTableFile);

% ========================== Helper functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end