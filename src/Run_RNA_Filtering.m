function Run_RNA_Filtering(Task_Number)
tic
%% Input
 Task_Number = 1;
maxSample = 8;
maxIm = 11;
Ych = [0 1 1 0 1 0 0 0 1];                                                  % YCY7 DAPI CY5 AF594 TMR YFP GFP/AF488 AF700 TRANS; defines which fluorescent dyes have been imaged; 
Ypos = 1;                                                                   %This is for the Gaussian filtering. Always have = 1
Ywin = 1;

chi = find(Ych == 1);
ch = size(chi,2);
file_type = 'multiple';                                                     % 'single' or 'multiple'

% globalpath = '/Daten/Microscopy/FISH-Folkert/';
% exppath = '110617'; % Day of experiments
% chamber = cellstr(char('HP_noprobe','HP_YPD','HP_YPA','HP_05h','HP_1h','HP_2h','HP_3h','HP_4h','DP_noprobe','DP_YPD','DP_YPA','DP_05h','DP_1h','DP_2h','DP_3h','DP_4h'));
% filepath1 = '/Daten/Projects/Folkert-ncRNA/Analysis/110617/Data/';

if Ywin == 0
globalpath = '/gpfs23/scratch/neuertg/BenK/';
else
  globalpath = 'C:\Users\keslerbk\Desktop\Microscopy\';  
end
% globalpath = 'C:\Users\lig\Music\'; % man directory
exppath = {'20190216','20190217','20190218','20190219','20190220'}; % Day of experiments
strain = 'F1-2-1';
osmo = '0.2MStep';
genepair = 'Xist-CY5-Tsix-TMR';
exp_names =  {'0d' '12hr' '12hr','1d','2d','5d','9d','9d'};
TMRname = 'Tsix-TMR-th';         %This is just for naming purposes of the final file
AF594name = '';    %This is just for naming purposes of the final file
CY5name = 'Xist-CY5-th';        %This is just for naming purposes of the final file

if Ywin == 0        
filepath1 = '/gpfs23/scratch/neuertg/BenK/2018-07-12 mESC Segmentation (20170607 timecourse) Output/'; %specifies where the segment files are
else
filepath1 =    'C:\Users\keslerbk\Dropbox (VU Basic Sciences)\Vanderbilt Computer\2018-07-12 mESC Segmentation (20170607 timecourse) Output\'
end
% filepath1 = 'C:\Users\lig\Music\'; % segmented cells are stored here

positions_touse = { ...
 [1:11];... %  0min 4                                                        %Change each number to the range of image numbers each timepoint has
 [1:11];... %  1min 4
 [1:11];... %  2min 
 [1:11];... %  04min
 [2:5];... %  06min
 [1:4];... %  08min
 [1:5];... %  10min
 [1:5];... %  15min
 [1:5];... %  20min
 [1:4];... %  25min
 [1:4];... %  30min
 [1:4];... %  35min
 [1:4];... %  40min
 [1:5];... %  45min
 [1:4];... %  50min
 [1:4];... %  55min
};
    

%% Constrant threshold
%% Constrant threshold
thCY7all = [125:25:225]
thAF700all = [125:25:225]
thCY5all = [31 42 50 65 75] % [5:5:25];
thAF594all = [34 47 55 63 69]
thTMRall = [98 115 130 145 170] %[50:10:90]
thYFPall = [100:25:200]
thGFPall = [100:25:200]

i = 1;
for i = 1:5;                                                                % Run different thresholds;
    thCY7 = thCY5all(i)
    thAF700 = thAF700all(i)
    thCY5 = thCY5all(i)
    thAF594 = thAF594all(i)
    thTMR = thTMRall(i)
    thYFP = thYFPall(i)
    thGFP = thGFPall(i)

    B_ClusterImageProcessingTMR_AF594_CY5_B(Task_Number,maxSample,maxIm,Ych,Ypos,globalpath,exppath,osmo,strain,genepair,exp_names,filepath1,positions_touse,...
        thCY7,thAF700,thCY5,thAF594,thTMR,thYFP,thGFP,chi,ch,file_type,Ywin,TMRname,AF594name,CY5name);
end;
toc

%%%%%%%%%%%%%%
