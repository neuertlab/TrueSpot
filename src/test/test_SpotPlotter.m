%

%%
%base_dir = 'D:\usr\bghos\labdat\imgproc';
base_dir = 'C:\Users\Blythe\labdata\imgproc';

%---------- Test inputs!
%Tsix 4d
SpotDat_Path = [base_dir '\data\preprocess\feb2018\Tsix_AF594\Tsix\all_3d\Tsix-AF594_IMG1_all_3d'];
Bkg_Path = [base_dir '\data\bkgmask\mESC_4d\Tsix-AF594_img_1_bkgmask.mat'];
Out_Dir = [base_dir '\data\spotplots\mESC_4d\Tsix-AF594_img_1'];
ProbeNames = [{'Tsix-AF594'};
             {'Bkg_CFilter'}; 
             {'SCR-AF594'}; 
              {'No Probe'}];
%Bkg_Direct_Path = [base_dir '\data\preprocess\feb2018\Tsix_AF594\bkg\bkg'];
Bkg_CFilter_Path = [SpotDat_Path '_bkgmasked'];
%load([Bkg_CFilter_Path '_coordTable'], 'coord_table');
SCR_Path = [base_dir '\data\preprocess\feb2018\Tsix_AF594\SCR\all_3d\SCR-AF594_all_3d']; %If applicable
NP_Path = [base_dir '\data\preprocess\feb2018\NoProbe\AF594\all_3d\NoProbe-AF594_all_3d']; %If applicable
ctrl_paths = [{Bkg_CFilter_Path};
             {SCR_Path};
             {NP_Path}];
          
%Xist 4d
%SpotDat_Path = [base_dir '\data\preprocess\feb2018\Xist_CY5\Xist'];
%Bkg_Path = [base_dir '\data\bkgmask\mESC_4d\Xist-CY5_img_1_bkgmask.mat'];
%Out_Dir = [base_dir '\data\spotplots\mESC_4d\Xist-CY5_img_1'];
%ProbeNames = [{'Xist-CY5'};
%              {'Bkg_CFilter'}; 
%              {'SCR-CY5'}; 
%              {'No Probe'}];
%ProbeNames = [{'Xist-CY5'};
%              {'Bkg_CFilter'}];
%Bkg_CFilter_Path = [SpotDat_Path '_bkgmasked'];
%SCR_Path = [base_dir '\data\preprocess\feb2018\Xist_CY5\SCR\all_3d\SCR-CY5_all_3d']; %If applicable
%NP_Path = [base_dir '\data\preprocess\feb2018\NoProbe\CY5\all_3d\NoProbe-CY5_all_3d']; %If applicable
%ctrl_paths = [{Bkg_CFilter_Path};
%              {SCR_Path};
%              {NP_Path}];
%ctrl_paths = [{Bkg_CFilter_Path}];
          

%Yeast E1R2
%SpotDat_Path = [base_dir '\data\preprocess\YeastFISH\E1R2\Ch1\all_3d\E1R2-CH1_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\YeastFISH\E1R2\Ch2\all_3d\E1R2-CH2_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\yeast\Exp1_rep2_10min_im4_bkgmask.mat'];
%Out_Dir = [base_dir '\data\spotplots\YeastFISH\E1R2\Ch1'];
%Out_Dir = [base_dir '\data\spotplots\YeastFISH\E1R2\Ch2'];
%ProbeNames = [{'STL1-TMR'};
%ProbeNames = [{'CTT1-CY5'};
%              {'Bkg_CFilter'}];
%Bkg_CFilter_Path = [SpotDat_Path '_bkgmasked'];
%load([Bkg_CFilter_Path '_coordTable'], 'coord_table');
%ctrl_paths = [{Bkg_CFilter_Path}];

%---------------- Code
addpath('..');  
addpath('./core');  

%Plotter
RNA_Threshold_Plotter.plotPreprocessedData(SpotDat_Path, ctrl_paths, ProbeNames, Out_Dir, true, 1, 69);