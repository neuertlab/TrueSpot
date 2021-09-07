function B_ClusterImageProcessingTMR_AF594_CY5_B_RepA(Task_Number,maxSample,maxIm,Ych,Ypos,globalpath,exppath,osmo,strain,genepair,exp_names,filepath1,...
    positions_touse,thCY7,thAF700,thCY5,thAF594,thTMR,thYFP,thGFP,chi,ch,file_type,Ywin,TMRname,AF594name,CY5name)

% maxFish = 25;
% maxDapi = 25;

maxx = maxSample * maxIm;
if Task_Number > maxx;
    return;
else;
end;


% files = char(strcat(filepath1,'mRNA_TMR-th',num2str(thTMR),'_AF-th',num2str(thAF594),'_CY5-th',num2str(thCY5),'_',exp_names(S1,:),'_im',h)); %strcat(exppath{S1},'_',ch,'_th',ths,'_p',ps,'_im',h)
% if exist(files, 'file')
%     return;
% else;
% end;

%% Determine task number
Task_Number
S1 = floor((Task_Number-1)/(maxIm))+1 % sample
Im1 = Task_Number - ((maxIm)*(S1-1)) %  Image
exp_name = [exppath{S1} '_' exp_names{S1} '_' strain '_' genepair '_img_'];
%                           1           2       3        4        5       6       7      8       9            10        11       12       13      14      15      16
%fish_dir = char(strcat(globalpath, exppath{S1}, '/', exp_name, num2str(S1)))

h=sprintf('%01d',Im1); 

'Check if file exists'
files = char(strcat(filepath1,'mRNA_',exppath{S1},'_STL1-TMR-th',num2str(thTMR),'_GPP2-Af594-th',num2str(thAF594),'_CTT1-GPD1-th',num2str(thCY5),'_',exp_names(:,S1),'_im',h,'.mat'))
%strcat(exppath{S1},'_',ch,'_th',ths,'_p',ps,'_im',h)
if exist(files, 'file') == 0
 
    if strcmp(file_type,'multiple');
    if Ywin == 0;;
        filename = [globalpath exppath{S1} '/' exp_name h '/' exp_name num2str(Im1) '_MMStack.ome.tif']
        else
         filename = [globalpath exppath{S1} '\' exp_name h '\' exp_name num2str(Im1) '_MMStack.ome.tif']
        end;
        [stack, img_read] = tiffread2(filename);
        ImSt = img_read/ch;
        
              
        if false;%Ych(8) == 1;                                                     % CY7
            ij = find(chi==8);
            CY7im = [ij:ch:img_read];
            for i = 1:ImSt;
                CY7_ims(:,:,i) = stack(1,CY7im(i)).data;
            end;
        else
            CY7_ims = NaN;
        end;

        if Ych(8) == 1;                                                     % AF700
            ij = find(chi==8);
            AF700im = [ij:ch:img_read];
            for i = 1:ImSt;
                AF700_ims(:,:,i) = stack(1,AF700im(i)).data;
            end;
        else
            AF700_ims = NaN;
        end;

        if Ych(3) == 1;                                                     % CY5
            ij = find(chi==3);
            CY5im = [ij:ch:img_read];
            for i = 1:ImSt;
                CY5_ims(:,:,i) = stack(1,CY5im(i)).data;
            end;
        else
            CY5_ims = NaN;
        end;

        if Ych(4) == 1;                                                     % AF594
            ij = find(chi==4);
            AF594im = [ij:ch:img_read];
            for i = 1:ImSt;
                AF594_ims(:,:,i) = stack(1,AF594im(i)).data;
            end;
        else
            AF594_ims = NaN;
        end;


        if Ych(5) == 1;                                                     % TMR
            ij = find(chi==5);
            TMRim = [ij:ch:img_read];
            for i = 1:ImSt;
                TMR_ims(:,:,i) = stack(1,TMRim(i)).data;
            end;
        else
            TMR_ims = NaN;
        end;


        if Ych(6) == 1;                                                     % YFP
            ij = find(chi==6);
            YFPim = [ij:ch:img_read];
            for i = 1:ImSt;
                YFP_ims(:,:,i) = stack(1,YFPim(i)).data;
            end;
        else
            YFP_ims = NaN;
        end;

        if Ych(7) == 1;                                                     % GFP / AF488
            ij = find(chi==7);
            GFPim = [ij:ch:img_read];
            for i = 1:ImSt;
                GFP_ims(:,:,i) = stack(1,GFPim(i)).data;
            end;
        else
            GFP_ims = NaN;
        end;
    %         
    %         
    %         if Ych(8) == 1;                                                     % DAPI                                 
    %             ij = find(chi==8);
    %             DAPIim = [ij:ch:img_read];
    %             for i = 1:ImSt;
    %                 DAPI_ims(:,:,i) = stack(1,DAPIim(i)).data;
    %             end;
    %         else
    %         end;
    %         
    %         if Ych(9) == 1;                                                     % TRANS
    %             ij = find(chi==9);
    %             TRANSim = [ij:ch:img_read];
    %             for i = 1:ImSt;
    %                 TRANS_ims(:,:,i) = stack(1,TRANSim(i)).data;
    %             end;
    %         else
    %         end;
    % 
    % 
    elseif strcmp(file_type,'single');
        if false;%Ych(8) == 1;                                                     % CY7
            filebase = 'img_000000000_8-CY7_';
            i = 1;
            for i = 1:50;
                h =sprintf('%03d',i-1); 
                if Ywin == 0;;
                fileN = [globalpath exppath{S1} '/' exp_name num2str(Im1) '/' filebase h '.tif'];
                else
                fileN = [globalpath exppath{S1} '\' exp_name num2str(Im1) '\' filebase h '.tif'];
                end;
                if exist(fileN) == 2;
                    ims = tiffread2(fileN);
                    CY7_ims(:,:,i) = ims.data;
                    ImSt = i;
                else;
                end;
            end;
        else
            CY7_ims = NaN;
        end;

        if Ych(8) == 1;                                                     % AF700
            filebase = 'img_000000000_9-AF700_';
            i = 1;
            for i = 1:50;
                h =sprintf('%03d',i-1); 
                if Ywin == 0;;
                fileN = [globalpath exppath{S1} '/' exp_name num2str(Im1) '/' filebase h '.tif'];
                else
                fileN = [globalpath exppath{S1} '\' exp_name num2str(Im1) '\' filebase h '.tif'];
                end;
                if exist(fileN) == 2;
                    ims = tiffread2(fileN);
                    AF700_ims(:,:,i) = ims.data;
                    ImSt = i;
                else;
                end;
            end;
        else
            AF700_ims = NaN;
        end;

        if Ych(3) == 1;                                                     % CY5
            filebase = 'img_000000000_2-CY5_';
            i = 1;
            for i = 1:50;
                h =sprintf('%03d',i-1); 
                if Ywin == 0;;
                fileN = [globalpath exppath{S1} '/' exp_name num2str(Im1) '/' filebase h '.tif'];
                else
                fileN = [globalpath exppath{S1} '\' exp_name num2str(Im1) '\' filebase h '.tif'];
                end;
                if exist(fileN) == 2;
                    ims = tiffread2(fileN);
                    CY5_ims(:,:,i) = ims.data;
                    ImSt = i;
                else;
                end;
            end;
        else
            CY5_ims = NaN;
        end;

        if Ych(4) == 1;                                                     % AF594
            filebase = 'img_000000000_3-AF594_';
            i = 1;
            for i = 1:50;
                h =sprintf('%03d',i-1); 
                if Ywin == 0;;
                fileN = [globalpath exppath{S1} '/' exp_name num2str(Im1) '/' filebase h '.tif'];
                else
                fileN = [globalpath exppath{S1} '\' exp_name num2str(Im1) '\' filebase h '.tif'];
                end;
                if exist(fileN) == 2;
                    ims = tiffread2(fileN);
                    AF594_ims(:,:,i) = ims.data;
                    ImSt = i;
                else;
                end;
            end;
        else
            AF594_ims = NaN;
        end;

        if Ych(5) == 1;                                                     % TMR
            filebase = 'img_000000000_4-TMR_';
            i = 1;
            for i = 1:50;
                h =sprintf('%03d',i-1); 
                if Ywin == 0;;
                fileN = [globalpath exppath{S1} '/' exp_name num2str(Im1) '/' filebase h '.tif'];
                else
                fileN = [globalpath exppath{S1} '\' exp_name num2str(Im1) '\' filebase h '.tif'];
                end;
                if exist(fileN) == 2;
                    ims = tiffread2(fileN);
                    TMR_ims(:,:,i) = ims.data;
                    ImSt = i;
                else;
                end;
            end;
        else
            TMR_ims = NaN;
        end;

        if Ych(6) == 1;                                                     % YFP
            filebase = 'img_000000000_5-YFP_';
            i = 1;
            for i = 1:50;
                h =sprintf('%03d',i-1); 
                if Ywin == 0;;
                fileN = [globalpath exppath{S1} '/' exp_name num2str(Im1) '/' filebase h '.tif'];
                else
                fileN = [globalpath exppath{S1} '\' exp_name num2str(Im1) '\' filebase h '.tif'];
                end;
                if exist(fileN) == 2;
                    ims = tiffread2(fileN);
                    YFP_ims(:,:,i) = ims.data;
                    ImSt = i;
                else;

                end;
            end;
        else
            YFP_ims = NaN;
        end;

        if Ych(7) == 1;                                                     % GFP / AF488
            filebase = 'img_000000000_6-GFP_';
            i = 1;
            for i = 1:50;
                h =sprintf('%03d',i-1); 
                if Ywin == 0;;
                fileN = [globalpath exppath{S1} '/' exp_name num2str(Im1) '/' filebase h '.tif'];
                else
                fileN = [globalpath exppath{S1} '\' exp_name num2str(Im1) '\' filebase h '.tif'];
                end;
                if exist(fileN) == 2;
                    ims = tiffread2(fileN);
                    GFP_ims(:,:,i) = ims.data;
                    ImSt = i;
                else;
                end;
            end;
        else
            GFP_ims = NaN;
        end;
    else
        'test'
        return;
    end;

    %% check if image is usefull
    if positions_touse{S1,1},Im1 == Im1;
    %% Load segmented cells',ch,'_im',num2str(m),'.mat'
      exp_name_1 = [exppath{S1} '_' exp_names{S1} '_' strain '_' genepair '_img']
      if Ywin == 0;
            fileLAB = char(strcat(filepath1, '/Lab_',exp_name_1,num2str(Im1),'.mat'))
        else;
            fileLAB = char(strcat(filepath1, '\Lab_',exp_name_1,num2str(Im1),'.mat'))
        end;  
        %load(fileLAB,'Lab','Nuc3D'); %'NucInfo','CellInfo',
        load(fileLAB,'cells'); %'NucInfo','CellInfo',
        Lab = cells;
        mm = max(Lab(:)); % maximum number of cells
        NoCell(Im1,S1) = mm;

    %% Load 3D nucleus    
        if Ywin == 0;
            fileNUC = char(strcat(filepath1, '/nuclei_',exp_name_1,num2str(Im1),'.mat'))
        else;
            fileNUC = char(strcat(filepath1, '\nuclei_',exp_name_1,num2str(Im1),'.mat'))
        end;
        %load(fileLAB,'Lab','Nuc3D'); %'NucInfo','CellInfo',
        load(fileNUC,'nuclei','Label_low','Label_mid','Label_hi'); %'NucInfo','CellInfo',
%         Nuc3D = nuclei;

    %% Determine offset for TMR,AF594,CY5 from beads relative to DAPI;
    %     [OFFset,C,TMRmaxFoff,TMR3D3off,TMR3Dorgoff,TMRmaxoff] = B3_ImOFF2D_Bead(TMRmaxF,AF594maxF,CY5maxF,TMR3D3,TMR3Dorg,TMRmax,maxFish);

        %[TMR3D3,TMRmax] = TMRprocess2D(fileTMR,w1,w2,th1,p,con,Y,maxFish);
        %[TMR3D3,TMRmax] = TMRprocess2Dim(fileTMR,w1,w2,th1,p,con,Y,maxFish);
    %     if (Y == 1);
    %         [TMR3D3immax,TMRmax] = TMRprocess3Dim_Test(fileTMR,w1,w2,th1,p,con,Y,maxFish,sample,x,y,Lab,st);
    %         return;
    %     elseif (Y == 0);

    %% Filter images
            % TMRmax:           Max projection original 3D stack
            % TMRmaxF:          Max projection filtered 3D stack
            % thallTMR          for each plane 1: std; 2: mean; 3: variance; 4: meadian
            % TMR3D3immax       Brightest pixel of each spot == maximum spot intensity
            % TMRmaxFimmax      Maximum projection of brightest pixel from each spot
            % TMR3Dorg:         Original 3D stack
            % TMR3filter:       Filtered 3D stack
            % TMR3Dback         Background image == filtered fluorencent image


        %% filter Cy7 images        
            if false;%Ych(8) == 1;
                f = 4;
                'filter Cy7 images'
                [CY7max,CY7maxF,thallCY7,CY73D3immax,CY7maxFimmax,CY73Dorg,CY73Dfilter,CY73Dback]...
                    = B2_RNAprocess3Dim(CY7_ims,thCY7,f);
        %         figure; imshow(CY5maxF,[0 1]);
            else;
                CY7max = NaN; CY7maxF = NaN; thallCY7 = NaN; CY73D3immax = NaN; CY7maxFimmax = NaN; CY73Dorg = NaN; CY73Dfilter = NaN; CY73Dback = NaN;
            end;

            %% filter AF700 images
            if Ych(8) == 1;
                f = 4;
                'filter AF700 images'
                [AF700max,AF700maxF,thallAF700,AF7003D3immax,AF700maxFimmax,AF7003Dorg,AF7003Dfilter,AF7003Dback]...
                = B2_RNAprocess3Dim(AF700_ims,thAF700,f);
        %         figure; imshow(AF700maxF,[0 1]);
            else;
                AF700max = NaN; AF700maxF = NaN; thallAF700 = NaN; AF7003D3immax = NaN; AF700maxFimmax = NaN; AF7003Dorg = NaN; AF7003Dfilter = NaN; AF7003Dback = NaN;
            end;

            %% filter Cy5 images        
            if Ych(3) == 1;
                f = 4;
                'filter Cy5 images'
                [CY5max,CY5maxF,thallCY5,CY53D3immax,CY5maxFimmax,CY53Dorg,CY53Dfilter,CY53Dback]...
                    = B2_RNAprocess3DimCloud(CY5_ims,thCY5,f);
               %  figure; imshow(CY5maxF,[0 1]);
            else;
                CY5max = NaN; CY5maxF = NaN; thallCY5 = NaN; CY53D3immax = NaN; CY5maxFimmax = NaN; CY53Dorg = NaN; CY53Dfilter = NaN; CY53Dback = NaN;
            end;

            %% filter AF594 images
            if Ych(4) == 1;
                f = 4;
                'filter AF594 images'
                [AF594max,AF594maxF,thallAF594,AF5943D3immax,AF594maxFimmax,AF5943Dorg,AF5943Dfilter,AF5943Dback]...
                    = B2_RNAprocess3Dim(AF594_ims,thAF594,f);
        %         figure; imshow(AF594maxF,[0 1]);
            else;
                AF594max = NaN; AF594maxF = NaN; thallAF594 = NaN; AF5943D3immax = NaN; AF594maxFimmax = NaN; AF5943Dorg = NaN; AF5943Dfilter = NaN; AF5943Dback = NaN;
            end;

            %% filter TMR images
            if Ych(5) == 1;
                f = 4;
                'filter TMR images'
                [TMRmax,TMRmaxF,thallTMR,TMR3D3immax,TMRmaxFimmax,TMR3Dorg,TMR3Dfilter,TMR3Dback]...
                    = B2_RNAprocess3DimCloud(TMR_ims,thTMR,f);
               %figure; imshow(TMRmaxF,[0 1]);
            else;
                TMRmax = NaN; TMRmaxF = NaN; thallTMR = NaN; TMR3D3immax = NaN; TMRmaxFimmax = NaN; TMR3Dorg = NaN; TMR3Dfilter = NaN; TMR3Dback = NaN;
            end;

            %% filter YFP images
            if Ych(6) == 1;
                f = 4;
                'filter YFP images'
                [YFPmax,YFPmaxF,thallYFP,YFP3D3immax,YFPmaxFimmax,YFP3Dorg,YFP3Dfilter,YFP3Dback]...
                = B2_RNAprocess3Dim(YFP_ims,thYFP,f);
        %         figure; imshow(YFPmaxF,[0 1]);
            else;
                YFPmax = NaN; YFPmaxF = NaN; thallYFP = NaN; YFP3D3immax = NaN; YFPmaxFimmax = NaN; YFP3Dorg = NaN; YFP3Dfilter = NaN; YFP3Dback = NaN;
            end;

            %% filter GFP images
            if Ych(7) == 1;
                f = 4;
                'filter GFP images'
                [GFPmax,GFPmaxF,thallGFP,GFP3D3immax,GFPmaxFimmax,GFP3Dorg,GFP3Dfilter,GFP3Dback]...
                = B2_RNAprocess3Dim(GFP_ims,thGFP,f);
        %         figure; imshow(GFPmaxF,[0 1]);
            else;
                GFPmax = NaN; GFPmaxF = NaN; thallGFP = NaN; GFP3D3immax = NaN; GFPmaxFimmax = NaN; GFP3Dorg = NaN; GFP3Dfilter = NaN; GFP3Dback = NaN;  
            end;

    %     maxFish = size(TMR3Dorg,3);

    %     %% Determine offset for TMR from images;
    %     if Ytmr == 1 & Ycy5 == 1 & Yoff == 1;
    %         'Determine offset for TMR from images'
    %         [OFFsetTMR,CTMR,TMRmaxFoff,TMR3Dfilteroff,TMR3Dorgoff,TMR3Dbackoff,TMRmaxoff,TMR3D3immaxoff,TMRmaxFimmaxoff]...
    %             = B3_ImOFF2D(TMRmaxF,CY5maxF,TMR3Dfilter,TMR3Dorg,TMR3Dback,TMRmax,maxFish);
    %         % OFFsetTMR         x-y offset
    %         % CTMR              correlation coefficent for determine the x-y OFFset
    %         % TMRmaxFoff        Max projection filtered 3D stack w/ x-y OFFset
    %         % TMR3Dfilteroff     Filtered 3D stack w/ x-y OFFset
    %         % TMR3Dorgoff       Original 3D stack w/ x-y OFFset
    %         % TMR3Dbackoff      Background image w/ x-y OFFset
    %         % TMRmaxoff         Max projection 3D stack w/ x-y OFFset
    %         % TMR3D3immaxoff    Imregional max stack w/ x-y OFFset
    %         % TMRmaxFimmaxoff   Imregional max stack projection w/ x-y OFFset
    %     else;
    %         'NO offset for TMR from images determined'
    %     end;
    % 
    %     %% Determine offset for AF594 from images;
    %     if Yaf594 == 1 & Ycy5 == 1 & Yoff == 1;
    %         'Determine offset for AF594 from images'
    %         [OFFsetAF594,CAF594,AF594maxFoff,AF5943Dfilteroff,AF5943Dorgoff,AF5943Dbackoff,AF594maxoff,AF5943D3immaxoff,AF594maxFimmaxoff]...
    %             = B3_ImOFF2D(AF594maxF,CY5maxF,AF5943Dfilter,AF5943Dorg,AF5943Dback,AF594max,maxFish);
    %     else;
    %         'NO offset for AF594 from images determined'
    %     end;

    %mm = 10;
        %% Determine CY7 mRNA spot positions and intensity
        %% Determine CY7 mRNA spot positions and intensity nongaussian fit
  tic
        if false;%Ych(8) == 1
            'Low threshold nuc'
            Nuc3D = Label_low;
            [CELLmaxRNAcy7_low,RNAposCY7]...
                =B3_nongaussRNApos(Lab,mm,CY73Dorg,CY73Dback,CY73D3immax,CY73Dfilter,Nuc3D);
            'Mid threshold nuc'
            Nuc3D = Label_mid;
            [CELLmaxRNAcy7_mid,RNAposCY7]...
                =B3_nongaussRNApos(Lab,mm,CY73Dorg,CY73Dback,CY73D3immax,CY73Dfilter,Nuc3D);
            'Hi threshold nuc'
            Nuc3D = Label_hi;
            [CELLmaxRNAcy7_hi, RNAposCY7]...
                =B3_nongaussRNApos(Lab,mm,CY73Dorg,CY73Dback,CY73D3immax,CY73Dfilter,Nuc3D);
        else;
            CELLmaxRNAcy7_low = NaN;
            CELLmaxRNAcy7_mid = NaN;
            CELLmaxRNAcy7_hi = NaN;
            RNAposCY7 = NaN;
        end;
        %if ypos == 1 then do gaussian fit
        if false;%Ych(8) == 1 & Ypos == 1
            Nuc3D = Label_mid;
            [PARcy7_mid] = B4_FitGaussians3(Lab, mm,CY73Dorg,CY73Dback,CY73D3immax,CY73Dfilter,RNAposCY7,Nuc3D);
             Nuc3D = Label_low;
            [PARcy7_low] = B5_RNAgaussTHRESH(CY73Dorg,Lab,PARcy7_mid,mm,Nuc3D);
            Nuc3D = Label_hi;
            [PARcy7_hi] = B5_RNAgaussTHRESH(CY73Dorg,Lab,PARcy7_mid,mm,Nuc3D);
        else;
            PARcy7_low = NaN;
            PARcy7_mid = NaN;
            PARcy7_hi = NaN;
        end;   		
        

        %% Determine AF700 mRNA spot positions and intensity nongaussian fit
        if Ych(8) == 1
            'Low threshold nuc'
            Nuc3D = Label_low;
            [CELLmaxRNAaf700_low,RNAposAF700]...
                =B3_nongaussRNApos(Lab,mm,AF7003Dorg,AF7003Dback,AF7003D3immax,AF7003Dfilter,Nuc3D);
            'Mid threshold nuc'
            Nuc3D = Label_mid;
            [CELLmaxRNAaf700_mid,RNAposAF700]...
                =B3_nongaussRNApos(Lab,mm,AF7003Dorg,AF7003Dback,AF7003D3immax,AF7003Dfilter,Nuc3D);
            'Hi threshold nuc'
            Nuc3D = Label_hi;
            [CELLmaxRNAaf700_hi, RNAposAF700]...
                =B3_nongaussRNApos(Lab,mm,AF7003Dorg,AF7003Dback,AF7003D3immax,AF7003Dfilter,Nuc3D);
        else;
            CELLmaxRNAaf700_low = NaN;
            CELLmaxRNAaf700_mid = NaN;
            CELLmaxRNAaf700_hi = NaN;
            RNAposAF700 = NaN;
        end;
        %if ypos == 1 then do gaussian fit
        if Ych(8) == 1 & Ypos == 1
            Nuc3D = Label_mid;
            [PARaf700_mid] = B4_FitGaussians3(Lab, mm,AF7003Dorg,AF7003Dback,AF7003D3immax,AF7003Dfilter,RNAposAF700,Nuc3D);
             Nuc3D = Label_low;
            [PARaf700_low] = B5_RNAgaussTHRESH(AF7003Dorg,Lab,PARaf700_mid,mm,Nuc3D);
            Nuc3D = Label_hi;
            [PARaf700_hi] = B5_RNAgaussTHRESH(AF7003Dorg,Lab,PARaf700_mid,mm,Nuc3D);
        else;
            PARaf700_low = NaN;
            PARaf700_mid = NaN;
            PARaf700_hi = NaN;
        end;   
        

 
		%% Determine CY5 mRNA spot positions and intensity nongaussian fit
        if Ych(3) == 1
            'Low threshold nuc'
            Nuc3D = Label_low;
            [CELLmaxRNAcy5_low,RNAposCY5,tran_cloud_CY5,clouds1_CY5,clouds2_CY5]...
                =B3_nongaussRNAposCloud(Lab,mm,CY53Dorg,CY53Dback,CY53D3immax,CY53Dfilter,Nuc3D);
            'Mid threshold nuc'
            Nuc3D = Label_mid;
            [CELLmaxRNAcy5_mid,RNAposCY5,tran_cloud_CY5,clouds1_CY5,clouds2_CY5]...
                =B3_nongaussRNAposCloud(Lab,mm,CY53Dorg,CY53Dback,CY53D3immax,CY53Dfilter,Nuc3D);
            'Hi threshold nuc'
            Nuc3D = Label_hi;
            [CELLmaxRNAcy5_hi, RNAposCY5,tran_cloud_CY5,clouds1_CY5,clouds2_CY5]...
                =B3_nongaussRNAposCloud(Lab,mm,CY53Dorg,CY53Dback,CY53D3immax,CY53Dfilter,Nuc3D);
        else;
            CELLmaxRNAcy5_low = NaN;
            CELLmaxRNAcy5_mid = NaN;
            CELLmaxRNAcy5_hi = NaN;
            RNAposCY5 = NaN;
        end;
        %if ypos == 1 then do gaussian fit
        if Ych(3) == 1 & Ypos == 1
            Nuc3D = Label_mid;
            [PARcy5_mid] = B4_FitGaussians3Cloud(Lab, mm,CY53Dorg,CY53Dback,CY53D3immax,CY53Dfilter,RNAposCY5,Nuc3D,tran_cloud_CY5,clouds1_CY5,clouds2_CY5);
             Nuc3D = Label_low;
            [PARcy5_low] = B5_RNAgaussTHRESH(CY53Dorg,Lab,PARcy5_mid,mm,Nuc3D);
            Nuc3D = Label_hi;
            [PARcy5_hi] = B5_RNAgaussTHRESH(CY53Dorg,Lab,PARcy5_mid,mm,Nuc3D);
        else;
            PARcy5_low = NaN;
            PARcy5_mid = NaN;
            PARcy5_hi = NaN;
        end;   
		
		
        %% Determine AF594 mRNA spot positions and intensity nongaussian fit
        if Ych(4) == 1
            'Low threshold nuc'
            Nuc3D = Label_low;
            [CELLmaxRNAaf594_low,RNAposAF594]...
                =B3_nongaussRNApos(Lab,mm,AF5943Dorg,AF5943Dback,AF5943D3immax,AF5943Dfilter,Nuc3D);
            'Mid threshold nuc'
            Nuc3D = Label_mid;
            [CELLmaxRNAaf594_mid,RNAposAF594]...
                =B3_nongaussRNApos(Lab,mm,AF5943Dorg,AF5943Dback,AF5943D3immax,AF5943Dfilter,Nuc3D);
            'Hi threshold nuc'
            Nuc3D = Label_hi;
            [CELLmaxRNAaf594_hi, RNAposAF594]...
                =B3_nongaussRNApos(Lab,mm,AF5943Dorg,AF5943Dback,AF5943D3immax,AF5943Dfilter,Nuc3D);
        else;
            CELLmaxRNAaf594_low = NaN;
            CELLmaxRNAaf594_mid = NaN;
            CELLmaxRNAaf594_hi = NaN;
            RNAposAF594 = NaN;
        end;
        %if ypos == 1 then do gaussian fit
        if Ych(4) == 1 & Ypos == 1
            Nuc3D = Label_mid;
            [PARaf594_mid] = B4_FitGaussians3(Lab, mm,AF5943Dorg,AF5943Dback,AF5943D3immax,AF5943Dfilter,RNAposAF594,Nuc3D);
             Nuc3D = Label_low;
            [PARaf594_low] = B5_RNAgaussTHRESH(AF5943Dorg,Lab,PARaf594_mid,mm,Nuc3D);
            Nuc3D = Label_hi;
            [PARaf594_hi] = B5_RNAgaussTHRESH(AF5943Dorg,Lab,PARaf594_mid,mm,Nuc3D);
        else;
            PARaf594_low = NaN;
            PARaf594_mid = NaN;
            PARaf594_hi = NaN;
        end;   
		
		
        %% Determine TMR mRNA spot positions and intensity nongaussian fit
        if Ych(5) == 1
            'Low threshold nuc'
            Nuc3D = Label_low;
            [CELLmaxRNAtmr_low,RNAposTMR,tran_cloud_TMR,clouds1_TMR,clouds2_TMR]...
                =B3_nongaussRNAposCloud(Lab,mm,TMR3Dorg,TMR3Dback,TMR3D3immax,TMR3Dfilter,Nuc3D);
            'Mid threshold nuc'
            Nuc3D = Label_mid;
            [CELLmaxRNAtmr_mid,RNAposTMR,tran_cloud_TMR,clouds1_TMR,clouds2_TMR]...
                =B3_nongaussRNAposCloud(Lab,mm,TMR3Dorg,TMR3Dback,TMR3D3immax,TMR3Dfilter,Nuc3D);
            'Hi threshold nuc'
            Nuc3D = Label_hi;
            [CELLmaxRNAtmr_hi, RNAposTMR,tran_cloud_TMR,clouds1_TMR,clouds2_TMR]...
                =B3_nongaussRNAposCloud(Lab,mm,TMR3Dorg,TMR3Dback,TMR3D3immax,TMR3Dfilter,Nuc3D);
        else;
            CELLmaxRNAtmr_low = NaN;
            CELLmaxRNAtmr_mid = NaN;
            CELLmaxRNAtmr_hi = NaN;
            RNAposTMR = NaN;
        end;
        %if ypos == 1 then do gaussian fit
        if Ych(5) == 1 & Ypos == 1
            Nuc3D = Label_mid;
            [PARtmr_mid] = B4_FitGaussians3Cloud(Lab, mm,TMR3Dorg,TMR3Dback,TMR3D3immax,TMR3Dfilter,RNAposTMR,Nuc3D,tran_cloud_CY5,clouds1_CY5,clouds2_CY5);
             Nuc3D = Label_low;
            [PARtmr_low] = B5_RNAgaussTHRESH(TMR3Dorg,Lab,PARtmr_mid,mm,Nuc3D);
            Nuc3D = Label_hi;
            [PARtmr_hi] = B5_RNAgaussTHRESH(TMR3Dorg,Lab,PARtmr_mid,mm,Nuc3D);
        else;
            PARtmr_low = NaN;
            PARtmr_mid = NaN;
            PARtmr_hi = NaN;
        end;
        

		
        %% Determine YFP mRNA spot positions and intensity nongaussian fit
        if Ych(6) == 1
            'Low threshold nuc'
            Nuc3D = Label_low;
            [CELLmaxRNAyfp_low,RNAposYFP]...
                =B3_nongaussRNApos(Lab,mm,YFP3Dorg,YFP3Dback,YFP3D3immax,YFP3Dfilter,Nuc3D);
            'Mid threshold nuc'
            Nuc3D = Label_mid;
            [CELLmaxRNAyfp_mid,RNAposYFP]...
                =B3_nongaussRNApos(Lab,mm,YFP3Dorg,YFP3Dback,YFP3D3immax,YFP3Dfilter,Nuc3D);
            'Hi threshold nuc'
            Nuc3D = Label_hi;
            [CELLmaxRNAyfp_hi, RNAposYFP]...
                =B3_nongaussRNApos(Lab,mm,YFP3Dorg,YFP3Dback,YFP3D3immax,YFP3Dfilter,Nuc3D);
        else;
            CELLmaxRNAyfp_low = NaN;
            CELLmaxRNAyfp_mid = NaN;
            CELLmaxRNAyfp_hi = NaN;
            RNAposYFP = NaN;
        end;
        %if ypos == 1 then do gaussian fit
        if Ych(6) == 1 & Ypos == 1
            Nuc3D = Label_mid;
            [PARyfp_mid] = B4_FitGaussians3(Lab, mm,YFP3Dorg,YFP3Dback,YFP3D3immax,YFP3Dfilter,RNAposYFP,Nuc3D);
             Nuc3D = Label_low;
            [PARyfp_low] = B5_RNAgaussTHRESH(YFP3Dorg,Lab,PARyfp_mid,mm,Nuc3D);
            Nuc3D = Label_hi;
            [PARyfp_hi] = B5_RNAgaussTHRESH(YFP3Dorg,Lab,PARyfp_mid,mm,Nuc3D);
        else;
            PARyfp_low = NaN;
            PARyfp_mid = NaN;
            PARyfp_hi = NaN;
        end;
        
		
        %% Determine GFP mRNA spot positions and intensity nongaussian fit
        if Ych(7) == 1
            'Low threshold nuc'
            Nuc3D = Label_low;
            [CELLmaxRNAgfp_low,RNAposGFP]...
                =B3_nongaussRNApos(Lab,mm,GFP3Dorg,GFP3Dback,GFP3D3immax,GFP3Dfilter,Nuc3D);
            'Mid threshold nuc'
            Nuc3D = Label_mid;
            [CELLmaxRNAgfp_mid,RNAposGFP]...
                =B3_nongaussRNApos(Lab,mm,GFP3Dorg,GFP3Dback,GFP3D3immax,GFP3Dfilter,Nuc3D);
            'Hi threshold nuc'
            Nuc3D = Label_hi;
            [CELLmaxRNAgfp_hi, RNAposGFP]...
                =B3_nongaussRNApos(Lab,mm,GFP3Dorg,GFP3Dback,GFP3D3immax,GFP3Dfilter,Nuc3D);
        else;
            CELLmaxRNAgfp_low = NaN;
            CELLmaxRNAgfp_mid = NaN;
            CELLmaxRNAgfp_hi = NaN;
            RNAposGFP = NaN;
        end;
        %if ypos == 1 then do gaussian fit
        if Ych(7) == 1 & Ypos == 1
            Nuc3D = Label_mid;
            [PARgfp_mid] = B4_FitGaussians3(Lab, mm,GFP3Dorg,GFP3Dback,GFP3D3immax,GFP3Dfilter,RNAposGFP,Nuc3D);
             Nuc3D = Label_low;
            [PARgfp_low] = B5_RNAgaussTHRESH(GFP3Dorg,Lab,PARgfp_mid,mm,Nuc3D);
            Nuc3D = Label_hi;
            [PARgfp_hi] = B5_RNAgaussTHRESH(GFP3Dorg,Lab,PARgfp_mid,mm,Nuc3D);
        else;
            PARgfp_low = NaN;
            PARgfp_mid = NaN;
            PARgfp_hi = NaN;
        end;
        
 toc
 

    else
    end;

    %% Save data
    'Save Data'
    
    %files = char(strcat(filepath1,exppath{S1},'/SegFiles/','mRNA_',exp_name,num2str(Im1),'.mat')) %strcat(exppath{S1},'_',ch,'_th',ths,'_p',ps,'_im',h)
    files = char(strcat(filepath1,'mRNA_',exppath{S1},'_',strain,'_',TMRname,num2str(thTMR),'_', AF594name,num2str(thAF594),'_',CY5name,num2str(thCY5),'_',exp_names(:,S1),'_im',num2str(Im1),'_','non_stringent_cloud','.mat')); %strcat(exppath{S1},'_',ch,'_th',ths,'_p',ps,'_im',h)
    save(char(files),...
    'PARcy7_low','CELLmaxRNAcy7_low','PARcy7_mid','CELLmaxRNAcy7_mid','PARcy7_hi','CELLmaxRNAcy7_hi',...
    'CY73D3immax','CY73Dfilter','CY7max','CY7maxF','CY7maxFimmax','thallCY7',... 
    'PARaf700_low','CELLmaxRNAaf700_low','PARaf700_mid','CELLmaxRNAaf700_mid','PARaf700_hi','CELLmaxRNAaf700_hi',...
    'AF7003D3immax','AF7003Dfilter','AF700max','AF700maxF','AF700maxFimmax','thallAF700',...
    'PARcy5_low','CELLmaxRNAcy5_low','PARcy5_mid','CELLmaxRNAcy5_mid','PARcy5_hi','CELLmaxRNAcy5_hi',...
    'CY53D3immax','CY53Dfilter','CY5max','CY5maxF','CY5maxFimmax','thallCY5',... 
    'PARaf594_low','CELLmaxRNAaf594_low','PARaf594_mid','CELLmaxRNAaf594_mid','PARaf594_hi','CELLmaxRNAaf594_hi',...
    'AF5943D3immax','AF5943Dfilter','AF594max','AF594maxF','AF594maxFimmax','thallAF594',...
    'PARtmr_low','CELLmaxRNAtmr_low','PARtmr_mid','CELLmaxRNAtmr_mid','PARtmr_hi','CELLmaxRNAtmr_hi',...
    'TMR3D3immax','TMR3Dfilter','TMRmax','TMRmaxF','TMRmaxFimmax','thallTMR',... 
		'PARyfp_low','CELLmaxRNAyfp_low','PARyfp_mid','CELLmaxRNAyfp_mid','PARyfp_hi','CELLmaxRNAyfp_hi',...
		'YFP3D3immax','YFP3Dfilter','YFPmax','YFPmaxF','YFPmaxFimmax','thallYFP',... 
		'PARgfp_low','CELLmaxRNAgfp_low','PARgfp_mid','CELLmaxRNAgfp_mid','PARgfp_hi','CELLmaxRNAgfp_hi',...
		'GFP3D3immax','GFP3Dfilter','GFPmax','GFPmaxF','GFPmaxFimmax','thallGFP','-v7.3');                               
else;
end;
