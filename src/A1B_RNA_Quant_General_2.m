function [CY5_AF594_TMR,CY5_AF594_TMR_tp,...
    CY5_in_CY5cloud,CY5_in_CY5cloud_tp,CY5_in_TMRcloud,CY5_in_TMRcloud_tp,TMR_in_CY5cloud,TMR_in_CY5cloud_tp,TMR_in_TMRcloud,TMR_in_TMRcloud_tp,...
    CY5_pos_init,CY5_pos_init_tp,TMR_pos_init,TMR_pos_init_tp,T_C_col,C_T_col,T_C_col_tp,C_T_col_tp,T_T_col,...
    T_T_col_tp,CY5_int,CY5_int_tp,TMR_int,TMR_int_tp,C_T_col_ind,C_T_col_ind_tp,T_C_col_ind,T_C_col_ind_tp,...
    CY5_Nuc_dist,CY5_Nuc_dist_tp,TMR_Nuc_dist,TMR_Nuc_dist_tp,Cl_Nuc_dist,Cl_Nuc_dist_tp,...
    C1_CY5_dist,C1_CY5_dist_tp,C2_CY5_dist,C2_CY5_dist_tp,C1_TMR_dist,C1_TMR_dist_tp,C2_TMR_dist,C2_TMR_dist_tp]...
    = A1B_RNA_Quant_General(...
    strain,dates0,times0,times1,max_tps,max_imnum,...
    thCY5all,thAF594all,thTMRall,CY5_name,AF594_name,TMR_name,thNum_CY5,thNum_AF594,thNum_TMR,...
    CY5_hascloud,AF594_hascloud,TMR_hascloud,CY5_cloudInt,AF594_cloudInt,TMR_cloudInt,nucTH,dist_col,CY5_tran_thres,TMR_tran_thres...
    );
%% The code below takes the mRNA files from the RNA analysis and processes them, primarily making the matrices: CY5_AF594_TMR and CY5_AF594_TMR_tp, which is the same but separated by timepoint.
% CY5_AF594_TMR_tp is a cell array with the cell array entry being the
% timepoint (example: CY5_AF594_TMR_tp{1} is the first timepoint).
% Otherwise, it is in the same format as CY5_AF594_TMR.
% CY5_AF594_TMR rows correspond to a single cell. Each column corresponds
% to:
%1: Number of transcripts from CY5
%2: Number of transcripts in AF594
%3: Number of transcripts in TMR
%4: timepoint
%5: date of experiment
%6: image number
%7: Cell number
%8: Number of CY5 transcription sites
%9: Number of AF594 transcription sites
%10: Number of TMR transcription sites
%11: Number of transcripts at the first transcription site/cloud CY5
%12: Number of transcripts at the second transcription site/cloud CY5
%13: Number of transcripts at the first transcription site/cloud TMR
%14: Number of transcripts at the second transcription site/cloud TMR
%15: Used for bins in boxplots (on A0 Analysis page)
%16: Volume of cloud1 CY5
%17: Volume of cloud2 CY5
%18: Volume of cloud1 TMR
%19: Volume of cloud2 TMR
%20: X position of Transcription Site 1 CY5
%21: Y position of Transcription Site 1 CY5
%22: Z position of Transcription Site 1 CY5
%23: X position of Transcription Site 2 CY5
%24: Y position of Transcription Site 2 CY5
%25: Z position of Transcription Site 2 CY5
%26: X position of Transcription Site 3 CY5
%27: Y position of Transcription Site 3 CY5
%28: Z position of Transcription Site 3 CY5
%29: X position of Transcription Site 4 CY5
%30: Y position of Transcription Site 4 CY5
%31: Z position of Transcription Site 4 CY5
%32: X position of Transcription Site 1 TMR
%33: Y position of Transcription Site 1 TMR
%34: Z position of Transcription Site 1 TMR
%35: X position of Transcription Site 2 TMR
%36: Y position of Transcription Site 2 TMR
%37: Z position of Transcription Site 2 TMR
%38: X position of Transcription Site 3 TMR
%39: Y position of Transcription Site 3 TMR
%40: Z position of Transcription Site 3 TMR
%41: X position of Transcription Site 4 TMR
%42: Y position of Transcription Site 4 TMR
%43: Z position of Transcription Site 4 TMR
%CY5_cloud CY5_cloud_tp TMR_cloud TMR_cloud_tp
    %These denote whether the indicated transcript was within the cloud or not
    %0: not a transcript
    %1: a transcript not part of either cloud
    %2: a transcript part of the first cloud
    %3: a transcript part of the second cloud
    %4: a transcript part of both clouds
%CY5_pos_init,CY5_pos_init_tp,TMR_pos_init,TMR_pos_init_tp 
    %are the initially determined positions
    %z=1: X Position of transcript (init)
    %z=2: Y Position of transcript (init)
    %z=3: Z Position of transcript (init)
%CY5_Nuc_dist = Distance of each CY5 transcript from the nucleus 
%CY5_Nuc_dist_tp = Distance of each CY5 transcript from the nucleus 
%TMR_Nuc_dist = Distance of each TMR transcript from the nucleus 
%TMR_Nuc_dist_tp = Distance of each TMR transcript from the nucleus 
%Cl_Nuc_dist = Distance of each cloud from the nucleus (first is cloud1,next is cloud2)
%Cl_Nuc_dist_tp = Distance of each cloud from the nucleus (first is cloud1,next is cloud2)
%CY5_C1_dist = Distance of each CY5 transcript from the first cloud (arranged like CY5_AF594_TMR)
%CY5_C1_dist_tp = Distance of each CY5 transcript from the first cloud (arranged like CY5_AF594_TMR_tp)
%CY5_C2_dist = Distance of each CY5 transcript from the 2nd cloud (arranged like CY5_AF594_TMR)
%CY5_C2_dist_tp = Distance of each CY5 transcript from the 2nd cloud (arranged like CY5_AF594_TMR_tp)
%TMR_C1_dist = Distance of each TMR transcript from the first cloud (arranged like CY5_AF594_TMR)
%TMR_C1_dist_tp = Distance of each TMR transcript from the first cloud (arranged like CY5_AF594_TMR_tp)
%TMR_C2_dist = Distance of each TMR transcript from the 2nd cloud (arranged like CY5_AF594_TMR)
%TMR_C2_dist_tp = Distance of each TMR transcript from the 2nd cloud (arranged like CY5_AF594_TMR_tp)

%% Sample variables
%%%example name:_STL1-TMR-th
% CY5_name = {'_Xist-CY5-th' '_Xist-CY5-th' '_Xist-CY5-th'};
% AF594_name = {'_Tsix-AF594-th' '_Tsix-AF594-th' '_Tsix-AF594-th'};
% TMR_name = {'_Jpx-TMR-th' '_Jpx-TMR-th' '_Jpx-TMR-th'};
% %%% copy over thresholds from experiment
% % thCY5all = [9 11 14 19 23] ;
% % thAF594all = [11 13 16 19 22];
% % thTMRall = [19 22 25 30 35] ;
% thCY5all = {[9 11 15 17 19],[9 11 14 19 23], [9 11 15 17 19]} ;
% thAF594all = {[11 13 16 18 20],[11 13 16 19 22],[10 11 13 16 18]};
% thTMRall = {[19 22 25 28 31] ,[19 22 25 30 35],[19 22 24 26 29]} ;
% %%% which threshold sets
% thNum_CY5 = 3; %Which set of the above thresholds to use (which column)
% thNum_TMR = 4; %Which set of the above thresholds to use (which column)
% thNum_AF594 = 2; %Which set of the above thresholds to use (which column)
% %%%
% dates0 = {'20161109','20161113','20161122'};
% times0 = {{'1hr','48hr','72hr'},{'96hr','120hr','168hr'},{'11d','13d'}}; %168 hour should be somewhere sometime
% times1 = [0,2,3,4,5,6,11,13];
% % dates0 = {'20161113'};
% % times0 = {{'96hr','120hr','168hr'}}; %168 hour should be somewhere sometime
% % times1 = [4,5,6];
% max_imnum = 6;
% max_tps = 3;
% %%% Whether to do Cloud quant or not
% CY5_hascloud = 1;
% AF594_hascloud = 0;
% TMR_hascloud = 0;
% %%%
% CY5_cloudInt = [1887,1887,1887];
% AF594_cloudInt = [];
% TMR_cloudInt = [];
% %%%
% nucTH = 'low'
%% Do analysis
perform_cell_functions = 0; %This is for determining if something is in a cloud or not, and must load the cell boundaries (which take time)
perform_nuc_functions = 0; %This is for determining the distance to the nucleus, and it must load the nucleus info, which takes time to load
if perform_cell_functions == 0
    perform_nuc_functions = 0;
end
perform_col_functions = 1;  %This is for whether to determine colocalization
cloud_correction = 5;
same_dist = 7;   %How close transcripts have to be before they are considered the same transcript
% dist_col = 5;
% CY5_tran_thres = 3; %How many transcripts needed to be determined as a transcription site or cloud
% TMR_tran_thres = 2; %How many transcripts needed to be determined as a transcription site or cloud
% AF594_cloudInt
z_adjust = 3;   %This adjusts the multiplier applied to the z axis to account for the fact that distance between slices will be different that distances between neighboring pixels
clear CY5_AF594_TMR_tp T_C_col_tp C_T_col_tp CY5_int_tp TMR_int_tp T_C_col_ind_tp C_T_col_ind_tp CY5_pos_init_tp TMR_pos_init_tp CY5_in_CY5cloud_tp TMR_in_CY5cloud_tp CY5_in_TMRcloud_tp TMR_in_TMRcloud_tp T_T_col_tp
T_T_col_tp = {};
T_T_col = zeros(2,3);       %shows if there is a colocalized spot in the CY5 channel for every spot in the CY5 channel (control) (0 if not a transcript, 1 if a transcript but no colocalized spot, 2 if there is a colocalized spot)
T_C_col = zeros(2,3);       %shows if there is a colocalized spot in the CY5 channel for every spot in the TMR channel (0 if not a transcript, 1 if a transcript but no colocalized spot, 2 if there is a colocalized spot)
C_T_col = zeros(2,3);       %shows if there is a colocalized spot in the CY5 channel for every spot in the TMR channel (0 if not a transcript, 1 if a transcript but no colocalized spot, 2 if there is a colocalized spot)
C_T_col_ind = zeros(2,3);   %stores the index of the colocalized TMR spot for each CY5 spot (0 if not a transcript, 1 if a transcript but no colocalized spot, 2 if there is a colocalized spot)
C_T_col_ind_tp = {};        %stores the index of the colocalized TMR spot for each CY5 spot  for each timepoint (0 if not a transcript, 1 if a transcript but no colocalized spot, 2 if there is a colocalized spot)
T_C_col_ind = zeros(2,3);   %stores the index of the colocalized CY5 spot for each TMR spot (0 if not a transcript, 1 if a transcript but no colocalized spot, 2 if there is a colocalized spot)
T_C_col_ind_tp = {};        %stores the index of the colocalized CY5 spot for each TMR spot for each timepoint (0 if not a transcript, 1 if a transcript but no colocalized spot, 2 if there is a colocalized spot)
CY5_AF594_TMR = zeros(2,43);         %matrix. Each row corresponds to a single cell, first column is how many CY5 (Xist), second column is how many AF594 (Tsix), third column is how many TMR (Jpx)
T_C_col_tp = {}; %(0 if not a transcript, 1 if a transcript but no colocalized spot, 2 if there is a colocalized spot)
C_T_col_tp = {}; %(0 if not a transcript, 1 if a transcript but no colocalized spot, 2 if there is a colocalized spot)
CY5_int = zeros(2,3);  %Has the fluorescence of every RNA spot (does not include clouds)
CY5_int_tp = {};  %Has the fluorescence of every RNA spot (does not include clouds)
TMR_int = zeros(2,3);  %Has the fluorescence of every RNA spot (does not include clouds)
TMR_int_tp = {};  %Has the fluorescence of every RNA spot (does not include clouds)
CY5_AF594_TMR_tp = {};
CY5_in_CY5cloud = zeros(2,3); %indicates whether CY5 transcripts exist inside the CY5 cloud or not
TMR_in_CY5cloud = zeros(2,3); %indicates whether TMR transcripts exist inside the CY5 cloud or not
CY5_in_CY5cloud_tp = {};  %indicates whether CY5 transcripts exist inside the CY5 cloud or not
TMR_in_CY5cloud_tp = {}; %indicates whether TMR transcripts exist inside the CY5 cloud or not
CY5_in_TMRcloud = zeros(2,3); %indicates whether CY5 transcripts exist inside the CY5 cloud or not
TMR_in_TMRcloud = zeros(2,3); %indicates whether TMR transcripts exist inside the CY5 cloud or not
CY5_in_TMRcloud_tp = {};  %indicates whether CY5 transcripts exist inside the CY5 cloud or not
TMR_in_TMRcloud_tp = {}; %indicates whether TMR transcripts exist inside the CY5 cloud or not
CY5_Nuc_dist = zeros(2,3); CY5_Nuc_dist_tp = {};    %Distance from CY5 spots to the nucleus
TMR_Nuc_dist = zeros(2,3); TMR_Nuc_dist_tp = {};    %Distance from TMR spots to the nucleus
Cl_Nuc_dist = zeros(2,3); Cl_Nuc_dist_tp = {};      %Distance from clouds to the nucleus
C1_CY5_dist = zeros(2,3); C1_CY5_dist_tp = {};      %Distance from cloud 1 to CY5 spots
C2_CY5_dist = zeros(2,3); C2_CY5_dist_tp = {};      %Distance from cloud 2 to CY5 spots
C1_TMR_dist = zeros(2,3); C1_TMR_dist_tp = {};      %Distance from cloud 1 to TMR spots
C2_TMR_dist = zeros(2,3); C2_TMR_dist_tp = {};      %Distance from cloud 1 to TMR spots
if isempty(CY5_name)
    thNum_CY5 = [];
end
if isempty(AF594_name)
    %thNum_AF594 = [];
end
if isempty(TMR_name)
    thNum_TMR = [];
end
for i = 1: size(unique(times1),2)
    CY5_AF594_TMR_tp{i} = zeros(1,43); %Same as CY5_AF594_TMR, but separated by timepoint (first dimension of cell is timepoint)
end
counter9 = 1;   %counter for cells in CY5_AF594_TMR (all timepoints)
counter11 = 1; %counter for which timepoint (third dimension of CY5_AF594_TMR_tp)
fail_count = 0; k = 0; temp_tpname = 'not_a_tp';
for i = 1:size(dates0,2)
%     i
    for j = 1:max_tps
%         j
        if fail_count <k & not(strcmp(temp_tpname,times0{i}{j}));
            counter11 = counter11+1;
        end
        fail_count = 0;   %Skips timepoint if all images fail
        if not(strcmp(temp_tpname,times0{i}{j}))
        counter10 = 1; %counter for CY5_AF594_TMR_tp (cells in a single timepoint)
        end
        for k = 1:max_imnum
             k
             try
                counter9_int = counter9;
                counter10_int = counter10;
                %mRNA_20190117_FemaleTsixInd_Tsix-TMR-th98_34_Xist-CY5-th31_WithDox_im7_stringent_cloud
                if perform_cell_functions
                %%%Load segmentation information (for using ccoordinates of
                %%%cloud)                
                %Lab_20180619_0d_F1-2-1_Xist-CY5-Maoa-TMR_img1
%                 fileLAB = char(['Lab_' dates0{i} '_' times0{i}{j} '_' strain  CY5_name{i}(1:size(CY5_name{i},2)-3) '-' TMR_name{i}(2:size(CY5_name{i},2)-3) '_img' num2str(k) '.mat']); % For regulat
                fileLAB = char(['Lab_' dates0{i} '_'  strain times0{i}{j} CY5_name{i}(1:size(CY5_name{i},2)-3) '-' TMR_name{i}(2:size(CY5_name{i},2)-3) '_img' num2str(k) '.mat']); %For Tsix ind experiments
                load(fileLAB,'cells'); %'NucInfo','CellInfo',
                Lab = cells;
                end
                if perform_nuc_functions
                fileNUC = char(['nuclei_' dates0{i} '_' times0{i}{j} '_' strain  CY5_name{i}(1:size(CY5_name{i},2)-3) '-' TMR_name{i}(2:size(CY5_name{i},2)-3) '_img' num2str(k) '.mat']);
                load(fileNUC,'Label_mid'); %'NucInfo','CellInfo',
                end
                if not(isempty(TMR_name)) & TMR_hascloud & perform_cell_functions   %Load TMR clouds
                    filename = ['mRNA_' dates0{i} '_' strain TMR_name{i} num2str(thTMRall{i}(thNum_TMR)) '_' num2str(thAF594all{i}(thNum_TMR)) CY5_name{i} num2str(thCY5all{i}(thNum_TMR)) '_' times0{i}{j} '_im' num2str(k) '_non_stringent_cloud.mat'];
%                     filename = ['mRNA_' dates0{i} '_' strain TMR_name{i} num2str(thTMRall{i}(thNum_TMR)) '_' num2str(thAF594all{i}(thNum_TMR)) CY5_name{i} num2str(thCY5all{i}(thNum_TMR)) '_' times0{i}{j} '_im' num2str(k)];
                    % AF594_name{i} num2str(thAF594all{i}(thNum_TMR))
                    load(filename,['PARtmr_' nucTH]);
                end
                if not(isempty(CY5_name))
                    filename = ['mRNA_' dates0{i} '_' strain TMR_name{i} num2str(thTMRall{i}(thNum_CY5)) '_' num2str(thAF594all{i}(thNum_CY5)) CY5_name{i} num2str(thCY5all{i}(thNum_CY5)) '_' times0{i}{j} '_im' num2str(k) '_non_stringent_cloud.mat'];
%                     filename = ['mRNA_' dates0{i} '_' strain TMR_name{i} num2str(thTMRall{i}(thNum_CY5)) '_' num2str(thAF594all{i}(thNum_CY5)) CY5_name{i} num2str(thCY5all{i}(thNum_CY5)) '_' times0{i}{j} '_im' num2str(k)];
                    %AF594_name{i} num2str(thAF594all{i}(thNum_CY5))
                    load(filename,['PARcy5_' nucTH]);
                    if size(PARcy5_mid.TotExpInt,2) == 1        %fives another column if there is only one
                        PARcy5_mid.TotExpInt(1:size(PARcy5_mid.TotExpInt,1),2) = NaN;
                    end   
                    %figure(2); clf; imshow(max(PARcy5_mid.clouds1,[],3),[])
                    for m = 1: size(PARcy5_mid.TotExpInt,1) %Go through every cell in this image
%                         m
                        CY5_AF594_TMR(counter9,5) = str2num(dates0{i})/10^6;  %5: date of experiment
                        CY5_AF594_TMR_tp{counter11}(counter10,5) = str2num(dates0{i})/10^6; %5: date of experiment
                        CY5_AF594_TMR(counter9,6) = k;                      %6: image number
                        CY5_AF594_TMR_tp{counter11}(counter10,6) = k;       %6: image number
                        CY5_AF594_TMR(counter9,7) = m;                      %7: cell number
                        CY5_AF594_TMR_tp{counter11}(counter10,7) = m;       %7: cell number                        
                        if perform_cell_functions
                        %%%Determine cell boundaries
                        A1 = size(cells,2);                            %limits of the segmented image
                        A2 = size(cells,1);                             %limits of the segmented image
                        k1 = Lab == m ; %== j;%sieve out dots in cell j.
                        %figure; imshow(k1,[0 1]);
                        k1=uint16(k1);
                        k3=regionprops(k1,'BoundingBox','Area');
                        k4 = k3.BoundingBox; %create the rectangular box around the cell.
                        X0=round(k4(1))-4;
                        Y0=round(k4(2))-4;
                        X1=round(k4(1)+ k4(3))+4;
                        Y1=round(k4(2)+ k4(4))+4;
                        if X0 < 1;  X0 = 1; end
                        if Y0 < 1;  Y0 = 1; end
                        if X1 > A1; X1 = A1; end
                        if Y1 > A2; Y1 = A2; end
                        if CY5_hascloud | TMR_hascloud
                            if CY5_hascloud
                                CY5_cloud1_cell = PARcy5_mid.clouds1(Y0:Y1,X0:X1,:);      %Generate just CY5 cloud 1 for this cell
                                CY5_cloud2_cell = PARcy5_mid.clouds2(Y0:Y1,X0:X1,:);      %Generate CY5 cloud 2 for this cell
                                CY5_AF594_TMR(counter9,16) = sum(CY5_cloud1_cell(:));            %Find volume of first CY5 cloud in pixels
                                CY5_AF594_TMR(counter9,17) = sum(CY5_cloud2_cell(:));            %Find volume of second CY5 cloud in pixels
                                CY5_AF594_TMR_tp{counter11}(counter10,16) = sum(CY5_cloud1_cell(:));            %Find volume of first CY5 cloud in pixels
                                CY5_AF594_TMR_tp{counter11}(counter10,17) = sum(CY5_cloud2_cell(:));            %Find volume of second CY5 cloud in pixels
                            end
                            if not(isempty(TMR_name)) & TMR_hascloud
                                TMR_cloud1_cell = PARtmr_mid.clouds1(Y0:Y1,X0:X1,:);      %Generate TMR cloud 1 for this cell
                                TMR_cloud2_cell = PARtmr_mid.clouds2(Y0:Y1,X0:X1,:);      %Generate TMR cloud 2 for this cell
                                CY5_AF594_TMR(counter9,18) = sum(TMR_cloud1_cell(:));            %Find volume of first CY5 cloud in pixels
                                CY5_AF594_TMR(counter9,19) = sum(TMR_cloud2_cell(:));            %Find volume of second CY5 cloud in pixels
                                CY5_AF594_TMR_tp{counter11}(counter10,18) = sum(TMR_cloud1_cell(:));            %Find volume of first CY5 cloud in pixels
                                CY5_AF594_TMR_tp{counter11}(counter10,19) = sum(TMR_cloud2_cell(:));            %Find volume of second CY5 cloud in pixels
                            end
                        end
                        %                                 if sum(CY5_cloud1_cell(:)) > 50                           %to check how mnay clouds in image
                        %                                     ['cell ' num2str(m) ' has a cloud']
                        %                                     figure; clf; imshow(max(CY5_cloud1_cell,[],3),[])
                        %                                 end
                        %                             end
                        %                         end
                        end
                        %%%Eliminate dead pixels by looking for many spots with the exact same x position (needs to be improved)
                        for p = 1:size(PARcy5_mid.xinit,2)
                            same_pos = find(PARcy5_mid.xinit(m,:) == PARcy5_mid.xinit(m,p));
                            if size(same_pos,2)>10
                                for q = 1:size(same_pos,2)
                                    PARcy5_mid.xinit(m,same_pos(q)) = NaN;
                                    PARcy5_mid.yinit(m,same_pos(q)) = NaN;
                                    PARcy5_mid.zinit(m,same_pos(q)) = NaN;
                                    PARcy5_mid.TotExpInt(m,same_pos(q)) = NaN;
                                end
                            end
                        end
                        %%% Below removes duplicates of the same
                        %%% transcript. Goes through every transcript,
                        %%% finds distance between it and other
                        %%% transcripts, and eliminates others
                        for p = 1:size(PARcy5_mid.xinit,2)
                            dists = zeros(1,size(size(PARcy5_mid.xinit,2),3)); %stores distance from transcript to others
                            if not(isnan(PARcy5_mid.xinit(m,p)))
                                for pp = 1:size(PARcy5_mid.xinit,2)    %find distance squared
                                    dists(1,pp,1) = (PARcy5_mid.xinit(m,p)- PARcy5_mid.xinit(m,pp))^2;
                                    dists(1,pp,2) = (PARcy5_mid.yinit(m,p)- PARcy5_mid.yinit(m,pp))^2;
                                    dists(1,pp,3) = ((PARcy5_mid.zinit(m,p)- PARcy5_mid.zinit(m,pp))*z_adjust)^2;
                                end
                                distF = sqrt(dists(1,:,1)+dists(1,:,2)+dists(1,:,3));
                                for pp = 1:size(PARcy5_mid.xinit,2)
                                    if p ~= pp & distF(pp) < same_dist
                                        if PARcy5_mid.TotExpInt(m,pp) < PARcy5_mid.TotExpInt(m,p) %Only keep the brighter spot
                                        PARcy5_mid.xinit(m,pp) = NaN;
                                        PARcy5_mid.yinit(m,pp) = NaN;
                                        PARcy5_mid.zinit(m,pp) = NaN;
                                        PARcy5_mid.TotExpInt(m,p) = PARcy5_mid.TotExpInt(m,p)+PARcy5_mid.TotExpInt(m,pp);
                                        PARcy5_mid.TotExpInt(m,pp) = NaN;
                                        else
                                        PARcy5_mid.xinit(m,p) = NaN;
                                        PARcy5_mid.yinit(m,p) = NaN;
                                        PARcy5_mid.zinit(m,p) = NaN;
                                        PARcy5_mid.TotExpInt(m,pp) = PARcy5_mid.TotExpInt(m,p)+PARcy5_mid.TotExpInt(m,pp);
                                        PARcy5_mid.TotExpInt(m,p) = NaN; 
                                        break
                                        end
                                    end
                                end
                            end
                        end
                        %%%store the positions of all transcripts (initial
                        %%%positions)
                        CY5_pos_init(counter9,1:size(PARcy5_mid.xinit,2),1) = PARcy5_mid.xinit(m,1:size(PARcy5_mid.xinit,2));
                        CY5_pos_init(counter9,1:size(PARcy5_mid.xinit,2),2) = PARcy5_mid.yinit(m,1:size(PARcy5_mid.yinit,2));
                        CY5_pos_init(counter9,1:size(PARcy5_mid.xinit,2),3) = PARcy5_mid.zinit(m,1:size(PARcy5_mid.zinit,2));
                        CY5_pos_init_tp{counter11}(counter10,1:size(PARcy5_mid.xinit,2),1) = PARcy5_mid.xinit(m,1:size(PARcy5_mid.xinit,2));
                        CY5_pos_init_tp{counter11}(counter10,1:size(PARcy5_mid.xinit,2),2) = PARcy5_mid.yinit(m,1:size(PARcy5_mid.yinit,2));
                        CY5_pos_init_tp{counter11}(counter10,1:size(PARcy5_mid.xinit,2),3) = PARcy5_mid.zinit(m,1:size(PARcy5_mid.zinit,2));
                        %% Determine distance to edge of nucleus
                        if perform_nuc_functions
                        Label_mid_cell = Label_mid(Y0:Y1,X0:X1,:);  %Get nucleus for just the cell
                        NucBorder1 = bwmorph3(Label_mid_cell,'remove');
                        %                         figure(300); clf;
                        %                         imshow(NucBorder1(:,:,round(size(NucBorder1,3)/2)),[]) %Visualize middle slice
                        %                         for zzz = 1:size(NucBorder1,3)    %All slices
                        %                             figure(300); clf; imshow(NucBorder1(:,:,zzz),[])
                        %                             pause(.5)
                        %                         end
                        lin_Nuc = find(NucBorder1 == 1);[xNuc,yNuc,zNuc] = ind2sub(size(NucBorder1),lin_Nuc(1:4:end)); clear lin_Nuc; %These are vertical vectors (x differs)
                        pNuc1 = xNuc; pNuc1(:,:,2) = yNuc; pNuc1(:,:,3) = zNuc; clear xNuc yNuc zNuc NucBorder1; %combine matrices
                        pNuc1 = int16(pNuc1);                        
                        pNuc2 = repmat(pNuc1,[1,size(PARcy5_mid.xinit,2),1]);  %replicate in y axis the matricies
                        pTran = PARcy5_mid.xinit(m,1:size(PARcy5_mid.xinit,2));
                        pTran(:,:,2) = PARcy5_mid.yinit(m,1:size(PARcy5_mid.yinit,2));
                        pTran(:,:,3) = PARcy5_mid.zinit(m,1:size(PARcy5_mid.zinit,2));
                        pTran = int16(pTran);
                        pTran1 = repmat(pTran,[size(pNuc2,1),1,1]);
                        pDiff = pNuc2-pTran1; clear pNuc2 pTran1; pDiff = int32(pDiff);
                        pDiff(:,:,3) = pDiff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
                        pDiff = sum(pDiff.^2,3); %pDiff = pDiff.^.5; %To keep in in integer form, do not square root
                        CY5_Nuc_dist(counter9,1:size(PARcy5_mid.xinit,2)) = min(pDiff,[],1);
                        CY5_Nuc_dist_tp{counter11}(counter10,1:size(PARcy5_mid.xinit,2)) = min(pDiff,[],1);
                        clear pDiff pNuc2
                        if CY5_hascloud
                            if not(isnan(PARcy5_mid.xinit(m,1)))
                                C1_bord = bwmorph3(CY5_cloud1_cell,'remove');
                                %                                 figure(301); clf;
                                %                                 imshow(C1_bord(:,:,round(size(C1_bord,3)/2)),[]) %Visualize middle slice
                                %                                 for zzz = 1:size(C1_bord,3)    %All slices
                                %                                     figure(301); clf; imshow(C1_bord(:,:,zzz),[])
                                %                                     title(['slice ' num2str(zzz)])
                                %                                     pause(.5)
                                %                                 end
                                lin_C1 = find(C1_bord == 1); [xC1,yC1,zC1] =  ind2sub(size(C1_bord),lin_C1(1:3:end)); clear lin_C1; %These are vertical vectors (x differs)
                                xC1 = xC1'; yC1 = yC1'; zC1 = zC1';     %Change x's to y's
                                pC1 = xC1; pC1(:,:,2) = yC1; pC1(:,:,3) = zC1; clear xC1 yC1 zC1 %populate one matrix
                                pC1 = int16(pC1);
                                %%% fast but memory-intensive way of determining distances
                                if size(pC1,2) < 1
                                    Cl_Nuc_dist(counter9,1) = NaN;
                                    Cl_Nuc_dist_tp{counter11}(counter10,1) = NaN;
                                    C1_CY5_dist(counter9,1:size(PARcy5_mid.xinit,2)) = nan(1,size(PARcy5_mid.xinit,2));
                                    C1_CY5_dist_tp{counter11}(counter10,1:size(PARcy5_mid.xinit,2)) = nan(1,size(PARcy5_mid.xinit,2));                                      
                                else
                                    pC1_2 = repmat(pC1,[size(pNuc1,1),1,1]); %make copies for comparison. Commented outfor memory
                                    pNuc2 = repmat(pNuc1,[1,size(pC1,2),1]);        %Make copies for comparison
                                    pDiff = pNuc2-pC1_2; clear pNuc2 pC1_2; pDiff = int32(pDiff);
                                    pDiff(:,:,3) = pDiff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
                                    pDiff = sum(pDiff.^2,3); %pDiff = pDiff.^.5;
                                    Cl_Nuc_dist(counter9,1) = min(pDiff(:));
                                    Cl_Nuc_dist_tp{counter11}(counter10,1) = min(pDiff(:));
                                    clear pDiff
                                    %%%determine distances to each transcript
                                    pTran = PARcy5_mid.xinit(m,1:size(PARcy5_mid.xinit,2))';
                                    pTran(:,:,2) = PARcy5_mid.yinit(m,1:size(PARcy5_mid.yinit,2))';
                                    pTran(:,:,3) = PARcy5_mid.zinit(m,1:size(PARcy5_mid.zinit,2))';
                                    pTran = int16(pTran);
                                    pTran1 =  repmat(pTran,[1,size(pC1,2),1]);        %Make copies for comparison
                                    pC1_2 = repmat(pC1,[size(pTran1,1),1,1]); %make copies for comparison.
                                    pDiff = pTran1-pC1_2; clear pC1_2 pTran1; pDiff = int32(pDiff);
                                    pDiff(:,:,3) = pDiff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
                                    pDiff = sum(pDiff.^2,3); %pDiff = pDiff.^.5; 
                                    C1_CY5_dist(counter9,1:size(PARcy5_mid.xinit,2)) = min(pDiff,[],2)';
                                    C1_CY5_dist_tp{counter11}(counter10,1:size(PARcy5_mid.xinit,2)) = min(pDiff,[],2)';
                                    clear pDiff
                                end
                                %%%
                                %%% Slow but memory-efficient way
%                                 if size(pC1,2) < 1
%                                     Cl_Nuc_dist(counter9,2) = NaN;
%                                     Cl_Nuc_dist_tp{counter11}(counter10,2) = NaN;
%                                 else
%                                     
%                                     pDiff = zeros(size(pC1,2),1);
%                                     for cl_dim = 1:size(pC1,2)  %less memory-intensive, but slower way for determining distances
%                                         pC1_temp = pC1(:,cl_dim,:);
%                                         pC1_temp = repmat(pC1_temp,[size(pNuc1,1),1,1]); %make copies for comparison.
%                                         temp_Diff = pC1_temp-pNuc1;
%                                         temp_Diff(:,:,3) = temp_Diff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
%                                         temp_Diff = sum(temp_Diff.^2,3); temp_Diff = temp_Diff.^.5;
%                                         pDiff(cl_dim) = min(temp_Diff(:));
%                                     end
%                                     Cl_Nuc_dist(counter9,1) = min(pDiff(:));
%                                     Cl_Nuc_dist_tp{counter11}(counter10,1) = min(pDiff(:));
%                                 end
                                %%%
                            else                                           %If no cloud, make NaN
                                Cl_Nuc_dist(counter9,1) = NaN;
                                Cl_Nuc_dist_tp{counter11}(counter10,1) = NaN;
                                C1_CY5_dist(counter9,1:size(PARcy5_mid.xinit,2)) = nan(1,size(PARcy5_mid.xinit,2));
                                C1_CY5_dist_tp{counter11}(counter10,1:size(PARcy5_mid.xinit,2)) = nan(1,size(PARcy5_mid.xinit,2));                                                            
                            end
                            if not(isnan(PARcy5_mid.xinit(m,2)))               %Do the same, but for the second cloud
                                C1_bord = bwmorph3(CY5_cloud2_cell,'remove');
                                %                                 figure(301); clf;
                                %                                 imshow(C1_bord(:,:,round(size(C1_bord,3)/2)),[]) %Visualize middle slice
                                %                                 for zzz = 1:size(C1_bord,3)    %All slices
                                %                                     figure(301); clf; imshow(C1_bord(:,:,zzz),[])
                                %                                     title(['slice ' num2str(zzz)])
                                %                                     pause(.5)
                                %                                 end
                                lin_C1 = find(C1_bord == 1); [xC1,yC1,zC1] =  ind2sub(size(C1_bord),lin_C1(1:2:end)); clear lin_C1; %These are vertical vectors (x differs)
                                xC1 = xC1'; yC1 = yC1'; zC1 = zC1';     %Change x's to y's
                                pC1 = xC1; pC1(:,:,2) = yC1; pC1(:,:,3) = zC1; clear xC1 yC1 zC1 %populate one matrix
                                pC1 = int16(pC1);
                                %%% fast but memory-intensive way of determining distances
                                if size(pC1,2) < 1
                                    Cl_Nuc_dist(counter9,2) = NaN;
                                    Cl_Nuc_dist_tp{counter11}(counter10,2) = NaN;
                                    C2_CY5_dist(counter9,1:size(PARcy5_mid.xinit,2)) = nan(1,size(PARcy5_mid.xinit,2));
                                    C2_CY5_dist_tp{counter11}(counter10,1:size(PARcy5_mid.xinit,2)) = nan(1,size(PARcy5_mid.xinit,2));                                     
                                else
                                    pC1_2 = repmat(pC1,[size(pNuc1,1),1,1]); %make copies for comparison. Commented outfor memory
                                    pNuc2 = repmat(pNuc1,[1,size(pC1,2),1]);        %Make copies for comparison
                                    pDiff = pNuc2-pC1_2; clear pNuc2 pC1_2; pDiff = int32(pDiff);
                                    pDiff(:,:,3) = pDiff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
                                    pDiff = sum(pDiff.^2,3);% pDiff = pDiff.^.5;
                                    Cl_Nuc_dist(counter9,2) = min(pDiff(:));
                                    Cl_Nuc_dist_tp{counter11}(counter10,2) = min(pDiff(:));
                                    clear pDiff
                                    %%%determine distances to each transcript
                                    pTran1 =  repmat(pTran,[1,size(pC1,2),1]);        %Make copies for comparison
                                    pC1_2 = repmat(pC1,[size(pTran1,1),1,1]); %make copies for comparison.
                                    pDiff = pTran1-pC1_2; clear pC1_2 pTran1; pDiff = int32(pDiff);
                                    pDiff(:,:,3) = pDiff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
                                    pDiff = sum(pDiff.^2,3); %pDiff = pDiff.^.5;
                                    C2_CY5_dist(counter9,1:size(PARcy5_mid.xinit,2)) = min(pDiff,[],2)';
                                    C2_CY5_dist_tp{counter11}(counter10,1:size(PARcy5_mid.xinit,2)) = min(pDiff,[],2)';                                    
                                    clear pDiff
                                end
                                %%%
                                %%% Slow but memory-efficient way
%                                 if size(pC1,2) < 1
%                                     Cl_Nuc_dist(counter9,2) = NaN;
%                                     Cl_Nuc_dist_tp{counter11}(counter10,2) = NaN;
%                                 else
%                                     pDiff = zeros(size(pC1,2),1);
%                                     for cl_dim = 1:size(pC1,2)  %less memory-intensive, but slower way for determining distances
%                                         pC1_temp = pC1(:,cl_dim,:);
%                                         pC1_temp = repmat(pC1_temp,[size(pNuc1,1),1,1]); %make copies for comparison.
%                                         temp_Diff = pC1_temp-pNuc1;
%                                         temp_Diff(:,:,3) = temp_Diff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
%                                         temp_Diff = sum(temp_Diff.^2,3); temp_Diff = temp_Diff.^.5;
%                                         pDiff(cl_dim) = min(temp_Diff(:));
%                                     end
%                                     Cl_Nuc_dist(counter9,2) = min(pDiff(:));
%                                     Cl_Nuc_dist_tp{counter11}(counter10,2) = min(pDiff(:));
%                                 end
                                %%%
                            else
                                Cl_Nuc_dist(counter9,2) = NaN;
                                Cl_Nuc_dist_tp{counter11}(counter10,2) = NaN;
                                C2_CY5_dist(counter9,1:size(PARcy5_mid.xinit,2)) = nan(1,size(PARcy5_mid.xinit,2));
                                C2_CY5_dist_tp{counter11}(counter10,1:size(PARcy5_mid.xinit,2)) = nan(1,size(PARcy5_mid.xinit,2));                                                           
                            end
                        end
                        end
                        %%
                        %CY5_C1_dist = Distance of each CY5 transcript from the first cloud (arranged like CY5_AF594_TMR)
                        %CY5_C1_dist_tp = Distance of each CY5 transcript from the first cloud (arranged like CY5_AF594_TMR_tp)
                        %CY5_C2_dist = Distance of each CY5 transcript from the 2nd cloud (arranged like CY5_AF594_TMR)
                        %CY5_C2_dist_tp = Distance of each CY5 transcript from the 2nd cloud (arranged like CY5_AF594_TMR_tp)
                        %TMR_C1_dist = Distance of each TMR transcript from the first cloud (arranged like CY5_AF594_TMR)
                        %TMR_C1_dist_tp = Distance of each TMR transcript from the first cloud (arranged like CY5_AF594_TMR_tp)
                        %TMR_C2_dist = Distance of each TMR transcript from the 2nd cloud (arranged like CY5_AF594_TMR)
                        %TMR_C2_dist_tp = Distance of each TMR transcript from the 2nd cloud (arranged like CY5_AF594_TMR_tp)
                        %%%
                        transc_num = size(find(PARcy5_mid.xinit(m,:)>0),2);
                        CY5_int(counter9,3:size(PARcy5_mid.TotExpInt,2)) = PARcy5_mid.TotExpInt(m,3:size(PARcy5_mid.TotExpInt,2));  %Store total experimental intensity
                        CY5_int(counter9,1:2) = PARcy5_mid.TotExpInt(m,1:2)/cloud_correction;
                        CY5_int_tp{counter11}(counter10,3:size(PARcy5_mid.TotExpInt,2)) = PARcy5_mid.TotExpInt(m,3:size(PARcy5_mid.TotExpInt,2)); %Store total experimental intensity
                        CY5_int_tp{counter11}(counter10,1:2) = PARcy5_mid.TotExpInt(m,1:2)/cloud_correction;
                        CY5_AF594_TMR(counter9,1) = transc_num; %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        CY5_AF594_TMR_tp{counter11}(counter10,1) = transc_num; %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        % CY5_AF594_TMR(counter9,1) = nansum([nansum(PARcy5_mid.cytoRNA(m,:)),nansum(PARcy5_mid.nucRNA(m))]); %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        %CY5_AF594_TMR_tp{counter11}(counter10,1) = nansum([nansum(PARcy5_mid.cytoRNA(m,:)),nansum(PARcy5_mid.nucRNA(m))]); %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        %CY5_cloud1 CY5_cloud1_tp CY5_cloud2 CY5_cloud2_tp TMR_cloud1 TMR_cloud1_tp
                        %TMR_cloud2 TMR_cloud2_tp
                        %% Determine whether transcripts are inside clouds
                        if perform_cell_functions
                        if CY5_hascloud
                            for temp = 1:size(PARcy5_mid.xinit,2)
                                trans_loc = [CY5_pos_init(counter9,temp,1) CY5_pos_init(counter9,temp,2) CY5_pos_init(counter9,temp,3)]; %temporary x,y,z positions
                                if not(isnan(trans_loc(1))) %Check that it is not NaN
                                    part_c1 = CY5_cloud1_cell(trans_loc(1),trans_loc(2),trans_loc(3));  %This needs to be adjusted. clouds1 refers to entire image while locations are for just the cell
                                    part_c2 = CY5_cloud2_cell(trans_loc(1),trans_loc(2),trans_loc(3));
                                    if part_c1 == 1 & part_c2 == 1
                                        CY5_in_CY5cloud(counter9,temp) = 4;
                                        CY5_in_CY5cloud_tp{counter11}(counter10,temp) = 4;
                                    elseif part_c1 == 1
                                        CY5_in_CY5cloud(counter9,temp) = 2;
                                        CY5_in_CY5cloud_tp{counter11}(counter10,temp) = 2;
                                    elseif part_c2 == 1
                                        CY5_in_CY5cloud(counter9,temp) = 3;
                                        CY5_in_CY5cloud_tp{counter11}(counter10,temp) = 3;
                                    else
                                        CY5_in_CY5cloud(counter9,temp) = 1;
                                        CY5_in_CY5cloud_tp{counter11}(counter10,temp) = 1;
                                    end
                                else
                                    CY5_in_CY5cloud(counter9,temp) = 0;
                                    CY5_in_CY5cloud_tp{counter11}(counter10,temp) = 0;
                                end
                            end
                        end
                        if TMR_hascloud
                            for temp = 1:size(PARcy5_mid.xinit,2)
                                trans_loc = [CY5_pos_init(counter9,temp,1) CY5_pos_init(counter9,temp,2) CY5_pos_init(counter9,temp,3)]; %temporary x,y,z positions
                                if not(isnan(trans_loc(1))) %Check that it is not NaN
                                    part_c1 = TMR_cloud1_cell(trans_loc(1),trans_loc(2),trans_loc(3));  %This needs to be adjusted. clouds1 refers to entire image while locations are for just the cell
                                    part_c2 = TMR_cloud2_cell(trans_loc(1),trans_loc(2),trans_loc(3));
                                    if part_c1 == 1 & part_c2 == 1
                                        CY5_in_TMRcloud(counter9,temp) = 4;
                                        CY5_in_TMRcloud_tp{counter11}(counter10,temp) = 4;
                                    elseif part_c1 == 1
                                        CY5_in_TMRcloud(counter9,temp) = 2;
                                        CY5_in_TMRcloud_tp{counter11}(counter10,temp) = 2;
                                    elseif part_c2 == 1
                                        CY5_in_TMRcloud(counter9,temp) = 3;
                                        CY5_in_TMRcloud_tp{counter11}(counter10,temp) = 3;
                                    else
                                        CY5_in_TMRcloud(counter9,temp) = 1;
                                        CY5_in_TMRcloud_tp{counter11}(counter10,temp) = 1;
                                    end
                                else
                                    CY5_in_TMRcloud(counter9,temp) = 0;
                                    CY5_in_TMRcloud_tp{counter11}(counter10,temp) = 0;
                                end
                            end
                        end
                        end
                        %%%
                        if CY5_hascloud && not(isnan(PARcy5_mid.TotExpInt(m,1))) && not(isnan(PARcy5_mid.TotExpInt(m,2)))
                            CY5_AF594_TMR(counter9,1) = CY5_AF594_TMR(counter9,1)-2 + PARcy5_mid.TotExpInt(m,1)/CY5_cloudInt(i)/cloud_correction+ PARcy5_mid.TotExpInt(m,2)/CY5_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR_tp{counter11}(counter10,1) =CY5_AF594_TMR_tp{counter11}(counter10,1) - 2 + PARcy5_mid.TotExpInt(m,1)/CY5_cloudInt(i)/cloud_correction + PARcy5_mid.TotExpInt(m,2)/CY5_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR(counter9,8) = sum(PARcy5_mid.TotExpInt(m,1:2) > CY5_tran_thres*CY5_cloudInt(i)*cloud_correction)+sum(PARcy5_mid.TotExpInt(m,3:size(PARcy5_mid.TotExpInt,2)) > CY5_tran_thres*CY5_cloudInt(i));        %Number of transcription sites/clouds in the cell
                            CY5_AF594_TMR_tp{counter11}(counter10,8) =CY5_AF594_TMR(counter9,8);      %Number of transcription sites/clouds in the cell
                            CY5_AF594_TMR(counter9,11) = PARcy5_mid.TotExpInt(m,1)/CY5_cloudInt(i)/cloud_correction;  %Number of transcripts inside the first cloud
                            CY5_AF594_TMR(counter9,12) = PARcy5_mid.TotExpInt(m,2)/CY5_cloudInt(i)/cloud_correction; %Number of transcripts inside the second cloud
                            CY5_AF594_TMR_tp{counter11}(counter10,11) = CY5_AF594_TMR(counter9,11);    %Number of transcripts inside the first cloud (timepoint sorted)
                            CY5_AF594_TMR_tp{counter11}(counter10,12) = CY5_AF594_TMR(counter9,12);     %Number of transcripts inside the second cloud (timepoint sorted)
                            %%%Below adds position information for the
                            %%%transcription sites                            
                            counter20 = 1;  %Counter for transcription sites
                            Tsite_ind = zeros(1,4); %Will have index for transcription sites                            
                            if CY5_AF594_TMR(counter9,11) >= CY5_tran_thres
                                Tsite_ind(counter20) = 1;
                                counter20 = counter20+1;
                            end
                            if CY5_AF594_TMR(counter9,12) >= CY5_tran_thres
                                Tsite_ind(counter20) = 2;
                                counter20 = counter20+1;
                            end
                            [temp_int,int_ind] = sort(CY5_int(counter9,3:size(PARcy5_mid.TotExpInt,2)),'descend');   %Sort the transcript intensities to find the transcription sites
                            int_ind = int_ind+2;    %Adjust for leaving out first two
                            Poss_tsites = find(temp_int >= CY5_tran_thres*CY5_cloudInt(i));
                            Tsites_left = 5-counter20;
                            if not(isempty(Poss_tsites)) & Tsites_left > 0
                                Tsites_add = min([size(Poss_tsites,2) Tsites_left]);    %See how many transcription sites to add
                                Tsite_ind(counter20:counter20+Tsites_add-1) = int_ind(Poss_tsites(1:Tsites_add)); %Add the indexes of the transcription sites
                            end
                            for index1 = 1:4
                                if Tsite_ind(index1) >0
                                    CY5_AF594_TMR(counter9,(20:22)+3*(index1-1)) = CY5_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                                    CY5_AF594_TMR_tp{counter11}(counter10,(20:22)+3*(index1-1)) = CY5_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                                end
                            end
                            %%%
                        elseif CY5_hascloud && not(isnan(PARcy5_mid.TotExpInt(m,1)))
                            CY5_AF594_TMR(counter9,1) = CY5_AF594_TMR(counter9,1)-1 + PARcy5_mid.TotExpInt(m,1)/CY5_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR_tp{counter11}(counter10,1) =CY5_AF594_TMR_tp{counter11}(counter10,1) - 1 + PARcy5_mid.TotExpInt(m,1)/CY5_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR(counter9,8) = sum(PARcy5_mid.TotExpInt(m,1) > CY5_tran_thres*CY5_cloudInt(i)*cloud_correction)+sum(PARcy5_mid.TotExpInt(m,3:size(PARcy5_mid.TotExpInt,2)) > CY5_tran_thres*CY5_cloudInt(i));        %Number of transcription sites/clouds in the cell
                            CY5_AF594_TMR_tp{counter11}(counter10,8) =CY5_AF594_TMR(counter9,8);      %Number of transcription sites/clouds in the cell
                            CY5_AF594_TMR(counter9,11) = PARcy5_mid.TotExpInt(m,1)/CY5_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR_tp{counter11}(counter10,11) = CY5_AF594_TMR(counter9,11);
                            %%%Below adds position information for the
                            %%%transcription sites
                            counter20 = 1;  %Counter for transcription sites
                            Tsite_ind = zeros(1,4); %Will have index for transcription sites                            
                            if CY5_AF594_TMR(counter9,11) >= CY5_tran_thres
                                Tsite_ind(counter20) = 1;
                                counter20 = counter20+1;
                            end
                            if CY5_AF594_TMR(counter9,12) >= CY5_tran_thres
                                Tsite_ind(counter20) = 2;
                                counter20 = counter20+1;
                            end
                            [temp_int,int_ind] = sort(CY5_int(counter9,3:size(PARcy5_mid.TotExpInt,2)),'descend');   %Sort the transcript intensities to find the transcription sites
                            int_ind = int_ind+2;    %Adjust for leaving out first two
                            Poss_tsites = find(temp_int >= CY5_tran_thres*CY5_cloudInt(i));
                            Tsites_left = 5-counter20;
                            if not(isempty(Poss_tsites)) & Tsites_left > 0
                                Tsites_add = min([size(Poss_tsites,2) Tsites_left]);    %See how many transcription sites to add
                                Tsite_ind(counter20:counter20+Tsites_add-1) = int_ind(Poss_tsites(1:Tsites_add)); %Add the indexes of the transcription sites
                            end
                             for index1 = 1:4
                                if Tsite_ind(index1) >0
                            CY5_AF594_TMR(counter9,(20:22)+3*(index1-1)) = CY5_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                            CY5_AF594_TMR_tp{counter11}(counter10,(20:22)+3*(index1-1)) = CY5_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                                end
                             end
                             %%%
                        else
                            CY5_AF594_TMR(counter9,8) = sum(PARcy5_mid.TotExpInt(m,:) > CY5_tran_thres*CY5_cloudInt(i));        %Number of transcription sites/clouds in the cell
                            CY5_AF594_TMR_tp{counter11}(counter10,8) =  CY5_AF594_TMR(counter9,8);        %Number of transcription sites/clouds in the cell
                            %%%Below adds position information for the
                            %%%transcription sites
                            Tsite_ind = zeros(1,4); %Will have index for transcription sites                            
                            [temp_int,int_ind] = sort(CY5_int(counter9,1:size(PARcy5_mid.TotExpInt,2)),'descend');   %Sort the transcript intensities to find the transcription sites
                            Poss_tsites = find(temp_int >= CY5_tran_thres*CY5_cloudInt(i));
                            Tsites_left = 4;
                            if not(isempty(Poss_tsites)) & Tsites_left > 0
                                Tsites_add = min([size(Poss_tsites,2) Tsites_left]);    %See how many transcription sites to add
                                Tsite_ind(1:Tsites_add) = int_ind(Poss_tsites(1:Tsites_add)); %Add the indexes of the transcription sites
                            end
                            for index1 = 1:4
                                if Tsite_ind(index1) >0
                                    CY5_AF594_TMR(counter9,(20:22)+3*(index1-1)) = CY5_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                                    CY5_AF594_TMR_tp{counter11}(counter10,(20:22)+3*(index1-1)) = CY5_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                                end
                            end
                            %%%
                        end
                        %CY5_transc_site TMR_transc_site
                        CY5_AF594_TMR(counter9,4) = times1(counter11);  %stores timepoint
                        CY5_AF594_TMR(counter9,5) = str2num(dates0{i})*.000001; %stores date of experiment
                        CY5_AF594_TMR(counter9,6) = k; %stores image number
                        counter9 = counter9+1;      %increase counter for cell (all)
                        counter10 = counter10+1;    %increase counter for cell (resets each timepoint)
                    end
                    counter9 = counter9_int;
                    counter10 = counter10_int;
                end
                %filename = ['mRNA_' dates0{i} '_STL1-TMR-th27_AF594-GPP2-th17_CTT1-GPD1-th21_' times0{i}{j} '_im' num2str(k)];
                if not(isempty(AF594_name))
                    filename = ['mRNA_' dates0{i} '_' strain TMR_name{i} num2str(thTMRall{i}(thNum_AF594)) AF594_name{i} num2str(thAF594all{i}(thNum_AF594)) CY5_name{i} num2str(thCY5all{i}(thNum_AF594)) '_' times0{i}{j} '_im' num2str(k) '_non_stringent_cloud.mat'];
%                     filename = ['mRNA_' dates0{i} '_' strain TMR_name{i} num2str(thTMRall{i}(thNum_AF594)) AF594_name{i} num2str(thAF594all{i}(thNum_AF594)) CY5_name{i} num2str(thCY5all{i}(thNum_AF594)) '_' times0{i}{j} '_im' num2str(k)];
                    load(filename,['PARaf594_' nucTH]);
                    for m = 1: size(PARaf594_mid.TotExpInt,1) %Go through every cell in this image
                        %%%Eliminate dead pixels  (needs to be improved)
                        for p = 1:size(PARaf594_mid.xinit,2)
                            same_pos = find(PARaf594_mid.xinit(m,:) == PARaf594_mid.xinit(m,p));
                            if size(same_pos,2)>10
                                for q = 1:size(same_pos,2)
                                    PARaf594_mid.xinit(m,same_pos(q)) = NaN;
                                end
                            end
                        end
                        %%% Below removes duplicates of the same
                        %%% transcript. Goes through every transcript,
                        %%% finds distance between it and other
                        %%% transcripts, and eliminates others
                        for p = 1:size(PARaf594_mid.xinit,2)
                            dists = zeros(1,size(size(PARaf594_mid.xinit,2),3)); %stores distance from transcript to others
                            if not(isnan(PARcy5_mid.xinit(m,p)))
                                for pp = 1:size(PARaf594_mid.xinit,2)    %find distance squared
                                    dists(1,pp,1) = (PARaf594_mid.xinit(m,p)- PARaf594_mid.xinit(m,pp))^2;
                                    dists(1,pp,2) = (PARaf594_mid.yinit(m,p)- PARaf594_mid.yinit(m,pp))^2;
                                    dists(1,pp,3) = ((PARaf594_mid.zinit(m,p)- PARaf594_mid.zinit(m,pp))*3)^2;
                                end
                                distF = sqrt(dists(1,:,1)+dists(1,:,2)+dists(1,:,3));
                                for pp = 1:size(PARaf594_mid.xinit,2)
                                    if p ~= pp & distF(pp) < same_dist
                                        PARaf594_mid.xinit(m,pp) = NaN;
                                        PARaf594_mid.yinit(m,pp) = NaN;
                                        PARaf594_mid.zinit(m,pp) = NaN;
                                        PARaf594_mid.TotExpInt(m,pp) = NaN;
                                    end
                                end
                            end
                        end
                        %%%
                        transc_num = size(find(PARaf594_mid.xinit(m,:)>0),2);
                        CY5_AF594_TMR(counter9,2) = transc_num; %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        CY5_AF594_TMR_tp{counter11}(counter10,2) = transc_num; %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        %                         CY5_AF594_TMR(counter9,2) = nansum([nansum(PARaf594_mid.cytoRNA(m,:)),nansum(PARaf594_mid.nucRNA(m))]); %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        %                         CY5_AF594_TMR_tp{counter11}(counter10,2) = nansum([nansum(PARaf594_mid.cytoRNA(m,:)),nansum(PARaf594_mid.nucRNA(m))]); %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        if AF594_hascloud & not(isnan(PARaf594_mid.TotExpIntCloud1(m,1))) & not(isnan(PARaf594_mid.TotExpIntCloud1(m,2)))
                            CY5_AF594_TMR(counter9,2) = CY5_AF594_TMR(counter9,2)-2 + PARaf594_mid.TotExpInt(m,1)/CY5_cloudInt(i)/cloud_correction+ PARaf594_mid.TotExpInt(m,2)/CY5_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR_tp{counter11}(counter10,2) = CY5_AF594_TMR_tp{counter11}(counter10,2) - 2 + PARaf594_mid.TotExpInt(m,1)/CY5_cloudInt(i)/cloud_correction + PARaf594_mid.TotExpInt(m,2)/CY5_cloudInt(i)/cloud_correction;
                        elseif AF594_hascloud & not(isnan(PARaf594_mid.TotExpIntCloud1(m,1)))
                            CY5_AF594_TMR(counter9,2) = CY5_AF594_TMR(counter9,2)-1 + PARaf594_mid.TotExpInt(m,1)/CY5_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR_tp{counter11}(counter10,2) = CY5_AF594_TMR_tp{counter11}(counter10,2) - 1 + PARaf594_mid.TotExpInt(m,1)/CY5_cloudInt(i)/cloud_correction;
                        end
                        CY5_AF594_TMR(counter9,4) = times1(counter11); %stores timepoint
                        CY5_AF594_TMR(counter9,5) = str2num(dates0{i})*.000001; %stores date of experiment
                        CY5_AF594_TMR(counter9,6) = k; %stores image number
                        counter9 = counter9+1;      %increase counter for cell (all)
                        counter10 = counter10+1;    %increase counter for cell (resets each timepoint)
                    end
                    counter9 = counter9_int;
                    counter10 = counter10_int;
                end
                %filename = ['mRNA_' dates0{i} '_STL1-TMR-th27_AF594-GPP2-th17_CTT1-GPD1-th21_' times0{i}{j} '_im' num2str(k)];
                
                if not(isempty(TMR_name))
                    filename = ['mRNA_' dates0{i} '_' strain TMR_name{i} num2str(thTMRall{i}(thNum_TMR)) '_' num2str(thAF594all{i}(thNum_TMR)) CY5_name{i} num2str(thCY5all{i}(thNum_TMR)) '_' times0{i}{j} '_im' num2str(k) '_non_stringent_cloud.mat'];
%                     filename = ['mRNA_' dates0{i} '_' strain TMR_name{i} num2str(thTMRall{i}(thNum_TMR)) '_' num2str(thAF594all{i}(thNum_TMR)) CY5_name{i} num2str(thCY5all{i}(thNum_TMR)) '_' times0{i}{j} '_im' num2str(k)];
                    % AF594_name{i} num2str(thAF594all{i}(thNum_TMR))
                    load(filename,['PARtmr_' nucTH]);
                    for m = 1: size(PARtmr_mid.TotExpInt,1) %Go through every cell in this image
                        if CY5_hascloud | TMR_hascloud
                            if perform_cell_functions
                            %%%Determine cell boundaries
                            A1 = size(cells,2);                            %limits of the segmented image
                            A2 = size(cells,1);                             %limits of the segmented image
                            k1 = Lab == m ; %== j;%sieve out dots in cell j.
                            %figure; imshow(k1,[0 1]);
                            k1=uint16(k1);
                            k3=regionprops(k1,'BoundingBox','Area');
                            k4 = k3.BoundingBox; %create the rectangular box around the cell.
                            X0=round(k4(1))-4;
                            Y0=round(k4(2))-4;
                            X1=round(k4(1)+ k4(3))+4;
                            Y1=round(k4(2)+ k4(4))+4;
                            if X0 < 1;  X0 = 1; end
                            if Y0 < 1;  Y0 = 1; end
                            if X1 > A1; X1 = A1; end
                            if Y1 > A2; Y1 = A2; end
                            if CY5_hascloud
                                CY5_cloud1_cell = PARcy5_mid.clouds1(Y0:Y1,X0:X1,:);      %Generate just CY5 cloud 1 for this cell
                                CY5_cloud2_cell = PARcy5_mid.clouds2(Y0:Y1,X0:X1,:);      %Generate CY5 cloud 2 for this cell
                            end
                            %                                 if sum(CY5_cloud1_cell(:)) > 50                           %to check how mnay clouds in image
                            %                                     ['cell ' num2str(m) ' has a cloud']
                            %                                     figure; clf; imshow(max(CY5_cloud1_cell,[],3),[])
                            %                                 end
                            %                             end
                            %                         end
                            %
                            if not(isempty(TMR_name)) & TMR_hascloud
                                TMR_cloud1_cell = PARtmr_mid.clouds1(Y0:Y1,X0:X1,:);      %Generate TMR cloud 1 for this cell
                                TMR_cloud2_cell = PARtmr_mid.clouds2(Y0:Y1,X0:X1,:);      %Generate TMR cloud 2 for this cell
                            end
                            end
                        end
                        %%%Eliminate dead pixels  (needs to be improved)
                        for p = 1:size(PARtmr_mid.xinit,2)
                            same_pos = find(PARtmr_mid.xinit(m,:) == PARtmr_mid.xinit(m,p));
                            if size(same_pos,2)>10
                                for q = 1:size(same_pos,2)
                                    PARtmr_mid.xinit(m,same_pos(q)) = NaN;
                                    PARtmr_mid.yinit(m,same_pos(q)) = NaN;
                                    PARtmr_mid.zinit(m,same_pos(q)) = NaN;
                                    PARtmr_mid.TotExpInt(m,same_pos(q)) = NaN;
                                end
                            end
                        end
                        %%% Below removes duplicates of the same
                        %%% transcript. Goes through every transcript,
                        %%% finds distance between it and other
                        %%% transcripts, and eliminates others
                        for p = 1:size(PARtmr_mid.xinit,2)
                            dists = zeros(1,size(size(PARtmr_mid.xinit,2),3)); %stores distance from transcript to others
                            if not(isnan(PARtmr_mid.xinit(m,p)))
                                for pp = 1:size(PARtmr_mid.xinit,2)    %find distance squared
                                    dists(1,pp,1) = (PARtmr_mid.xinit(m,p)- PARtmr_mid.xinit(m,pp))^2;
                                    dists(1,pp,2) = (PARtmr_mid.yinit(m,p)- PARtmr_mid.yinit(m,pp))^2;
                                    dists(1,pp,3) = ((PARtmr_mid.zinit(m,p)- PARtmr_mid.zinit(m,pp))*3)^2;
                                end
                                distF = sqrt(dists(1,:,1)+dists(1,:,2)+dists(1,:,3));
                                for pp = 1:size(PARtmr_mid.xinit,2)
                                    if p ~= pp & distF(pp) < same_dist
                                        PARtmr_mid.xinit(m,pp) = NaN;
                                        PARtmr_mid.yinit(m,pp) = NaN;
                                        PARtmr_mid.zinit(m,pp) = NaN;
                                        PARtmr_mid.TotExpInt(m,p) = PARtmr_mid.TotExpInt(m,p)+PARtmr_mid.TotExpInt(m,pp);
                                        PARtmr_mid.TotExpInt(m,pp) = NaN;
                                    end
                                end
                            end
                        end
                         %% Determine distance to edge of nucleus
                        if perform_nuc_functions
                         NucBorder1 = bwmorph3(Label_mid_cell,'remove');
                        %                         figure(300); clf;
                        %                         imshow(NucBorder1(:,:,round(size(NucBorder1,3)/2)),[]) %Visualize middle slice
                        %                         for zzz = 1:size(NucBorder1,3)    %All slices
                        %                             figure(300); clf; imshow(NucBorder1(:,:,zzz),[])
                        %                             pause(.5)
                        %                         end
                        lin_Nuc = find(NucBorder1 == 1);[xNuc,yNuc,zNuc] = ind2sub(size(NucBorder1),lin_Nuc(1:4:end)); clear lin_Nuc; %These are vertical vectors (x differs)
                        pNuc1 = xNuc; pNuc1(:,:,2) = yNuc; pNuc1(:,:,3) = zNuc; clear xNuc yNuc zNuc NucBorder1; %combine matrices
                        pNuc1 = int16(pNuc1); 
                        pNuc2 = repmat(pNuc1,[1,size(PARtmr_mid.xinit,2),1]);  %replicate in y axis the matricies
                        %pNuc2 = int16(pNuc2);
                        pTran = PARtmr_mid.xinit(m,1:size(PARtmr_mid.xinit,2));
                        pTran(:,:,2) = PARtmr_mid.yinit(m,1:size(PARtmr_mid.yinit,2));
                        pTran(:,:,3) = PARtmr_mid.zinit(m,1:size(PARtmr_mid.zinit,2));
                        pTran = int16(pTran);
                        pTran1 = repmat(pTran,[size(pNuc2,1),1,1]);
                        pDiff = pNuc2-pTran1; clear pNuc2 pTran1; pDiff = int32(pDiff);
                        pDiff(:,:,3) = pDiff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
                        pDiff = sum(pDiff.^2,3); %pDiff = pDiff.^.5; %To keep in in integer form, do not square root
                        TMR_Nuc_dist(counter9,1:size(PARtmr_mid.xinit,2)) = min(pDiff,[],1);
                        TMR_Nuc_dist_tp{counter11}(counter10,1:size(PARtmr_mid.xinit,2)) = min(pDiff,[],1);
                        clear pDiff
                        end
                        %%
                        %%%store the positions of all transcripts
                        TMR_pos_init(counter9,1:size(PARtmr_mid.xinit,2),1) = PARtmr_mid.xinit(m,1:size(PARtmr_mid.xinit,2));
                        TMR_pos_init(counter9,1:size(PARtmr_mid.xinit,2),2) = PARtmr_mid.yinit(m,1:size(PARtmr_mid.yinit,2));
                        TMR_pos_init(counter9,1:size(PARtmr_mid.xinit,2),3) = PARtmr_mid.zinit(m,1:size(PARtmr_mid.zinit,2));
                        TMR_pos_init_tp{counter11}(counter10,1:size(PARtmr_mid.xinit,2),1) = PARtmr_mid.xinit(m,1:size(PARtmr_mid.xinit,2));
                        TMR_pos_init_tp{counter11}(counter10,1:size(PARtmr_mid.xinit,2),2) = PARtmr_mid.yinit(m,1:size(PARtmr_mid.yinit,2));
                        TMR_pos_init_tp{counter11}(counter10,1:size(PARtmr_mid.xinit,2),3) = PARtmr_mid.zinit(m,1:size(PARtmr_mid.zinit,2));
                        %%% Determine transcript number
                        transc_num = size(find(PARtmr_mid.xinit(m,:)>0),2);
                        CY5_AF594_TMR(counter9,3) = transc_num; %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        CY5_AF594_TMR_tp{counter11}(counter10,3) = transc_num; %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        TMR_int(counter9,3:size(PARtmr_mid.TotExpInt,2)) = PARtmr_mid.TotExpInt(m,3:size(PARtmr_mid.TotExpInt,2));  %Store total experimental intensity
                        TMR_int_tp{counter11}(counter10,3:size(PARtmr_mid.TotExpInt,2)) = PARtmr_mid.TotExpInt(m,3:size(PARtmr_mid.TotExpInt,2)); %Store total experimental intensity
                        CY5_AF594_TMR(counter9,13) = PARtmr_mid.TotExpInt(m,1)/TMR_cloudInt(i);
                        CY5_AF594_TMR(counter9,14) = PARtmr_mid.TotExpInt(m,2)/TMR_cloudInt(i);
                        CY5_AF594_TMR_tp{counter11}(counter10,13) = CY5_AF594_TMR(counter9,13);
                        CY5_AF594_TMR_tp{counter11}(counter10,14) = CY5_AF594_TMR(counter9,14);
                        if TMR_hascloud
                            TMR_int_tp{counter11}(counter10,1:2) = PARtmr_mid.TotExpInt(m,1:2)/cloud_correction; %Store total experimental intensity
                            TMR_int(counter9,1:2) = PARtmr_mid.TotExpInt(m,1:2)/cloud_correction; %Store total experimental intensity
                        else
                            TMR_int_tp{counter11}(counter10,1:2) = PARtmr_mid.TotExpInt(m,1:2); %Store total experimental intensity
                            TMR_int(counter9,1:2) = PARtmr_mid.TotExpInt(m,1:2); %Store total experimental intensity
                        end
                        %CY5_AF594_TMR(counter9,3) = nansum([nansum(PARtmr_mid.cytoRNA(m,:)),nansum(PARtmr_mid.nucRNA(m))]); %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        %                         CY5_AF594_TMR_tp{counter11}(counter10,3) = nansum([nansum(PARtmr_mid.cytoRNA(m,:)),nansum(PARtmr_mid.nucRNA(m))]); %Set Jpx value to nuclear + Cytoplasmic RNA molecules
                        %%% Determine whether transcripts are inside clouds
                        if perform_cell_functions
                        if CY5_hascloud
                            for temp = 1:size(PARtmr_mid.xinit,2)
                                trans_loc = [TMR_pos_init(counter9,temp,1) TMR_pos_init(counter9,temp,2) TMR_pos_init(counter9,temp,3)]; %temporary x,y,z positions
                                if not(isnan(trans_loc(1))) %Check that it is not NaN
                                    part_c1 = CY5_cloud1_cell(trans_loc(1),trans_loc(2),trans_loc(3));  %This needs to be adjusted. clouds1 refers to entire image while locations are for just the cell
                                    part_c2 = CY5_cloud2_cell(trans_loc(1),trans_loc(2),trans_loc(3));
                                    if part_c1 == 1 & part_c2 == 1
                                        TMR_in_CY5cloud(counter9,temp) = 4;
                                        TMR_in_CY5cloud_tp{counter11}(counter10,temp) = 4;
                                    elseif part_c1 == 1
                                        TMR_in_CY5cloud(counter9,temp) = 2;
                                        TMR_in_CY5cloud_tp{counter11}(counter10,temp) = 2;
                                    elseif part_c2 == 1
                                        TMR_in_CY5cloud(counter9,temp) = 3;
                                        TMR_in_CY5cloud_tp{counter11}(counter10,temp) = 3;
                                    else
                                        TMR_in_CY5cloud(counter9,temp) = 1;
                                        TMR_in_CY5cloud_tp{counter11}(counter10,temp) = 1;
                                    end
                                else
                                    TMR_in_CY5cloud(counter9,temp) = 0;
                                    TMR_in_CY5cloud_tp{counter11}(counter10,temp) = 0;
                                end
                            end
                            %% Determine distances to the clouds 
                            pTran = PARtmr_mid.xinit(m,1:size(PARtmr_mid.xinit,2));
                            pTran(:,:,2) = PARtmr_mid.yinit(m,1:size(PARtmr_mid.yinit,2));
                            pTran(:,:,3) = PARtmr_mid.zinit(m,1:size(PARtmr_mid.zinit,2));
                            pTran = int16(pTran);
                            if not(isnan(PARcy5_mid.xinit(m,1)))
                                C1_bord = bwmorph3(CY5_cloud1_cell,'remove');
                                %                                 figure(301); clf;
                                %                                 imshow(C1_bord(:,:,round(size(C1_bord,3)/2)),[]) %Visualize middle slice
                                %                                 for zzz = 1:size(C1_bord,3)    %All slices
                                %                                     figure(301); clf; imshow(C1_bord(:,:,zzz),[])
                                %                                     title(['slice ' num2str(zzz)])
                                %                                     pause(.5)
                                %                                 end
                                lin_C1 = find(C1_bord == 1); [xC1,yC1,zC1] =  ind2sub(size(C1_bord),lin_C1(1:2:end)); clear lin_C1; %These are vertical vectors (x differs)
%                                 xC1 = xC1'; yC1 = yC1'; zC1 = zC1';     %Change x's to y's
                                pC1 = xC1; pC1(:,:,2) = yC1; pC1(:,:,3) = zC1; clear xC1 yC1 zC1 %populate one matrix
                                pC1 = int16(pC1);
                                %%% fast but memory-intensive way of determining distances
                                if size(pC1,1) < 1
                                    C1_TMR_dist(counter9,1:size(PARtmr_mid.xinit,2)) = nan(1,size(PARtmr_mid.xinit,2));
                                    C1_TMR_dist_tp{counter11}(counter10,1:size(PARtmr_mid.xinit,2)) = nan(1,size(PARtmr_mid.xinit,2));
                                else
                                    pTran1 = repmat(pTran,[size(pC1,1),1,1]);
                                    pC1_2 = repmat(pC1,[1,size(pTran1,2),1]); %make copies for comparison. Commented outfor memory
                                    pDiff = pTran1-pC1_2; clear pC1_2 pTran1; pDiff = int32(pDiff);
                                    pDiff(:,:,3) = pDiff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
                                    pDiff = sum(pDiff.^2,3); %pDiff = pDiff.^.5;
                                    C1_TMR_dist(counter9,1:size(PARtmr_mid.xinit,2)) = min(pDiff,[],1);
                                    C1_TMR_dist_tp{counter11}(counter10,1:size(PARtmr_mid.xinit,2)) = min(pDiff,[],1);
                                    clear pDiff
                                end
                                %%%
                                %%% Slow but memory-efficient way
%                                 if size(pC1,2) < 1
%                                     Cl_Nuc_dist(counter9,2) = NaN;
%                                     Cl_Nuc_dist_tp{counter11}(counter10,2) = NaN;
%                                 else
%                                     
%                                     pDiff = zeros(size(pC1,2),1);
%                                     for cl_dim = 1:size(pC1,2)  %less memory-intensive, but slower way for determining distances
%                                         pC1_temp = pC1(:,cl_dim,:);
%                                         pC1_temp = repmat(pC1_temp,[size(pNuc1,1),1,1]); %make copies for comparison.
%                                         temp_Diff = pC1_temp-pNuc1;
%                                         temp_Diff(:,:,3) = temp_Diff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
%                                         temp_Diff = sum(temp_Diff.^2,3); temp_Diff = temp_Diff.^.5;
%                                         pDiff(cl_dim) = min(temp_Diff(:));
%                                     end
%                                     Cl_Nuc_dist(counter9,1) = min(pDiff(:));
%                                     Cl_Nuc_dist_tp{counter11}(counter10,1) = min(pDiff(:));
%                                 end
                                %%%
                            else                                           %If no cloud, make NaN
                                C1_TMR_dist(counter9,1:size(PARtmr_mid.xinit,2)) = nan(1,size(PARtmr_mid.xinit,2));
                                C1_TMR_dist_tp{counter11}(counter10,1:size(PARtmr_mid.xinit,2)) = nan(1,size(PARtmr_mid.xinit,2));
                            end
                            if not(isnan(PARcy5_mid.xinit(m,1)))
                                C1_bord = bwmorph3(CY5_cloud2_cell,'remove');
                                %                                 figure(301); clf;
                                %                                 imshow(C1_bord(:,:,round(size(C1_bord,3)/2)),[]) %Visualize middle slice
                                %                                 for zzz = 1:size(C1_bord,3)    %All slices
                                %                                     figure(301); clf; imshow(C1_bord(:,:,zzz),[])
                                %                                     title(['slice ' num2str(zzz)])
                                %                                     pause(.5)
                                %                                 end
                                lin_C1 = find(C1_bord == 1); [xC1,yC1,zC1] =  ind2sub(size(C1_bord),lin_C1(1:2:end)); clear lin_C1; %These are vertical vectors (x differs)
%                                 xC1 = xC1'; yC1 = yC1'; zC1 = zC1';     %Change x's to y's
                                pC1 = xC1; pC1(:,:,2) = yC1; pC1(:,:,3) = zC1; clear xC1 yC1 zC1 %populate one matrix
                                pC1 = int16(pC1);
                                %%% fast but memory-intensive way of determining distances
                                if size(pC1,1) < 1
                                    C2_TMR_dist(counter9,1:size(PARtmr_mid.xinit,2)) = nan(1,size(PARtmr_mid.xinit,2));
                                    C2_TMR_dist_tp{counter11}(counter10,1:size(PARtmr_mid.xinit,2)) = nan(1,size(PARtmr_mid.xinit,2));
                                else
                                    pTran1 = repmat(pTran,[size(pC1,1),1,1]);
                                    pC1_2 = repmat(pC1,[1,size(pTran1,2),1]); %make copies for comparison. Commented outfor memory
                                    pDiff = pTran1-pC1_2; clear pC1_2 pTran1; pDiff = int32(pDiff);
                                    pDiff(:,:,3) = pDiff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
                                    pDiff = sum(pDiff.^2,3); %pDiff = pDiff.^.5;
                                    C2_TMR_dist(counter9,1:size(PARtmr_mid.xinit,2)) = min(pDiff,[],1);
                                    C2_TMR_dist_tp{counter11}(counter10,1:size(PARtmr_mid.xinit,2)) = min(pDiff,[],1);
                                    clear pDiff
                                end
                                %%%
                                %%% Slow but memory-efficient way
%                                 if size(pC1,2) < 1
%                                     Cl_Nuc_dist(counter9,2) = NaN;
%                                     Cl_Nuc_dist_tp{counter11}(counter10,2) = NaN;
%                                 else
%                                     
%                                     pDiff = zeros(size(pC1,2),1);
%                                     for cl_dim = 1:size(pC1,2)  %less memory-intensive, but slower way for determining distances
%                                         pC1_temp = pC1(:,cl_dim,:);
%                                         pC1_temp = repmat(pC1_temp,[size(pNuc1,1),1,1]); %make copies for comparison.
%                                         temp_Diff = pC1_temp-pNuc1;
%                                         temp_Diff(:,:,3) = temp_Diff(:,:,3)*z_adjust;   %Adjust distance determined by z axis based on how large
%                                         temp_Diff = sum(temp_Diff.^2,3); temp_Diff = temp_Diff.^.5;
%                                         pDiff(cl_dim) = min(temp_Diff(:));
%                                     end
%                                     Cl_Nuc_dist(counter9,1) = min(pDiff(:));
%                                     Cl_Nuc_dist_tp{counter11}(counter10,1) = min(pDiff(:));
%                                 end
                                %%%
                            else                                           %If no cloud, make NaN
                                C2_TMR_dist(counter9,1:size(PARtmr_mid.xinit,2)) = nan(1,size(PARtmr_mid.xinit,2));
                                C2_TMR_dist_tp{counter11}(counter10,1:size(PARtmr_mid.xinit,2)) = nan(1,size(PARtmr_mid.xinit,2));
                            end
                        end
                        if TMR_hascloud
                            for temp = 1:size(PARtmr_mid.xinit,2)
                                trans_loc = [TMR_pos_init(counter9,temp,1) TMR_pos_init(counter9,temp,2) TMR_pos_init(counter9,temp,3)]; %temporary x,y,z positions
                                if not(isnan(trans_loc(1))) %Check that it is not NaN
                                    part_c1 = TMR_cloud1_cell(trans_loc(1),trans_loc(2),trans_loc(3));  %This needs to be adjusted. clouds1 refers to entire image while locations are for just the cell
                                    part_c2 = TMR_cloud2_cell(trans_loc(1),trans_loc(2),trans_loc(3));
                                    if part_c1 == 1 & part_c2 == 1
                                        TMR_in_TMRcloud(counter9,temp) = 4;
                                        TMR_in_TMRcloud_tp{counter11}(counter10,temp) = 4;
                                    elseif part_c1 == 1
                                        TMR_in_TMRcloud(counter9,temp) = 2;
                                        TMR_in_TMRcloud_tp{counter11}(counter10,temp) = 2;
                                    elseif part_c2 == 1
                                        TMR_in_TMRcloud(counter9,temp) = 3;
                                        TMR_in_TMRcloud_tp{counter11}(counter10,temp) = 3;
                                    else
                                        TMR_in_TMRcloud(counter9,temp) = 1;
                                        TMR_in_TMRcloud_tp{counter11}(counter10,temp) = 1;
                                    end
                                else
                                    TMR_in_TMRcloud(counter9,temp) = 0;
                                    TMR_in_TMRcloud_tp{counter11}(counter10,temp) = 0;
                                end
                            end
                        end
                        end
                        %%%
                        if TMR_hascloud & not(isnan(PARtmr_mid.TotExpInt(m,1))) & not(isnan(PARtmr_mid.TotExpInt(m,2)))
                            CY5_AF594_TMR(counter9,3) = CY5_AF594_TMR(counter9,3)-2 + PARtmr_mid.TotExpInt(m,1)/TMR_cloudInt(i)/cloud_correction+ PARtmr_mid.TotExpInt(m,2)/TMR_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR_tp{counter11}(counter10,3) = CY5_AF594_TMR_tp{counter11}(counter10,3) - 2 + PARtmr_mid.TotExpInt(m,1)/TMR_cloudInt(i)/cloud_correction + PARtmr_mid.TotExpInt(m,2)/TMR_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR(counter9,10) = sum(PARtmr_mid.TotExpInt(m,1:2) > TMR_tran_thres*TMR_cloudInt(i)*cloud_correction)+sum(PARtmr_mid.TotExpInt(m,3:size(PARtmr_mid.TotExpInt,2)) > TMR_tran_thres*TMR_cloudInt(i));        %Number of transcription sites/clouds in the cell
                            CY5_AF594_TMR_tp{counter11}(counter10,10) =CY5_AF594_TMR(counter9,10);      %Number of transcription sites/clouds in the cell
                            CY5_AF594_TMR(counter9,13) = PARtmr_mid.TotExpInt(m,1)/TMR_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR(counter9,14) = PARtmr_mid.TotExpInt(m,2)/TMR_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR_tp{counter11}(counter10,13) = CY5_AF594_TMR(counter9,13);
                            CY5_AF594_TMR_tp{counter11}(counter10,14) = CY5_AF594_TMR(counter9,14);
                            %%%Below adds position information for the
                            %%%transcription sites                            
                            counter20 = 1;  %Counter for transcription sites
                            Tsite_ind = zeros(1,4); %Will have index for transcription sites                            
                            if CY5_AF594_TMR(counter9,11) >= TMR_tran_thres
                                Tsite_ind(counter20) = 1;
                                counter20 = counter20+1;
                            end
                            if CY5_AF594_TMR(counter9,12) >= TMR_tran_thres
                                Tsite_ind(counter20) = 2;
                                counter20 = counter20+1;
                            end
                            [temp_int,int_ind] = sort(TMR_int(counter9,3:size(PARtmr_mid.TotExpInt,2)),'descend');   %Sort the transcript intensities to find the transcription sites
                            int_ind = int_ind+2;    %Adjust for leaving out first two
                            Poss_tsites = find(temp_int >= TMR_tran_thres*TMR_cloudInt(i));
                            Tsites_left = 5-counter20;
                            if not(isempty(Poss_tsites)) & Tsites_left > 0
                                Tsites_add = min([size(Poss_tsites,2) Tsites_left]);    %See how many transcription sites to add
                                Tsite_ind(counter20:counter20+Tsites_add-1) = int_ind(Poss_tsites(1:Tsites_add)); %Add the indexes of the transcription sites
                            end
                            for index1 = 1:4
                                if Tsite_ind(index1) >0
                                    CY5_AF594_TMR(counter9,(32:34)+3*(index1-1)) = TMR_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                                    CY5_AF594_TMR_tp{counter11}(counter10,(32:34)+3*(index1-1)) = TMR_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                                end
                            end
                            %%%
                        elseif TMR_hascloud & not(isnan(PARtmr_mid.TotExpIntCloud1(m,1)))
                            CY5_AF594_TMR(counter9,3) = CY5_AF594_TMR(counter9,3)-1 + PARtmr_mid.TotExpInt(m,1)/TMR_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR_tp{counter11}(counter10,3) = CY5_AF594_TMR_tp{counter11}(counter10,3) - 1 + PARtmr_mid.TotExpInt(m,1)/TMR_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR(counter9,10) = sum(PARtmr_mid.TotExpInt(m,1) > TMR_tran_thres*TMR_cloudInt(i)*cloud_correction)+sum(PARtmr_mid.TotExpInt(m,3:size(PARtmr_mid.TotExpInt,2)) > TMR_tran_thres*TMR_cloudInt(i));        %Number of transcription sites/clouds in the cell
                            CY5_AF594_TMR_tp{counter11}(counter10,10) =CY5_AF594_TMR(counter9,10);      %Number of transcription sites/clouds in the cell
                            CY5_AF594_TMR(counter9,13) = PARtmr_mid.TotExpInt(m,1)/TMR_cloudInt(i)/cloud_correction;
                            CY5_AF594_TMR_tp{counter11}(counter10,13) = CY5_AF594_TMR(counter9,13);
                            %%%Below adds position information for the
                            %%%transcription sites   
                            counter20 = 1;  %Counter for transcription sites
                            Tsite_ind = zeros(1,4); %Will have index for transcription sites                            
                            if CY5_AF594_TMR(counter9,11) >= TMR_tran_thres
                                Tsite_ind(counter20) = 1;
                                counter20 = counter20+1;
                            end
                            if CY5_AF594_TMR(counter9,12) >= TMR_tran_thres
                                Tsite_ind(counter20) = 2;
                                counter20 = counter20+1;
                            end
                            [temp_int,int_ind] = sort(TMR_int(counter9,3:size(PARtmr_mid.TotExpInt,2)),'descend');   %Sort the transcript intensities to find the transcription sites
                            int_ind = int_ind+2;    %Adjust for leaving out first two
                            Poss_tsites = find(temp_int >= TMR_tran_thres*TMR_cloudInt(i));
                            Tsites_left = 5-counter20;
                            if not(isempty(Poss_tsites)) & Tsites_left > 0
                                Tsites_add = min([size(Poss_tsites,2) Tsites_left]);    %See how many transcription sites to add
                                Tsite_ind(counter20:counter20+Tsites_add-1) = int_ind(Poss_tsites(1:Tsites_add)); %Add the indexes of the transcription sites
                            end
                            for index1 = 1:4
                                if Tsite_ind(index1) >0
                                    CY5_AF594_TMR(counter9,(32:34)+3*(index1-1)) = TMR_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                                    CY5_AF594_TMR_tp{counter11}(counter10,(32:34)+3*(index1-1)) = TMR_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                                end
                            end
                            %%%
                        else
                            CY5_AF594_TMR(counter9,10) = sum(PARtmr_mid.TotExpInt(m,:) > TMR_tran_thres*TMR_cloudInt(i));        %Number of transcription sites/clouds in the cell
                            CY5_AF594_TMR_tp{counter11}(counter10,10) =  CY5_AF594_TMR(counter9,10);        %Number of transcription sites/clouds in the cell
                            %%%Below adds position information for the
                            %%%transcription sites
                            Tsite_ind = zeros(1,4); %Will have index for transcription sites                            
                            [temp_int,int_ind] = sort(TMR_int(counter9,1:size(PARtmr_mid.TotExpInt,2)),'descend');   %Sort the transcript intensities to find the transcription sites
                            Poss_tsites = find(temp_int >= TMR_tran_thres*TMR_cloudInt(i));
                            Tsites_left = 4;
                            if not(isempty(Poss_tsites)) & Tsites_left > 0
                                Tsites_add = min([size(Poss_tsites,2) Tsites_left]);    %See how many transcription sites to add
                                Tsite_ind(1:Tsites_add) = int_ind(Poss_tsites(1:Tsites_add)); %Add the indexes of the transcription sites
                            end
                            for index1 = 1:4
                                if Tsite_ind(index1) >0
                                    CY5_AF594_TMR(counter9,(32:34)+3*(index1-1)) = TMR_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                                    CY5_AF594_TMR_tp{counter11}(counter10,(32:34)+3*(index1-1)) = TMR_pos_init(counter9,Tsite_ind(index1),:);  %Add position information for first transcription site
                                end
                            end
                            %%%
                        end
                        CY5_AF594_TMR(counter9,4) = times1(counter11);  %stores timepoint
                        CY5_AF594_TMR(counter9,5) = str2num(dates0{i})*.000001; %stores date of experiment
                        CY5_AF594_TMR(counter9,6) = k; %stores image number
                        %%% Determine if there is a colocalized spot in CY5
                        %%% for TMR spots
                        if perform_col_functions
                        for p = 1:size(PARtmr_mid.xinit,2)
                            dists = zeros(1,size(size(PARcy5_mid.xinit,2),3)); %stores distance from transcript to others
                            if not(isnan(PARtmr_mid.xinit(m,p)))
                                for pp = 1:size(PARcy5_mid.xinit,2)    %find distance squared
                                    dists(1,pp,1) = (PARtmr_mid.xinit(m,p)- PARcy5_mid.xinit(m,pp))^2;
                                    dists(1,pp,2) = (PARtmr_mid.yinit(m,p)- PARcy5_mid.yinit(m,pp))^2;
                                    dists(1,pp,3) = ((PARtmr_mid.zinit(m,p)- PARcy5_mid.zinit(m,pp))*z_adjust)^2;
                                end
                                distF = sqrt(dists(1,:,1)+dists(1,:,2)+dists(1,:,3));
                                if size(find(distF < dist_col),2) > 0
                                    T_C_col(counter9,p) = 2;
                                    T_C_col_tp{counter11}(counter10,p) = 2;
                                    T_C_col_ind(counter9,p) = find(distF < dist_col,1,'last');
                                    T_C_col_ind_tp{counter11}(counter10,p) = find(distF < dist_col,1,'last');
                                else
                                    T_C_col(counter9,p) = 1;
                                    T_C_col_tp{counter11}(counter10,p) = 1;
                                end
                            else
                                T_C_col(counter9,p) = 0;
                                T_C_col_tp{counter11}(counter10,p) = 0;
                            end
                        end
                        %%% Determine if there is a colocalized spot in TMR
                        %%% for TMR spots (control)
                        for p = 1:size(PARtmr_mid.xinit,2)
                            dists = zeros(1,size(size(PARtmr_mid.xinit,2),3)); %stores distance from transcript to others
                            if not(isnan(PARtmr_mid.xinit(m,p)))
                                for pp = 1:size(PARtmr_mid.xinit,2)    %find distance squared
                                    dists(1,pp,1) = (PARtmr_mid.xinit(m,p)- PARtmr_mid.xinit(m,pp))^2;
                                    dists(1,pp,2) = (PARtmr_mid.yinit(m,p)- PARtmr_mid.yinit(m,pp))^2;
                                    dists(1,pp,3) = ((PARtmr_mid.zinit(m,p)- PARtmr_mid.zinit(m,pp))*z_adjust)^2;
                                end
                                distF = sqrt(dists(1,:,1)+dists(1,:,2)+dists(1,:,3));
                                if size(find(distF < dist_col),2) > 0
                                    T_T_col(counter9,p) = 2;
                                    T_T_col_tp{counter11}(counter10,p) = 2;
                                else
                                    T_T_col(counter9,p) = 1;
                                    T_T_col_tp{counter11}(counter10,p) = 1;
                                end
                            else
                                T_T_col(counter9,p) = 0;
                                T_T_col_tp{counter11}(counter10,p) = 0;
                            end
                        end
                        %%%
                                          T_C_col_tp{counter11}(T_C_col_tp{counter11} == 0) = NaN;
                        end
                          counter9 = counter9+1;      %increase counter for cell (all)
                        counter10 = counter10+1;    %increase counter for cell (resets each timepoint)
                    end
                    counter9 = counter9_int;
                    counter10 = counter10_int;
                    %%% Determine if there is a colocalized spot in CY5
                    %%% for TMR spots
                    if perform_col_functions                    
                    if not(isempty(CY5_name))
                        counter9 = counter9_int;
                        counter10 = counter10_int;
                        for m = 1: size(PARcy5_mid.TotExpInt,1)
                            for p = 1:size(PARcy5_mid.xinit,2)
                                dists = zeros(1,size(size(PARcy5_mid.xinit,2),3)); %stores distance from transcript to others
                                if not(isnan(PARcy5_mid.xinit(m,p)))
                                    for pp = 1:size(PARtmr_mid.xinit,2)    %find distance squared
                                        dists(1,pp,1) = (PARcy5_mid.xinit(m,p)- PARtmr_mid.xinit(m,pp))^2;
                                        dists(1,pp,2) = (PARcy5_mid.yinit(m,p)- PARtmr_mid.yinit(m,pp))^2;
                                        dists(1,pp,3) = ((PARcy5_mid.zinit(m,p)- PARtmr_mid.zinit(m,pp))*z_adjust)^2;
                                    end
                                    distF = sqrt(dists(1,:,1)+dists(1,:,2)+dists(1,:,3));
                                    if size(find(distF < dist_col),2) > 0
                                        C_T_col(counter9,p) = 2;
                                        C_T_col_tp{counter11}(counter10,p) = 2;
                                        C_T_col_ind(counter9,p) = find(distF < dist_col,1,'last');
                                        C_T_col_ind_tp{counter11}(counter10,p) = find(distF < dist_col,1,'last');
                                    else
                                        C_T_col(counter9,p) = 1;
                                        C_T_col_tp{counter11}(counter10,p) = 1;
                                    end
                                else
                                    C_T_col(counter9,p) = 0;
                                    C_T_col_tp{counter11}(counter10,p) = 0;
                                end
                            end
                            counter9 = counter9+1;      %increase counter for cell (all)
                            counter10 = counter10+1;    %increase counter for cell (resets each timepoint)
                        end
                    end
                    end
                end                   
            catch
                fail_count = fail_count+1;
                ['processing ' filename ' failed']
            end
        end
        temp_tpname = times0{i}{j};
    end
end
T_C_col(T_C_col == 0) = NaN;
T_T_col(T_T_col == 0) = NaN;
C_T_col(T_C_col == 0) = NaN;

%% Analysis

% figure(20); clf; gscatter(CY5_AF594_TMR(:,2),CY5_AF594_TMR(:,1),CY5_AF594_TMR(:,4),[],'o',10);
% xlabel('Tsix','Fontsize',16)
% ylabel('Xist','Fontsize',16)
% legend('0 day','2 day','3 day','4 day','5 day','6 day','11 day','13 day','FontSize',20)
%
% figure(40); clf; gscatter(CY5_AF594_TMR(:,2),CY5_AF594_TMR(:,3),CY5_AF594_TMR(:,4),[],'o',10);
% xlabel('Tsix','Fontsize',16)
% ylabel('Jpx','Fontsize',16)
% legend('0 day','2 day','3 day','4 day','5 day','6 day','11 day','13 day','FontSize',20)
%
% figure(50); clf; gscatter(CY5_AF594_TMR(:,3),CY5_AF594_TMR(:,1),CY5_AF594_TMR(:,4),[],'o',10);
% xlabel('Jpx','Fontsize',16)
% ylabel('Xist','Fontsize',16)
% legend('0 day','2 day','3 day','4 day','5 day','6 day','11 day','13 day','FontSize',20)
% %Populate distributions for Jpx

% dist0 = zeros(1,1); %Xist Distribution with zero Jpx
% dist1 = zeros(1,1);  %Xist Distribution with one Jpx
% dist2 = zeros(1,1);  %Xist distribution with 2 Jpx
% dist3 = zeros(1,1);
% dist4 = zeros(1,1);
% dist5 = zeros(1,1);
% counter_dist0 = 1;
% counter_dist1 = 1;
% counter_dist2 = 1;
% counter_dist3 = 1;
% counter_dist4 = 1;
% counter_dist5 = 1;
% CY5_AF594_TMR2 = CY5_AF594_TMR;       %will be used for boxplot (4th column has which category for grouping)
% for i = 1:size(CY5_AF594_TMR,1)
%     if CY5_AF594_TMR(i,2) == 0
%         dist0(counter_dist0,1) = CY5_AF594_TMR(i,1);
%         counter_dist0 = counter_dist0+1;
%         CY5_AF594_TMR2(i,4) = 0;
%     elseif CY5_AF594_TMR(i,2) == 1
%         dist1(counter_dist1,1) = CY5_AF594_TMR(i,1);
%         counter_dist1 = counter_dist1+1;
%         CY5_AF594_TMR2(i,4) = 1;
%     elseif CY5_AF594_TMR(i,2) == 2
%         dist2(counter_dist2,1) = CY5_AF594_TMR(i,1);
%         counter_dist2 = counter_dist2+1;
%         CY5_AF594_TMR2(i,4) = 2;
%     elseif CY5_AF594_TMR(i,2) == 3
%         dist3(counter_dist3,1) = CY5_AF594_TMR(i,1);
%         counter_dist3 = counter_dist3+1;
%         CY5_AF594_TMR2(i,4) = 3;
%     elseif CY5_AF594_TMR(i,2) >= 4
%         dist4(counter_dist4,1) = CY5_AF594_TMR(i,1);
%         counter_dist4 = counter_dist4+1;
%         CY5_AF594_TMR2(i,4) = 4;
% %     elseif CY5_AF594_TMR(i,2) >= 5
% %         dist5(counter_dist5,1) = CY5_AF594_TMR(i,1);
% %         counter_dist5 = counter_dist5+1;
%     end
% end
% Jpx_labels = {'0 Jpx','1 Jpx','2 Jpx','3 Jpx','>=4 Jpx'};
% %Number of cells with different amounts of Jpx
% Numcell_Jpx = [size(dist0,1), size(dist1,1),size(dist2,1),size(dist3,1),size(dist4,1)];
% figure(21); clf; bar([1:size(Numcell_Jpx,2)],Numcell_Jpx,'r');
% set(gca,'xticklabel',Jpx_labels)
% ylabel('Number of Cells')
%
% %Mean Xist for cells with different numbers of Jpx
% Xist_means = [mean(dist0),mean(dist1),mean(dist2),mean(dist3),mean(dist4),];
% figure(22); clf; bar([1:size(Numcell_Jpx,2)],Xist_means);
% set(gca,'xticklabel',Jpx_labels)
% ylabel('Mean Tsix Molecules')
%
% % boxplot of Tsix for different numbers of Jpx
%
% figure(25); clf;
% boxplot(CY5_AF594_TMR2(:,1),CY5_AF594_TMR2(:,4),'Labels',Jpx_labels)
% %xlabel('Jpx per Cell','Fontsize',10)
% ylabel('Tsix per Cell','Fontsize',10)
% hold on
% plot(0:6,mean(CY5_AF594_TMR2(:,1)):.0001:mean(CY5_AF594_TMR2(:,1))+.0006,'--r')
% 
% %%Distribution of Tsix
% bin_size = (max(CY5_AF594_TMR(:,1))-min(CY5_AF594_TMR(:,1)))/10;
% binsX = 0:bin_size:max(CY5_AF594_TMR(:,1));
% labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
% binsall = zeros(1,size(binsX,2)-1);
% for i = 1:size(binsall,2)
%     binsall(1,i) = size(find(CY5_AF594_TMR(:,1) >= binsX(i)),1)-size(find(CY5_AF594_TMR(:,1) >= binsX(i+1)),1);
% end
% binall = binsall/sum(binsall);
% figure(24); clf
% plot(labelsX,binall,'r','LineWidth',2); hold on
% xlabel('Tsix Molecules')
% ylabel('Probability')
% 
% %%Distribution of Xist for each Jpx
% binnum= 20;
% max_bin = 100;
% %max_bin = max(CY5_AF594_TMR(:,1));
% bin_size = (max_bin-min(CY5_AF594_TMR(:,1)))/binnum;
% binsX = 0:bin_size:max_bin;
% labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
% bins0 = zeros(1,size(binsX,2)-1);
% bins1 = zeros(1,size(binsX,2)-1);
% bins2 = zeros(1,size(binsX,2)-1);
% bins3 = zeros(1,size(binsX,2)-1);
% bins4 = zeros(1,size(binsX,2)-1);
% for i = 1:size(bins0,2)
%     bins0(1,i) = size(find(dist0 >= binsX(i)),1)-size(find(dist0 >= binsX(i+1)),1);
%     bins1(1,i) = size(find(dist1 >= binsX(i)),1)-size(find(dist1 >= binsX(i+1)),1);
%     bins2(1,i) = size(find(dist2 >= binsX(i)),1)-size(find(dist2 >= binsX(i+1)),1);
%     bins3(1,i) = size(find(dist3 >= binsX(i)),1)-size(find(dist3 >= binsX(i+1)),1);
%     bins4(1,i) = size(find(dist4 >= binsX(i)),1)-size(find(dist4 >= binsX(i+1)),1);
% end
% bin0 = bins0/sum(bins0);
% bin1 = bins1/sum(bins1);
% bin2 = bins2/sum(bins2);
% bin3 = bins3/sum(bins3);
% bin4 = bins4/sum(bins4);
% figure(23); clf
% plot(labelsX,bin0,'r','LineWidth',2); hold on
% plot(labelsX,bin1,'k','LineWidth',2);
% plot(labelsX,bin2,'g','LineWidth',2);
% plot(labelsX,bin3,'b','LineWidth',2);
% plot(labelsX,bin4,'m','LineWidth',2);
% xlabel('Tsix Molecules')
% ylabel('Probability')
% legend(Jpx_labels)
% %%% Mean Xist per Jpx, per timepoint
% CY5_AF594_TMR2_tp = zeros(4,size(CY5_AF594_TMR_tp,2));
% for i = 1: size(CY5_AF594_TMR_tp,2)
% for j = 1:4
%     CY5_AF594_TMR2_tp(j,i) = mean(CY5_AF594_TMR_tp{i}(find(CY5_AF594_TMR_tp{i}(:,2) == j),1));
% end
% end
% figure(26); clf; mesh(times1,1:4,CY5_AF594_TMR2_tp)
% xlabel('Time (Days)')
% ylabel('Jpx Molecules')
% zlabel('Mean Xist')
% %hold on; testmesh = CY5_AF594_TMR2_tp + 5; mesh(times1,1:4,testmesh) 
% figure(27); clf; mesh(1:4,times1,CY5_AF594_TMR2_tp')
% ylabel('Time (Days)')
% xlabel('Jpx Molecules')
% zlabel('Mean Xist')
% %%% Correlation Plots
% % load Data_Canada
% % corrplot(DataTable)
