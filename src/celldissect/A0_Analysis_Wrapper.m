%% Define Variables
%%For 20170525 and 20170531
tic
dates0 = {'20180619','20180620','20180621','20180622','20180623','20180626','20180628','20180629'};
times0 = {{'0d'} {'6hr'} {'12hr'},{'1d'},{'2d'},{'5d'},{'9d'},{'9d'}}; 
times1 = [0,.25,.5,1,2,5,9,9];
%%%
max_imnum = 11;
max_tps = 1;
%%% Channel Names. Make empty if channel not used
CY5_name = {'_Xist-CY5-th' '_Xist-CY5-th' '_Xist-CY5-th' ...
    '_Xist-CY5-th' '_Xist-CY5-th' '_Xist-CY5-th' '_Xist-CY5-th' '_Xist-CY5-th'};
AF594_name = {};
%TMR_name = {'_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th'};
%TMR_name = {'_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th'};
%TMR_name = {'_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' };
TMR_name = {'_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' };  
%% copy over thresholds from experiment
thCY5all = {[31 42 50 65 75],[31 42 50 65 75], [31 42 50 65 75],...
    [31 42 50 65 75],[31 42 50 65 75],...
    [31 42 50 65 75],[31 42 50 65 75],[31 42 50 65 75]} ;
thAF594all = {[34 47 55 63 69],[34 47 55 63 69],[34 47 55 63 69],...
    [34 47 55 63 69],[34 47 55 63 69],...
    [34 47 55 63 69],[34 47 55 63 69],[34 47 55 63 69]};
thTMRall = {[98 115 130 145 170] ,[98 115 130 145 170] ,[98 115 130 145 170] ,...
    [98 115 130 145 170] ,[98 115 130 145 170] ,...
    [98 115 130 145 170] ,[98 115 130 145 170],[98 115 130 145 170] } ;
%%% which threshold sets
thNum_CY5 = 3; %Which set of the above thresholds to use (which column)
thNum_TMR = 2; %Which set of the above thresholds to use (which column) Maoa and Tsix should probably be 2
thNum_AF594 = 3; %Which set of the above thresholds to use (which column)

%%% Whether to do Cloud quant or not
CY5_hascloud = 1;
AF594_hascloud = 0;
TMR_hascloud = 1;
%%%
CY5_cloudInt = [10000,10000,10000,...
    10000,10000,10000,10000,10000];
AF594_cloudInt = [];
 TMR_cloudInt = [12550,12550,12550,12550,12550,12550,12550,12550,12550];   %Xist/RepA
 %TMR_cloudInt = [6750,6750,6750,6750,6750,6750,6750,6750,6750];    %Maoa
%TMR_cloudInt = [8650,8650,8650,8650,8650,8650,8650,8650];  %Tsix
%TMR_cloudInt = [8950,8950,8950,8950,8950,8950,8950,8950,8950,8950,8950];    %Pdk3
%%%
nucTH = 'low'
same_dist = 6;   %How close transcripts have to be before they are considered the same transcript
dist_col = 5;
CY5_tran_thres = 3; %How many transcripts needed to be determined as a transcription site or cloud
TMR_tran_thres = 3; %How many transcripts needed to be determined as a transcription site or cloud
AF594_cloudInt
%% Visualize Segmentation Files
%A0A_Segmentation_Checker(dates0,times0,max_tps,max_imnum,CY5_name,AF594_name,TMR_name,nucTH)

%% Gather transcript information
[CY5_AF594_TMR,CY5_AF594_TMR_tp,...
    CY5_in_CY5cloud,CY5_in_CY5cloud_tp,CY5_in_TMRcloud,CY5_in_TMRcloud_tp,TMR_in_CY5cloud,TMR_in_CY5cloud_tp,TMR_in_TMRcloud,TMR_in_TMRcloud_tp,...
    CY5_pos_init,CY5_pos_init_tp,TMR_pos_init,TMR_pos_init_tp,...
    T_C_col,C_T_col,T_C_col_tp,C_T_col_tp,T_T_col,T_T_col_tp...
    ,CY5_int,CY5_int_tp,TMR_int,TMR_int_tp,C_T_col_ind,C_T_col_ind_tp,T_C_col_ind,T_C_col_ind_tp,...
    CY5_Nuc_dist,CY5_Nuc_dist_tp,TMR_Nuc_dist,TMR_Nuc_dist_tp,Cl_Nuc_dist,Cl_Nuc_dist_tp,...
    C1_CY5_dist,C1_CY5_dist_tp,C2_CY5_dist,C2_CY5_dist_tp,C1_TMR_dist,C1_TMR_dist_tp,C2_TMR_dist,C2_TMR_dist_tp]...
    = A1B_RNA_Quant_General(...
    dates0,times0,times1,max_tps,max_imnum,...
    thCY5all,thAF594all,thTMRall,CY5_name,AF594_name,TMR_name,thNum_CY5,thNum_AF594,thNum_TMR,...
    CY5_hascloud,AF594_hascloud,TMR_hascloud,CY5_cloudInt,AF594_cloudInt,TMR_cloudInt,nucTH,...
    dist_col,CY5_tran_thres,TMR_tran_thres); toc

save([CY5_name{1} TMR_name{1} '_' '20181030_Xist_RepA_allcalc_0zadjust'],'CY5_AF594_TMR','CY5_AF594_TMR_tp',...
    'CY5_in_CY5cloud','CY5_in_CY5cloud_tp','CY5_in_TMRcloud','CY5_in_TMRcloud_tp','TMR_in_CY5cloud','TMR_in_CY5cloud_tp','TMR_in_TMRcloud','TMR_in_TMRcloud_tp',...
    'CY5_pos_init','CY5_pos_init_tp','TMR_pos_init','TMR_pos_init_tp',...
    'T_C_col','C_T_col','T_C_col_tp','C_T_col_tp','T_T_col','T_T_col_tp'...
    ,'CY5_int','CY5_int_tp','TMR_int','TMR_int_tp','C_T_col_ind','C_T_col_ind_tp','T_C_col_ind','T_C_col_ind_tp',...
        'CY5_Nuc_dist','CY5_Nuc_dist_tp','TMR_Nuc_dist','TMR_Nuc_dist_tp','Cl_Nuc_dist','Cl_Nuc_dist_tp',...
    'C1_CY5_dist','C1_CY5_dist_tp','C2_CY5_dist','C2_CY5_dist_tp','C1_TMR_dist','C1_TMR_dist_tp','C2_TMR_dist','C2_TMR_dist_tp')
%% Correlate Xist cloud to transcripts and transcription sites

%%%Correlate largest cloud size through all timepoints
figure(101); clf;
CY5_label = 'Xist';
TMR_label = 'Pdk3';
clims=[0,.3];
large_cloud = zeros(size(CY5_AF594_TMR,1),1);        %This will store the largest cloud for each cell in the timecourse
for i = 1:size(CY5_AF594_TMR,1)
    large_cloud(i,1) = max(CY5_AF594_TMR(i,11:12));
end
%CY5_AF594_TMR(find(CY5_AF594_TMR(:,1) == max(CY5_AF594_TMR(:,1))),2:6)
bin_num = 10;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 50;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 0;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
temp1 = large_cloud;
temp2 = temp1;%(temp1 >= minV);
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
binsall = zeros(1,size(binsX,2)-1);
transc_cloud = zeros(3, size(binsall,2));
temp_TMR = CY5_AF594_TMR(:,10);   %Number of transcription sites for TMR
temp_TMR(temp_TMR>2) = 2;     %Make maximum sites = 2
%binsall = imhistc(temp2,bin_num,0,bin_num);
for j = 1:size(CY5_AF594_TMR,1)
    for i = 1:size(binsall,2)
        if i  == size(binsall,2) && temp2(j,1) >= binsX(i)
            transc_cloud(temp_TMR(j,1)+1,i) = transc_cloud(temp_TMR(j,1)+1,i)+1;
        elseif temp2(j,1) >= binsX(i) && temp2(j,1) < binsX(i+1)
            transc_cloud(temp_TMR(j,1)+1,i) = transc_cloud(temp_TMR(j,1)+1,i)+1;
        end
    end
end
    transc_cloud = flipud(transc_cloud);
    transc_cloud = transc_cloud/sum(transc_cloud(:));       %normalize to probability
    imagesc(transc_cloud,clims)
    %title(times0{i}{1})
    if CY5_x
        xlabel([CY5_label ' Transcripts in Cloud'])
        ylabel([TMR_label ' Transcription Sites'])
    else
        xlabel([TMR_label ' Transcription Sites'])
        ylabel([CY5_label ' Transcripts in Cloud'])
    end
    yticks([1 2 3])
    yticklabels({'2','1','0'})
    xticklabels(labelsX)
    colorbar
    h = colorbar;
    title(h, 'Probability') 

%%%Separate by timepoint

figure(102); clf;
CY5_label = 'Xist';
TMR_label = 'Pdk3';
clims=[0,.3];
for k = 1:size(CY5_AF594_TMR_tp,2)
large_cloud = zeros(size(CY5_AF594_TMR_tp{k},1),1);        %This will store the largest cloud for each cell in the timecourse
for i = 1:size(CY5_AF594_TMR_tp{k},1)
    large_cloud(i,1) = max(CY5_AF594_TMR_tp{k}(i,11:12));
end
%CY5_AF594_TMR_tp{k}(find(CY5_AF594_TMR_tp{k}(:,1) == max(CY5_AF594_TMR_tp{k}(:,1))),2:6)
bin_num = 15;
%maxV = max(CY5_AF594_TMR_tp{k}(:,1));
maxV = 150;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 0;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
temp1 = large_cloud;
temp2 = temp1;%(temp1 >= minV);
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
binsall = zeros(1,size(binsX,2)-1);
transc_cloud = zeros(3, size(binsall,2));
temp_TMR = CY5_AF594_TMR_tp{k}(:,10);   %Number of transcription sites for TMR
temp_TMR(temp_TMR>2) = 2;     %Make maximum sites = 2
%binsall = imhistc(temp2,bin_num,0,bin_num);
for j = 1:size(CY5_AF594_TMR_tp{k},1)
    for i = 1:size(binsall,2)
        if i  == size(binsall,2) && temp2(j,1) >= binsX(i)
            transc_cloud(temp_TMR(j,1)+1,i) = transc_cloud(temp_TMR(j,1)+1,i)+1;
        elseif temp2(j,1) >= binsX(i) && temp2(j,1) < binsX(i+1)
            transc_cloud(temp_TMR(j,1)+1,i) = transc_cloud(temp_TMR(j,1)+1,i)+1;
        end
    end
end
    transc_cloud = flipud(transc_cloud);
    transc_cloud = transc_cloud/sum(transc_cloud(:));       %normalize to probability
    subplot(2,4,k)
    imagesc(transc_cloud,clims)
    %title(times0{i}{1})
    if CY5_x
        xlabel([CY5_label ' Transcripts in Cloud'])
        ylabel([TMR_label ' Transcription Sites'])
    else
        xlabel([TMR_label ' Transcription Sites'])
        ylabel([CY5_label ' Transcripts in Cloud'])
    end
    yticks([1 2 3])
    yticklabels({'2','1','0'})
    xticks(1:2:size(binsall,2))
    xticklabels(labelsX(1:2:size(binsall,2)))
    colorbar
    h = colorbar;
    title(h, 'Probability') 
end


%% Scatterplots for every transcript (pooled timepoints)

figure(20); clf; gscatter(CY5_AF594_TMR(:,2),CY5_AF594_TMR(:,1),CY5_AF594_TMR(:,4),[],'o',10);
xlabel('Tsix','Fontsize',16)
ylabel('Xist','Fontsize',16)
legend('0 day','2 day','3 day','4 day','5 day','6 day','11 day','13 day','FontSize',20)

figure(40); clf; gscatter(CY5_AF594_TMR(:,2),CY5_AF594_TMR(:,3),CY5_AF594_TMR(:,4),[],'o',10);
xlabel('Tsix','Fontsize',16)
ylabel('Jpx','Fontsize',16)
legend('0 day','2 day','3 day','4 day','5 day','6 day','11 day','13 day','FontSize',20)

figure(50); clf; gscatter(CY5_AF594_TMR(:,3),CY5_AF594_TMR(:,1),CY5_AF594_TMR(:,4),[],'o',10);
xlabel('Tsix','Fontsize',16)
ylabel('Xist','Fontsize',16)
legend({'day 0','6 hours','12 hours','1 day','2 days','5 days','9 days'},'FontSize',14)

figure(51); clf; scatplot(CY5_AF594_TMR(:,3),CY5_AF594_TMR(:,1),'circles',15,100,5,3,15);
ylim([0 500])
xlabel('Pdk3','Fontsize',16)
ylabel('Xist','Fontsize',16)

figure(51); clf; scatplot(CY5_AF594_TMR(:,1),CY5_AF594_TMR(:,3),'circles',15,100,5,3,15);
xlim([0 500])
xlabel('Xist','Fontsize',16)
ylabel('Maoa','Fontsize',16)
%legend({'day 0','6 hours','12 hours','1 day','2 days','5 days','9 days'},'FontSize',14)
%% Transcription site distributions
figure(90); clf
%CY5_AF594_TMR(find(CY5_AF594_TMR(:,1) == max(CY5_AF594_TMR(:,1))),2:6)
bin_num = 10;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 10;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 0;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
temp1 = CY5_AF594_TMR(:,10);
temp2 = temp1(temp1 >= minV);
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
binsall = zeros(1,size(binsX,2)-1);
%binsall = imhistc(temp2,bin_num,0,bin_num);
for i = 1:size(binsall,2)
    binsall(1,i) = size(find(temp2 >= binsX(i)),1)-size(find(temp2 >= binsX(i+1)),1);
end
binall = binsall/sum(binsall);  %normalize
bar(labelsX,binsall)
%plot(labelsX,binall,'r','LineWidth',3); hold on
title('Distribution for TMR Transcription Sites')
xlabel('TMR Transcription Sites')
ylabel('Probability')




%% Transcription Site correlations
CY5_x = 1;      %Make 1 if you want the CY5 channel on the x axis
CY5_label = 'Xist';
TMR_label = 'RepA';
clims=[0,.5]
figure(91); clf;

transc_corr = zeros(3,3,size(CY5_AF594_TMR_tp,2));     %will store how many cells are in each bin 
for i = 1:size(CY5_AF594_TMR_tp,2)
    temp_CY = CY5_AF594_TMR_tp{i}(:,8);     %Number of transcription sites for CY5
    temp_TMR = CY5_AF594_TMR_tp{i}(:,10);   %Number of transcription sites for TMR
    temp_CY(temp_CY>2) = 2;     %Make maximum sites = 2
    temp_TMR(temp_TMR>2) = 2;     %Make maximum sites = 2
    for j = 1:size(CY5_AF594_TMR_tp{i},1)
        if CY5_x
            transc_corr(temp_TMR(j)+1,temp_CY(j)+1,i) = transc_corr(temp_TMR(j)+1,temp_CY(j)+1,i)+1;
        else
            transc_corr(temp_CY(j)+1,temp_TMR(j)+1,i) = transc_corr(temp_CY(j)+1,temp_TMR(j)+1,i)+1;
        end
    end
    transc_corr(:,:,i) = flipud(transc_corr(:,:,i));
    temp_tc = transc_corr(:,:,i);       %correlation for this timepoint
    temp_tc = temp_tc/sum(temp_tc(:));       %normalize to probability
    subplot(2,4,i)
    imagesc(temp_tc,clims)
    title(times0{i}{1})
    if CY5_x
        xlabel([CY5_label ' Transcription Sites'])
        ylabel([TMR_label ' Transcription Sites'])
    else
        xlabel([TMR_label ' Transcription Sites'])
        ylabel([CY5_label ' Transcription Sites'])
    end
    xticks([1 2 3])
    xticklabels({'0','1','2'})
    yticks([1 2 3])
    yticklabels({'2','1','0'})
    colorbar
    h = colorbar;
    title(h, 'Probability') 
end
transc_corr 

figure(92); clf
transc_corr_all = zeros(3);     %will store how many cells are in each bin 
    temp_CY = CY5_AF594_TMR(:,8);     %Number of transcription sites for CY5
    temp_TMR = CY5_AF594_TMR(:,10);   %Number of transcription sites for TMR
    temp_CY(temp_CY>2) = 2;     %Make maximum sites = 2
    temp_TMR(temp_TMR>2) = 2;     %Make maximum sites = 2
    for j = 1:size(CY5_AF594_TMR,1)
        if CY5_x
            transc_corr_all(temp_TMR(j)+1,temp_CY(j)+1) = transc_corr_all(temp_TMR(j)+1,temp_CY(j)+1)+1;
        else
            transc_corr_all(temp_CY(j)+1,temp_TMR(j)+1) = transc_corr_all(temp_CY(j)+1,temp_TMR(j)+1)+1;
        end
    end
transc_corr_all 
%transc_corr_all = fliplr(transc_corr_all);
transc_corr_all = flipud(transc_corr_all);
%HeatMap(transc_corr_all,'Symmetric','False','Colormap','redbluecmap','RowLabels',{'0','1','2'},'ColumnLabels',{'0','1','2'})
imagesc(transc_corr_all)
    if CY5_x
        xlabel([CY5_label ' Transcription Sites'])
        ylabel([TMR_label ' Transcription Sites'])
    else
        xlabel([TMR_label ' Transcription Sites'])
        ylabel([CY5_label ' Transcription Sites'])
    end
xticks([1 2 3])
xticklabels({'0','1','2'})
yticks([1 2 3])
yticklabels({'2','1','0'})
colorbar
h = colorbar;
title(h, 'Number of Cells')

%% Distributions of each transcript (for positive cells)
%%%Xist
figure(21); clf;
%CY5_AF594_TMR(find(CY5_AF594_TMR(:,1) == max(CY5_AF594_TMR(:,1))),2:6)
bin_num = 15;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 300;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 1;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
temp1 = CY5_AF594_TMR(:,1);
temp2 = temp1(temp1 >= minV);
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
binsall = zeros(1,size(binsX,2)-1);
%binsall = imhistc(temp2,bin_num,0,bin_num);
for i = 1:size(binsall,2)
    binsall(1,i) = size(find(temp2 >= binsX(i)),1)-size(find(temp2 >= binsX(i+1)),1);
end
binall = binsall/sum(binsall);  %normalize
bar(labelsX,binsall)
%plot(labelsX,binall,'r','LineWidth',3); hold on
title('Distribution for Xist Positive Cells')
xlabel('Xist Molecules')
ylabel('Number of Cells')

%%%Tsix
figure(22); clf;
%CY5_AF594_TMR(find(CY5_AF594_TMR(:,1) == max(CY5_AF594_TMR(:,1))),2:6)
bin_num = 10;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 30;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 1;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
temp1 = CY5_AF594_TMR(:,2);
temp2 = temp1(temp1 >= minV);
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
binsall = zeros(1,size(binsX,2)-1);
%binsall = imhistc(temp2,bin_num,0,maxV);
for i = 1:size(binsall,2)
     binsall(1,i) = size(find(temp2 >= binsX(i)),1)-size(find(temp2 >= binsX(i+1)),1);
end
binall = binsall/sum(binsall);  %normalize
bar(labelsX,binsall)
%plot(labelsX,binall,'b','LineWidth',3); hold on
title('Distribution for Tsix Positive Cells')
xlabel('Tsix Molecules')
ylabel('Number of Cells')

%%%Jpx
% for i = 1:5
%     thNum_TMR = i;
%     [CY5_AF594_TMR,CY5_AF594_TMR_tp] = A1B_RNA_Quant_General(dates0,times0,times1,max_tps,max_imnum,...
%     thCY5all,thAF594all,thTMRall,CY5_name,AF594_name,TMR_name,thNum_CY5,thNum_AF594,thNum_TMR,...
%     CY5_hascloud,AF594_hascloud,TMR_hascloud,CY5_cloudInt,AF594_cloudInt,TMR_cloudInt,nucTH)
% figure(23+i); clf;
%CY5_AF594_TMR(find(CY5_AF594_TMR(:,3) == max(CY5_AF594_TMR(:,3))),2:6)
% CY5_AF594_TMR(find(CY5_AF594_TMR(:,3) > 60),:)
% figure(100); clf; imshow(cells,[])
% figure(101); clf; imshow(cells == 25,[])
figure(23); clf;
bin_num = 20;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 100;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 1;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
temp1 = CY5_AF594_TMR(:,3);
temp2 = temp1(temp1 >= minV);
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1);%+bin_size/2;
binsall = zeros(1,size(binsX,2)-1);
%binsall = imhistc(temp2,bin_num,0,maxV);
for i = 1:size(binsall,2)
     binsall(1,i) = size(find(temp2 >= binsX(i)),1)-size(find(temp2 >= binsX(i+1)),1);
end
%binall = binsall/sum(binsall);  %normalize
bar(labelsX,binsall)
%plot(labelsX,binall,'g','LineWidth',3); hold on
title('Distribution for Pdk3 Positive Cells')
xlabel('Pdk3 Molecules')
ylabel('Number of Cells')
%end
%% Distributions over time

figure(24); clf;
colors = {'k' 'b' 'm' 'r' 'g' 'c' [255 159 0]/(255+159)}
bin_num = 15;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 149;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 0;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
temp1 = CY5_AF594_TMR(:,3);
temp2 = temp1(temp1 >= minV);
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1);%+bin_size/2;
binsall = {};
for j = 1:7
temp1 = CY5_AF594_TMR_tp{j}(:,3);
temp2 = temp1(temp1 >= minV);
binsall{j} = zeros(1,size(binsX,2)-1);
%binsall = imhistc(temp2,bin_num,0,maxV);
for i = 1:size(binsall{j},2)
     binsall{j}(1,i) = size(find(temp2 >= binsX(i)),1)-size(find(temp2 >= binsX(i+1)),1);
end
binall = binsall{j}/sum(binsall{j});  %normalize
plot(labelsX,binall,'Color',colors{j},'LineWidth',3); hold on
end
%binall = binsall/sum(binsall);  %normalize
%bar(labelsX,binsall)
%plot(labelsX,binall,'g','LineWidth',3); hold on
legend('0d', '6hr', '12hr', '1d', '2d', '5d','9d')
%title('Distribution for Maoa Positive Cells')
xlabel('Pdk3 Molecules')
%ylabel('Number of Cells')
ylabel('Probability')
%end

figure(25); clf;
colors = {'k' 'b' 'm' 'r' 'g' 'c' [255 159 0]/(255+159)}
bin_num = 15;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 299;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 1;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
temp1 = CY5_AF594_TMR(:,1);
temp2 = temp1(temp1 >= minV);
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1);%+bin_size/2;
binsall = {};
for j = 1:7
temp1 = CY5_AF594_TMR_tp{j}(:,1);
temp2 = temp1(temp1 >= minV);
binsall{j} = zeros(1,size(binsX,2)-1);
%binsall = imhistc(temp2,bin_num,0,maxV);
for i = 1:size(binsall{j},2)
     binsall{j}(1,i) = size(find(temp2 >= binsX(i)),1)-size(find(temp2 >= binsX(i+1)),1);
end
binall = binsall{j}/sum(binsall{j});  %normalize
plot(labelsX,binall,'Color',colors{j},'LineWidth',3); hold on
end
%binall = binsall/sum(binsall);  %normalize
%bar(labelsX,binsall)
%plot(labelsX,binall,'g','LineWidth',3); hold on
legend('0d', '6hr', '12hr', '1d', '2d', '5d','9d')
%title('Distribution for Maoa Positive Cells')
xlabel('Xist Molecules')
%ylabel('Number of Cells')
ylabel('Probability')
%end



%% Boxplots:
%%%Boxplots with Xist on x axis (if box_col = 1) or tmr on x axis (box_col
%%%=3)
clear dists1 X_labels
x_label = 'Xist';
y_label = 'Maoa';
y_col = 3;          %column used for values on y axis
% CY5_AF594_TMR(:,18) = log(CY5_AF594_TMR(:,1));
% for i = 1:size(CY5_AF594_TMR_tp,2)
%     CY5_AF594_TMR_tp{i}(:,18) = log(CY5_AF594_TMR_tp{i}(:,y_col));
% end
% y_col = 18;          %column used for values on y axis
box_col = 1;        %column used for values on x axis
dists1 = {}; %store into individual matrices the cells for each bin (not used)
bin_num = 10;
maxV = 200;         %cutoff for max value shown (x axis). Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 0;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1); %+bin_size/2;
CY5_AF594_TMR2 = CY5_AF594_TMR;
X_labels = {};
for i = 1:bin_num+1
    if i == bin_num+1
        X_labels{i} = [num2str(binsX(i)) '+'];                  %If last label, make it show that it can be anything above as well
    else
        if binsX(i+1)-1-binsX(i) > 0
            X_labels{i} = [num2str(binsX(i)) '-' num2str(binsX(i+1)-1)];      %Make labels for the bins
        else
            X_labels{i} = num2str(binsX(i));
        end
    end
end
%maxV = max(CY5_AF594_TMR(:,1));
for i = 1:bin_num+1
    dists1{i} = zeros(2,size(CY5_AF594_TMR,2)-1);
    counter1 = 1;
    for j = 1:size(CY5_AF594_TMR,1)
        if i == bin_num+1
            if CY5_AF594_TMR(j,box_col) >= binsX(i)
           CY5_AF594_TMR2(j,8) = i-1;
           counter1 = counter1+1;            
            end
       elseif CY5_AF594_TMR(j,box_col) >= binsX(i) & CY5_AF594_TMR(j,box_col) < binsX(i+1)
           CY5_AF594_TMR2(j,8) = i-1;
           counter1 = counter1+1;
        end
    end
end
labels_here = unique(CY5_AF594_TMR2(:,8));
true_end = size(CY5_AF594_TMR2,1);    %The size of the matrix before fake values are added
for m = 0:size(X_labels,2)-1
    if sum(labels_here == m) == 0
        CY5_AF594_TMR2(size(CY5_AF594_TMR2,1)+1,8) = m;
    end
end
%%%Boxplot of TMR on y axis   
figure(24); clf; 
boxplot(CY5_AF594_TMR2(:,y_col),CY5_AF594_TMR2(:,8),'Labels',X_labels)
%ylim([0 100])
xlabel([x_label ' per Cell'],'Fontsize',10)
ylabel([y_label ' per Cell'],'Fontsize',10)
hold on
plot(0:bin_num+2,mean(CY5_AF594_TMR2(1:true_end,y_col)):.0001:mean(CY5_AF594_TMR2(1:true_end,y_col))+.0001*(bin_num+2),'--r','LineWidth',2)


%%% Same as above, but changing timepoints
figure(25); clf;
clear dists1 X_labels
% y_col = 3;          %column used for values on y axis
% box_col = 1;        %column used for values on x axis
dists1 = {}; %store into individual matrices the cells for each bin (not used)
% bin_num = 10;
% maxV =60;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
% minV = 0;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1); %+bin_size/2;
CY5_AF594_TMR_tp2 = CY5_AF594_TMR_tp;
X_labels = {};
for i = 1:bin_num+1
    if i == bin_num+1
        X_labels{i} = [num2str(binsX(i)) '+'];                  %If last label, make it show that it can be anything above as well
    else
        if binsX(i+1)-1-binsX(i) > 0
            X_labels{i} = [num2str(binsX(i)) '-' num2str(binsX(i+1)-1)];      %Make labels for the bins
        else
            X_labels{i} = num2str(binsX(i));
        end
    end
end
%maxV = max(CY5_AF594_TMR(:,1));
for k = 1:size(CY5_AF594_TMR_tp2,2)
for i = 1:bin_num+1
    counter1 = 1;
    for j = 1:size(CY5_AF594_TMR_tp{k},1)
        if i == bin_num+1
            if CY5_AF594_TMR_tp2{k}(j,box_col) >= binsX(i)
           CY5_AF594_TMR_tp2{k}(j,8) = i-1;
           counter1 = counter1+1;            
            end
       elseif CY5_AF594_TMR_tp2{k}(j,box_col) >= binsX(i) & CY5_AF594_TMR_tp2{k}(j,box_col) < binsX(i+1)
           CY5_AF594_TMR_tp2{k}(j,8) = i-1;
           counter1 = counter1+1;
        end
    end
end
labels_here = unique(CY5_AF594_TMR_tp2{k}(:,8));
true_end = size(CY5_AF594_TMR_tp2{k},1);    %The size of the matrix before fake values are added
for m = 0:size(X_labels,2)-1
    if sum(labels_here == m) == 0
        CY5_AF594_TMR_tp2{k}(size(CY5_AF594_TMR_tp2{k},1)+1,8) = m;
    end
end
        
%%%Boxplot 
subplot(2,4,k) 
h = boxplot(CY5_AF594_TMR_tp2{k}(:,y_col),CY5_AF594_TMR_tp2{k}(:,8),'Labels',X_labels)
xlabel([x_label ' per Cell'],'Fontsize',10)
ylabel([y_label ' per Cell'],'Fontsize',10)
title(times0{k})
ylim([0 max(CY5_AF594_TMR(:,y_col))])
%ylim([0 100])
xtickangle(45)
hold on
plot(0:bin_num+2,mean(CY5_AF594_TMR_tp2{k}(1:true_end,y_col)):.0001:mean(CY5_AF594_TMR_tp2{k}(1:true_end,y_col))+.0001*(bin_num+2),'--r','LineWidth',2)
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
set(h,'LineWidth',2)
end
 

%%%Old boxplot code
% boxplot of Tsix for different numbers of Jpx
%%%This part first sets bins for the TMR channel and then makes a column
%%%where the bin number is specified
clear dists1
dists1 = {}; %store into individual matrices the cells for each bin (not used)
bin_num = 15;
maxV = 149;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 0;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1); %+bin_size/2;
CY5_AF594_TMR2 = CY5_AF594_TMR;
%maxV = max(CY5_AF594_TMR(:,1));
for i = 1:bin_num+1
    dists1{i} = zeros(2,size(CY5_AF594_TMR,2)-1);
    counter1 = 1;
    for j = 1:size(CY5_AF594_TMR,1)
        if i == bin_num+1
            if CY5_AF594_TMR(j,3) >= binsX(i)
           dists1{i}(counter1,:)= CY5_AF594_TMR(j,[1,2,4:size(CY5_AF594_TMR,2)])
           CY5_AF594_TMR2(j,7) = i-1;
           counter1 = counter1+1;            
            end
       elseif CY5_AF594_TMR(j,3) >= binsX(i) & CY5_AF594_TMR(j,3) < binsX(i+1)
           dists1{i}(counter1,:)= CY5_AF594_TMR(j,[1,2,4:size(CY5_AF594_TMR,2)])
           CY5_AF594_TMR2(j,7) = i-1;
           counter1 = counter1+1;
        end
    end
end


%%%Boxplot of Xist
figure(25); clf; 
boxplot(CY5_AF594_TMR2(:,1),CY5_AF594_TMR2(:,7));%,'Labels',Jpx_labels)
xlabel('Jpx per Cell','Fontsize',10)
ylabel('Xist per Cell','Fontsize',10)
hold on
plot(0:7,mean(CY5_AF594_TMR2(:,1)):.0001:mean(CY5_AF594_TMR2(:,1))+.0007,'--r')

times2 = [0,1,2,3,4,5,6,11,13];
%%%Boxplots of Tsix 
    figure(25); clf; hold on
    figure(26); clf; hold on
for i = 1:size(times2,2)
    tempT =  CY5_AF594_TMR2(find(CY5_AF594_TMR2(:,4) == times2(i)),:)
    %%%Plot Tsix
    figure(25)
    if i <= 5
        subplot(2,5,i)
    else
        subplot(2,4,i-1)
    end
    boxplot(tempT(:,2),tempT(:,7));%,'Labels',Jpx_labels)
    xlabel('Jpx per Cell','Fontsize',10)
    ylabel('Tsix per Cell','Fontsize',10)
    title([num2str(times2(i)) ' Days'])
    hold on
    plot(0:7,mean(tempT(:,2)):.0001:mean(tempT(:,2))+.0007,'--r')
    %%%Plot Xist
    figure(26)
    if i <= 5
        subplot(2,5,i)
    else
        subplot(2,4,i-1)
    end
    boxplot(tempT(:,1),tempT(:,7));%,'Labels',Jpx_labels)
    xlabel('Jpx per Cell','Fontsize',10)
    ylabel('Xist per Cell','Fontsize',10)
    title([num2str(times2(i)) ' Days'])
    hold on
    plot(0:7,mean(tempT(:,2)):.0001:mean(tempT(:,2))+.0007,'--r')
end
times3 = [0,2,3,4]

%% Average Transcripts per timepoint

times = [0,.25,.5,1,2,5,9];
avgs = zeros(5,1);
for i = 1:7
    avgs(i) = mean(CY5_AF594_TMR_tp{i}(:,3));
end
avgs = avgs/max(avgs(:));
figure(44); clf; 
plot(times,avgs,'r-o','LineWidth',3)
%title('Average Pdk3 Over Time')
xlabel('Time (Days)')
ylabel('Number of Tsix Molecules')
ylim([0 max(avgs)+5])
hold on
figure(44); clf; 
times = [0,.25,.5,1,2,5,9];
avgs = zeros(5,1);
for i = 1:7
    avgs(i) = mean(CY5_AF594_TMR_tp{i}(:,1));
end
%figure(44); clf; 
avgs = avgs/max(avgs(:));
plot(times,avgs,'b-o','LineWidth',3); hold on
%title('Normalized Average Expression Over Time')
xlabel('Time (Days)')
ylabel('Number of Molecules')
ylim([0 max(avgs)+5])
ylim([0 1])
legend('Tsix','Xist')
%%%Plot on fraction
avgs = zeros(5,1);
thresh1=10;
for i = 1:7
    avgs(i) = sum(CY5_AF594_TMR_tp{i}(:,1)>thresh1)/size(CY5_AF594_TMR_tp{i},1);
end
%figure(44); clf; 
%avgs = avgs/max(avgs(:));
plot(times,avgs,'r-o','LineWidth',3); hold on
%%%Plot cells with one cloud
avgs = zeros(5,1);
thresh1=10;
for i = 1:7
    avgs(i) = sum(CY5_AF594_TMR_tp{i}(:,11)>thresh1)/size(CY5_AF594_TMR_tp{i},1);
end
%figure(44); clf; 
%avgs = avgs/max(avgs(:));
plot(times,avgs,'k-o','LineWidth',3); hold on
%%%Plot cells with two clouds
avgs = zeros(5,1);
thresh1=10;
for i = 1:7
    both = zeros(1,size(CY5_AF594_TMR_tp{i},1));
    for m = 1:size(CY5_AF594_TMR_tp{i},1)
        both(m) = CY5_AF594_TMR_tp{i}(m,12)>thresh1 & CY5_AF594_TMR_tp{i}(m,11)>thresh1;
    end
    avgs(i) = sum(both)/size(CY5_AF594_TMR_tp{i},1);
end
%figure(44); clf; 
%avgs = avgs/max(avgs(:));
plot(times,avgs,'g-o','LineWidth',3); hold on
%title('Normalized Average Expression Over Time')
xlabel('Time (Days)')
ylabel('Fraction')
ylim([0 max(avgs)+5])
ylim([0 1.3])
set(gca,'FontSize',14)
legend('Normalized Average Xist',['On Fraction Xist (above ' num2str(thresh1) ')'],['Fraction of Cells with at least One Cloud (' num2str(thresh1) ' transcripts)'],['Fraction of Cells with at least Two Clouds (' num2str(thresh1) ' transcripts)']  )

%% Joint probability
T1_T3 = zeros(1,2);
T1_T3(1:size(CY5_AF594_TMR(:,1),1),1)= CY5_AF594_TMR(:,1);
T1_T3(1:size(CY5_AF594_TMR(:,1),1),2)= CY5_AF594_TMR(:,3);
figure(80); clf
hist3(T1_T3)

%%%by timepoint
T1_T3 = zeros(1,2);
T1_T3(1:size(CY5_AF594_TMR_tp{7}(:,1),1),1)= CY5_AF594_TMR_tp{7}(:,1);
T1_T3(1:size(CY5_AF594_TMR_tp{7}(:,1),1),2)= CY5_AF594_TMR_tp{7}(:,3);
figure(80); clf
hist3(T1_T3)


%% Repeat of on fraction cells
figure(81); clf
colors = {'k' 'b' 'g' 'm' 'r' 'y' 'c'}
times2 = [0,.25,.5,1,2,5,9];
Xist_thres = 3;  %equal or greater is positive
Tsix_thres = 3; %equal or greater is positive
Tsix_only = zeros(size(times2));
Xist_only = zeros(size(times2));
Tsix_Xist = zeros(size(times2));
Neither = zeros(size(times2));
for i = 1:size(times2,2)
    for j = 1:size(CY5_AF594_TMR_tp{i},1)
        if CY5_AF594_TMR_tp{i}(j,1) >=  Xist_thres & CY5_AF594_TMR_tp{i}(j,3) >= Tsix_thres
            Tsix_Xist(i) = Tsix_Xist(i)+1;
        elseif CY5_AF594_TMR_tp{i}(j,1) <  Xist_thres & CY5_AF594_TMR_tp{i}(j,3) < Tsix_thres
            Neither(i) = Neither(i)+1;
        elseif CY5_AF594_TMR_tp{i}(j,1) >=  Xist_thres
            Xist_only(i) = Xist_only(i) + 1;
        elseif CY5_AF594_TMR_tp{i}(j,3) >= Tsix_thres
            Tsix_only(i) = Tsix_only(i) + 1;
        end
    end
    Tsix_Xist(i)+Neither(i)+Xist_only(i)+Tsix_only(i)
    total = size(CY5_AF594_TMR_tp{i},1)
    Tsix_Xist(i) = Tsix_Xist(i)/total;
    Neither(i) = Neither(i)/total;
    Xist_only(i) = Xist_only(i)/total;
    Tsix_only(i) = Tsix_only(i)/total;
end
plot(times2,Xist_only,'r-o','LineWidth',3)
hold on
plot(times2,Tsix_only,'b-o','LineWidth',3)
plot(times2,Tsix_Xist,'m-o','LineWidth',3)
plot(times2,Neither,'k-o','LineWidth',3)
legend('Xist Only','Maoa Only', 'Both','Neither')
xlabel('Time (Days)')
ylabel('Probability')
ylim([0,1])

%% On fraction (including both)
figure(82); clf
colors = {'k' 'b' 'g' 'm' 'r' 'y' 'c'}
times2 = [0,.25,.5,1,2,5,9];
Xist_thres = 3;  %equal or greater is positive
Tsix_thres = 3; %equal or greater is positive
Tsix_on = zeros(size(times2));
Xist_on = zeros(size(times2));
Tsix_Xist = zeros(size(times2));
Neither = zeros(size(times2));
for i = 1:size(times2,2)
    Xist_pos_temp = CY5_AF594_TMR_tp{i}(:,1)>Xist_thres;    %Matrix indicating positive CY5
    Tsix_pos_temp = CY5_AF594_TMR_tp{i}(:,3)>Tsix_thres;    %Matrix indicating positive TMR
    Tsix_on(i) = sum(Tsix_pos_temp);    
    Xist_on(i) = sum(Xist_pos_temp);
    Tsix_Xist(i) = sum(immultiply(Tsix_pos_temp,Xist_pos_temp));
    Xist_or_Tsix = Tsix_pos_temp+Xist_pos_temp;
    Xist_or_Tsix = Xist_or_Tsix>0;
    total = size(CY5_AF594_TMR_tp{i},1);
    Neither(i) = total-sum(Xist_or_Tsix(:));
    Tsix_Xist(i) = Tsix_Xist(i)/total;
    Neither(i) = Neither(i)/total;
    Xist_on(i) = Xist_on(i)/total;
    Tsix_on(i) = Tsix_on(i)/total;
    Neither(i)+ Xist_on(i) +Tsix_on(i)-Tsix_Xist(i);
end
plot(times2,Xist_on,'r-o','LineWidth',3)
hold on
plot(times2,Tsix_on,'b-o','LineWidth',3)
plot(times2,Tsix_Xist,'m-o','LineWidth',3)
plot(times2,Neither,'k-o','LineWidth',3)
legend('Xist On','Maoa On', 'Both','Neither')
xlabel('Time (Days)')
ylabel('Probability')
ylim([0,1])

%%% Plot all transcript dynamics
figure(83); clf       %For on fractions
figure(84); clf        %for averages
figure(85); clf
colors = {'k' 'b' 'g' 'm' 'r' 'y' 'c'}
times2 = [0,.25,.5,1,2,5,9];
TMR_names = {{'_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th'},...
{'_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th'},...
{'_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' },...
{'_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' }};  
for j = 1:size(TMR_names,2)
Xist_thres = 3;  %equal or greater is positive
Tsix_thres = 3; %equal or greater is positive
if j == 4
    TMR_hascloud = 1;
else
    TMR_hascloud = 0;
end

TMR_name = TMR_names{j};
[CY5_AF594_TMR,CY5_AF594_TMR_tp,T_C_col,C_T_col,T_C_col_tp,C_T_col_tp,T_T_col] = A1B_RNA_Quant_General(...
    dates0,times0,times1,max_tps,max_imnum,...
    thCY5all,thAF594all,thTMRall,CY5_name,AF594_name,TMR_name,thNum_CY5,thNum_AF594,thNum_TMR,...
    CY5_hascloud,AF594_hascloud,TMR_hascloud,CY5_cloudInt,AF594_cloudInt,TMR_cloudInt,nucTH,...
    dist_col,CY5_tran_thres,TMR_tran_thres);
Tsix_on = zeros(size(times2));
Xist_on = zeros(size(times2));
Tsix_Xist = zeros(size(times2));
Neither = zeros(size(times2));
Xist_avg =  zeros(size(times2));
Tsix_avg =  zeros(size(times2));
for i = 1:size(times2,2)
    Xist_pos_temp = CY5_AF594_TMR_tp{i}(:,1)>Xist_thres;    %Matrix indicating positive CY5
    Tsix_pos_temp = CY5_AF594_TMR_tp{i}(:,3)>Tsix_thres;    %Matrix indicating positive TMR
    Tsix_on(i) = sum(Tsix_pos_temp);    
    Xist_on(i) = sum(Xist_pos_temp);
    Tsix_Xist(i) = sum(immultiply(Tsix_pos_temp,Xist_pos_temp));
    Xist_or_Tsix = Tsix_pos_temp+Xist_pos_temp;
    Xist_or_Tsix = Xist_or_Tsix>0;
    Xist_avg(i) = nanmean(CY5_AF594_TMR_tp{i}(:,1));
    Tsix_avg(i) = nanmean(CY5_AF594_TMR_tp{i}(:,3));
    total = size(CY5_AF594_TMR_tp{i},1);
    Neither(i) = total-sum(Xist_or_Tsix(:));
    Tsix_Xist(i) = Tsix_Xist(i)/total;
    Neither(i) = Neither(i)/total;
    Xist_on(i) = Xist_on(i)/total;
    Tsix_on(i) = Tsix_on(i)/total;
    Neither(i)+ Xist_on(i) +Tsix_on(i)-Tsix_Xist(i);
end
Xist_avg_norm = Xist_avg/max(Xist_avg(:));
Tsix_avg_norm = Tsix_avg/max(Tsix_avg(:));
figure(83)
hold on
if j == 1
    plot(times2,Xist_on,'k-o','LineWidth',3)
end
plot(times2,Tsix_on,[colors{j+1} '-o'],'LineWidth',3)
%plot(times2,Tsix_Xist,'m-o','LineWidth',3)
%plot(times2,Neither,'k-o','LineWidth',3)
figure(84)
hold on
if j == 1
    plot(times2,Xist_avg,'k-o','LineWidth',3)
end
plot(times2,Tsix_avg,[colors{j+1} '-o'],'LineWidth',3)
%plot(times2,Tsix_Xist,'m-o','LineWidth',3)
%plot(times2,Neither,'k-o','LineWidth',3)
figure(85)
hold on
if j == 1
    plot(times2,Xist_avg_norm,'k-o','LineWidth',3)
end
plot(times2,Tsix_avg_norm,[colors{j+1} '-o'],'LineWidth',3)
%plot(times2,Tsix_Xist,'m-o','LineWidth',3)
%plot(times2,Neither,'k-o','LineWidth',3)
end
figure(83)
legend('Xist On','Pdk3 On','Maoa On','Tsix On','RepA On', 'Both','Neither')
xlabel('Time (Days)')
ylabel('Probability')
ylim([0,1])
figure(84)
legend('Xist','Pdk3','Maoa','Tsix','RepA', 'Both','Neither')
xlabel('Time (Days)')
ylabel('Average Expression')
figure(85)
legend('Xist','Pdk3','Maoa','Tsix','RepA', 'Both','Neither')
xlabel('Time (Days)')
ylabel('Normalized Average Expression')
ylim([0,1])

%% colocalization
figure(94); clf
xs = 1:20;
ys = zeros(2,size(xs,2));
for i = xs
    dist_col = i;
  [CY5_AF594_TMR,CY5_AF594_TMR_tp,T_C_col,C_T_col,T_C_col_tp,C_T_col_tp,T_T_col] = A1B_RNA_Quant_General(...
    dates0,times0,times1,max_tps,max_imnum,...
    thCY5all,thAF594all,thTMRall,CY5_name,AF594_name,TMR_name,thNum_CY5,thNum_AF594,thNum_TMR,...
    CY5_hascloud,AF594_hascloud,TMR_hascloud,CY5_cloudInt,AF594_cloudInt,TMR_cloudInt,nucTH,...
    dist_col,CY5_tran_thres,TMR_tran_thres);  
% sum(T_C_col(:)>0)
% sum(T_C_col(:)>1)
% sum(C_T_col(:)>0)
% sum(C_T_col(:)>1)

ys(1,i) = sum(T_C_col(:)>1)/sum(T_C_col(:)>0);      %how much TMR has a colocalized CY5 spot
ys(2,i) = sum(C_T_col(:)>1)/sum(C_T_col(:)>0);      %how much CY5 correlates with TMR
end
plot(xs,ys(1,:),'LineWidth',3,'Color','r')
hold on
plot(xs,ys(2,:),'LineWidth',3,'Color','b')
xlabel('Pixels Apart')
ylabel('Fraction of Spots')
ylim([0 1])
legend({'Fraction RepA with Colocalized Xist','Fraction Xist with Colocalized RepA'})


sum(T_T_col(:)==1)
sum(T_T_col(:)==2)
for i = 1:7
    i
%     sum(T_C_col_tp{i}(:)>1)
%     sum(T_C_col_tp{i}(:)>0)
sum(T_C_col_tp{i}(:)>1)/sum(T_C_col_tp{i}(:)>0)
sum(C_T_col_tp{i}(:)>1)/sum(C_T_col_tp{i}(:)>0)
end

%% colocalization by timepoint

T_C_col_tp
C_T_col_tp
figure(94); clf
xs = 1:20;
ys = zeros(2,size(xs,2),7);
for i = xs
    dist_col = i;
  [CY5_AF594_TMR,CY5_AF594_TMR_tp,T_C_col,C_T_col,T_C_col_tp,C_T_col_tp,T_T_col] = A1B_RNA_Quant_General(...
    dates0,times0,times1,max_tps,max_imnum,...
    thCY5all,thAF594all,thTMRall,CY5_name,AF594_name,TMR_name,thNum_CY5,thNum_AF594,thNum_TMR,...
    CY5_hascloud,AF594_hascloud,TMR_hascloud,CY5_cloudInt,AF594_cloudInt,TMR_cloudInt,nucTH,...
    dist_col,CY5_tran_thres,TMR_tran_thres);  
% sum(T_C_col(:)>0)
% sum(T_C_col(:)>1)
% sum(C_T_col(:)>0)
% sum(C_T_col(:)>1)
for j = 1:size(T_C_col_tp,2)
ys(1,i,j) = sum(T_C_col_tp{j}(:)>1)/sum(T_C_col_tp{j}(:)>0);      %how much TMR has a colocalized CY5 spot
ys(2,i,j) = sum(C_T_col_tp{j}(:)>1)/sum(C_T_col_tp{j}(:)>0);      %how much CY5 correlates with TMR
end
end
for j = 1:size(T_C_col_tp,2)
    subplot(2,4,j)

plot(xs,ys(1,:,j),'LineWidth',3,'Color','r')
hold on
plot(xs,ys(2,:,j),'LineWidth',3,'Color','b')
xlabel('Pixels Apart')
ylabel('Fraction of Spots')
ylim([0 1])
legend({'Fraction RepA with Colocalized Xist','Fraction Xist with Colocalized RepA'})
title(times0{j})
end

sum(T_T_col(:)==1)
sum(T_T_col(:)==2)
for i = 1:7
    i
%     sum(T_C_col_tp{i}(:)>1)
%     sum(T_C_col_tp{i}(:)>0)
sum(T_C_col_tp{i}(:)>1)/sum(T_C_col_tp{i}(:)>0)
sum(C_T_col_tp{i}(:)>1)/sum(C_T_col_tp{i}(:)>0)
end

%% Number of Xist Clouds over time

figure(71);clf
colors = {'k' 'b' 'g' 'm' 'r' 'y' 'c'}
times2 = [0,.25,.5,1,2,5,9];
threshes = [5:10:95];
for mm = 1:size(threshes,2)
    subplot(2,size(threshes,2)/2,mm)
cloud_thresh = threshes(mm);  %Number of transcripts to be determined an Xist cloud
cloud_num = 0;      %Matrix with cloud numebrs. First row is 0 clouds, second row is 1, third is two
for i = 1:size(times2,2)
    num_cloud_temp = [0];
    for j = 1:size(CY5_AF594_TMR_tp{i},1)
        num_cloud_temp(j,1) = sum(CY5_AF594_TMR_tp{i}(j,11:12)>=cloud_thresh);
    end
    for k = 0:2
        cloud_num(k+1,i) = sum(num_cloud_temp==k)/size(CY5_AF594_TMR_tp{i},1);
    end
end
for i = 1:3
    plot(times2,cloud_num(i,:),[colors{i} '-o'],'LineWidth',3)
    hold on
end
legend('Zero Clouds','One Cloud','Two Clouds')
xlabel('Time (Days)')
ylabel('Fraction of Cells')
set(gca,'FontSize',14)
xlim([0 9])
title([num2str(cloud_thresh) ' Transcripts'])
end

%% Distribution of Xist Clouds Over time
figure(72);clf
colors = {'k' 'b' 'g' 'm' 'r' 'y' 'c'}
times2 = [0,.25,.5,1,2,5,9];
threshes = 20;
bin_num = 15;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 300;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)


%plot(labelsX,binall,'r','LineWidth',3); hold on
title('Distribution for Xist Positive Cells')
xlabel('Xist Molecules')
ylabel('Number of Cells')
for mm = 1:size(threshes,2)

cloud_thresh = threshes(mm);  %Number of transcripts to be determined an Xist cloud
cloud_num = 0;      %Matrix with cloud numebrs. First row is 0 clouds, second row is 1, third is two
for i = 1:size(times2,2)
        subplot(2,4,i)
    temp1 = CY5_AF594_TMR_tp{i}(:,11:12);
    minV = cloud_thresh;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
    temp2 = temp1(temp1 >= minV);
    bin_size = (maxV-minV)/bin_num;
    binsX = minV:bin_size:maxV;
    labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
    binsall = zeros(1,size(binsX,2)-1);
%binsall = imhistc(temp2,bin_num,0,bin_num);
for l = 1:size(binsall,2)
    if l == size(binsall,2)
         binsall(1,l) = size(find(temp2 >= binsX(l)),1);
    else        
        binsall(1,l) = size(find(temp2 >= binsX(l)),1)-size(find(temp2 >= binsX(l+1)),1);
    end
end
binall = binsall/sum(binsall);  %normalize
bar(labelsX,binall)
xlabel('Time (Days)')
ylabel('Fraction of Cells')
set(gca,'FontSize',14)
ylim([0 1])
title(times0{i})
end
end
%% Determine fluorescence values of colocalized versus noncolocalized spots
C_T_col(:,1:2) = 0;             %Eliminate Clouds 
T_C_col(:,1:2) = 0;             %Eliminate Clouds
T_T_col(:,1:2) = 0;             %Eliminate Clouds
CY5_colpos = 0;         %How many cells have CY5 spots colocalized with TMR
CY5_colneg = 0;         %How many cells have CY5 spots not colocalized with TMR
TMR_colpos = 0;         %How many cells have TMR spots colocalized with CY5
TMR_colneg = 0;         %How many cells have TMR spots not colocalized with CY5
CY5_neg_only = 0;       %Number of cells that only have non-colocalized CY5 spots
TMR_neg_only = 0;       %Number of cells that only have non-colocalized TMR spots
CY5_TMR_colneg = 0;     %Number of cells that have both non-colocalized CY5 and non-colocalized TMR
for i = 1:size(C_T_col,1)           %Determine numbers of cells in each category
    if sum(C_T_col(i,:) == 2) > 0
        CY5_colpos = CY5_colpos+1;
    elseif sum(C_T_col(i,:) == 1) > 0 & sum(T_C_col(i,:) == 1) == 0 %if no colocalized spots (elseif), and with non-colocalized CY5, and no non-colocalized TMR
        CY5_neg_only = CY5_neg_only+1;
    elseif sum(C_T_col(i,:) == 1) > 0 & sum(T_C_col(i,:) == 1) > 0 %if no colocalized spots (elseif), with non colocalized CY5 and TMR
        CY5_TMR_colneg = CY5_TMR_colneg+1;
    end
    if sum(C_T_col(i,:) == 1) > 0
        CY5_colneg = CY5_colneg+1;
  %      if sum(
    end
    if sum(T_C_col(i,:) == 2) > 0
        TMR_colpos = TMR_colpos+1;
    elseif sum(T_C_col(i,:) == 1) > 0 & sum(C_T_col(i,:) == 1) == 0 %if no colocalized spots (elseif), and with non-colocalized TMR, and no non-colocalized CY5
        TMR_neg_only = TMR_neg_only+1;
    end
    if sum(T_C_col(i,:) == 1) > 0
        TMR_colneg = TMR_colneg+1;
    end
end
col_br_CY5 = immultiply(C_T_col == 2,CY5_int);
col_br_TMR = immultiply(T_C_col == 2,TMR_int);
col_br_TMR_TMR = immultiply(T_T_col == 2,TMR_int);
noncol_br_CY5 = immultiply(C_T_col == 1,CY5_int);
noncol_br_TMR = immultiply(T_C_col == 1,TMR_int);
col_br_CY5(col_br_CY5 == 0) = NaN;
col_br_TMR(col_br_TMR == 0) = NaN;
col_br_TMR_TMR(col_br_TMR_TMR == 0) = NaN;
noncol_br_CY5(noncol_br_CY5 == 0) = NaN;
noncol_br_TMR(noncol_br_TMR == 0) = NaN;
clear col_br_CY5_tp col_br_TMR_tp noncol_br_CY5_tp noncol_br_TMR_tp
col_br_CY5_tp = {};
col_br_TMR_tp = {};
for i = 1:size(CY5_AF594_TMR_tp,2)
    col_br_CY5_tp{i} = immultiply(C_T_col_tp{i} == 2,CY5_int_tp{i});
    col_br_TMR_tp{i} = immultiply(T_C_col_tp{i} == 2,TMR_int_tp{i});
    noncol_br_CY5_tp{i} = immultiply(C_T_col_tp{i} == 1,CY5_int_tp{i});
    noncol_br_TMR_tp{i} = immultiply(T_C_col_tp{i} == 1,TMR_int_tp{i});
    col_br_CY5_tp{i}(col_br_CY5_tp{i} == 0) = NaN;
    col_br_TMR_tp{i}(col_br_TMR_tp{i} == 0) = NaN;
    noncol_br_CY5_tp{i}(noncol_br_CY5_tp{i} == 0) = NaN;
    noncol_br_TMR_tp{i}(noncol_br_TMR_tp{i} == 0) = NaN;
end
figure(91); clf
colors = {'k' 'b' 'g' 'm' 'r' 'y' 'c'}
times2 = [0,.25,.5,1,2,5,9];
bin_num = 80;
maxV = 40000;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 0;
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
binsall = zeros(1,size(binsX,2)-1);
%binsall = imhistc(temp2,bin_num,0,bin_num);
for l = 1:size(binsall,2)
    if l == size(binsall,2)
        binsall(1,l) = sum(sum(col_br_CY5 >= binsX(l)));
        binsall(2,l) = sum(sum(col_br_TMR >= binsX(l)));
        binsall(3,l) = sum(sum(noncol_br_CY5 >= binsX(l)));
        binsall(4,l) = sum(sum(noncol_br_TMR >= binsX(l)));
        binsall(5,l) = sum(sum(col_br_TMR_TMR >= binsX(l)));
    else
        binsall(1,l) = sum(sum(col_br_CY5 >= binsX(l))-sum(col_br_CY5 >= binsX(l+1)));
        binsall(2,l) = sum(sum(col_br_TMR >= binsX(l))-sum(col_br_TMR >= binsX(l+1)));
        binsall(3,l) = sum(sum(noncol_br_CY5 >= binsX(l))-sum(noncol_br_CY5 >= binsX(l+1)));
        binsall(4,l) = sum(sum(noncol_br_TMR >= binsX(l))-sum(noncol_br_TMR >= binsX(l+1)));
        binsall(5,l) = sum(sum(col_br_TMR_TMR >= binsX(l))-sum(col_br_TMR_TMR >= binsX(l+1)));    
    end
end
for m = 1:size(binsall,1)
    binall(m,1:size(binsall,2)) = binsall(m,:)/sum(binsall(m,:));  %normalize
    plot(labelsX,binall(m,:),[colors{m} '-o'],'LineWidth',3)
    hold on
end    
xlabel('Integrated Fluorescence')
ylabel('Fraction of Cells')
set(gca,'FontSize',14)
%ylim([0 1])
legend(['Colocalized CY5 spots (' num2str(sum(C_T_col(:)==2)) ' spots), ' num2str(CY5_colpos) ' cells'],...
    ['Colocalized TMR Spots (' num2str(sum(T_C_col(:)==2)) ' spots), ' num2str(TMR_colpos) ' cells'],...
    ['Non-colocalized CY5 Spots (' num2str(sum(C_T_col(:)==1)) ' spots), ' num2str(CY5_colneg) ' cells'],...
    ['Non-colocalized TMR Spots (' num2str(sum(T_C_col(:)==1)) ' spots), ' num2str(TMR_colneg) ' cells'],'TMR Spots')

%%%Also show intensity of colocalized versus noncolocalized spots, but by timepoint
clear binsall binall
binsall = zeros(1,size(binsX,2)-1);
figure(92); clf
for i = 1:size(times2,2)
    subplot(2,4,i)
    C_T_col_tp{i}(:,1:2) = 0;             %Eliminate Clouds 
T_C_col_tp{i}(:,1:2) = 0;             %Eliminate Clouds
T_T_col_tp{i}(:,1:2) = 0;             %Eliminate Clouds
    CY5_colpos = 0;         %How many cells have CY5 spots colocalized with TMR
CY5_colneg = 0;         %How many cells have CY5 spots not colocalized with TMR
TMR_colpos = 0;
TMR_colneg = 0;
for j = 1:max([size(C_T_col_tp{i},1) size(T_C_col_tp{i},1)])
    if sum(C_T_col_tp{i}(j,:) == 2) > 0
        CY5_colpos = CY5_colpos+1;
    end
    if sum(C_T_col_tp{i}(j,:) == 1) > 0
        CY5_colneg = CY5_colneg+1;
    end
    if sum(T_C_col_tp{i}(j,:) == 2) > 0
        TMR_colpos = TMR_colpos+1;
    end
    if sum(T_C_col_tp{i}(j,:) == 1) > 0
        TMR_colneg = TMR_colneg+1;
    end
end
for l = 1:size(binsall,2)
    if l == size(binsall,2)
        binsall(1,l) = sum(sum(col_br_CY5_tp{i} >= binsX(l)));
        binsall(2,l) = sum(sum(col_br_TMR_tp{i} >= binsX(l)));
        binsall(3,l) = sum(sum(noncol_br_CY5_tp{i} >= binsX(l)));
        binsall(4,l) = sum(sum(noncol_br_TMR_tp{i} >= binsX(l)));
    else
        binsall(1,l) = sum(sum(col_br_CY5_tp{i} >= binsX(l))-sum(col_br_CY5_tp{i} >= binsX(l+1)));
        binsall(2,l) = sum(sum(col_br_TMR_tp{i} >= binsX(l))-sum(col_br_TMR_tp{i} >= binsX(l+1)));
        binsall(3,l) = sum(sum(noncol_br_CY5_tp{i} >= binsX(l))-sum(noncol_br_CY5_tp{i} >= binsX(l+1)));
        binsall(4,l) = sum(sum(noncol_br_TMR_tp{i} >= binsX(l))-sum(noncol_br_TMR_tp{i} >= binsX(l+1)));
    end
end
for m = 1:size(binsall,1)
    binall(m,1:size(binsall,2)) = binsall(m,:)/sum(binsall(m,:));  %normalize
    plot(labelsX,binall(m,:),[colors{m} '-o'],'LineWidth',2)
    hold on
end    
xlabel('Integrated Fluorescence')
ylabel('Fraction of Cells')
set(gca,'FontSize',11)
%ylim([0 1])
legend(['Colocalized CY5 spots (' num2str(sum(C_T_col_tp{i}(:)==2)) ' spots), ' num2str(CY5_colpos) ' cells'],...
    ['Colocalized TMR Spots (' num2str(sum(T_C_col_tp{i}(:)==2)) ' spots), ' num2str(TMR_colpos) ' cells'],...
    ['Non-colocalized CY5 Spots (' num2str(sum(C_T_col_tp{i}(:)==1)) ' spots), ' num2str(CY5_colneg) ' cells'],...
    ['Non-colocalized TMR Spots (' num2str(sum(T_C_col_tp{i}(:)==1)) ' spots), ' num2str(TMR_colneg) ' cells'],'TMR Spots')
title(times0{i})
end
clear binsall binall
%Scatterplot of intensities of colocalized spots
C_T_corr_int = zeros(2,3);          %This will have the intensity of the corresponding TMR spot for every colocalized CY5 spot (same index placement as rnas in other matrices)
CY5_TMR_corr_int_only = zeros(2,1);         %This will only have the intensities of CY5 spots in first column and corresponding intensities in TMR in the second column
counter = 1;            %counter for CY5_TMR_corr_int_only
for i = 1:size(CY5_AF594_TMR,1)
    if sum(C_T_col(i,:)==2)>0       %search if it has a colocalized spot
        init_ind = find(C_T_col(i,:)==2);   %Find indices of CY5 spots colocalized   
        for j = 1:size(init_ind,2)
            if C_T_col_ind(i,init_ind(j)) > 2       %Avoid clouds
                CY5_TMR_corr_int_only(counter,1) =  col_br_CY5(i,init_ind(j));
                CY5_TMR_corr_int_only(counter,2) =  col_br_TMR(i,C_T_col_ind(i,init_ind(j)));
                counter = counter+1;
            end
        end
    end
end
counter = 1
%for i = 25:50:275
i = 75;
figure(120+counter); clf; scatplot(CY5_TMR_corr_int_only(:,1),CY5_TMR_corr_int_only(:,2),'circles',i,100,5,3,15);
counter = counter+1;
ylim([0 60000])
xlim([0 60000])
xlabel('Xist intensity','Fontsize',16)
ylabel('RepA intensity','Fontsize',16)   
title('Integrated Intensities of Colocalized Spots','FontSize',16)
%end
CY5_TMR_corr_int_only(CY5_TMR_corr_int_only <1) = 1;
log_temp = log2(CY5_TMR_corr_int_only);

 figure(121); clf; scatplot(log_temp(:,1),log_temp(:,2),'circles',.7,100,5,3,15);   %changing the value after 'circles' is important for dynamic range
ylim([10 20])
xlim([10 20])
xlabel('log2(Xist intensity)','Fontsize',16)
ylabel('log2(RepA intensity)','Fontsize',16) 
title('Integrated Intensities of Colocalized Spots','FontSize',16)
%% Determine transcription sites within clouds (over time)
%This will look at the intensities in  CY5_in_CY5cloud,CY5_in_CY5cloud_tp,
%CY5_in_TMRcloud,CY5_in_TMRcloud_tp,TMR_in_CY5cloud,TMR_in_CY5cloud_tp,TMR_in_TMRcloud,and TMR_in_TMRcloud_tp 
%in order to determine transcription sites. Then it will use
%CY5_pos_init and related variables to determine the position

%This will look at all clouds with greater than or equal to the cloud
%thres and determine for cells with at least one cloud that size whether a
%transcription site exists for the TMR channel inside a cloud or not 
%(has to be as bright as tsite_thres*intensity of individual RNA molecule)
%%% Relevant numbers in CY5_AF594_TMR
%11: Number of transcripts at the first transcription site/cloud CY5
%12: Number of transcripts at the second transcription site/cloud CY5
%13: Number of transcripts at the first transcription site/cloud TMR
%14: Number of transcripts at the second transcription site/cloud TMR
cloud_thres = 3;   %Threshold for how big a cloud is
tsite_thres = 3;    %Threshold for brightness / brightness of individual transcript to be considered a tsite
cl_size_w_tsite = [0,0];    %Will have number of transcripts in CY5 clouds with TMR tsites
cl_size_no_tsite = [0,0];   %Will have number of transcripts in CY5 clouds with no TMR tsites
counter = 1;    %Counter for storing size of cloud that has tsite
counter2 = 1;   %counter for storing size of cloud that does not have tsite
num_both = [0,0];   %Will have the number 
num_tsite_cl = zeros(3,size(TMR_in_CY5cloud_tp,2)); %Number of tsites in each cloud per timepoint. First row is no cloud, second row is cloud 1, third row is cloud 2     
cl_cells = zeros(3,size(TMR_in_CY5cloud_tp,2));  %Number of cells with number of clouds per timepoint (first row is 0 clouds, second row is 1 cloud, first column is first timepoint, etc)
for i = 1:size(TMR_in_CY5cloud_tp,2)        %Go to each timepoint
    for m = 1:size(TMR_in_CY5cloud_tp{i},1) %Go to each cell in timepoint
        num_clouds_temp = sum(CY5_AF594_TMR_tp{i}(m,11:12) >= cloud_thres); %Number of clouds in this cell
        cl_cells(num_clouds_temp+1,i) = cl_cells(num_clouds_temp+1,i) + 1;  %Add cell to appropriate indicator depending on cloud num          
        %if num_clouds_temp>0    %only look at cells with at least one cloud (comment out if you want to look at all of them)
        cl_tsite_temp = zeros(3,1);  %Number of TMR tsites in no clouds (row1), cloud 1 (row2), cloud2 (row3)
        for k = 1:size(TMR_in_CY5cloud_tp{i},2) %Go through each transcript in cell          
            if TMR_int_tp{i}(m,k) >= TMR_cloudInt(i)*tsite_thres  %If a TMR transcription site
                
                if TMR_in_CY5cloud_tp{i}(m,k) == 4      %If transcription site is in both clouds
                    if cl_tsite_temp(2,1) == 0 & cl_tsite_temp(3,1) == 0  %Only add cloud transcript number if there were not already other clouds
                        cl_size_w_tsite(counter:counter+1) = CY5_AF594_TMR_tp{i}(m,11:12);  %Add both clouds as clouds with a transcription site
                        counter = counter+2;
                    elseif cl_tsite_temp(2,1) == 0  %Only add if another transcription site wasn't already added to the first
                        cl_size_w_tsite(counter) = CY5_AF594_TMR_tp{i}(m,11);   
                        counter = counter+1;
                    elseif cl_tsite_temp(3,1) == 0
                        cl_size_w_tsite(counter) = CY5_AF594_TMR_tp{i}(m,12);
                        counter = counter+1;
                    end
                    cl_tsite_temp(2,1) = cl_tsite_temp(2,1)+1;  %Add number of transcription sites in both clouds
                    cl_tsite_temp(3,1) = cl_tsite_temp(3,1)+1;  %Add number of transcription sites in both clouds
                elseif TMR_in_CY5cloud_tp{i}(m,k) == 3          %If tsite in second cloud
                    if cl_tsite_temp(3,1) == 0
                        cl_size_w_tsite(counter) = CY5_AF594_TMR_tp{i}(m,12);
                        counter = counter+1;
                    end                    
                    cl_tsite_temp(3,1) = cl_tsite_temp(3,1)+1;
                elseif TMR_in_CY5cloud_tp{i}(m,k) == 2  %If tsite in first cloud
                    if cl_tsite_temp(2,1) == 0 & cl_tsite_temp(3,1) == 0
                        cl_size_w_tsite(counter) = CY5_AF594_TMR_tp{i}(m,11);
                        counter = counter+1;
                    end                   
                    cl_tsite_temp(2,1) = cl_tsite_temp(2,1)+1;
                elseif TMR_in_CY5cloud_tp{i}(m,k) == 1
                    cl_tsite_temp(1,1) = cl_tsite_temp(1,1)+1;  %Number of TMR tsites in no clouds 
                end
            end
        end
        if CY5_AF594_TMR_tp{i}(m,11) < cloud_thres      %If the size of the first cloud is not big enough
            cl_tsite_temp(1,1) = cl_tsite_temp(1,1)+ cl_tsite_temp(2,1);   %Add it to no clouds
            cl_tsite_temp(2,1) = 0; %set tsites in first cloud to 0
        end
        if CY5_AF594_TMR_tp{i}(m,12) < cloud_thres      %If size of second cloud is not big enough
            cl_tsite_temp(1,1) = cl_tsite_temp(1,1)+ cl_tsite_temp(3,1);    %Add it to no clouds
            cl_tsite_temp(3,1) = 0; %set tsites in second cloud to 0
        end
        num_tsite_cl(1:3,i) = num_tsite_cl(1:3,i)+cl_tsite_temp;
        if cl_tsite_temp(2,1) == 0
            cl_size_no_tsite(1,counter2) = CY5_AF594_TMR_tp{i}(m,11);    %add cloud 1 to clouds with no tsite
            counter2 = counter2+1;
        end
        if cl_tsite_temp(3,1) == 0
            cl_size_no_tsite(1,counter2) = CY5_AF594_TMR_tp{i}(m,12);    %add cloud 2 to clouds with no tsite
            counter2 = counter2+1;
        end        
        if sum(cl_tsite_temp(2:3)>1) > 0
            ['Warning: there is a cloud with ' num2str(max(cl_tsite_temp(2:3))) ' TMR transcription sites inside']
        end
        %end  %comment out if you want to look at all cells (not just ones with clouds)
    end
end
num_tsite_cl
figure(126); clf
colors = {'k' 'b' 'm' 'r' 'g' 'c' [255 159 0]/(255+159)}
times2 = [0,.25,.5,1,2,5,9];
plot(times2,sum(num_tsite_cl,1),colors{1},'LineWidth',3)   %Plot total number of transcription sites
hold on
plot(times2,sum(num_tsite_cl(2:3,:),1),colors{2},'LineWidth',3)   %Plot total number of transcription sites inside clouds over time
plot(times2,sum(cl_cells(2:3,:),1),colors{3},'LineWidth',3)   %Plot total number of cells with clouds over time
plot(times2,sum(cl_cells(1,:),1),colors{4},'LineWidth',3)   %Plot total number of cells without clouds over time
legend('Total number of transcription sites','Number of transcription sites in clouds','Number of cells with clouds','Number of cells without clouds')
title('Pdk3 Transcription Sites in Clouds')
set(gca,'FontSize',14)
xlim([0 9])
xlabel('Time (Days)')
%title('Fraction of Pdk3 transcription sites over time within clouds, given at least one cloud')

%% Look at distance to the nucleus over time
times2 = [0,.25,.5,1,2,5,9];
mean_nuc_dist = zeros(1,size(Nuc_dist_tp,2));
mean_nuc_dist_c1 = zeros(1,size(Nuc_dist_tp,2));
mean_nuc_dist_c2 = zeros(1,size(Nuc_dist_tp,2));
mean_nuc_dist_noncloud = zeros(1,size(Nuc_dist_tp,2));
median_nuc_dist = zeros(1,size(Nuc_dist_tp,2));
median_nuc_dist_c1 = zeros(1,size(Nuc_dist_tp,2));
median_nuc_dist_c2 = zeros(1,size(Nuc_dist_tp,2));
median_nuc_dist_noncloud = zeros(1,size(Nuc_dist_tp,2));
for i = 1:size(Nuc_dist_tp,2)
Nuc_dist_tp{i}(Nuc_dist_tp{i} == 0) = NaN;
lin_dist = Nuc_dist_tp{i}(:,1:2);   %Extracts nuclear distance of both clouds so it can be linearized later
lin_dist1 = Nuc_dist_tp{i}(:,3:end);   %Extracts nuclear distance of all transcripts so it can be linearized later
mean_nuc_dist(1,i) = nanmean(lin_dist(:));    %Mean distance of both clouds
mean_nuc_dist_c1(1,i) = nanmean(Nuc_dist_tp{i}(:,1));  %Mean distance of first cloud
mean_nuc_dist_c2(1,i) = nanmean(Nuc_dist_tp{i}(:,2));   %Mean distance of second cloud  
mean_nuc_dist_noncloud(1,i) = nanmean(lin_dist1(:));
median_nuc_dist(1,i) = nanmedian(lin_dist(:));  %Median distance of both clouds
median_nuc_dist_c1(1,i) = nanmedian(Nuc_dist_tp{i}(:,1));  %Median distance of first cloud
median_nuc_dist_c2(1,i) = nanmedian(Nuc_dist_tp{i}(:,2));  %Median distance of second cloud
median_nuc_dist_noncloud(1,i) = nanmedian(lin_dist1(:));
end
figure(127); clf; hold on
plot(times2,mean_nuc_dist,colors{1},'LineWidth',3)
plot(times2,mean_nuc_dist_c1,colors{2},'LineWidth',3)
plot(times2,mean_nuc_dist_c2,colors{3},'LineWidth',3)
plot(times2,mean_nuc_dist_noncloud,colors{4},'LineWidth',3)
legend('Mean distance of clouds to the nucleus','Mean distance of first cloud to the nucleus','Mean distance of second cloud to the nucleus','Mean distance of non clouds to the nucleus')
title('Xist Clouds in Pdk3 dataset judged by brightest pixel')
xlabel('Time (Days)')
ylabel('Distance (pixels, z is adjusted by factor of 3')
set(gca,'FontSize',13)
xlim([0 9])
% lin_dist = Nuc_dist(:,1:2);
% ylim([0 max(Nuc_dist(:))])
figure(128); clf; hold on
plot(times2,median_nuc_dist,colors{1},'LineWidth',3)
plot(times2,median_nuc_dist_c1,colors{2},'LineWidth',3)
plot(times2,median_nuc_dist_c2,colors{3},'LineWidth',3)
plot(times2,median_nuc_dist_noncloud,colors{4},'LineWidth',3)
legend('median distance of clouds to the nucleus','median distance of first cloud to the nucleus','median distance of second cloud to the nucleus','median distance of non clouds to the nucleus')
title('Xist Clouds in Pdk3 dataset judged by brightest pixel')
xlabel('Time (Days)')
ylabel('Distance (pixels, z is adjusted by factor of 3')
set(gca,'FontSize',13)
xlim([0 9])

%% Look at cloud size (transcript number) and pdk3 transcripts versus distance to nucleus
lin_Nuc = Nuc_dist(:,1:2);
lin_cloud = CY5_AF594_TMR(:,11:12);
lin_Nuc_mindist = nanmin(Nuc_dist(:,1:2),[],2)
lin_Nuc_mindist(lin_Nuc_mindist == 0) = NaN;
figure(128);clf; scatplot(lin_Nuc(:),lin_cloud(:),'circles',15,100,5,3,15);
xlabel('Distance to Nucleus from Brightest Point (pixels, z-adjusted)')
ylabel('Number of Transcripts in Cloud')
set(gca,'FontSize',13)

figure(128);clf; scatplot(lin_Nuc_mindist',CY5_AF594_TMR(:,3)','circles',15,100,5,3,15);
xlabel('Min Distance to Nucleus from Brightest Point of Cloud (pixels, z-adjusted)')
ylabel('Number of Pdk3 Transcripts')
set(gca,'FontSize',13)
%scatter(lin_Nuc(:),lin_cloud(:))

%% Looking for transcription sites of TMR (RepA) that are not colocalized with an Xist spot
cell_nums = zeros(1,size(TMR_int_tp,2)); %Stores the number of cells per timepoint
tsite_nums = zeros(1,size(TMR_int_tp,2)); %Stores the number of transcription sites per timepoint
tsite_cells = zeros(1000,size(TMR_int_tp,2)); %Stores number of cells with number of of transcription sites. Row is how many transcription sites (+1), column is timepoint
tsite_col = zeros(1,size(TMR_int_tp,2)); %Stores the number of transcription sites that colocalizes with an Xist spot per timepoint
tsite_noncol = zeros(1,size(TMR_int_tp,2)); %Stores the number of transcription sites that does not colocalize with an Xist spot per timepoint
tsite_thres = 2.5;
for i = 1:size(TMR_int_tp,2)    %searches through every timepoint
    cell_nums(i) = size(TMR_int_tp{i},1);   %Stores how many cells per timepoint
    tsite_counter_tp = 0;   %will count how many transcription sites per timepoint
    tsite_col_temp = 0;     %will count how many transcription sites colocalize with CY5 per timepoint
    tsite_noncol_temp = 0;  %will count how many transcription sites DO NOT colocalize with CY5 per timepoint
    for j = 1:size(TMR_int_tp{i},1) %searches through every cell
        tsite_temp = sum(TMR_int_tp{i}(j,:) >= tsite_thres*TMR_cloudInt(i));
        tsite_cells(tsite_temp+1,i) = tsite_cells(tsite_temp+1,i)+1;
        tsite_counter_tp = tsite_counter_tp+tsite_temp;
        for k = 1:size(TMR_int_tp{i},2) %searches through every transcript
            if TMR_int_tp{i}(j,k) >= tsite_thres*TMR_cloudInt(i)    %if it is a transcription site
                %tsite_counter_tp = tsite_counter_tp +1; %Add to transcription site counter
                if T_C_col_tp{i}(j,k) == 2              %If there is a colocalized Xist
                    tsite_col_temp = tsite_col_temp+1;  %Add to counter of colocalized transcription sites
                elseif T_C_col_tp{i}(j,k) == 1              %If there is not a colocalized Xist
                    tsite_noncol_temp = tsite_noncol_temp+1; %Add to counter of noncolocalized transcription sites
                else
                    'There was an error. Transcription site determined where colocalization matrix said there was no transcript'
                end
            end
        end
    end
    tsite_col(i) = tsite_col_temp;
    tsite_noncol(i) = tsite_noncol_temp;
    tsite_nums(i) = tsite_counter_tp;
end
i = 5
j = 1
TMR_pos_init_tp_ts = immultiply(TMR_pos_init_tp{i},repmat(TMR_int_tp{i}>tsite_thres*TMR_cloudInt(i),[1,1,3]));
figure(900); clf;
scatter3(TMR_pos_init_tp{i}(j,:,2),TMR_pos_init_tp{i}(j,:,1),TMR_pos_init_tp{i}(j,:,3));
title('Localization of RepA spots in a single cell')
figure(901); clf; 
scatter3(TMR_pos_init_tp_ts(j,:,2),TMR_pos_init_tp_ts(j,:,1),TMR_pos_init_tp_ts(j,:,3));
title('Localization of RepA transcription sites in a single cell')
figure(902); clf;
scatter(TMR_pos_init_tp{i}(j,:,2),-1*TMR_pos_init_tp{i}(j,:,1));
title('Localization of RepA spots in a single cell')
xlim([0 400])
ylim([-400 0])
figure(903); clf; 
scatter(TMR_pos_init_tp_ts(j,:,2),-1*TMR_pos_init_tp_ts(j,:,1));
title('Localization of RepA transcription sites in a single cell')
xlim([0 400])
ylim([-400 0])
['timepoint ' num2str(CY5_AF594_TMR_tp{i}(j,4))]
['Exp date ' num2str(CY5_AF594_TMR_tp{i}(j,5))]
['image number ' num2str(CY5_AF594_TMR_tp{i}(j,6))]
 i = 5
 
        
                
                    
                
        
        