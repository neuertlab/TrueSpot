%% Define Variables
%%For 20170525 and 20170531
for tran_thresh = 4 %1.5:.5:4
tic
dates0 = {'20190216','20190217','20190218','20190219','20190220','20190221'};
times0 = {{'0d'} {'.5d'} {'1d'} {'2d'} {'2d'},{'3d'}}; 
strain = 'F1-2-1';
times1 = [0,.5,1,2,2,3];
%%%
max_imnum = 50;
max_tps = 1;
%%% Channel Names. Make empty if channel not used
CY5_name = {'_Xist-CY5-th' '_Xist-CY5-th' '_Xist-CY5-th' ...
    '_Xist-CY5-th' '_Xist-CY5-th' '_Xist-CY5-th' '_Xist-CY5-th' '_Xist-CY5-th'};
AF594_name = {};
% TMR_name = {'_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th' '_Pdk3-TMR-th'};
% TMR_name = {'_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th' '_Maoa-TMR-th'};
TMR_name = {'_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' '_Tsix-TMR-th' };
% TMR_name = {'_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' '_RepA-TMR-th' };  
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
thNum_CY5 = 4; %Which set of the above thresholds to use (which column) %for 20190125 timecourse, 4
thNum_TMR = 4; %Which set of the above thresholds to use (which column) %for 20190125 timecourse, 4
thNum_AF594 = 3; %Which set of the above thresholds to use (which column)

%%% Whether to do Cloud quant or not
CY5_hascloud = 1;
AF594_hascloud = 0;
TMR_hascloud = 0;
%%%
CY5_cloudInt = [10000,10000,10000,...
    10000,10000,10000,10000,10000];
AF594_cloudInt = [];
%  TMR_cloudInt = [12550,12550,12550,12550,12550,12550,12550,12550,12550];   %Xist/RepA
%  TMR_cloudInt = [6750,6750,6750,6750,6750,6750,6750,6750,6750];    %Maoa
TMR_cloudInt = [8650,8650,8650,8650,8650,8650,8650,8650];  %Tsix
% TMR_cloudInt = [8950,8950,8950,8950,8950,8950,8950,8950,8950,8950,8950];    %Pdk3
%%%
nucTH = 'mid'
same_dist = 6;   %How close transcripts have to be before they are considered the same transcript
dist_col = 7;
CY5_tran_thres = tran_thresh ; %How many transcripts needed to be determined as a transcription site or cloud
TMR_tran_thres = tran_thresh ; %How many transcripts needed to be determined as a transcription site or cloud
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
    = A1B_RNA_Quant_General_2(...
    strain,dates0,times0,times1,max_tps,max_imnum,...
    thCY5all,thAF594all,thTMRall,CY5_name,AF594_name,TMR_name,thNum_CY5,thNum_AF594,thNum_TMR,...
    CY5_hascloud,AF594_hascloud,TMR_hascloud,CY5_cloudInt,AF594_cloudInt,TMR_cloudInt,nucTH,...
    dist_col,CY5_tran_thres,TMR_tran_thres); toc

save([CY5_name{1} TMR_name{1} '_' date '_3zadjust_6pixsame' num2str(TMR_tran_thres) 'transcriptsInTSite' '_CY5th' num2str(thNum_CY5)  '_TMRth' num2str(thNum_TMR)  ],'CY5_AF594_TMR','CY5_AF594_TMR_tp',...
    'CY5_in_CY5cloud','CY5_in_CY5cloud_tp','CY5_in_TMRcloud','CY5_in_TMRcloud_tp','TMR_in_CY5cloud','TMR_in_CY5cloud_tp','TMR_in_TMRcloud','TMR_in_TMRcloud_tp',...
    'CY5_pos_init','CY5_pos_init_tp','TMR_pos_init','TMR_pos_init_tp',...
    'T_C_col','C_T_col','T_C_col_tp','C_T_col_tp','T_T_col','T_T_col_tp'...
    ,'CY5_int','CY5_int_tp','TMR_int','TMR_int_tp','C_T_col_ind','C_T_col_ind_tp','T_C_col_ind','T_C_col_ind_tp',...
        'CY5_Nuc_dist','CY5_Nuc_dist_tp','TMR_Nuc_dist','TMR_Nuc_dist_tp','Cl_Nuc_dist','Cl_Nuc_dist_tp',...
    'C1_CY5_dist','C1_CY5_dist_tp','C2_CY5_dist','C2_CY5_dist_tp','C1_TMR_dist','C1_TMR_dist_tp','C2_TMR_dist','C2_TMR_dist_tp',...
    'thNum_CY5','thNum_AF594','thNum_TMR','CY5_tran_thres','TMR_tran_thres')
end
%% Correlate Xist cloud to transcripts and transcription sites

%%%Correlate largest cloud size through all timepoints
figure(101); clf;
CY5_x = 1;
CY5_label = 'Xist';
TMR_label = 'Tsix';
clims=[0,.3];
large_cloud = zeros(size(CY5_AF594_TMR,1),1);        %This will store the largest cloud for each cell in the timecourse
for i = 1:size(CY5_AF594_TMR,1)
    large_cloud(i,1) = max(CY5_AF594_TMR(i,11:12));
end
large_cloud = log2(large_cloud);    %Use the second log
%CY5_AF594_TMR(find(CY5_AF594_TMR(:,1) == max(CY5_AF594_TMR(:,1))),2:6)
bin_num = 5;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 9;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
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
        xlabel(['log2 ' CY5_label ' Transcripts in Cloud'])
        ylabel([TMR_label ' Transcription Sites'])
    else
        xlabel([TMR_label ' Transcription Sites'])
        ylabel(['log2 ' CY5_label ' Transcripts in Cloud'])
    end
    yticks([1 2 3])
    yticklabels({'2','1','0'})
%     xticklabels(labelsX)
    colorbar
    h = colorbar;
    title(h, 'Probability') 

%%%Separate by timepoint

figure(102); clf;
CY5_label = 'Xist';
TMR_label = 'Tsix';
clims=[0,.3];
for k = 1:size(CY5_AF594_TMR_tp,2)
large_cloud = zeros(size(CY5_AF594_TMR_tp{k},1),1);        %This will store the largest cloud for each cell in the timecourse
for i = 1:size(CY5_AF594_TMR_tp{k},1)
    large_cloud(i,1) = max(CY5_AF594_TMR_tp{k}(i,11:12));
end
large_cloud = log2(large_cloud);
%CY5_AF594_TMR_tp{k}(find(CY5_AF594_TMR_tp{k}(:,1) == max(CY5_AF594_TMR_tp{k}(:,1))),2:6)
bin_num = 5;
%maxV = max(CY5_AF594_TMR_tp{k}(:,1));
maxV = 10;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
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

figure(50); clf; gscatter(CY5_AF594_TMR(:,3),CY5_AF594_TMR(:,1),CY5_AF594_TMR(:,5),[],'o',10);
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
CY5_x = 0;      %Make 1 if you want the CY5 channel on the x axis
CY5_label = 'Xist';
TMR_label = 'Tsix';
clims=[0,.5]
figure(2001); clf;

p_vals = zeros(1,size(CY5_AF594_TMR_tp,2))
transc_corr = zeros(3,3,size(CY5_AF594_TMR_tp,2));     %will store how many cells are in each bin 
for i = 1:4%size(CY5_AF594_TMR_tp,2)
    temp_CY = CY5_AF594_TMR_tp{i}(:,8);     %Number of transcription sites for CY5
    temp_TMR = CY5_AF594_TMR_tp{i}(:,10);   %Number of transcription sites for TMR
    temp_CY(temp_CY>2) = 2;     %Make maximum sites = 2
    temp_TMR(temp_TMR>2) = 2;     %Make maximum sites = 2
    [blah,p_val] = corrcoef(CY5_AF594_TMR_tp{i}(:,8),CY5_AF594_TMR_tp{i}(:,10))
    p_vals(i) = p_val(1,2);
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
    subplot(1,4,i)
    imagesc(temp_tc,clims)
    title([times0{i}{1} ' (' num2str(size(CY5_AF594_TMR_tp{i},1)) ' cells)'])
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
    for k = 1:3
        for l = 1:3
            text(k,l,num2str(round(temp_tc(l,k)*100)/100),'HorizontalAlignment','center','Color','w')
        end
    end
    if false;%i == 7
    h = colorbar;
    title(h, 'Probability') 
    end
end
transc_corr 

figure(2002); clf
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
temp1 = log2(CY5_AF594_TMR(:,1));
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
bin_num = 30;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 30;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 0;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
temp1 = CY5_AF594_TMR(:,3);
temp2 = temp1(temp1 >= minV);
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1);%+bin_size/2;
binsall = {};
for j = 1:4
temp1 = CY5_AF594_TMR_tp{j}(:,3);
temp2 = temp1(temp1 >= minV);
binsall{j} = zeros(1,size(binsX,2)-1);
%binsall = imhistc(temp2,bin_num,0,maxV);
for i = 1:size(binsall{j},2)
    if i == size(binsall{j},2)
     binsall{j}(1,i) = size(find(temp2 >= binsX(i)),1)
    else
     binsall{j}(1,i) = size(find(temp2 >= binsX(i)),1)-size(find(temp2 >= binsX(i+1)),1);
    end
end
binall = binsall{j}/sum(binsall{j});  %normalize
plot(labelsX,binall,'Color',colors{j},'LineWidth',3); hold on
end
%binall = binsall/sum(binsall);  %normalize
%bar(labelsX,binsall)
%plot(labelsX,binall,'g','LineWidth',3); hold on
legend('0d', '12hr', '1d', '2d', '5d','9d')
%title('Distribution for Maoa Positive Cells')
xlabel('Tsix Molecules')
%ylabel('Number of Cells')
ylabel('Probability')
%end
set(gca,'FontSize',14)
ylim([0 .3])
xlim([0 20])

figure(25); clf;
colors = {'k' 'b' 'm' 'r' 'g' 'c' [255 159 0]/(255+159)}
bin_num = 15;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 9;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = .5;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
temp1 = CY5_AF594_TMR(:,1);
% maxV = max(temp1);
temp1(temp1 < 1) = 1;
temp2 = log2(temp1);
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
xlabel('log2 Xist Molecules')
%ylabel('Number of Cells')
ylabel('Probability')
%end



%% Boxplots:
%%%Boxplots with Xist on x axis (if box_col = 1) or tmr on x axis (box_col
%%%=3)
clear dists1 X_labels
x_label = 'Tsix';
y_label = 'Xist';
y_col = 1;          %column used for values on y axis
% CY5_AF594_TMR(:,18) = log(CY5_AF594_TMR(:,1));
% for i = 1:size(CY5_AF594_TMR_tp,2)
%     CY5_AF594_TMR_tp{i}(:,18) = log(CY5_AF594_TMR_tp{i}(:,y_col));
% end
% y_col = 18;          %column used for values on y axis
box_col = 3;        %column used for values on x axis
dists1 = {}; %store into individual matrices the cells for each bin (not used)
bin_num = 9;
maxV = 9;         %cutoff for max value shown (x axis). Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
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
%%%Populate 8th column with which bin each cell falls into
for i = 1:bin_num+1
    dists1{i} = zeros(2,size(CY5_AF594_TMR,2)-1);
    counter1 = 1;
    for j = 1:size(CY5_AF594_TMR,1)
        if i == bin_num+1
            if CY5_AF594_TMR(j,box_col) >= binsX(i)
           CY5_AF594_TMR2(j,15) = i-1;
           counter1 = counter1+1;            
            end
       elseif CY5_AF594_TMR(j,box_col) >= binsX(i) & CY5_AF594_TMR(j,box_col) < binsX(i+1)
           CY5_AF594_TMR2(j,15) = i-1;
           counter1 = counter1+1;
        end
    end
end
labels_here = unique(CY5_AF594_TMR2(:,15));
true_end = size(CY5_AF594_TMR2,1);    %The size of the matrix before fake values are added
%%% Add fake cells as placeholders for bins where there are no cells and determine how many cells are in each bin
cell_nums1 = zeros(1,size(X_labels,2));
for m = 0:size(X_labels,2)-1
    cell_nums1(1,m+1) = sum(CY5_AF594_TMR2(:,15) == m);  %Determine cell number
    if sum(labels_here == m) == 0
        CY5_AF594_TMR2(size(CY5_AF594_TMR2,1)+1,15) = m;
    end
end
%%%Boxplot of TMR on y axis   
figure(24); clf; 
h = boxplot(CY5_AF594_TMR2(:,y_col),CY5_AF594_TMR2(:,15),'Labels',X_labels)
%ylim([0 100])
xlabel([x_label ' per Cell'],'Fontsize',10)
ylabel([y_label ' per Cell'],'Fontsize',10)
hold on
plot(0:bin_num+2,mean(CY5_AF594_TMR2(1:true_end,y_col)):.0001:mean(CY5_AF594_TMR2(1:true_end,y_col))+.0001*(bin_num+2),'--g','LineWidth',2)
xtickangle(45)
yyaxis right
plot(1:bin_num+1,cell_nums1,'--','LineWidth',2);
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
set(h,'LineWidth',2)
ylabel('Cell Number')

figure(28);clf
h = boxplot(CY5_AF594_TMR2(:,10),CY5_AF594_TMR2(:,15),'Labels',X_labels)
xlabel([x_label ' per Cell'],'Fontsize',10)
ylabel(['Tsix Transcription Sites per Cell'],'Fontsize',10)
ylim([0 6])   
xtickangle(45)
hold on
% yyaxis right
% plot(1:bin_num+1,cell_nums1,'--','LineWidth',2);
% %ylabel('Cell Number')
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
set(h,'LineWidth',2)
yyaxis right
plot(1:bin_num+1,cell_nums1,'--','LineWidth',2);
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
set(h,'LineWidth',2)
ylabel('Cell Number')

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
% for i = 1:size(CY5_AF594_TMR_tp2,2)
%     CY5_AF594_TMR_tp2{i}(:,1) = log(CY5_AF594_TMR_tp2{i}(:,1));
% end
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
           CY5_AF594_TMR_tp2{k}(j,15) = i-1;
           counter1 = counter1+1;            
            end
       elseif CY5_AF594_TMR_tp2{k}(j,box_col) >= binsX(i) & CY5_AF594_TMR_tp2{k}(j,box_col) < binsX(i+1)
           CY5_AF594_TMR_tp2{k}(j,15) = i-1;
           counter1 = counter1+1;
        end
    end
end
labels_here = unique(CY5_AF594_TMR_tp2{k}(:,15));
%%% Add fake cells as placeholders for bins where there are no cells and determine how many cells are in each bin
true_end = size(CY5_AF594_TMR_tp2{k},1);    %The size of the matrix before fake values are added
cell_nums1 = zeros(1,size(X_labels,2));
for m = 0:size(X_labels,2)-1
    cell_nums1(1,m+1) = sum(CY5_AF594_TMR_tp2{k}(:,15) == m);  %Determine cell number    
    if sum(labels_here == m) == 0
        CY5_AF594_TMR_tp2{k}(size(CY5_AF594_TMR_tp2{k},1)+1,15) = m;
    end
end
        
%%%Boxplot 
figure(25)
subplot(1,2,3-k)
%subplot('Position',[left bottom width height])
% figure(30+k); clf

h = boxplot(CY5_AF594_TMR_tp2{k}(:,y_col),CY5_AF594_TMR_tp2{k}(:,15),'Labels',X_labels)
xlabel([x_label ' per Cell'],'Fontsize',10)
ylabel(['log(' y_label ') per Cell'],'Fontsize',10)
title(times0{k})
if k <=4
ylim([0 20])   
else
ylim([0 max(CY5_AF594_TMR(:,y_col))])
end
%ylim([0 100])
xtickangle(45)
hold on
plot(0:bin_num+2,mean(CY5_AF594_TMR_tp2{k}(1:true_end,y_col)):.0001:mean(CY5_AF594_TMR_tp2{k}(1:true_end,y_col))+.0001*(bin_num+2),':m','LineWidth',2)
yyaxis right
plot(1:bin_num+1,cell_nums1,'--','LineWidth',2);
ylabel('Cell Number')
set(gca,'FontSize',13)
set(gca,'LineWidth',2)
set(h,'LineWidth',2)
%%% boxplots for number of transcription sites
figure(26);
subplot(1,2,k)
%subplot('Position',[left bottom width height])
h = boxplot(CY5_AF594_TMR_tp2{k}(:,10),CY5_AF594_TMR_tp2{k}(:,15),'Labels',X_labels)
xlabel([x_label ' per Cell'],'Fontsize',10)
ylabel(['Tsix Transcription Sites per Cell'],'Fontsize',10)
title(times0{k})
% if k <=4
 ylim([0 6])   
% else
% ylim([0 max(CY5_AF594_TMR(:,y_col))])
% end
%ylim([0 100])
xtickangle(45)
hold on
plot(0:bin_num+2,mean(CY5_AF594_TMR_tp2{k}(1:true_end,y_col)):.0001:mean(CY5_AF594_TMR_tp2{k}(1:true_end,y_col))+.0001*(bin_num+2),'--m','LineWidth',2)
yyaxis right
plot(1:bin_num+1,cell_nums1,'--m','LineWidth',2);
ylabel('Cell Number')
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
set(h,'LineWidth',2)
%%%Boxplot for transcription sites (on x axis)
figure(27);
subplot(1,2,k)
%subplot('Position',[left bottom width height])
h = boxplot(CY5_AF594_TMR_tp2{k}(:,3),CY5_AF594_TMR_tp2{k}(:,10));%'Labels',X_labels)
ylabel([x_label ' per Cell'],'Fontsize',10)
xlabel(['Tsix Transcription Sites per Cell'],'Fontsize',10)
title(times0{k})
 xlim([0 6]) 
% if k <=4
% ylim([0 20])   
% else
% ylim([0 max(CY5_AF594_TMR(:,y_col))])
% end
xtickangle(45)
hold on
% plot(0:bin_num+2,mean(CY5_AF594_TMR_tp2{k}(1:true_end,y_col)):.0001:mean(CY5_AF594_TMR_tp2{k}(1:true_end,y_col))+.0001*(bin_num+2),'--r','LineWidth',2)
% yyaxis right
% plot(1:bin_num+1,cell_nums1,'--','LineWidth',2);
% ylabel('Cell Number')
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

% times = [0,.25,.5,1,2,5,9];
times = [0,.5,1,2];
avgs = zeros(4,1);
for i = 1:4
    avgs(i) = mean(CY5_AF594_TMR_tp{i}(:,3));
end
% avgs = avgs/max(avgs(:));
figure(44); clf; 
yyaxis left
plot(times,avgs,'-o','LineWidth',3)
% title('Average Pdk3 Over Time')
xlabel('Time (Days)')
ylabel('Mean Tsix')
% xlim([0 max(times)])
ylim([0 max(avgs)+3])
set(gca,'FontSize',18)
hold on

% figure(44); clf; 
yyaxis right
% times = [0,.25,.5,1,2,5,9];
times = [0,.5,1,2];
avgs = zeros(4,1);
for i = 1:4
    avgs(i) = mean(CY5_AF594_TMR_tp{i}(:,1));
end
%figure(44); clf; 
% avgs = avgs/max(avgs(:));
plot(times,avgs,'-o','LineWidth',3); hold on
% title('Normalized Average Expression Over Time','FontSize',18)
xlabel('Time (Days)')
ylabel('Mean Xist')
ylim([0 max(avgs)+5])
% ylim([0 1])
% legend('Tsix','Xist')
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
figure(44); clf; 
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
ylim([0 1])
set(gca,'FontSize',14)
% legend('Normalized Average Xist',['On Fraction Xist (above ' num2str(thresh1) ')'],['Fraction of Cells with at least One Cloud (' num2str(thresh1) ' transcripts)'],['Fraction of Cells with at least Two Clouds (' num2str(thresh1) ' transcripts)']  )
legend(['Fraction of Cells with at least One Cloud (' num2str(thresh1) ' transcripts)'],['Fraction of Cells with at least Two Clouds (' num2str(thresh1) ' transcripts)'])
%% Joint probability
%Throughout timecourse
T1_T3 = zeros(1,2);
x_col = 3;
y_col = 1;
log_y = 1;
log_x = 0
T1_T3(1:size(CY5_AF594_TMR(:,1),1),1)= CY5_AF594_TMR(:,x_col);
T1_T3(1:size(CY5_AF594_TMR(:,1),1),2)= CY5_AF594_TMR(:,y_col);
for i = 1:size(T1_T3,1)
    if log_y
        if T1_T3(i,2) < 1
            T1_T3(i,2) = 1;
        end
        T1_T3(i,2) = log2(T1_T3(i,2));
    end
    if log_x
        if T1_T3(i,1) < 1
            T1_T3(i,1) = 1;
        end
        T1_T3(i,1) = log2(T1_T3(i,2));
    end
end
figure(80); clf
hist3(T1_T3,'CdataMode','auto')
xlabel('Tsix')
ylabel('log2 Xist')
h = colorbar;
title(h, 'Probability') 
view(2)

%For specific timepoint
figure(81); clf
for j = 1:4%size(CY5_AF594_TMR_tp,2)
T1_T3 = zeros(1,2);
x_col = 3;
y_col = 1;
log_y = 1;
log_x = 0;
T1_T3(1:size(CY5_AF594_TMR_tp{j}(:,1),1),1)= CY5_AF594_TMR_tp{j}(:,x_col);
T1_T3(1:size(CY5_AF594_TMR_tp{j}(:,1),1),2)= CY5_AF594_TMR_tp{j}(:,y_col);
for i = 1:size(T1_T3,1)
    if log_y
        if T1_T3(i,2) < 1
            T1_T3(i,2) = 1;
        end
        T1_T3(i,2) = log2(T1_T3(i,2));
    end
    if log_x
        if T1_T3(i,1) < 1
            T1_T3(i,1) = 1;
        end
        T1_T3(i,1) = log2(T1_T3(i,2));
    end
end
subplot(1,4,j)
x_edges = [0:1:8];
y_edges = [0:8];
h1 = hist3(T1_T3,'CdataMode','auto','Edges',{x_edges y_edges});
 h1 = log(h1); h1(h1<0) = 0;
h1 = h1/sum(h1(:))
g = imagesc(flipud(h1'),[0 .1])
 yticklabels(fliplr(0:8));
 xticks(1:9)
 xticklabels(0:8);
 
initx_siz = (max(x_edges)-min(x_edges))/(size(x_edges,2)-1);
inity_siz = (max(y_edges)-min(y_edges))/(size(y_edges,2)-1);
% yticks1 =-.5:inity_siz*2:size(y_edges,2)-.5;
%  yticks(yticks1)
%  xticks1 =-.5:initx_siz*2:size(x_edges,2)-.5;
%  xticks(xticks1)
%  ylabels1 = {}
%  for k = 1:size(yticks1,2)-1
%      ylabels1{k} = num2str(size(y_edges,2)+1-(2*k));
%  end
%     yticklabels(ylabels1)
%      xlabels1 = {}
%      for k = 1:size(xticks1,2)-1
%      xlabels1{k} = num2str(size(x_edges,2)+1-(2*k));
%      end
%      xticklabels(xlabels1)
    colorbar
    h = colorbar;
    title(h, 'Log(Probability)') 
%   hN1 = hist3(mRNA_C1,'Edges',{RNA1 RNA2});                               % Compute joint probability distribution for RNA1 and RNA2
%     JointOut.low.nuc(:,:,i) = hN1./sum(hN1(:));                             % Compute normalized joint probability distribution for RNA1 and RNA2
%     SingleOut1.low.nuc(:,i) = sum(JointOut.low.nuc(:,:,i),2);               % Compute normalized marginal probability distribution for RNA1 from the joint probability distribution
%     SingleOut2.low.nuc(:,i) = sum(JointOut.low.nuc(:,:,i),1);               % Compute normalized marginal probability distribution for RNA2 from the joint probability distribution
%     Cor.low.nuc(i,1) = corr2(mRNA_C1(1:NoCells1(i),1),mRNA_C1(1:NoCells1(i),2)); % Compute pearson correlation coefficent of expression for RNA1 and RNA2 
%     hN2 = hist(mRNA_C1(:,1),RNA1);                                   % Compute normalized marginal probability distribution for RNA1 
%     SingleOut3.low.nuc(:,i) = hN2/sum(hN2(:));
%     hN2 = hist(mRNA_C1(:,2),RNA2);                                   % Compute normalized marginal probability distribution for RNA1
%     SingleOut4.low.nuc(:,i) = hN2/sum(hN2(:));                              

% h1 = histogram2(T1_T3(:,1),T1_T3(:,2),'DisplayStyle','tile','ShowEmptyBins','on','Normalization', 'probability')
% Xedges = h1.XBinEdges;
% Yedges = h1.YBinEdges;
% probs = h1.Values;
% probs(probs == 0) = 1;
% probs = log2(probs);
% 
% histogram2('XBinEdges',Xedges,'YBinEdges',Yedges,'BinCounts',probs)


title(times0{j})
xlabel('Tsix','FontSize',14)
ylabel('log2 Xist','FontSize',14)
h = colorbar;
title(h, 'Log(Probability)') 
set(gca,'FontSize',14)
view(2)
end


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
set(gca,'FontSize',14)


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
% C_T_col(:,1:2) = 0;             %Eliminate Clouds 
% T_C_col(:,1:2) = 0;             %Eliminate Clouds
% T_T_col(:,1:2) = 0;             %Eliminate Clouds
if TMR_hascloud 
T_C_col(TMR_in_CY5cloud == 2) = 0;
for i = 1:4
    T_C_col_tp{i}(TMR_in_CY5cloud_tp{i} == 2) = 0;
end
end
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
for i = 1:4%size(CY5_AF594_TMR_tp,2)
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
bin_num = 100;
maxV = 120000;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 0;
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
binsall = zeros(1,size(binsX,2)-1);
binall = zeros(1,size(binsX,2)-1);
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
legend(['Colocalized Xist spots (' num2str(sum(C_T_col(:)==2)) ' spots), ' num2str(CY5_colpos) ' cells'],...
    ['Colocalized RepA Spots (' num2str(sum(T_C_col(:)==2)) ' spots), ' num2str(TMR_colpos) ' cells'],...
    ['Non-colocalized Xist Spots (' num2str(sum(C_T_col(:)==1)) ' spots), ' num2str(CY5_colneg) ' cells'],...
    ['Non-colocalized RepA Spots (' num2str(sum(T_C_col(:)==1)) ' spots), ' num2str(TMR_colneg) ' cells'],'TMR Spots')

%%%Also show intensity of colocalized versus noncolocalized spots, but by timepoint
clear binsall binall
binsall = zeros(1,size(binsX,2)-1);
figure(92); clf
bin_num = bin_num/5;
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
binsall = zeros(1,size(binsX,2)-1);
binall = zeros(1,size(binsX,2)-1);
%binsall = imhistc(temp2,bin_num,0,bin_num);
for i = 1:4%size(times2,2)
    subplot(2,4,i)
%     C_T_col_tp{i}(:,1:2) = 0;             %Eliminate Clouds 
% T_C_col_tp{i}(:,1:2) = 0;             %Eliminate Clouds
% T_T_col_tp{i}(:,1:2) = 0;             %Eliminate Clouds
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
legend(['Col Xist (' num2str(sum(C_T_col_tp{i}(:)==2)) ' spots), ' num2str(CY5_colpos) ' cells'],...
    ['Col RepA (' num2str(sum(T_C_col_tp{i}(:)==2)) ' spots), ' num2str(TMR_colpos) ' cells'],...
    ['Non-col Xist(' num2str(sum(C_T_col_tp{i}(:)==1)) ' spots), ' num2str(CY5_colneg) ' cells'],...
    ['Non-col RepA (' num2str(sum(T_C_col_tp{i}(:)==1)) ' spots), ' num2str(TMR_colneg) ' cells'],'TMR Spots')
title(times0{i})
end
clear binsall binall
%% Scatterplot of intensities of colocalized spots
C_T_corr_int = zeros(2,3);          %This will have the intensity of the corresponding TMR spot for every colocalized CY5 spot (same index placement as rnas in other matrices)
CY5_TMR_corr_int_only = zeros(2,1);         %This will only have the intensities of CY5 spots in first column and corresponding intensities in TMR in the second column
counter = 1;            %counter for CY5_TMR_corr_int_only
CY5_int_noncloud = CY5_int(:,3:end);
for i = 1:size(CY5_AF594_TMR,1)
    if sum(C_T_col(i,:)==2)>0       %search if it has a colocalized spot
        init_ind = find(C_T_col(i,:)==2);   %Find indices of CY5 spots colocalized  
        noncol_ind = find(C_T_col(i,:)==1);   %Find indices of CY5 spots not colocalized  
        noncol_ind1 = find(T_C_col(i,:)==1);   %Find indices of CY5 spots not colocalized  
        for j = 1:size(init_ind,2)
%             if C_T_col_ind(i,init_ind(j)) > 2       %Avoid clouds
                CY5_TMR_corr_int_only(counter,1) =  CY5_int(i,init_ind(j));  %intensity of CY5
                CY5_TMR_corr_int_only(counter,2) =  TMR_int(i,C_T_col_ind(i,init_ind(j)));
                counter = counter+1;
%             end
        end
%         %%%Noncolocalized spots
%         for j = 1:size(noncol_ind,2)                %noncolocalized CY5 spots
% %             if C_T_col_ind(i,init_ind(j)) > 2       %Avoid clouds
%                 CY5_TMR_corr_int_only(counter,1) =  CY5_int(i,noncol_ind(j));  %intensity of CY5
%                 CY5_TMR_corr_int_only(counter,2) =  min(TMR_int(TMR_int>0));    %Set equal to minimum TMR intensity
%                 counter = counter+1;
% %            end
%         end
%         for j = 1:size(noncol_ind1,2)             %noncolocalized TMR spots
% %             if C_T_col_ind(i,init_ind(j)) > 2       %Avoid clouds
%                 CY5_TMR_corr_int_only(counter,2) =  TMR_int(i,noncol_ind1(j));  %intensity of TMR
%                 CY5_TMR_corr_int_only(counter,1) =  min(CY5_int_noncloud(CY5_int_noncloud>0));    %Set equal to minimum CY5 intensity
%                 counter = counter+1;
% %            end
%         end
%         %%%
%         %%%Do this to include transcription sites not colocalized        
%         tssite_ind = find(TMR_int(i,:) >= 3*TMR_cloudInt(1)); %threshold for number of transcripts (second is greater than)
%         for j = 1:size(tssite_ind,2)
%             CY5_TMR_corr_int_only(counter,1) =  1;
%             CY5_TMR_corr_int_only(counter,2) =  TMR_int(i,tssite_ind(j));
%             counter = counter+1;
%         end
    end
end
counter = 1
for i = 4000/9000
% i = 75;
% figure(122); clf; scatplot(log(CY5_TMR_corr_int_only(:,1)/CY5_cloudInt(1)),log(CY5_TMR_corr_int_only(:,2)/TMR_cloudInt(1)),'circles',i,10,10,3,15);
figure(122); clf; scatplot(CY5_TMR_corr_int_only(:,1)/CY5_cloudInt(1),CY5_TMR_corr_int_only(:,2)/TMR_cloudInt(1),'circles',i,10,10,3,15);

counter = counter+1;
% ylim([0 60000])
% xlim([0 100000])
% xlabel('log(Xist Transcripts (from tsite intensity))','Fontsize',16)
% ylabel('log(Tsix Transcripts (from tsite intensity))','Fontsize',16) 
xlabel('Xist Transcripts (from tsite intensity)','Fontsize',16)
ylabel('Tsix Transcripts (from tsite intensity)','Fontsize',16) 
title('Integrated Intensities of Colocalized Spots','FontSize',16)
pause(.5)


end
%%% kstest of tsite intensities
save('Colocalization_intensities_second timecourse','CY5_TMR_corr_int_only')
CY5_coloc = CY5_TMR_corr_int_only(:,1)/CY5_cloudInt(1);
TMR_coloc = CY5_TMR_corr_int_only(:,2)/TMR_cloudInt(1);
[h,pval_cy5_thresh] = kstest2(TMR_coloc(CY5_coloc>=6),TMR_coloc(CY5_coloc<6))
figure(131); clf
subplot(1,2,1); histogram(TMR_coloc(CY5_coloc>=6),9); title('TMR Transc site above CY5 threshold')
subplot(1,2,2); histogram(TMR_coloc(CY5_coloc<6),9); title('TMR Transc site below CY5 threshold')
[h,pval_tmr_thresh] = kstest2(CY5_coloc(TMR_coloc>=6),CY5_coloc(TMR_coloc<6))
figure(132); clf
subplot(1,2,1); histogram(CY5_coloc(TMR_coloc>=6),9); title('CY5 Transc site above TMR threshold')
subplot(1,2,2); histogram(CY5_coloc(TMR_coloc<6),9); title('CY5 Transc site below TMR threshold')
%%% Spearman correlation
[rho,pval] = corr(CY5_coloc,TMR_coloc,'Type','Spearman')
[rho,pval] = corr(CY5_coloc,TMR_coloc,'Type','Kendall')
%%%
cov

figure(123); clf
CY5_TMR_corr_int_only_all = zeros(2,3);          %This will have the intensity of the corresponding TMR spot for every colocalized CY5 spot (same index placement as rnas in other matrices), but for all timepoints
counter2 = 1;
for k = 1:4
    C_T_corr_int = zeros(2,3);          %This will have the intensity of the corresponding TMR spot for every colocalized CY5 spot (same index placement as rnas in other matrices)
CY5_TMR_corr_int_only = zeros(2,1);         %This will only have the intensities of CY5 spots in first column and corresponding intensities in TMR in the second column
counter = 1;            %counter for CY5_TMR_corr_int_only
    for i = 1:size(CY5_AF594_TMR_tp{k},1)
%     if sum(C_T_col_tp{k}(i,:)==2)>0       %search if it has a colocalized spot
        init_ind = find(C_T_col_tp{k}(i,:)==2);   %Find indices of CY5 spots colocalized   
        noncol_ind = find(C_T_col_tp{k}(i,:)==1);   %Find indices of CY5 spots colocalized 
        noncol_ind1 = find(T_C_col_tp{k}(i,:)==1);   %Find indices of CY5 spots colocalized 
        for j = 1:size(init_ind,2)
%             if C_T_col_ind_tp{k}(i,init_ind(j)) > 2       %Avoid clouds
%                 CY5_TMR_corr_int_only(counter,1) =  col_br_CY5_tp{k}(i,init_ind(j));  %intensity of CY5           
%                 CY5_TMR_corr_int_only(counter,2) =  col_br_TMR_tp{k}(i,C_T_col_ind_tp{k}(i,init_ind(j)));
                CY5_TMR_corr_int_only(counter,1) =  CY5_int_tp{k}(i,init_ind(j));  %intensity of CY5                   
                CY5_TMR_corr_int_only(counter,2) =  TMR_int_tp{k}(i,C_T_col_ind_tp{k}(i,init_ind(j)));
                CY5_TMR_corr_int_only_all(counter2,1) =  CY5_int_tp{k}(i,init_ind(j));  %intensity of CY5                   
                CY5_TMR_corr_int_only_all(counter2,2) =  TMR_int_tp{k}(i,C_T_col_ind_tp{k}(i,init_ind(j)));
                CY5_TMR_corr_int_only_all(counter2,3) = k;
                counter2 = counter2+1
                counter = counter+1;
                %             end
        end
%         for j = 1:size(noncol_ind,2)                    %find noncolocalized spots
%             %             if C_T_col_ind(i,init_ind(j)) > 2       %Avoid clouds
%             CY5_TMR_corr_int_only(counter,1) =  CY5_int_tp{k}(i,noncol_ind(j));  %intensity of CY5
%             CY5_TMR_corr_int_only(counter,2) =  min(TMR_int(TMR_int>0));    %Set equal to minimum TMR intensity
%             counter = counter+1;
%             %             end
%         end
%         for j = 1:size(noncol_ind1,2)                    %find noncolocalized spots
%             %             if C_T_col_ind(i,init_ind(j)) > 2       %Avoid clouds
%             CY5_TMR_corr_int_only(counter,1) =  TMR_int_tp{k}(i,noncol_ind1(j));  %intensity of CY5
%             CY5_TMR_corr_int_only(counter,2) =  min(CY5_int(CY5_int>0));    %Set equal to minimum TMR intensity
%             counter = counter+1;
%             %             end
%         end
%         %%%Do this to include transcription sites not colocalized        
%         tssite_ind = find(TMR_int(i,:) >= 3*TMR_cloudInt(1)); %threshold for number of transcripts (second is greater than)
%         for j = 1:size(tssite_ind,2)
%             CY5_TMR_corr_int_only(counter,1) =  1;
%             CY5_TMR_corr_int_only(counter,2) =  TMR_int(i,tssite_ind(j));
%             counter = counter+1;
%         end
%     end
    end
    subplot(1,4,k)
 scatplot(CY5_TMR_corr_int_only(:,1)/CY5_cloudInt(1),CY5_TMR_corr_int_only(:,2)/TMR_cloudInt(1),'circles',4000/9000,10,10,3,15);
% scatter(CY5_TMR_corr_int_only(:,1)/CY5_cloudInt(1),CY5_TMR_corr_int_only(:,2)/TMR_cloudInt(1),25);
  
ylim([0 35])
 xlim([0 30])
if k >4 
    xlabel('Xist intensity (in transcripts)','Fontsize',16)
end
if k ==1 | k ==5
    ylabel('Tsix intensity (in transcripts)','Fontsize',16)   
end
title(times0{k}(1),'FontSize',12)
pause(.5)    
CY5_TMR_corr_int_only(CY5_TMR_corr_int_only <1) = 1;
log_temp = log2(CY5_TMR_corr_int_only);

CY5_TMR_corr_int_only(:,2)
end
figure(124); clf
gscatter(CY5_TMR_corr_int_only_all(:,1)/CY5_cloudInt(1),CY5_TMR_corr_int_only_all(:,2)/TMR_cloudInt(1),CY5_TMR_corr_int_only_all(:,3),...
    'MarkerSize',3);
  

 figure(121); clf; scatplot(log_temp(:,1),log_temp(:,2),'circles',.7,100,5,3,15);   %changing the value after 'circles' is important for dynamic range
ylim([10 20])
xlim([10 20])
xlabel('log2(Xist intensity)','Fontsize',16)
ylabel('log2(Tsix intensity)','Fontsize',16) 
title('Integrated Intensities of Colocalized Spots','FontSize',16)
counter1 = 1;
counter2 = 1;
hists_Xlow = [0;0];
hists_Xhi = [0;0];
for i = 1:size(log_temp,1)
    if log_temp(i,2) > 15
        log_temp(i,3) = 1;
        hists_Xhi(counter2,1) = log_temp(i,1);
        counter2=counter2+1
    else
        log_temp(i,3) = 2;
        hists_Xlow(counter1,1) = log_temp(i,1);
        counter1 = counter1+1
    end
end
figure(122); clf
h = boxplot(log_temp(:,1),log_temp(:,3),'Labels',{'Low colocalized Tsix intensity','High colocalized Tsix intensity'})
% xtickangle(90)
ylabel('log2 Xist intensity')
set(gca,'FontSize',14)
set(h,'LineWidth',2)
figure(125); clf;
h1 = histogram(hists_Xlow);
h1.Normalization = 'probability';
h1.BinWidth = 1;
hold on
% figure(124); clf;
h2 = histogram(hists_Xhi);
h2.Normalization = 'probability';
h2.BinWidth = 1;
legend('Low colocalized Tsix intensity','High colocalized Tsix intensity')
title('Brightness of Xist spots','FontSize',13)
xlabel('Fluorescence Intensity','FontSize',13)
ylabel('Probability','FontSize',13)

%Show histograms of different intensities
TMR_int(TMR_int <= 0) = NaN;
TMR_int_temp = TMR_int;
TMR_int_temp(isnan(TMR_int_temp)) = [];
TMR_int_transc = immultiply(TMR_int >= 3*TMR_cloudInt(1), TMR_int);
TMR_int_transc(TMR_int_transc <= 0) = [];
figure(123); clf;
h1 = histogram(CY5_TMR_corr_int_only(:,2));
h1.Normalization = 'probability';
h1.BinWidth = 5000;
hold on
% figure(124); clf;
h2 = histogram(TMR_int_transc(:));
h2.Normalization = 'probability';
h2.BinWidth = 5000;
% figure(125); clf;
h3 =  histogram(TMR_int_temp(:));
h3.Normalization = 'probability';
sum(h1.Values(:))
h3.BinWidth = 5000;
legend('Tsix colocalized with Xist','Tsix Transcription Sites (by Intensity)','All Tsix Spots')
xlabel('Fluorescence Intensity','FontSize',13)
ylabel('Probability','FontSize',13)
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

%% Looking for transcription sites of TMR (RepA) that are not colocalized with an Xist spot and vice-versa
TMR_name1 = 'Tsix'
times2 = [0,.25,.5,1,2,5,9];
colors = {'k' 'b' 'm' 'r' 'g' 'c' [255 159 0]/(255+159)}
cell_nums = zeros(1,size(TMR_int_tp,2)); %Stores the number of cells per timepoint
tsite_nums_TMR = zeros(1,size(TMR_int_tp,2)); %Stores the number of transcription sites per timepoint
tsite_cells_TMR = zeros(50,size(TMR_int_tp,2)); %Stores number of cells with number of of transcription sites. Row is how many transcription sites (+1), column is timepoint
tsite_col_TMR = zeros(1,size(TMR_int_tp,2)); %Stores the number of transcription sites that colocalizes with an Xist spot per timepoint
tsite_noncol_TMR = zeros(1,size(TMR_int_tp,2)); %Stores the number of transcription sites that does not colocalize with an Xist spot per timepoint
tsite_thres_TMR = 3;
tsite_nums_CY5 = zeros(1,size(TMR_int_tp,2)); %Stores the number of transcription sites per timepoint
tsite_cells_CY5 = zeros(50,size(TMR_int_tp,2)); %Stores number of cells with number of of transcription sites. Row is how many transcription sites (+1), column is timepoint
tsite_col_CY5 = zeros(1,size(TMR_int_tp,2)); %Stores the number of transcription sites that colocalizes with an Xist spot per timepoint
tsite_noncol_CY5 = zeros(1,size(TMR_int_tp,2)); %Stores the number of transcription sites that does not colocalize with an Xist spot per timepoint
tsite_thres_CY5 = 3;
for i = 1:size(TMR_int_tp,2)    %searches through every timepoint
    cell_nums(i) = size(TMR_int_tp{i},1);   %Stores how many cells per timepoint
    tsite_counter_tp_TMR = 0;   %will count how many transcription sites per timepoint
    tsite_col_temp_TMR = 0;     %will count how many transcription sites colocalize with CY5 per timepoint
    tsite_noncol_temp_TMR = 0;  %will count how many transcription sites DO NOT colocalize with CY5 per timepoint
    tsite_counter_tp_CY5 = 0;   %will count how many transcription sites per timepoint
    tsite_col_temp_CY5 = 0;     %will count how many transcription sites colocalize with CY5 per timepoint
    tsite_noncol_temp_CY5 = 0;  %will count how many transcription sites DO NOT colocalize with CY5 per timepoint
    for j = 1:size(TMR_int_tp{i},1) %searches through every cell
        
        tsite_temp_TMR = sum(TMR_int_tp{i}(j,:) >= tsite_thres_TMR*TMR_cloudInt(i));
        tsite_cells_TMR(tsite_temp_TMR+1,i) = tsite_cells_TMR(tsite_temp_TMR+1,i)+1;
        tsite_counter_tp_TMR = tsite_counter_tp_TMR+tsite_temp_TMR;
        
        tsite_temp_CY5 = sum(CY5_int_tp{i}(j,:) >= tsite_thres_CY5*CY5_cloudInt(i));
        tsite_cells_CY5(tsite_temp_CY5+1,i) = tsite_cells_CY5(tsite_temp_CY5+1,i)+1;
        tsite_counter_tp_CY5 = tsite_counter_tp_CY5+tsite_temp_CY5;
        %TMR
        for k = 1:size(TMR_int_tp{i},2) %searches through every transcript (TMR)
            if TMR_int_tp{i}(j,k) >= tsite_thres_TMR*TMR_cloudInt(i)    %if it is a transcription site
                %tsite_counter_tp = tsite_counter_tp +1; %Add to transcription site counter
                if T_C_col_tp{i}(j,k) == 2              %If there is a colocalized Xist
                    tsite_col_temp_TMR = tsite_col_temp_TMR+1;  %Add to counter of colocalized transcription sites
                elseif T_C_col_tp{i}(j,k) == 1              %If there is not a colocalized Xist
                    tsite_noncol_temp_TMR = tsite_noncol_temp_TMR+1; %Add to counter of noncolocalized transcription sites
                else
                    'There was an error. Transcription site determined where colocalization matrix said there was no transcript'
                end
            end
        end
        %CY5
       for k = 1:size(CY5_int_tp{i},2) %searches through every transcript (CY5)
            if CY5_int_tp{i}(j,k) >= tsite_thres_CY5*CY5_cloudInt(i)    %if it is a transcription site
                %tsite_counter_tp = tsite_counter_tp +1; %Add to transcription site counter
                if C_T_col_tp{i}(j,k) == 2              %If there is a colocalized Xist
                    tsite_col_temp_CY5 = tsite_col_temp_CY5+1;  %Add to counter of colocalized transcription sites
                elseif C_T_col_tp{i}(j,k) == 1              %If there is not a colocalized Xist
                    tsite_noncol_temp_CY5 = tsite_noncol_temp_CY5+1; %Add to counter of noncolocalized transcription sites
                else
                    'There was an error. Transcription site determined where colocalization matrix said there was no transcript'
                end
            end
        end
        
    end
    tsite_col_TMR(i) = tsite_col_temp_TMR;
    tsite_noncol_TMR(i) = tsite_noncol_temp_TMR;
    tsite_nums_TMR(i) = tsite_counter_tp_TMR;
    
    tsite_col_CY5(i) = tsite_col_temp_CY5;
    tsite_noncol_CY5(i) = tsite_noncol_temp_CY5;
    tsite_nums_CY5(i) = tsite_counter_tp_CY5;
end
norm1 = cell_nums.^-1; %Multiplying by this will normalize by number of cells
figure(800); clf; hold on
plot(times2,immultiply(norm1,tsite_nums_TMR),colors{4},'LineWidth',3)
plot(times2,immultiply(norm1,tsite_col_TMR),colors{3},'LineWidth',3)
plot(times2,immultiply(norm1,tsite_noncol_TMR),'Color',colors{7},'LineWidth',3)
plot(times2,immultiply(norm1,tsite_nums_CY5),colors{2},'LineWidth',3)
plot(times2,immultiply(norm1,tsite_col_CY5),colors{5},'LineWidth',3)
plot(times2,immultiply(norm1,tsite_noncol_CY5),colors{6},'LineWidth',3)
legend('Tsix Transcription Sites per Cell','Tsix Transcription Sites colocalized with Xist per Cell','Tsix Transcription Sites not Colocalized with Xist per Cell',...
    'Xist Transcription Sites per Cell','Xist Transcription Sites colocalized with Tsix per Cell','Xist Transcription Sites not Colocalized with Tsix per Cell')
xlabel('Time (Days)')
ylabel('Number per Cell')
% xlim([0 1])
% ylim([0 .5])
set(gca,'FontSize',13)

%%% Visualize transcription sites for TMR
TMR_pos_init_tp_ts = immultiply(TMR_pos_init_tp{i},repmat(TMR_int_tp{i}>tsite_thres_TMR*TMR_cloudInt(i),[1,1,3]));
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

%%% Visualize transcription sites for CY5
CY5_pos_init_tp_ts = immultiply(CY5_pos_init_tp{i},repmat(CY5_int_tp{i}>tsite_thres_CY5*CY5_cloudInt(i),[1,1,3]));
figure(900); clf;
scatter3(CY5_pos_init_tp{i}(j,:,2),CY5_pos_init_tp{i}(j,:,1),CY5_pos_init_tp{i}(j,:,3));
title('Localization of Xist spots in a single cell')
figure(901); clf; 
scatter3(CY5_pos_init_tp_ts(j,:,2),CY5_pos_init_tp_ts(j,:,1),CY5_pos_init_tp_ts(j,:,3));
title('Localization of Xist transcription sites in a single cell')
figure(902); clf;
scatter(CY5_pos_init_tp{i}(j,:,2),-1*CY5_pos_init_tp{i}(j,:,1));
title('Localization of Xist spots in a single cell')
xlim([0 400])
ylim([-400 0])
figure(903); clf; 
scatter(CY5_pos_init_tp_ts(j,:,2),-1*CY5_pos_init_tp_ts(j,:,1));
title('Localization of Xist transcription sites in a single cell')
xlim([0 400])
ylim([-400 0])
['timepoint ' num2str(CY5_AF594_TMR_tp{i}(j,4))]
['Exp date ' num2str(CY5_AF594_TMR_tp{i}(j,5))]
['image number ' num2str(CY5_AF594_TMR_tp{i}(j,6))]

%% Find number of trancription sites in cells with certain numbers of transcripts
tran_thres = 6; %threshold for number of transcripts (second is greater than)
tran_nums = zeros(2,50);    %Will store number fo transcription sites. The first row are cells below the threshold, second is above
tran_nums_tp = zeros(2,50,size(CY5_AF594_TMR_tp,2));    %Will store number fo transcription sites. The first row are cells below the threshold, second is above
tran_test = 3;      %Specifies which transcript is being thresholded for (1 is CY5, 2 is AF594, 3 is TMR)
tsite_test = 3;     %Specifies which transcription sites are being looked for (1 is CY5, 2 is AF594, 3 is TMR)
for i = 1:size(CY5_AF594_TMR_tp,2)  %Go through every timepoint
    for j = 1: size(CY5_AF594_TMR_tp{i},1) %Go through every cell
        if CY5_AF594_TMR_tp{i}(j,tran_test) > tran_thres   %Checks if there are more transcripts than threshold
            tran_nums_tp(2,CY5_AF594_TMR_tp{i}(j,7+tsite_test)+1,i) = tran_nums_tp(2,CY5_AF594_TMR_tp{i}(j,7+tsite_test)+1,i)+1;
            tran_nums(2,CY5_AF594_TMR_tp{i}(j,7+tsite_test)+1) = tran_nums(2,CY5_AF594_TMR_tp{i}(j,7+tsite_test)+1)+1;
        else
            tran_nums_tp(1,CY5_AF594_TMR_tp{i}(j,7+tsite_test)+1,i) = tran_nums_tp(1,CY5_AF594_TMR_tp{i}(j,7+tsite_test)+1,i)+1;
            tran_nums(1,CY5_AF594_TMR_tp{i}(j,7+tsite_test)+1) = tran_nums(1,CY5_AF594_TMR_tp{i}(j,7+tsite_test)+1)+1;
        end
    end
end
figure(1100); clf
bar(0:4,[tran_nums(1,1:5)',tran_nums(2,1:5)']); 
xlabel('Number of Transcription Sites')
ylabel('Number of Cells')
legend(['Below ' num2str(tran_thres) ' transcripts'],['Above ' num2str(tran_thres) ' transcripts'])

figure(1101); scatplot(CY5_AF594_TMR(:,10),CY5_AF594_TMR(:,3),'circles',15,100,5,3,15);
xlabel('Number of Transcription Sites')
ylabel('Number of Transcripts')
title('Tsix')

figure(1101); scatter(CY5_AF594_TMR(:,10),CY5_AF594_TMR(:,3));

figure(1102); clf; 
x = randn(2000,1);
y = 1 + randn(5000,1);
h1 = histogram(x);
hold on
h2 = histogram(y);

%% Find the distributions over time for the TMR transcripts in cells that are either Xist positive or negative
colors = {'k' 'b' 'm' 'r' 'g' 'c' [255 159 0]/(255+159)}
for mm = 1:15;%:2:12; %2:2:10 %Go through different thresholds
figure(20+mm); clf; hold on
figure(40+mm); clf; hold on
figure(60+mm); clf; hold on
figure(80+mm); clf; hold on
CY5_thres = mm; %The threshold that CY5 (Xist) needs to be equal or higher to be part of the larger group
bin_num = 10;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 20;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
 minV = 0;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
% temp1 = CY5_AF594_TMR(:,3);
% temp2 = temp1(temp1 >= minV);
x_col = 3;  %column of data from CY5_AF594_TMR that will be on the X axis (16 is log2 of column 1)
x_name = 'Tsix';
log_x = 1;  %Set to 1 to use log2 scale for x axis
thres_col = 1;  %column of data from CY5_AF594_TMR that will be thresholded
thres_name = 'Xist';
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1);%+bin_size/2;
%%%Make labels
clear X_labels
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
%%%
clear binsall
binsall = {};
for j = 1:7
    counter1 = 1;   %counter for low CY5 value matrix
    counter2 = 1; %counter for high CY5 value matrix
    temp1 = zeros(2,1);      %matrix for low CY5 value
    temp2 = zeros(2,1);      %matrix for high CY5 values
    temp3 = CY5_AF594_TMR_tp{j}(:,1);   %Used for log values
    temp3(temp3 < minV) = minV;   
    CY5_AF594_TMR_tp{j}(:,16) = log2(temp3);
    temp4 = CY5_AF594_TMR_tp{j};
    %%%Find cells below or above the threshold
    for k = 1:size(CY5_AF594_TMR_tp{j},1)
        if CY5_AF594_TMR_tp{j}(k,thres_col) < CY5_thres
            temp1(counter1) = CY5_AF594_TMR_tp{j}(k,x_col);
            temp4(k,17) = 1;    %1 is low Xist expression
            counter1 = counter1+1;
        else
            temp2(counter2) = CY5_AF594_TMR_tp{j}(k,x_col);
            counter2 = counter2+1;
            temp4(k,17) = 2;    %1 is high Xist expression
        end
    end    
    %%%
    [something1,p_val] = kstest2(temp1,temp2);
    figure(80+mm); hold on
    subplot(2,4,j)
    h = boxplot(temp4(:,3),temp4(:,17),'Labels',{'Xist-','Xist+'})
% xtickangle(90)
ylabel('Tsix Molecules')
title(times0{j})
set(gca,'FontSize',14)
set(h,'LineWidth',2)
sgtitle(['Threshold = ' num2str(CY5_thres)])
binsall{j} = zeros(2,size(binsX,2));      %Has number of cells in each bin. first row is lower CY5, second is higher
binall = zeros(2,size(binsX,2));      %Has number of cells in each bin. first row is lower CY5, second is higher
%binsall = imhistc(temp2,bin_num,0,maxV);
for i = 1:size(X_labels,2)
    if i == size(X_labels,2)
     binsall{j}(1,i) = size(find(temp1 >= binsX(i)),1)
     binsall{j}(1,i) = size(find(temp2 >= binsX(i)),1)
    else
     binsall{j}(1,i) = size(find(temp1 >= binsX(i)),1)-size(find(temp1 >= binsX(i+1)),1);
     binsall{j}(2,i) = size(find(temp2 >= binsX(i)),1)-size(find(temp2 >= binsX(i+1)),1);
    end
end
binall(1,:) = binsall{j}(1,:)/sum(binsall{j}(1,:));  %normalize
binall(2,:) = binsall{j}(2,:)/sum(binsall{j}(2,:));  %normalize
figure(20+mm)
subplot(1,2,1)
plot(1:size(X_labels,2),binall(1,:),'Color',colors{j},'LineWidth',3); hold on
xticks(1:size(X_labels,2))
xticklabels(X_labels)
xtickangle(45)
subplot(1,2,2)
plot(size(X_labels,2),binall(2,:),'Color',colors{j},'LineWidth',3); hold on
xticks(1:size(X_labels,2))
xticklabels(X_labels)
xtickangle(45)
figure(40+mm)
subplot(2,4,j)
plot(1:size(X_labels,2),binall(1,:),'Color',colors{j},'LineWidth',3); hold on
plot(1:size(X_labels,2),binall(2,:),'Color',colors{j},'LineWidth',3,'LineStyle','--'); hold on
xticks(1:size(X_labels,2))
xticklabels(X_labels)
xtickangle(45)
legend([thres_name '- (' num2str(counter1-1) ' cells)'],[thres_name '+ (' num2str(counter2-1) ' cells)'],'Location','northeast')
xlabel([x_name ' Molecules'])
%ylabel('Number of Cells')
ylabel('Probability')
title(times0{j})
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
sgtitle(['Threshold = ' num2str(CY5_thres)])
figure(60+mm)
subplot(2,4,j)
binall_cum = zeros(size(binall));
for ab = 1:size(binall,2)
    binall_cum(1,ab) = sum(binall(1,1:ab));
    binall_cum(2,ab) = sum(binall(2,1:ab));
end 
plot(1:size(X_labels,2),binall_cum(1,:),'Color',colors{j},'LineWidth',3); hold on
plot(1:size(X_labels,2),binall_cum(2,:),'Color',colors{j},'LineWidth',3,'LineStyle','--'); hold on
xticks(1:size(X_labels,2))
xticklabels(X_labels)
xtickangle(45)
ylim([0 1])
text(5,.5,['p = ' num2str(p_val)])
legend([thres_name '- (' num2str(counter1-1) ' cells)'],[thres_name '+ (' num2str(counter2-1) ' cells)'],'Location','southeast')
xlabel([x_name ' Molecules'])
%ylabel('Number of Cells')
ylabel('Cumulative Probability')
title(times0{j})
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
sgtitle(['Threshold = ' num2str(CY5_thres)])
end
%binall = binsall/sum(binsall);  %normalize
%bar(labelsX,binsall)
%plot(labelsX,binall,'g','LineWidth',3); hold on
figure(20+mm)
subplot(1,2,1)
legend('0d', '6hr', '12hr', '1d', '2d', '5d','9d')
title(['Distribution for ' thres_name ' Negative Cells'])
xlabel([x_name ' Molecules'])
%ylabel('Number of Cells')
ylabel('Probability')
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
subplot(1,2,2)
legend('0d', '6hr', '12hr', '1d', '2d', '5d','9d')
title(['Distribution for ' thres_name ' Positive Cells'])
xlabel(['log2 ' x_name ' Molecules'])
%ylabel('Number of Cells')
ylabel('Probability')
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
end
for i = 1:15
    figure(60+i)
    pause(5)
end
%% On fractions of cells with different amounts of one transcript (binned) over time
colors = {'k' 'b' 'm' 'r' 'g' 'c' [255 159 0]/(255+159)}
for mm = 1:10:101   %Go through different thresholds for being considered "on"
figure(40+mm); clf; hold on
bin_num = 10;
%maxV = max(CY5_AF594_TMR(:,1));
maxV = 30;         %cutoff for max value shown. Hist is made with vale of final bin = inside bin + anything above (no cells are lost)
minV = 0;           %cutoff for min. THIS REMOVES CELLS BELOW THE MIN FROM THE HISTOGRAM
thres_on = mm;
% temp1 = CY5_AF594_TMR(:,3);
% temp2 = temp1(temp1 >= minV);
x_col = 3;  %column of data from CY5_AF594_TMR that will be on the X axis (binned)
x_name = 'Tsix';
log_x = 0;  %Set to 1 to use log2 scale for x axis
y_col = 1;  %column of data from CY5_AF594_TMR that will be thresholded
y_name = 'Xist';
bin_size = (maxV-minV)/bin_num;
binsX = minV:bin_size:maxV;
labelsX = binsX(1:size(binsX,2)-1);%+bin_size/2;
clear binsall
binsall = {};
%%%Make labels
clear X_labels
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
%%%
for j = 1:7
    temp1 = CY5_AF594_TMR_tp{j}(:,x_col);      %matrix for binned x axis 
%     CY5_AF594_TMR_tp{j}(:,16) = CY5_AF594_TMR_tp{j}(:,y_col)>=thres_on;    
binsall{j} = zeros(2,size(binsX,2));      %Has number of cells in each bin. first row is the number of cells, second is on fraction
binall = zeros(2,size(binsX,2));      %Has number of cells in each bin. first row is the number of cells, second is on fraction
CY5_AF594_TMR_tp{j}(:,15) = 0;
%binsall = imhistc(temp2,bin_num,0,maxV);
for i = 1:size(X_labels,2)  % Go through each bin
    %%%Makes index for bin number in 15th column
        for k = 1:size(CY5_AF594_TMR_tp{j},1)   %Go through each cell
        if i == bin_num+1
            if CY5_AF594_TMR(k,x_col) >= binsX(i)
                CY5_AF594_TMR_tp{j}(k,15) = i;            
            end
       elseif CY5_AF594_TMR_tp{j}(k,x_col) >= binsX(i) & CY5_AF594_TMR_tp{j}(k,x_col) < binsX(i+1)
           CY5_AF594_TMR_tp{j}(k,15) = i;
        end
        end
    %%% Stores On fraction of cells (in y_col) for each bin
    temp2 = CY5_AF594_TMR_tp{j}(:,y_col);
    temp3 = temp2(find(CY5_AF594_TMR_tp{j}(:,15) == i));
    binsall{j}(2,i) = sum(temp3 > thres_on)/size(temp3,1);
    %%%Stores the number of cells in each bin
    if i == size(binsall{j},2)
     binsall{j}(1,i) = size(find(temp1 >= binsX(i)),1)
    else
     binsall{j}(1,i) = size(find(temp1 >= binsX(i)),1)-size(find(temp1 >= binsX(i+1)),1);
    end
    %%%
end
%  binall(1,:) = binsall{j}(1,:)/sum(binsall{j}(1,:));  %normalize
binall = binsall{j};
subplot(2,4,j)
plot(1:size(X_labels,2),binall(2,:),'Color',colors{j},'LineWidth',3); hold on
xticks(1:size(X_labels,2))
xticklabels(X_labels)
xtickangle(45)
ylabel([ y_name ' On Fraction'])
ylim([0 1])
yyaxis right
% plot(1:size(X_labels,2),binall(1,:),'LineWidth',3); hold on
scatter(1:size(X_labels,2),binall(1,:)); hold on
 xlabel([x_name ' Molecules'])
ylabel('Number of Cells')
title(times0{j})
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
end
end
%% Look at TMR distances to the cloud over time
%%%Note: It is impossible to store NaN into integer matrices, so be sure to
%%%remove a lot of the distances based on whether there is a NaN in the
%%%position matrices.
%%%Also, the position matrices have 0's that should be NaN. Be sure to
%%%change that
colors = {'k' 'b' 'm' 'r' 'g' 'c' [255 159 0]/(255+159)}
CY5_int1 = CY5_int;
CY5_int1(CY5_int <= 0) = NaN;    %Make all 0's NaN
TMR_int1 = TMR_int;
TMR_int1(TMR_int <= 0) = NaN;    %Make all 0's NaN
C1_TMR_dist1 = C1_TMR_dist;
C1_TMR_dist1(isnan(TMR_int1)) = NaN; %Make everything NaN that should be in the distances
% C1_TMR_dist1(C1_TMR_dist1<=2) = NaN; %Make everything NaN that should be in the distances
times2 = {'0d', '6hr', '12hr','1d','2d','5d','9d'}; 
figure(700); clf; histogram(C1_TMR_dist1);
bin_num = 10;
figure(701); clf    %all spots
figure(702); clf    %mean for each cell
figure(703); clf    %median for each cell
for i = 1:size(C1_TMR_dist_tp,2)
    TMR_int1 = TMR_int_tp{i};
    TMR_int1(TMR_int1 <= 0) = NaN;    %Make all 0's NaN
    C1_TMR_dist1 = C1_TMR_dist_tp{i};
    C1_TMR_dist1(isnan(TMR_int1)) = NaN; %Make everything NaN that should be in the distances
    C1_TMR_dist1(C1_TMR_dist1<=2) = NaN; %Make everything NaN that should be in the distances
    C1_TMR_dist1 = C1_TMR_dist1.^.5;
    [counts1, binCenters1] = hist(C1_TMR_dist1(:), bin_num);
    counts1 = counts1/sum(counts1);
    figure(701)
    plot(binCenters1, counts1, 'Color',colors{i});
    hold on;   
    grid on;
    %%%Plot the histogram of mean distances per cell
    mean_dists = nanmean(C1_TMR_dist1,2);  %Find the mean distance for each cell
    [counts1, binCenters1] = hist(mean_dists, bin_num);
    counts1 = counts1/sum(counts1);
    figure(702)
    plot(binCenters1, counts1, 'Color',colors{i});
    hold on
    median_dists = nanmedian(C1_TMR_dist1,2);  %Find the mean distance for each cell
    [counts1, binCenters1] = hist(median_dists, bin_num);
    counts1 = counts1/sum(counts1);
    figure(703)
    plot(binCenters1, counts1, 'Color',colors{i});
    hold on
end
    figure(701)
    title('Distances from TMR spots to the first Cloud')
    xlabel('Distance')
    ylabel('Probability')
    legend(times2); 
    figure(702)
    title('Mean Distances per Cell from TMR spots to the first Cloud')
    xlabel('Distance')
    ylabel('Probability')
    legend(times2);    
    figure(703)
    title('Median Distances per Cell from TMR spots to the first Cloud')
    xlabel('Distance')
    ylabel('Probability')
    legend(times2);    
%% Look at the distance of TMR transcription sites away from xist transcription sites (also based on t site intensity)
colors = {'k' 'b' 'm' 'r' 'g' 'c' [255 159 0]/(255+159)}
CY5_int1 = immultiply(CY5_int,CY5_int>CY5_cloudInt(1));
CY5_int1(CY5_int <= 0) = NaN;    %Make all 0's NaN
TMR_int1 = immultiply(TMR_int,TMR_int>TMR_cloudInt(1)); %Find TMR transcription sites
TMR_int1(TMR_int <= 0) = NaN;    %Make all 0's NaN
C1_TMR_dist1 = C1_TMR_dist;
C1_TMR_dist1(isnan(TMR_int1)) = NaN; %Make everything NaN that should be in the distances
% C1_TMR_dist1(C1_TMR_dist1<=2) = NaN; %Make everything NaN that should be in the distances
times2 = {'0d', '6hr', '12hr','1d','2d','5d','9d'}; 
figure(700); clf; histogram(C1_TMR_dist1);
bin_num = 10;
figure(701); clf    %all spots
figure(702); clf    %mean for each cell
figure(703); clf    %median for each cell
for i = 1:size(C1_TMR_dist_tp,2)
    TMR_int1 = TMR_int_tp{i};
    TMR_int1(TMR_int1 <= 0) = NaN;    %Make all 0's NaN
    C1_TMR_dist1 = C1_TMR_dist_tp{i};
    C1_TMR_dist1(isnan(TMR_int1)) = NaN; %Make everything NaN that should be in the distances
    C1_TMR_dist1(C1_TMR_dist1<=2) = NaN; %Make everything NaN that should be in the distances
    C1_TMR_dist1 = C1_TMR_dist1.^.5;
    [counts1, binCenters1] = hist(C1_TMR_dist1(:), bin_num);
    counts1 = counts1/sum(counts1);
    figure(701)
    plot(binCenters1, counts1, 'Color',colors{i});
    hold on;   
    grid on;
    %%%Plot the histogram of mean distances per cell
    mean_dists = nanmean(C1_TMR_dist1,2);  %Find the mean distance for each cell
    [counts1, binCenters1] = hist(mean_dists, bin_num);
    counts1 = counts1/sum(counts1);
    figure(702)
    plot(binCenters1, counts1, 'Color',colors{i});
    hold on
    median_dists = nanmedian(C1_TMR_dist1,2);  %Find the mean distance for each cell
    [counts1, binCenters1] = hist(median_dists, bin_num);
    counts1 = counts1/sum(counts1);
    figure(703)
    plot(binCenters1, counts1, 'Color',colors{i});
    hold on
end
    figure(701)
    title('Distances from TMR spots to the first Cloud')
    xlabel('Distance')
    ylabel('Probability')
    legend(times2); 
    figure(702)
    title('Mean Distances per Cell from TMR spots to the first Cloud')
    xlabel('Distance')
    ylabel('Probability')
    legend(times2);    
    figure(703)
    title('Median Distances per Cell from TMR spots to the first Cloud')
    xlabel('Distance')
    ylabel('Probability')
    legend(times2);    
%% Look for number of transcription sites of Tsix with colocalized Xist
col_tp = zeros(1,7);
noncol_tp = zeros(1,7);
for j = 1:7;
    number_col = 0;
    number_noncol = 0;
    for m = 1:size(TMR_int_tp{j},1)
        TSite_ind = find(TMR_int_tp{j}(m,:) >= 3*TMR_cloudInt(1));
        if not(isempty(TSite_ind))
            number_col = number_col+sum(T_C_col_tp{j}(m,TSite_ind)==2);
            number_noncol = number_noncol+sum(T_C_col_tp{j}(m,TSite_ind)==1);
        end
    end
    col_tp(j) = number_col;
    noncol_tp(j) = number_noncol;
end
 figure(99); clf; plot(times1(1:7),col_tp./(col_tp+noncol_tp),'b','LineWidth',3); hold on
  plot(times1(1:7),noncol_tp./(col_tp+noncol_tp),'r','LineWidth',3)
  xlabel('Time (Days)')
  ylabel('Fraction of Transcription Sites')
  legend('Colocalized with Xist','Not Colocalized with Xist')
  title(['Tsix Transcription Sites within ' num2str(dist_col) ' pixels of Xist'])
  set(gca,'Fontsize',14)
 
  %% Colocalization Visualization         
 tp1 = 1;
 dist_bins = 1;
 im_num = 1;
  imnum_cell_ind = find(CY5_AF594_TMR_tp{tp1}(:,6) == im_num); %Find each cell i image
%%% Load files
 filename = ['mRNA_' dates0{tp1} '_' strain TMR_name{tp1} num2str(thTMRall{tp1}(thNum_CY5)) '_' num2str(thAF594all{tp1}(thNum_CY5)) CY5_name{tp1} num2str(thCY5all{tp1}(thNum_CY5)) '_' times0{tp1}{1} '_im' num2str(im_num) '_non_stringent_cloud'];
%   filename = ['mRNA_' dates0{tp1} '_F1-2-1' TMR_name{tp1} num2str(thTMRall{tp1}(thNum_CY5)) '_' num2str(thAF594all{tp1}(thNum_CY5)) CY5_name{tp1} num2str(thCY5all{tp1}(thNum_CY5)) '_' times0{tp1}{1} '_im' num2str(im_num)];
 load(filename,['PARcy5_' nucTH],['PARtmr_' nucTH]); 
im_name = [dates0{tp1} '_' times0{tp1}{1} '_' strain CY5_name{tp1}(1:size(CY5_name{tp1},2)-3) '-' TMR_name{tp1}(2:size(CY5_name{tp1},2)-3) '_img_' num2str(im_num)];
% im_name_diff = [im_name(1:size(im_name,2)-2) im_name(size(im_name,2))]
% filename = ['D:\' dates0{tp1} '\' im_name '\' im_name '_MMStack.ome.tif'];
filename = ['C:\Users\keslerbk\Desktop\Microscopy\' dates0{tp1} '\' im_name '\' im_name '_MMStack.ome.tif'];
[stack, img_read] = tiffread2(filename);
Ych = [0 1 1 0 1 0 0 0 1];
%     chi = find(handles.Ych == 1);                                                       % Determines how many channels have been imaged
%     ch = size(chi,2);
%in Ych, first is the total number of channels, second is the marker slice, and third is the boundary slice
chi = find(Ych == 1);
ch = size(chi,2);
ImSt = img_read/ch

        if Ych(3) == 1;                                                     % CY5
            ij = find(chi==3);           
            CY5im = [ij:ch:img_read];
            CY5_ims = zeros([size(stack(1,CY5im(1)).data) size(CY5im,2)]);
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
            TMR_ims = zeros([size(stack(1,TMRim(1)).data) size(TMRim,2)]);
            for i = 1:ImSt;
                TMR_ims(:,:,i) = stack(1,TMRim(i)).data;
            end;
        else
        end
     fileLAB = char(['Lab_' dates0{tp1} '_' times0{tp1}{1} '_F1-2-1'  CY5_name{tp1}(1:size(CY5_name{tp1},2)-3) '-' TMR_name{tp1}(2:size(CY5_name{tp1},2)-3) '_img_' num2str(im_num) '.mat']);
 load(fileLAB,'cells'); %'NucInfo','CellInfo',
 Lab = cells;
%  imnum_cell_ind = 1:max(cells(:));  %Determines the indexes of cells (only works if looking at images1)
%  Mid_DAPI = max(DAPI_ims(:,:,round(ImSt/2)),[],3);
%%%Visualize Clouds
Clouds0 = max(PARcy5_mid.clouds1+PARcy5_mid.clouds2,[],3);  %Get both CY5 clouds
Clouds1 = bwmorph(Clouds0,'remove')*100000;
Clouds1 = bwmorph(Clouds1,'thicken',2);
Clouds1 = double(Clouds1);
% Clouds1 = zeros(size(cells));
if TMR_hascloud
    TClouds0 = max(PARtmr_mid.clouds1+PARtmr_mid.clouds2,[],3); %Get both TMR clouds
    TClouds1 = bwmorph(TClouds0,'remove')*100000;
%     TClouds1 = bwmorph(TClouds1,'thicken',2);
else
    TClouds1 = zeros(size(cells));
end
TClouds1 = double(TClouds1);
CY50 = max(CY5_ims(:,:,5:end),[],3);
 CY51 = CY50/max(CY50(:));
 low_CY = max([0.0 median(CY51(:))-.01])
 hi_CY = min([1.0 median(CY51(:))+.1])
 CY52 = imadjust(CY51,[low_CY hi_CY]);
TMR0 = max(TMR_ims(:,:,5:end),[],3);
  TMR1 = TMR0/max(TMR0(:));
   low_TMR = max([0 median(TMR1(:))+.01])
 hi_TMR = min([1 median(TMR1(:))+.1])
 TMR2 = imadjust(TMR1,[low_TMR hi_TMR]);
 %          Nuc1 = Mid_DAPI/max(DAPI_ims(:));
 %             Nuc2 = imadjust(Nuc1,[0 1]);
CellBorder1 = bwmorph(cells,'remove')*100000;
CellBorder1 = bwmorph(CellBorder1,'thicken',2);
CellBorder1 = double(CellBorder1);
%  Trans1 = max(TRANS_ims(:,:,round(ImSt/2)),[],3);
%  Trans2 = Trans1/max(TRANS_ims(:));
%  Trans3 = imadjust(Trans2,[0 1]);
 R = TMR2+CellBorder1+TClouds1;
 B = CellBorder1+Clouds1+TClouds1;
 G = CY52+CellBorder1;
 RGB = cat(3,R,G,B);
 figure(3000); clf
%   imshow(CY50,[]);
 imshow(RGB,[]); hold on
  A = size(CY52,1);
  for m = imnum_cell_ind' %Go through each cell in image
      k1 = Lab == m-imnum_cell_ind(1)+1 ; %== j;%sieve out dots in cell j.
      %figure; imshow(k1,[0 1]);
      k1=uint16(k1);
      k3=regionprops(k1,'BoundingBox','Area');
      k4 = k3.BoundingBox; %create the rectangular box around the cell.
      k5 = regionprops(k1,'Centroid');
      if CY5_AF594_TMR_tp{tp1}(m,44) == 0 
          color1 = 'magenta';
      elseif CY5_AF594_TMR_tp{tp1}(m,44) == 1 
          color1 = 'cyan';
      else
          color1 = 'magenta';
      end    
      text(k5.Centroid(1),k5.Centroid(2),num2str(CY5_AF594_TMR_tp{tp1}(m,44)),'Color',color1)
      X0=round(k4(1))-4;
      Y0=round(k4(2))-4;
      X1=round(k4(1)+ k4(3))+4;
      Y1=round(k4(2)+ k4(4))+4;
      if X0 < 1; X0 = 1; end;
      if Y0 < 1; Y0 = 1; end;
      if X1 > A; X1 = A; end;
      if Y1 > A; Y1 = A; end;
      %%%Using TSite location in CY5_AF594_TMR_tp
%       for i = 1:size(CY5_AF594_TMR_tp{tp1},2)    %Go through each transcript
          shape1 = 'x'; %Make square if transcription site
          for j = 1:4
              if CY5_AF594_TMR_tp{tp1}(m,20+(3*(j-1))) > 0
                  plot(CY5_AF594_TMR_tp{tp1}(m,21+(3*(j-1)))+X0,CY5_AF594_TMR_tp{tp1}(m,20+(3*(j-1)))+Y0,['r' shape1]);
                  %                   'There is a determined transcription site'
                  if CY5_AF594_TMR_tp{tp1}(m,21+(3*(j-1))) == CY5_AF594_TMR_tp{tp1}(m,20+(3*(j-1)))
                      ['They are the same and equal to' num2str(CY5_AF594_TMR_tp{tp1}(m,21+(3*(j-1))))]
                  end
              end
              if CY5_AF594_TMR_tp{tp1}(m,33+(3*(j-1))) > 0
                  plot(CY5_AF594_TMR_tp{tp1}(m,33+(3*(j-1)))+X0,CY5_AF594_TMR_tp{tp1}(m,32+(3*(j-1)))+Y0,['g' shape1]);
                  %                   'There is a determined transcription site'
                  if CY5_AF594_TMR_tp{tp1}(m,33+(3*(j-1))) == CY5_AF594_TMR_tp{tp1}(m,32+(3*(j-1)))
                      ['They are the same and equal to' num2str(CY5_AF594_TMR_tp{tp1}(m,33+(3*(j-1))))];
                  end
              end
          end
%       end
      %     k2 = k1(Y0:Y1,X0:X1);
      %%%Using colocalization information
      for i = 1:size(T_C_col_tp{tp1},2)    %Go through each transcript
          if TMR_int_tp{tp1}(m,i) > TMR_tran_thres*TMR_cloudInt(1)
              shape1 = 's'; %Make square if transcription site
          else
              shape1 = 'o' ; %Make circle if not transcription site
          end
          if T_C_col_tp{tp1}(m,i) == 1
              plot(TMR_pos_init_tp{tp1}(m,i,2)+X0,TMR_pos_init_tp{tp1}(m,i,1)+Y0,['g' shape1])
          elseif T_C_col_tp{tp1}(m,i) == 2
              plot(TMR_pos_init_tp{tp1}(m,i,2)+X0,TMR_pos_init_tp{tp1}(m,i,1)+Y0,['y' shape1])
          end
      end
      for i = 1:size(C_T_col_tp{tp1},2)    %Go through each transcript
          if CY5_int_tp{tp1}(m,i) > CY5_tran_thres*CY5_cloudInt(1)
              shape1 = 's';
          else
              shape1 = 'o';
          end
          if C_T_col_tp{tp1}(m,i) == 1
              plot(CY5_pos_init_tp{tp1}(m,i,2)+X0,CY5_pos_init_tp{tp1}(m,i,1)+Y0,['r' shape1])
          elseif C_T_col_tp{tp1}(m,i) == 2
              plot(CY5_pos_init_tp{tp1}(m,i,2)+X0,CY5_pos_init_tp{tp1}(m,i,1)+Y0,['y' shape1])
          end
      end
  end
 savefig(['Colocalization plot_' date CY5_name{1} TMR_name{1}  'tp' times0{tp1}{1} '_im' num2str(im_num) '_CY5th' num2str(thNum_CY5) '_TMRth' num2str(thNum_TMR)]);

%% Quantify X Chromosome pairing
figure(400); clf
figure(401); clf
figure(402); clf
figure(403); clf
figure(404); clf
figure(405); clf
clear pairwise_dists min_pairwise_dists xticklabel max_pairwise_dists
minDist = 14;   %distance has to be above this to be considered a real distance (and not colocalized signal)
edges1 = 0:2000/65:10*2000/65;
edges2 = edges1; edges2(1) = minDist;
edges3 = 0:1000/65:10*2000/65;
xticklabel1 = {};
for i = 1:size(edges1,2)-1
    if false;%i == 1 | i == size(edges1,2)-1
        xticklabel1{i} = ['<' num2str(edges1(i+1)*65/1000)];
    else
        xticklabel1{i} = [num2str(edges1(i+1)*65/1000)];
    end
end
pairwise_dists = {};
pairwise_dists_maxint = {};
min_pairwise_dists = {};
max_pairwise_dists = {};
max_pairwise_dists_cy5tmr = {};
min_pairwise_dists_maxint = {};
for j = 1:7 %Timepoint
    pairwise_dists{j} = zeros(size(CY5_AF594_TMR_tp{j},1),1);
    min_pairwise_dists{j} = zeros(size(CY5_AF594_TMR_tp{j},1),1);
    max_pairwise_dists{j} = zeros(size(CY5_AF594_TMR_tp{j},1),1);
    max_pairwise_dists_cy5tmr{j} = zeros(size(CY5_AF594_TMR_tp{j},1),1);
    for m = 1:size(CY5_AF594_TMR_tp{j},1)
        pos_cy51 = zeros(1,1,3);
        pos_tmr1 = zeros(1,1,3);
        pos_both1 = zeros(1,1,3);
        countercy5 = 1;
        countertmr = 1;
        counterboth = 1;
        for i = 20:3:29
            if CY5_AF594_TMR_tp{j}(m,i) > 0
                pos_cy51(1,countercy5,1:3) = CY5_AF594_TMR_tp{j}(m,i:i+2);
                pos_both1(1,counterboth,1:3) = CY5_AF594_TMR_tp{j}(m,i:i+2);
                countercy5 = countercy5+1;
                counterboth = counterboth+1;
            else
                pos_cy51(1,countercy5,1:3) = NaN;
                countercy5 = countercy5+1;
            end
        end
        for i = 32:3:41
            if CY5_AF594_TMR_tp{j}(m,i) > 0
                pos_tmr1(1,countertmr,1:3) = CY5_AF594_TMR_tp{j}(m,i:i+2);
                pos_both1(1,counterboth,1:3) = CY5_AF594_TMR_tp{j}(m,i:i+2);
                counterboth = counterboth+1;
                countertmr = countertmr+1;
            else
                pos_tmr1(1,countertmr,1:3) = NaN;
                countertmr = countertmr+1;                
            end
        end
        %%%Distances between TMR and CY5 transcription sites
        pos_tmr2 = repmat(permute(pos_tmr1,[2,1,3]),[1, size(pos_cy51,2),1]);
        pos_cy52 = repmat(pos_cy51,[size(pos_tmr1,2),1,1]);
        diffs_sq = (pos_tmr2-pos_cy52).^2;
        dists = sum(diffs_sq,3);
        dists = dists.^.5;
        dists = dists(:);
        pairwise_dists{j}(m,1:size(dists,1)) = dists;
        max_pairwise_dists_cy5tmr{j}(m,1) = max(dists(:));
        %%%
        pos_both2 = repmat(permute(pos_both1,[2,1,3]),[1, size(pos_both1,2),1]);
        pos_both3 = repmat(pos_both1,[size(pos_both1,2),1,1]);
        pos_both2_maxint = pos_both2(:,:,1:2)
        pos_both3_maxint = pos_both3(:,:,1:2)
        diffs_sq = (pos_both2-pos_both3).^2;
        dists = sum(diffs_sq,3);
        dists = dists.^.5;
        dists = dists(:);
        pairwise_dists{j}(m,1:size(dists,1)) = dists;
        if not(isempty(dists(dists>minDist)))
            min_pairwise_dists{j}(m,1) = min(dists(dists>minDist));
        end
        if not(isempty(dists(dists>minDist)))
            max_pairwise_dists{j}(m,1) = max(dists(dists>minDist));
        end
        for i = 1:size(edges2,2)-1
                if max_pairwise_dists{j}(m,1) < edges2(i+1) & max_pairwise_dists{j}(m,1) >= edges2(i)
                CY5_AF594_TMR_tp{j}(m,44) = i;
                end
        end
        %%% for maximum intensity projection
        diffs_sq = (pos_both2_maxint-pos_both3_maxint).^2;
        dists = sum(diffs_sq,3);
        dists = dists.^.5;
        dists = dists(:);
        pairwise_dists_maxint{j}(m,1:size(dists,1)) = dists;
        if not(isempty(dists(dists>minDist)))
            min_pairwise_dists_maxint{j}(m,1) = min(dists(dists>minDist));
        end
    end
    %%% For minimum distance
    temp_min_pairwise_dists = min_pairwise_dists{j}(min_pairwise_dists{j}>0);
    max(temp_min_pairwise_dists)
    cellnum = size(temp_min_pairwise_dists,1);
    figure(400);
    subplot(2,4,j)
    histogram(temp_min_pairwise_dists,edges1,'Normalization','probability')
    xticks(edges1(1:1:end))
     xticklabels(xticklabel1(1:1:end))
     if j >4
     xlabel('Distance (\mum)')
     end
     xtickangle(45)
     if j == 1 | j == 5
         ylabel('Probability')
     end
     ylim([0 .5])
    title([times0{j}{1} ' n=' num2str(cellnum) '/' num2str(size(CY5_AF594_TMR_tp{j},1))])
    set(gca,'Fontsize',14)
    %%% for maximum intensity projection of min distance
    temp_min_pairwise_dists = min_pairwise_dists_maxint{j}(min_pairwise_dists_maxint{j}>0);
    max(temp_min_pairwise_dists)
    cellnum = size(temp_min_pairwise_dists,1);
    figure(401);
    subplot(2,4,j)
    histogram(temp_min_pairwise_dists,edges1,'Normalization','probability')
    xticks(edges1(1:1:end))
     xticklabels(xticklabel1(1:1:end))
     if j >4
     xlabel('Distance (\mum)')
     end
     xtickangle(45)
     if j == 1 | j == 5
         ylabel('Probability')
     end
     ylim([0 .5])
    title([times0{j}{1} ' n=' num2str(cellnum) '/' num2str(size(CY5_AF594_TMR_tp{j},1))])
    set(gca,'Fontsize',14)
    %%% For max distance
        temp_max_pairwise_dists = max_pairwise_dists{j}(max_pairwise_dists{j}>0);
%     max(temp_max_pairwise_dists)
    cellnum = size(temp_max_pairwise_dists,1);
    figure(402);
    subplot(2,4,j)
    histogram(temp_max_pairwise_dists,edges1,'Normalization','probability')
    xticks(edges1(1:1:end))
     xticklabels(xticklabel1(1:1:end))
     if j >4
     xlabel('Distance (\mum)')
     end
     xtickangle(45)
     if j == 1 | j == 5
         ylabel('Probability')
     end
     ylim([0 .5])
    title([times0{j}{1} ' n=' num2str(cellnum) '/' num2str(size(CY5_AF594_TMR_tp{j},1))])
    set(gca,'Fontsize',14)
    %%%
    figure(403);
      subplot(2,4,j)
    temp_cells = CY5_AF594_TMR_tp{j}(CY5_AF594_TMR_tp{j}(:,44)>0,:)
    bins12 = 1:size(xticklabel1,2)
    for lm = bins12
        if sum(temp_cells(:,44)==lm) == 0
            temp_cells(size(temp_cells,1)+1,44) = lm;
        end
    end           
    h = boxplot(temp_cells(:,3),temp_cells(:,44),'Labels',xticklabel1);
%ylim([0 100])
xlabel(['Distance (\mum)'],'Fontsize',10)
ylabel(['Tsix Transcripts per Cell'],'Fontsize',10)
ylim([0 max(CY5_AF594_TMR(:,3))])
hold on
xtickangle(45)
figure(404);
      subplot(2,4,j)
    temp_cells = CY5_AF594_TMR_tp{j}(CY5_AF594_TMR_tp{j}(:,44)>0,:)
    bins12 = 1:size(xticklabel1,2)
    for lm = bins12
        if sum(temp_cells(:,44)==lm) == 0
            temp_cells(size(temp_cells,1)+1,44) = lm;
        end
    end           
    h = boxplot(temp_cells(:,10),temp_cells(:,44),'Labels',xticklabel1);
%ylim([0 100])
xlabel(['Distance (\mum)'],'Fontsize',10)
ylabel(['Tsix Transcription Sites per Cell'],'Fontsize',10)
ylim([0 max(CY5_AF594_TMR(:,10))])
hold on
xtickangle(45)
figure(405);
      subplot(2,4,j)
    temp_cells = CY5_AF594_TMR_tp{j}(CY5_AF594_TMR_tp{j}(:,44)>0,:)
    bins12 = 1:size(xticklabel1,2)
    for lm = bins12
        if sum(temp_cells(:,44)==lm) == 0
            temp_cells(size(temp_cells,1)+1,44) = lm;
        end
    end           
    h = boxplot(temp_cells(:,1),temp_cells(:,44),'Labels',xticklabel1);
%ylim([0 100])
xlabel(['Distance (\mum)'],'Fontsize',10)
ylabel(['Xist Transcripts per Cell'],'Fontsize',10)
if j <= 4
ylim([0 100])    
else
ylim([0 max(CY5_AF594_TMR(:,1))])
end
hold on
xtickangle(45)
figure(406);
      subplot(2,4,j)
    temp_cells = CY5_AF594_TMR_tp{j}(CY5_AF594_TMR_tp{j}(:,44)>0,:)
    bins12 = 1:size(xticklabel1,2)
    for lm = bins12
        if sum(temp_cells(:,44)==lm) == 0
            temp_cells(size(temp_cells,1)+1,44) = lm;
        end
    end           
    h = boxplot(temp_cells(:,8),temp_cells(:,44),'Labels',xticklabel1);
%ylim([0 100])
xlabel(['Distance (\mum)'],'Fontsize',10)
ylabel(['Xist Transcription Sites per Cell'],'Fontsize',10)
if j <=4
ylim([0 5])    
else
ylim([0 max(CY5_AF594_TMR(:,8))])
end
hold on
xtickangle(45)
%%% For max distance
        temp_max_pairwise_dists = max_pairwise_dists_cy5tmr{j}(max_pairwise_dists_cy5tmr{j}>0);
%     max(temp_max_pairwise_dists)
    cellnum = size(temp_max_pairwise_dists,1);
    figure(407);
    subplot(2,4,j)
    histogram(temp_max_pairwise_dists,edges3,'Normalization','probability')
    xticks(edges1(1:1:end))
     xticklabels(xticklabel1(1:1:end))
     if j >4
     xlabel('Distance (\mum)')
     end
     xtickangle(45)
     if j == 1 | j == 5
         ylabel('Probability')
     end
     ylim([0 .5])
    title([times0{j}{1} ' n=' num2str(cellnum) '/' num2str(size(CY5_AF594_TMR_tp{j},1))])
    set(gca,'Fontsize',14)
end
%% Distribution of transcription site intensity


%% Bar plot with sliding window
x_axis1 = 3;
y_axis1 = 1;
num_bins = 20;
start_valx = min(CY5_AF594_TMR(:,x_axis1));
 start_valy = min(CY5_AF594_TMR(:,y_axis1));
% end_valx = max(CY5_AF594_TMR(:,x_axis1));
end_valx = 20;
 end_valy = max(CY5_AF594_TMR(:,y_axis1));
bins1 = start_valx:(end_valx-start_valx)/num_bins:end_valx
window_radius = 2;
% figure(15); clf
% find(CY5_AF594_TMR_tp{i}(:,3) == 6)
% CY5_AF594_TMR_tp{i}(find(CY5_AF594_TMR_tp{i}(:,3) == 5),1)
%%% Median
% for i = 1:size(CY5_AF594_TMR_tp,2)
%     cell_bins = zeros(size(bins1)); %How many cells are in each bin
%     cell_bins = cell_bins+.0001;
%     yinbin1 = zeros(size(bins1)); %How much of other transcript in each bin
%     for j = 1:size(bins1,2)    
%     yinbin1(j) = median(CY5_AF594_TMR_tp{i}(find(CY5_AF594_TMR_tp{i}(:,3) == j),y_axis1))
%     cell_bins(j) = size(find(CY5_AF594_TMR_tp{i}(:,3) == j),1);
%     end
%     subplot(1,size(CY5_AF594_TMR_tp,2),i)
%     bar(bins1,yinbin1)
%     hold on
%     title(times0{i})
%     xlabel('Tsix')
%     ylabel('Median Xist')
%     yyaxis right
%     plot(bins1,cell_bins,'--','LineWidth',2);
%     ylabel('Cell Number')
% end
%%%Mean depending on timepoint
figure(16); clf
for i = 1:4%size(CY5_AF594_TMR_tp,2)
    cell_bins = zeros(size(bins1)); %How many cells are in each bin
    cell_bins = cell_bins+.0001;
    yinbin1 = zeros(size(bins1)); %How much of other transcript in each bin
    pvals1 = zeros(size(yinbin1));
    for j = 1:size(bins1,2)    
    yinbin1(j) = mean(CY5_AF594_TMR_tp{i}(find(CY5_AF594_TMR_tp{i}(:,x_axis1) >= j-window_radius & CY5_AF594_TMR_tp{i}(:,x_axis1) <= min([j+window_radius,end_valx])),y_axis1))
    cell_bins(j) = size(find(CY5_AF594_TMR_tp{i}(:,x_axis1) >= j-window_radius & CY5_AF594_TMR_tp{i}(:,x_axis1) <= min([j+window_radius,end_valx])),1);   
    end
    subplot(1,4,i)
    bar(bins1,yinbin1)
    hold on
    title(times0{i})
    xlabel('Tsix')
    ylabel('Mean Xist')
    yyaxis right
    plot(bins1,cell_bins,'--','LineWidth',2);
    ylabel('Cell Number')
    set(gca,'FontSize',15)
end

%% Randomization of data and the bar plot
figure(18); clf
num_rand = 10000;
x_axis1 = 3;
y_axis1 = 1;
for i = 1:4
    yinbin2 = zeros(size(bins1,1),size(bins1,2),num_rand); %How much of other transcript in each bin
        cell_bins = zeros(size(bins1)); %How many cells are in each bin
for k = 1:num_rand
    a1 = CY5_AF594_TMR_tp{i}(:,x_axis1);
    b1 = CY5_AF594_TMR_tp{i}(:,y_axis1);
    a1_rand = a1(randperm(length(a1)));
    for j = 1:size(bins1,2)    
    yinbin2(1,j,k) = mean(CY5_AF594_TMR_tp{i}(find(a1_rand >= j-window_radius & a1_rand <= min([j+window_radius,end_valx])),y_axis1));
    cell_bins(1,j,k) = size(find(a1_rand == j),1);
    end
end  
    cell_bins = cell_bins+.0001;
    subplot(1,4,i)
    bar(bins1,mean(yinbin2,3));     hold on
    errorbar(bins1, mean(yinbin2,3), std(yinbin2,0,3)./sqrt(num_rand), '.k');
    title(times0{i})
    xlabel('Tsix')
    ylabel('Mean Xist')
    yyaxis right
    plot(bins1,mean(cell_bins,3),'--','LineWidth',2);
    ylabel('Cell Number')
    prob1 = zeros(1,size(bins1,2));   %probability based off number of randoms
    for k = 1:size(bins1,2)
        prob1(k) = sum(yinbin2(1,k,:)<=yinbin1(1,k))/size(yinbin2,3);
    end
    prob1
end
    


 %% Nascent Transcription of TMR versus transcript numbers of CY5
 % NOTE: DO NOT USE TO DETERMINE A CORRELATION. RANDOMIZING LEADS TO
 % SIMILAR PLOTS
 nasc_tran = zeros(size(CY5_AF594_TMR,1),1); %,sum(CY5_AF594_TMR(:,13:14),2)
 for i = 1:size(CY5_AF594_TMR,1)
     temp_TMR = TMR_int(i,:);
     nasc_tran(i,1) = sum(temp_TMR(temp_TMR>TMR_tran_thres*TMR_cloudInt(1)))/TMR_cloudInt(1);
 end
 tsite_filter = CY5_AF594_TMR(:,10)<=4;
 nasc_tran1 = nasc_tran(tsite_filter);
 cy5_tran = CY5_AF594_TMR(:,1);
 cy5_tran1 = cy5_tran(tsite_filter);
 tsites1 = CY5_AF594_TMR(:,10);
  tsites2 = tsites1(tsite_filter);
 
%  figure(2000); clf; gscatter(nasc_tran,CY5_AF594_TMR(:,1),CY5_AF594_TMR(:,10))
 figure(2000); clf; gscatter(nasc_tran1,cy5_tran1, tsites2,'kbgrc','d',5)
 xlabel('Total Nascent Transcription of Tsix')
 ylabel('Number of Xist transcripts')
 xlim([-1 20])
 ylim([-1 80])
 set(gca,'FontSize',14)
 
 %%% Randomize
 num_rand = 4
 for i = 1: num_rand
     a1 = cy5_tran1;
    a1_rand = a1(randperm(length(a1)));
     figure(2000+i); clf; gscatter(nasc_tran1,a1_rand, tsites2,'kbgrc','d',5)
 xlabel('Total Nascent Transcription of Tsix')
 ylabel('Number of Xist transcripts')
 xlim([-1 20])
 ylim([-1 80])
 set(gca,'FontSize',14)
 end
    
 %% Bar plot nascent transcription versus other transcript
 nasc_tran
  nasc_tran = zeros(size(CY5_AF594_TMR,1),1); %,sum(CY5_AF594_TMR(:,13:14),2)
 for i = 1:size(CY5_AF594_TMR,1)
     temp_TMR = TMR_int(i,:);
     nasc_tran(i,1) = sum(temp_TMR(temp_TMR>TMR_tran_thres*TMR_cloudInt(1)))/TMR_cloudInt(1);
 end
 tsite_filter = CY5_AF594_TMR(:,10)<=4;
 nasc_tran1 = nasc_tran(tsite_filter);
 cy5_tran = CY5_AF594_TMR(:,1);
 cy5_tran1 = cy5_tran(tsite_filter);
 tsites1 = CY5_AF594_TMR(:,10);
  tsites2 = tsites1(tsite_filter);
 bin_num = 10;
 min_tmr_nasc = 0;
 max_tmr_nasc = 20;
 bins1 = min_tmr_nasc:(max_tmr_nasc-min_tmr_nasc)/bin_num:max_tmr_nasc;
 bins2 = (bins1(2:end)+bins1(1:size(bins1,2)-1))/2;
 mean_cy5 = zeros(size(bins2));
 for i = 1:size(bins2,2)
     mean_cy5(i) = mean(cy5_tran1(nasc_tran1 >= bins1(i) & nasc_tran1 < bins1(i+1)));
 end
 figure(2010); clf
 bar(bins2,mean_cy5)
 xlabel('Nascent Transcription')
 ylabel('Mean Xist Value')
 ylim([0 6])
 set(gca,'FontSize',14)
 
 %%% Randomize
 num_rand = 4
 for j = 1: num_rand
     a1 = cy5_tran1;
    a1_rand = a1(randperm(length(a1)));
     figure(2010+j); clf;
 mean_cy5 = zeros(size(bins2));
 for i = 1:size(bins2,2)
     mean_cy5(i) = mean(a1_rand(nasc_tran1 >= bins1(i) & nasc_tran1 < bins1(i+1)));
 end
 bar(bins2,mean_cy5)
 xlabel('Nascent Transcription')
 ylabel('Mean Xist Value')
 set(gca,'FontSize',14)
 ylim([0 6])
 end
 
 %% plot distributions before and after apparent threshold and do KS test
 pvals1 = zeros(10,4);
% figure(4000); clf
for j = 1:10
for i = 1:4
    cy5_trans = CY5_AF594_TMR_tp{i}(:,1);
    tmr_trans = CY5_AF594_TMR_tp{i}(:,3);
    [h,pvals1(j,i)] = kstest2(cy5_trans(CY5_AF594_TMR_tp{i}(:,3)<j),cy5_trans(CY5_AF594_TMR_tp{i}(:,3)>=j));
end
end
 pvals1
  pvals_cy5thresh = zeros(10,4);
  pvals_tmrthresh = zeros(10,4);
 figure(4000); clf
for j = 1:40
for i = 1%:4
    cy5_trans = CY5_AF594_TMR_tp{i}(:,1);
    tmr_trans = CY5_AF594_TMR_tp{i}(:,3);
%     [h,pvals_tmrthresh(j,i)] = kstest2(cy5_trans(tmr_trans<j),cy5_trans(tmr_trans>=j));
    [h,pvals_cy5thresh(j,i)] = kstest2(tmr_trans(cy5_trans<j),tmr_trans(cy5_trans>=j));
    figure(4001); clf;
    subplot(1,2,1);
    histogram(tmr_trans(cy5_trans<j))
    xlim([-1 20])
    subplot(1,2,2);
    histogram(tmr_trans(cy5_trans>=j))
    xlim([-1 20])
    title(['Threshold = ' num2str(j)])
    pause(1)
end
end
 pvals1
 pvals_cy5thresh
 pvals_trans_spear = zeros(2,7);
pvals_trans_Kend = zeros(2,7);
for i = 1:7
        cy5_trans = CY5_AF594_TMR_tp{i}(:,1);
    tmr_trans = CY5_AF594_TMR_tp{i}(:,3);
 %%% Spearman correlation
[pvals_trans_spear(1,i),pvals_trans_spear(2,i)] = corr(cy5_trans ,tmr_trans,'Type','Spearman')
[pvals_trans_Kend(1,i),pvals_trans_Kend(2,i)] = corr(cy5_trans,tmr_trans,'Type','Kendall')
%%%                   
end      
%% Look at how many transcription sites determined for each threshold of transcription site determination
% figure(4000); clf
% figure(4001); clf
threshes = {'1.5','2','2.5','3','3.5'};
threshes0 = .5:.5:3.5;
CY5_trans = zeros(size(threshes0));
TMR_trans = zeros(size(threshes0));
for i = 1:size(threshes0,1)
% load(['_Xist-CY5-th_Tsix-TMR-th_04-Mar-2019_3zadjust_6pixsame' threshes{i} 'transcriptsInTSite_CY5th1_TMRth1.mat'])
% figure(4000);
% subplot(1,size(size(threshes0,2),2),i)
Cloud_Tsite = CY5_int(:,1:2)>CY5_cloudInt(1)*5*threshes0(i);
Non_cloud_Tsite = CY5_int(:,3:end)>CY5_cloudInt(1)*threshes0(i);
CY5_trans(i) = (sum(Cloud_Tsite(:))+sum(Non_cloud_Tsite(:)))/size(find(CY5_int>0),1);
% histogram(CY5_AF594_TMR(:,8))
% figure(4001);
% subplot(1,size(threshes,2),i)
% histogram(CY5_AF594_TMR(:,10))
TMR_Tsites = TMR_int>TMR_cloudInt(1)*threshes0(i);
TMR_trans(i) = sum(TMR_Tsites(:))/size(find(TMR_int>0),1);
end
figure(4002); clf; plot(threshes0,CY5_trans)
figure(4003); clf; plot(threshes0,TMR_trans)