function [CY5_AF594_TMR,CY5_AF594_TMR_tp] = A1A_Cloud_Quant(dates0,times0,times1,max_tps,max_imnum,thCY5all,thAF594all,thTMRall,CY5_name,AF594_name,TMR_name,thNum_CY5,thNum_AF594,thNum_TMR,CY5_hascloud,AF594_hascloud,TMR_hascloud)
%%examples name:_STL1-TMR-th
% CY5_name = '_Xist-CY5-th';
% AF594_name = '_Tsix-AF594-th';
% TMR_name = '_Jpx-TMR-th';
% %%% copy over thresholds from experiment
% thCY5all = [12 15 18 21 24] 
% thAF594all = [11 13 15 17 19]
% thTMRall = [20 23 25 27 30] 
% %%% which threshold sets
% thNum_CY5 = 3; %Which set of the above thresholds to use (which column)
% thNum_TMR = 5; %Which set of the above thresholds to use (which column)
% thNum_AF594 = 3; %Which set of the above thresholds to use (which column)
% %%%
% dates0 = {'20161109','20161113','20161122'};
% times0 = {{'1hr','48hr','72hr'},{'96hr','120hr'},{'11d','13d'}}; %168 hour should be somewhere sometime
% times1 = [0,2,3,4,5,11,13];
% max_imnum = 6;
% max_tps = 3;
% %%%
% nucTH = 'low'
%% Processing data


avg_cloud_size = [0];
max_cloud_size = [0];
avg_trans_numb = [0];
avg_trans_numb_pos = [0]; %average transcript number in positive cells
median_trans_fl = [0];
percent_cloud = [0];
percent_Xist = [0];
cloud_distributions = {};
All_spots = NaN(200,1);
All_non_cloud = NaN(200,1);
All_non_cloud_tp = NaN(500,200,8);
All_cloud = NaN(200,1);
All_cloud2 = NaN(1,2);
counter = 1; %counter for timepoint
counter2 = 1; %counter for all spot intensities
counter3 = 1; %counter for noncloud spot intensity
counter4 = 1; %counter for cloud spots
ints = zeros(8,1);

for i = 1:size(dates0,2)
    for j = 1:max_tps
        counter1 = 1; %counter for image
        temp_trans_numb = NaN(100,5); %Matrix with transcript number for each cell in each image for this timepoint
        temp_cloud_size = NaN(100,5); %Matrix with cloud size for each cell in each image for this timepoint
        %temp_trans_fl = NaN(100,5); %Matrix with non-cloud spot intensities
        temp_cell_count = 0; %keeps track of total number of cells
        counter5 = 1; %counter for spot intensities to keep all images in timepoint in same z
        fail_count = 0;   %Skips timepoint if all images fail
        for k = 1:max_imnum
            try
                %filename = ['mRNA_' dates0{i} '_STL1-TMR-th25_AF594-GPP2-th15_CTT1-GPD1-th18_' times0{i}{j} '_im' num2str(k)];
                %filename = ['mRNA_' dates0{i} TMR_name num2str(thTMRall(thNum_CY5)) AF594_name num2str(thAF594all(thNum_CY5)) CY5_name num2str(thCY5all(thNum_CY5)) '_' times0{i}{j} '_im' num2str(k)];
                filename = ['mRNA_' dates0{i} '_F1-2-1' TMR_name{i} num2str(thTMRall{i}(thNum_CY5)) '_' num2str(thAF594all{i}(thNum_CY5)) CY5_name{i} num2str(thCY5all{i}(thNum_CY5)) '_' times0{i}{j} '_im' num2str(k)];
                load(filename,['PARtmr_' nucTH]);
                         for m = 1: size(PARtmr_low.TotExpInt,1) %Go through every cell in this image 
                        %%%Eliminate dead pixels  (needs to be improved)                      
                        for p = 1:size(PARtmr_low.xinit,2)
                                same_pos = find(PARtmr_low.xinit(m,:) == PARtmr_low.xinit(m,p));
                                if size(same_pos,2)>5
                                    for q = 1:size(same_pos,2)
                                        PARtmr_low.xinit(m,same_pos(q)) = NaN;
                                        PARtmr_low.yinit(m,same_pos(q)) = NaN;
                                        PARtmr_low.zinit(m,same_pos(q)) = NaN;
                                        PARtmr_low.TotExpInt(m,same_pos(q)) = NaN;
                                    end
                                end 
                        end
                        %%% Below removes duplicates of the same
                        %%% transcript. Goes through every transcript,
                        %%% finds distance between it and other
                        %%% transcripts, and eliminates others 
                        dists = zeros(1,size(size(PARtmr_low.xinit,2),3)); %stores distance from transcript to others
                        for p = 1:size(PARtmr_low.xinit,2)
                            if not(isnan(PARtmr_low.xinit(m,p)))
                            for pp = 1:size(PARtmr_low.xinit,2)    %find distance squared
                                dists(1,pp,1) = (PARtmr_low.xinit(m,p)- PARtmr_low.xinit(m,pp))^2;
                                dists(1,pp,2) = (PARtmr_low.yinit(m,p)- PARtmr_low.yinit(m,pp))^2;
                                dists(1,pp,3) = ((PARtmr_low.zinit(m,p)- PARtmr_low.zinit(m,pp))*3)^2;
                            end
                            distF = sqrt(dists(1,:,1)+dists(1,:,2)+dists(1,:,3));
                            for pp = 1:size(PARtmr_low.xinit,2)
                                if p ~= pp & distF(pp) < 3
                                    PARtmr_low.xinit(m,pp) = NaN;
                                    PARtmr_low.yinit(m,pp) = NaN;
                                    PARtmr_low.zinit(m,pp) = NaN;
                                    PARtmr_low.TotExpInt(m,pp) = NaN;
                                end
                            end
                            end
                        end   
                        %%%
                         end
                ints(i,1) = ints(i,1)+ size(find(PARtmr_low.TotExpInt >8000),1);
                All_spots(1:size(PARtmr_low.TotExpInt,1),counter2:counter2-1+size(PARtmr_low.TotExpInt,2)) = PARtmr_low.TotExpInt;
                All_non_cloud(1:size(PARtmr_low.TotExpInt,1),counter3:counter3-3+size(PARtmr_low.TotExpInt,2)) = PARtmr_low.TotExpInt(:,3:end);          %RNAs without the cloud
                All_non_cloud_tp(counter5:counter5+size(PARtmr_low.TotExpInt,1)-1,1:size(PARtmr_low.TotExpInt,2)-2,counter) = PARtmr_low.TotExpInt(:,3:end); 
                All_cloud(1:size(PARtmr_low.TotExpInt,1),counter4) = PARtmr_low.TotExpInt(:,1);
                All_cloud(size(PARtmr_low.TotExpInt,1)+1:size(PARtmr_low.TotExpInt,1)*2,counter4) = PARtmr_low.TotExpInt(:,2);
                counter2 =counter2+size(PARtmr_low.TotExpInt,2);
                counter3 = counter3+size(PARtmr_low.TotExpInt,2)-1;
                counter4 = counter4+1;
                counter5 = counter5+size(PARtmr_low.TotExpInt,1);
                if size(PARtmr_low.TotExpInt,2) >1 & size(find(PARtmr_low.TotExpInt(:,3:end) <6000),1) > 0  %If there is a cell with more than one RNA spot and they are less than 4000
                    non_cloud = PARtmr_low.TotExpInt(:,3:end);          %RNAs without the cloud
                    RNA_int = nanmedian(non_cloud(find(PARtmr_low.TotExpInt(:,3:end))))  %RNA intensity for one spot is median of non-clod spots under 4000 intensity
                else
                    RNA_int = 1887;   %Default intensity from day 13 image 1
                end
                for m = 1:size(PARtmr_low.TotExpInt,1)      %loop for number of cells in image
                    if PARtmr_low.TotExpInt(m,1) > 10*RNA_int
                        temp_cloud_size(m,counter1) = PARtmr_low.TotExpInt(m,1)/RNA_int;
                        temp_trans_numb(m,counter1) = PARtmr_low.TotExpInt(m,1)/RNA_int+size(find(PARtmr_low.TotExpInt(m,3:end)),2);
                       % temp_trans_fl(m,counter1) = PARtmr_low.TotExpInt(m,3:end)
                    else
                        temp_trans_numb(m,counter1) = size(find(PARtmr_low.TotExpInt(m,:)<6000),2);
                    end
                end
                counter1 = counter1+1;
                temp_cell_count = temp_cell_count+size(PARtmr_low.TotExpInt,1);
            catch
                fail_count = fail_count+1;
            end
        end
        if fail_count <k
            temp_trans_numb;
            avg_cloud_size(counter) = nanmean(temp_cloud_size(:));
            max_cloud_size(counter) = nanmax(temp_cloud_size(:));
            avg_trans_numb(counter) = nansum(temp_trans_numb(:))/temp_cell_count;
            percent_cloud(counter) = size(find(temp_cloud_size(:)>0),1)/temp_cell_count;
            percent_Xist(counter) = size(find(temp_trans_numb(:)>0),1)/temp_cell_count;
            cloud_distributions{counter} = temp_cloud_size(temp_cloud_size>0);
            counter = counter+1;
        else
            ['Timepoint ' times0{i}{j} ' failed']
        end
    end
end
xs1 = times1;       
avg_cloud_size
max_cloud_size
avg_trans_numb
percent_cloud
percent_Xist
cloud_distributions
All_non_cloud_tp(find(All_non_cloud_tp == 0)) = NaN;

%%

xs1 = [0,.25,.5,1,2,5,9];
figure(2);clf;plot(xs1,max_cloud_size,'Color','b');hold on; plot(xs1,avg_cloud_size,'Color','g'); legend('Maximum Cloud Size','Average Cloud Size')
xlabel('Day')
ylabel('Number of Xist transcripts')
%plot cloud distribution
figure(1); clf; hold on
colors = {'k','r','g','b','m','c',[1,.5,.5],[.75,.75,.25]};
min1 = floor(nanmin(max_cloud_size)/100)*100;
max1 = ceil(nanmax(max_cloud_size)/100)*100;
bins = min1:50:max1;
cloud_hists = {};
xs = zeros(1,size(bins,2)-1);
ys = zeros(1,size(bins,2)-1);
for i = 1:size(avg_cloud_size,2)    %go through every timepoint
   %cloud_distributions{i}
    for j = 1:(size(bins,2)-1)
        xs(j) = mean(bins(j:j+1));
        if size(cloud_distributions{i},1)>0
            ys(j) = (size(find(cloud_distributions{i} < bins(j+1)),1)-size(find(cloud_distributions{i}<=bins(j)),1))/size(cloud_distributions{i},1);
        else
            ys(j) = 0;
        end
    end
    plot(xs,ys,'Color',colors{i},'LineWidth',2)
end
legend('0 days','2 days','3 days','4 days','5 days','11 days','13 days')

%plot all spot distribution
figure(3); clf; hold on
colors = {'k','r','g','b','m','c',[1,.5,.5],[.75,.75,.25]};
for i = 1:size(All_spots,1)
    for j = 1:size(All_spots,2)
        if All_spots(i,j) == 0 | All_spots(i,j) <0
            All_spots(i,j) = NaN;
        end
    end
end
min1 = floor(nanmin(All_spots(:))/100)*100;
max1 = ceil(nanmax(All_spots(:))/100)*100;
%min1 = 0;
%max1 = 20000;
bins = min1:1000:max1;
xs = zeros(1,size(bins,2)-1);
ys = zeros(1,size(bins,2)-1);
total1 = size(find(All_spots(:) >0),1);
for j = 1:(size(bins,2)-1)  %go through every bin
    xs(j) = mean(bins(j:j+1));  %x value is the middle of the bin
    
    ys(j) = (size(find(All_spots(:) < bins(j+1)),1)-size(find(All_spots(:)<=bins(j)),1));
    
end
sum(ys)
total1
ys = ys/sum(ys)/10;
plot(xs,ys,'LineWidth',3,'Color','b')

%legend('All Spots','2 days','3 days','4 days','5 days','11 days','13 days')

%% plot all non cloud spot distribution
figure(4); clf; hold on
for i = 1:size(All_non_cloud,1)
    for j = 1:size(All_non_cloud,2)
        if All_non_cloud(i,j) == 0 | All_non_cloud(i,j) <0
            All_non_cloud(i,j) = NaN;
        end
    end
end
min1 = floor(nanmin(All_non_cloud(:))/100)*100;
max1 = 50000;
%max1 = ceil(nanmax(All_non_cloud(:))/100)*100;
bins = min1:100:max1;
xs = zeros(1,size(bins,2)-1);
ys = zeros(1,size(bins,2)-1);
total2 = size(find(All_non_cloud(:) >0),1);
for j = 1:(size(bins,2)-1)  %go through every bin
    xs(j) = mean(bins(j:j+1));  %x value is the middle of the bin
    
    ys(j) = (size(find(All_non_cloud(:) < bins(j+1)),1)-size(find(All_non_cloud(:)<=bins(j)),1));
    
end
sum(ys)
total2
ys = ys/sum(ys);
plot(xs,ys,'LineWidth',3,'Color','r')

legend('Non-cloud spots','3 days','4 days','5 days','11 days','13 days')

%plot all cloud spot fluorescence distribution
%figure(5); clf; hold on
for i = 1:size(All_cloud,1)
    for j = 1:size(All_cloud,2)
        if All_cloud(i,j) == 0 | All_cloud(i,j) <1887*10
            All_cloud(i,j) = NaN;
        end
    end
end
min1 = floor(nanmin(All_cloud(:))/100)*100;
max1 = ceil(nanmax(All_cloud(:))/100)*100;
bins = min1:1000:max1;
xs = zeros(1,size(bins,2)-1);
ys = zeros(1,size(bins,2)-1);
total2 = size(find(All_cloud(:) >0),1);
for j = 1:(size(bins,2)-1)  %go through every bin
    xs(j) = mean(bins(j:j+1));  %x value is the middle of the bin
    ys(j) = (size(find(All_cloud(:) < bins(j+1)),1)-size(find(All_cloud(:)<=bins(j)),1)); 
end
sum(ys)
total2
ys = ys/sum(ys)/10;
plot(xs,ys,'LineWidth',3,'Color','b')
xlabel('Fluorescence')
ylabel('Probability')
legend('Non-cloud spots','Clouds','3 days','4 days','5 days','11 days','13 days')

%plot all cloud transcript distribution
figure(6); clf; hold on
bin_size = 50;
counter2 = 1;
for i = 1:size(cloud_distributions,2)
    for j = 1:size(cloud_distributions{i},1)
        All_cloud2(counter2,1) = cloud_distributions{i}(j);
        counter2 = counter2+1;
    end
end
min1 = floor(nanmin(All_cloud2(:))/bin_size)*bin_size;
max1 = ceil(nanmax(All_cloud2(:))/bin_size)*bin_size;
bins = min1:bin_size:max1;
xs = zeros(1,size(bins,2)-1);
ys = zeros(1,size(bins,2)-1);
total2 = size(find(All_cloud2(:) >0),1);
for j = 1:(size(bins,2)-1)  %go through every bin
    xs(j) = mean(bins(j:j+1));  %x value is the middle of the bin
    ys(j) = (size(find(All_cloud2(:) < bins(j+1)),1)-size(find(All_cloud2(:)<=bins(j)),1)); 
end
sum(ys)
total2
ys = ys/sum(ys);
plot(xs,ys,'LineWidth',3,'Color','b')
xlabel('Xist Transcripts')
ylabel('Probability')
legend('Clouds','Non-cloud spots','3 days','4 days','5 days','11 days','13 days')

%Plot median of non-cloud spot intensities over time

figure(7); clf; hold on
bin_size = 50;
xs = [0,.25,.5,1,2,5,9,9];
ys = [0];
ys1 = [0];
counter2 = 1;
for i = 1:size(All_non_cloud_tp,3)
    temp_mat = All_non_cloud_tp(:,:,i);
    ys(i) = nanmedian(temp_mat(:));
    ys1(i) = nanmean(temp_mat(:));
end
plot(xs,ys,'LineWidth',3,'Color','b')
hold on
plot(xs,ys1,'LineWidth',3,'Color','r')
xlabel('Time (Days)')
ylim([0 30000])
ylabel('Fluorescence')
legend('Median Spot intensity','Average Spot Intensity','3 days','4 days','5 days','11 days','13 days')
