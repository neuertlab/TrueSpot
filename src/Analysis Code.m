%% This is to plot comparisons of threshold vs how many spots would be detected
Ych2 = 2;
exp_date = '20161109';                                                        % Date of the experiment and folder name on the hard disk
strain = 'mESC';
osmo = '0.2Mstep';
genepairs = {'no-probes','SCR-CY5-SCR-AF594-SCR-TMR','STL1-TMR-GPD1-CY5-GPP2-AF594'};
file_end = '_MMStack.ome.tif';
exp_names = {'1hr','48hr' '72hr'}; %'2min' '4min' '6min' '8min' '10min' '12min' '14min',...
    %'16min' '18min' '20min' '22min' '25min' '40min' '50min' '120min' '122min' '129min'};
Ychs = {'3','4','5'};
Ych_names = {'Cy5','AF594','TMR'};
colors = {'r','b','g','m','c'};
styles = {'--',':','-','-.','..'};
figure(101);clf
hold on
Ych = 2;
pairs = [1 1;1 2; 1 3; 2 1;2 2;2 3; 3 3]                   %timepoint, genepair 
counter = 1;                        %counter for legends
legends = {};                       %Where legend entries will be recorded.
for m = 2                       %for Ych
    m
    for k = 1:2         %for image num
        k
        Ych2 = str2num(Ychs{m});
        for i = [2 5]                  %Uses the pairs matrix to look at each
            i
            exp_name = [exp_date '_' exp_names{pairs(i,1)} '_' strain  '_' genepairs{pairs(i,2)} '_img']; 
            load(['SpecRNA ' exp_name num2str(k) '_Ych2_' num2str(Ych2) '.mat'],'SpecRNA')
            SpecRNA = SpecRNA(:,1:100);
            thresholds = SpecRNA(1,1:100);
            %         RNA_counts = zeros(1);            %for total RNA counts in image
            %         for j = 1:100
            %             RNA_counts(1,j) = sum(SpecRNA(2:size(SpecRNA,1),j));
            %         end
            spotpcell = zeros(1);
            for j = 1:100
                spotpcell(1,j) = sum(SpecRNA(2:size(SpecRNA,1),j))/(size(SpecRNA,1)-1);
            end
            %spotpcell = asinh(spotpcell);
            %       spotpcell = spotpcell/sum(spotpcell);               %Normalizes to total number of spots in all thresholds
            plot(thresholds,log(spotpcell),'Color',colors{i}) %,'Linestyle',styles{k}
            legends{counter} = [exp_names{pairs(i,1)} ' ' genepairs{pairs(i,2)} ' channel ' Ych_names{m} 'image num ' num2str(k)];
            counter = counter+1;
        end
    end
end
legend(legends)

%% This is to look at the pixel distributions in the filtered or unfiltered (change line 74) images

Ych2 = 2;
exp_date = '20161109';                                                        % Date of the experiment and folder name on the hard disk
strain = 'mESC';
osmo = '0.2Mstep';
genepairs = {'no-probes','SCR-CY5-SCR-AF594-SCR-TMR','Xist-CY5-Tsix-AF594-Jpx-TMR'};
genepair_names = {'no probes','scrambled probes','gene probes'};
probe_names = {'Xist CY5','Tsix AF594','Jpx TMR'};
file_end = '_MMStack.ome.tif';
exp_names = {'1hr','48hr' '72hr'}; %'2min' '4min' '6min' '8min' '10min' '12min' '14min',...
Ychs = {'3','4','5'};
Ych_names = {'Cy5','AF594','TMR'};
colors = {'r','b','m','g','c'};
styles = {':','--','-','-.','..'};
figure(101);clf
hold on
Ych = 2;
pairs = [1 1;1 2; 1 3; 2 1;2 2;2 3; 3 3]                   %timepoint, genepair
counter = 1;                        %counter for legends
legends = {};                       %Where legend entries will be recorded.
for m = 1:3                       %for Ych
    subplot(1,3,m)
    hold on
    counter = 1;                        %counter for legends
    Ych2 = str2num(Ychs{m});
    for i = [1:3]                  %Uses the pairs matrix to look at each
        for k = 1 %img_num
            exp_name = [exp_date '_' exp_names{pairs(i,1)} '_' strain  '_' genepairs{pairs(i,2)} '_img'];
            load(['PixRNA_unfiltered_' exp_name num2str(k) '_Ych2_' num2str(Ych2) '.mat'],'PixRNA')
            temp12 = PixRNA(1:2048*2048);      %linear version of image
            minp = min(temp12);
            maxp = max(temp12);
            if m == 1
                xval = 100:5:200;   %Sets x values used for bins
                xvalp = 102.5:5:197.5;    %Sets the graphed x values in the middle of the bins
            elseif m == 2
                xval = 100:2:200;   %Sets x values used for bins
                xvalp = 101:2:199;    %Sets the graphed x values in the middle of the bins
            elseif m == 3
                xval = 100:5:350;   %Sets x values used for bins
                xvalp = 102.5:5:347.5;    %Sets the graphed x values in the middle of the bins
            end
            yvalp = [0];
            for j = 1:size(xval,2)-1
                yvalp(1,j) = size(find(temp12 <= xval(j+1)),2)-size(find(temp12 <= xval(j)),2);     %sets y value for bin to anything that is greater than to the lower number and less than or equal to the higher number
            end
            %       yvalp(1,j) = yvalp(1,j) + size(find(temp12 == xval(j)),1);          %Adds the last values
            yvalp = yvalp/sum(yvalp);
            plot(xvalp,yvalp,'Color',colors{m},'Linestyle',styles{pairs(i,2)},'Linewidth',2)   %
            legends{counter} = [exp_names{pairs(i,1)} ' ' genepair_names{pairs(i,2)}];% ' channel ' Ych_names{m}];% 'image num ' num2str(k)];
            counter = counter+1;
            mean12 = mean(temp12(find(temp12>0)))
            stdv12 = std(double(temp12(find(temp12>0))))
        end
    end
    xlabel('Pixel Intensity');ylabel('Probability within a Cell');
    legend(legends)
    legends = {};                       %Where legend entries will be recorded.
    title(probe_names{m});
end


Ych2 = 2;
exp_date = '20161109';                                                        % Date of the experiment and folder name on the hard disk
strain = 'mESC';
osmo = '0.2Mstep';
genepairs = {'no-probes','SCR-CY5-SCR-AF594-SCR-TMR','Xist-CY5-Tsix-AF594-Jpx-TMR'};
genepair_names = {'no probes','scrambled probes','gene probes'};
probe_names = {'Xist CY5','Tsix AF594','Jpx TMR'};
file_end = '_MMStack.ome.tif';
exp_names = {'1hr','48hr' '72hr'}; %'2min' '4min' '6min' '8min' '10min' '12min' '14min',...
    %'16min' '18min' '20min' '22min' '25min' '40min' '50min' '120min' '122min' '129min'};
Ychs = {'3','4','5'};
Ych_names = {'Cy5','AF594','TMR'};
colors = {'r','b','m','g','c'};
styles = {':','--','-','-.','..'};
figure(101);clf
hold on
Ych = 2;
pairs = [1 1;1 2; 1 3; 2 1;2 2;2 3; 3 3]                   %timepoint, genepair 
counter = 1;                        %counter for legends
legends = {};                       %Where legend entries will be recorded.
for m = 2                       %for Ych
        Ych2 = str2num(Ychs{m});
    for i = [2:3]                  %Uses the pairs matrix to look at each
        for k = 1       
        exp_name = [exp_date '_' exp_names{pairs(i,1)} '_' strain  '_' genepairs{pairs(i,2)} '_img']; 
        load(['PixRNA_unfiltered_' exp_name num2str(k) '_Ych2_' num2str(Ych2) '.mat'],'PixRNA')
        temp12 = PixRNA(1:2048*2048);      %linear version of image
        minp = min(temp12);
        maxp = max(temp12);
%         binnum12 = 10;
        if m == 1
        xval = 100:5:200;   %Sets x values used for bins
        xvalp = 102.5:5:197.5;    %Sets the graphed x values in the middle of the bins
        elseif m == 2
        xval = 100:5:200;   %Sets x values used for bins
        xvalp = 102.5:5:197.5;    %Sets the graphed x values in the middle of the bins
        elseif m == 3
        xval = 100:5:350;   %Sets x values used for bins
        xvalp = 102.5:5:347.5;    %Sets the graphed x values in the middle of the bins
        end
        yvalp = [0];
%         xval = minp:(maxp-minp)/binnum12:maxp;          %Sets x values used for bins according to bin size
%         xvalp = minp+(maxp-minp)/binnum12/2:(maxp-minp)/binnum12:maxp-(maxp-minp)/binnum12/2;   %Sets the graphed x values in the middle of the bins
        for j = 1:size(xval,2)-1
            yvalp(1,j) = size(find(temp12 <= xval(j+1)),2)-size(find(temp12 <= xval(j)),2);     %sets y value for bin to anything that is greater than to the lower number and less than or equal to the higher number
        end
        %       yvalp(1,j) = yvalp(1,j) + size(find(temp12 == xval(j)),1);          %Adds the last values
        yvalp = yvalp/sum(yvalp);
        yvalp1(i-1,1:size(yvalp,2)) = yvalp;
        plot(xvalp,yvalp,'Color',colors{m},'Linestyle',styles{pairs(i,2)},'Linewidth',2)   % 
        legends{counter} = [exp_names{pairs(i,1)} ' ' genepair_names{pairs(i,2)}];% ' channel ' Ych_names{m}];% 'image num ' num2str(k)];
        counter = counter+1;
        mean12 = mean(temp12(find(temp12>0)))
        stdv12 = std(double(temp12(find(temp12>0))))
        end
    end
end
diff1 = yvalp1(2,:)-yvalp1(1,:);
xvalp
xvalp(10+find(diff1(10:20) == max(diff1(10:20))))









%% This is to generate a histogram depicting how high a threshold images needed to be in order to have a percent of their cells with 0 RNA spots
%load('RNA_thresholds_20130914')
timmins = [0,15];
%no probes 0 min
thress = [0,0];                 %one row that has the lowest thresholds that result in the percentage of cells having 0 spots. Each value corresponds to one image
thres0 = zeros(2,1);            %first row is threshold, next rows are each image's percent of cells with transcripts
prct = .95;                     %cutoff percent of cells with 0 spots
counter19 = 2;                  %counter for image number (thres0)
for i = [1:5]       %Task ids 1-3, timepoints 0,1,2
    i
    SpecRNA = importdata(strcat('SpecRNA 20160206_0.4Mstep_WT_00min_no_probes_img',num2str(i)),'SpecRNA');            %SpecRNA seems to have a list of how many RNA spots were determined for each cell in each image for each threshold.
    %First row of Spec RNA has the threshold, every subsequent row is for a
    %cell in the image. Any value on (x+1,y) is the number of RNA spots
    %for cell x for threshold y(1).
    SpecRNA = SpecRNA(:,1:100);
    thres0 = zeros(2,size(SpecRNA,2));                              %thres0 will have the percent of cells without transcripts for each threshold. First row is thresholds and second row is the decimal fraction
    thres0(1,:) = SpecRNA(1,:);
    for j = 1:size(SpecRNA,2)
        j
        thres0(2,j) = numel(find(SpecRNA(2:size(SpecRNA,1),j)==0))/...
            (size(SpecRNA,1)-1);                                      %the number of cells that have 0 spots divided by the total for each threshold          
    end
    if ~isempty(thres0(1,find(thres0(2,:)>=prct,1,'first')))            %Make sure there is a threshold at which there are enough cells with 0 spots
        thress(counter19-1) = thres0(1,find(thres0(2,:)>=prct,1,'first'));  %If there is a threshold that results in enough off cells, use lowest such threshold
    else
        thress(counter19-1) = max(thres0(1,:));                             %If no threshold with enough off cells, use maximum threshold
    end
    counter19 = counter19+1
end
figure(16);clf;hist(thress,6); xlabel('Threshold');ylabel('Image Number');
axis([0 max(thress)+(min(thress)-0) 0 numel(find(thress == mode(thress)))+3])
title(strcat('Thresholds for ', num2str(prct*100), ' percent cells with 0 transcripts at 0 min timepoint with no probes'))
mean1 = mean(thress)
stdv1 = std(thress)
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',['mean = ' num2str(mean1) ', stdv = ' num2str(stdv1)],'FitBoxToText','on');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','black')


%no probes 15 min

thress = [0,0];                 %one row that has the lowest thresholds that result in the percentage of cells having 0 spots. Each value corresponds to one image
thres0 = zeros(2,1);            %first row is threshold, next rows are each image's percent of cells with transcripts
counter19 = 2;                  %counter for image number (thres0)
for i = [6:10]       %Task ids 1-3, timepoints 0,1,2
    i
    SpecRNA = importdata(strcat('SpecRNA 20160206_0.4Mstep_WT_15min_no_probes_img',num2str(i)),'SpecRNA');            %SpecRNA seems to have a list of how many RNA spots were determined for each cell in each image for each threshold.
    %First row of Spec RNA has the threshold, every subsequent row is for a
    %cell in the image. Any value on (x+1,y) is the number of RNA spots
    %for cell x for threshold y(1).
    SpecRNA = SpecRNA(:,1:100);
    thres0 = zeros(2,size(SpecRNA,2));                              %thres0 will have the percent of cells without transcripts for each threshold. First row is thresholds and second row is the decimal fraction
    thres0(1,:) = SpecRNA(1,:);
    for j = 1:size(SpecRNA,2)
        j
        thres0(2,j) = numel(find(SpecRNA(2:size(SpecRNA,1),j)==0))/...
            (size(SpecRNA,1)-1);                                      %the number of cells that have 0 spots divided by the total for each threshold          
    end
    if ~isempty(thres0(1,find(thres0(2,:)>=prct,1,'first')))            %Make sure there is a threshold at which there are enough cells with 0 spots
        thress(counter19-1) = thres0(1,find(thres0(2,:)>=prct,1,'first'));  %If there is a threshold that results in enough off cells, use lowest such threshold
    else
        thress(counter19-1) = max(thres0(1,:));                             %If no threshold with enough off cells, use maximum threshold
    end
    counter19 = counter19+1
end
figure(17);clf;hist(thress,6); xlabel('Threshold');ylabel('Image Number');
axis([0 max(thress)+(min(thress)-0) 0 numel(find(thress == mode(thress)))+3])
title(strcat('Thresholds for ', num2str(prct*100), ' percent cells with 0 transcripts at 15 min timepoint with no probes'))
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','black')
mean1 = mean(thress)
stdv1 = std(thress)
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',['mean = ' num2str(mean1) ', stdv = ' num2str(stdv1)],'FitBoxToText','on');

% Both together

thress = [0,0];                 %one row that has the lowest thresholds that result in the percentage of cells having 0 spots. Each value corresponds to one image
thres0 = zeros(2,1);            %first row is threshold, next rows are each image's percent of cells with transcripts
counter19 = 2;                  %counter for image number (thres0)
for i = [1:10]       %Task ids 1-3, timepoints 0,1,2
    i
    if i < 6
        SpecRNA = importdata(strcat('SpecRNA 20160206_0.4Mstep_WT_00min_no_probes_img',num2str(i)),'SpecRNA');            %SpecRNA seems to have a list of how many RNA spots were determined for each cell in each image for each threshold.
    else
        SpecRNA = importdata(strcat('SpecRNA 20160206_0.4Mstep_WT_15min_no_probes_img',num2str(i)),'SpecRNA');            %SpecRNA seems to have a list of how many RNA spots were determined for each cell in each image for each threshold.
    end
        %First row of Spec RNA has the threshold, every subsequent row is for a
    %cell in the image. Any value on (x+1,y) is the number of RNA spots
    %for cell x for threshold y(1).
    SpecRNA = SpecRNA(:,1:100);
    thres0 = zeros(2,size(SpecRNA,2));                              %thres0 will have the percent of cells without transcripts for each threshold. First row is thresholds and second row is the decimal fraction
    thres0(1,:) = SpecRNA(1,:);
    for j = 1:size(SpecRNA,2)
        j
        thres0(2,j) = numel(find(SpecRNA(2:size(SpecRNA,1),j)==0))/...
            (size(SpecRNA,1)-1);                                      %the number of cells that have 0 spots divided by the total for each threshold          
    end
    if ~isempty(thres0(1,find(thres0(2,:)>=prct,1,'first')))            %Make sure there is a threshold at which there are enough cells with 0 spots
        thress(counter19-1) = thres0(1,find(thres0(2,:)>=prct,1,'first'));  %If there is a threshold that results in enough off cells, use lowest such threshold
    else
        thress(counter19-1) = max(thres0(1,:));                             %If no threshold with enough off cells, use maximum threshold
    end
    counter19 = counter19+1
end
figure(18);clf;hist(thress,6); xlabel('Threshold');ylabel('Image Number');
axis([0 max(thress)+(min(thress)-0) 0 numel(find(thress == mode(thress)))+5])
title(strcat('Thresholds for ', num2str(prct*100), ' percent cells with 0 transcripts at 15 min timepoint with no probes'))
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','black')
mean1 = mean(thress)
stdv1 = std(thress)
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',['mean = ' num2str(mean1) ', stdv = ' num2str(stdv1)],'FitBoxToText','on');

%% Using scrambled Probes

thress = [0,0];                 %one row that has the lowest thresholds that result in the percentage of cells having 0 spots. Each value corresponds to one image
thres0 = zeros(2,1);            %first row is threshold, next rows are each image's percent of cells with transcripts
counter19 = 2;                  %counter for image number (thres0)
for i = [1:4]       %Task ids 1-3, timepoints 0,1,2
    i
    SpecRNA = importdata(strcat('SpecRNA 20160206_0.4Mstep_WT_00min_SCR-TMR-3500_img',num2str(i)),'SpecRNA');            %SpecRNA seems to have a list of how many RNA spots were determined for each cell in each image for each threshold.
    %First row of Spec RNA has the threshold, every subsequent row is for a
    %cell in the image. Any value on (x+1,y) is the number of RNA spots
    %for cell x for threshold y(1).
    SpecRNA = SpecRNA(:,1:100);
    thres0 = zeros(2,size(SpecRNA,2));                              %thres0 will have the percent of cells without transcripts for each threshold. First row is thresholds and second row is the decimal fraction
    thres0(1,:) = SpecRNA(1,:);
    for j = 1:size(SpecRNA,2)
        j
        thres0(2,j) = numel(find(SpecRNA(2:size(SpecRNA,1),j)==0))/...
            (size(SpecRNA,1)-1);                                      %the number of cells that have 0 spots divided by the total for each threshold          
    end
    if ~isempty(thres0(1,find(thres0(2,:)>=prct,1,'first')))            %Make sure there is a threshold at which there are enough cells with 0 spots
        thress(counter19-1) = thres0(1,find(thres0(2,:)>=prct,1,'first'));  %If there is a threshold that results in enough off cells, use lowest such threshold
    else
        thress(counter19-1) = max(thres0(1,:));                             %If no threshold with enough off cells, use maximum threshold
    end
    counter19 = counter19+1
end
figure(19);clf;hist(thress,6); xlabel('Threshold');ylabel('Image Number');
axis([0 max(thress)+(min(thress)-0) 0 numel(find(thress == mode(thress)))+3])
title(strcat('Thresholds for ', num2str(prct*100), ' percent cells with 0 transcripts at 0 min timepoint with scrambled probes'))
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','black')
mean1 = mean(thress)
stdv1 = std(thress)
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',['mean = ' num2str(mean1) ', stdv = ' num2str(stdv1)],'FitBoxToText','on');

%no probes 15 min

thress = [0,0];                 %one row that has the lowest thresholds that result in the percentage of cells having 0 spots. Each value corresponds to one image
thres0 = zeros(2,1);            %first row is threshold, next rows are each image's percent of cells with transcripts
counter19 = 2;                  %counter for image number (thres0)
for i = [5:9]       %Task ids 1-3, timepoints 0,1,2
    i
    SpecRNA = importdata(strcat('SpecRNA 20160206_0.4Mstep_WT_15min_SCR-TMR-3500_img',num2str(i)),'SpecRNA');            %SpecRNA seems to have a list of how many RNA spots were determined for each cell in each image for each threshold.
    %First row of Spec RNA has the threshold, every subsequent row is for a
    %cell in the image. Any value on (x+1,y) is the number of RNA spots
    %for cell x for threshold y(1).
    SpecRNA = SpecRNA(:,1:100);
    thres0 = zeros(2,size(SpecRNA,2));                              %thres0 will have the percent of cells without transcripts for each threshold. First row is thresholds and second row is the decimal fraction
    thres0(1,:) = SpecRNA(1,:);
    for j = 1:size(SpecRNA,2)
        j
        thres0(2,j) = numel(find(SpecRNA(2:size(SpecRNA,1),j)==0))/...
            (size(SpecRNA,1)-1);                                      %the number of cells that have 0 spots divided by the total for each threshold          
    end
    if ~isempty(thres0(1,find(thres0(2,:)>=prct,1,'first')))            %Make sure there is a threshold at which there are enough cells with 0 spots
        thress(counter19-1) = thres0(1,find(thres0(2,:)>=prct,1,'first'));  %If there is a threshold that results in enough off cells, use lowest such threshold
    else
        thress(counter19-1) = max(thres0(1,:));                             %If no threshold with enough off cells, use maximum threshold
    end
    counter19 = counter19+1
end
figure(20);clf;hist(thress,6); xlabel('Threshold');ylabel('Image Number');
axis([0 max(thress)+(min(thress)-0) 0 numel(find(thress == mode(thress)))+3])
title(strcat('Thresholds for ', num2str(prct*100), ' percent cells with 0 transcripts at 15 min timepoint with scrambled probes'))
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','black')
mean1 = mean(thress)
stdv1 = std(thress)
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',['mean = ' num2str(mean1) ', stdv = ' num2str(stdv1)],'FitBoxToText','on');

% Both together

thress = [0,0];                 %one row that has the lowest thresholds that result in the percentage of cells having 0 spots. Each value corresponds to one image
thres0 = zeros(2,1);            %first row is threshold, next rows are each image's percent of cells with transcripts
counter19 = 2;                  %counter for image number (thres0)
for i = [1:9]       %Task ids 1-3, timepoints 0,1,2
    i
    if i < 5
        SpecRNA = importdata(strcat('SpecRNA 20160206_0.4Mstep_WT_00min_SCR-TMR-3500_img',num2str(i)),'SpecRNA');            %SpecRNA seems to have a list of how many RNA spots were determined for each cell in each image for each threshold.
    else
        SpecRNA = importdata(strcat('SpecRNA 20160206_0.4Mstep_WT_15min_SCR-TMR-3500_img',num2str(i)),'SpecRNA');            %SpecRNA seems to have a list of how many RNA spots were determined for each cell in each image for each threshold.
    end
        %First row of Spec RNA has the threshold, every subsequent row is for a
    %cell in the image. Any value on (x+1,y) is the number of RNA spots
    %for cell x for threshold y(1).
    SpecRNA = SpecRNA(:,1:100);
    thres0 = zeros(2,size(SpecRNA,2));                              %thres0 will have the percent of cells without transcripts for each threshold. First row is thresholds and second row is the decimal fraction
    thres0(1,:) = SpecRNA(1,:);
    for j = 1:size(SpecRNA,2)
        j
        thres0(2,j) = numel(find(SpecRNA(2:size(SpecRNA,1),j)==0))/...
            (size(SpecRNA,1)-1);                                      %the number of cells that have 0 spots divided by the total for each threshold          
    end
    if ~isempty(thres0(1,find(thres0(2,:)>=prct,1,'first')))            %Make sure there is a threshold at which there are enough cells with 0 spots
        thress(counter19-1) = thres0(1,find(thres0(2,:)>=prct,1,'first'));  %If there is a threshold that results in enough off cells, use lowest such threshold
    else
        thress(counter19-1) = max(thres0(1,:));                             %If no threshold with enough off cells, use maximum threshold
    end
    counter19 = counter19+1
end
figure(21);clf;hist(thress,6); xlabel('Threshold');ylabel('Image Number');
axis([0 max(thress)+(min(thress)-0) 0 numel(find(thress == mode(thress)))+5])
title(strcat('Thresholds for ', num2str(prct*100), ' percent cells with 0 transcripts at both timepoints with scrambled probes'))
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','black')
mean1 = mean(thress)
stdv1 = std(thress)
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',['mean = ' num2str(mean1) ', stdv = ' num2str(stdv1)],'FitBoxToText','on');

%% Look at number of RNA spots per cell over time using cutoff

%CY5
cutoff = mean1+stdv1;                                                        %cutoff using mean and stdv
thresss = abs(thres0(1,:)-cutoff);                                           %find closest threshold to the cutoff
cutoff = thres0(1,find(thresss == min(thresss)));                             %new cutoff
counter19 = 1;
counter20 = 1;      %counter for timspot_temp (images of same timepoint to be averaged)
counter21 = 1;
timspot = timmins;
timspot_temp = [0,0];    %images of same timepoint to be averaged
for i = [1:2:131]       %All Task ids
    i;
    load(strcat('SpecRNA ',num2str(i)));
    SpecRNA = SpecRNA(:,1:100);
    thres_index = find(SpecRNA(1,:) == cutoff);
    spotnum = sum(SpecRNA(2:size(SpecRNA,1),thres_index));
    cellnum = size(SpecRNA,1)-1;
    timspot1(1,counter19) = spotnum/cellnum;
    if i > 2 & RNA_thresholds(2,i) < RNA_thresholds(2,i-2)
        counter20 = 1;
        timspot(2,counter21) = mean(timspot_temp);
        timspot(3,counter21) = std(timspot_temp);
        counter21 = counter21+1;
        timspot_temp = [0,0];
    elseif i > size(RNA_thresholds,2)-2
        timspot_temp(1,counter20) = spotnum/cellnum;
        timspot(2,counter21) = mean(timspot_temp);
        timspot(3,counter21) = std(timspot_temp);
    end
    timspot_temp(1,counter20) = spotnum/cellnum;
    counter20 = counter20+1;
    counter19 = counter19+1;
end
figure(17);clf;errorbar(timspot(1,:),timspot(2,:),timspot(3,:))
xlabel('Time After Input');ylabel('Spots per Cell');
title('CY5 (CTT1)')
%% Look at percent of on cells for each timepoint

%CY5

thresss = abs(thres0(1,:)-cutoff);                                           %find closest threshold
cutoff = thres0(1,find(thresss == min(thresss)));                             %new cutoff
counter19 = 1;
counter20 = 1;      %counter for timspot_temp (images of same timepoint to be averaged)
counter21 = 1;
timspot = timmins;
timspot_temp = [0,0];    %images of same timepoint to be averaged
for i = [1:2:131]       %All Task ids
    i;
    load(strcat('SpecRNA ',num2str(i)));
    SpecRNA = SpecRNA(:,1:100);
    cellon = numel(find(SpecRNA(2:size(SpecRNA,1),thres_index) > 0));
    celloff = numel(find(SpecRNA(2:size(SpecRNA,1),thres_index) == 0));
    timspot1(1,counter19) = cellon/(cellon+celloff);
    if i > 2 & RNA_thresholds(2,i) < RNA_thresholds(2,i-2)
        counter20 = 1;
        timspot(2,counter21) = mean(timspot_temp);
        timspot(3,counter21) = std(timspot_temp);
        counter21 = counter21+1;
        timspot_temp = [0,0];
    elseif i > size(RNA_thresholds,2)-2
        timspot_temp(1,counter20) = cellon/(cellon+celloff);
        timspot(2,counter21) = mean(timspot_temp);
        timspot(3,counter21) = std(timspot_temp);
    end
    timspot_temp(1,counter20) = cellon/(cellon+celloff);
    counter20 = counter20+1;
    counter19 = counter19+1;
end
figure(18);clf;errorbar(timspot(1,:),timspot(2,:),timspot(3,:))
xlabel('Time After Input');ylabel('Fraction of on Cells');
title('CY5 (CTT1)')

%% This is to generate a histogram depicting how high a threshold images needed to be in order to have a percent of their cells with 0 RNA spots

%TMR
for prct = [.95];                     %cutoff percent of cells with 0 spots
thress = [0,0];                 %one row that has the lowest thresholds that result in the percentage of cells having 0 spots
thres0 = zeros(2,1);            %first row is threshold, next rows are each image's percent of cells with transcripts
counter19 = 2;                  %counter for image number (thres0)
for i = [2:2:24,96:2:124]       %Task ids 1,2,3,13-16, timepoints 0,1,2,40,45,50,55,
    i
    load(strcat('SpecRNA ',num2str(i)));
    SpecRNA = SpecRNA(:,1:100);
    thres0 = zeros(2,size(SpecRNA,2));
    thres0(1,:) = SpecRNA(1,:);
    for j = 1:size(SpecRNA,2)
        j;
        thres0(2,j) = numel(find(SpecRNA(2:size(SpecRNA,1),j)==0))/...
            (size(SpecRNA,1)-1);                                      %the number of cells that have 0 spots divided by the total for each threshold          
    end
    if ~isempty(thres0(1,find(thres0(2,:)>=prct,1,'first')))
        thress(counter19-1) = thres0(1,find(thres0(2,:)>=prct,1,'first'));
    else
        thress(counter19-1) = max(thres0(1,:));
    end
    counter19 = counter19+1
end
figure(19);clf;hist(thress,6); xlabel('Threshold');ylabel('Image Number');
title(['TMR (STL1) Thresholds for ' num2str(prct*100) ' percent cells with 0 transcripts at off timepoints'])
axis([0 max(thress)+(min(thress)-0) 0 numel(find(thress == mode(thress)))+6])
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','black')
mean2 = mean(thress)
stdv2 = std(thress)
pause(1)
end


%% Look at number of RNA spots per cell over time using cutoff
%TMR
cutoff = mean2+stdv2;                                                        %cutoff using mean and stdv
thresss = abs(thres0(1,:)-cutoff);                                           %find closest threshold
cutoff = thres0(1,find(thresss == min(thresss)));                             %new cutoff
counter19 = 1;
counter20 = 1;      %counter for timspot_temp (images of same timepoint to be averaged)
counter21 = 1;
timspot = timmins;
timspot_temp = [0,0];    %images of same timepoint to be averaged
for i = [2:2:132]       %All Task ids
    i;
    load(strcat('SpecRNA ',num2str(i)));
    SpecRNA = SpecRNA(:,1:100);
    thres_index = find(SpecRNA(1,:) == cutoff);
    spotnum = sum(SpecRNA(2:size(SpecRNA,1),thres_index));
    cellnum = size(SpecRNA,1)-1;
    timspot1(1,counter19) = spotnum/cellnum;
    if i > 2 & RNA_thresholds(2,i) < RNA_thresholds(2,i-2)
        counter20 = 1;
        timspot(2,counter21) = mean(timspot_temp);
        timspot(3,counter21) = std(timspot_temp);
        counter21 = counter21+1;
        timspot_temp = [0,0];
    elseif i > size(RNA_thresholds,2)-2
        timspot_temp(1,counter20) = spotnum/cellnum;
        timspot(2,counter21) = mean(timspot_temp);
        timspot(3,counter21) = std(timspot_temp);
    end
    timspot_temp(1,counter20) = spotnum/cellnum;
    counter20 = counter20+1;
    counter19 = counter19+1;
end
figure(20);clf;errorbar(timspot(1,:),timspot(2,:),timspot(3,:))
xlabel('Time After Input');ylabel('Spots per Cell');
title('TMR (STL1)')

%% Look at percent of on cells for each timepoint
%TMR
thresss = abs(thres0(1,:)-cutoff);                                           %find closest threshold
cutoff = thres0(1,find(thresss == min(thresss)));                             %new cutoff
timspot = timmins;
counter19 = 1;
counter20 = 1;      %counter for timspot_temp (images of same timepoint to be averaged)
counter21 = 1;
timspot_temp = [0,0];    %images of same timepoint to be averaged
for i = [2:2:132]       %All Task ids
    i;
    load(strcat('SpecRNA ',num2str(i)));
    SpecRNA = SpecRNA(:,1:100);
    thres_index = find(SpecRNA(1,:) == cutoff);
    cellon = numel(find(SpecRNA(2:size(SpecRNA,1),thres_index) > 0));
    celloff = numel(find(SpecRNA(2:size(SpecRNA,1),thres_index) == 0));
    timspot1(1,counter19) = cellon/(cellon+celloff);
    if i > 2 & RNA_thresholds(2,i) < RNA_thresholds(2,i-2)
        counter20 = 1;
        timspot(2,counter21) = mean(timspot_temp);
        timspot(3,counter21) = std(timspot_temp);
        counter21 = counter21+1;
        timspot_temp = [0,0];
    elseif i > size(RNA_thresholds,2)-2
        timspot_temp(1,counter20) = cellon/(cellon+celloff);
        timspot(2,counter21) = mean(timspot_temp);
        timspot(3,counter21) = std(timspot_temp);
    end
    timspot_temp(1,counter20) = cellon/(cellon+celloff);
    counter20 = counter20+1;
    counter19 = counter19+1;
end
figure(21);clf;errorbar(timspot(1,:),timspot(2,:),timspot(3,:))
xlabel('Time After Input');ylabel('Fraction of on Cells');
title('TMR (STL1)')
