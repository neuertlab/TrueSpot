function [thCY7Man,thAF700Man,thCY5Man,thAF594Man,thTMRMan,thYFPMan,thGFPMan] = AB_FindThreshold_TMR_AF594_CY5_B...
    (im_size,Im1,cells,trans_plane,S,CY7_ims,AF700_ims,CY5_ims,AF594_ims,TMR_ims,YFP_ims,GFP_ims,thA,ths,Ych,max_int_thres,max_int_spot_th);
%% This code will find recurring pixels (dead pixels) that should be filtered out

%if Im1 == 1 & S == 1
     testmat = zeros(3,1);
     counter_test = 1;
%    for zeb =  [60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,99.9]
%    for seb = [2,3,4,5,6,7,8] 
%          for jeb = [1/4,3/8,1/2,5/8,3/4,7/8,1]
%            percentile1 = zeb;
%      slice_cut = jeb;            %cutoff for how many slices spot is bright
%            testmat(1,counter_test) = jeb;
%             testmat(5,counter_test) = seb;
            slice_cut = 4/8;            %cutoff for how many slices spot is bright
            slice_check = 6;            %number of slices from bottom to check
            counter99 = 0;
            hi_pixels_all = zeros(1,1);             %This will have all the highest intensity pixels for each stack
            'determining recurring pixels'
            %%% First Strategy Using a Percentile cutoff to determine bright spots, and
%%% Strategy using a percentile cutoff for determining bright spots seeing if consistent images
%             percentile1 = 99.99;         %percentile for bright spots
%             if size(CY5_ims,1) > 1;
%                 for i = 1:slice_check
%                     counter99 = counter99+1;
%                     slice = CY5_ims(:,:,i);
%                     cutoff = prctile(slice(:),percentile1);
%                     temp_pix = find(slice > cutoff);
% %                     figure(23); clf; imshow(slice,[]); hold on
% %                     [xs,ys] = ind2sub(size(CY5_ims(:,:,1)),temp_pix); 
% %                     plot(ys,xs,'or','markersize',10)
%                     hi_pixels_all(counter99,1:size(temp_pix,1)) = temp_pix;
%                 end
%             end
%             if size(AF594_ims,1) > 1;
%                 for i = 1:slice_check
%                     counter99 = counter99+1;
%                     slice = AF594_ims(:,:,i);
%                     cutoff = prctile(slice(:),percentile1);
%                     temp_pix = find(slice > cutoff);
%                     hi_pixels_all(counter99,1:size(temp_pix,1)) = temp_pix;
%                 end
%             end
%             if size(TMR_ims,1) > 1;
%                 for i = 1:slice_check
%                     counter99 = counter99+1;
%                     slice = TMR_ims(:,:,i);
%                     cutoff = prctile(slice(:),percentile1);
%                     temp_pix = find(slice > cutoff);
%                     hi_pixels_all(counter99,1:size(temp_pix,1)) = temp_pix;
%                 end
%             end
%             slice_num = size(hi_pixels_all,1);
%%% Second Strategy using a filter
w2 = [-1 -1 -1;...   % edge detection filter
      -1 +8 -1;...
      -1 -1 -1];
% w2 = [-1 -1 -1;...   % edge detection filter (has to be 4x mean intensity of surrounding pixels)
%       -1 +2 -1;...
%       -1 -1 -1];
            if size(CY5_ims,1) > 1;
                for i = 1:slice_check
                    counter99 = counter99+1;
                    slice = imfilter(CY5_ims(:,:,i),w2);
                     cutoff = mean(slice(:))+3*std(slice(:))
%                     cutoff = 1;
                    temp_pix = find(slice > cutoff);
%                     temp_pix_not = slice(:);
%                     temp_pix_not(temp_pix) = []; 
%                     slice2 = CY5_ims(:,:,i);
%                     figure(33); clf; imshow(slice,[]); hold on
%                     [xs,ys] = ind2sub(size(CY5_ims(:,:,1)),temp_pix); 
%                      plot(ys,xs,'or','markersize',10)
                    hi_pixels_all(counter99,1:size(temp_pix,1)) = temp_pix;
                end
            end
%             if size(AF594_ims,1) > 1;
%                 for i = 1:slice_check
%                     counter99 = counter99+1;
%                     slice = imfilter(AF594_ims(:,:,i),w2);
%                     cutoff = mean(slice(:))+10*std(slice(:));
%                     temp_pix = find(slice > cutoff);
% %                     figure(23); clf; imshow(slice,[]); hold on
% %                     [xs,ys] = ind2sub(size(CY5_ims(:,:,1)),temp_pix); 
% %                      plot(ys,xs,'or','markersize',10)
%                     hi_pixels_all(counter99,1:size(temp_pix,1)) = temp_pix;
%                 end
%             end
%             if size(TMR_ims,1) > 1;
%                 for i = 1:slice_check
%                     counter99 = counter99+1;
%                     slice = imfilter(TMR_ims(:,:,i),w2);
%                     cutoff = mean(slice(:))+10*std(slice(:));
%                     temp_pix = find(slice > cutoff);
%                     figure(23); clf; imshow(slice,[]); hold on
%                     [xs,ys] = ind2sub(size(CY5_ims(:,:,1)),temp_pix); 
%                      plot(ys,xs,'or','markersize',10)
%                     hi_pixels_all(counter99,1:size(temp_pix,1)) = temp_pix;
%                 end
%             end
            %%%
            slice_num = size(hi_pixels_all,1);
            %%%for graphing number of recurring pixels depending on cutoff  for slices
            %     total = size(hi_pixels_all,1)*size(hi_pixels_all,2);                        %go through every pixel
            %     counter20 = 0;              %counter for recurring pixels
            %     recurring_pixels = [0];
            %     recnum(1,1:slice_num) = 1:slice_num;                                      %a matrix. the first row specifies how many images the pixels recur on, and the second row has the number of pixels
            %     for j = 1:slice_num
            %         counter20 = 0;              %counter for recurring pixels
            %         recurring_pixels = [0];
            %         for i = 1:total
            %             if hi_pixels_all(i) ~= 0 & ...
            %                     size(find(hi_pixels_all == hi_pixels_all(i)),1) >= j & ...
            %                     size(find(recurring_pixels == hi_pixels_all(i)),1) == 0
            %                 counter20 = counter20+1;
            %                 recurring_pixels(counter20,1) = hi_pixels_all(i);
            %             end
            %         end
            %         if j == slice_num*slice_cut
            %             save('recurring pixels 2 out of 4','recurring_pixels')
            %         end
            %         recnum(2,j) = counter20;
            %     end
            %%%
%             %%% Determine how many recurring pixels would occur at random
            rand_hi= zeros(size(hi_pixels_all));
            for j = 1:slice_num
                temp_hi =hi_pixels_all(j,:);                                            %extract current slice
                temp_hi1 = temp_hi;
                temp_hi1(temp_hi1 ==0) = [];
                rand_hi_temp = randsample(im_size(1)*im_size(1),size(temp_hi1,2)); %generates random positions equal to toal number of positions in hi_pixels_all
                rand_hi_temp = rand_hi_temp';
                rand_hi(j,1:size(rand_hi_temp,2))=rand_hi_temp; %generates random positions equal to toal number of positions in hi_pixels_all
             end
             hist_pix_rand = imhistc(double(rand_hi(:)),double(im_size(1)*im_size(1)),0,double(im_size(1)*im_size(1)));      %histogram for number of times each position is randomly selected
             rand_recur_pix = find(hist_pix_rand >= slice_num*slice_cut);
             testmat(2,counter_test) = size(rand_recur_pix,1);
             [num2str(size(rand_recur_pix,1)) ' randomly recurring pixels' ]
% %             %%%
            hi_pixels_all(hi_pixels_all == 0) = [];         %remove all entries equal to 0
            hist_pix = imhistc(double(hi_pixels_all(:)),double(im_size(1)*im_size(1)),0,double(im_size(1)*im_size(1)));     %make a histogram for number of times each position is determined as a high pixel
            recurring_pixels = find(hist_pix >= slice_num*slice_cut);                   %Find all the positions that are above the cutoff for percent of times the pixel is bright
            [num2str(size(recurring_pixels,1)) ' recurring pixels' ]
             testmat(3,counter_test) = size(recurring_pixels,1);
            counter_test = counter_test+1;
            save('recurring pixels 3 out of 6','recurring_pixels')

%          end
%       end
% %         end
%%%Look at empirically detemined cs random recurring pixels in a 2D plot
 testmat(4,:) = testmat(3,:)-testmat(2,:);
     testmat(6,:)  = testmat(3,:)./testmat(2,:);
%     figure(15);clf; plot(testmat(1,:),testmat(3,:),'r'); hold on
%     plot(testmat(1,:),testmat(2,:),'b');
%     plot(testmat(1,:),testmat(4,:),'g');
%     legend({'Empirically Recurring Pixels','Randomly Recurring Pixels','Difference between Empirical and Random'})
%     %xlabel('Percentile Cutoff')
%     xlabel('Image Number Cutoff')
%     xlim([0 1])
%     ylabel('Number of Pixels')
%%% Looking at optimization with multiple values for both percentile
%      %%% and image cutoff
%     %Need to run lines defining zeb and jeb (entire vector)
%     index1 = find(testmat(4,:)==max(testmat(4,:)));
%     [num2str(testmat(1,index1)) ' optimal proportion images cutoff']
% %    [num2str(testmat(5,index1)) ' optimal percentile pixels cutoff']
%     [num2str(testmat(5,index1)) ' optimal standard deviation multiplier for threshold']
%     [num2str(testmat(3,index1)) ' empirical recurring pixels']
%     [num2str(testmat(2,index1)) ' random recurring pixels']
%     [num2str(testmat(4,index1)) ' more empirical recurring pixels than random']
%   %difference between empirical and random spots  
%     zmesh1 = zeros(size(seb,2),size(jeb,2));
%     for j = 1:size(testmat,2)
%         zmesh1(find(seb == testmat(5,j)),find(jeb == testmat(1,j)))= testmat(4,j);
%     end
%     figure(17); mesh(jeb,seb,zmesh1);
%     xlabel('Proportion of Images to be Determined Recurring')
% %    ylabel('Percentile Cutoff to be Determined Bright')
%     ylabel('Standard Deviation Multiplier for Threshold')
%     zlabel('Difference between Empirically Determined Bright spots and Random') 
%     % ratio if empirical and random
%         zmesh2 = zeros(size(seb,2),size(jeb,2));
%     for j = 1:size(testmat,2)
%         zmesh2(find(seb == testmat(5,j)),find(jeb == testmat(1,j)))= testmat(6,j);
%     end
%     figure(18); mesh(jeb,seb,zmesh2);
%     xlabel('Proportion of Images to be Determined Recurring')
%     ylabel('Standard Deviation Multiplier for Threshold')
% %    ylabel('Percentile Cutoff to be Determined Bright')
%     zlabel('Ratio between Empirically Determined Bright spots and Random') 
% %end
%     %empirical spots determined
%     zmesh2 = zeros(size(seb,2),size(jeb,2));
%     for j = 1:size(testmat,2)
%         zmesh2(find(seb == testmat(5,j)),find(jeb == testmat(1,j)))= testmat(3,j);
%     end
%     figure(19); mesh(jeb,seb,zmesh2);
%     xlabel('Proportion of Images to be Determined Recurring')
% %    ylabel('Percentile Cutoff to be Determined Bright')
%     ylabel('Standard Deviation Multiplier for Threshold')
%     zlabel('Empirically Determined Bright spots') 

% figure(22); clf; imshow(CY5_ims(:,:,1),[]); hold on
% [xs,ys] = ind2sub(size(CY5_ims(:,:,1)),recurring_pixels); 
% plot(ys,xs,'or','markersize',10)
%figure(); plot(recnum(1,2:12),recnum(2,2:12)); xlabel('Number of images'); ylabel('Number of Times Pixel is Present')   
h =sprintf('%03d',Im1);  % defines label
%% Shows Trans image
% figure(10000); clf; imshow(cells,[]); 
%     title(['S: ' num2str(S) ', Im: ' h ' Click the upper left-hand corner of the viewing window']); impixelinfo;  hold on; 
%% Show Fluorescence images
        CY51 = max(CY5_ims,[],3);
         CY51 = CY51/max(CY51(:));
         LminF = min(CY51(:)); %median(CY51(:))-round(0*std(CY51(:))) % min(min(TMR3Dfilter(:,:,I)))
        LmaxF = 3*(median(CY51(:))+round(10*std(CY51(:)))) %75  % max(max(TMR3Dfilter(:,:,I)))/1
            CY51 = imadjust(CY51,[LminF LmaxF]);
        AF5941 = max(AF594_ims,[],3);
         AF5941 = AF5941/max(AF5941(:));
       LminF = min(AF5941(:)); %median(AF5941(:))-round(0*std(AF5941(:))) % min(min(TMR3Dfilter(:,:,I)))
        LmaxF = 3*(median(AF5941(:))+round(10*std(AF5941(:)))) %75  % max(max(TMR3Dfilter(:,:,I)))/1
            AF5941 = imadjust(AF5941,[LminF LmaxF]);   
        TMR1 = max(TMR_ims,[],3);
         TMR1 = TMR1/max(TMR1(:));
        LminF = min(TMR1(:)); %median(TMR1(:))-round(0*std(TMR1(:))) % min(min(TMR3Dfilter(:,:,I)))
        LmaxF = 3*(median(TMR1(:))+round(10*std(TMR1(:)))) %75  % max(max(TMR3Dfilter(:,:,I)))/1
            TMR1 = imadjust(TMR1,[LminF  LmaxF ]);
            if size(CY51) == im_size; R = CY51; else; R = zeros(im_size);  end
            if size(AF5941) == im_size; B = AF5941; else; B = zeros(im_size);end
            if size(TMR1) == im_size; G = TMR1; else; G = zeros(im_size); end
            RGB = cat(3,R,G,B); 
            figure(10000); clf; imshow(RGB,[]); 
    title(['S: ' num2str(S) ', Im: ' h ' Click the upper left-hand corner of the viewing window']); impixelinfo;  hold on; 
clear CY51 AF5941 TMR1
    %% Shows Segmentation File
% mm = max(cells(:)); % maximum number of cells
% figure(100); clf;subplot(1,2,1); imshow(trans_plane,[]); subplot(1,2,2); imshow(cells,[0 1]); ...
%     title(['S: ' num2str(S) ', Im: ' h ', Cells: ' num2str(mm)]); impixelinfo;  hold on; 
%%

[xi,yi,but] = ginput(1); % reads the mouse curser position
if but == 1;
    plot(xi,yi,'ro'); % plots an red circle around the choosen spot
    xy(2) = round(xi); % get the x and y coordinates
    xy(1) = round(yi);

    %% filter CY7 images
    if false; %Ych(1) == 1;
        T = 'filter CY7 images'
        [thCY7Man] = AB2_Threshold3Dim(CY7_ims,xy,thA,T,max_int_thres,max_int_spot_th);
    else;
        thCY7Man = NaN;
    end;

    %% filter AF700 images
    if false; %Ych(2) == 1;
        T = 'filter AF700 images'
        [thAF700Man] = AB2_Threshold3Dim(AF700_ims,xy,thA,T,max_int_thres,max_int_spot_th);
    else;
        thAF700Man = NaN;
    end;


    %% filter CY5 images
    if Ych(4);%Ych(3) == 1;
        T = 'filter CY5 images'
        Ych2 = 3;
        [thCY5Man] = AB2_Threshold3Dim(CY5_ims,xy,thA,T,max_int_thres,max_int_spot_th);
    else;
        thCY5Man = NaN;
    end;

    %% filter AF594 images
    if Ych(5);%Ych(4) == 1;
        T = 'filter AF594 images'
        [thAF594Man] = AB2_Threshold3Dim(AF594_ims,xy,thA,T,max_int_thres,max_int_spot_th);
    else;
        thAF594Man = NaN;
    end;

    %% filter TMR images
    if Ych(6);%Ych(5) == 1;
        T = 'filter TMR images'
        Ych2 = 5;
        [thTMRMan] = AB2_Threshold3Dim(TMR_ims,xy,thA,T,max_int_thres,max_int_spot_th);
    else;
        thTMRMan = NaN;
    end;

    %% filter YFP images
    if Ych(6) == 1;
        T = 'filter YFP images'
        [thYFPMan] = AB2_Threshold3Dim(YFP_ims,xy,thA,T,max_int_thres,max_int_spot_th);
    else;
        thYFPMan = NaN;
    end;

    %% filter GFP images
    if Ych(7) == 1;
        T = 'filter GFP images'
        [thGFPMan] = AB2_Threshold3Dim(GFP_ims,xy,thA,T,max_int_thres,max_int_spot_th);
    else;
        thGFPMan = NaN;
    end;
else;
end;