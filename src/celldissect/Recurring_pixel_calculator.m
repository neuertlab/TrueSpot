%This function uses a blank image (images taken with no sample or slide on the microscope) to determine the recurring pixels in our
%images

[file,path] = uigetfile; %choose a blank image
  [stack, img_read] = tiffread2([path '/' file]);
            im_size = size(stack(1,1).data);
            testmat = zeros(7,1);
     counter_test = 1;
     for k = 1:size(stack,2)
            testmat(1,k) = k;
            slice_cut = k/size(stack,2);            %cutoff for how many slices spot is bright
            slice_check = img_read;            %number of slices from bottom to check
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
            if size(stack,2) > 1;
                for i = 1:slice_check
                    counter99 = counter99+1;
                    slice = imfilter(stack(i).data,w2);
%                     slice = CY5_ims(:,:,i);
                     cutoff = mean(slice(:))+3*std(double(slice(:)));
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
%              hist_pix_rand = imhistc(double(rand_hi(:)),double(im_size(1)*im_size(1)),0,double(im_size(1)*im_size(1)));      %histogram for number of times each position is randomly selected
            figure; 
            h = histogram(double(rand_hi(:)),double(im_size(1)*im_size(1)));
             hist_pix_rand = h.Values;     %histogram for number of times each position is randomly selected
             rand_recur_pix = find(hist_pix_rand >= slice_num*slice_cut);
             testmat(2,counter_test) = size(rand_recur_pix,2);
             [num2str(size(rand_recur_pix,2)) ' randomly recurring pixels' ]
% %             %%%
            hi_pixels_all(hi_pixels_all == 0) = [];         %remove all entries equal to 0
%             hist_pix = imhistc(double(hi_pixels_all(:)),double(im_size(1)*im_size(1)),0,double(im_size(1)*im_size(1)));     %make a histogram for number of times each position is determined as a high pixel
            h = histogram(double(hi_pixels_all(:)),double(im_size(1)*im_size(1)));     %make a histogram for number of times each position is determined as a high pixel
            hist_pix = h.Values;     %make a histogram for number of times each position is determined as a high pixel
                        close
            recurring_pixels = find(hist_pix >= slice_num*slice_cut);                   %Find all the positions that are above the cutoff for percent of times the pixel is bright
            [num2str(size(recurring_pixels,2)) ' recurring pixels' ]
             testmat(3,counter_test) = size(recurring_pixels,2);
            counter_test = counter_test+1;
            recurring_pixels = recurring_pixels';
            save('recurring pixels 3 out of 6','recurring_pixels')

          end
%       end
% %         end
%%%Look at empirically detemined cs random recurring pixels in a 2D plot
 testmat(4,:) = testmat(3,:)-testmat(2,:);
     testmat(6,:)  = testmat(3,:)./testmat(2,:);
    figure(15);clf; plot(testmat(1,:),testmat(3,:),'r','LineWidth',2); hold on
    plot(testmat(1,:),testmat(2,:),'b','LineWidth',2);
    plot(testmat(1,:),testmat(4,:),'g','LineWidth',2);
    legend({'Empirically Recurring Pixels','Randomly Recurring Pixels','Difference between Empirical and Random'})
    %xlabel('Percentile Cutoff')
    xlabel('Image Number Cutoff')
    xlim([0 i])
    ylim([0 max(testmat(3,:))])
    ylabel('Number of Pixels')
    %% Find minimum coefficient of variation
    for i = 2:size(testmat,2)-1
        testmat(7,i) = testmat(3,i)/std(testmat(3,i-1:i+1));
    end
    testmat(7,1) = testmat(7,2); testmat(7,size(testmat,2)) = testmat(7,size(testmat,2)-1); %Add ends
    yyaxis right
    plot(testmat(1,:),testmat(7,:),'m','LineWidth',2);
    ylabel('Inverse Coefficient of Variation')
    set(gca,'FontSize',14)
    legend({'Empirically Recurring Pixels','Randomly Recurring Pixels','Difference between Empirical and Random','Inverse Coefficient of Variation'})
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
%% Determine cutoff
%     cutoff1 = find(testmat(7,:) == max(testmat(7,:))) %Based on the maximum inverse coefficient of variation
     cutoff1 = find(testmat(2,:) <= 1,1) %Based on when randomly recurring pixels is 1 or 0
%% Recalculate recurring pixels
 im_size = size(stack(1,1).data);
     counter_test = 1;
     for k = cutoff1
            testmat(1,k) = k;
            slice_cut = k/size(stack,2);            %cutoff for how many slices spot is bright
            slice_check = img_read;            %number of slices from bottom to check
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
            if size(stack,2) > 1;
                for i = 1:slice_check
                    counter99 = counter99+1;
                    slice = imfilter(stack(i).data,w2);
%                     slice = CY5_ims(:,:,i);
                     cutoff = mean(slice(:))+3*std(double(slice(:)));
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
%              hist_pix_rand = imhistc(double(rand_hi(:)),double(im_size(1)*im_size(1)),0,double(im_size(1)*im_size(1)));      %histogram for number of times each position is randomly selected
            figure; 
            h = histogram(double(rand_hi(:)),double(im_size(1)*im_size(1)));
             hist_pix_rand = h.Values;     %histogram for number of times each position is randomly selected
             rand_recur_pix = find(hist_pix_rand >= slice_num*slice_cut);
             testmat(2,counter_test) = size(rand_recur_pix,2);
             [num2str(size(rand_recur_pix,2)) ' randomly recurring pixels' ]
% %             %%%
            hi_pixels_all(hi_pixels_all == 0) = [];         %remove all entries equal to 0
%             hist_pix = imhistc(double(hi_pixels_all(:)),double(im_size(1)*im_size(1)),0,double(im_size(1)*im_size(1)));     %make a histogram for number of times each position is determined as a high pixel
            h = histogram(double(hi_pixels_all(:)),double(im_size(1)*im_size(1)));     %make a histogram for number of times each position is determined as a high pixel
            hist_pix = h.Values;     %make a histogram for number of times each position is determined as a high pixel
                        close
            recurring_pixels = find(hist_pix >= slice_num*slice_cut);                   %Find all the positions that are above the cutoff for percent of times the pixel is bright
            [num2str(size(recurring_pixels,2)) ' recurring pixels' ]
             testmat(3,counter_test) = size(recurring_pixels,2);
            counter_test = counter_test+1;
            save('recurring pixels 3 out of 6','recurring_pixels')

          end
%       end
% %         end
%% Visualizing the recurring pixels
 figure(22); clf; imshow(stack(1).data,[median(stack(1).data(:))+std(double(stack(1).data(:))) cutoff]); hold on
[xs,ys] = ind2sub(size(stack(1).data),recurring_pixels); 
plot(ys,xs,'or','markersize',10)
%figure(); plot(recnum(1,2:12),recnum(2,2:12)); xlabel('Number of images'); ylabel('Number of Times Pixel is Present')   
