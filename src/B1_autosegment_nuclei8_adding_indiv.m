function [dapi_threshold,dapi_label,counter2,max_nuclei,DAPI_ims_added,dapi_label_low1] = B1_autosegment_nuclei8_adding(DAPI_ims,Yim,ManTh,counter2,max_nuclei,min_nucleus_size,max_nucleus_size,Ywin,dxy1);
%% This segments the cells by adding together binary images at every threshold, then removing pixels that are not present in most, and then going to each cell to segment it further.
%might need to change line 33 to "if false" in certain systems
%Line 150 is a kind of threshold that might need to be changed
%Yim = 1;
% min_nucleus_size = 20000;                                                     % minimum nuclear size
% max_nucleus_size = 10000;  %Ben Kesler 2/14/15 Changed from 20000 to 10000   % maximum nuclear size

a = size(DAPI_ims,1);                                                       % image size in number of pixels
z = size(DAPI_ims,3);                                                       % stack size in number of images

STD2D = NaN(1,z);                                                           %BK 5/16, needed this so timepoints with more stacks wouldn't be left over 
for i = 1:z;                                                                % Find the Image with the strongest DAPI signal (largest STD)
    STD2D(i) = std2(DAPI_ims(:,:,i));
end;
[p,ip] = max(STD2D);
range = 3                                                                   %use the 3 images around the max image
while ip + range > size(DAPI_ims,3) | ip - range < 1                        %reduce the range if this z stack is near the edges
    range = range - 1
end   
dapi_max = max(DAPI_ims(:,:,ip-range:ip+range),[],3);                               % maximum intensity projection in z-direction
if Yim == 1;
    figure(1); clf; imshow(dapi_max,[]); impixelinfo;  
    imwrite(uint16(dapi_max),['DAPI_max.tif']); 
else;
end;
Dapi_Mean = mean(dapi_max(:))
Dapi_Min = min(dapi_max(:))
Dapi_Max = max(dapi_max(:))
Dapi_Median = median(dapi_max(:))
if false;%Ywin
    threshold_sampling = 100
else
    threshold_sampling =200
end
%% find the nuclei-maximizing dapi threshold
if 10*Dapi_Median < Dapi_Max                                         %BK 4/27/2016
    dd = round((10*Dapi_Median-Dapi_Min)/threshold_sampling)
    dapi_threshold2 = Dapi_Min:dd:10*Dapi_Median;
else
    dd = round((Dapi_Max-Dapi_Min)/threshold_sampling)
    dapi_threshold2 = Dapi_Min:dd:Dapi_Max;
end
% for j = 1:size(dapi_threshold2,2)
%     dapi_bw = DAPI_ims > dapi_threshold2(j);                                        % only take DAPI intensities above the identified threshold for 3D stack
%     dapi_bw_max2(:,:,j) = max(dapi_bw,[],3);                                            % maxium z-direction
% end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % Ben Kesler 2/10/15 This goes through thresholds and determines how many
% % nuclei would result from each. The threshold at which the max number of
% % nuclei is obtained will be used
nuclei_num = zeros(1,size(dapi_threshold2,2));  %added 9/7 BK
DAPI_ims_added = zeros(size(DAPI_ims,1),size(DAPI_ims,2)); %added 4/29 BK
f = waitbar(1/size(dapi_threshold2,2),'Finding Nuclei')
for i = threshold_sampling/10:size(dapi_threshold2,2);                          %see how many nuclei for every threshold %
    waitbar(i/size(dapi_threshold2,2),f,'Finding Nuclei')
    i
    dapi_threshold = dapi_threshold2(i);
    dapi_bw = DAPI_ims > dapi_threshold;                                        % only take DAPI intensities above the identified threshold for 3D stack
    dapi_bw_max = max(dapi_bw(:,:,13:end),[],3);                                           % Maximum projection for the binary pixels above the threshold
    if Yim == 1;
        figure(2); clf; imshow(dapi_bw_max,[]); 
        imwrite(uint16(dapi_bw_max),['DAPI_THA_' num2str(i) '.tif']); 
    else;
    end;
    %% remove nuclei > than max_nucleus and < than min_nucleus
    dapi_normal = bwareaopen(dapi_bw_max, min_nucleus_size);                        % remove DAPI signal that are too small
    dapi_huge = bwareaopen(dapi_bw_max, max_nucleus_size);                          % remove DAPI signal that are too large
    dapi_bw2 = dapi_normal - dapi_huge;                                         % DAPI signal of the right size
    dapi_bw_max = max(dapi_bw2,[],3);                                           % Maximum projection for the binary pixels above the threshold
    dapi_OK = bwareaopen(dapi_bw_max,4);                                       % segment all DAPI spots
    dapi_label = bwlabeln(dapi_OK,8);                                          % label all segmented DAPI spots
    nuclei_num(i) = max(dapi_label(:));
    if Yim == 1;
        figure(3); clf; imshow(dapi_bw_max,[]); 
        imwrite(uint16(dapi_bw_max),['DAPI_THB_' num2str(i) '.tif']); 
    else;
    end;
    %added 4-29 BK
    DAPI_ims_added = DAPI_ims_added + dapi_OK;
    
end
close(f)
%% These are all code associated with testing and seeing what happens. Not needed to actual segmentation
% figure(60); clf; imshow(log(DAPI_ims_added),[])                                       % show what all thresholds added together looks like
%figure(7); clf; plot(dapi_threshold2,nuclei_num);                               % plot the number of nuclei depending on the threshold
%xlabel('Threshold');ylabel('Nuclei Number')
%filename = ['nuclei vs threshold samp ' num2str(threshold_sampling) ' task_id ' num2str(task_id) ' img ' num2str(counter) ' exp date ' exp_date '.fig' ]; savefig(filename);
%  figure(20); clf; imshow(DAPI_ims_final)
%  filename = ['Spot image ' num2str(threshold_sampling) ' task_id ' num2str(task_id) ' img ' num2str(counter) ' exp date ' exp_date '.fig' ]; savefig(filename);
% max_nuclei(1,counter2) = task_id;
% max_nuclei(2,counter2) = threshold_sampling;                                                     %record the sampling number
% max_nuclei(3,counter2) = dapi_threshold2(find(nuclei_num == max(nuclei_num),1,'last')); %threshold determined by max nuclei number
% max_nuclei(4,counter2) = nuclei_num(find(nuclei_num == max(nuclei_num),1,'last')); %find last max nuclei number
%counter2 = counter2+1

%% This is code for Manual thresholding based on choosing a threshold
%ik = 20;                              %set the threshold to the largest threshold that results in the maximum cells
% if ManTh == 1
% figure(7); clf; plot(dapi_threshold2,nuclei_num);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%  buu = 2
% while buu == 2;
% DAPI_ims_final = DAPI_ims_added;
% DAPI_ims_final(find(DAPI_ims_final < max(DAPI_ims_final(:))*ik/100)) = 0;
% dapi_normal = bwareaopen(DAPI_ims_final, min_nucleus_size);                        % remove DAPI signal that are too small
% dapi_huge = bwareaopen(DAPI_ims_final, max_nucleus_size);                          % remove DAPI signal that are too large
% DAPI_ims_final = dapi_normal - dapi_huge; 
% dapi_label = bwlabeln(DAPI_ims_final,8);
% xs = zeros(size(dapi_label(:)));
% ys = zeros(size(dapi_label(:)));
% for j = 1:max(dapi_label(:))
%     [x_es,y_es] = find(dapi_label == j);
%     xs(j,1) = mean(x_es);
%     ys(j,1) = mean(y_es);
% end
% figure(20); clf; imshow(DAPI_ims_final); hold on; 
% plot(ys,xs,'or','markersize',10);
% title([', cell num: ' num2str(max(dapi_label(:))) ' percent cutoff ' num2str(ik)]);
% 
% 
% [xi,yi,but] = ginput(1);
%     if but == 46; % >
% %        pause(1);
%         ik = ik + 1;
% %        figure(2); clf; imshow((uint16(dapi_bw_max2(:,:,ik))),[]); hold on; title([', th: ' num2str(dapi_threshold2(ik))]); % plot the image 500:1500,500:1500
% %        pause(1);                                                              %Ben Kesler 2/10/15 It was redundant to have these at
% %                                                                                 each if statement since it could instead be placed in the while loop
%     elseif but == 44; % <
% %        pause(1)
%         ik = ik - 1;
% %        figure(2); clf; imshow((uint16(dapi_bw_max2(:,:,ik))),[]); hold on; title([', th: ' num2str(dapi_threshold2(ik))]); % plot the image
% %        pause(1);
%     elseif but == 47; % >>
% %        pause(1);
%         ik = ik + 10;
% %        figure(2); clf; imshow((uint16(dapi_bw_max2(:,:,ik))),[]); hold on; title([', th: ' num2str(dapi_threshold2(ik))]); % plot the image 500:1500,500:1500
% %        pause(1);
%     elseif but == 109; % <<
% %        pause(1)
%         ik = ik - 10;
% %        figure(2); clf; imshow((uint16(dapi_bw_max2(:,:,ik))),[]); hold on; title([', th: ' num2str(dapi_threshold2(ik))]); % plot the image
% %        pause(1);
%     elseif but == 1;
%         buu = 1;
%     end;
% end;
% 
% end
% figure(20); clf; imshow(DAPI_ims_final); hold on; 
% plot(ys,xs,'or','markersize',10);
% title([', cell num: ' num2str(max(dapi_label(:))) ' percent cutoff ' num2str(ik)]);
% filename = ['Spot image ' num2str(threshold_sampling) ' task_id ' num2str(task_id) ' img ' num2str(counter) ' exp date ' exp_date '.fig' ]; savefig(filename);
% dapi_label = bwlabeln(DAPI_ims_final,8);
% max_nuclei(5,counter2) = max(dapi_label(:)); %find last max nuclei number
% max_nuclei(6,counter2) = ik;
% counter2 = counter2+1
%% Apply cutoff for added image
ik = 5                                                                     %This is a kind of cutoff. It's a percent of the maximum number of times a pixel is present 
DAPI_ims_final = DAPI_ims_added;
if ik > 0
    DAPI_ims_final(find(DAPI_ims_final < max(DAPI_ims_final(:))*ik/100)) = 0;
    dapi_normal = bwareaopen(DAPI_ims_final, min_nucleus_size);                        % remove DAPI signal that are too small
    dapi_huge = bwareaopen(DAPI_ims_final, max_nucleus_size);                          % remove DAPI signal that are too large
    DAPI_ims_cut = dapi_normal - dapi_huge;
    DAPI_ims_final = immultiply(DAPI_ims_cut,DAPI_ims_final);
end
dapi_label = bwlabeln(DAPI_ims_final,8);
xs = zeros(size(dapi_label(:)));
ys = zeros(size(dapi_label(:)));
for j = 1:max(dapi_label(:))                                                %collect data about each spot. Find the center (not adjusted by intensity right now)
    [x_es,y_es] = find(dapi_label == j);
    xs(j,1) = mean(x_es);
    ys(j,1) = mean(y_es);
end
% if Yim == 1
% figure(30); clf; imshow(DAPI_ims_added); hold on;
% filename = ['Added image ' num2str(threshold_sampling) ' task_id ' num2str(task_id) ' img ' num2str(counter) ' exp date ' exp_date '.fig' ]; savefig(filename);
% figure(20); clf; imshow(DAPI_ims_final); hold on; 
% end
% max_nuclei(5,counter2) = max(dapi_label(:));                                %the number of spots before investigating each one by one
% if Yim == 1
%     plot(ys,xs,'or','markersize',10);                                           %draw circles around every spot
%     title([', cell num: ' num2str(max(dapi_label(:))) ' percent cutoff ' num2str(ik)]);
%     filename = ['Spot image ' num2str(threshold_sampling) ' task_id ' num2str(task_id) ' img ' num2str(counter) ' exp date ' exp_date '.fig' ]; savefig(filename);
% end
%% Separate individual DAPI spots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dxy1 = 150;                                                                  % 7 ; define number of pixels from the max DAPI pixel'
dxy2 = dxy1*2.5;
m22 = max(dapi_label(:));                                                   % determine maximum number of DAPI stained nuclei
counter16 = 1;                                                               %counter for how many times the box had more than one nucleus
f = waitbar(1/m22,'Separating Joined Spots')
dapi_label_low1 = DAPI_ims_final;
for i = 1:m22;                                                              
    waitbar(i/m22,f,'Separating Joined Spots')
    i
    k1 = dapi_label == i;                                                  % Single cell per image
    k2a = immultiply(k1,DAPI_ims_final);
    m111 = max(k2a(:));
%     figure(20);clf;imshow(k1)
%     figure(21);clf;imshow(k2a,[])
%     numel(find(k1 > 0))
%     numel(find(k2a > 0))
%     numel(find(DAPI_ims_cut > 0))
%     numel(find(DAPI_ims_final > 0))
%     figure(20);clf;imshow(DAPI_ims_cut); title('DAPI ims cut')
%     figure(21);clf;imshow(DAPI_ims_final); title('DAPI ims final')
%     numel(find(dapi_label > 0))
%     numel(find(DAPI_ims_final > 0))
%     figure(20);clf;imshow(dapi_label); title('dapi label')
%     figure(21);clf;imshow(DAPI_ims_final); title('DAPI ims final')
    xA = round(xs(i));                                                      %use the center found earlier as the starting point for drawing the box
    yA = round(ys(i));
    A = xA-dxy2;                                                            % generate corner points with are dxy1 pixels away from the nuclear center
    B = xA+dxy2;
    C = yA-dxy2;
    D = yA+dxy2;
    if  A < 1;                                                              % if corner A is below 1 pixel or negative set equals 1
        A = 1;
    else;
    end;
    if B > a;                                                               % if corner B is above a pixels (1024 or 2048) set equals a
        B = a;
    else;
    end;
    if C < 1;                                                               % if corner C is below 1 pixel or negative set equals 1
        C = 1;
    else;
    end;
    if D > a;                                                                % if corner D is above a pixels (1024 or 2048) set equals a
        D = a;
    else;
    end;
    sqim = k2a(A:B,C:D);                                         % makes an image that consists only of the box around the nucleus
%     if Yim == 1;
%         figure(99); clf; imshow(sqim,[]);
%         figure(100); clf; imshow(dapi_normal,[]);
%     else;
%     end;
    sqim_lbl = bwlabeln(sqim,8);                                            %label the image
    thresmat = 1:max(k2a(:));                                                % The array that will have the cell numb for diff thresholds
    if true%max(sqim_lbl(:)) == 1;                                               %if there is more or less than one cell in the box, don't do the next steps
        for y = 1:max(k2a(:))
            sqim_temp = sqim > 0;
            sqim_temp(find(sqim < y)) = 0;
            dapi_normal = bwareaopen(sqim_temp, min_nucleus_size);                        % remove DAPI signal that are too small
            dapi_huge = bwareaopen(sqim_temp, max_nucleus_size);                          % remove DAPI signal that are too large
            sqim_temp = dapi_normal - dapi_huge; 
            sqim_lbl1 = bwlabeln(sqim_temp,8);
            thresmat(2,y) = max(sqim_lbl1(:));
%             k2a_temp = k2a;
%             k2a_temp(find(k2a < y)) = 0;
%             dapi_normal = bwareaopen(k2a_temp, min_nucleus_size);                        % remove DAPI signal that are too small
%             dapi_huge = bwareaopen(k2a_temp, max_nucleus_size);                          % remove DAPI signal that are too large
%             k2a_temp = dapi_normal - dapi_huge; 
%             k2a_lbl1 = bwlabeln(k2a_temp,8);
%             thresmat(2,y) = max(k2a_lbl1(:));
        end
        thres = thresmat(1,find(thresmat(2,:) == max(thresmat(2,:)),1,'last'));  %find the threshold that results in the most number of nuclei (maximum threhsold)
        thres2 = thresmat(1,find(thresmat(2,:) == max(thresmat(2,:)),1,'first'));  %find the threshold that results in the most number of nuclei (minimum threshold)      
        k2a(k2a < thres) = 0;                                        %Change the scaled image so it is at the new threshold
        if Yim
            sqim(sqim < thres) = 0;
            figure(101); clf; imshow(sqim,[])
        end
        k2a_bin = k2a > 0;                                                   %Make binary of new scaled image
        subtr_im = k1 - k2a_bin;                                            %Determine pixels that need to be subtracted from scaled image
        DAPI_ims_final(find(subtr_im > 0)) = 0;                                     %Change final image to subtract other pixels      
        k2a = immultiply(k1,dapi_label_low1);                                %Reset the image
        k2a(k2a < thres2) = 0;                                        %Change the scaled image so it is at the new threshold (for minimum stacks removed)
        if Yim
            sqim = k2a(A:B,C:D);                                         % makes an image that consists only of the box around the nucleus
            sqim(find(sqim < thres)) = 0;
            figure(101); clf; imshow(sqim,[])
        end
        k2a_bin = k2a > 0;                                                   %Make binary of new scaled image
        subtr_im = k1 - k2a_bin;                                            %Determine pixels that need to be subtracted from scaled image
        dapi_label_low1(find(subtr_im > 0)) = 0;                                     %Change final image to subtract other pixels
    else
        counter16 = counter16+1;
        thresmat(2,1) = 0;
    end
    if max(thresmat(2,:)) > 1
        thresmat;                                                           %If the semicolon is removed this will show the thresmat for whenever there was more than one nucleus
    end
end
close(f)
counter16
if Yim == 1;
    figure(5); clf; imshow(dapi_bw_max,[]); 
    figure(8); clf; imshow(log(DAPI_ims_final),[]); 
    figure(9); clf; imshow(log(dapi_label_low1),[]); 
    imwrite(uint16(DAPI_ims_final),['DAPI_FINAL1_.tif']); 
else;
end;

dapi_normal = bwareaopen(DAPI_ims_final, min_nucleus_size);                        % remove DAPI signal that are too small
dapi_huge = bwareaopen(DAPI_ims_final, max_nucleus_size);                          % remove DAPI signal that are too large
DAPI_ims_final = dapi_normal - dapi_huge; 
dapi_normal = bwareaopen(dapi_label_low1, min_nucleus_size);                        % remove DAPI signal that are too small
dapi_huge = bwareaopen(dapi_label_low1, max_nucleus_size);                          % remove DAPI signal that are too large
dapi_label_low1 = dapi_normal - dapi_huge; 
if Yim == 1;
    figure(6); clf; imshow(log(dapi_label_low1),[]); 
    figure(7); clf; imshow(dapi_bw_max,[]); 
    imwrite(uint16(DAPI_ims_final),['DAPI_FINAL2_.tif']); 
else;
end;
dapi_label = bwlabeln(DAPI_ims_final,8);
dapi_label_low1 = bwlabeln(dapi_label_low1,8);
xs = zeros(size(dapi_label(:)));
ys = zeros(size(dapi_label(:)));
for j = 1:max(dapi_label(:))
    [x_es,y_es] = find(dapi_label == j);
    xs(j,1) = mean(x_es);
    ys(j,1) = mean(y_es);
end
if Yim == 1;
    figure(7); clf; imshow(dapi_bw_max,[]); 
    imwrite(uint16(dapi_label),['DAPI_FINAL3_.tif']); 
else;
end;

% if Yim == 1
%     figure(22); clf; imshow(DAPI_ims_final); hold on;
%     plot(ys,xs,'or','markersize',10);
%     title([', cell num: ' num2str(max(dapi_label(:))) ' percent cutoff ' num2str(ik)]);
%     filename = ['Spot image with removal ' num2str(threshold_sampling) ' task_id ' num2str(task_id) ' img ' num2str(counter) ' exp date ' exp_date '.fig' ]; savefig(filename);
% end
% max_nuclei(6,counter2) = max(dapi_label(:)); %find last max nuclei number
% max_nuclei(7,counter2) = ik;
counter2 = counter2+1
% %% remove nuclei > than max_nucleus and < than min_nucleus
% dapi_normal = bwareaopen(dapi_bw, min_nucleus_size);                        % remove DAPI signal that are too small
% dapi_huge = bwareaopen(dapi_bw, max_nucleus_size);                          % remove DAPI signal that are too large
% dapi_bw2 = dapi_normal - dapi_huge;                                         % DAPI signal of the right size
% dapi_bw_max = max(dapi_bw2,[],3);                                           % Maximum projection for the binary pixels above the threshold
% if Yim == 1;
%     figure(4); clf; imshow(dapi_bw_max,[]);
% else;
% end;
% 
% %% Determine DAPI threshold for each individual cell
% % segment the DAPI signals in the image
% dapi_OK = bwareaopen(dapi_bw_max,4);                                       % segment all DAPI spots
% dapi_label = bwlabeln(dapi_OK,8);                                          % label all segmented DAPI spots
% m22 = max(dapi_label(:))                                                   % determine maximum nuber of DAPI stained nuclei
% if Yim == 1
% figure(5); clf; imshow(dapi_label,[]); impixelinfo;
% end
