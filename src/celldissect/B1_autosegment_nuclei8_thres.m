function [dapi_threshold,dapi_label] = B1_autosegment_nuclei8_thres(DAPI_ims,Yim,ManTh,min_nucleus_size,max_nucleus_size);
%Yim = 1;
% min_nucleus_size = 200;                                                     % minimum nuclear size
% max_nucleus_size = 10000;  %Ben Kesler 2/14/15 Changed from 20000 to 10000   % maximum nuclear size

a = size(DAPI_ims,1);                                                       % image size in number of pixels
z = size(DAPI_ims,3);                                                       % stack size in number of images

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
    figure(1); clf; imshow(dapi_max,[]); impixelinfo;  imwrite(uint16(dapi_max),['dapi_max.tif']);                      
else;
end;
Dapi_Mean = mean(dapi_max(:))
Dapi_Min = min(dapi_max(:))
Dapi_Max = max(dapi_max(:))
Dapi_Median = median(dapi_max(:))
if ManTh
    threshold_sampling = 100
else
    threshold_sampling = 200
end
%% find the nuclei-maximizing dapi threshold
if 10*Dapi_Median < Dapi_Max 
    dd = round((10*Dapi_Median-Dapi_Min)/threshold_sampling)
    dapi_threshold2 = Dapi_Min:dd:10*Dapi_Median;
else
    dd = round((Dapi_Max-Dapi_Min)/threshold_sampling)
    dapi_threshold2 = Dapi_Min:dd:Dapi_Max;
end
for j = 1:size(dapi_threshold2,2)
    dapi_bw = DAPI_ims > dapi_threshold2(j);                                        % only take DAPI intensities above the identified threshold for 3D stack
    dapi_bw_max2(:,:,j) = max(dapi_bw,[],3);                                            % maxium z-direction
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % Ben Kesler 2/10/15 This goes through thresholds and determines how many
% % nuclei would result from each. The threshold at which the max number of
% % nuclei is obtained will be used
% Check subsets of the image and 
for i = threshold_sampling/10:size(dapi_threshold2,2);                          %see how many nuclei for every threshold
    i
    dapi_threshold = dapi_threshold2(i);
    dapi_bw = DAPI_ims > dapi_threshold;                                        % only take DAPI intensities above the identified threshold for 3D stack
    dapi_bw_max = max(dapi_bw,[],3);                                           % Maximum projection for the binary pixels above the threshold
    if Yim == 1;
        figure(3); clf; imshow(dapi_bw_max,[]); %imwrite(uint16(dapi_bw_max),['dapi_bw_max.tif']);      
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
        figure(4); clf; imshow(dapi_bw_max,[]); title(['th = ' num2str(dapi_threshold2(i))]); %imwrite(uint16(dapi_bw_max),['dapi_bw_max.tif']);      
    else;
    end;
end


ik = find(nuclei_num == max(nuclei_num),1,'last');                              %set the threshold to the largest threshold that results in the maximum cells
if ManTh == 1
figure(7); clf; plot(dapi_threshold2,nuclei_num);     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
buu = 2
while buu == 2;
%     dapi_threshold = ik;
%     dapi_bw = DAPI_ims > dapi_threshold;
%     %%%Ben Kesler 2/10/15 This removes the large and small spots for the
%     %%%image that the user looks at to better inform their decision
     dapi_threshold = dapi_threshold2(ik);                                      
     dapi_bw = DAPI_ims > dapi_threshold;                                        % only take DAPI intensities above the identified threshold for 3D stack
     dapi_bw_max = max(dapi_bw,[],3);                                           % Maximum projection for the binary pixels above the threshold
     dapi_normal = bwareaopen(dapi_bw_max, min_nucleus_size);                        % remove DAPI signal that are too small
     dapi_huge = bwareaopen(dapi_bw_max, max_nucleus_size);                          % remove DAPI signal that are too large
     dapi_bw2 = dapi_normal - dapi_huge;                                         % DAPI signal of the right size
%      dapi_bw_max3 = max(dapi_bw2,[],3);
%     figure(3); clf; imshow(dapi_bw_max); hold on; title([', th: ' num2str(dapi_threshold2(ik))]); %plot the image after removing the small and large nuclei
    %subplot(1,2,2); imshow(dapi_bw2);
    figure(4); clf; subplot(1,2,1); imshow(dapi_bw_max);  hold on; subplot(1,2,2); imshow(bwlabeln(dapi_bw2)); title([', th: ' num2str(dapi_threshold2(ik)) 'cell num ' num2str(nuclei_num(ik))] ); %plot the image after removing the small and large nuclei
    [xi,yi,but] = ginput(1);
    %%%Figure 2 is the version that includes the large and small spots
    %figure(2); clf; imshow((uint16(dapi_bw_max2(:,:,ik))),[]); hold on; title([', th: ' num2str(dapi_threshold2(ik))]); % plot the image 500:1500,500:1500
    %
    if but == 46; % >
%        pause(1);
        ik = ik + 1;
%        figure(2); clf; imshow((uint16(dapi_bw_max2(:,:,ik))),[]); hold on; title([', th: ' num2str(dapi_threshold2(ik))]); % plot the image 500:1500,500:1500
%        pause(1);                                                              %Ben Kesler 2/10/15 It was redundant to have these at
%                                                                                 each if statement since it could instead be placed in the while loop
    elseif but == 44; % <
%        pause(1)
        ik = ik - 1;
%        figure(2); clf; imshow((uint16(dapi_bw_max2(:,:,ik))),[]); hold on; title([', th: ' num2str(dapi_threshold2(ik))]); % plot the image
%        pause(1);
    elseif but == 47; % >>
%        pause(1);
        ik = ik + 10;
%        figure(2); clf; imshow((uint16(dapi_bw_max2(:,:,ik))),[]); hold on; title([', th: ' num2str(dapi_threshold2(ik))]); % plot the image 500:1500,500:1500
%        pause(1);
    elseif but == 109; % <<
%        pause(1)
        ik = ik - 10;
%        figure(2); clf; imshow((uint16(dapi_bw_max2(:,:,ik))),[]); hold on; title([', th: ' num2str(dapi_threshold2(ik))]); % plot the image
%        pause(1);
    elseif but == 1;
        buu = 1;
    end;
end;

end
dapi_threshold = dapi_threshold2(ik)

dapi_bw = DAPI_ims > dapi_threshold;                                        % only take DAPI intensities above the identified threshold for 3D stack
dapi_bw_max = max(dapi_bw,[],3);                                           % Maximum projection for the binary pixels above the threshold
if Yim == 1;
    figure(3); clf; imshow(dapi_bw_max,[]); imwrite(uint16(dapi_bw_max),['dapi_bw_max.tif']);    
else;
end;

%% remove nuclei > than max_nucleus and < than min_nucleus
dapi_normal = bwareaopen(dapi_bw_max, min_nucleus_size);                        % remove DAPI signal that are too small
dapi_huge = bwareaopen(dapi_bw_max, max_nucleus_size);                          % remove DAPI signal that are too large
dapi_bw2 = dapi_normal - dapi_huge;                                         % DAPI signal of the right size
dapi_bw_max = max(dapi_bw2,[],3);                                           % Maximum projection for the binary pixels above the threshold
if Yim == 1;
    figure(4); clf; imshow(dapi_bw_max,[]);  imwrite(uint16(dapi_bw_max),['dapi_bw_max2.tif']);    
else;
end;

%% Determine DAPI threshold for each individual cell
% segment the DAPI signals in the image
 dapi_OK = bwareaopen(dapi_bw_max,4);                                       % segment all DAPI spots
dapi_label = bwlabeln(dapi_OK,8);                                          % label all segmented DAPI spots
m22 = max(dapi_label(:));                                                   % determine maximum nuber of DAPI stained nuclei
if Yim == 1
figure(5); clf; imshow(dapi_label,[]); impixelinfo;  imwrite(uint16(dapi_label),['dapi_label.tif']);    
end
