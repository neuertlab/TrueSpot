function [TMRmax,TMRmaxF,thall,TMR3D3immax,TMRmaxFimmax,TMR3D,TMR3Dfilter,TMR3Dback]...
    = B2_TMRprocess3Dim(TMR_ims,th1,f); % TMR3D,TMR3Dfilter,TMR3Dback

% f =4;
% TMR_ims = CY5_ims;
%th1 = thTMR;
%th1 = thCY5;
%% filter sets up for gaussian filtering
mu1 = 7;
mu2 = 7;
s1 = 2; % 2s1 = FWHM of a diffraction limited RNA spot at TMR fitted with a gaussian
fs = 13;
[X1,Y1] = meshgrid(1:1:fs);
w1 = 10*exp(-((X1-mu1).^2)./(2*s1^2)-((Y1-mu2).^2)./(2*s1^2));
w1 = w1./sum(w1(:));

w2 = [-1 -1 -1;...
      -1 +8 -1;...
      -1 -1 -1];
  
%% load TMR images
%filepath = 'C:\Daten\Code\GN_FISH_Yeast_Image_Processing\data\';
%strTMR = char(strcat('TMR010'));
% fileTMR = [filepath strTMR];
% tmr1=tiffread1(fileTMR); % read raw data and take the pixels from the active cell label
% TMR3D = NaN(A,A,Z);
% for i = 1:Z;
%     tmr = tmr1(i).data;
%     TMR3D(:,:,i) = uint16(tmr);
% end;
TMR3D = uint16(TMR_ims);
TMRmax = max(TMR3D,[],3); %% Max projection original 3D stack
A = size(TMR3D,1);
Z = size(TMR3D,3);
%% Make recurring bright pixels an average of the pixel intensities surrounding it (added from AB2_Threshold3DIM_BK_auto on 1/10/2017)
load('recurring pixels 3 out of 6')
exmat(1:size(TMR3D,1),1) = 1:size(TMR3D,1);                                 %exmat is a matrix that contains the border indices, which will not be subject filtered for recurring pixels                                                  
exmat(1:size(TMR3D,1),4) = (size(TMR3D,1)^2-(size(TMR3D,1)-1)):size(TMR3D,1)^2;
exmat(1:2,2:3) = 1;
exmat(3:size(TMR3D,1),2) = size(TMR3D,1)*(2:(size(TMR3D,1)-1));
exmat(3:size(TMR3D,1),3) = size(TMR3D,1)*(1:(size(TMR3D,1)-2))+1;
for kl=1:Z;
    TMR = TMR3D(:,:,kl);
    for j = 1:size(recurring_pixels,1)          %go through each documented recurring pixel
        a1 = recurring_pixels(j);               %load the recurring pixel                  
        a2 = exmat == a1;                       %check if a1 is a border index
        if max(a2(:)) == 0;                                             %excludes pixels on the border
            surr_pix = [a1-1,a1+1,a1-2048, a1-2049,a1-2047,a1+2048,a1+2047,a1+2049]; %surrounding pixels CHANGE IF IMAGE RESOLUTION CHANGES
            TMR(a1) = mean(TMR(surr_pix));                                    %replace pixel with mean of surrounding pixels
        end
    end
    TMR3D(:,:,kl) = TMR;
end
%figure(55);imshow(TMRmax,[102 250])
%% Filtering images with a laplacian of a gaussian
TMR3Dfilter = uint16(NaN(A,A,Z));
p = 3;
con = 8;
for kk=1:Z;
    for nn = 1:f;
        kk;
        TMR = TMR3D(:,:,kk);
        TMR1 = imfilter(TMR, w1); % gauss smoothing
        TMR2 = imfilter(TMR1,w2); % spot filtering
        TMR3 = TMR2 > th1; % remove pixels below th1
        TMR4 = bwareaopen(TMR3,p,con); % remove objects thare larger as p - pixels and count connected regions with con pixels as one object
        TMR5 = immultiply(TMR4,TMR); 
        thall(kk,1) = std2(TMR2(:)); % determine threshold in each plane from the STD of background pixels;
        thall(kk,2) = mean2(TMR2(:)); % determine threshold in each plane from the STD of background pixels;
        thall(kk,3) = var(double(TMR2(:))); % determine threshold in each plane from the STD of background pixels;
        thall(kk,4) = median(double(TMR2(:))); % determine threshold in each plane from the STD of background pixels;
%        TMR6 = imfilter(TMR5, w1); % gauss smoothing
        TMR3Dfilter(:,:,kk) = TMR5;
%         if (kk == 8) % visulization
%             figure(1);clf;subplot(1,2,1);imshow(TMR,[150 300]);subplot(1,2,2);imshow(TMR1,[150 300]);impixelinfo; % CY5: [150 300]
%             figure(2);clf;subplot(1,2,1);imshow(TMR,[150 300]);subplot(1,2,2);imshow(TMR2,[0 50]);impixelinfo;    % CY5: [0 50]
%             figure(3);clf;subplot(1,2,1);imshow(TMR,[150 300]);subplot(1,2,2);imshow(TMR5,[0 300]);impixelinfo;    % CY5: [0 50]
%         else
%         end;
        
    end;
end;
% %th1 = round(1.5*std2(double(TMR3Dfilter(:)))) % determine threshold in each plance from the STD of background pixels;
% %th1 = max(thall); % determine threshold in each plance from the STD of background pixels;
% for kk=1:Z;
%     TMR = TMR3Dfilter(:,:,kk);
%     TMR3 = TMR > th1; % remove pixels below th1
%     TMR4 = bwareaopen(TMR3,p,con); % remove objects thare larger as p - pixels and count connected regions with con pixels as one object
%     TMR5 = immultiply(TMR4,TMR); 
%     if (kk == 8) % visulization
%         figure(3);clf;subplot(1,2,1);imshow(TMR,[150 300]);subplot(1,2,2);imshow(TMR5,[0 20]);
%         %figure;subplot(1,2,1);imshow(log(double(TMR)),[]);subplot(1,2,2);imshow(TMR5,[0 20]);
%     else
%     end;
%     for k = 1:3; % repead filtering three more times
%         TMR1 = imfilter(TMR5,w1);
%         TMR2 = imfilter(TMR1,w2); 
%         TMR3 = TMR2 > th1; % remove pixels below th1
%         TMR4 = bwareaopen(TMR3,p,con); 
%         TMR5 = immultiply(TMR4,TMR2); 
%         if (kk == 8); % visulization
%             figure(3+k);subplot(1,2,1);imshow(TMR,[150 300]);subplot(1,2,2);imshow(TMR5,[0 20]);
%             %figure;subplot(1,2,1);imshow(log(double(TMR)),[]);subplot(1,2,2);imshow(TMR5,[0 20]);
%         else
%         end;
%     end
%     TMR3Dfilter(:,:,kk) = TMR5;
% end;
TMR3Dfilter(1:mu1*2,:,:) = 0;
TMR3Dfilter(A-mu1*2:A,:,:) = 0;
TMR3Dfilter(:,1:mu1*2,:) = 0;
TMR3Dfilter(:,A-mu1*2:A,:) = 0;
TMRmaxF = max(TMR3Dfilter,[],3);                %% Max projection filtered 3D stack

%% Determine brightest pixel in 3D for each mRNA spot
TMR3D3bw = imregionalmax(TMR3Dfilter,26);       % 6 18 26        
TMR3D3immax = immultiply(TMR3Dfilter, TMR3D3bw);     %% Filtered 3D stack
TMRmaxFimmax = max(TMR3D3immax,[],3);                %% Max projection filtered 3D stack
% IM3 = TMR3D3immax;
% try clear xall yall; end;
% [yall,xall] = find(IM3 > th1);
%%%%%%%%%%%%%%%%
% %% START here
% for i = 1:23;
%     clear x1 y1;
%     [y1,x1] = find(IM3(:,:,i) > th1);        % find all pixel coordianates for pixels above previosly selected threshold 
%     figure(1112+i); clf; 
%     imshow(TMR3D(:,:,i),[150 300]); hold on; plot(x1,y1,'or','markersize',10); impixelinfo; % TMR: 1000 3000
% end;
%%%%%

% %% Deterime background for each image
%TMR3Dback = NaN;
TMR3Dback = uint16(NaN(A,A,Z));
for kk=1:Z;
    kk;
    TMR = TMR3D(:,:,kk);
    TMR3Dback(:,:,kk) = medfilt2(TMR, [20 20]); % gauss smoothing
%     TMR2 = imsubtract(TMR,TMR1);
%     figure(kk); clf; 
%     subplot(1,2,1); imshow(TMR1,[2000 6000]); hold on;
%     subplot(1,2,2); imshow(TMR2,[0 5000]); hold on;
end;


% %% Deterime x-y positions for initial fits
% clear xall yall
% xall = NaN(10000,size(TMR3D,3));
% yall = NaN(10000,size(TMR3D,3));
% for kk = 1:Z;
%     [y1,x1] = find(TMR3D3immax(:,:,kk) > th1);
%     xall(1:size(x1,1),kk) = x1;
%     yall(1:size(y1,1),kk) = y1;
%     x2 = find(x1 == 2);
%     y2 = find(y1 == 2);
%     xall(x2,kk) = NaN;
%     yall(x2,kk) = NaN;
%     xall(y2,kk) = NaN;
%     yall(y2,kk) = NaN;
%     x2 = find(x1 == 1023);
%     y2 = find(y1 == 1023);
%     xall(x2,kk) = NaN;
%     yall(x2,kk) = NaN;
%     xall(y2,kk) = NaN;
%     yall(y2,kk) = NaN;
% end;