function [th] = AB2_Threshold3Dim(TMR_ims,xy,thA,T,max_int_thres,max_int_spot_th)
TMR3D = TMR_ims;

%% filter sets up for gaussian filtering based on experimental PSF
mu1 = 7;
mu2 = 7;
s1 = 2; % 2s1 = FWHM of a diffraction limited RNA spot at TMR fitted with a gaussian
fs = 13;
[X1,Y1] = meshgrid(1:1:fs);
w1 = 10*exp(-((X1-mu1).^2)./(2*s1^2)-((Y1-mu2).^2)./(2*s1^2));
w1 = w1./sum(w1(:));
%figure(1010); clf; subplot(1,2,1); imshow(g,[]);impixelinfo; subplot(1,2,2); plot(g(mu1,:)); ;


w2 = [-1 -1 -1;...   % edge detection filter
      -1 +8 -1;...
      -1 -1 -1];

A = size(TMR3D,1);
Z = size(TMR3D,3);

TMR3D = uint16(TMR3D); % convert 3D image sections into 16bit format to save memory
TMRmax = max(TMR3D,[],3); %% Max projection original 3D stack
% figure; imshow(TMRmax,[])

%% Make recurring bright pixels an average of the pixel intensities surrounding it (added from AB2_Threshold3DIM_BK_auto on 1/10/2017)
load('recurring pixels 3 out of 6')
exmat(1:size(TMR3D,1),1) = 1:size(TMR3D,1);                                 %exmat is a matrix that contains the border indices, which will not be subject filtered for recurring pixels                                                  
exmat(1:size(TMR3D,1),4) = (size(TMR3D,1)^2-(size(TMR3D,1)-1)):size(TMR3D,1)^2;
exmat(1:2,2:3) = 1;
exmat(3:size(TMR3D,1),2) = size(TMR3D,1)*(2:(size(TMR3D,1)-1));
exmat(3:size(TMR3D,1),3) = size(TMR3D,1)*(1:(size(TMR3D,1)-2))+1;
'Averaging recurring pixels'
tic
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
toc
clear recurring_pixels
%% Filtering images with a laplacian of a gaussian for a 150 x 150 px field of view
dxy = 700;
TMR3Dsmall = uint16(NaN(dxy+1,dxy+1,Z));
TMR3Dfilter = uint16(NaN(dxy+1,dxy+1,Z));
p = 3;
con = 8;
kk = 8;
for kk=1:Z; % go through each plane in the fluorecent images
    kk;
    TMR = TMR3D(xy(1):xy(1)+dxy,xy(2):xy(2)+dxy,kk); % take only small fraction of the image
    TMR3Dsmall(:,:,kk) = TMR; % load into 3D stack
    TMR1 = imfilter(TMR, w1); % gauss smoothing
    TMR2 = imfilter(TMR1,w2); % spot filtering
    TMR2(1:5,:) = 0; % set corner pixels to zero
    TMR2(:,1:5) = 0;
    TMR2(dxy-3:dxy+1,:) = 0;
    TMR2(:,dxy-3:dxy+1) = 0;
    TMR3Dfilter(:,:,kk) = TMR2; % collect all filtered images into 3D stack
%    figure(100+kk); clf; subplot(1,2,1); imshow(TMR,[]); subplot(1,2,2); imshow(TMR2,[]); impixelinfo;
end;
TMR3Dfilter(1:fs,:,:) = 0; % set corners and borders to zero
TMR3Dfilter(:,1:fs,:) = 0;
TMR3Dfilter(dxy-fs:dxy+1,:,:) = 0;
TMR3Dfilter(:,dxy-fs:dxy+1,:) = 0;
TMRmaxf = max(TMR3Dfilter,[],3); %% Max projection original 3D stack
figure(200); clf; subplot(1,2,1); imshow(max(TMR3Dsmall,[],3),[]); subplot(1,2,2); imshow(TMRmaxf,[]); impixelinfo; title(T);% show results imshow TMRmaxf changed from [0 500] to []

for kk=1:Z;
    thall(kk,1) = std2(TMR3Dfilter(:,:,kk));
    thall(kk,2) = mean2(TMR3Dfilter(:,:,kk)); % determine threshold in each plane from the STD of background pixels;
end;% determine threshold in each plane from the STD of background pixels;

[M,I] = nanmax(thall(:,2)) % determine maimum threshold in each image
% figure(1); clf; subplot(1,2,1); imshow(TMR3D(minI:maxI,minI:maxI,I),[2000 30000]); subplot(1,2,2);imshow(TMR3Dfilter(:,:,I),[0 2000]); title(['STD: ' num2str(thallstd(I,1))]); impixelinfo ;
clear num
%th1 = 200:200:4000;

%% define variables
ii = 2;
xAll = NaN(1000,5);
yAll = NaN(1000,5);
th1 = thA(1);
dth = thA(2);
th2 = thA(3);
thh = [th1:dth:th2];
xx = NaN(100000,size(thh,2));
yy = NaN(100000,size(thh,2));
IMall = NaN(dxy+1,dxy+1,size(thh,2));
ii = 5;
I
size(thh,2);

%% Go through the different thresholds
'generating threshold spots'
for th5 = 5:size(thh,2)
    th = thh(th5); % pick threshold
    if max_int_thres
        if max_int_spot_th
             IM = immultiply(max(TMR3Dfilter,[],3),max(TMR3Dfilter,[],3) > th); % generate FISH image with only pixels above the threshold
             %    figure; imshow(IM,[]);
             if sum(IM(:)) > 0; % test that the images has spots
                 IM2 = imregionalmax(IM,26);       % 6 18 26     % determine the brightes pixel for each RNA spot in 3D
                 IM3 = immultiply(max(TMR3Dfilter,[],3), IM2);      % generate image including only the brightes pixel and their pixel values
                 m1(th5) = min(IM3(:));                          % determine minimum intentisy pixel
                 m2(th5) = max(IM3(:))./2;                       % determine maximum intensity pixel
                 [y1,x1] = find(IM3 > round(std2(max(TMR3Dfilter,[],3))));        % find all pixel coordianates for pixels above previously selected threshold
                 %[y1,x1] = find(IM3 > round(thall(j,2)));        % find all pixel coordianates for pixels above previously selected threshold
                 %    [y1,x1] = find(IM3 > round(thall(th5,2)));
                 xx(1:size(x1),th5) = x1;                        % save all x coordiantes for each threshold
                 yy(1:size(x1),th5) = y1;                %y coordiantes for each threshold
             end
        else                     
            index1 = 1;
               IM = immultiply(TMR3Dfilter(:,:,5:end),TMR3Dfilter(:,:,5:end) > th); % generate FISH image with only pixels above the threshold
               if sum(IM(:)) > 0; % test that the images has spots
                   IM2 = imregionalmax(IM,26);       % 6 18 26     % determine the brightes pixel for each RNA spot in 3D
                   IM3 = immultiply(TMR3Dfilter(:,:,5:end), IM2);      % generate image including only the brightes pixel and their pixel values
                   m1(th5) = min(IM3(:));                          % determine minimum intentisy pixel
                   m2(th5) = max(IM3(:))./2;                       % determine maximum intensity pixel
                   [y1,x1,z1] = find(IM3 > round(min(thall(:,2))));% round(thall(j,2)));        % find all pixel coordianates for pixels above previously selected threshold
                   %    [y1,x1] = find(IM3 > round(thall(th5,2)));
                   xx(index1:(index1+size(x1)-1),th5) = x1;                        % save all x coordiantes for each threshold
                   yy(index1:(index1+size(x1)-1),th5) = y1;                % y coordiantes for each threshold
                   %               IMall(:,:,th5) = IM3;                                         % (didn't adapt to max intensity) save all intensties for each threshold
                   index1 = index1+size(x1);
               end
            %% Old way of doing slice by slice
%             for j = 5:(Z-5)
%                 IM = immultiply(TMR3Dfilter(:,:,j),TMR3Dfilter(:,:,j) > th); % generate FISH image with only pixels above the threshold
%                 %    figure; imshow(IM,[]);
%                 if sum(IM(:)) > 0; % test that the images has spots
%                     IM2 = imregionalmax(IM,26);       % 6 18 26     % determine the brightes pixel for each RNA spot in 3D
%                     IM3 = immultiply(TMR3Dfilter(:,:,j), IM2);      % generate image including only the brightes pixel and their pixel values
%                     m1(th5) = min(IM3(:));                          % determine minimum intentisy pixel
%                     m2(th5) = max(IM3(:))./2;                       % determine maximum intensity pixel
%                     [y1,x1] = find(IM3 > round(thall(j,2)));        % find all pixel coordianates for pixels above previously selected threshold
%                     %    [y1,x1] = find(IM3 > round(thall(th5,2)));
%                     xx(index1:(index1+size(x1)-1),th5) = x1;                        % save all x coordiantes for each threshold
%                     yy(index1:(index1+size(x1)-1),th5) = y1;                % y coordiantes for each threshold
%                     %               IMall(:,:,th5) = IM3;                                         % (didn't adapt to max intensity) save all intensties for each threshold
%                     index1 = index1+size(x1);
%                     
%                 end
%             end
            %%
        end
    else
        IM = immultiply(TMR3Dfilter(:,:,I),TMR3Dfilter(:,:,I) > th); % generate FISH image with only pixels above the threshold
        %    figure; imshow(IM,[]);
        if sum(IM(:)) > 0; % test that the images has spots
            IM2 = imregionalmax(IM,26);       % 6 18 26     % determine the brightes pixel for each RNA spot in 3D
            IM3 = immultiply(TMR3Dfilter(:,:,I), IM2);      % generate image including only the brightes pixel and their pixel values
            m1(th5) = min(IM3(:));                          % determine minimum intentisy pixel
            m2(th5) = max(IM3(:))./2;                       % determine maximum intensity pixel
            [y1,x1] = find(IM3 > round(thall(I,2)));        % find all pixel coordianates for pixels above previosly selected threshold
            %    [y1,x1] = find(IM3 > round(thall(th5,2)));
            xx(1:size(x1),th5) = x1;                        % save all x coordiantes for each threshold
            yy(1:size(x1),th5) = y1;                        % save all y coordiantes for each threshold
            IMall(:,:,th5) = IM3;                           % save all intensties for each threshold
        else;
        end;
        %     figure(ii); clf; imshow(IMall(:,:,th5),[0 1]);impixelinfo;
        ii = ii + 1;
    end
end;
'showing memory used'
whos 

% clear RGB;
%i = round(size(thh,2)/2);

%% plot results
%[i3,i] = find(thh == ths);
%figure(111); clf; imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,2),yy(:,2),'or','markersize',10); title(strTMR); impixelinfo;
% m1 = round(min(min(TMR3Dfilter(:,:,I)))*2)
% m2 = round(max(max(TMR3Dfilter(:,:,I)))/2)
I
% figure(1111); clf; 
% subplot(1,2,1);
% imshow(log(double(TMR3Dfilter(:,:,I))),[1 9]); hold on; plot(xx(:,I),yy(:,I),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(I)) ', ' T]); impixelinfo;        
% subplot(1,2,2);
% imshow(log(double(TMR3Dsmall(:,:,I))),[]); hold on; plot(xx(:,I),yy(:,I),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(I)) ', ' T]); impixelinfo;        

%%% This below section commented out 1/12/2017 (BK) and replaced with code
%%% further below. if you want to pick the visualization thresholds
%%% manually do this
% if Ych2 == 3;
% LminF = 0  % min(min(TMR3Dfilter(:,:,I)))
% LmaxF = 200 %75  % max(max(TMR3Dfilter(:,:,I)))/1
% Lmin = 100 %100 % min(min(TMR3Dsmall(:,:,I)))
% Lmax = 6000 %600 % max(max(TMR3Dsmall(:,:,I)))/1
% elseif Ych2 == 5;
% LminF = 0 % min(min(TMR3Dfilter(:,:,I)))
% LmaxF = 500 % max(max(TMR3Dfilter(:,:,I)))*8
% Lmin = 1000 %1000 % min(min(TMR3Dsmall(:,:,I)))
% Lmax = 10000 % max(max(TMR3Dsmall(:,:,I)))*8
% end;
%%%
TMR3Dsmaller = TMR3Dsmall(:,:,5:(size(TMR3Dsmall,3)-5));
max_proj = double(max(TMR3Dsmaller,[],3));                                    %generate maximum intensity projection for original stack
max_proj_f = double(max(TMR3Dfilter,[],3));                                 %generate maximum intensity projection for filtered stack


LminF = median(max_proj_f(:))-round(0*std(max_proj_f(:))); % min(min(TMR3Dfilter(:,:,I)))
LmaxF = median(max_proj_f(:))+round(10*std(max_proj_f(:))); %75  % max(max(TMR3Dfilter(:,:,I)))/1
%Lmin = median(max_proj(:))-round(10*std(max_proj_f(:))); %100 % min(min(TMR3Dsmall(:,:,I)))
Lmin = min(max_proj(:)) %100 % min(min(TMR3Dsmall(:,:,I)))
Lmax = median(max_proj(:))+round(10*std(max_proj(:))); %600 % max(max(TMR3Dsmall(:,:,I)))/1


figure(1112); clf; 
subplot(1,2,1);
imshow(max_proj_f,[LminF LmaxF]); hold on; plot(xx(:,I),yy(:,I),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(I)) ', ' T]); impixelinfo;       
%imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,I),yy(:,I),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(I)) ', ' T]); impixelinfo;       
subplot(1,2,2);
imshow(max_proj,[Lmin Lmax]); hold on; plot(xx(:,I),yy(:,I),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(I)) ', ' T]); impixelinfo;        
%imshow(max_proj,[Lmin Lmax]); hold on; plot(xx(:,I),yy(:,I),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(I)) ', ' T]); impixelinfo;        


%imshow((double(TMR3Dsmall(:,:,I))),[]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%imshow(log(double(TMR3Dfilter(:,:,I))),[3 10]); hold on; plot(xx(:,2),yy(:,2),'or','markersize',10); title(strTMR); impixelinfo;
% figure(111); clf; imshow(IMall(:,:,2),[0 1]); hold on; plot(xx(:,2),yy(:,2),'o'); title(strTMR); impixelinfo;

I
buu = 5;
% m1(I)
% m2(I)
%% select the best threshold by decreasing (<) or increasing (>) the threshold
ik = I
sgtitle(['You can store low, mid, and high thresholds with 1,2, and 3. Middle mouse button exits', char(10) 'Change Threshold: z <<, x <<, c <, v >, b >>, n >>>'])
while buu == 5
    [xi,yi,but] = ginput(1)
    if but == 118; % > 1 th-step up
        ik = min(thA(3),ik + 1);
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo; 
        if max_int_thres
                  subplot(1,2,1);
        imshow(max_proj_f,[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
        subplot(1,2,2);
        imshow(max_proj,[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
        else
        subplot(1,2,1);
        imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
        subplot(1,2,2);
%        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
        imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
%         imshow(IMall(:,:,i),[0 1]); hold on; plot(xx(:,i),yy(:,i),'o'); title([strTMR ,' plane: ' num2str(i)]); impixelinfo; 
        end
       pause(1)
    elseif but == 99; % < 1 th-step down
        ik = max(thA(1),ik - 1);
        if max_int_thres
                  subplot(1,2,1);
        imshow(max_proj_f,[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
        subplot(1,2,2);
        imshow(max_proj,[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
        else
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
%        imshow(log(double(TMR3Dfilter(:,:,I))),[2 11]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
        subplot(1,2,1);
        imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
        subplot(1,2,2);
%        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
        imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
        end
        pause(1)
    elseif but == 98; % >>  10 th-steps up
        ik = min(thA(3),ik + 10);
        if max_int_thres
        subplot(1,2,1);
        imshow(max_proj_f,[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
        subplot(1,2,2);
        imshow(max_proj,[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
        else
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
        subplot(1,2,1);
        imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
        subplot(1,2,2);
%        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
        imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
%         imshow(IMall(:,:,i),[0 1]); hold on; plot(xx(:,i),yy(:,i),'o'); title([strTMR ,' plane: ' num2str(i)]); impixelinfo;  
        end
        pause(1)
    elseif but == 120; % << 10 th-steps down
        ik = max(thA(1),ik - 10);
        if max_int_thres
                  subplot(1,2,1);
        imshow(max_proj_f,[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
        subplot(1,2,2);
        imshow(max_proj,[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
        else
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
%        imshow(log(double(TMR3Dfilter(:,:,I))),[2 11]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
        subplot(1,2,1);
        imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
        subplot(1,2,2);
%        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
        imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
        end
        pause(1)
   elseif but == 110; % >>  100 th-steps up
        ik = min(thA(3),ik + 100);
        if max_int_thres
        subplot(1,2,1);
        imshow(max_proj_f,[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
        subplot(1,2,2);
        imshow(max_proj,[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
        else
%        imsho
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
        subplot(1,2,1);
        imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
        subplot(1,2,2);
%        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
        imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
%         imshow(IMall(:,:,i),[0 1]); hold on; plot(xx(:,i),yy(:,i),'o'); title([strTMR ,' plane: ' num2str(i)]); impixelinfo;
        end
        pause(1)
    elseif but == 122; % << 100 th-steps down
        ik = max(thA(1),ik - 100);
        if max_int_thres
        subplot(1,2,1);
        imshow(max_proj_f,[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
        subplot(1,2,2);
        imshow(max_proj,[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
        else
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
%        imshow(log(double(TMR3Dfilter(:,:,I))),[2 11]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
        subplot(1,2,1);
        imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
        subplot(1,2,2);
%        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
        imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
        end
        pause(1)
    elseif but == 49; % number '1' bottom low threshold
        th(1) = thh(ik)
        pause(1)
    elseif but == 50; % number '2' bottom medium threshold
        th(2) = thh(ik)
        pause(1)
    elseif but == 51; % number '3' bottom high threshold
        th(3) = thh(ik)
        pause(1)
    elseif but == 2; % middle mouse finish        
         buu = 2;

%   elseif but == 1; % left mouse bottom low threshold
%        th(1) = thh(ik)
%       pause(1)
%   elseif but == 3; % right mouse bottom high threshold
%       th(2) = thh(ik)
%       pause(1)
%   elseif but == 2; % middle mouse finish        
%         buu = 2;

    end;
end;

% th = thh(ik)
% th = input('Best Threshold:');

% %% Plot all planes with this threshold
% for kk=1:size(tmr1,2);
%     figure(10+kk);
%     imshow(log(double(TMR3Dfilter(:,:,I))),[2 11]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
% end;
