function [th,RNA_thresholds,counter6] = AB2_Threshold3Dim_auto(TMR_ims,xy,thA,T,Ych2,counter,S,RNA_thresholds,counter6,cells,exp_name,outfile_prefix_RNA)
TMR3D = TMR_ims;

%% filter sets up for gaussian filtering based on experimental PSF
load('recurring pixels 8 out of 12')
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

%% Filtering images with a laplacian of a gaussian for a 150 x 150 px field of view
% dxy = 300;
% TMR3Dsmall = uint16(NaN(dxy+1,dxy+1,Z));
% TMR3Dfilter = uint16(NaN(dxy+1,dxy+1,Z));
TMR3Dsmall = uint16(NaN(size(TMR3D,1),size(TMR3D,2),Z));
TMR3Dfilter = uint16(NaN(size(TMR3D,1),size(TMR3D,2),Z));
p = 3;
con = 8;
kk = 1;
exmat(1:size(TMR3D,1),1) = 1:size(TMR3D,1);                                 %exmat is a matrix that contains the border indices, which will not be filtered                                                  
exmat(1:size(TMR3D,1),4) = (size(TMR3D,1)^2-(size(TMR3D,1)-1)):size(TMR3D,1)^2;
exmat(1:2,2:3) = 1;
exmat(3:size(TMR3D,1),2) = size(TMR3D,1)*(2:(size(TMR3D,1)-1));
exmat(3:size(TMR3D,1),3) = size(TMR3D,1)*(1:(size(TMR3D,1)-2))+1;
for kk=1:Z; % go through each plane in the fluorecent images
    kk;
    %TMR = TMR3D(xy(1):xy(1)+dxy,xy(2):xy(2)+dxy,kk); % take only small fraction of the image
    TMR = TMR3D(:,:,kk);
%    figure; imshow(TMR,[0 2000]);impixelinfo;
%     TMR0 = imfilter(TMR, w3); % finding dead pixels
%     hi_pixels = find(TMR > prctile(TMR(:),99.95));
%     hi_pixels12 = hi_pixels;
%     for j = 1:size(hi_pixels,1)
%         a1 = hi_pixels(j);                      
%         a2 = exmat == a1;
%         if max(a2(:)) == 0;                                             %excludes pixels on the border
%             surr_pix = [a1-1,a1+1,a1-2048, a1-2049,a1-2047,a1+2048,a1+2047,a1+2049]; %surrounding pixels
%             TMR(a1) = mean(TMR(surr_pix));                                    %replace pixel with mean of surrounding pixels
%         end
%     end
    for j = 1:size(recurring_pixels,1)
        a1 = recurring_pixels(j);                      
        a2 = exmat == a1;
        if max(a2(:)) == 0;                                             %excludes pixels on the border
            surr_pix = [a1-1,a1+1,a1-2048, a1-2049,a1-2047,a1+2048,a1+2047,a1+2049]; %surrounding pixels
            TMR(a1) = mean(TMR(surr_pix));                                    %replace pixel with mean of surrounding pixels
        end
    end
%    figure; imshow(TMR,[0 2000]);impixelinfo;
%     TMR_bin = bwareaopen(TMR,1);                                            %BK 5/21/15 remove dots that are only single pixels
%     TMR = immultiply(TMR_bin,TMR);                                          %BK 5/21/15 multiply the binary by original image
    TMR3Dsmall(:,:,kk) = TMR; % load into 3D stack
    TMR1 = imfilter(TMR, w1); % gauss smoothing
    %figure; imshow(TMR1,[]);impixelinfo;
    TMR2 = imfilter(TMR1,w2); % spot filtering
    %figure; imshow(TMR2,[]);impixelinfo;
    TMR2(1:5,:) = 0; % set corner pixels to zero
    TMR2(:,1:5) = 0;
%     TMR2(dxy-3:dxy+1,:) = 0;
%     TMR2(:,dxy-3:dxy+1) = 0;
    TMR2(size(TMR,1)-3:size(TMR,1),:) = 0;
    TMR2(:,size(TMR,2)-3:size(TMR,2)) = 0;
    TMR3Dfilter(:,:,kk) = TMR2; % collect all filtered images into 3D stack
%    figure(100+kk); clf; subplot(1,2,1); imshow(TMR,[]); subplot(1,2,2); imshow(TMR2,[]); impixelinfo;
end;
TMR3Dfilter(1:fs,:,:) = 0; % set corners and borders to zero
TMR3Dfilter(:,1:fs,:) = 0;
TMR3Dfilter(size(TMR,1)-fs:size(TMR,1),:,:) = 0;
TMR3Dfilter(:,size(TMR,1)-fs:size(TMR,1),:) = 0;
TMRmaxf = max(TMR3Dfilter,[],3); %% Max projection original 3D stack
%figure(200); clf; subplot(1,2,1); imshow(max(TMR3Dsmall,[],3),[]); subplot(1,2,2); imshow(TMRmaxf,[0 500]); impixelinfo; title(T);% show results

for kk=1:Z;
    thall(kk,1) = std2(TMR3Dfilter(:,:,kk));
    thall(kk,2) = mean2(TMR3Dfilter(:,:,kk)); % determine threshold in each plane from the STD of background pixels;
end;% determine threshold in each plane from the STD of background pixels;

[M,I] = nanmax(thall(:,2)) % determine maimum threshold in each image
% figure(1); clf; subplot(1,2,1); imshow(TMR3D(minI:maxI,minI:maxI,I),[2000 30000]); subplot(1,2,2);imshow(TMR3Dfilter(:,:,I),[0 2000]); title(['STD: ' num2str(thallstd(I,1))]); impixelinfo ;
clear num
%th1 = 200:200:4000;
%% Code to Store Pixel intensities %%
filename = ['PixRNA_' exp_name num2str(counter) '_Ych2_' num2str(Ych2) '.mat']
%filename = ['PixRNA_' exp_name '_' num2str(counter) '_Ych2_' num2str(Ych2) '.mat']
if exist([outfile_prefix_RNA filename])
    load([outfile_prefix_RNA filename])
else
PixRNA = zeros(2048,2048);                    %Stores max pixels for each image: x axis, y-axis, timepoint, image num, channel)
PixRNA1d = zeros(2048*2048,1);
end

TMRmax = uint16(TMRmaxf); % convert 3D image sections into 16bit format to save memory
CellPix = immultiply(TMRmax,uint16(cells > 0));
PixRNA(1:2048,1:2048) = CellPix;
PixRNA1d(1:2048*2048) = CellPix(:);
save([outfile_prefix_RNA filename],'PixRNA','PixRNA1d','-v7.3');
clear PixRNA
clear PixRNAd
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
IMall = NaN(size(TMR,1),size(TMR,1),size(thh,2));

ii = 5;
I
size(thh,2);
%% Go through the different thresholds
counter5 = 1;
SpecRNA = zeros(max(cells(:))+1,size(thh,2)+1);                             %This will store how many spots are in each cell in the image
for th5 = 1:size(thh,2)
    th = thh(th5); % pick threshold
    IM = immultiply(TMR3Dfilter(:,:,I),TMR3Dfilter(:,:,I) > th); % generate FISH image with only pixels above the threshold
%    figure; imshow(IM,[]);
    if sum(IM(:)) > 0; % test that the images has spots
    IM2 = imregionalmax(IM,26);       % 6 18 26     % determine the brightes pixel for each RNA spot in 3D    
    IM3 = immultiply(TMR3Dfilter(:,:,I), IM2);      % generate iamge including only the brightes pixel and their pixel values
    m1(th5) = min(IM3(:));                          % determine minimum intentisy pixel
    m2(th5) = max(IM3(:))./2;                       % determine maximum intensity pixel
    [y1,x1] = find(IM3 > round(thall(I,2)));        % find all pixel coordianates for pixels above previosly selected threshold 
%    [y1,x1] = find(IM3 > round(thall(th5,2)));
    xx(1:size(x1),th5) = x1;                        % save all x coordiantes for each threshold
    yy(1:size(x1),th5) = y1;                        % save all y coordiantes for each threshold
    IMall(:,:,th5) = IM3;                           % save all intensties for each threshold
    else;
    IM2 = zeros(2048,2048);
    end;
%     figure(ii); clf; imshow(IMall(:,:,th5),[0 1]);impixelinfo;
    ii = ii + 1;
    RNA_num(counter5) = numel(find(IM2 == 1));
    SpecRNA(1,counter5) = th;
    for j = 1:max(cells(:))
        temp_cell = cells == j;
        temp_RNA = immultiply(uint16(IM2),uint16(temp_cell));
        SpecRNA(j+1,counter5) = numel(find(temp_RNA == 1));
    end 
    counter5 = counter5 + 1
end;
save(['SpecRNA ' exp_name num2str(counter) '_Ych2_' num2str(Ych2) '.mat'],'SpecRNA','-v7.3')
thh_RNAnum = NaN(2,size(thh,2));
thh_RNAnum(1,:) = thh;
thh_RNAnum(2,:) = RNA_num;
save(strcat('IMall ',exp_name, '_', num2str(counter), '_Ych2_', num2str(Ych2),'.mat'),'IMall')
save(strcat('RNA spot vs thres ', exp_name, '_', num2str(counter), '_Ych2_', num2str(Ych2),'.mat'), 'thh_RNAnum') 
figure(7);clf;plot(thh,RNA_num)
filename = ['RNA num vs thres task_id ' num2str(S) ' img ' num2str(counter) ' channel number ' num2str(Ych2) '.fig']
savefig(filename);
for stdv = .01%*2.^(1:10)
%%%%%% BK Plot the Inverse CV %%%%%%%%%%%% 
cv_width = 10
 inv_cv = zeros(size(thh)-cv_width);
 inv_cv_x = zeros(size(thh)-cv_width);
 
for k = [1:(size(thh,2)-cv_width)]                                   % populate variables for inv cv 
     std_cv = std(RNA_num(k:k+cv_width));
    if std_cv < stdv;
        std_cv = stdv;
     end
     inv_cv(k) = mean(RNA_num(k:k+cv_width))/std_cv;    % find the inverse coefficient of variation (inverse cv) across a certain number of timepoints
     inv_cv_x(k) = mean([thh(k),thh(k+cv_width)]);         % the mean of the points in the calculation will be the x value
     
end 
figure(8);clf;plot(inv_cv_x,inv_cv)
filename = ['inv cv width ' num2str(cv_width) ' task_id ' num2str(S) ' img ' num2str(counter) ' max threshold ' num2str(thA(3)) ' channel number ' num2str(Ych2) ' low std ' num2str(stdv) '.fig']
savefig(filename);
RNA_thresholds(1,counter6) = S;                         %task id number
RNA_thresholds(2,counter6) = counter;                   %image number
RNA_thresholds(3,counter6) = Ych2;                      %channel number
RNA_thresholds(4,counter6) = stdv;                      %min stdv
RNA_thresholds(5,counter6) = thA(3);                    %threshold cutoff
RNA_thresholds(6,counter6) = thh(find(inv_cv == max(inv_cv),1,'first'));    %threshold determined by inv cv
RNA_thresholds(7,counter6) = RNA_num(find(inv_cv == max(inv_cv),1,'first')); %RNA number determined by inv cv
counter6 = counter6 + 1;                                %add a counter for RNA_thresholds to move to next spot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
if Ych2 == 3;
LminF = 0  % min(min(TMR3Dfilter(:,:,I)))
LmaxF = 400 %100 %75  % max(max(TMR3Dfilter(:,:,I)))/1
Lmin = 100 %30 %100 % min(min(TMR3Dsmall(:,:,I)))
Lmax = 4000 %1000 %600 % max(max(TMR3Dsmall(:,:,I)))/1
elseif Ych2 == 5;
LminF = 0 % min(min(TMR3Dfilter(:,:,I)))
LmaxF = 400 % max(max(TMR3Dfilter(:,:,I)))*8
Lmin = 100 %1000 % min(min(TMR3Dsmall(:,:,I)))
Lmax = 4000 % max(max(TMR3Dsmall(:,:,I)))*8
else
LminF = 0 % min(min(TMR3Dfilter(:,:,I)))
LmaxF = 400 % max(max(TMR3Dfilter(:,:,I)))*8
Lmin = 100 %1000 % min(min(TMR3Dsmall(:,:,I)))
Lmax = 4000 % max(max(TMR3Dsmall(:,:,I)))*8    
end;

% % figure(1112); clf; 
% % subplot(1,2,1);
% % imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,I),yy(:,I),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(I)) ', ' T]); impixelinfo;        
% % subplot(1,2,2);
% % imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,I),yy(:,I),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(I)) ', ' T]); impixelinfo;        
%imshow((double(TMR3Dsmall(:,:,I))),[]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%imshow(log(double(TMR3Dfilter(:,:,I))),[3 10]); hold on; plot(xx(:,2),yy(:,2),'or','markersize',10); title(strTMR); impixelinfo;
% figure(111); clf; imshow(IMall(:,:,2),[0 1]); hold on; plot(xx(:,2),yy(:,2),'o'); title(strTMR); impixelinfo;
I
buu = 5;
% m1(I)
% m2(I)
%% select the best threshold by decreasing (<) or increasing (>) the threshold
ik = find(inv_cv == max(inv_cv),1,'first')
%  imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%         subplot(1,2,1);
%         imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
%         subplot(1,2,2);
% %        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%         imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
% %         imshow(IMall(:,:,i),[0 1]); hold on; plot(xx(:,i),yy(:,i),'o'); title([strTMR ,' plane: ' num2str(i)]); impixelinfo;
% filename = ['RNA image width ' num2str(cv_width) ' task_id ' num2str(S) ' img ' num2str(counter) ' max threshold ' num2str(thA(3)) ' channel number ' num2str(Ych2) ' low std ' num2str(stdv) '.fig']
% savefig(filename);
th = thh(ik);
end

% while buu == 5
%     [xi,yi,but] = ginput(1)
%     if but == 118 & ik < max(thA); % > 1 th-step up
%         pause(1)
%         ik = ik + 1;
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%         subplot(1,2,1);
%         imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
%         subplot(1,2,2);
% %        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%         imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
% %         imshow(IMall(:,:,i),[0 1]); hold on; plot(xx(:,i),yy(:,i),'o'); title([strTMR ,' plane: ' num2str(i)]); impixelinfo;        
%         pause(1)
%     elseif but == 99 & ik > min(thA); % < 1 th-step down
%         pause(1)
%         ik = ik - 1;
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
% %        imshow(log(double(TMR3Dfilter(:,:,I))),[2 11]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
%         subplot(1,2,1);
%         imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
%         subplot(1,2,2);
% %        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%         imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
%     elseif but == 98 & ik <= (max(thA)-10); % >>  10 th-steps up
%         pause(1)
%         ik = ik + 10;
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%         subplot(1,2,1);
%         imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
%         subplot(1,2,2);
% %        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%         imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
% %         imshow(IMall(:,:,i),[0 1]); hold on; plot(xx(:,i),yy(:,i),'o'); title([strTMR ,' plane: ' num2str(i)]); impixelinfo;        
%         pause(1)
%     elseif but == 120 & ik >= (min(thA)+10); % << 10 th-steps down
%         pause(1)
%         ik = ik - 10;
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
% %        imshow(log(double(TMR3Dfilter(:,:,I))),[2 11]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
%         subplot(1,2,1);
%         imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
%         subplot(1,2,2);
% %        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%         imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
%     elseif but == 110 & ik <= (max(thA)-100); % >>  100 th-steps up
%         pause(1)
%         ik = ik + 100;
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%         subplot(1,2,1);
%         imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo; % plot the image
%         subplot(1,2,2);
% %        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%         imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;      % plot the identified spots in color on top
% %         imshow(IMall(:,:,i),[0 1]); hold on; plot(xx(:,i),yy(:,i),'o'); title([strTMR ,' plane: ' num2str(i)]); impixelinfo;        
%         pause(1)
%     elseif but == 122 & ik >= (min(thA)+100); % << 100 th-steps down
%         pause(1)
%         ik = ik - 100;
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
% %        imshow(log(double(TMR3Dfilter(:,:,I))),[2 11]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
%         subplot(1,2,1);
%         imshow((double(TMR3Dfilter(:,:,I))),[LminF LmaxF]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
%         subplot(1,2,2);
% %        imshow(log(double(TMR3Dfilter(:,:,I))),[4 8]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
%         imshow((double(TMR3Dsmall(:,:,I))),[Lmin Lmax]); hold on; plot(xx(:,ik),yy(:,ik),'or','markersize',10); title([' plane: ' num2str(I) ', th: ' num2str(thh(ik)) ', ' T]); impixelinfo;        
% %        imshow(log(double(TMRmaxf)),[3 10]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;    
% 
%     elseif but == 49; % number '1' bottom low threshold
%         th(1) = thh(ik)
%         pause(1)
%     elseif but == 50; % number '2' bottom medium threshold
%         th(2) = thh(ik)
%         pause(1)
%     elseif but == 51; % number '3' bottom high threshold
%         th(3) = thh(ik)
%         pause(1)
%     elseif but == 2; % middle mouse finish        
%          buu = 2;
% 
% %   elseif but == 1; % left mouse bottom low threshold
% %        th(1) = thh(ik)
% %       pause(1)
% %   elseif but == 3; % right mouse bottom high threshold
% %       th(2) = thh(ik)
% %       pause(1)
% %   elseif but == 2; % middle mouse finish        
% %         buu = 2;
% 
%     end;
% end;

% th = thh(ik)
% th = input('Best Threshold:');

% %% Plot all planes with this threshold
% for kk=1:size(tmr1,2);
%     figure(10+kk);
%     imshow(log(double(TMR3Dfilter(:,:,I))),[2 11]); hold on; plot(xx(:,i),yy(:,i),'or','markersize',10); title([strTMR ,' plane: ' num2str(i) ', th: ' num2str(thh(i))]); impixelinfo;        
% end;
