function [Label_low,Label_mid,Label_hi,nuclei] = B1_autosegment_nuclei5(DAPI_ims,dapi_label,Yim,dapi_threshold,min_nucleus_size,max_nucleus_size);
% Yim = 1;
% min_nucleus_size = 50;                                                     % minimum nuclear size
% max_nucleus_size = 20000;                                                   % maximum nuclear size

a = size(DAPI_ims,1);                                                       % image size in number of pixels
z = size(DAPI_ims,3);                                                       % stack size in number of images

i = 100;
dxy1 = 11;                                                                  % 7 ; define number of pixels from the max DAPI pixel'
dxy2 = dxy1;
dxy = dxy1;
dapi_cell2 = NaN(1+2*dxy1,1+2*dxy2,z);                                      % generate NaN 3D matrix with 2*dxy in x and y direction
Label_low = zeros(a,a,z);                                                   % generate zero 3D matrix of the size of the image stack
Label_mid = zeros(a,a,z);
Label_hi = zeros(a,a,z);
Label_index = zeros(a,a);
DAPI_ims_max = max(DAPI_ims,[],3);
% go through all the DAPI nucleus labels one by one
m22 = max(dapi_label(:))                                                   % determine maximum number of DAPI stained nuclei
for i = 1:m22;                                                              
    i
    k1 = dapi_label == i;                                                  % Single cell per image
    k2=uint16(k1);                                                          % convert to 16bit image
    k2a = immultiply(k1,DAPI_ims_max);
    m111 = max(k2a(:));
    sumposx = 0;                                                             %will be numerator for calculating weighted average for x
    sumposy = 0;                                                            %will be numerator for calculating weighted average for y
    sumfluor = 0;                                                            %will be denominator for calculating weighted average
    %k2atest = zeros(1,1024*1024);
    counter = 1;
    for j = 1:size(k2a,1)
        for k = 1:size(k2a,2)
            %k2atest(counter) = k1(j,k);
            counter = counter+1;
            sumposx = sumposx + (k1(j,k)*j);
            sumposy = sumposy + (k1(j,k)*k);
            sumfluor = sumfluor + k1(j,k);
        end
    end
    xA = round(sumposx/sumfluor);
    yA = round(sumposy/sumfluor);
    
    %[xA,yA] = find(k2a == m111);
%     k3=regionprops(k2,'Centroid');                                          % center of maximun in 3D called centroid
%     k4 = k3.Centroid;                                                      
%     yA = round(k4(1,1));                                                    % round yA to full pixel
%     xA = round(k4(1,2));                                                    % round xA to full pixel
    Imzero = uint16(zeros(a,a));
    Imzero(xA,yA) = 1;
    area1 = zeros(a,a,z);                                                   % generate zero 3D matrix of the size of the image stack
    A = xA-dxy1;                                                            % generate corner points with are dxy1 pixels away from the nuclear center
    B = xA+dxy1;
    C = yA-dxy2;
    D = yA+dxy2;
    
    for j = [1,numel(A)];                                                   % if corner A is below 1 pixel or negative set equals 1
        if A(j) < 1;
            A(j) = 1;
        end
    end
    for j = [1,numel(B)];                                                   % if corner B is above a pixels (1024 or 2048) set equals a
        if B(j) > a;
            B(j) = a;
        end
    end
    for j = [1,numel(C)];                                                   % if corner C is below 1 pixel or negative set equals 1
        if C(j) < 1;
            C(j) = 1;
        end
     end
    for j = [1,numel(D)];                                                   % if corner D is above a pixels (1024 or 2048) set equals a
        if D(j) > a;
            D(j) = a;
        end
    end
%     if  A < 1;                                                              % if corner A is below 1 pixel or negative set equals 1
%         A = 1;
%     else;
%     end;
%     if B > a;                                                               % if corner B is above a pixels (1024 or 2048) set equals a
%         B = a;
%     else;
%     end;
%     if C < 1;                                                               % if corner C is below 1 pixel or negative set equals 1
%         C = 1;
%     else;
%     end;
%     if D > a;                                                                % if corner D is above a pixels (1024 or 2048) set equals a
%         D = a;
%     else;
%     end;
    area1(A:B,C:D,:) = 1;                                                   % set pixels within A,B,C,D equals 1
%     figure(10); clf; imshow(max(area1,[],3),[]); impixelinfo;
    dapi_cell = immultiply(uint16(DAPI_ims),uint16(area1));                         % generate dapi image stack inside area1 and set the rest of the image to zero
    dapi_cell_max = max(dapi_cell,[],3);                                    % determine maximum projection
%     figure(12); clf; imshow(dapi_cell_max,[]); impixelinfo;
%     dapi_cell2 = dapi_cell(A:B,C:D,:);
    dapi_cell2max = dapi_cell_max(A:B,C:D);                                 % cut out the maxiumum projection
%     figure(13); clf; imshow(max(dapi_cell2,[],3),[]); impixelinfo;
%     figure(14); clf; imshow(dapi_cell2max,[]); impixelinfo;
    mm1 = min(dapi_cell2max(:));                                            % determine minimum dapi signal
    mm2 = max(dapi_cell2max(:));                                            % determine maximum dapi signal
    mmd = mm2-mm1;                                                          % determine the difference between maximum and minimum
    m1 = mm1+mmd*0.4;                                                       % set "LOW" threshold to 40% of difference
    m2 = mm1+mmd*0.5;                                                       % set "MID" threshold to 50% of difference
    m3 = mm1+mmd*0.6;                                                       % set "HI" threshold to 60% of difference
    dapi_label3A = dapi_cell > (m1);                                        % "LOW" generate 3D binary image above the thresholds
    dapi_label3B = dapi_cell > (m2);                                        % "MID" generate 3D binary image above the thresholds
    dapi_label3C = dapi_cell > (m3);                                        % "HI" generate 3D binary image above the thresholds
    Label_low = Label_low + dapi_label3A;                                   % Add to label matrix 'LOW'
    Label_mid = Label_mid + dapi_label3B;                                   % Add to label matrix 'MID'
    Label_hi = Label_hi + dapi_label3C;                                     % Add to label matrix 'HI'
    Label_index = uint16(Label_index) + Imzero;
end;
if Yim == 1;
    figure(21); clf; imshow(max(Label_low,[],3),[]); impixelinfo; imwrite(uint16(max(Label_low,[],3)),['Label_low.tif']); 
    figure(22); clf; imshow(max(Label_mid,[],3),[]); impixelinfo; imwrite(uint16(max(Label_mid,[],3)),['Label_mid.tif']); 
    figure(23); clf; imshow(max(Label_hi,[],3),[]); impixelinfo; imwrite(uint16(max(Label_hi,[],3)),['Label_hi.tif']); 
    figure(24); clf; imshow(Label_index,[]); impixelinfo; imwrite(uint16(Label_index),['Label_index.tif']); 
else;
end;
nuclei = Label_index;                                                       % max projection of the segmented dapi signals using the 50% threshold
