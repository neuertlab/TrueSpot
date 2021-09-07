function [Lab,plane,CellInfo,trans_plane,CellImage] = B2_autosegment_cells_new(TRANS_ims,nuclei,auto_mode,Yim,min_cell_size,max_cell_size,dapi_label,img_stacks);  
% auto_mode = 'last5';                                                        % segmentation mode 
% Yim = 1;                                                                    % Show images
% min_cell_size = 1000;                                                       % 250; % minimum cell size
% max_cell_size = 10000;                                                      % 1300;% maximum cell size
A = size(TRANS_ims,1);                                                      % image size in number of pixels
Z = size(TRANS_ims,3);                                                      % stack size in number of images
% f = msgbox('Segmenting Cell Boundaries');
%{ 
   auto_mode is either 'midplane' or 'max_cells'
   The plane for segmentation is automatically chosen.  One method for this
   is to segment each stack and choose the one that results the most cells.
   This type of auto-segment is done by invoking this function with:
    auto-mode= max_cells

   The alternative is to find the lowest-contrast plane, which in yeast is the
   very middle of the cell, and then use the max-projection of the 5 planes
   above the middle (plane +3 through +7 above the middle).
%}
%nuclei_max = max(nuclei,[],3);                                              % Maximum intensity projection of the nulei
nuclei_max = dapi_label;
%nuclei_max = DAPI_ims_added;
if Yim == 1;
    figure(51); clf; imshow(nuclei_max,[]); impixelinfo; imwrite(uint16(nuclei_max),['nuclei_max.tif']); 
else;
end;

%% Select segmentation mode

% segment cell using maxium contrast
if strcmp(auto_mode,'max_cells')                                            
    plane =C2_find_trans_bestplane(TRANS_ims,nuclei_max);                   % identifies the best plane
    plane;
    trans_plane = TRANS_ims(:,:,plane);                                     % choose the trans image to find cell boundary
if Yim == 1;
        figure(52); clf; imshow(trans_plane,[]); title('max'); impixelinfo; imwrite(uint16(trans_plane),['trans_plane.tif']);     % comment out for analysis
else;
end;
    trans2 = imopen(trans_plane,strel('disk',10)); % 25                     % smooth trans image with a disk filter of 50px to generate an background image
if Yim == 1;
        figure(53); clf; imshow(trans2,[]);impixelinfo; imwrite(uint16(trans2),['trans2.tif']);
else;
end;
    trans3 = imsubtract(trans_plane,trans2);                                % subtract background for the real image to enhance signal
if Yim == 1;
        figure(54); clf; imshow(trans3,[]);impixelinfo; imwrite(uint16(trans3),['trans3.tif']);
else;
end;
    trans_plane = double(trans3);                                           % convert image into double format

% segment cell using the mid plane with the lowest contrast
elseif strcmp(auto_mode,'midplane')
    plane= C3_find_trans_midplane(TRANS_ims);                               % find plane with the lowest contrast
    plane
if Yim == 1;
        figure(52); clf; imshow(TRANS_ims(:,:,plane),[]);impixelinfo; imwrite(uint16(TRANS_ims(:,:,plane),[]),['TRANS_ims.tif']);
else;
end;
    start_plane = plane + 3;                                                % consider 3 planes above the focus plane to see the cell boundary
    end_plane = plane + 7;                                                  % consider 7 planes above the focus plane to see the cell boundary
    if end_plane> size(TRANS_ims,3)
        end_plane = size(TRANS_ims,3);
    end
    if start_plane > size(TRANS_ims,3)
        start_plane = size(TRANS_ims,3);
    end
    trans_ring_planes = TRANS_ims(:,:,...
                           start_plane:end_plane);
    trans_plane = max(trans_ring_planes,[],3);                              % Maximum intensity z-projection
    
if Yim == 1;
        figure(53); clf; imshow(trans_plane,[]);  title('mid');  impixelinfo; imwrite(uint16(trans_plane),['trans_plane.tif']);% comment out for analysis
else;
end;
    trans2 = imopen(trans_plane,strel('disk',10)); % 25                     % smooth trans image with a disk filter of 50px to generate an background image
if Yim == 1;
        figure(54); clf; imshow(trans2,[]);impixelinfo; imwrite(uint16(trans2),['trans2.tif']);
else;
end;
    trans3 = imsubtract(trans_plane,trans2);                                % subtract background for the real image to enhace signal
if Yim == 1;
        figure(55); clf; imshow(trans3,[]);impixelinfo; imwrite(uint16(trans3),['trans3.tif']);
else;
end;
    trans_plane = double(trans3);                                           % convert image into double format

elseif strcmp(auto_mode,'midplane2')
    [ims2,plane] = C4_find_trans_plane(TRANS_ims);
    plane
if Yim == 1;
        figure(52); clf; imshow(TRANS_ims(:,:,plane),[]);impixelinfo; imwrite(uint16(TRANS_ims(:,:,plane),[]),['TRANS_ims.tif']);
else;
end;
    start_plane = plane + 3;                              % consider 3 planes above the focus plane to see the cell boundary
    end_plane = plane + 7;                                % consider 7 planes above the focus plane to see the cell boundary
    if end_plane> size(TRANS_ims,3)
        end_plane = size(TRANS_ims,3);
    end
    if start_plane > size(TRANS_ims,3)
        start_plane = size(TRANS_ims,3);
    end
    trans_ring_planes = TRANS_ims(:,:,...
                           start_plane:end_plane);
    trans_plane = max(trans_ring_planes,[],3);            % Maximum intensity z-projection
    trans_plane = ims2;
if Yim == 1;
        figure(53); clf; imshow(trans_plane,[]);  title('mid2');  impixelinfo; imwrite(uint16(trans_plane),['trans_plane.tif']);  % comment out for analysis
else;
end;
    trans2 = imopen(trans_plane,strel('disk',10)); % 25                                   % smooth trans image with a disk filter of 50px to generate an background image
if Yim == 1;
        figure(54); clf; imshow(trans2,[]);impixelinfo; imwrite(uint16(trans2),['trans2.tif']);
else;
end;
    trans3 = imsubtract(trans_plane,trans2);                                                  % subtract background for the real image to enhace signal
if Yim == 1;
        figure(55); clf; imshow(trans3,[]);impixelinfo; imwrite(uint16(trans3),['trans3.tif']);
else;
end;
    trans_plane = double(trans3);                                                        % convert image into double format
    
elseif strcmp(auto_mode,'first5')
    plane = 1;
    start_plane = 1;                              % consider 3 planes above the focus plane to see the cell boundary
    end_plane = 5;                                % consider 7 planes above the focus plane to see the cell boundary
    if end_plane> size(TRANS_ims,3)
        end_plane = size(TRANS_ims,3);
    end
    if start_plane > size(TRANS_ims,3)
        start_plane = size(TRANS_ims,3);
    end
    trans_ring_planes = TRANS_ims(:,:,...
                           start_plane:end_plane);
    trans_plane = max(trans_ring_planes,[],3);            % Maximum intensity z-projection
if Yim == 1;    
        figure(52); clf; imshow(trans_plane,[]);   title('first5');   impixelinfo; imwrite(uint16(trans_plane),['trans_plane.tif']);% comment out for analysis
else;
end;
    trans2 = imopen(trans_plane,strel('disk',10)); % 25                                   % smooth trans image with a disk filter of 50px to generate an background image
if Yim == 1;
        figure(53); clf; imshow(trans2,[]);impixelinfo; imwrite(uint16(trans2),['trans2.tif']);
else;
end;
    trans3 = imsubtract(trans_plane,trans2);                                                  % subtract background for the real image to enhace signal
if Yim == 1;
        figure(54); clf; imshow(trans3,[]);impixelinfo; imwrite(uint16(trans3),['trans3.tif']);
else;
end;
    trans_plane = double(trans3);                                                        % convert image into double format

    
    
elseif strcmp(auto_mode,'last5')
    plane = Z;
    start_plane = img_stacks(1);                              % consider 3 planes above the focus plane to see the cell boundary
                                                    %Changed to 9 for
                                                    %mammalian cells BK
                                                    %4/7/16
    end_plane = img_stacks(2);                                % consider 7 planes above the focus plane to see the cell boundary
    if end_plane> size(TRANS_ims,3)
        end_plane = size(TRANS_ims,3);
    end
    if start_plane > size(TRANS_ims,3)
        start_plane = size(TRANS_ims,3);
    end
    trans_ring_planes = TRANS_ims(:,:,...
                           start_plane:end_plane);
    trans_plane = max(trans_ring_planes,[],3);            % Maximum intensity z-projection
else%     trans1 = trans_plane; %imopen(trans_plane,strel('disk',5)); % 25 %Added Gaussian smoothing
if Yim == 1;
        figure(52); clf; imshow(trans1,[]);   title('last5');   impixelinfo; imwrite(uint16(trans1),['trans1.tif']); % comment out for analysis
else;
end;
    trans2 = imopen(trans_plane,strel('disk',round(sqrt((min_cell_size)/(pi))/2))); % 25 (when yeast?) % Changed from 40 (when doing mESCs) to being dependent on cell size
%     trans2 = imopen(trans_plane,strel('disk',100)); % 25 (when yeast?) % Changed from 40 (when doing mESCs) to being dependent on cell size
if Yim == 1;
        figure(53); clf; imshow(trans2,[]);impixelinfo; imwrite(uint16(trans2),['trans2.tif']);
else;
end;
    trans2 = imopen(trans2,strel('disk',10)); % 25                                   % smooth YFP image with a disk filter of 50px to generate an background image
if Yim == 1;
        figure(99); clf; imshow(trans2,[]);impixelinfo; imwrite(uint16(trans2),['trans2b.tif']);
else;
end;
    trans3 = imsubtract(trans1,trans2);                                                  % subtract background for the real image to enhace signal
if Yim == 1;
        figure(54); clf; imshow(trans3,[]);impixelinfo; imwrite(uint16(trans3),['trans3.tif']);
else;
end;
    trans_plane = double(trans3);                                                        % convert image into double format

    fprintf(['invalid segmentation mode;' ...
             'must be either "midplane" or "max_cells"\n']);
    return;
end


%% segment using Matlab's watershed, relying on DAPI nuclei
'size trans_plane'
size(trans_plane)
'size nuclei_max'
size(nuclei_max)
gradmag2 = imimposemin(trans_plane, nuclei_max,26);                            % generate superimposed image of the trans signal and the DAPI stained nucleus       
if Yim == 1;
    figure(56); clf; imshow(gradmag2,[]); impixelinfo; imwrite(uint16(gradmag2),['gradmag2.tif']) % comment out for analysis
else;
end;
bck = mode(gradmag2(:));
bckgrnd = gradmag2 == bck;
if Yim == 1;
    figure(97); clf; imshow(bckgrnd,[]); impixelinfo; imwrite(uint16(bckgrnd),['bckgrnd.tif']); % comment out for analysis
else;
end;
cells = watershed(gradmag2);                                                % segment cell using water shed algorithem 
if Yim == 1;
    figure(57); clf; imshow(cells,[]); impixelinfo; imwrite(uint16(cells),['cells.tif']); % comment out for analysis
else;
end;

cells_normal = bwareaopen(cells, min_cell_size);                            %filter out cells that are too small
cells_huge = bwareaopen(cells, max_cell_size);                              %filter out cells that are too large
cells_filtered = cells_normal - cells_huge;                                 % remaining cells have the right size
if Yim == 1;
    figure(58); clf; imshow(cells_filtered,[]); impixelinfo; imwrite(uint16(cells_filtered),['cells_filtered.tif']); % comment out for analysis
else;
end;

cells= bwlabeln(cells_filtered);                                            % labels cells

Label1 = uint16(cells);                                                     % convert semented cell image to 16bit
if Yim == 1;
    figure(img_stacks(2)); clf; imshow(Label1,[]); impixelinfo;  imwrite(uint16(Label1),['Label1.tif']);
else;
end;
m = max(Label1(:))                                                          % determine the maximum number of segmented cells 

%% Go through each cell and analyze the shape, remove cells and get information
%Label1 = Label
for j = 1:m                                                                 % m number of cells
    j;
    k1 = Label1 == j;                                                       % Single cell per image
    k11 = double(k1>0);
    %kkk = k1.*j;
    k2=uint16(k1);
    k3=regionprops(k2,'BoundingBox','Area','Centroid');                     % cut out rectange with cell / increases post processing
    k4 = k3.BoundingBox;                                                    %create the rectangular box around the cell 
    k5(j) = k3.Area;
    X0=round(k4(1));
    Y0=round(k4(2));
    X1=round(k4(1)+ k4(3))-1;
    Y1=round(k4(2)+ k4(4))-1;
    
    %% remove cells that are too small or big
    if (k5(j)> max_cell_size)|(sum(k5(j)< min_cell_size ))
            kkk = uint16(k11) * j;
            Label1 = Label1 - kkk;   
    else
            Label1 = Label1;
    end

    %% remove cells on the border of the image
    if (X0 <4)|(X1>A-4)|(Y0 <4)|(Y1>A-4)
        Label1 = Label1 - uint16(k11)*j;
    end
    
    %% analyze cell shape and remove cells that are too big
%     k6=regionprops(k2,'MajorAxisLength','MinorAxisLength');
%     k7 = k6.MajorAxisLength;
%     k8 = k6.MinorAxisLength;
%     if (k7/k8 > 1.5); % | (k7/k8 < 0.5);
%            Label1= Label1 - uint16(k11)*j;
%     end
    clear k1
    clear k2
    clear k3
    clear k4
    clear k
end

%% generate cell label
Labelf = zeros(A,A);
Label1 = Label1 > 0;                                                        % generate a binary image of the segmented cells
Label1 = imfill(Label1, 'holes');                                           % closes holes in the cells
Labelf = imadd(double(Labelf), double(Label1));                             % generate segmented cells with orginial intenity profile
if Yim == 1;
    figure(61);clf; imshow(Labelf,[]); impixelinfo; imwrite(uint16(Labelf),['Labelf.tif']); % comment out for analysis
else;
end;
Lab = Labelf > 0;                                                           % make binary image
if Yim == 1;
    figure(62);clf; imshow(Lab,[0 2]); imwrite(uint16(Lab),['Lab.tif']);
else;
end;
[Lab,mm] =bwlabel(Label1,8);                                                % index cells
Lab = uint16(Lab);                                                          % convert image to 16bit format to reduce image size and save memory
mm
j = 3;

%% Analyze cell shape for further analysis
for j = 1 : mm;                                                              % organise cell information
    j;
    % Determine cell properties
    L1 = Lab == j;%sieve out dots in cell j
    L2 = double(L1>0);
    L3 = uint16(L2);
    %figure; imshow(L3,[]); pixval;
    try L4 = regionprops(L3,'Centroid','MajorAxisLength','MinorAxisLength','FilledArea','Image'); end; % erxtract several properties from each cell
    try CellInfo(j,[1:2]) = L4.Centroid; catch CellInfo(j,[1:2]) = NaN; end;                               % Cell Center
    try CellInfo(j,3) = L4.MajorAxisLength; catch CellInfo(j,3) = NaN; end;                            % Cell long axis
    try CellInfo(j,4) = L4.MinorAxisLength; catch CellInfo(j,4) = NaN; end;                            % Cell short Axis
    try CellInfo(j,5) = L4.FilledArea; catch CellInfo(j,5) = NaN; end;                                 % Cell Area
    try CellImage(1:size(L4.Image,1),1:size(L4.Image,2),j) = L4.Image; end; % Images of the cells; can be used for cell cycle analysis
    
end;


%-------------------------- End of function
%------------------ random analysis/visualization code below

%{
figure(1);
subplot(1,2,1);
imshow(max(dapi4,[],3),[0 1000]);

subplot(1,2,2);
imshow(max(nuclei,[],3),[0 1]);

figsize=get(figure(1),'Position');
figsize(3)=800;
figsize(4)=800;
set(figure(1),'Position',figsize);

%}

%{

%analyze sizes of the nuclei, plotting a histogram

nuclei_lab = bwlabeln(nuclei);
figure,imshow(max(nuclei_lab,[],3),[0 1])

%for i=1:max(nuclei_lab(:))

stats = regionsprops(nuclei_lab,'Area');

areas = [stats.Area];

[area_counts,area_xaxis]=hist(areas,0:100:3000);


bar(area_xaxis,area_counts);
title_str=sprintf('DAPI Areas, total: %g',size(areas,2));
title(title_str);
xlabel('Area');
ylabel('Total DAPI Count');

%set(gca,'yscale','log');

%}

