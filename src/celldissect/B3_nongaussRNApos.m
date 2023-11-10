function [CELLmaxRNA,RNApos]...
    =B3_nongaussRNApos(Lab,mm,TMR3Dorg,TMR3Dback,TMR3D3immax,TMR3Dfilter,Nuc3D)

%This function places the non gaussian fits of RNA spots into matrix RNApos
clear CELLmaxRNA RNApos CellRNAindex2  CellRNAorg2 Cellback2 Nuc
b = double(ones(9,9)); 
b(2:8,2:8) = 0; 
b1= (1-b)*0.5;
b2 = b1;
b2(3:7,3:7) = 2;
b2 = b2 ./2;
clear x3y3 g1 g2 g3 XY
zz = size(TMR3Dorg,3);
nz = size(Nuc3D,3);
PAR = NaN(mm,16,zz,25); %% (Number of cells,number of samples,number of planes, number of variables) 
%NN = round([2:zz+1]*nz/zz);
CELLmaxRNA = NaN(mm,3);
A = size(TMR3Dorg,1);
mm
j = 3; %76;


for j = 1:mm;
    %j = 18
    %tic; 
    clear k1 k2 Nuc Cyto NucTMR2d CellTMR1g CellTMR2d L CellTMR2dorg NucTMR2dorg CellRNA CellRNA2 CellRNAorg Nuc Cyto Cell CytoRNA2 NucRNA2 rLabeltcell NUMcell1 rLabeltcyto NUMcyto1 rLabeltnuc NUMnuc1;

    %% single cell
    k1 = Lab == j ; %== j;%sieve out dots in cell j.
    %figure; imshow(k1,[0 1]);
    k1=uint16(k1);
    k3=regionprops(k1,'BoundingBox','Area');
    k4 = k3.BoundingBox; %create the rectangular box around the cell. 
    X0=round(k4(1))-4;
    Y0=round(k4(2))-4;
    X1=round(k4(1)+ k4(3))+4;
    Y1=round(k4(2)+ k4(4))+4;
    if X0 < 1;
        X0 = 1;
    else;
    end;
    if Y0 < 1;
        Y0 = 1;
    else;
        
    end;
    if X1 > A;
        X1 = A;
    else;
    end;
    if Y1 > A;
        Y1 = A;
    else;
    end;
    k2 = k1(Y0:Y1,X0:X1);
% figure(100); clf; imshow(k1,[]); hold on; plot(X0,Y0,'o');hold on; plot(X1,Y1,'or'); title(['Cell: ' num2str(j)]); pixval;

    %% Separate Cell, Nuclear and cyto plasmic RNA
    i = 1;
    clear CellRNAfilter CellRNAindex CellRNAorg CellRNAback Nuc Cyto Cell 
    
    for i=1:zz;
        Tf = TMR3Dfilter(Y0:Y1,X0:X1,i);
        CellRNAfilter(:,:,i) = immultiply(Tf,k2); % Filtered image of RNA in the cell        
        CellRNAindex(:,:,i) = immultiply(TMR3D3immax(Y0:Y1,X0:X1,i),k2);  % Filtered image of RNA in the cell
        CellRNAorg(:,:,i) = immultiply(TMR3Dorg(Y0:Y1,X0:X1,i),k2); % Original RNA image of the cell
        CellRNAback(:,:,i) = immultiply(TMR3Dback(Y0:Y1,X0:X1,i),k2); % Original RNA image of the cell
        Nuc(:,:,i) = immultiply(double(Nuc3D(Y0:Y1,X0:X1,i)),double(k2)); % Nucleus object
        Cyto(:,:,i) = imsubtract(double(k1(Y0:Y1,X0:X1)),double(Nuc(:,:,i))); % Cytoplasm object
        Cell(:,:,i) = double(k1(Y0:Y1,X0:X1)); 
    end;

%% Cell, Cytoplasmic, Nuclear RNA
    clear CellRNAorg2 CellRNAindex2 Cellback2 CytoRNA2 NucRNA2
    CellRNAorg2 = immultiply(double(CellRNAorg),Cell);     
    CellRNAindex2 = immultiply(double(CellRNAindex),Cell); 
    Cellback2 = immultiply(double(CellRNAback),Cell); 
    CytoRNA2 = immultiply(double(CellRNAindex),double(Cyto)); 
    NucRNA2 = immultiply(double(CellRNAindex),double(Nuc)); 
      %figure(102); clf; imshow(k2 + uint16(max(CellRNAindex2,[],3)),[]); title(['Cell: ' num2str(j)]); pixval;

    %% RNA in the cell
    clear rLabeltcell NUMcell1 rLabeltcyto NUMcyto1 rLabeltnuc NUMnuc1
    [rLabeltcell,NUMcell1] = bwlabeln(CellRNAindex2 > 0,26); 
    [rLabeltcyto,NUMcyto1] = bwlabeln(CytoRNA2 > 0,26); 
    [rLabeltnuc,NUMnuc1] = bwlabeln(NucRNA2 > 0,26); 
    % figure(103); clf; imshow(k2 + uint16(max(rLabeltcell,[],3)),[]); title(['Cell: ' num2str(j)]); %pixval;
    
    %location of the RNA in three dimentions non-gaussian
   clear len
   if NUMcell1 ~= 0
    for i=1:NUMcell1;
       [row,col,vec]=ind2sub(size(rLabeltcell),find(rLabeltcell == i));
       if size(row,1)==1;
           y = row;
           x = col;
           z = vec;
       elseif size(row,1)> 1;  %if more than one pixel give the brightest one the position vecttor
           len=(size(row,1));
           for foo=1:len;
               kp(foo) = CellRNAindex2(row(foo),col(foo),vec(foo));
           end;
           [val, indx] = max(kp);
           y = row(indx);
           x = col(indx);
           z = vec(indx);
           clear kp
       end;
        RNApos(j,i,1) = x;
        RNApos(j,i,2) = y;
        RNApos(j,i,3) = z;
    end;
   else
        RNApos(j,:,1) = 0;
        RNApos(j,:,2) = 0;
        RNApos(j,:,3) = 0;
   end;
       

    %% Count RNA in the cell, cytoplasm, nucleus
    CELLmaxRNA(j,1) = NUMcell1;
    CELLmaxRNA(j,2) = NUMcyto1;
    CELLmaxRNA(j,3) = NUMnuc1;
    
end;
end
