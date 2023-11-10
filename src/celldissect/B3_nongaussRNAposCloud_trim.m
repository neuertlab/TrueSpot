function [CELLmaxRNA,RNApos,tran_cloud,clouds1,clouds2]...
    =B3_nongaussRNAposCloud(Lab,mm,TMR3Dorg,TMR3Dback,TMR3D3immax,TMR3Dfilter,Nuc3D)

%This function places the non gaussian fits of RNA spots into matrix RNApos
clear CELLmaxRNA RNApos CellRNAindex2  CellRNAorg2 Cellback2 Nuc tran_cloud clouds1 clouds2
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
j = 2; %76;
tic
counter56 = 1;
tran_cloud = zeros(4,mm);            %1st row is pixel size of sphere, 2nd row is original transcripts determined per cloud
clouds1 = zeros(size(TMR3Dorg));     %binary 3D image with all clouds
clouds2 = zeros(size(TMR3Dorg));     %binary 3D image with all clouds

    for j = 1:mm %Looping thru cells
        cloud2 = 0;         %Changes to 1 if there is a second cloud
        %j = 18
        j
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
        
        if X0 < 1
            X0 = 1;
        end
        if Y0 < 1
            Y0 = 1;
        end
        if X1 > A
            X1 = A;
        end
        if Y1 > A
            Y1 = A;
        end
        k2 = k1(Y0:Y1,X0:X1);

        %% Separate Cell, Nuclear and cyto plasmic RNA
        i = 1;
        clear CellRNAfilter CellRNAindex CellRNAorg CellRNAback Nuc Cyto Cell 
    
        for i=1:zz
            Tf = TMR3Dfilter(Y0:Y1,X0:X1,i);
            CellRNAfilter(:,:,i) = immultiply(Tf,k2); % Filtered image of RNA in the cell        
            CellRNAindex(:,:,i) = immultiply(TMR3D3immax(Y0:Y1,X0:X1,i),k2);  % Filtered image of RNA in the cell
            CellRNAorg(:,:,i) = immultiply(TMR3Dorg(Y0:Y1,X0:X1,i),k2); % Original RNA image of the cell
            CellRNAback(:,:,i) = immultiply(TMR3Dback(Y0:Y1,X0:X1,i),k2); % Original RNA image of the cell
            Nuc(:,:,i) = immultiply(double(Nuc3D(Y0:Y1,X0:X1,i)),double(k2)); % Nucleus object
            Cyto(:,:,i) = imsubtract(double(k1(Y0:Y1,X0:X1)),double(Nuc(:,:,i))); % Cytoplasm object
            Cell(:,:,i) = double(k1(Y0:Y1,X0:X1)); 
        end

        %% Cell, Cytoplasmic, Nuclear RNA
        clear CellRNAorg2 CellRNAindex2 Cellback2 CytoRNA2 NucRNA2
        CellRNAorg2 = immultiply(double(CellRNAorg),Cell);     
        CellRNAindex2 = immultiply(double(CellRNAindex),Cell); 
        Cellback2 = immultiply(double(CellRNAback),Cell); 
        CytoRNA2 = immultiply(double(CellRNAindex),double(Cyto));
        NucRNA2 = immultiply(double(CellRNAindex),double(Nuc));
        CytoRNA3 = immultiply(double(CellRNAorg),double(Cyto));
        NucRNA3 = immultiply(double(CellRNAorg),double(Nuc));
        CellRNAback = double(CellRNAback);
        CellRNAback(CellRNAback == 0) = NaN;
        CytoRNA3(CytoRNA3 == 0) = NaN;
        NucRNA3(NucRNA3 == 0) = NaN;


        %% RNA in the cell
        clear rLabeltcell NUMcell1 rLabeltcyto NUMcyto1 rLabeltnuc NUMnuc1
        [rLabeltcell,NUMcell1] = bwlabeln(CellRNAindex2 > 0,26); 
        [rLabeltcyto,NUMcyto1] = bwlabeln(CytoRNA2 > 0,26); 
        [rLabeltnuc,NUMnuc1] = bwlabeln(NucRNA2 > 0,26); 

        %%% Second try, start from maximum point and expand
        if NUMcell1 ~= 0
            z_adjust = 1;   %adjusts viewing window for z direction (due to slice size being different from pixel size)
            max_dxy = 50;
            fract_thres = .97;
            thres_devs = 3;        %This makes the threshold for what is determined to be the cloud to be the median+this many standard deviations
            [row,col,vec]=ind2sub(size(CellRNAorg),find(CellRNAorg == max(CellRNAorg(:))));  %Find brightest RNA pixel in original RNA image
            surrpix = CellRNAorg(row-1:row+1,col-1:col+1,vec);  %makes 3 x 3 matrix of surrounding pixels and brightest
            sumtot = sum(surrpix(:));   %sum of all pixls
            while sumtot < 11*CellRNAorg(row,col,vec)/3   %if average pixel intensity of surrounding pixels is less than 1/3 brightest pixel
                CellRNAorg(row,col,vec) = floor((sumtot-CellRNAorg(row,col,vec))/8);   %replaces pixel with average of surrounding pixels on same plane
                [row,col,vec]=ind2sub(size(CellRNAorg),find(CellRNAorg == max(CellRNAorg(:))));  %Find brightest RNA pixel in original RNA image
                surrpix = CellRNAorg(row-1:row+1,col-1:col+1,vec);  %makes 3 x 3 matrix of surrounding pixels and brightest
                sumtot = sum(surrpix(:));   %sum of all pixls
            end      
            if size(row,1)==1
                y = row;
                x = col;
                z = vec;
            elseif size(row,1)> 1  %if more than one pixel give the brightest one the position vecttor
                %i
                len=(size(row,1));
                for foo=1:len
                    kp(foo) = CellRNAindex2(row(foo),col(foo),vec(foo));
                end
                [val, indx] = max(kp);
                y = row(indx);
                x = col(indx);
                z = vec(indx);
                clear kp
            end
            counter76 = 1;
            fluor_per =[0;0];
            new = 1;
            for dxy1=2:2:max_dxy
                dxy1
                area1 = zeros(size(CellRNAorg));
                [ny,nx,nz]=size(CellRNAorg);
                [xx1,yy1,zz1] = meshgrid((1:nx)-x,(1:ny)-y,((1:nz)-z)*z_adjust);
                area1_dist=sqrt(xx1.^2 + yy1.^2 + zz1.^2);
                area1 = area1_dist <= dxy1;

                CellRNAorg3 = immultiply(double(CellRNAorg),area1);
                CellRNAorg3(CellRNAorg3 == 0) = NaN;

                new_fluor_per = nansum(CellRNAorg3(:))/size(find(not(isnan(CellRNAorg3(:)))),1);
                if new ~= 1
                    if new_fluor_per >= (old_fluor_per*fract_thres) & new_fluor_per < (4 * max(Cellback2(:)));
                        break
                    end
                end
                new = 0;    
                old_fluor_per = nansum(CellRNAorg3(:))/size(find(not(isnan(CellRNAorg3(:)))),1);
                fluor_per(1,counter76) = dxy1;
                fluor_per(2,counter76) = old_fluor_per;
                counter76 = counter76+1;
            end

            CellRNAorgN = double(CellRNAorg);
            CellRNAorgN(CellRNAorgN == 0) = NaN;
            mm4 = nanmedian(CellRNAorgN(:))+thres_devs*nanstd(CellRNAorgN(:));  %Sets threshold to median plus 10 standard deviations
            CellRNAorg4 = CellRNAorg3 > mm4;

%%% Below is for filling holes int he image and only keeping the largest
%%% connected neighborhood for each slice
            CellRNAorg4A = zeros(size(CellRNAorg4));
            for k = 1:zz
                CellRNAorg4A(:,:,k) = imfill(CellRNAorg4(:,:,k),'holes');                                        % fill holes in image
            end
            dapi_label_low = bwlabeln(CellRNAorg4A(:,:,:),26);                                          % label different connected neighborhoods

            CellRNAorg4AA = dapi_label_low;                                        % "LOW"
            CellRNAorg4AA(CellRNAorg4AA==0)=[];
            CellRNAorg4A(:,:,:) = dapi_label_low == mode(CellRNAorg4AA);                                        % only keep biggest connected neighborhood 

            clouds1(Y0:Y1,X0:X1,:) = clouds1(Y0:Y1,X0:X1,:)+double(CellRNAorg4A);
            tran_cloud(1,counter56) = dxy1;                             %the radius of sphere
            trans_cl = immultiply(CellRNAorg4A,CellRNAindex2);
            tran_cloud(3,counter56) = size(find(trans_cl(:)>0),1);
            CellRNAorg5 = zeros(size(CellRNAorg3));
            CellRNAorg5A = zeros(size(CellRNAorg5));
            mask = CellRNAorg4A == 0;    %All areas outside first cloud
            temp_Labeltcell= immultiply(rLabeltcell,mask);
        end
        
        if NUMcell1 > 1 & max(temp_Labeltcell(:)) > 0 %Look for another cloud
            cloud2 = 1;
            temp_CellRNAorg = immultiply(CellRNAorg,mask);  %only include areas outside first cloud       
            [row,col,vec]=ind2sub(size(temp_CellRNAorg),find(temp_CellRNAorg == max(temp_CellRNAorg(:))));  %Find brightest RNA pixel in original RNA image
            surrpix = temp_CellRNAorg(row-1:row+1,col-1:col+1,vec);  %makes 3 x 3 matrix of surrounding pixels and brightest
            sumtot = sum(surrpix(:));   %sum of all pixls
            while sumtot < 11*temp_CellRNAorg(row,col,vec)/3   %if average pixel intensity of surrounding pixels is less than 1/3 brightest pixel
                temp_CellRNAorg(row,col,vec) = floor((sumtot-temp_CellRNAorg(row,col,vec))/8);   %replaces pixel with average of surrounding pixels on same plane
                [row,col,vec]=ind2sub(size(temp_CellRNAorg),find(temp_CellRNAorg == max(temp_CellRNAorg(:))));  %Find brightest RNA pixel in original RNA image
                surrpix = temp_CellRNAorg(row-1:row+1,col-1:col+1,vec);  %makes 3 x 3 matrix of surrounding pixels and brightest
                sumtot = sum(surrpix(:));   %sum of all pixls
            end      
            if size(row,1)==1
                y1 = row;
                x1 = col;
                z1 = vec;
            elseif size(row,1)> 1  %if more than one pixel give the brightest one the position vecttor
                %i
                len=(size(row,1));
                for foo=1:len
                    kp(foo) = CellRNAindex2(row(foo),col(foo),vec(foo));
                end
                [val, indx] = max(kp);
                y1 = row(indx);
                x1 = col(indx);
                z1 = vec(indx);
                clear kp
            end
            counter76 = 1;
            fluor_per =[0;0];
            new = 1;
            for dxy1=2:2:max_dxy
                dxy1
                area1 = zeros(size(CellRNAorg));
                [ny,nx,nz]=size(CellRNAorg);
                [xx1,yy1,zz1] = meshgrid((1:nx)-x1,(1:ny)-y1,((1:nz)-z1)*z_adjust);
                area1_dist=sqrt(xx1.^2 + yy1.^2 + zz1.^2);
                area1 = area1_dist <= dxy1;

                CellRNAorg3 = immultiply(double(CellRNAorg),area1);
                CellRNAorg3(CellRNAorg3 == 0) = NaN;

                new_fluor_per = nansum(CellRNAorg3(:))/size(find(not(isnan(CellRNAorg3(:)))),1);
                if new ~= 1
                    if new_fluor_per >= (old_fluor_per*fract_thres) & new_fluor_per < (4 * max(Cellback2(:)))
                        break
                    end
                end
                new = 0;    
                old_fluor_per = nansum(CellRNAorg3(:))/size(find(not(isnan(CellRNAorg3(:)))),1);
                fluor_per(1,counter76) = dxy1;
                fluor_per(2,counter76) = old_fluor_per;
                counter76 = counter76+1;
            end

            CellRNAorgN = double(CellRNAorg);
            CellRNAorgN(CellRNAorgN == 0) = NaN;
            mm4 = nanmedian(CellRNAorgN(:))+thres_devs*nanstd(CellRNAorgN(:));  %Sets threshold to median plus 10 standard deviations
            CellRNAorg5 = CellRNAorg3 > mm4;
%%% Below is for filling holes int he image and only keeping the largest
%%% connected neighborhood for each slice
            CellRNAorg5A = zeros(size(CellRNAorg5));
            for k = 1:zz
                CellRNAorg5A(:,:,k) = imfill(CellRNAorg5(:,:,k),'holes');                                        % fill holes in image
            end
            dapi_label_low = bwlabeln(CellRNAorg5A(:,:,:),26);                                          % label different connected neighborhoods
            CellRNAorg5AA = dapi_label_low;                                        % "LOW"
            CellRNAorg5AA(CellRNAorg5AA==0)=[];
            CellRNAorg5A(:,:,:) = dapi_label_low == mode(CellRNAorg5AA);                                        % only keep biggest connected neighborhood
%%%
            bot = z1-dxy1/2-2;
            top = z1+dxy1/2+2;
            if bot < 1
                bot = 1;
            end
            if top > zz
                top = zz;
            end
 
            clouds2(Y0:Y1,X0:X1,:) = clouds2(Y0:Y1,X0:X1,:)+double(CellRNAorg5A);
            tran_cloud(2,counter56) = dxy1;                             %the radius of sphere
            trans_cl = immultiply(CellRNAorg5A,CellRNAindex2);
            tran_cloud(4,counter56) = size(find(trans_cl(:)>0),1);
            CellRNAorg4A = CellRNAorg4A+CellRNAorg5A;
        end  
        counter56 = counter56+1;

    %location of the RNA in three dimentions non-gaussian
        clear len
        counter80 = 1;
        if NUMcell1 ~= 0
            RNApos(j,1,1) = x;                                                 %Add the brightest intensity pixel as the first spot (first column
            RNApos(j,1,2) = y;
            RNApos(j,1,3) = z;
            counter80 = 2;                                                     %This counter will be for the RNA in the cell, so RNA within clouds are eliminated 
            if cloud2
                RNApos(j,2,1) = x1;                                                 %Add the brightest intensity pixel of second cloud as second entry
                RNApos(j,2,2) = y1;
                RNApos(j,2,3) = z1;
                counter80 = 3; 
            end
            for i=1:NUMcell1
                [row,col,vec]=ind2sub(size(rLabeltcell),find(rLabeltcell == i));   %This still goes through every determined spot
                    if size(row,1)==1
                        y = row;
                        x = col;
                        z = vec;
                    elseif size(row,1)> 1  %if more than one pixel give the brightest one the position vector
                        %i
                        len=(size(row,1));
                        for foo=1:len
                            kp(foo) = CellRNAindex2(row(foo),col(foo),vec(foo));
                        end
                        [val, indx] = max(kp);
                        y = row(indx);
                        x = col(indx);
                        z = vec(indx); 
                        clear kp 
                    end
                    if CellRNAorg4A(y,x,z) == 0 & CellRNAorg5A(y,x,z) == 0          %check if RNA is part of cloud before adding it (RNAs within cloud are not added other than the first)
                        RNApos(j,counter80,1) = x;
                        RNApos(j,counter80,2) = y;
                        RNApos(j,counter80,3) = z;
                        counter80 = counter80+1;
                    else
                        CellRNAindex2(rLabeltcell == i) = 0;
                        CytoRNA2(CytoRNA2 == i) = 0;
                        NucRNA2(NucRNA2 == i) = 0;
                    end
            end
        else
            RNApos(j,:,1) = 0;
            RNApos(j,:,2) = 0;
            RNApos(j,:,3) = 0;
        end
        if NUMcell1 ~= 0
            CellRNAindex2(RNApos(j,1,2),RNApos(j,1,1),RNApos(j,1,3)) = 1;
            if cloud2
                CellRNAindex2(RNApos(j,2,2),RNApos(j,2,1),RNApos(j,2,3)) = 1;
            end
        end
       %% RNA in the cell
        clear rLabeltcell NUMcell1 rLabeltcyto NUMcyto1 rLabeltnuc NUMnuc1
        [rLabeltcell,NUMcell1] = bwlabeln(CellRNAindex2 > 0,26);
        CytoRNA2 = immultiply(double(CellRNAindex2),double(Cyto));
        NucRNA2 = immultiply(double(CellRNAindex2),double(Nuc));
        [rLabeltcyto,NUMcyto1] = bwlabeln(CytoRNA2 > 0,26); 
        [rLabeltnuc,NUMnuc1] = bwlabeln(NucRNA2 > 0,26); 
       

    %% Count RNA in the cell, cytoplasm, nucleus
        CELLmaxRNA(j,1) = NUMcell1;
        CELLmaxRNA(j,2) = NUMcyto1;
        CELLmaxRNA(j,3) = NUMnuc1;
    
    end

toc
end
