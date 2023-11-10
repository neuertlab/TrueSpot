function [thCY7Man,thAF700Man,thCY5Man,thAF594Man,thTMRMan,thYFPMan,thGFPMan,RNA_thresholds,counter6] = AB_FindThreshold_TMR_AF594_CY5_B_auto...
    (counter,cells,trans_plane,S,CY7_ims,AF700_ims,CY5_ims,AF594_ims,TMR_ims,YFP_ims,GFP_ims,thA,ths,Ych,RNA_thresholds,counter6,exp_name,outfile_prefix_RNA);%(Im1,cells,trans_plane,S,CY7_ims,AF700_ims,CY5_ims,AF594_ims,TMR_ims,YFP_ims,GFP_ims,thA,ths,Ych);


h =sprintf('%03d',counter);  % defines label
mm = max(cells(:)); % maximum number of cells
%figure(100); clf;subplot(1,2,1); imshow(trans_plane,[]); subplot(1,2,2); imshow(cells,[0 1]); ...
 %   title(['S: ' num2str(S) ', Im: ' h ', Cells: ' num2str(mm)]); impixelinfo;  hold on; 

% [xi,yi,but] = ginput(1); % reads the mouse curser position
xi = 1;yi = 1;but = 1; 
if but == 1;
    plot(xi,yi,'ro'); % plots an red circle around the choosen spot
    xy(2) = round(xi); % get the x and y coordinates
    xy(1) = round(yi);

    %% filter CY7 images
    if false; %Ych(1) == 1;
        T = 'filter CY7 images'
        Ych2 = 1;
        [thCY7Man,RNA_thresholds,counter6] = AB2_Threshold3Dim_BK(CY7_ims,xy,thA,T,Ych2,counter,S,RNA_thresholds,counter6,cells,exp_name,outfile_prefix_RNA);
    else;
        thCY7Man = NaN;
    end;

    %% filter AF700 images
    if false; %Ych(2) == 1;
        T = 'filter AF700 images'
        Ych2 = 2;
        [thAF700Man,RNA_thresholds,counter6] = AB2_Threshold3Dim_BK_auto(AF700_ims,xy,thA,T,Ych2,counter,S,RNA_thresholds,counter6,cells,exp_name,outfile_prefix_RNA);
    else;
        thAF700Man = NaN;
    end;


    %% filter CY5 images
    if Ych(3) == 1; 
        T = 'filter CY5 images'
        Ych2 = 3;
        [thCY5Man,RNA_thresholds,counter6] = AB2_Threshold3Dim_BK_auto(CY5_ims,xy,thA,T,Ych2,counter,S,RNA_thresholds,counter6,cells,exp_name,outfile_prefix_RNA);
    else;
        thCY5Man = NaN;
    end;

    %% filter AF594 images
    if Ych(4) == 1;
        T = 'filter AF594 images'
        Ych2 = 4;
        [thAF594Man,RNA_thresholds,counter6] = AB2_Threshold3Dim_BK_auto(AF594_ims,xy,thA,T,Ych2,counter,S,RNA_thresholds,counter6,cells,exp_name,outfile_prefix_RNA);
    else;
        thAF594Man = NaN;
    end;

    %% filter TMR images
    if Ych(5) == 1;
        T = 'filter TMR images'
        Ych2 = 5;
        [thTMRMan,RNA_thresholds,counter6] = AB2_Threshold3Dim_BK_auto(TMR_ims,xy,thA,T,Ych2,counter,S,RNA_thresholds,counter6,cells,exp_name,outfile_prefix_RNA);
    else;
        thTMRMan = NaN;
    end;

    %% filter YFP images
    if Ych(6) == 1;
        T = 'filter YFP images'
        Ych2 = 6;
        [thYFPMan,RNA_thresholds,counter6] = AB2_Threshold3Dim_BK_auto(YFP_ims,xy,thA,T,Ych2,counter,S,RNA_thresholds,counter6,cells,exp_name,outfile_prefix_RNA);
    else;
        thYFPMan = NaN;
    end;

    %% filter GFP images
    if Ych(7) == 1;
        T = 'filter GFP images'
        Ych2 = 7;
        [thGFPMan,RNA_thresholds,counter6] = AB2_Threshold3Dim_BK_auto(GFP_ims,xy,thA,T,Ych2,counter,S,RNA_thresholds,counter6,cells,exp_name,outfile_prefix_RNA);
    else;
        thGFPMan = NaN;
    end;
else;
end;