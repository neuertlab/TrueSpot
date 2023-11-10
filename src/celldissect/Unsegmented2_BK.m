function [F_C_all,F_N_all,FN_C_all,FN_N_all,FP_C_all,FP_N_all,Prec_C_all,Prec_N_all,Sens_C_all,Sens_N_all,TP_C_all,TP_N_all] ...
    = Unsegmented2(trans_plane,cells,Label_mid,DAPI,CellInfo,NN,d1,Sc,LT,Font_size1,all_name,filenum,...
    F_C_all,F_N_all,FN_C_all,FN_N_all,FP_C_all,FP_N_all,Prec_C_all,Prec_N_all,Sens_C_all,Sens_N_all,TP_C_all,TP_N_all,nuc_num)

No = max(cells(:));
%d1 = 20;
DAPI1 = DAPI(d1*2-1).data;
fig = 1;

% Segmentation border
xx = 10;
px = size(DAPI1,1);
SB = zeros(px,px);
SB(xx:px-xx,xx:px-xx) = 1;
%figure(fig); imshow(SB,[]); impixelinfo; fig = fig + 1;
%%%scale
x1 = (px/2) - (px*Sc)/2 + 1
x2 = (px/2) + (px*Sc)/2
y1 = (px/2) - (px*Sc)/2 + 1
y2 = (px/2) + (px*Sc)/2
%%%
trans_plane = immultiply(trans_plane,SB);
TRANS2 = im2double(trans_plane./max(trans_plane(:)));
%figure(fig); clf; imshow(TRANS2,[]); title('trans'); impixelinfo; fig = fig + 1;
n = 3; m = 0;
Tm = mean2(TRANS2(:));
Ts = std2(TRANS2(:));
TRANS3 = imadjust(TRANS2,[Tm-m*Ts Tm+n*Ts],[]); %Tm-n*Ts
%TRANS3 = imadjust(TRANS2); % [0 1] add auto adjusting. 
figure(fig); clf; imshow(TRANS3(x1:x2,y1:y2,:)); title('trans'); impixelinfo; fig = fig + 1;


%figure(fig); clf; imshow(DAPI1,[]); title('DAPI'); impixelinfo; fig = fig + 1;
DAPI2 = im2double(DAPI1); % ./max(DAPI1(:))
DAPI2 = immultiply(DAPI2,SB);
figure(fig); clf; imshow(DAPI2(x1:x2,y1:y2,:),[]); title('DAPI'); impixelinfo; fig = fig + 1;
DAPI3 = imadjust(DAPI2); %[0 0.2] add auto adjusting. 
n = 3; m = 0;
Dm = mean2(DAPI2(:));
Ds = std2(DAPI2(:));
DAPI3 = imadjust(DAPI2,[Dm-m*Ds Dm+n*Ds],[]); %Dm-n*Ds
%figure(fig); clf; imshow(DAPI3); title('DAPI'); impixelinfo; fig = fig + 1;

cells = immultiply(cells,uint16(SB));
SEG1 = cells > 0; %.marker;
SEG2 = bwmorph(SEG1,'remove')*100000;
SEG3 = imadjust(SEG2,[])*100000;
SEG4 = bwmorph(SEG3,'thicken');
SEG4 = bwmorph(SEG4,'close');
    SEG4 = bwmorph(SEG4,'thicken');
    SEG4 = bwmorph(SEG4,'close');

% if LT > 1;
%     SEG4 = bwmorph(SEG4,'thicken');
%     SEG4 = bwmorph(SEG4,'close');
% else;
% end;
figure(fig); clf; imshow(SEG4(x1:x2,y1:y2,:)); title('Seg Cells'); impixelinfo; fig = fig + 1;

NUC0 = Label_mid(:,:,d1);
NUC0 = immultiply(NUC0,SB);
NUC1 = (NUC0) > 0; %.marker;
NUC2 = bwmorph(NUC1,'remove')*100000;
NUC3 = imadjust(NUC2,[])*100000;
NUC4 = bwmorph(NUC3,'thicken');
NUC4 = bwmorph(NUC4,'close');
NUC4 = bwmorph(NUC4,'thicken');
NUC4 = bwmorph(NUC4,'close');
% if LT > 1;
%     NUC4 = bwmorph(NUC4,'thicken');
%     NUC4 = bwmorph(NUC4,'close');
% else;
% end;
figure(fig); clf; imshow(NUC4(x1:x2,y1:y2,:)); title('Seg Nuc'); impixelinfo; fig = fig + 1;

%% Overlays + RGB display
R = TRANS3 + NUC4;
B = DAPI3 + NUC4 + TRANS3;
G = SEG4 + TRANS3;
RGB = cat(3,R,G,B);
x1 = (px/2) - (px*Sc)/2 + 1
x2 = (px/2) + (px*Sc)/2
y1 = (px/2) - (px*Sc)/2 + 1
y2 = (px/2) + (px*Sc)/2
figure(fig); clf; imshow(RGB(x1:x2,y1:y2,:),[]); impixelinfo; hold on; % display data
pos = find(CellInfo(:,1)-CellInfo(:,3)/2 > x1 & CellInfo(:,1)+CellInfo(:,3)/2 < x2 & CellInfo(:,2)-CellInfo(:,3)/2 > y1 & CellInfo(:,2)+CellInfo(:,3)/2 < y2);
NoS = size(pos,1);
for i = 1: size(pos,1); %No;
    i
    text(CellInfo(pos(i,1),1)-x1,CellInfo(pos(i,1),2)-y1,num2str(i),'Color','white','FontSize',Font_size1);
%    plot(CellInfo(i,1),CellInfo(i,2),[num2str(i)]); hold on;
end;
title({['\fontsize{8} 1 (cyan) - False Positive Nuc (:  automatic segmented area but not manual  ||  2 (magenta) - False Negative Nuc: manual segmented area but not automatic'],['\fontsize{8} 3 (yellow) - False Positive Cell (:  automatic segmented area but not manual  ||  4 (green) - False Negative Cell: manual segmented area but not automatic || c: Cells (not on border) Q: quit']});

%% Mark single cells;
n = 0;
FP_N = 0;    %False positive nucleus
FN_N = 0;   %False negative nucleus
FP_C = 0;    %false positive cell
FN_C = 0;   %false negative cell
TP_C = 0;   %true positive cell
while true; %n < NoS;
    [x,y,button] = ginput(1);
    n = n + 1;
    %     if button == 1;             % Left
    %         FP = FP + 1;
    %         text(x,y,[num2str(FP) ' FP'],'Color','red','FontSize',18) ; %num2str(n)
    %     elseif button == 3;
    %         FN = FN + 1;
    %         text(x,y,[num2str(FN) ' FN'],'Color','green','FontSize',18) ; %num2str(n)
    %     elseif button == 113;
    %         break; %return;
    %     end;
    if button == 49;             % number 1
        FP_N = FP_N + 1;
        text(x,y,[num2str(FP_N) ' FP_N'],'Color','cyan','FontSize',Font_size1) ; %num2str(n)
    elseif button == 50;        %number 2
        FN_N = FN_N + 1;
        text(x,y,[num2str(FN_N) ' FN_N'],'Color','magenta','FontSize',Font_size1) ; %num2str(n)
    elseif button == 51;        %number 3
        FP_C = FP_C + 1;
        text(x,y,[num2str(FP_C) ' FP_C'],'Color','yellow','FontSize',Font_size1) ; %num2str(n)
    elseif button == 52;        %number 4
        FN_C = FN_C + 1;
        text(x,y,[num2str(FN_C) ' FN_C'],'Color','green','FontSize',Font_size1) ; %num2str(n)
%     elseif button == 99;        %lowercase c
%         TP_C = TP_C + 1;
%         text(x,y,['C' num2str(TP_C)],'Color','white','FontSize',Font_size1) ; %num2str(n)
    elseif button == 113;       %Q
        break; %return;
    end;
end;
%msgbox('Done collecting points');

% %% Mark single cells;
% n = 0;
% [x,y,button] = ginput;
% n = n + 1;
% FPi = find(button == 1);
% FNi = find(button == 3);
% for j = 1:size(FNi,1);
%     text(x(FNi(j)),y(FNi(j)),[num2str(FNi(j)) ' FN'],'Color','green','FontSize',18) ; %num2str(n)
% end;
% for k = 1:size(FPi,1);
%     text(x(FPi(k)),y(FPi(k)),[num2str(FPi(k)) ' FP'],'Color','red','FontSize',18) ; %num2str(n)
% end;
% 
%% Quantification
% FP = size(FPi,1)
% FN = size(FNi,1)
%%% Old quantification with TP_C being every labeled cell
% NoS = TP_C;     %Set number of cells equal to user-defined function
%%% Updated
TP_N = nuc_num-FP_N;     %subtract false positives
TP_C = NoS-FP_C;     %subtract false positives
%%%
Prec_N = double(TP_N) ./double(TP_N + FP_N);
Sens_N = double(TP_N) ./double(TP_N + FN_N);
Prec_C = double(TP_C) ./double(TP_C + FP_C);
Sens_C = double(TP_C) ./double(TP_C + FN_C);
F_C = 2*(Prec_C*Sens_C./(Prec_C+Sens_C));
F_N = 2*(Prec_N*Sens_C./(Prec_N+Sens_N));
title(['\fontsize{12} SegCells: ' num2str(NoS) ', FP_C: ' num2str(FP_C) ', FN_C: ' num2str(FN_C) ', Precision_C: ' num2str(Prec_C) ', Sensitivity_C: ' num2str(Sens_C) ', F-score: ' num2str(F_C) ', ' NN]);
saveas(gcf,[NN '.eps'],'epsc')
save([NN '.mat'],'TP_C','TP_N','FP_N','FN_N','Prec_N','Sens_N','F_N','F_C','FP_C','FN_C','Prec_C','Sens_C','TRANS3','DAPI3','SEG4','NUC4','RGB','CellInfo','NN');
%% Store variables for all images
F_C_all(1,filenum) = F_C;
F_N_all(1,filenum) = F_N;
FN_C_all(1,filenum) = FN_C;
FN_N_all(1,filenum) = FN_N;
FP_C_all(1,filenum) = FP_C;
FP_N_all(1,filenum) = FP_N;
Prec_C_all(1,filenum) = Prec_C;
Prec_N_all(1,filenum) = Prec_N;
Sens_C_all(1,filenum) = Sens_C;
Sens_N_all(1,filenum) = Sens_N;
TP_C_all(1,filenum) = TP_C;
TP_N_all(1,filenum) = TP_N;
% elseif (button == 49)
%     text(x(j),y(j),num2str(j),'Color','white','FontSize',30) ; %num2str(n)

