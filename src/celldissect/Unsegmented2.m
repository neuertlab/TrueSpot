function Unsegmented2(trans_plane,cells,Label_mid,DAPI,CellInfo,NN,d1,Sc,LT)

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

trans_plane = immultiply(trans_plane,SB);
TRANS2 = im2double(trans_plane./max(trans_plane(:)));
%figure(fig); clf; imshow(TRANS2,[]); title('trans'); impixelinfo; fig = fig + 1;
n = 3; m = 0;
Tm = mean2(TRANS2(:));
Ts = std2(TRANS2(:));
TRANS3 = imadjust(TRANS2,[Tm-m*Ts Tm+n*Ts],[]); %Tm-n*Ts
%TRANS3 = imadjust(TRANS2); % [0 1] add auto adjusting. 
figure(fig); clf; imshow(TRANS3); title('trans'); impixelinfo; fig = fig + 1;


%figure(fig); clf; imshow(DAPI1,[]); title('DAPI'); impixelinfo; fig = fig + 1;
DAPI2 = im2double(DAPI1); % ./max(DAPI1(:))
DAPI2 = immultiply(DAPI2,SB);
figure(fig); clf; imshow(DAPI2,[]); title('DAPI'); impixelinfo; fig = fig + 1;
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
figure(fig); clf; imshow(SEG4); title('Seg Cells'); impixelinfo; fig = fig + 1;

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
figure(fig); clf; imshow(NUC4); title('Seg Nuc'); impixelinfo; fig = fig + 1;

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
pos = find(CellInfo(:,1) > x1 & CellInfo(:,1) < x2 & CellInfo(:,2) > y1 & CellInfo(:,2) < y2);
NoS = size(pos,1);
for i = 1: size(pos,1); %No;
    i
    text(CellInfo(pos(i,1),1)-x1,CellInfo(pos(i,1),2)-y1,num2str(i),'Color','white','FontSize',18);
%    plot(CellInfo(i,1),CellInfo(i,2),[num2str(i)]); hold on;
end;
title(['\fontsize{14} Left (red) - False Positive (:  automatic segmented area but not manual  ||  Right (green) - False Negative: manual segmented area but not automatic || Q: quit' ]);

%% Mark single cells;
n = 0;
FP = 0
FN = 0;
while n < NoS;
    [x,y,button] = ginput(1);
    n = n + 1;
    if button == 1;             % Left
        FP = FP + 1;
        text(x,y,[num2str(FP) ' FP'],'Color','red','FontSize',18) ; %num2str(n)
    elseif button == 3;
        FN = FN + 1;
        text(x,y,[num2str(FN) ' FN'],'Color','green','FontSize',18) ; %num2str(n)
    elseif button == 113;
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
Prec = double(NoS) ./double(NoS + FP)
Sens = double(NoS) ./double(NoS + FN)
F = 2*(Prec*Sens./(Prec+Sens))
title(['\fontsize{16} SegCells: ' num2str(NoS) ', FP: ' num2str(FP) ', FN: ' num2str(FN) ', Precision: ' num2str(Prec) ', Sensitivity: ' num2str(Sens) ', F-score: ' num2str(F) ', ' NN]);
saveas(gcf,[NN '.eps'],'epsc')
save([NN '.mat'],'F','FP','FN','Prec','Sens','TRANS3','DAPI3','SEG4','NUC4','RGB','CellInfo','NN');


% elseif (button == 49)
%     text(x(j),y(j),num2str(j),'Color','white','FontSize',30) ; %num2str(n)

