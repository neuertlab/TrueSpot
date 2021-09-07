function [Label_low,Label_mid,Label_hi,nuclei,Nuc_int,Nuc_vol,Maj_axis,Min_axis,dapi_label] = B1_autosegment_nuclei5(DAPI_ims,dapi_label,Yim,dapi_threshold,min_nucleus_size,max_nucleus_size,dxy1,dapi_label_low1);
% Yim = 1;
% min_nucleus_size = 50;                                                     % minimum nuclear size
% max_nucleus_size = 10000;                                                   % maximum nuclear size

a = size(DAPI_ims,1);                                                       % image size in number of pixels
z = size(DAPI_ims,3);                                                       % stack size in number of images

dxy3 = 10; %size of box to find max
i = 110;
%dxy1 = 120;                                                                  % 7 ; define number of pixels from the max DAPI pixel' CHANGED TO 100 from 11 on 1/9/2017
% dxy1 = round(dxy1*1.3);    %Use a large circle
dxy2 = round(dxy1*1.3);
dxy = dxy1;
dapi_cell2 = NaN(1+2*dxy2,1+2*dxy2,z);                                      % generate NaN 3D matrix with 2*dxy in x and y direction
Label_low = zeros(a,a,z);                                                   % generate zero 3D matrix of the size of the image stack
Label_mid = zeros(a,a,z);
Label_hi = zeros(a,a,z);
Label_index = zeros(a,a);
%figure(802); clf; imshow(dapi_label_added,[])
% figure(802); clf; imshow(DAPI_ims_added,[])
DAPI_ims_max = max(DAPI_ims(:,:,13:end),[],3);
m22 = max(dapi_label_low1(:))                                                   % determine maximum number of DAPI stained nuclei
% figure(801); clf; imshow(dapi_label,[])
Nuc_int = zeros(3,1);                                                     %Will store the integrated intensity of each nucleus
Maj_axis = zeros(3,1);                                                     %Will store the length of the major axis of each nucleus
Min_axis = zeros(3,1);                                                     %Will store the length of the minor axis of each nucleus
Nuc_vol = zeros(3,1);                                                     %Will store the volume of each nucleus
Maj_axis_MaxInt = zeros(3,1);                                             %Will store the major axis of the maximum intensity projection excluding the top and bottom 5 slices
Min_axis_MaxInt = zeros(3,1);                                              %Will store the minor axis of the maximum intensity projection excluding the top and bottom 5 slices
% go through all the DAPI nucleus labels one by one
bar1 = waitbar(1/m22,'Segmenting Nuclei')
if Yim
figure(101); clf
figure(103); clf
end
%% Thresholding by individual cells
for i = 1:m22;                                                              
    waitbar(i/m22,bar1,'Segmenting Nuclei')
    i
    k1 = dapi_label_low1 == i;                                                  % Single cell per image
%     figure(805); clf; imshow(k1,[])
    %    k2=uint16(k1);                                                          % convert to 16bit image
    k2a = immultiply(k1,DAPI_ims_max);
     %% Look at distribution of pixels inside originally determined nuclei
%     [h.Values,h.BinEdges] = histcounts(k2a(k2a>0),100);
%     ysss = h.Values;
%     xsss = zeros(size(ysss));
%     for j = 1:size(ysss,2)
%     xsss(j) = (h.BinEdges(j)+h.BinEdges(j+1))/2;
%     end
%     if Yim
%     figure(101); scatter(xsss,ysss); %hold on                                     %plot cumulative distribution function
%     saveas(gcf,['DAPI_CDF_' num2str(i) '.fig']);
%     %     figure(102); clf ; plot(xssss,yssss)
%     end   
%     %%%Smooth
%     num_maxs = 100; %Initial value for loop
% window_num = size(xsss,2)/2; %Number of windows for smoothing purposes
% fail_count = 0;
% while num_maxs > 10
% window_num = window_num-2;
% if window_num < 2
%     window_num = window_num+3;
%     num_maxs = 6;
%     fail_count = fail_count+1;
%     if fail_count>5
%         break
%     else
%     continue
%     end
% end
% % cdf_slopes1 = smooth(cdf_slopes1,'sgolay',1);
% cdf_slopes1 = smooth(xsss,ysss,1/window_num);
% %%%Alternate way of finding maxima
% DataInv = (cdf_slopes1);
% [Maxima,MaxIdx] = findpeaks(DataInv);
% % MinIdx
% %%%
% num_maxs = size(MaxIdx,1);
% if num_maxs < 1
%     window_num = window_num+3;
%     num_maxs = 6;
%     fail_count = fail_count+1;
%     if fail_count>5
%         break
%     else
%     continue
%     end
% end
% if Yim
% maxs_scaty1 = cdf_slopes1(MaxIdx);
% maxs_scatx1 = xsss(MaxIdx);
% figure(103); clf; plot(xsss',cdf_slopes1'); hold on; 
% scatter(maxs_scatx1,maxs_scaty1,'b')
% title(['windows ' num2str(window_num) ' number of maxes '  num2str(num_maxs)]);
% pause(.3)
% end
% end
% %%% Fit to log normal
% fun2 = @(x)a1*exp(-((x-b1)/c1)^2);
% xsss1 = log(xsss);
% % f = fit(xsss1',cdf_slopes1,'gauss2');
% g = fit(xsss1',cdf_slopes1,'gauss1','Upper',[10^10 10^10 10^10],'Lower',[0 min(xsss1) 0]);
% mm5 = exp(g.b1-2.25*g.c1); 
% mm4 = .9*mm5;
% mm6 = 1.1*mm5;
% %%%Using the integral
% % fun3 = @(x)g.a1*exp(-((x-g.b1)/g.c1)^2);
% % tests1 = 0:(max(xsss1)-0)/100:max(xsss1);
% % blah = integral(fun3,0,min(xsss1),tests1(50))
% %%%
% % lsqcurvefit(fun2,[f.a1 f.b1 f.c1],xsss1,cdf_slopes1) %Does not work for some reason at the moment
% if Yim
% figure(107); clf; plot(g,xsss1,cdf_slopes1); hold on
% end
% % if f.a1> f.a2
% %     f = fit(xsss1',cdf_slopes1,'gauss2','Lower', [0 0 0 0 min([f.b1+1.5*f.c1 max(xsss1(:))]) f.c1/2],'Upper',[Inf max(xsss1(:)) max(xsss1(:)) Inf max(xsss1(:)) max(xsss1(:))]) %f.c1/2
% % else
% %     f = fit(xsss1',cdf_slopes1,'gauss2','Lower', [0 0 0 0 min([f.b2+1.5*f.c2 max(xsss1(:))]) f.c2/2],'Upper',[Inf max(xsss1(:)) max(xsss1(:)) Inf max(xsss1(:)) max(xsss1(:))]) %f.c2/2
% % end
% % fun1 = @(x) a1/(x*c1*sqrt(2*pi))*exp(-((log(x)-b1)^2/2*c1^2)^2) + a2*exp(-((x-b2)/c2)^2); %Add lognormal to gaussian
% % ft = fittype('a1/(x*c1*sqrt(2*pi))*exp(-((log(x)-b1)^2/2*c1^2)^2) + a2*exp(-((x-b2)/c2)^2)')
% % f = fit(xsss',cdf_slopes1,ft)
% 
% % f = fit(xsss',cdf_slopes1,fun1,'Lower', [0 0 0 0 f.b2+1*f.c2 f.c2/2])
% % max_th = max([f.b1+2*f.c1 f.b2+2*f.c2]);
% % if integrate(f,max_th,0)/integrate(f,max(xsss(:)),0) > .95    %Look at integral of function
% %     f = fit(xsss',cdf_slopes1,'gauss2')
% % end 
% if Yim
%     figure(107); clf; plot(g,xsss1,cdf_slopes1); hold on
% end
% %%%assume threshold is the Gaussian with the highest amplitude
% % if f.a1> f.a2
% % mm_temp = exp(f.b1-1.67*f.c1)     %
% % else
% %   mm_temp = exp(f.b2-1.67*f.c2)     %  
% % end
% %   mm4 = mm5*.9;
% % mm6 = mm5 * 1.1;
% % scatter(log(mm5),f(log(mm5)));
% pause(.1)
    %%
    m111 = max(k2a(:));
    sumposx = 0;                                                             %will be numerator for calculating weighted average for x
    sumposy = 0;                                                            %will be numerator for calculating weighted average for y
    sumfluor = 0;                                                            %will be denominator for calculating weighted average
    %k2atest = zeros(1,1024*1024);
    counter = 1;
    [xs,ys] = ind2sub(size(k1),find(k1 == 1));
    sumposx = sum(xs(:)); sumposy = sum(ys(:));
    sumfluor = sum(k1(:));
    %%%Old slow way of finding weighted average
%     for j = 1:size(k2a,1)
%         for k = 1:size(k2a,2)
%             %k2atest(counter) = k1(j,k);
%             counter = counter+1;
%             sumposx = sumposx + (k1(j,k)*j);
%             sumposy = sumposy + (k1(j,k)*k);
%             sumfluor = sumfluor + k1(j,k);
%         end
%     end
    %%%
    xA = round(sumposx/sumfluor);
    yA = round(sumposy/sumfluor);

    %[xA,yA] = find(k2a == m111);
%     k3=regionprops(k2,'Centroid');                                          % center of maximun in 3D called centroid
%     k4 = k3.Centroid;                                                      
%     yA = round(k4(1,1));                                                    % round yA to full pixel
%     xA = round(k4(1,2));                                                    % round xA to full pixel
      clear k1 k2 k2a sumposx sumposy  sumfluor xs ys
    Imzero = uint16(zeros(a,a));
    Imzero(xA,yA) = 1;
%     area1 = zeros(a,a,z);                                                   % generate zero 3D matrix of the size of the image stack
    A = xA-dxy2;                                                            % generate corner points with are dxy2 pixels away from the nuclear center
    B = xA+dxy2;
    C = yA-dxy2;
    D = yA+dxy2;
    A1 = xA-dxy3;                                                            % generate corner points for determining max which are dxy3 pixels away from the nuclear center
    B1 = xA+dxy3;
    C1 = yA-dxy3;
    D1 = yA+dxy3;
    stop1 = 0;
    for j = [1,numel(A)];                                                   % if corner A is below 1 pixel or negative set equals 1
        if A(j) < 1;
            A(j) = 1;
            stop1 = 1;
        end
    end
    for j = [1,numel(B)];                                                   % if corner B is above a pixels (1024 or 2048) set equals a
        if B(j) > a;
            B(j) = a;
            stop1 = 1;
        end
    end
    for j = [1,numel(C)];                                                   % if corner C is below 1 pixel or negative set equals 1
        if C(j) < 1;
            C(j) = 1;
            stop1 = 1;
        end
     end
    for j = [1,numel(D)];                                                   % if corner D is above a pixels (1024 or 2048) set equals a
        if D(j) > a;
            D(j) = a;
            stop1 = 1;
        end
    end
    if stop1 == 1
%         continue
    end
    if  A1 < 1;                                                              % if corner A is below 1 pixel or negative set equals 1
        A1 = 1;
%         continue
    else;
    end;
    if B1 > a;                                                               % if corner B is above a pixels (1024 or 2048) set equals a
%         continue
        B1 = a;
    else;
    end;
    if C1 < 1;                                                               % if corner C is below 1 pixel or negative set equals 1
%         continue
        C1 = 1;
    else;
    end;
    if D1 > a;                                                                % if corner D is above a pixels (1024 or 2048) set equals a
%         continue
        D1 = a;
    else;
    end;
    
    %%%Fast way of generating circle
%     figure(6); imshow(max(area1,[],3),[])    
    area1 = zeros(a,a,1);
%     tic
    for j = 1:a                                                                 %This sets the area by making a circle around the center
        for k = 1:a
            if sqrt((j-xA)^2+(k-yA)^2) <= dxy2
                area1(j,k,1) = 1;
            end
        end
    end
%     toc
%     figure(6); imshow(max(area1,[],3),[]) 
    %%%
    %%% Slower way of generating circle
%     area1 = zeros(a,a);
%     tic
%     xs = zeros([a,a]);
%     ys = zeros([a,a]);
%     for j = 1:a                                                                 %This sets the area by making a circle around the center
%         xs(j,:) = j;
%     end
%     for k = 1:a  
%         ys(:,k) = k;
%     end
%     dists_sq = (xs-xA).^2+(ys-yA).^2;
%     area1 = dists_sq <= dxy2^2;
%     area1 = repmat(area1,[1 1 z]);
%     toc
%     figure(6); imshow(max(area1,[],3),[]) 
    %%%
%     area1(A:B,C:D,:) = 1;                                                   % set pixels within A,B,C,D equals 1
%     figure(10); clf; imshow(max(area1,[],3),[]); impixelinfo;
    %% Determine how much area inside circle (area1) is taken up by the dapi_label_low1 image   
    total_area = sum(area1(:))/size(area1,3);       %Total area inside area1
    nuc_area = immultiply(area1(:,:,1),dapi_label_low1==i); 
    if Yim 
        nuc_area1 = nuc_area;
    end
    nuc_area = sum(nuc_area(:));
    cumul_cutoff = 1-(nuc_area/total_area); %Tentative cutoff for the cumulative distribution
    %%
 dapi_cell = immultiply(uint16(DAPI_ims),uint16(repmat(area1,[1,1,z])));                         % generate dapi image stack inside area1 and set the rest of the image to zero
other_nuc =  immultiply(dapi_label_low1 ~= i,dapi_label_low1 >0);     %  Find where other nuclei are
other_nuc =  immultiply(other_nuc,area1(:,:,1));        %find where other nuclei are in area
 dapi_cell_max = immultiply(DAPI_ims_max,area1(:,:,1));                                    % determine maximum projection
 temp1 = dapi_cell_max(dapi_cell_max>0);
 dapi_cell(repmat(other_nuc,[1,1,z])==1) =  min(temp1(:)); 
%  dapi_cell(repmat(other_nuc,[1,1,z])==1) =  min(dapi_cell_max(dapi_cell_max>0),[],'all'); 
temp1 = dapi_cell_max(dapi_cell_max>0);
% dapi_cell_max(other_nuc==1) = min(dapi_cell_max(dapi_cell_max>0),[],'all');                                   %Set the intensity of the other nuclei to the min  
dapi_cell_max(other_nuc==1) = min(temp1(:));                                   %Set the intensity of the other nuclei to the min      
if Yim
        Nuc0 = dapi_cell_max(A:B,C:D);
         Nuc1 = Nuc0/max(Nuc0(:));
        Nuc2 = imadjust(Nuc1,[min(Nuc1(:)) 1]); 
        NucBorder1 = bwmorph(nuc_area1(A:B,C:D),'remove')*100000;
         NucBorder1 = bwmorph(NucBorder1,'thicken',6);
        se = strel('disk',50);
        NucBorder1 = imclose(NucBorder1,se);
        AreaBorder1 = bwmorph(area1(A:B,C:D),'remove')*100000;
%         AreaBorder1 = bwmorph(AreaBorder1,'thicken',4);
        se = strel('disk',50);
        AreaBorder1 = imclose(AreaBorder1,se);
        R = Nuc2; G = Nuc2+AreaBorder1; Bl = Nuc2+NucBorder1; %dark green rgb(34,139,34)
         RGB = cat(3,R,G,Bl);
        figure(12); clf;
        imshow(RGB); impixelinfo;
        saveas(gcf,['DAPI_Max_circle' num2str(i) '.fig']);
    else;
    end;
      dapi_cell2 = dapi_cell(A:B,C:D,:);
    dapi_cell2max = double(dapi_cell_max(A:B,C:D));                                 % cut out the maxiumum projection
%     figure(12); clf; imshow(dapi_cell2max,[]); impixelinfo;
    %     dapi_cell3max = double(dapi_cell_max(A1:B1,C1:D1));                                 % cut out the maxiumum projection (smaller for max)
    dapi_cell2max(dapi_cell2max == 0) = NaN;                                    %This makes all zero elements NaN instead
%     figure(13); clf; imshow(max(dapi_cell2,[],3),[]); impixelinfo;
%     figure(14); clf; imshow(dapi_cell2max,[]); impixelinfo;
%     D1max = dapi_cell2max(:);
%     counter74 = 1;
%     for j = nanmin(dapi_cell2max(:)):10:nanmax(dapi_cell2max(:))            %This does a cumulative distribution function for total fluorescence
%         xssss(counter74) = j;
%         yssss(counter74) = nansum(D1max(find(D1max < j)))/nansum(D1max);
%         counter74 = counter74+1;
%     end
        %% Determine thresholds based on cutoffs in cdf
    [ysss,xsss] = ecdf(dapi_cell2max(:));                                       %determine cumulative distribution 
%     zsss = abs(ysss-(cumul_cutoff-.1));                                                    %Find values closest to certain threshold in cumulative distribution function
    z2sss = abs(ysss-cumul_cutoff); 
%     z3sss = abs(ysss-(cumul_cutoff+.1)); 
%     mm4 = xsss(find(zsss == min(zsss),1,'first'));                          %threshold at cutoff for cumulative distribution function
    mm5 = xsss(find(z2sss == min(z2sss),1,'first'));                          %threshold at cutoff for cumulative distribution function
%     mm6 = xsss(find(z3sss == min(z3sss),1,'first'));                          %threshold at cutoff for cumulative distribution function
    mm4 = mm5*.9;
    mm6 = mm5*1.1;   
    if Yim
%         scatter([mm4 mm5 mm6],[cumul_cutoff-.1 cumul_cutoff cumul_cutoff+.1]) %Cutoffs when based on different parts of cumulative distribution
%         scatter([mm5], [cumul_cutoff])
    end
    if Yim
        figure(201); clf; plot(xsss,ysss,'k','LineWidth',2); hold on
        first_xs = min(xsss(:)):10:mm5;
        plot(first_xs,repmat(cumul_cutoff,[1,size(first_xs,2)]),'b','LineWidth',2);
        first_ys = 0:.01:cumul_cutoff;
        plot(repmat(mm5,[1 size(first_ys,2)]),first_ys,'b','LineWidth',2);
        plot(xsss,repmat(1,[size(xsss,1),1]),'g','LineWidth',2);
        xlabel('DAPI Intensity')
        ylabel('Cumulative Probability')
        ylim([0 1.1])
        xlim([min(xsss(:)) max(xsss(:))])
    end
    %% Histogram of entire nucleus
%     h = struct;
%     [h.Values,h.BinEdges] = histcounts(dapi_cell2max(dapi_cell2max>0),100);
%     ysss = h.Values;
%     xsss = zeros(size(ysss));
%     for j = 1:size(ysss,2)
%     xsss(j) = (h.BinEdges(j)+h.BinEdges(j+1))/2;
%     end
%     if Yim
%     figure(101); scatter(xsss,ysss); %hold on                                     %plot cumulative distribution function
%     saveas(gcf,['DAPI_CDF_' num2str(i) '.fig']);
%     %     figure(102); clf ; plot(xssss,yssss)
%     end
% %% Determine thresold based on slope of cumulative distribution function (finding first local minimum)
% num_mins = 100; %Initial value for loop
% window_num = size(xsss,2)/2; %Number of windows for smoothing purposes
% fail_count = 0;
% while num_mins > 10
% window_num = window_num-2;
% if window_num == 0
%     window_num = window_num+3;
%     num_mins = 100;
%     fail_count = fail_count+1;
%     if fail_count>5
%         break
%     else
%     continue
%     end
% end
% % % %%% using the ecdf
% % % % window_rad1 = round(size(ysss,1)/window_num);                               %radius of the window for determining the slope of the cdf
% % % % cdf_slopes1 = zeros(size(ysss,1),1);                                        %Will have the slope of the cdf at each point
% % % % for j = (window_rad1+1):(size(ysss,1)-(1+window_rad1))
% % % %     cdf_slopes1(j,1) = (ysss(j+window_rad1)-ysss(j-window_rad1))/window_rad1;
% % % % end
% % % % cdf_slopes1 = cdf_slopes1(round(size(cdf_slopes1,1)/50):end);  %cut off the very beginning
% % % % xsss1 = xsss(round(size(xsss,1)/50):end);  %cut off the very beginning
% %%% Using histogram
% if Yim
% figure(101); clf
% end
% %%%
% % cdf_slopes1 = smooth(cdf_slopes1,'sgolay',1);
% cdf_slopes1 = smooth(xsss,ysss,1/window_num);
% [mins1,Ps1] = islocalmin(cdf_slopes1);
% %%%Alternate way of finding minima
% DataInv = 1.01*max(cdf_slopes1) - cdf_slopes1;
% [Minima,MinIdx] = findpeaks(DataInv);
% % MinIdx
% %%%
% num_mins = size(MinIdx,1);
% if num_mins == 0
%     window_num = window_num+3;
%     num_mins = 100;
%     fail_count = fail_count+1;
%     if fail_count>5
%         break
%     else
%     continue
%     end
% end
% if Yim
% mins_scaty = cdf_slopes1(mins1);
% mins_scatx = xsss(mins1);
% mins_scaty1 = cdf_slopes1(MinIdx);
% mins_scatx1 = xsss(MinIdx);
% figure(103); clf; plot(xsss',cdf_slopes1'); hold on; 
% scatter(mins_scatx,mins_scaty)
% scatter(mins_scatx1,mins_scaty1,'b')
% title(['windows ' num2str(window_num) ' number of mins '  num2str(num_mins)]);
% pause(.3)
% end
% end
%%
% % cdf_slopes2 = zeros(size(ysss,1),1);                                        %Will have the 2nd derivative of the cdf at each point
% % for j = (window_rad1+1):(size(ysss,1)-(1+window_rad1))
% %     cdf_slopes2(j,1) = (cdf_slopes1(j+window_rad1)-cdf_slopes1(j-window_rad1))/window_rad1;
% % end
% % cdf_slopes2(window_rad1:(size(ysss,1)-window_rad1+1)) = ((ysss(window_rad1+1:end)-ysss(1:(size(ysss,1)-window_rad1)))/window_rad1);
% % if Yim
% % figure(103); clf; plot(xsss',cdf_slopes1'); hold on
% % % figure(104); clf; plot(xsss',cdf_slopes2')
% % end
% %     mm4 = min(xsss(mins1))*.9;           %threshold at cutoff for cumulative distribution function
% %     mm5 = min(xsss(mins1));             %threshold at cutoff for cumulative distribution function
% %     mm6 = min(xsss(mins1))*1.1;           %threshold at cutoff for cumulative distribution function
% %% Determine threshold based on fitting two Gaussians 
% % cdf_slopes1 = cdf_slopes1/sum(cdf_slopes1(:));  %Normalize
% % cdf_slopes1 = ysss; %Do this if you do not want smoothing
% xsss1 = log(xsss);
% f = fit(xsss1',cdf_slopes1,'gauss2','Lower', [0 0 0 0 g.b1 g.c1/2],'Upper',[10^10 g.b1-g.c1/3 10^10 10^10 g.b1 g.c1*2]);
% if Yim
% figure(107); clf; plot(f,xsss1,cdf_slopes1); hold on
% end
% % if f.a1> f.a2
% %     f = fit(xsss1',cdf_slopes1,'gauss2','Lower', [0 0 0 0 min([f.b1+1.5*f.c1 max(xsss1(:))]) f.c1/2],'Upper',[Inf max(xsss1(:)) max(xsss1(:)) Inf max(xsss1(:)) max(xsss1(:))]) %f.c1/2
% % else
% %     f = fit(xsss1',cdf_slopes1,'gauss2','Lower', [0 0 0 0 min([f.b2+1.5*f.c2 max(xsss1(:))]) f.c2/2],'Upper',[Inf max(xsss1(:)) max(xsss1(:)) Inf max(xsss1(:)) max(xsss1(:))]) %f.c2/2
% % end
% % fun1 = @(x) a1/(x*c1*sqrt(2*pi))*exp(-((log(x)-b1)^2/2*c1^2)^2) + a2*exp(-((x-b2)/c2)^2); %Add lognormal to gaussian
% % ft = fittype('a1/(x*c1*sqrt(2*pi))*exp(-((log(x)-b1)^2/2*c1^2)^2) + a2*exp(-((x-b2)/c2)^2)')
% % f = fit(xsss',cdf_slopes1,ft)
% 
% % f = fit(xsss',cdf_slopes1,fun1,'Lower', [0 0 0 0 f.b2+1*f.c2 f.c2/2])
% % max_th = max([f.b1+2*f.c1 f.b2+2*f.c2]);
% % if integrate(f,max_th,0)/integrate(f,max(xsss(:)),0) > .95    %Look at integral of function
% %     f = fit(xsss',cdf_slopes1,'gauss2')
% % end 
% if Yim
%     figure(107); clf; plot(f,xsss1,cdf_slopes1); hold on
% end
% %%%assume threshold is the Gaussian with the highest amplitude
% if f.b1< f.b2
% mm_temp = exp(f.b1+1.28*f.c1)     %
% else
%   mm_temp = exp(f.b2+1.28*f.c2)     %  
% end
% % mm5 = mean([mm5 mm_temp]);
% % mm5 = max([mm5 exp(min([f.b1 f.b2]))]);  %Make threshold at least mean of background
% %   mm4 = mm5*.9;
% % mm6 = mm5 * 1.1;
% if Yim
% scatter(log(mm5),f(log(mm5)));
% end
% pause(.1)
% %%% Determine thresholds based on arbitrary cutoffs in cdf
% %     zsss = abs(ysss-.2);                                                    %Find values closest to certain threshold in cumulative distribution function
% %     z2sss = abs(ysss-.3);                                                    %Find values closest to certain threshold in cumulative distribution function
% %     z3sss = abs(ysss-.4);                                                    %Find values closest to certain threshold in cumulative distribution function
% %     mm4 = xsss(find(zsss == min(zsss),1,'first'));                          %threshold at cutoff for cumulative distribution function
% %     mm5 = xsss(find(z2sss == min(z2sss),1,'first'));                          %threshold at cutoff for cumulative distribution function
% %     mm6 = xsss(find(z3sss == min(z3sss),1,'first'));                          %threshold at cutoff for cumulative distribution function
% %%%
%     %     mm1 = nanmin(dapi_cell2max(:));                                            % determine minimum dapi signal
% %     mm2 = nanmax(dapi_cell3max(:));                                            % determine maximum dapi signal
% %     mm3 = nanmedian(dapi_cell2max(:));                                            % determine median dapi signal
% %     mmd = mm2-mm1;                                                          % determine the difference between maximum and minimum
% %     m1 = mm1+mmd*0.4;                                                       % set "LOW" threshold to 40% of difference
% %     m2 = mm1+mmd*0.5;                                                       % set "MID" threshold to 50% of difference
% %     m3 = mm1+mmd*0.6;                                                       % set "HI" threshold to 60% of difference
% %     m1 = mm3*0.9;                                                       % set "LOW" threshold to 90% of median
% %     m2 = mm3*1.0;                                                       % set "MID" threshold to 100% of median
% %     m3 = mm3*1.1;                                                       % set "HI" threshold to 110% of median
    m1 = mm4;                                                       % set "LOW" threshold to first cumulative distribution cutoff
    m2 = mm5;                                                       % set "MID" threshold to second of cumulative distribution cutoff
    m3 = mm6;                                                       % set "HI" threshold to third of cumulative distribution cutoff
    dapi_label3A = dapi_cell2 > (m1);                                        % "LOW" generate 3D binary image above the thresholds
    dapi_label3B = dapi_cell2 > (m2);                                        % "MID" generate 3D binary image above the thresholds
    dapi_label3C = dapi_cell2 > (m3);                                        % "HI" generate 3D binary image above the thresholds 
% if Yim
%      figure(902); clf; imshow(max(dapi_label3B,[],3),[])
% end
%% This fills holes in the image in 3D, but it is very computationally intensive
%     dapi_label3A = imfill(dapi_label3A,'holes');                                        % "LOW" fill holes in image
%     dapi_label3B = imfill(dapi_label3B,'holes');                                        % "MID" fill holes in image
%     dapi_label3C = imfill(dapi_label3C,'holes');                                        % "HI" fill holes in image
%%%
%%% Below fills in holes and takes the most connected area in 2D
    for k = 1:z
    dapi_label3A(:,:,k) = imfill(dapi_label3A(:,:,k),'holes');                                        % "LOW" fill holes in image
    dapi_label3B(:,:,k) = imfill(dapi_label3B(:,:,k),'holes');                                        % "MID" fill holes in image
    dapi_label3C(:,:,k) = imfill(dapi_label3C(:,:,k),'holes');                                        % "HI" fill holes in image 
   % figure(14); clf; imshow(dapi_label3A(:,:,k),[])
%     dapi_label_low = bwlabeln(dapi_label3A(:,:,k),8);                                          % "LOW" label different connected neighborhoods
%     dapi_label_mid = bwlabeln(dapi_label3B(:,:,k),8);                                          % "MID" label different connected neighborhoods
%     dapi_label_hi = bwlabeln(dapi_label3C(:,:,k),8);                                          % "HI" label different connected neighborhoods
%     %figure(13); clf; imshow(dapi_label_low,[])
%     dapi_label3A(:,:,k) = dapi_label_low == mode(dapi_label_low(dapi_label_low>0));                                        % "LOW" fill holes in image
%     dapi_label3B(:,:,k) = dapi_label_mid == mode(dapi_label_mid(dapi_label_mid>0));                                        % "MID" fill holes in image
%     dapi_label3C(:,:,k) = dapi_label_hi == mode(dapi_label_hi(dapi_label_hi>0));                                        % "HI" fill holes in image 
    end
    %Find 3D connected area in 3D
    dapi_label_low = bwlabeln(dapi_label3A,6);                                          % "LOW" label different connected neighborhoods
    dapi_label_mid = bwlabeln(dapi_label3B,6);                                          % "MID" label different connected neighborhoods
    dapi_label_hi = bwlabeln(dapi_label3C,6);                                          % "HI" label different connected neighborhoods
    %figure(13); clf; imshow(dapi_label_low,[])
    dapi_label3A = dapi_label_low == mode(dapi_label_low(dapi_label_low>0));                                        % "LOW" fill holes in image
    dapi_label3B = dapi_label_mid == mode(dapi_label_mid(dapi_label_mid>0));                                        % "MID" fill holes in image
    dapi_label3C = dapi_label_hi == mode(dapi_label_hi(dapi_label_hi>0)); 
    if Yim
     figure(903); clf; imshow(max(dapi_label3B,[],3),[])
    end
        %% Checking each slice for correct segmentation (comment out unless testing)
% % 
% for j = 1:z %Go through each stack
%     figure(1101); clf; imshow(dapi_label3B(:,:,j),[]);
%     Mid_DAPI = dapi_cell2(:,:,j);
%          Nuc1 = Mid_DAPI/max(dapi_cell2(:));
%             Nuc2 = imadjust(Nuc1); %,[0 1]
%                 NucBorder1 = bwmorph(dapi_label3B(:,:,j),'remove')*100000;
%                 NucBorder1 = bwmorph(NucBorder1,'thicken',2);
%                 se = strel('disk',3);
%                 NucBorder1 = imclose(NucBorder1,se);
%             R = NucBorder1; B = Nuc2; G = NucBorder1;
%             RGB = cat(3,R,G,B);
%             figure(1100); clf;
%             imshow(RGB,[]); 
%             title(['slice ' num2str(j)])
% %             pause(.1)
% end
% %     min3 = 200;
% %     max3 = 300;
% %     [xs3,ys3,zs3] = ind2sub(size(dapi_cell2),intersect(find(dapi_cell2>min3),find(dapi_cell2<max3)));
% %     figure(106); clf;
% %     scatter3(xs3,ys3,zs3)
    %% Eliminate Cells at the border
    temp_props = regionprops3(dapi_label3A,'BoundingBox') ;              %Obtain axis lengths for nucleus with low threshold
%     figure;volshow(dapi_label3A)
    k4 = temp_props.BoundingBox;                                                    %create the rectangular box around the cell 
    stop1 = 0;
    for shapes1 = 1:size(k4,1)  %Go through each shape
    X0=round(k4(shapes1,1));
    Y0=round(k4(shapes1,2));
    X1=round(k4(shapes1,1)+ k4(shapes1,4))-1;
    Y1=round(k4(shapes1,2)+ k4(shapes1,5))-1;
    if (X0+C <4)|(X1+C>a-4)|(Y0+A <4)|(Y1+A>a-4)
%         Label1 = Label1 - uint16(k11)*j;
        stop1 = 1;
    end
%     Label_low(X0+A:X1+A,[Y0 Y1]+C,:) = 1;
%     Label_mid(X0+A:X1+A,[Y0 Y1]+C,:) = 1;
%     Label_hi(X0+A:X1+A,[Y0 Y1]+C,:) = 1;
%     Label_low([X0 X1]+A,Y0+C:Y1+C,:) = 1;
%     Label_mid([X0 X1]+A,Y0+C:Y1+C,:) = 1;
%     Label_hi([X0 X1]+A,Y0+C:Y1+C,:) = 1;
    end
    if stop1 == 1
%         [num2str(i) ' was at the edge']
dapi_label(dapi_label_low1 == i) = 0; 
        continue
    end
%% Store aspects of nucleus
Nuc_int(1,i) = sum(sum(sum(immultiply(dapi_label3A,dapi_cell2))));          %Store integrated intensity of nucleus with lower threshold
Nuc_int(2,i) = sum(sum(sum(immultiply(dapi_label3B,dapi_cell2))));          %Store integrated intensity of nucleus with middle threshold
Nuc_int(3,i) = sum(sum(sum(immultiply(dapi_label3C,dapi_cell2))));          %Store integrated intensity of nucleus with lower threshold
Nuc_vol(1,i) = sum(dapi_label3A(:));                                        %Store volume of nucleus with low threshold
Nuc_vol(2,i) = sum(dapi_label3B(:));                                        %Store volume of nucleus with medium threshold
Nuc_vol(3,i) = sum(dapi_label3C(:));                                        %Store volume of nucleus with high threshold
% temp_props = regionprops3(dapi_label3A,'PrincipalAxisLength') ;              %Obtain axis lengths for nucleus with low threshold
% try Maj_axis(1,i) = temp_props.PrincipalAxisLength(1);  catch   Maj_axis(1,i) = NaN; end       %Store the Principal Axis length for nucleus with low threshold
% try Min_axis(1,i) = temp_props.PrincipalAxisLength(2);  catch  Min_axis(1,i) = NaN; end       %Store the Principal Axis length for nucleus with low threshold
% temp_props = regionprops3(dapi_label3B,'PrincipalAxisLength') ;              %Obtain axis lengths for nucleus with medium threshold
% try Maj_axis(2,i) = temp_props.PrincipalAxisLength(1); catch Maj_axis(2,i) = NaN; end          %Store the Principal Axis length for nucleus with medium threshold
% try Min_axis(2,i) = temp_props.PrincipalAxisLength(2); catch Min_axis(2,i) =  NaN; end        %Store the Principal Axis length for nucleus with medium threshold
% temp_props = regionprops3(dapi_label3C,'PrincipalAxisLength')  ;             %Obtain axis lengths for nucleus with high threshold
% try Maj_axis(3,i) = temp_props.PrincipalAxisLength(1); catch  Maj_axis(3,i) = NaN; end        %Store the Principal Axis length for nucleus with high threshold
% try Min_axis(3,i) = temp_props.PrincipalAxisLength(2); catch  Min_axis(3,i) = NaN; end         %Store the Principal Axis length for nucleus with high thresholdtemp_props = regionprops3(dapi_label3C,"PrincipalAxisLength")  ;
% temp_props = regionprops(max(dapi_label3A(:,:,6:size(dapi_label3A,3)-6),[],3),'MajorAxisLength','MinorAxisLength');
% try Maj_axis_MaxInt(1,i) = temp_props.MajorAxisLength; catch Maj_axis_MaxInt(1,i) = NaN; end;
% try Min_axis_MaxInt(1,i) = temp_props.MinorAxisLength; catch Min_axis_MaxInt(1,i) = NaN;  end
% temp_props = regionprops(max(dapi_label3B(:,:,6:size(dapi_label3A,3)-6),[],3),'MajorAxisLength','MinorAxisLength');
% try Maj_axis_MaxInt(2,i) = temp_props.MajorAxisLength; catch Maj_axis_MaxInt(2,i) = NaN;  end
% try Min_axis_MaxInt(2,i) = temp_props.MinorAxisLength; catch Min_axis_MaxInt(2,i) = NaN;  end
% temp_props = regionprops(max(dapi_label3C(:,:,6:size(dapi_label3A,3)-6),[],3),'MajorAxisLength','MinorAxisLength');
% try Maj_axis_MaxInt(3,i) = temp_props.MajorAxisLength; catch Maj_axis_MaxInt(3,i) = NaN;  end
% try Min_axis_MaxInt(3,i) = temp_props.MinorAxisLength; catch Min_axis_MaxInt(3,i) = NaN;  end
%%%
Label_low(A:B,C:D,:) = Label_low(A:B,C:D,:) + dapi_label3A;                                   % Add to label matrix 'LOW'
Label_mid(A:B,C:D,:) = Label_mid(A:B,C:D,:) + dapi_label3B;                                   % Add to label matrix 'MID'
Label_hi(A:B,C:D,:) = Label_hi(A:B,C:D,:) + dapi_label3C;                                     % Add to label matrix 'HI'

%%% Below is one way to avoid extra spots in the image (, but it takes up a lot of RAM
%     dapi_label_low = bwlabeln(dapi_label3A,26);                                          % "LOW" label different connected neighborhoods
%     dapi_label_mid = bwlabeln(dapi_label3B,26);                                          % "MID" label different connected neighborhoods
%     dapi_label_hi = bwlabeln(dapi_label3C,26);                                          % "HI" label different connected neighborhoods
%     dapi_label3A = dapi_label_low;                                        % "LOW" generate 3D binary image above the thresholds
%     dapi_label3B = dapi_label_mid;                                        % "MID" generate 3D binary image above the thresholds
%     dapi_label3C = dapi_label_hi;                                        % "HI" generate 3D binary image above the thresholds  
%     dapi_label3A(dapi_label3A==0)=[];
%     dapi_label3B(dapi_label3B==0)=[];
%     dapi_label3C(dapi_label3C==0)=[];
%     Label_low = Label_low + dapi_label_low == mode(dapi_label3A);                                   % Add to label matrix 'LOW'
%     Label_mid = Label_mid + dapi_label_mid == mode(dapi_label3B);                                   % Add to label matrix 'MID'
%     Label_hi = Label_hi + dapi_label_hi == mode(dapi_label3C);                                     % Add to label matrix 'HI'
 %%%  For looking at different slices of the image 
%     for k = [10:10:70]
%     figure(15); clf; imshow(Label_low(:,:,k)); impixelinfo; title(['slice ' num2str(k)])
%     figure(16); clf; imshow(Label_mid(:,:,k)); impixelinfo; title(['slice' num2str(k)])
%     figure(17); clf; imshow(Label_hi(:,:,k)); impixelinfo; title(['slice ' num2str(k)])
%     pause(10)
%     end
%     figure(16); clf; imshow(dapi_cell2max,[]); impixelinfo;
%%%
    Label_index = uint16(Label_index) + Imzero;
end;
close(bar1)
% figure(22); clf; imshow(max(Label_mid,[],3),[]);
% for i = 1:z
%     figure(22); clf; imshow(Label_mid(:,:,i)); pause(.1)
% end
% for i = 1:z
%     
%     figure(22); clf; imshow(Label_mid(:,:,i)); pause(.1)
% end
if Yim == 1;
    figure(21); clf; imshow(max(Label_low,[],3),[]); impixelinfo; imwrite(uint16(max(Label_low,[],3)),['Label_low.tif']); 
    figure(22); clf; imshow(max(Label_mid,[],3),[]); impixelinfo; imwrite(uint16(max(Label_mid,[],3)),['Label_mid.tif']); 
    figure(23); clf; imshow(max(Label_hi,[],3),[]); impixelinfo; imwrite(uint16(max(Label_hi,[],3)),['Label_hi.tif']); 
    figure(24); clf; imshow(Label_index,[]); impixelinfo; imwrite(uint16(Label_index),['Label_index.tif']); 
else;
end;
% figure(25); clf;  volshow(Label_low); 
% figure(26); clf;  volshow(Label_mid); 
% figure(27); clf;  volshow(Label_hi); 
    Label_low = Label_low > 0;                                   % Make sure it is binary image
    Label_mid = Label_mid > 0;                                   % Make sure it is binary image
    Label_hi = Label_hi > 0;                                     % Make sure it is binary image
% for i = 1:z
%     figure(1); clf; imshow(Label_low(:,:,i),[])
%     pause(.5)
% end
% figure(2); volshow(Label_low)
size(Label_low)
    nuclei = Label_index;                                                       % max projection of the segmented dapi signals using the 50% threshold
%% Calculate aspects of the nuclei such as volume, integrated intensity, and 
