function [PAR]...
    = B4_FitGaussians3(Lab, mm,TMR3Dorg,TMR3Dback,TMR3D3immax,TMR3Dfilter,RNAposTMR,Nuc3D)  

d = 4; %Defines the image size to be 4 pixels greater than nucleus boundary
% generate masks
b = double(ones(9,9)); 
b(2:8,2:8) = 0; 
b1= (1-b)*0.5;
b2 = b1;
b2(3:7,3:7) = 2;
b2 = b2 ./2;
clear x3y3 g1 g2 g3 XY len
zz = size(TMR3Dorg,3);
A = size(TMR3Dorg,1);
q = 47;
cnt = 5;
%RNAposTMR = RNAposCY5
% figure(104); clf; imshow(CellNuc + max(CellRNAindex2,[],3),[0 3]); hold on; plot(RNAposTMR(,yall,'o'); title(['Cell: ' num2str(j)]); pixval;
        %% Fit gaussians to mRNA spots in each cell
for i=1:mm
    len(i) = nnz(RNAposTMR(i,:,1));
end
nmbr = max(len);
%declare parameter matrix
PAR.xfit = NaN(mm,nmbr);            %x position from gaussian fit.  Row!
PAR.yfit = NaN(mm,nmbr);            %y postion from gaussian fit. Column!
PAR.xinit = NaN(mm,nmbr);           %x initial position from nongaussian fit.
PAR.yinit = NaN(mm,nmbr);           %y initial position from nongaussian fit.
PAR.xgw = NaN(mm,nmbr);             %x gaussian width
PAR.ygw = NaN(mm,nmbr);             %y gaussian width
PAR.xFWHM = NaN(mm,nmbr);           %x Full Width at Half Max of gaussian fit
PAR.yFWHM = NaN(mm,nmbr);           %y Full Width at Half Max of gaussian fit
PAR.expMInt = NaN(mm,nmbr);         %Maximum RNA pixel intensity experiment
PAR.fitMInt = NaN(mm,nmbr);         %Maximum RNA spot intensity from fit
PAR.TotExpInt = NaN(mm,nmbr);       %Integrated Pixel Intensity for experiment (with background subtracted out)
PAR.TotFitInt = NaN(mm,nmbr);       %Integrated gaussian intensity fit (with background subtraction)
PAR.r = NaN(mm,nmbr);               %Ratio of long to short axis of experimental spot
PAR.rFit = NaN(mm,nmbr);            %Ratio of long to short axis of Fitted gaussian spot
PAR.back = NaN(mm,nmbr);            %Background intensity of the RNA spot
PAR.xsem = NaN(mm,nmbr);            %x fit Standard Error of the Mean in pixels found from the number of photons at each pixel
PAR.ysem = NaN(mm,nmbr);            %y fit Standard Error of the Mean in pixels found from the number of photons at each pixel
PAR.zxfit = NaN(mm,nmbr);           %z postion from x fit
PAR.zyfit = NaN(mm,nmbr);           %z positoin from y fit
PAR.zinit = NaN(mm,nmbr);           %z initial postion from the nongaussian fit.  Integer number
PAR.zint = NaN(mm,nmbr);            %z position from intensity variance through stacks
PAR.zrel = NaN(mm,nmbr);            %Relative z position from the mean of the three z fits
PAR.zabs = NaN(mm,nmbr);            %absolute position in the z stack from the fits
PAR.zstd = NaN(mm,nmbr);            %Nan standard deviation of the z position from the three estimates
PAR.zqFit = NaN(mm,nmbr);           %Quality of zfit.  If <0 pass if >0 curvature is positive at some point and fit is fail
PAR.nucRNA = NaN(mm,nmbr);          %is one if nuclear RNA, 0 if not.  Number changes depending on the nuclear threshold
PAR.distRNA = NaN(mm,nmbr);         %distance in nm of the nuclear RNA to the nuclear membrane
PAR.szNUC = NaN(mm,nmbr);           %size in nm^2 of the nucleus in the z plane of the RNA spot
PAR.cytoRNA = NaN(mm,nmbr);         %is one is cytoplasmic RNA, 0 if not
PAR.xabsloc = NaN(mm,nmbr);         %x Absolute location in the image.  Column!
PAR.yabsloc = NaN(mm,nmbr);         %y aboslute location in image. Row! Plot x y on RNA image is in FunPlotGen.m function.  Must plot(yabsloc,xabsloc) on image!
PAR.normdist = NaN(mm,nmbr);        %Distance of the nuclear RNA to the nuclear membrane normalized by the radius of the nucleus
PAR.nucR = NaN(mm,nmbr);            %Ratio of long axis to short axis in nucleaus.
       
        for q = 1:mm
            len = nnz(RNAposTMR(q,:,1));  %the number of RNA spots in each cell
            
            clear k1 k2 Nuc Cyto NucTMR2d CellTMR1g CellTMR2d L CellTMR2dorg NucTMR2dorg CellRNA CellRNA2 CellRNAorg Nuc Cyto Cell CytoRNA2 NucRNA2 rLabeltcell NUMcell1 rLabeltcyto NUMcyto1 rLabeltnuc NUMnuc1;

            %% single cell
            k1 = Lab == q ; %== q;%sieve out dots in cell q.
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
            % figure(100); clf; imshow(k1,[]); hold on; plot(X0,Y0,'o');hold on; plot(X1,Y1,'or');plot(X0,Y1,'go');plot(X1,Y0,'mo'); title(['Cell: ' num2str(q)]); pixval;

            %% Separate Cell, Nuclear and cyto plasmic RNA
            i = 1;
            clear CellRNAfilter CellRNAindex CellRNAorg CellRNAback Nuc Cyto Cell 

            for i=1:zz;
                i;
                Tf = TMR3Dfilter(Y0:Y1,X0:X1,i);
                CellRNAfilter(:,:,i) = immultiply(Tf,k2); % Filtered image of RNA in the cell        
                CellRNAorg(:,:,i) = immultiply(TMR3Dorg(Y0:Y1,X0:X1,i),k2); % Original RNA image of the cell
                CellRNAindex(:,:,i) = immultiply(TMR3D3immax(Y0:Y1,X0:X1,i),k2);  % Filtered image of RNA in the cell
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
              %figure(102); clf; imshow(k2 + uint16(max(CellRNAindex2,[],3)),[]); title(['Cell: ' num2str(q)]); pixval;

            %% RNA in the cell
            clear rLabeltcell NUMcell1 rLabeltcyto NUMcyto1 rLabeltnuc NUMnuc1
            [rLabeltcell,NUMcell1] = bwlabeln(CellRNAindex2 > 0,26); 
            [rLabeltcyto,NUMcyto1] = bwlabeln(CytoRNA2 > 0,26); 
            [rLabeltnuc,NUMnuc1] = bwlabeln(NucRNA2 > 0,26); 
             %figure(103); clf; imshow(k2 + uint16(max(rLabeltcell,[],3)),[]); title(['Cell: ' num2str(q)]); %pixval;

            len = nnz(RNAposTMR(q,:,1));  %the number of RNA spots in each cell
            
            if len ~= 0 %only fit if there is an RNA spot
                for cnt = 1:len

                k = RNAposTMR(q,cnt,3);    %The z-stack containing the RNA position
                
                %% set lower and upper limits in the stack
                    if k < 2;
                        k3 = k-1;
                    else;
                        k3 = 2;
                    end;
                    if k > zz - 2;
                        k4 = 1;
                    else;
                        k4 = 2;
                    end;
                    %% include limits for x and y direction
                    clear TMRfree TMRim1 yA xA TMRim
               
                    xA = RNAposTMR(q,cnt,2);
                    yA =  RNAposTMR(q,cnt,1);
                    
 %                   yA = xall(xynan(i),k);
 %                   xA = yall(xynan(i),k);
                    if k-k3 == 0;
                        TMRfree(:,:,1) = zeros(size(CellRNAorg2,1),size(CellRNAorg2,2));
                    else
                        TMRfree(:,:,1) = CellRNAorg2(:,:,k-k3) - Cellback2(:,:,k-k3); % generate image stack -2:+2
                    end;
                    TMRfree(:,:,2) = CellRNAorg2(:,:,k-k3+1) - Cellback2(:,:,k-k3+1);
                    TMRfree(:,:,3) = CellRNAorg2(:,:,k) - Cellback2(:,:,k);
                    TMRfree(:,:,4) = CellRNAorg2(:,:,k+k4-1) - Cellback2(:,:,k+k4-1);
                    if k + k4 == zz+1;
                        TMRfree(:,:,5) = zeros(size(CellRNAorg2,1),size(CellRNAorg2,2));
                    else;
                        k + k4;
                        TMRfree(:,:,5) = CellRNAorg2(:,:,k+k4)  - Cellback2(:,:,k+k4);
                    end;
                    q
                    cnt
                    TMRim1(:,:,1) = TMRfree(xA-d:xA+d,yA-d:yA+d,1); % generate image for single pixcel 1:9 x 1:9
                    TMRim1(:,:,2) = TMRfree(xA-d:xA+d,yA-d:yA+d,2);
                    TMRim1(:,:,3) = TMRfree(xA-d:xA+d,yA-d:yA+d,3);
                    TMRim1(:,:,4) = TMRfree(xA-d:xA+d,yA-d:yA+d,4);
                    TMRim1(:,:,5) = TMRfree(xA-d:xA+d,yA-d:yA+d,5);
    %                 logTMRfree1 = TMRfree(:,:,3) >= 0;
    %                 logTMRfree2 = log(immultiply(TMRfree(:,:,3),logTMRfree1) + 1);
%     figure(105); clf; subplot(1,2,1); imshow(TMRfree(:,:,3),[]); hold on; plot(yA,xA,'o'); hold on; title(['Cell: ' num2str(q) ', Plane: ' num2str(k)]); ...
%     subplot(1,2,2); imshow(TMRim1(:,:,3),[]); title(['RNA: ' num2str(i) ', maxRNA: ' num2str(size(xynan,1))]); impixelinfo;

                    m =1;
                    for m = 1:5;
                        TMRb(m,1) = sum(sum(immultiply(TMRim1(:,:,m),double(b)),1))./sum(b(:)); % deterimine background from border pixvels;
                        TMRim(:,:,m) = immultiply(TMRim1(:,:,m),double(b1)); % set boundary pxvels to zero and reduce the second line pixel intensity by 1/2
                        TMRim(:,:,m) = round(immultiply(TMRim1(:,:,m),double(b2))./1);
                    end;
%     figure(106); clf; subplot(1,2,1); imshow(max(TMRim1,[],3),[]); title(['Cell: ' num2str(q) ', Plane: ' num2str(k)]); ...
%     subplot(1,2,2); imshow(max(TMRim,[],3),[]); title(['RNA: ' num2str(i)  ', maxRNA: ' num2str(size(xynan,1))]); impixelinfo;

                    TMRint(:,1) = sum(sum(TMRim,1),2);  % Determine total fluorecnes under the peak
        %            [X1,Y1] = meshgrid(1:1:2*d+1);
                    par0(1) = 2; % initial condition g1 x-position
                    par0(2) = 3; % initial condition g1 y-position
                    par0(3) = 1; % initial condition s1 width of gaussian
                    par0(4) = double(max(TMRim(:))); % initial condition A1 intensity of gaussian
                    mp = size(TMRim,1); % size of image to fit
                    [X1,Y1] = meshgrid(1:1:mp); 
                    
                    options = optimset('Display', 'none','MaxFunEvals', 200,'MaxIter', 200,'TolX', 1E-9); % fit conditions
                    for jj = 1:5;
                        TMRimm = double(TMRim(:,:,jj));% Spot image
                        TMRimm2(:,:,jj) = TMRimm;
                        TMRimm1 = TMRimm./max(TMRimm(:)); 
                        TMRimm1 = im2bw(TMRimm1,0.5); % reduce pixel below 1/2 of max 
                        IR = regionprops(double(TMRimm1),'MinorAxisLength','MajorAxisLength'); % determine shortest / longest axis of the spot
                        try r(jj) = IR.MajorAxisLength ./ IR.MinorAxisLength; end; % deterine ratio of the two axis => circular = 1 single spot; >1.3 multiple spots
                        lb = [2 2 par0(3) par0(4)]; % lower bounds for finding the position in x and y
                        ub = [mp-2 mp-2  par0(3) 5*par0(4)] + 0.001; % upper bounds for finding the position in x and y
                        [par,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Gauss3DFit,par0,mp,TMRimm,lb,ub,options); % least square fit for the position      
                        lb = [par(1) par(2) 0 0]; % lower bounds for finding the gaussian width
                        ub = [par(1) par(2) 5 2*par(4)] + 0.001; % upper bounds for finding the gaussian width
                        [par2,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Gauss3DFit,par0,mp,TMRimm,lb,ub,options); % least square fit for the gaussian width and amplitude  
                        lb = [par2(1) par2(2) 0 0 0]; % lower bounds for finding the gaussian width
                        ub = [par2(1) par2(2) par2(3)*3 par2(3)*3 par2(4)*3] + 0.001; % upper bounds for finding the gaussian width
                        par3([1 2 3 5]) = par2;
                        par3(4) = par2(3);
                        [par4,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Gauss3DFitS1S2,par3,mp,TMRimm,lb,ub,options); % least square fit for the gaussian width in x and y and amplitude  
                        g3 = par4(5)*exp(-((X1-par4(1)).^2)./(2*par4(3)^2)-((Y1-par4(2)).^2)./(2*par4(4)^2)); % compute gaussain after fit
                        gAll(:,:,jj) = g3;
                        par4;
                        par4(1) = par4(1)+xA-d-1; %fitted x-position
                        par4(2) = par4(2)+yA-d-1; %fitted y-position
                        par5(1:5,jj) = par4; % fitted parameter 3: gauss width; 4: Gauss amplitude
                    end;
                    try PAR.xfit(q,cnt) = par5(1,3); end; % fitted x-parameter
                    try PAR.yfit(q,cnt) = par5(2,3); end; % fitted y-parameter
                    try PAR.xinit(q,cnt) = xA; end; %initial x-position
                    try PAR.yinit(q,cnt) = yA; end; % initial; y-position
                    try PAR.xgw(q,cnt) = par5(3,3); end; % fitted gaussian width x
                    try PAR.ygw(q,cnt) = par5(4,3); end; % fitted gaussian width y
                    try PAR.xFWHM(q,cnt) = 2*sqrt(2*log(2)).*par5(3,3); end; % fitted FWHM x
                    try PAR.yFWHM(q,cnt) = 2*sqrt(2*log(2)).*par5(4,3); end; % fitted FWHM y
                    try PAR.expMInt(q,cnt) = double(max(TMRimm2(:))); end;  % experimental maximum spot intensity
                    try PAR.fitMInt(q,cnt) = par5(5,3); end; % fitted maximum spot intensity
                    try PAR.TotExpInt(q,cnt) = nansum(nansum(TMRimm2(:,:,3))); end; % experimental integrated spot intensity
                    try PAR.TotFitInt(q,cnt) = nansum(nansum(gAll(:,:,3))); end; % fitted integrated spot intensity
                    try PAR.r(q,cnt) = r(3); end; % ratio of shortest / longest spot axis; based on measured RNA spot properties
                    try 
                        if PAR.xFWHM(q,cnt) > PAR.yFWHM(q,cnt)
                            PAR.rFit(q,cnt) = PAR.xFWHM(q,cnt)./PAR.yFWHM(q,cnt);
                        else;
                            PAR.rFit(q,cnt) = PAR.yFWHM(q,cnt)./PAR.xFWHM(q,cnt);
                        end
                    end
                    try PAR.back(q,cnt) = TMRb(3,1); end; % background intensity
                    ps = 65; % Hamamatsu Orca Flash 4.2 ps = 65nm 
                    %See Thompson Biophysical Journal 82(5) 2775-2783 eq. 17  
                    %s=standard deviation of point spread function
                    %(PAR.xgw); a = pixel size (ps);  b = background noise (PAR.back)   N = number of photons collected (PAR.TotExpInt)
                    %This result is in pixels.  to convert to nm multiply by 65nm
                    try PAR.xsem(q,cnt) = sqrt(((PAR.xgw(q,cnt)*ps)^2)/PAR.TotExpInt(q,cnt) + (ps^2/(12*PAR.TotExpInt(q,cnt)) + (8*pi()*(PAR.xgw(q,cnt)*ps)^4)*PAR.back(q,cnt).^2)./(ps^2*PAR.TotExpInt(q,cnt).^2))./ps; end; % standard error of the mean for x-position
                    try PAR.ysem(q,cnt) = sqrt(((PAR.ygw(q,cnt)*ps)^2)/PAR.TotExpInt(q,cnt) + (ps^2/(12*PAR.TotExpInt(q,cnt)) + (8*pi()*(PAR.ygw(q,cnt)*ps)^4)*PAR.back(q,cnt).^2)./(ps^2*PAR.TotExpInt(q,cnt).^2))./ps; end; % standard error of the mean for y-position

                    % Determine z-position from gaussian width in x-direction
                    ydata1 = double(par5(3,:));
                    [s,I] = min(ydata1);
                     if (I > 2 && I < 5);
                        pardf = double([1 3 2 0.67 0.3]);
                        lb = [0 0 0 0 0]; % lower bounds for finding the gaussian width
                        ub = [10 10 10 10 10] + 0.001; % upper bounds for finding the gaussian width
                        z = 1:5;
                        options = optimset('Display', 'none','MaxFunEvals', 200,'MaxIter', 200,'TolX', 1E-9); % fit conditions
                        [dfpar] = lsqcurvefit(@Defoc4order,pardf,z,ydata1,lb,ub,options); % least square fit for the gaussian width and amplitude   ,resnorm,residual,exitflag,output,lambda,jacobian
                        z1 = [1:0.01:5];
                        df1 = dfpar(1).*sqrt(1 + ((z1-dfpar(2))./dfpar(3)).^2 + dfpar(4).*((z1-dfpar(2))./dfpar(3)).^3 + dfpar(5).*((z1-dfpar(2))./dfpar(3)).^4 );
 %                        figure; plot(ydata1,'o'); hold on; plot(z1,df1); hold on;
                        [e,gm1] = min(df1);
                        gm1 = z1(gm1);
                    else;
                        gm1 = NaN;
                    end;
                    try PAR.zxfit(q,cnt) = gm1; end; % fitted z-position from x-gaussian width

                    % Determine z-position from gaussian width in y-direction
                    ydata2 = double(par5(4,:));
                    [s,I] = min(ydata2);
                     if (I > 2 && I < 5);
                        pardf = double([1 3 2 0.67 0.3]);
                        lb = [0 0 0 0 0]; % lower bounds for finding the gaussian width
                        ub = [10 10 10 10 10] + 0.001; % upper bounds for finding the gaussian width
                        z = 1:5;
                        options = optimset('Display', 'none','MaxFunEvals', 200,'MaxIter', 200,'TolX', 1E-9); % fit conditions
                        [dfpar] = lsqcurvefit(@Defoc4order,pardf,z,ydata2,lb,ub,options); % least square fit for the gaussian width and amplitude   ,resnorm,residual,exitflag,output,lambda,jacobian
                        z1 = [1:0.01:5];
                        df2 = dfpar(1).*sqrt(1 + ((z1-dfpar(2))./dfpar(3)).^2 + dfpar(4).*((z1-dfpar(2))./dfpar(3)).^3 + dfpar(5).*((z1-dfpar(2))./dfpar(3)).^4 );
%                        plot(ydata2,'or'); hold on; plot(z1,df2,'r');
                        [e,gm2] = min(df2);
                        gm2 = z1(gm2);
                    else;
                        gm2 = NaN;
                    end;
                    try PAR.zyfit(q,cnt) = gm2; end; % fitted z-position from y-gaussian width

                    % Determine z-position from gaussian intensity
                    ydata3 = par5(5,:);
                    [s,I] = max(ydata3);
                    if (I > 2 && I < 5);
                        p3 = polyfit([1:5],ydata3,4);
                        x2 = 1:0.01:5;
                        [e,gm3] = max(polyval(p3,x2));
 %                       crvtur = diff(diff(polyval(p3,x2(gm3-10:gm3+10)),1));
 %                       if crvtur(gm3) < 0
 %                           r1 = min(diff(diff(polyval(p3,x2),2)));
 %                       else;
 %                           r1 = 0;
 %                       end
                        r1 = max(diff(diff(polyval(p3,x2),2)));
                        gm3 = x2(gm3);
                        
                        %% compute mean and std of z position from gaussoan width in x and y and the gauss intensity
                        gm = nanmean([gm1 gm2 gm3]);
                        gms = nanstd([gm1 gm2 gm3]);
                        gf = r1;
                    else;
                        gm3 = NaN;
                        gm = NaN;
                        gms = NaN;
                        gf = NaN;
                    end;
                    try PAR.zinit(q,cnt) = k; end;%initial frame height from nongaussian fit
                    try PAR.zint(q,cnt) = gm3; end; % fitted z-position from intensity
                    try PAR.zrel(q,cnt) = gm; end;% relative mean from three fitted z-position
                    try PAR.zabs(q,cnt) = gm + k-3; end;% absolute mean from three fitted z-position
                    try PAR.zstd(q,cnt) = gms; end;% STD from three fitted z-position
                    try PAR.zqFit(q,cnt) = r1; end;% Quality of fit through the curviture of the fit, >=0 failed; <0 curvature is negative everywhere (passed)
                    try PAR.nucRNA(q,cnt) = Nuc(round(par5(1,3)),round(par5(2,3)),k) > 0;% Nuclear RNA
                        if Nuc(round(par5(1,3)),round(par5(2,3)),k) > 0;
                            for i=1:(round(abs(X1-X0)))
                                if Nuc(round(par5(1,3)+i),round(par5(2,3)),k) == 0
                                    xloc = i;       %Number of pixels to nuclear membrane
                                    break;
                                elseif Nuc(round(par5(1,3)-i),round(par5(2,3)),k) == 0
                                    xloc = i;
                                    break;
                                end
                            end
                            for i =1:(round(abs(Y1-Y0)))
                                if Nuc(round(par5(1,3)),round(par5(2,3)+i),k) == 0
                                    yloc = i;
                                    break;
                                elseif Nuc(round(par5(1,3)),round(par5(2,3)-i),k) == 0
                                    yloc = i;
                                    break;
                                end
                            end
                            if xloc < yloc
                                PAR.distRNA(q,cnt) = xloc*ps;  %the result of the distance in nm of all nuclear RNA
                            elseif yloc < xloc
                                PAR.distRNA(q,cnt) = yloc*ps;
                            end
                            PAR.szNUC(q,cnt) = nnz(Nuc(:,:,k))*ps.^2;
                            PAR.normdist(q,cnt) = PAR.distRNA(q,cnt)./sqrt(PAR.szNUC(q,cnt)./pi()); %normalized to radius
                        end
                    end;
                    %walk in four directions from center of image until the
                    %edge is reached to find the approximate ellipticity
                    try 
                        clear len
                        len = size(Nuc(:,:,k),1);
                        len = round(len./2);
                        for i = 1:(round(abs(X1-X0)))
                            if Nuc(len+i,len,k) == 0
                                mxloc = abs(i);
                                break;
                            end
                        end
                        
                        for i = 1:(round(abs(X1-X0)))
                            if Nuc(len-i,len,k) == 0
                                pxloc = abs(i);
                                break;
                            end
                        end
                        
                        for i = 1:(round(abs(Y1-Y0)))
                            if Nuc(len,len+i,k) == 0
                                myloc = abs(i);
                                break;
                            end
                        end
                        
                        for i = 1:(round(abs(Y1-Y0)))
                            if Nuc(len,len-i,k) == 0
                                pyloc = abs(i);
                                break
                            end
                        end
                        
                        PAR.nucR(q,cnt) = (mxloc+pxloc)./(myloc+pyloc);
                    end
                    
                    try PAR.cytoRNA(q,cnt) = Nuc(round(par5(1,3)),round(par5(2,3)),k) == 0; end;% Cytoplasmic
                    try PAR.xabsloc(q,cnt) = Y0+PAR.xfit(q,cnt)-1; end; %Row!
                    try PAR.yabsloc(q,cnt) = X0+PAR.yfit(q,cnt)-1; end; %Colunm! 

%                     figure(110); clf; 
%                     subplot(4,6,1);imshow(TMRim(:,:,1),[]); hold on;title(['x-fit:' num2str(PAR.xfit(q,cnt))]);
%                     subplot(4,6,2);imshow(TMRim(:,:,2),[]);title(['y-fit:' num2str(PAR.yfit(q,cnt))]);
%                     subplot(4,6,3);imshow(TMRim(:,:,3),[]);title(['x-init:' num2str(PAR.xinit(q,cnt))]);
%                     subplot(4,6,4);imshow(TMRim(:,:,4),[]);title(['y-init:' num2str(PAR.yinit(q,cnt))]);
%                     subplot(4,6,5);imshow(TMRim(:,:,5),[]);title(['x-gw:' num2str(PAR.xgw(q,cnt))]);
%                     subplot(4,6,6);plot([1:5],ydata1,'o'); hold on; try plot(z1,df1); end; title(['y-gw:' num2str(PAR.ygw(q,cnt))]);
%                     subplot(4,6,7);imshow(gAll(:,:,1),[]);title(['x-FWHM fit:' num2str(PAR.xFWHM(q,cnt))]);
%                     subplot(4,6,8);imshow(gAll(:,:,2),[]);title(['y-FWHM fit:' num2str(PAR.yFWHM(q,cnt))]);
%                     subplot(4,6,9);imshow(gAll(:,:,3),[]);title(['maxint:' num2str(PAR.expMInt(q,cnt))]);
%                     subplot(4,6,10);imshow(gAll(:,:,4),[]);title(['fit int:' num2str(PAR.fitMint(q,cnt))]);
%                     subplot(4,6,11);imshow(gAll(:,:,5),[]);title(['SumInt:' num2str(PAR.TotExpInt(q,cnt))]);
%                     subplot(4,6,12);plot([1:5],ydata2,'o'); hold on; try plot(z1,df2); end;title(['fitSumInt:' num2str(PAR.TotFitInt(q,cnt))]);
%                     subplot(4,6,13);surf(X1,Y1,double(TMRim(:,:,1))./TMRint(1,1),'FaceColor','interp','EdgeColor','black'); axis([1 9 1 9 0 0.2]);title(['R:' num2str(PAR.r(q,cnt))]);
%                     subplot(4,6,14);surf(X1,Y1,double(TMRim(:,:,2))./TMRint(2,1),'FaceColor','interp','EdgeColor','black','FaceLighting','phong'); axis([1 9 1 9 0 0.2]);title(['bInt:' num2str(PAR.back(q,cnt))]);
%                     subplot(4,6,15);surf(X1,Y1,double(TMRim(:,:,3))./TMRint(3,1),'FaceColor','interp','EdgeColor','black','FaceLighting','phong'); axis([1 9 1 9 0 0.2]);title(['SEM x:' num2str(PAR.xsem(q,cnt))]);
%                     subplot(4,6,16);surf(X1,Y1,double(TMRim(:,:,4))./TMRint(4,1),'FaceColor','interp','EdgeColor','black','FaceLighting','phong'); axis([1 9 1 9 0 0.2]);title(['SEM y:' num2str(PAR.ysem(q,cnt))]);
%                     subplot(4,6,17);surf(X1,Y1,double(TMRim(:,:,5))./TMRint(5,1),'FaceColor','interp','EdgeColor','black','FaceLighting','phong'); axis([1 9 1 9 0 0.2]);title(['z-pos-x-gw:' num2str(PAR.zxfit(q,cnt))]);
%                     subplot(4,6,18);plot([1:5],ydata3,'o'); hold on; try plot(x2,polyval(p3,x2)); end; title(['z-pos-y-gw:' num2str(PAR.zyfit(q,cnt))]);
%                     subplot(4,6,19);surf(X1,Y1,double(gAll(:,:,1))./PAR.TotFitInt(q,cnt),'FaceColor','interp','EdgeColor','black'); axis([1 9 1 9 0 0.2]); title(['z-pos-Int:' num2str(PAR.zint(q,cnt))]);
%                     subplot(4,6,20);surf(X1,Y1,double(gAll(:,:,2))./PAR.TotFitInt(q,cnt),'FaceColor','interp','EdgeColor','black','FaceLighting','phong'); axis([1 9 1 9 0 0.2]); title(['z-rel:' num2str(PAR.zrel(q,cnt))]);
%                     subplot(4,6,21);surf(X1,Y1,double(gAll(:,:,3))./PAR.TotFitInt(q,cnt),'FaceColor','interp','EdgeColor','black','FaceLighting','phong'); axis([1 9 1 9 0 0.2]); title(['z-abs:' num2str(PAR.zabs(q,cnt)) ', nRNA:' num2str(PAR.nucRNA(q,cnt))]);
%                     subplot(4,6,22);surf(X1,Y1,double(gAll(:,:,4))./PAR.TotFitInt(q,cnt),'FaceColor','interp','EdgeColor','black','FaceLighting','phong'); axis([1 9 1 9 0 0.2]); title(['z-std:' num2str(PAR.zstd(q,cnt)) ', cRNA:' num2str(PAR.cytoRNA(q,cnt))]);
%                     subplot(4,6,23);surf(X1,Y1,double(gAll(:,:,5))./PAR.TotFitInt(q,cnt),'FaceColor','interp','EdgeColor','black','FaceLighting','phong'); axis([1 9 1 9 0 0.2]); title(['fit-Q:' num2str(PAR.zqFit(q,cnt))]);
%                     subplot(4,6,24); imshow(k2 + uint16(max(CellRNAorg,[],3)),[]); hold on; plot(yA,xA,'o'); hold on;title(['Cell: ' num2str(j) ', Plane:' num2str(k) ', RNA: ' num2str(i)]); impixelinfo;
%                 pause(5)
%                

                end;
       
            end;
        end;        
end







%code now in B3_nongaussRNApos.m function
% for j = 1:mm;
%     j
%     %tic; 
%     clear k1 k2 Nuc Cyto NucTMR2d CellTMR1g CellTMR2d L CellTMR2dorg NucTMR2dorg CellRNA CellRNA2 CellRNAorg Nuc Cyto Cell CytoRNA2 NucRNA2 rLabeltcell NUMcell1 rLabeltcyto NUMcyto1 rLabeltnuc NUMnuc1;
% 
%     %% single cell
%     k1 = Lab == j ; %== j;%sieve out dots in cell j.
%     %figure; imshow(k1,[0 1]);
%     k1=uint16(k1);
%     k3=regionprops(k1,'BoundingBox','Area');
%     k4 = k3.BoundingBox; %create the rectangular box around the cell. 
%     X0=round(k4(1))-4;
%     Y0=round(k4(2))-4;
%     X1=round(k4(1)+ k4(3))+4;
%     Y1=round(k4(2)+ k4(4))+4;
%     if X0 < 1;
%         X0 = 1;
%     else;
%     end;
%     if Y0 < 1;
%         Y0 = 1;
%     else;
%         
%     end;
%     if X1 > A;
%         X1 = A;
%     else;
%     end;
%     if Y1 > A;
%         Y1 = A;
%     else;
%     end;
%     k2 = k1(Y0:Y1,X0:X1);
% % figure(100); clf; imshow(k1,[]); hold on; plot(X0,Y0,'o');hold on; plot(X1,Y1,'or'); title(['Cell: ' num2str(j)]); pixval;
% 
%     %% Separate Cell, Nuclear and cyto plasmic RNA
%     i = 1;
%     clear CellRNAfilter CellRNAindex CellRNAorg CellRNAback Nuc Cyto Cell 
%     
%     for i=1:zz;
%         Tf = TMR3Dfilter(Y0:Y1,X0:X1,i);
%         CellRNAfilter(:,:,i) = immultiply(Tf,k2); % Filtered image of RNA in the cell        
%         CellRNAindex(:,:,i) = immultiply(TMR3D3immax(Y0:Y1,X0:X1,i),k2);  % Filtered image of RNA in the cell
%         CellRNAorg(:,:,i) = immultiply(TMR3Dorg(Y0:Y1,X0:X1,i),k2); % Original RNA image of the cell
%         CellRNAback(:,:,i) = immultiply(TMR3Dback(Y0:Y1,X0:X1,i),k2); % Original RNA image of the cell
%         Nuc(:,:,i) = immultiply(double(Nuc3D(Y0:Y1,X0:X1,i)),double(k2)); % Nucleus object
%         Cyto(:,:,i) = imsubtract(double(k1(Y0:Y1,X0:X1)),double(Nuc(:,:,i))); % Cytoplasm object
%         Cell(:,:,i) = double(k1(Y0:Y1,X0:X1)); 
%     end;
% 
% %% Cell, Cytoplasmic, Nuclear RNA
%     clear CellRNAorg2 CellRNAindex2 Cellback2 CytoRNA2 NucRNA2
%     CellRNAorg2 = immultiply(double(CellRNAorg),Cell);     
%     CellRNAindex2 = immultiply(double(CellRNAindex),Cell); 
%     Cellback2 = immultiply(double(CellRNAback),Cell); 
%     CytoRNA2 = immultiply(double(CellRNAindex),double(Cyto)); 
%     NucRNA2 = immultiply(double(CellRNAindex),double(Nuc)); 
%     %  figure(102); clf; imshow(k2 + uint16(max(CellRNAindex2,[],3)),[]); title(['Cell: ' num2str(j)]); pixval;
% 
%     %% RNA in the cell
%     clear rLabeltcell NUMcell1 rLabeltcyto NUMcyto1 rLabeltnuc NUMnuc1
%     [rLabeltcell,NUMcell1] = bwlabeln(CellRNAindex2 > 0,26); 
%     [rLabeltcyto,NUMcyto1] = bwlabeln(CytoRNA2 > 0,26); 
%     [rLabeltnuc,NUMnuc1] = bwlabeln(NucRNA2 > 0,26); 
%      % figure(103); clf; imshow(k2 + uint16(max(rLabeltcell,[],3)),[]); title(['Cell: ' num2str(j)]); pixval;
%     
%     %% Count RNA in the cell, cytoplasm, nucleus
%     CELLmaxRNA(j,1) = NUMcell1;
%     CELLmaxRNA(j,2) = NUMcyto1;
%     CELLmaxRNA(j,3) = NUMnuc1;
%     
%     %% Deterime x-y positions in each cell for initial fits
%     clear xall yall
%     th1 = 0;
%     xall = NaN(100,zz);
%     yall = NaN(100,zz);
%     try
%     for kk = 1:zz;
%         [y1,x1] = find(CellRNAindex2(:,:,kk) > th1);
%         xall(1:size(x1,1),kk) = x1;
%         yall(1:size(y1,1),kk) = y1;
%         x2 = find(x1 == 2);
%         y2 = find(y1 == 2);
%         xall(x2,kk) = NaN;
%         yall(x2,kk) = NaN;
%         xall(y2,kk) = NaN;
%         yall(y2,kk) = NaN;
%         x2 = find(x1 == 1023);
%         y2 = find(y1 == 1023);
%         xall(x2,kk) = NaN;
%         yall(x2,kk) = NaN;
%         xall(y2,kk) = NaN;
%         yall(y2,kk) = NaN;
%     end;
%     end;
%     
% end;