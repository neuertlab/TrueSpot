%Transcription plot
%Must have already run one threshold in
%B_ClusterImageProcessingTMR_AF594_CY5_B.m with all the variables in the
%work space.

figure(101)
clf;
%defined outlier as being more than 3 std deviations outside of mean value
imshow(max(CY53Dfilter,[],3),[]);
zz = size(CY53Dorg,3);
%hold on;
for i=1:mm
    len = sum(~isnan(PARcy5_mid.TotExpInt(i,:)));
    if len ~= 0
          clear k1 k2 k3 k4 CellRNAorg Cell CellRNAorg2
        A = size(CY53Dorg,1);
        k1 = Lab == i ; %== q;%sieve out dots in cell q
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

        for j = 1:zz
            CellRNAorg(:,:,j)=immultiply(CY53Dorg(Y0:Y1,X0:X1,j),k2);
            Cell(:,:,j) = double(k1(Y0:Y1,X0:X1));
        end
 
        
%imshow(k2 + uint16(max(CellRNAorg,[],3)),[]); 
       hold on; title(['cell: ' num2str(i) ' r=outlier; b=25-75%; g=other']);
%plot(X0,Y0,'r*'); plot(X1,Y1,'b*');plot(X1,Y0,'g*');plot(X0,Y1,'m*');  %corners of each cell
        qntl = quantile(PARcy5_mid.TotExpInt(i,:),[.25 .75]);

        clear outliers
        outliers = nan;
        for j=1:len
            if abs(PARcy5_mid.TotExpInt(i,j)-nanmean(PARcy5_mid.TotExpInt(i,:))) > 3*nanstd(PARcy5_mid.TotExpInt(i,:))
                outliers = [outliers,j];
            end
        end

       for j=1:len
           if PARcy5_mid.TotExpInt(i,j) > qntl(1) & PARcy5_mid.TotExpInt(i,j) < qntl(2)
               x = Y0+PARcy5_mid.xfit(i,j)-1;
               y = X0+PARcy5_mid.yfit(i,j)-1;
               plot(y,x,'bo')
           else
               x = Y0+PARcy5_mid.xfit(i,j)-1;
               y = X0+PARcy5_mid.yfit(i,j)-1;
               plot(y,x,'go')
       end

       if size(outliers,2)>1

           for j=2:size(outliers,2)
               x = Y0+PARcy5_mid.xfit(i,outliers(j))-1;
               y = X0+PARcy5_mid.yfit(i,outliers(j))-1;
               plot(y,x,'ro')
           end
       end

       end
    end
end

    