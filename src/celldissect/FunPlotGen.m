function FunPlotGen(Lab, PAR, Ych, mm, IMorg)

IMorg = CY53Dorg;
PAR = PARtmr_mid;
ttl = 'TMR'
figure()
clf;
%defined outlier as being more than 3 std deviations outside of mean value
imshow(max(CY53Dfilter,[],3),[]);
%hold on;
for i=1:mm
    len = sum(~isnan(PAR.TotExpInt(i,:)));
    if len ~= 0
          clear k1 k2 k3 k4 CellRNAorg Cell CellRNAorg2
        A = size(IMorg,1);
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
            CellRNAorg(:,:,j)=immultiply(IMorg(Y0:Y1,X0:X1,j),k2);
            Cell(:,:,j) = double(k1(Y0:Y1,X0:X1));
        end
 
        
%       imshow(k2 + uint16(max(CellRNAorg,[],3)),[]); 
       hold on; title(['cell: ' num2str(i) ' r=outlier; b=25-75%; g=other']);
        %plot(X0,Y0,'r*'); plot(X1,Y1,'b*');plot(X1,Y0,'g*');plot(X0,Y1,'m*');
        qntl = quantile(PAR.TotExpInt(i,:),[.25 .75]);

        clear outliers
        outliers = nan;
        for j=1:len
            if abs(PAR.TotExpInt(i,j)-nanmean(PAR.TotExpInt(i,:))) > 3*nanstd(PAR.TotExpInt(i,:))
                outliers = [outliers,j];
            end
        end

       for j=1:len
           if PAR.TotExpInt(i,j) > qntl(1) & PAR.TotExpInt(i,j) < qntl(2)
               x = Y0+PAR.xfit(i,j)-1;
               y = X0+PAR.yfit(i,j)-1;
               plot(y,x,'bo')
           else
               x = Y0+PAR.xfit(i,j)-1;
               y = X0+PAR.yfit(i,j)-1;
               plot(y,x,'go')
       end

       if size(outliers,2)>1

           for j=2:size(outliers,2)
               x = Y0+PAR.xfit(i,outliers(j))-1;
               y = X0+PAR.yfit(i,outliers(j))-1;
               plot(y,x,'ro')
           end
       end

       end
end
    


end






%junk code
%         clear k1 k2 k3 k4 CellRNAorg Cell CellRNAorg2
%         k1 = Lab == q ; %== q;%sieve out dots in cell q
%         k1=uint16(k1);
%         k3=regionprops(k1,'BoundingBox','Area');
%         k4 = k3.BoundingBox; %create the rectangular box around the cell. 
%         X0=round(k4(1))-4;
%         Y0=round(k4(2))-4;
%         X1=round(k4(1)+ k4(3))+4;
%         Y1=round(k4(2)+ k4(4))+4;
%         if X0 < 1;
%         X0 = 1;
%         else;
%         end;
%         if Y0 < 1;
%         Y0 = 1;
%         else;
%         end;
%         if X1 > A;
%         X1 = A;
%         else;
%         end;
%         if Y1 > A;
%         Y1 = A;
%         else;
%         end;
%         k2 = k1(Y0:Y1,X0:X1);

%         for j = 1:zz
%             CellRNAorg(:,:,j)=immultiply(IMorg(Y0:Y1,X0:X1,j),k2);
%              Cell(:,:,j) = double(k1(Y0:Y1,X0:X1));
%         end
%         clf;
%         CellRNAorg2 = immultiply(double(CellRNAorg),Cell);
%         imshow(k2 + uint16(max(CellRNAorg2,[],3)),[]); hold on; 