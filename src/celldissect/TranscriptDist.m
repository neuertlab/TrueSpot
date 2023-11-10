figure(2)
clf; hold on;
for i=1:mm
    len = sum(~isnan(PARcy5_mid.TotExpInt(i,:)));
    if len ~= 0
        for j = 1:len
            plot(PARcy5_mid.zabs(i,j),PARcy5_mid.szNUC(i,j),'r*', 'linewidth', 2)
        end
    end
end

PARcy5_mid.distRNA = PARcy5_mid.distRNA./65;


figure(3)
clf; hold on;
for i=1:mm
    len = sum(~isnan(PARcy5_mid.TotExpInt(i,:)));
    if len ~= 0
        for j = 1:len
            plot(PARcy5_mid.zabs(i,j),PARcy5_mid.szNUC(i,j),'r*', 'linewidth', 2)
        end
    end
end


histogram(PARcy5_mid.distRNA,10);
figure(4)
clf; hold on;
histogram(PARcy5_mid.normdist,5);

%Transcription distribution
m=1
for i=1:mm
    len = sum(~isnan(PARcy5_mid.TotExpInt(i,:)));
    if len ~= 0
        for z = 1:len
            if isnan(PARcy5_mid.TotFitInt(i,z))
                cnt = cnt+1;
            else
                foo(m,1)= PARcy5_mid.TotExpInt(i,z);
                foo(m,2)= PARcy5_mid.TotFitInt(i,z);
                m = m+1;

            end 
        end

    end
end

clear int celtot 
outcell=0; m = 1;

%for finding the number of cells and rna as a function of time
szCell = size(PARcy5_mid.xfit,1);
szRNA = size(PARcy5_mid.xfit,2);

%for finding the experimental mean and median
if ~isnan(mean(nanmean(PARcy5_mid.TotExpInt)));
    ttmen = mean(nanmean(PARcy5_mid.TotExpInt));
    ttmed = median(nanmedian(PARcy5_mid.TotExpInt));
end
outliers = NaN(szCell*szRNA,1);

for i=1:mm
    clear val
    len = sum(~isnan(PARcy5_mid.TotExpInt(i,:)));
    if len > 1  %have at least 2 RNA in the cell
        val = PARcy5_mid.TotExpInt(i,~isnan(PARcy5_mid.TotExpInt(i,:)));
%            qntl = quantile(val,[.25 .75]);
      for z =1:len
        if abs(val(z)-mean(foo(:,1))) > 3*std(foo(:,1))
            outliers(i+z,1) = PARcy5_mid.normdist(i,z);
        end
      end
    end 
end




figure(5)
histogram(outliers,20);

figure(6)
histogram(PARcy5_mid.szNUC,10);
    