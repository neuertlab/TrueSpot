sum(sum(sqrt((PARtmr_mid.yfit(:,:)-PARtmr_mid.yinit).^2+(PARtmr_mid.xfit(:,:)-PARtmr_mid.xinit(:,:)).^2)))./(size(PARtmr_mid.xfit,1)*size(PARtmr_mid.xfit,2))


sum(sum(PARtmr_mid.zabs(:,:)-PARtmr_mid.zinit(:,:)))/(size(PARtmr_mid.xfit,1)*size(PARtmr_mid.xfit,2))


PAR.TotExpInt(q,cnt) = nansum(nansum(TMRimm2(:,:,3))); % experimental integrated spot intensity
PAR.TotFitInt(q,cnt) = nansum(nansum(gAll(:,:,3))); % fitted integrated spot intensity


len = (size(PARtmr_mid.xfit,1)*size(PARtmr_mid.xfit,2))
sum(sum((PARtmr_mid.TotExpInt(:,:)-PARtmr_mid.TotFitInt(:,:))))/len

mean((PARtmr_mid.TotExpInt(:,:)))
mean((PARtmr_mid.TotFitInt(:,:)))

max(PARtmr_mid.TotFitInt(:,:))

max(max(PARtmr_mid.

clf;

hold on;
for i=1:mm
    plot(i,(PARtmr_mid.TotFitInt(i,:)),'bo')
    
end;

figure();
clf;

hold on;
for i=1:mm
    plot(i,(PARtmr_mid.TotExpInt(i,:)),'bo')
    
end;

max(nonzeros((PARtmr_mid.TotExpInt(:,:))))
mean(nonzeros((PARtmr_mid.TotExpInt(:,:))))

figure();
hold on;
for i=1:mm
x=std(nonzeros((PARtmr_mid.TotExpInt(i,:))));
plot(i,x,'bo');
end

figure();
hold on;
for i=1:mm
x=mean(nonzeros((PARtmr_mid.TotFitInt(i,:))));
plot(i,x,'bo');
pixval;
end
title('Integrated Intensity for each cell');
xlabel('cell index');
ylabel('integrated cell intensity (fit)');

figure();
hold on;
for i=1:mm
x=mean(nonzeros((PARtmr_mid.TotFitInt(i,:))));plot(i,x,'bo');
end
title('Integrated Intensity for each cell');
xlabel('cell index');
ylabel('integrated cell intensity (fit)');



figure();
hold on;       
for i=1:mm
    len = nnz(PARtmr_mid.xgw(i,:))
    i
    if len ~= 0
        for z = 1:len
        x = PARtmr_mid.xgw(i,z);
        y = PARtmr_mid.ygw(i,z);
        foo(i,z) = pi()*x./2*y./2;
        end
    end
end

plot(foo(:,:),'bo')
title('Area of each fit per cell');
xlabel('Cell index');
ylabel('area of gaussian fit');

figure();
hold on;       
for i=1:mm
    len = nnz(PARtmr_mid.r(i,:))
    i
    if len ~= 0
        for z = 1:len
        foo(i,z) = PARtmr_mid.r(i,z);
        end
    end
end

plot(foo(:,:),'bo')
title('elongation of axis of each fit per cell');
xlabel('Cell index');
ylabel('PAR.r');

figure();
hold on;
for i=1:mm
    len = nnz(PARtmr_mid.zinit(i,:));
    if len ~= 0
        for z = 1:len
        foo(i,z)= abs(PARtmr_mid.zinit(i,z)-PARtmr_mid.zabs(i,z));
        end
    end
end
plot(foo(:,:),'bo')
title('absolute diffence in z position initial and fitted');
xlabel('Cell index');
ylabel('z-stack');


figure();
hold on;
for i=1:mm
    len(i,1) = nnz(PARtmr_mid.nucRNA(i,:));
    len(i,2) = nnz(PARtmr_mid.cytoRNA(i,:));
end
bar(len);
len = nnz(PARtmr_mid.nucRNA(:,:));
title('Nuclear vs. Cytoplasmic RNA per cell');
xlabel('Cell index');
ylabel('frequency');
legend('Nuclear RNA','Cytoplasmic');

%for cy5

figure(12);
hold on;
for i=1:mm
x = nonzeros((PARcy5_mid.TotFitInt(i,:)));
vect= zeros(size(x))+i;
plot(vect,x,'bo'); impixelinfo;
end
title('Integrated Intensity for each cell (CY5)');

xlabel('cell index');
ylabel('integrated cell intensity (fit)');



figure();
hold on;       
for i=1:mm
    len = nnz(PARcy5_mid.xgw(i,:))
    i
    if len ~= 0
        for z = 1:len
        x = PARcy5_mid.xgw(i,z);
        y = PARcy5_mid.ygw(i,z);
        foo(i,z) = pi()*x./2*y./2;
        end
    end
end

plot(foo(:,:),'bo')
title('Area of each fit per cell (CY5)');
xlabel('Cell index');
ylabel('area of gaussian fit');

figure();
hold on;       
for i=1:mm
    len = nnz(PARcy5_mid.r(i,:))
    i
    if len ~= 0
        for z = 1:len
        foo(i,z) = PARcy5_mid.r(i,z);
        end
    end
end

plot(foo(:,:),'bo')
title('elongation of axis of each fit per cell (CY5)');
xlabel('Cell index');
ylabel('PAR.r');

figure();
hold on;
for i=1:mm
    len = nnz(PARcy5_mid.zinit(i,:));
    if len ~= 0
        for z = 1:len
        foo(i,z)= abs(PARcy5_mid.zinit(i,z)-PARcy5_mid.zabs(i,z));
        end
    end
end
plot(foo(:,:),'bo')
title('absolute diffence in z position initial and fitted (CY5)');
xlabel('Cell index');
ylabel('z-stack');


figure();
hold on;
for i=1:mm
    len(i,1) = nnz(PARcy5_mid.nucRNA(i,:));
    len(i,2) = nnz(PARcy5_mid.cytoRNA(i,:));
end
bar(len);
len = nnz(PARcy5_mid.nucRNA(:,:));
title('Nuclear vs. Cytoplasmic RNA per cell (CY5)');
xlabel('Cell index');
ylabel('frequency');
legend('Nuclear RNA','Cytoplasmic');



imshow(Lab,[])
impixelinfo

figure(9);
hold on;
for i=1:mm
    len = nnz(PARtmr_mid.TotExpInt(i,:));
    if len ~= 0
        for z = 1:len
        foo(i,z,1)= PARtmr_mid.TotExpInt(i,z);
        foo(i,z,2)= PARtmr_mid.TotFitInt(i,z);
        plot(foo(i,z,1),foo(i,z,2),'bo');
        end
    end
    len1 = nnz(PARtmr_mid.xinit(i,:));
    if len1 ~= len
        foo(i,z,1)= PARtmr_mid.TotExpInt(i,z);
        foo(i,z,2)= PARtmr_mid.TotFitInt(i,z);
        plot(foo(i,z,1),foo(i,z,2),'ro');
    end
end
title('correlation of fit to experimental integrated intensities (TMR)');
xlabel('experiment');
ylabel('fit');



figure(1);
hold on;
for i=1:mm
    len = nnz(PARcy5_mid.TotExpInt(i,:));
    if len ~= 0
        for z = 1:len
        foo(i,z,1)= PARcy5_mid.TotExpInt(i,z);
        foo(i,z,2)= PARcy5_mid.TotFitInt(i,z);
        plot(foo(i,z,1),foo(i,z,2),'bo');
        end
    end
    len1 = nnz(PARcy5_mid.xinit(i,:));
    if len1 ~= len
        foo(i,z,1)= PARcy5_mid.TotExpInt(i,z);
        foo(i,z,2)= PARcy5_mid.TotFitInt(i,z);
        plot(foo(i,z,1),foo(i,z,2),'ro');
    end
end
cc = corr2(PARcy5_mid.TotExpInt(~isnan(PARcy5_mid.TotExpInt)),PARcy5_mid.TotFitInt(~isnan(PARcy5_mid.TotFitInt)));
title(['correlation of fit to experimental integrated intensities (CY5), cc = ' num2str(cc)]);
xlabel('experiment');
ylabel('fit');

figure();
hold on;
for i=1:mm
    i
x = nonzeros((PARcy5_mid.TotFitInt(i,:)));
if size(x,1) ~= 0
vect = zeros(size(x,1),1);
for j=1:size(vect,1)
vect(j) = vect(j) + i;
end
plot(vect,x,'bo'); 
end
end
title('Integrated Intensity for each cell (cy5)');

xlabel('cell index');
ylabel('integrated cell intensity (fit)');

figure(2);
clf;
boxplot(PARcy5_mid.TotFitInt')
title('Integrated intensity (Fit) (CY5)')
xlabel('cell index')
ylabel('integrated spot intensity (Fit)')
figure(3)
boxplot(PARcy5_mid.TotExpInt')
title('Integrated intensity (Experiment) CY5')
xlabel('cell index')
ylabel('integrated spot intensity (Experiment)')

figure(8);
clf;
boxplot(PARtmr_mid.TotFitInt')
title('Integrated intensity (Fit) TMR')
xlabel('cell index')
ylabel('integrated spot intensity (Fit)')
figure(10)
boxplot(PARtmr_mid.TotExpInt')
title('Integrated intensity (Experiment) TMR')
xlabel('cell index')
ylabel('integrated spot intensity (Experiment)')



figure();
hold on;
for i=1:mm
plot(PARtmr_mid.r(i,:),PARtmr_mid.rFit(i,:),'ro')
end
title('correlation of ellipticity (r) of fit to experiment (TMR)');
xlabel('experiment');
ylabel('fit');

figure();
hold on;
for i=1:mm
plot(PARcy5_mid.r(i,:),PARcy5_mid.rFit(i,:),'ro')
end
title('correlation of ellipticity (r) of fit to experiment (CY5)');
xlabel('experiment');
ylabel('fit');


corr2(PARcy5_mid.r(~isnan(PARcy5_mid.r)),PARcy5_mid.rFit(~isnan(PARcy5_mid.rFit)))
corr2(PARcy5_mid.TotExpInt(~isnan(PARcy5_mid.TotExpInt)),PARcy5_mid.TotFitInt(~isnan(PARcy5_mid.TotFitInt)))
corr2(PARtmr_mid.TotExpInt(~isnan(PARtmr_mid.TotExpInt)),PARtmr_mid.TotFitInt(~isnan(PARtmr_mid.TotFitInt)))
corr2(PARtmr_mid.r(~isnan(PARtmr_mid.r)),PARtmr_mid.rFit(~isnan(PARtmr_mid.rFit)))
