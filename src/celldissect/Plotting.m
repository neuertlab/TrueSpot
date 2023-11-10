xsx = [0];
ysy = [0];
all1 = [0];
counter123 = 1;
for i = 1:size(PAR.TotExpInt,1)
    for j = 1:size(PAR.TotExpInt,2)
        if not(isnan(PAR.TotExpInt(i,j)))
            all1(counter123) = PAR.TotExpInt(i,j);
            counter123 = counter123+1;
        end
    end
end
all2 = log10(all1);
edges = min(all2):(max(all2)-min(all2))/30:max(all2);
%xsx = min(all1):1000:max(all1)
figure(20);clf;plot(edges,(hist(all2,edges)/sum(all2)))
    xlabel('log10(Fluorescence)')
    ylabel('Probability')
    title('Fluorescence intensity per spot')

d1hist = hist(datas{j},edges);%create histogram of data
    d1Nhist = d1hist/sum(d1hist(:));%normalize histogram
    [F,goF] = fit(edges',d1Nhist',type,specs);%fit to function
    variables(1,j) = F.mu;
    variables(2,j) = F.sigma;
    variables(3,j) = F.lambda;
    variables(4,j) = binnum;
    variables(5,j) = goF.rsquare
    figure(j)
    clf;
    %plot(edges,d1Nhist)
    bar(edges, d1Nhist, 'Facecolor','w')
    hold on
    edges = 1:.2:(max(datas{j})+10);
    d1plot = zeros(size(edges,2),1);
    for i = 1:size(edges,2)
        d1plot(i,1) = F.lambda/2*exp(F.lambda/2*(2*F.mu+F.lambda*(F.sigma)^2-2*edges(1,i)))*(erfc((F.mu+F.lambda*F.sigma^2-edges(1,i))/(sqrt(2)*F.sigma)));
    end
    plot(edges,d1plot);
    title(labels{j});
    TeXString = texlabel(strcat('rsquare:',num2str(goF.rsquare)));
    text(50,0.045, TeXString)
    TeXString=texlabel(strcat('lambda:',num2str(F.lambda)));
    text(50,0.035, TeXString)
    TeXString=texlabel(strcat('mu:',num2str(F.mu)));
    text(50,0.025, TeXString)
    TeXString=texlabel(strcat('sigma:',num2str(F.sigma)));
    text(50,0.015, TeXString)
    legend('Data','EMG Fit','NorthEast')
    xlabel('Intermitotic Time')
    ylabel('Probability')
