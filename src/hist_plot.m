bin_num = 10;
bin_size = (max(Xist_Jpx(:,1))-min(Xist_Jpx(:,1)))/bin_num;
binsX = 0:bin_size:max(Xist_Jpx(:,1));
labelsX = binsX(1:size(binsX,2)-1)+bin_size/2;
binsall = zeros(1,size(binsX,2)-1);
for i = 1:size(binsall,2)
    binsall(1,i) = size(find(Xist_Jpx(:,1) >= binsX(i)),1)-size(find(Xist_Jpx(:,1) >= binsX(i+1)),1);
end
binall = binsall/sum(binsall);
figure(24); clf
plot(labelsX,binall,'r','LineWidth',2); hold on
xlabel('Xist Molecules')
ylabel('Probability'