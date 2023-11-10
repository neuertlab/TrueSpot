clear all

%// Data (example):
X = randn(1,1e5); %// random variables.
Y = randn(1,1e5);

x_axis = -3:.2:3; %// Define edges of bins for x axis. Column vector
y_axis = -3:.2:3; %// Same for y axis

figure(4);clf
tps = {'0d', '6hr', '12hr', '1d', '2d', '5d','9d'};
%//Actual data
for i = 1:7;
    
X = CY5_AF594_TMR_tp{i}(:,1);
Y = CY5_AF594_TMR_tp{i}(:,3);
% X = CY5_AF594_TMR(:,1);
% Y = CY5_AF594_TMR(:,3);

X = X';
Y = Y'

x_axis = 0:50:500
y_axis = 0:50:500

%// Compute corners of 2D-bins:
[x_mesh_upper,y_mesh_upper] = meshgrid(x_axis(2:end),y_axis(2:end));
[x_mesh_lower,y_mesh_lower] = meshgrid(x_axis(1:end-1),y_axis(1:end-1));

%// Compute centers of 1D-bins:
x_centers = (x_axis(2:end)+x_axis(1:end-1))/2;
y_centers = (y_axis(2:end)+y_axis(1:end-1))/2;

%// Compute pdf:
pdf = mean( bsxfun(@le, X(:), x_mesh_upper(:).') ...
    & bsxfun(@gt, X(:), x_mesh_lower(:).') ...
    & bsxfun(@le, Y(:), y_mesh_upper(:).') ...
    & bsxfun(@gt, Y(:), y_mesh_lower(:).') );
pdf = reshape(pdf,length(x_axis)-1,length(y_axis)-1); %// pdf values at the
%// grid points defined by x_centers, y_centers
pdf = pdf ./ (y_mesh_upper-y_mesh_lower) ./ (x_mesh_upper-x_mesh_lower);
%// normalize pdf to unit integral

%// Compute cdf:
cdf = mean( bsxfun(@le, X(:), x_mesh_upper(:).') ...
    & bsxfun(@le, Y(:), y_mesh_upper(:).') );
cdf = reshape(cdf,length(x_axis)-1,length(y_axis)-1);

%// Plot pdf

%%% Plot with diff number on each row
% if i <4
% subplot(2,3,i)
% else
%     subplot(2,4,i+1)
% end
%%%
subplot(3,3,i)
imagesc(x_centers,y_centers,pdf)
axis xy
axis equal
colorbar
title(tps{i})
xlabel('Xist Molecules')
ylabel('Pdk3 Molecules')
xlim([min(x_axis) max(x_axis)])
ylim([min(y_axis) max(y_axis)])
axis square
end
% %// Plot cdf
% figure(5);clf
% imagesc(x_centers,y_centers,cdf)
% axis xy
% axis equal
% colorbar
% title 'cdf'
% xlabel('Xist Molecules')
% ylabel('Pdk3 Molecules')
% xlim([min(x_axis) max(x_axis)])
% ylim([min(y_axis) max(y_axis)])