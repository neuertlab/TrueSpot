
function fhandle = MakeRedBlueScatterplot(x_data, y_data, red_data, blue_data, x_label, y_label, figno)
    
%TODO, add jitter if needed?

    if nargin < 7
        figno = round(rand() * 10000);
    end

    %Scale red and blue axes
    red_max = max(red_data, [], 'all');
    blue_max = max(blue_data, [], 'all');

    point_count = size(x_data,2);

    %For coloring, do points one by one through a for loop for now
    fhandle = figure(figno);
    hold on;
    for i = 1:point_count
        red_lvl = red_data(i)/red_max;
        blue_lvl = blue_data(i)/blue_max;
        plot(x_data(i), y_data(i), 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none',...
            'MarkerFaceColor', [red_lvl 0.0 blue_lvl], 'MarkerSize', 5);
    end

    %Labels
    xlabel(x_label);
    ylabel(y_label);

end