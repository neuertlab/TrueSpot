
function fhandle = MakeRedBlueScatterplot(x_data, y_data, red_data, blue_data, x_label, y_label, figno)
    
%TODO, add jitter if needed?

    if nargin < 7
        figno = round(rand() * 10000);
    end

    %Remove cases where any values are NaN
    is_bad = isnan(x_data) | isnan(y_data) | isnan(red_data) | isnan(blue_data);
    keepidx = find(~is_bad);
    x_data = x_data(keepidx);
    y_data = y_data(keepidx);
    red_data = red_data(keepidx);
    blue_data = blue_data(keepidx);


    %Scale red and blue axes
    red_max = max(red_data, [], 'all');
    blue_max = max(blue_data, [], 'all');

    rdim = size(x_data, 1);
    cdim = size(x_data, 2);
    if rdim > cdim
        x_data = x_data.';
        y_data = y_data.';
        red_data = red_data.';
        blue_data = blue_data.';
    end
    
    point_count = size(x_data,2);

    %For coloring, do points one by one through a for loop for now
    fhandle = figure(figno);
    hold on;
    for i = 1:point_count
        red_lvl = red_data(i)/red_max;
        blue_lvl = blue_data(i)/blue_max;
        plot(x_data(i), y_data(i), 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none',...
            'MarkerFaceColor', [red_lvl 0.0 blue_lvl], 'MarkerSize', 5);
        hold on;
    end

    %Labels
    xlabel(x_label);
    ylabel(y_label);

end