%
%%

function fig_handle = GenMultiTool_SpotPlot(img_results_obj, figno)

if nargin < 2
    figno = round(rand() * 10000);
end

fig_handle = figure(figno);
subplot(2,2,1);
fig_handle_temp = img_results_obj.renderSpotCountPlot('homebrew', [1.0 0.0 0.0], figno, fig_handle);
fig_handle = handleReturnHandle(fig_handle, fig_handle_temp);
title('Neuert Lab');

subplot(2,2,2);
fig_handle_temp = img_results_obj.renderSpotCountPlot('bigfish', [0.0 0.0 1.0], figno, fig_handle);
fig_handle = handleReturnHandle(fig_handle, fig_handle_temp);
title('BigFISH');

subplot(2,2,3);
fig_handle_temp = img_results_obj.renderSpotCountPlot('rsfish', [0.0 0.8 0.0], figno, fig_handle);
fig_handle = handleReturnHandle(fig_handle, fig_handle_temp);
title('RS-FISH');

subplot(2,2,4);
fig_handle_temp = img_results_obj.renderSpotCountPlot('deepblink', [0.7 0.7 0.0], figno, fig_handle);
fig_handle = handleReturnHandle(fig_handle, fig_handle_temp);
title('DeepBlink');


end

function fighandle_out = handleReturnHandle(fighandle_og, fighandle_ret)
    if isempty(fighandle_ret)
        fighandle_out = fighandle_og;
    else
        fighandle_out = fighandle_ret;
    end
end