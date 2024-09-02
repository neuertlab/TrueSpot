%
%%
classdef VisInteractive

    %%
    properties
        spotsRun;

        imageCh; %Cache for loaded image channel.
        figHandle;

        useFilt = false; %Toggle between filtered and raw image channel
        maxProj = false; %Toggle between single slice view and max projection

        %Layer visibility
        imageLayerOn = true;
        cellLayerOn = false;
        nucLayerOn = false;
        spotCircleLayerOn = false;
        quantFitLayerOn = false;
        cloudLayerOn = false;

    end

    %%
    methods
    end

    %%
    methods(Static)
    end

end