%
%%
classdef GUIUtils

    methods(Static)

        %%------------- Logging/Messaging -------------
        function sendMessage(listener, origin, message)
            if ~isempty(listener)
                listener.showMessage(message);
                listener.logMessage(['[' origin ']' message]);
            end
        end

    end

end