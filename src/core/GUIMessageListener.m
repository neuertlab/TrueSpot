%
%%
classdef GUIMessageListener

    %%
    properties
        logFileHandle = [];
        dialogHandle = [];
        printToCmd = true;
    end

    %%
    methods

        function obj = openLog(obj, logPath)
            if ~isempty(obj.logFileHandle)
                fclose(obj.logFileHandle);
            end
            [logDir, ~, ~] = fileparts(logPath);
            if ~isfolder(logDir)
                mkdir(logDir);
            end
            obj.logFileHandle = fopen(logPath, 'w');
        end

        function obj = closeLog(obj)
            if ~isempty(obj.logFileHandle)
                fclose(obj.logFileHandle);
            end
            obj.logFileHandle = [];
        end

        function showMessage(obj, message)
            %Updates dialog box to show message
            if ~isempty(obj.dialogHandle)
                obj.dialogHandle.updateMessage(message);
            end
        end

        function logMessage(obj, message)
            %Writes message to log stream
            dt = datetime;
            if ~isempty(obj.logFileHandle)
                fprintf(obj.logFileHandle, '[%s]%s\n', dt, message);
            end
            if obj.printToCmd
                fprintf('[%s]%s\n', dt, message);
            end
        end
    end

end