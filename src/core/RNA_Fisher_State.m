%
%%

classdef RNA_Fisher_State
    
    properties
        rna_fisher_gui_mode; %bool
        tif_verbose; %bool
        
    end
    
    methods
        
        function [obj, initres] = initialize(obj, gui_mode)
            if gui_mode
                obj.rna_fisher_gui_mode = true;
            else
                obj.rna_fisher_gui_mode = false;
            end
            
            obj.tif_verbose = true;
            initres = 0;
        end
        
        function obj = outputMessageLine(obj, message, send_to_gui)
            %Either way, print to console.
            fprintf("%s\n", message);
            
            %If GUI, send message to GUI as well
            if send_to_gui
                if obj.rna_fisher_gui_mode
                    %TODO
                end
            end
        end
        
    end
    
    methods(Static)
        
        function rf_state = getStaticState()
            persistent rna_fisher_state;
            if isempty(rna_fisher_state)
                rna_fisher_state = RNA_Fisher_State;
            end
            rf_state = rna_fisher_state;
        end
        
        function clearStaticState()
            clear RNA_Fisher_State.getStaticState
        end
        
        function setGUIMode(mybool)
            rf_state = RNA_Fisher_State.getStaticState();
            rf_state.rna_fisher_gui_mode = mybool;
        end
        
        function outputMessageLineStatic(message, send_to_gui)
            rf_state = RNA_Fisher_State.getStaticState();
            rf_state.outputMessageLine(message, send_to_gui);
        end
        
    end
    
end
