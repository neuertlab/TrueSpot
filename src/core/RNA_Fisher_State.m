%
%%

classdef RNA_Fisher_State
    
    properties
        rna_fisher_gui_mode; %bool
        
    end
    
    methods
        
        function [obj, initres] = initialize(obj, gui_mode)
            if gui_mode
                obj.rna_fisher_gui_mode = true;
            else
                obj.rna_fisher_gui_mode = false;
            end
            
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
        
    end
    
end
