%
%%

classdef RNAThreshold
    
    methods (Static)
        
        function threshold_results = runDefaultParameters(spot_count_table, ctrl_count_table)
            param_struct = RNA_Threshold_Common.genEmptyThresholdParamStruct();
            param_struct.sample_spot_table = spot_count_table;
            
            if nargin > 1
                param_struct.control_spot_table = ctrl_count_table;
            end
            
            threshold_results = RNA_Threshold_Common.estimateThreshold(param_struct);
        end
        
        function threshold_results = runSavedParameters(rnaspots_run, verbosity, spot_table, ctrl_table)
            if nargin < 2
                verbosity = 0;
            end
            
            param_struct = RNA_Threshold_Common.genEmptyThresholdParamStruct();
            
            if nargin >= 3 & ~isempty(spot_table)
                param_struct.sample_spot_table = spot_table;
            else
                [rnaspots_run, param_struct.sample_spot_table, ~] = rnaspots_run.loadZTrimmedTables_Sample();
            end
            
            if nargin >= 4 & ~isempty(ctrl_table)
                param_struct.control_spot_table = ctrl_table;
            else
                [rnaspots_run, param_struct.control_spot_table, ~] = rnaspots_run.loadZTrimmedTables_Control();
            end
            
            param_struct.verbosity = verbosity;
            
            winmin = rnaspots_run.ttune_winsz_min;
            if winmin < 2; winmin = 2; end
            winincr = rnaspots_run.ttune_winsz_incr;
            if winincr < 1; winincr = 1; end
            winmax = rnaspots_run.ttune_winsz_max;
            if winmax < winmin; winmax = winmin+winincr; end
            param_struct.window_sizes = [winmin:winincr:winmax];
            
            madmin = rnaspots_run.ttune_madf_min;
            if isnan(madmin); madmin = -1.0; end
            madmax = rnaspots_run.ttune_madf_max;
            if isnan(madmax); madmax = madmin+0.5; end
            param_struct.mad_factor_min = madmin;
            param_struct.mad_factor_max = madmax;
            
            spline_itr = rnaspots_run.ttune_spline_itr;
            if spline_itr < 0; spline_itr = 0; end
            param_struct.spline_iterations = spline_itr;
            
            threshold_results = RNA_Threshold_Common.estimateThreshold(param_struct);
        end
        
    end
    
end