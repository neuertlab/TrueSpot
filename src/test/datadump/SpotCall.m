%
%%

classdef SpotCall
    
    properties
        %Image relative, base pixel call
        init_x = 0;
        init_y = 0;
        init_z = 0;
        init_peak = NaN; %Peak intensity
        init_total = NaN; %Total intensity
        
        %Any sort of fitting
        fit_x = NaN;
        fit_y = NaN;
        fit_z = NaN;
        fit_peak = NaN; %Peak intensity
        fit_total = NaN; 
        
        fit_xFWHM = NaN;
        fit_yFWHM = NaN;

        nearest_trues; %SpotCall or vector of SpotCalls
        nearest_true_distance;
        tfcall; %0 true +, 1 false +, -1 false -, 2 unset
    end
    
    methods

        function obj = initializeMe(obj)
            obj.nearest_trues = SpotCall.empty();
            obj.nearest_true_distance = [];
            obj.tfcall = [];
        end
        
        function obj = importFromQuantSpot(obj, rna_spot, x_off, y_off, z_off)
            gf = rna_spot.gauss_fit;
            obj.fit_x = gf.xfit + x_off;
            obj.fit_y = gf.yfit + y_off;
            obj.fit_z = gf.zabs + z_off;
            obj.fit_peak = gf.fitMInt;
            obj.fit_total = gf.TotFitInt;
            obj.fit_xFWHM = gf.xFWHM;
            obj.fit_yFWHM = gf.yFWHM;
            
            obj.init_x = gf.xinit + x_off;
            obj.init_y = gf.yinit + y_off;
            obj.init_z = gf.zinit + z_off;
            obj.init_peak = gf.expMInt;
            obj.init_total = gf.TotExpInt;
        end
        
    end
    
end