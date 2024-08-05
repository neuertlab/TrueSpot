%%
classdef VisCommon
    
    %%
    properties
        colortbl_red;
        colortbl_magenta;
        colortbl_yellow;
        colortbl_cyan;
        colortbl_white;

        colorBelowFar;
        colorBelowNear;
        colorAboveFar;
        colorAboveNear;

        idims;

        isInitialized = false;
    end

    %%
    methods
        
        %%
        function obj = initializeMe(obj, idims)
            if nargin < 2; idims = []; end

            if isempty(idims)
                idims = struct('x', 0, 'y', 0, 'z', 0);
            end

            obj = obj.generateColorTables(idims.z + 1);
            obj.idims = idims;

            obj.colorBelowFar = [0.60, 0.50, 0.0];
            obj.colorBelowNear = [1.0, 1.0, 0.0];
            obj.colorAboveFar = [0.0, 0.50, 0.20];
            obj.colorAboveNear = [0.0, 1.0, 0.0];

            obj.isInitialized = true;
        end

        %%
        function obj = generateColorTables(obj, levels)
            obj.colortbl_red = VisCommon.generateDarkColors(levels, [1.0, 0.0, 0.0]);
            obj.colortbl_magenta = VisCommon.generateDarkColors(levels, [1.0, 0.0, 1.0]);
            obj.colortbl_yellow = VisCommon.generateDarkColors(levels, [1.0, 1.0, 0.0]);
            obj.colortbl_cyan = VisCommon.generateDarkColors(levels, [0.0, 1.0, 1.0]);
            obj.colortbl_white = VisCommon.generateDarkColors(levels, [1.0, 1.0, 1.0]);
        end


    end

    %%
    methods(Static)

         %%
        function clrmtx = gen_z_color_matrix(color_tbl, coord_tbl, rel_z, range)
            maxidx = size(color_tbl,1);
            incr = 1;
            
            if range > 0
                incr = maxidx./range;
            end
            
            color_idxs = min((abs(rel_z - coord_tbl(:,3)) + 1) * incr, maxidx);
            clrmtx = color_tbl(color_idxs,:);
        end

        %%
        function color_tbl = generateDarkColors(levels, base_color)

            levels = double(levels+1);
            color_tbl = NaN(levels,3);
            konst = 1/log10(levels);

            r = base_color(1,1);
            g = base_color(1,2);
            b = base_color(1,3);

            color_tbl(1,1) = r;
            color_tbl(1,2) = g;
            color_tbl(1,3) = b;

            %r_div = r/levels;
            %g_div = g/levels;
            %b_div = b/levels;

            for l = 2:levels
                %r = r-r_div;
                r = (1 - (log10(l) * konst)) * base_color(1,1);
                color_tbl(l,1) = r;

                %g = g-g_div;
                g = (1 - (log10(l) * konst)) * base_color(1,2);
                color_tbl(l,2) = g;

                %b = b-b_div;
                b = (1 - (log10(l) * konst)) * base_color(1,3);
                color_tbl(l,3) = b;
            end

        end

        %%
        function color_tbl = generateLightColors(levels, base_color)

            levels = double(levels+1);
            color_tbl = NaN(levels,3);
            konst = 1/log10(levels);

            r = base_color(1,1);
            g = base_color(1,2);
            b = base_color(1,3);

            color_tbl(1,1) = r;
            color_tbl(1,2) = g;
            color_tbl(1,3) = b;

            for l = 2:levels
                %r = r-r_div;
                r = (log10(l) * konst) + base_color(1,1);
                color_tbl(l,1) = min(r, 1.0);

                %g = g-g_div;
                g = (log10(l) * konst) + base_color(1,2);
                color_tbl(l,2) = min(g, 1.0);

                %b = b-b_div;
                b = (log10(l) * konst) + base_color(1,3);
                color_tbl(l,3) = min(b, 1.0);
            end

        end

    end
end