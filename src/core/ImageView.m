%
%%

classdef ImageView

    methods(Static)

        function fighandle = renderImageView(imgdat, disprange, figno, z)

            if nargin < 2
                disprange = []; %Adjust later
            end
            if nargin < 3
                figno = 1;
            end
            if nargin < 4; z = 0; end

            %Resolve 3D...
            if ndims(imgdat) > 2
                if z > 0
                    imgdat = double(imgdat(:,:,z));
                else
                    imgdat = double(max(imgdat,[],3));
                end
            else
                imgdat = double(imgdat);
            end

            if ~isempty(disprange)
                Lmin = disprange(1);
                Lmax = disprange(2);
            else
                Lmin = min(imgdat(:));
                Lmax = median(imgdat(:)) + round(10 * std(imgdat(:)));
            end

            fighandle = figure(figno);
            clf;
            imshow(imgdat, [Lmin, Lmax]);
            hold on;
            impixelinfo;

        end

    end

end