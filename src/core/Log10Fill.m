%
%%
function output = Log10Fill(input)
    output = [];
    if isempty(input)
        return;
    end
    
    if ismatrix(input) & (size(input,1) > size(input,2))
        input = transpose(input);
    end
    
    if ~isfloat(input)
        input = double(input);
    end
    
    output = log10(input);
    output(output == -Inf) = NaN;
    
    %Slow forloop for now...
    P = size(output,2);
    for i = 1:P
        if ~isnan(output(i)); continue; end
        l = -1;
        r = -1;
        left = NaN;
        right = NaN;
        for j = i-1:-1:1
            if ~isnan(output(j))
                l = j;
                left = output(j);
                break;
            end
        end
        for j = i+1:P
            if ~isnan(output(j))
                r = j;
                right = output(j);
                break;
            end
        end
        
        if isnan(left)
            %No valid left found
            if ~isnan(right)
                %Extend right value to here.
                for j = 1:i
                    output(j) = right;
                end
                continue;
            end
            %Neither are good. Leave as NaN.
            continue;
        end
        
        if isnan(right)
            %Extend to right end.
            for j = i:P
                output(j) = left;
            end
            continue;
        end
        
        %Try to draw a line.
        m = (right - left) / (r - l);
        b = left - (m * l);
        for j = l+1:r-1
            output(j) = (m * j) + b;
        end
        
    end
end