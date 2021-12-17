%
%%

function boolean = Force2Bool(value)
    boolean = false;
    if ~islogical(value)
        if ischar(value)
            boolean = startsWith(value, "true", 'IgnoreCase', true);
        else
            if isnumeric(value)
                if value == 0
                    boolean = false;
                else
                    boolean = true;
                end
            end
        end
    else
        boolean = value;
    end
end