%
%%

function boolean = Force2Bool(value)
    boolean = FALSE;
    if ~islogical(value)
        if ischar(value)
            boolean = startsWith(value, "true", 'IgnoreCase', true);
        else
            if isnumeric(value)
                if value == 0
                    boolean = FALSE;
                else
                    boolean = TRUE;
                end
            end
        end
    else
        boolean = value;
    end
end