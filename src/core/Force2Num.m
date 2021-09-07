%

%%

function number = Force2Num(value)

    number = NaN;
    if isnumeric(value)
        number = value;
    end
    if ischar(value)
        number = str2double(value);
    end

end