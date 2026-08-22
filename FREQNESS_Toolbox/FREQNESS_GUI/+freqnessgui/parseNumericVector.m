function values = parseNumericVector(textValue,allowEmpty,fieldName)
%PARSENUMERICVECTOR Parse numbers and colon ranges without evaluating code.
%
% Supported examples: "1 2 3", "1, 2, 3", "1:5", and "1:0.5:5".

if nargin < 2
    allowEmpty = false;
end
if nargin < 3
    fieldName = 'Value';
end

if isnumeric(textValue)
    values = textValue(:)';
    return
end

if ~(ischar(textValue) || (isstring(textValue) && isscalar(textValue)))
    error('FREQNESS:GUI:InvalidNumericText','%s must be numeric text.',fieldName);
end

textValue = strtrim(char(textValue));
if strcmp(textValue,'[]')
    textValue = '';
end
if numel(textValue) >= 2 && textValue(1) == '[' && textValue(end) == ']'
    textValue = strtrim(textValue(2:end-1));
end

if isempty(textValue)
    if allowEmpty
        values = [];
        return
    end
    error('FREQNESS:GUI:MissingNumericValue','%s cannot be empty.',fieldName);
end

tokens = regexp(textValue,'[,;\s]+','split');
tokens = tokens(~cellfun('isempty',tokens));
values = [];

for tokeni = 1:numel(tokens)
    token = tokens{tokeni};
    colonParts = strsplit(token,':');
    numbers = cellfun(@str2double,colonParts);
    if any(~isfinite(numbers))
        error('FREQNESS:GUI:InvalidNumericText', ...
            '%s contains invalid numeric text: "%s".',fieldName,token);
    end

    if isscalar(numbers)
        expanded = numbers;
    elseif numel(numbers) == 2
        expanded = numbers(1):numbers(2);
    elseif numel(numbers) == 3
        if numbers(2) == 0
            error('FREQNESS:GUI:InvalidNumericText', ...
                '%s contains a zero range step.',fieldName);
        end
        expanded = numbers(1):numbers(2):numbers(3);
    else
        error('FREQNESS:GUI:InvalidNumericText', ...
            '%s contains an invalid colon expression: "%s".',fieldName,token);
    end

    if isempty(expanded)
        error('FREQNESS:GUI:InvalidNumericText', ...
            '%s contains an empty range: "%s".',fieldName,token);
    end
    values = [values expanded]; %#ok<AGROW>
end

end
