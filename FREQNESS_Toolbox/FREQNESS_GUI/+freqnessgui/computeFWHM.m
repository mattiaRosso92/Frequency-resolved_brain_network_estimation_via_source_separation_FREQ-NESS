function fwhm = computeFWHM(frequencies,fwidth,filterType)
%COMPUTEFWHM Reproduce the backend FWHM vector for GUI preview.

if ~isnumeric(frequencies) || ~isreal(frequencies) || ...
        isempty(frequencies) || ~isvector(frequencies) || ...
        any(~isfinite(frequencies)) || any(frequencies <= 0)
    error('FREQNESS:GUI:InvalidFrequencies', ...
        'Frequencies must be a positive finite numeric vector.');
end
frequencies = reshape(sort(frequencies),1,[]);
nFrequencies = numel(frequencies);
if ~(ischar(filterType) || (isstring(filterType) && isscalar(filterType))) || ...
        ~ismember(lower(char(filterType)),{'logarithmic','linear'})
    error('FREQNESS:GUI:InvalidFilterType', ...
        'Filter progression must be logarithmic or linear.');
end

if isempty(fwidth)
    fwhm = FREQNESS_ComputeFilterWidths(frequencies,filterType);
elseif isscalar(fwidth)
    fwhm = repmat(fwidth,1,nFrequencies);
elseif isvector(fwidth) && numel(fwidth) == nFrequencies
    fwhm = fwidth(:)';
else
    error('FREQNESS:GUI:InvalidFWHM', ...
        'FWHM must be empty, scalar, or contain one value per frequency.');
end

if ~isnumeric(fwhm) || ~isreal(fwhm) || any(~isfinite(fwhm)) || ...
        any(fwhm <= 0)
    error('FREQNESS:GUI:InvalidFWHM', ...
        'FWHM values must be positive and finite.');
end

end
