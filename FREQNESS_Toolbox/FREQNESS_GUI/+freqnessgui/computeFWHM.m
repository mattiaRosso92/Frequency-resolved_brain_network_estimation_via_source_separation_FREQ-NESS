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
    referenceFrequency = 2.439;
    referenceFWHM = 0.35;
    nAbove = 80;
    nBelow = 6;
    frequenciesAbove = linspace(referenceFrequency, ...
        referenceFrequency*(nAbove/2),nAbove);
    frequenciesBelow = linspace(referenceFrequency/2, ...
        referenceFrequency/(2*nBelow),nBelow);
    referenceFrequencies = [frequenciesBelow(end:-1:1),frequenciesAbove];
    fwhmAbove = logspace(log10(referenceFWHM), ...
        log10(referenceFWHM*nAbove),nAbove);
    fwhmBelow = logspace(log10(referenceFWHM), ...
        log10(referenceFWHM/nBelow),nBelow+1);
    referenceWidths = [fwhmBelow(end:-1:2),fwhmAbove];
    [~,closestIndex] = min(abs(referenceFrequencies-frequencies(1)));
    initialWidth = referenceWidths(closestIndex);

    if strcmpi(filterType,'logarithmic')
        fwhm = logspace(log10(initialWidth), ...
            log10(initialWidth*nFrequencies),nFrequencies);
    elseif strcmpi(filterType,'linear')
        fwhm = linspace(initialWidth, ...
            initialWidth*nFrequencies,nFrequencies);
    end
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
