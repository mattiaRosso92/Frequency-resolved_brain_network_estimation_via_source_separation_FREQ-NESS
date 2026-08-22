function [frequencyAxis,responses] = computeFilterResponses( ...
        frequencies,fwhm,plotLimits,nPoints)
%COMPUTEFILTERRESPONSES Build filterFGx-style Gaussian filter responses.

if ~isnumeric(frequencies) || ~isreal(frequencies) || ...
        isempty(frequencies) || ~isvector(frequencies) || ...
        any(~isfinite(frequencies)) || any(frequencies <= 0)
    error('FREQNESS:GUI:InvalidFrequencies', ...
        'Filter frequencies must be a positive finite vector.');
end
frequencies = reshape(frequencies,[],1);

if ~isnumeric(fwhm) || ~isreal(fwhm) || ~isvector(fwhm) || ...
        numel(fwhm) ~= numel(frequencies) || ...
        any(~isfinite(fwhm)) || any(fwhm <= 0)
    error('FREQNESS:GUI:InvalidFWHM', ...
        'FWHM must contain one positive finite value per frequency.');
end
fwhm = reshape(fwhm,[],1);

if ~isnumeric(plotLimits) || ~isreal(plotLimits) || ...
        numel(plotLimits) ~= 2 || any(~isfinite(plotLimits)) || ...
        plotLimits(1) < 0 || plotLimits(1) >= plotLimits(2)
    error('FREQNESS:GUI:InvalidFilterPlotLimits', ...
        'Filter plot limits must contain two ascending nonnegative values.');
end
plotLimits = reshape(plotLimits,1,[]);
if any(frequencies < plotLimits(1)) || any(frequencies > plotLimits(2))
    error('FREQNESS:GUI:InvalidFilterPlotLimits', ...
        'Every filter centre must lie inside the plot limits.');
end

if nargin < 4 || isempty(nPoints)
    nPoints = 1200;
end
if ~isnumeric(nPoints) || ~isscalar(nPoints) || ~isfinite(nPoints) || ...
        nPoints < 2 || nPoints ~= round(nPoints)
    error('FREQNESS:GUI:InvalidFilterPlotResolution', ...
        'Filter plot resolution must be an integer greater than one.');
end

frequencyAxis = linspace(plotLimits(1),plotLimits(2),nPoints);
normalizedWidths = fwhm*(2*pi-1)/(4*pi);
offsets = frequencyAxis-frequencies;
responses = exp(-0.5*(offsets./normalizedWidths).^2);
responses = responses./max(responses,[],2);

end
