function [indices,selectedFrequencies] = mapFrequencySelection(availableFrequencies,requested)
%MAPFREQUENCYSELECTION Map scalar/range requests to the imported FREQ grid.

if ~isnumeric(availableFrequencies) || ~isreal(availableFrequencies) || ...
        isempty(availableFrequencies) || ~isvector(availableFrequencies) || ...
        any(~isfinite(availableFrequencies)) || any(availableFrequencies <= 0)
    error('FREQNESS:GUI:InvalidAvailableFrequencies', ...
        'Available frequencies must be a positive finite vector.');
end
availableFrequencies = availableFrequencies(:)';
if any(diff(availableFrequencies) <= 0)
    error('FREQNESS:GUI:InvalidAvailableFrequencies', ...
        'Available frequencies must be strictly increasing.');
end
if ~isnumeric(requested) || ~isreal(requested) || ...
        ~ismember(numel(requested),[1 2]) || any(~isfinite(requested))
    error('FREQNESS:GUI:InvalidFrequencySelection', ...
        'A frequency selection must contain one value or two range endpoints.');
end

requested = reshape(requested,1,[]);
requested = sort(requested);
indices = zeros(size(requested));
for requesti = 1:numel(requested)
    [~,indices(requesti)] = min(abs(availableFrequencies-requested(requesti)));
end
if isscalar(indices)
    selectedFrequencies = availableFrequencies(indices);
else
    indices = sort(indices);
    selectedFrequencies = availableFrequencies(indices(1):indices(2));
end

end
