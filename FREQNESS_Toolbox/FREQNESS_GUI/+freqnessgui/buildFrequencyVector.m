function [frequencies,snappedRange] = buildFrequencyVector(requestedRange,step)
%BUILDFREQUENCYVECTOR Snap two endpoints to a fixed positive frequency grid.

if ~isnumeric(requestedRange) || ~isreal(requestedRange) || ...
        numel(requestedRange) ~= 2 || any(~isfinite(requestedRange))
    error('FREQNESS:GUI:InvalidFrequencyRange', ...
        'The frequency range must contain two finite numeric endpoints.');
end
if ~isnumeric(step) || ~isreal(step) || ~isscalar(step) || ...
        ~isfinite(step) || step <= 0
    error('FREQNESS:GUI:InvalidFrequencyStep', ...
        'Frequency resolution must be one positive finite value.');
end

requestedRange = reshape(sort(requestedRange),1,[]);
snappedRange = round(requestedRange./step).*step;
snappedRange(1) = max(step,snappedRange(1));
snappedRange(2) = max(snappedRange(1),snappedRange(2));

% Suppress floating-point tails while retaining the user-selected resolution.
decimalScale = 10^max(0,min(12,ceil(-log10(step))+10));
snappedRange = round(snappedRange.*decimalScale)./decimalScale;
frequencies = snappedRange(1):step:snappedRange(2);
frequencies = round(frequencies.*decimalScale)./decimalScale;
if isempty(frequencies) || frequencies(end) < snappedRange(2)
    frequencies(end+1) = snappedRange(2);
end

end
