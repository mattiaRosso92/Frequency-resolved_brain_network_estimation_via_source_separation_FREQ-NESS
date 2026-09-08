function fwidth = FREQNESS_ComputeFilterWidths(frex,filter)

% ========================================================================
%
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you use this toolbox, please cite:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain
%  Networks During Auditory Stimulation. Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%  This function computes physically calibrated Gaussian-filter widths
%  from the analyzed frequencies. The schedule depends on frequency in Hz
%  and its numerical range, not on the number or spacing of frequency bins.
%
%  The logarithmic schedule (default) decreases spectral selectivity
%  smoothly from Q = 7 at the lowest frequency to Q = 3.5 at the highest:
%
%      t(f) = log(f/fmin) / log(fmax/fmin)
%      Q(f) = 7 * (3.5/7)^t(f)
%      FWHM(f) = f / Q(f)
%
%  The linear schedule uses constant Q = 5, so FWHM grows linearly with
%  centre frequency. For a single unique frequency, the logarithmic mode
%  uses the geometric midpoint Q = sqrt(7*3.5).
%
% ------------------------------------------------------------------------
%  INPUTS:
% ------------------------------------------------------------------------
%  - frex   : Non-empty vector of positive finite frequencies in Hz.
%  - filter : 'logarithmic' (default) or 'linear'.
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%  - fwidth : Row vector containing one positive FWHM value in Hz for each
%             input frequency, in the same element order as frex.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso & Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%
% ========================================================================

if nargin < 1 || ~isnumeric(frex) || ~isreal(frex) || ...
        ~isvector(frex) || isempty(frex) || any(~isfinite(frex)) || ...
        any(frex <= 0)
    error('frex must be a non-empty vector of positive finite frequencies.');
end

if nargin < 2 || isempty(filter)
    filter = 'logarithmic';
end

if ~(ischar(filter) || (isstring(filter) && isscalar(filter))) || ...
        ~(strcmpi(filter,'logarithmic') || strcmpi(filter,'linear'))
    error('filter must be either ''logarithmic'' or ''linear''.');
end

frex = double(frex(:)');

if strcmpi(filter,'linear')
    q_values = 5 * ones(size(frex));
else
    q_low = 7;
    q_high = 3.5;
    fmin = min(frex);
    fmax = max(frex);

    if fmin == fmax
        q_values = sqrt(q_low*q_high) * ones(size(frex));
    else
        log_position = log(frex./fmin) ./ log(fmax/fmin);
        q_values = q_low .* (q_high/q_low).^log_position;
    end
end

fwidth = frex ./ q_values;

end
