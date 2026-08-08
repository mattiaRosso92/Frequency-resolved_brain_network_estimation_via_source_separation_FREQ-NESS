function analytic_signal = FREQNESS_AnalyticSignal(data)

% ========================================================================
%
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you use this toolbox, please cite:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ-NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved
%  Brain Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%  This function computes the discrete-time analytic signal through an FFT.
%  The operation is performed along the first dimension (time), with one
%  signal per column. It provides the analytic-signal calculation required
%  by FREQNESS without depending on an external Hilbert-transform function.
%
% ------------------------------------------------------------------------
%  INPUT:
%  - data : real numeric time-by-signals vector or matrix
%
%  OUTPUT:
%  - analytic_signal : complex analytic signal with the same size as data
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso & Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 08/08/2026
%
% ========================================================================

if ~isnumeric(data) || ~isreal(data) || isempty(data) || any(~isfinite(data(:)))
    error('data must be a non-empty real numeric vector or matrix containing finite values.');
end

if ndims(data) > 2
    error('data must be a vector or a 2D matrix with time along the first dimension.');
end

input_was_row = isrow(data);
if input_was_row
    data = data(:);
end

ntime = size(data,1);
if ntime < 2
    error('data must contain at least two time samples.');
end

analytic_filter = zeros(ntime,1);
analytic_filter(1) = 1;

if mod(ntime,2) == 0
    analytic_filter(ntime/2+1) = 1;
    analytic_filter(2:ntime/2) = 2;
else
    analytic_filter(2:(ntime+1)/2) = 2;
end

analytic_signal = ifft(fft(data,[],1).*analytic_filter,[],1);

if input_was_row
    analytic_signal = analytic_signal.';
end

end
