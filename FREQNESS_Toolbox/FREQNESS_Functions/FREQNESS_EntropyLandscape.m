function [ED,H2] = FREQNESS_EntropyLandscape(FREQ)

% ========================================================================
%
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you find this function useful, please cite the first FREQNESS paper:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain
%  Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%  This function estimates effective dimensionality (ED) and quadratic  
%  Rényi entropy (H2) across frequencies, given the eigenspectrum of a 
%  covariance matrix. It uses the entropy-based index derived by Pirk et al. 
%  (2012) to estimate the effective number of uncorrelated measurements, as  
%  described in Del Giudice (2020).
%
%  It can be applied to a 1-D vector of eigenvalues, as well as directly to
%  the FREQ.evals output produced by FREQNESS_NetworkEstimation().
%  If FREQ.evals is given as an input, this can consist of either a 2D or
%  3D matrix, depending on whether FREQNESS was run on individual
%  participants or at the group level. 
%
%  When FREQ.evals is provided as input, the ED output will reflect the
%  effective number of uncorrelated measurements across the frequency
%  spectrum. If FREQ.evals contains participants as a 3rd dimension, then ED
%  will also include an additional dimension for multiple participants.
%  This will allow to visualize the grand-average entropy landscape and
%  eventually carry out statistical testing.
% ========================================================================

% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - FREQ : structure with fields
%           • FREQ.evals -> eigenvalue array [nComp x nFrex x (nSubs)]
%           • FREQ.frex  -> frequency vector [nFrex x 1] (optional but
%                            used for the x-axis in the visualization)
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - ED    : Scalar or array containing the estimate of effective
%            dimensionality.
%
%  - H2    : Quadratic Rényi entropy of the eigenvalue distribution. 
%            Provides a log-scaled measure of the evenness of variance 
%            across components, without conversion to an equivalent 
%            number of dimensions.
%
% ------------------------------------------------------------------------
%  REFERENCES:
% ------------------------------------------------------------------------
%  - Del Giudice, M. (2021). Effective dimensionality: A tutorial.
%    Multivariate Behavioral Research, 56(3), 527–542.
%
%  - Pirk, R. J., Remley, K. A., & Patané, C. S. L. (2012).
%    Reverberation chamber measurement correlation.
%    IEEE Transactions on Electromagnetic Compatibility, 54(3), 533–545.
%
%  - Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%    Kringelbach, M. L., & Bonetti, L. (2025).
%    FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain
%    Networks During Auditory Stimulation. Advanced Science, 2413195.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso & Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 30/09/2025
%
% ========================================================================


%% Map inputs from FREQ 

if ~isstruct(FREQ)
    error('Input must be a struct with fields FREQ.evals (and optionally FREQ.frex).');
end
if ~isfield(FREQ,'evals') || isempty(FREQ.evals)
    error('FREQ.evals is missing or empty.');
end
eigenspectrum = FREQ.evals;

% Optional frequency axis for plotting
if isfield(FREQ,'frex') && ~isempty(FREQ.frex)
    frex = FREQ.frex(:);
end

%% Normalize input shape to [nComp x nFrex x nSubs]

if isvector(eigenspectrum)
    eigenspectrum = eigenspectrum(:);               % [nComp x 1]
end

% Get the number of dimensions in the input
nd = ndims(eigenspectrum);
if nd == 2
    % [nComp x nFrex] -> [nComp x nFrex x 1]
    eigenspectrum = reshape(eigenspectrum, size(eigenspectrum,1), size(eigenspectrum,2), 1);
elseif nd > 3
    error('eigenspectrum must be 1D, 2D, or 3D (nComp x nFrex x nSubs).');
end

% Define Ns
ncomps = size(eigenspectrum,1);
nfrex  = size(eigenspectrum,2);
nsubs  = size(eigenspectrum,3);

%% Compute effective dimensionality (ED)

num   = sum(eigenspectrum, 1).^2;       % [1 x nFrex x nSubs]
denom = sum(eigenspectrum.^2, 1);       % [1 x nFrex x nSubs]

% Guard against division by zero
denom = max(denom, realmin);

% Compute (sum(lambda))^2 / sum(lambda.^2) along component dimension
ED = squeeze(num ./ denom);             % [nFrex x nSubs]

%% Compute the quadratic Rényi entropy

% Normalize across components (dim 1) to get probabilities p
sumtemp = sum(eigenspectrum,1);                         
sumtemp = max(sumtemp, realmin);                         
p = eigenspectrum ./ sumtemp;                           

% Quadratic Rényi entropy H2 = -log(sum p_i^2)
H2 = squeeze( -log( sum(p.^2, 1) ) );


%% Visualize entropy landscapes: ED and H2 

% X-axis: use 'frex' if available, otherwise frequency index
if exist('frex','var') && numel(frex)==nfrex
    x = frex(:);
    xlab = 'Frequency';
else
    x = (1:nfrex).';
    xlab = 'Frequency index';
end

% Mean ± SEM across subjects
if nsubs > 1
    avg_ED  = mean(ED,  2, 'omitnan');
    avg_H2  = mean(H2,  2, 'omitnan');
    sem_ED  = std(ED,  0, 2, 'omitnan') ./ sqrt(nsubs);
    sem_H2  = std(H2,  0, 2, 'omitnan') ./ sqrt(nsubs);
end

figure('Color','w', 'Units','normalized', 'Position',[0.3 0.3 0.35 0.5]);

% --- ED ---
subplot(2,1,1); hold on
if nsubs > 1
    fill([x; flipud(x)], [avg_ED - sem_ED; flipud(avg_ED + sem_ED)], ...
        [0 0 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
plot(x, ED, 'k', 'LineWidth', 2);
ylabel('ED', 'FontSize', 13, 'FontWeight', 'bold');
title('Effective Dimensionality', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'off');
grid on; grid minor; xlim([min(x) max(x)]);

% --- H2 ---
subplot(2,1,2); hold on
if nsubs > 1
    fill([x; flipud(x)], [avg_H2 - sem_H2; flipud(avg_H2 + sem_H2)], ...
        [0 0 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
plot(x, H2, 'k', 'LineWidth', 2);
xlabel(xlab, 'FontSize', 13, 'FontWeight', 'bold');
ylabel('H_2', 'FontSize', 13, 'FontWeight', 'bold');
title('Quadratic Rényi Entropy', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'off');
grid on; grid minor; xlim([min(x) max(x)]);


%% Adjust outputs 

% If single frequency & single subject, return scalar
if isempty(ED)
    ED = NaN;
elseif isscalar(ED)
    % leave as scalar
else
    % ensure 2D output [nFrex x nSubs]
    ED = reshape(ED, nfrex, nsubs);
end

% If single frequency & single subject, return scalar
if isempty(H2)
    H2 = NaN;
elseif isscalar(H2)
    % leave as scalar
else
    % ensure 2D output [nFrex x nSubs]
    H2 = reshape(H2, nfrex, nsubs);
end

end
