function [H2,ED] = FREQNESS_EntropyLandscape(FREQ,varargin)

% ========================================================================
%
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you use this toolbox, please cite:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain
%  Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%  This function estimates quadratic Rényi entropy (H2) and effective
%  dimensionality (ED) across frequencies from the eigenspectrum of a
%  covariance matrix. It uses the entropy-based index derived by Pirk et al.
%  (2012) to estimate the effective number of uncorrelated measurements, as
%  described in Del Giudice (2021).
%
%  Both H2 and ED are returned for numerical analyses. The visualization
%  deliberately shows only the H2 entropy landscape.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - FREQ : structure with fields
%           • FREQ.evals -> complete eigenvalue array
%                          [nEigenvalues x nFrex x (nSubs)]
%           • FREQ.frex  -> frequency vector [nFrex x 1] (optional)
%
%  - Optional name-value pair:
%      - 'plot_all' : when true, additionally plot the H2 landscape for
%                     every participant. Default: false.
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - H2 : Quadratic Rényi entropy of the eigenvalue distribution.
%
%  - ED : Effective dimensionality. This output is retained for numerical
%         analysis and backwards compatibility but is not visualized.
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


%% Handle optional arguments

opts = struct('plot_all',false);
opts = parse_name_value_pairs(opts,varargin{:});
plot_all = opts.plot_all;

if ~islogical(plot_all) || ~isscalar(plot_all)
    error('plot_all must be one logical value (true or false).');
end


%% Map and validate inputs from FREQ

if ~isstruct(FREQ)
    error('Input must be a FREQ structure.');
end

if ~isfield(FREQ,'evals') || isempty(FREQ.evals)
    error('FREQ.evals is missing or empty.');
end

eigenspectrum = FREQ.evals;
if ~isnumeric(eigenspectrum) || ~isreal(eigenspectrum) || ...
        ndims(eigenspectrum) > 3 || any(isinf(eigenspectrum(:))) || ...
        any(eigenspectrum(:) < 0)
    error(['FREQ.evals must be a real numeric array of non-negative, ' ...
        'non-infinite eigenvalues in ' ...
        '[nEigenvalues x nFrex x (nSubs)] format.']);
end

if isvector(eigenspectrum)
    eigenspectrum = eigenspectrum(:);
end

nfrex = size(eigenspectrum,2);
nsubs = size(eigenspectrum,3);

% Use the analyzed frequencies when valid; otherwise use frequency indices
if isfield(FREQ,'frex') && isnumeric(FREQ.frex) && isreal(FREQ.frex) && ...
        isvector(FREQ.frex) && numel(FREQ.frex)==nfrex && ...
        all(isfinite(FREQ.frex))
    x = FREQ.frex(:);
    xlab = 'Frequency (Hz)';
else
    x = (1:nfrex)';
    xlab = 'Frequency index';
end

fprintf(['\nFREQNESS Entropy Landscape: computing Rényi entropy (H2) ' ...
    'and effective dimensionality (ED) for %d participants.\n'],nsubs);


%% Compute quadratic Rényi entropy and effective dimensionality

sumtemp = sum(eigenspectrum,1);
sumtemp = max(sumtemp,realmin);
p = eigenspectrum./sumtemp;

% Quadratic Rényi entropy: H2 = -log(sum(p_i^2))
H2 = reshape(-log(sum(p.^2,1)),nfrex,nsubs);

% Effective dimensionality: ED = sum(lambda)^2/sum(lambda^2)
num = sum(eigenspectrum,1).^2;
denom = max(sum(eigenspectrum.^2,1),realmin);
ED = reshape(num./denom,nfrex,nsubs);


%% Visualize H2 entropy landscape

avg_H2 = mean(H2,2,'omitnan');
if nsubs > 1
    sem_H2 = std(H2,0,2,'omitnan')/sqrt(nsubs);
else
    sem_H2 = [];
end

cmap = parula(256);
col_line = cmap(1,:);

figure('Color','w','Units','normalized','Position',[0.3 0.3 0.35 0.5]);
hold on
if nsubs > 1 && ~all(isnan(sem_H2))
    fill([x; flipud(x)],[avg_H2-sem_H2; flipud(avg_H2+sem_H2)], ...
        [0 0 0],'FaceAlpha',0.1,'EdgeColor','none');
end
plot(x,avg_H2,'Color',col_line,'LineWidth',2);
xlabel(xlab,'FontSize',13,'FontWeight','bold');
ylabel('H_2','FontSize',13,'FontWeight','bold');
title('Quadratic Rényi Entropy - Grand Average', ...
    'FontSize',14,'FontWeight','bold');
set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
grid on
grid minor
set_frequency_limits(x)


%% Optional individual-participant plots

if plot_all
    for subi = 1:nsubs
        figure('Color','w','Units','normalized','Position',[0.3 0.3 0.35 0.5]);
        plot(x,H2(:,subi),'Color',col_line,'LineWidth',1.5);
        xlabel(xlab,'FontSize',13,'FontWeight','bold');
        ylabel('H_2','FontSize',13,'FontWeight','bold');
        title(sprintf('Quadratic Rényi Entropy - Participant #%d',subi), ...
            'FontSize',14,'FontWeight','bold');
        set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
        grid on
        grid minor
        set_frequency_limits(x)
    end
end


%% Adjust scalar outputs

if numel(H2) == 1
    H2 = H2(1);
    ED = ED(1);
end

end


%% Helper function: Set frequency limits
function set_frequency_limits(frex)

if numel(frex) == 1
    padding = max(abs(frex(1))*0.05,0.5);
    xlim([frex(1)-padding frex(1)+padding])
else
    xlim([min(frex) max(frex)])
end

end


%% Helper function: Parse name-value pairs
function opts = parse_name_value_pairs(opts,varargin)

if mod(length(varargin),2) ~= 0
    error('Arguments must be given as name-value pairs.');
end

for i = 1:2:length(varargin)
    if ~(ischar(varargin{i}) || (isstring(varargin{i}) && isscalar(varargin{i})))
        error('Optional argument names must be text.');
    end
    name = lower(char(varargin{i}));
    if isfield(opts,name)
        opts.(name) = varargin{i+1};
    else
        error(['Unrecognized argument: ' name]);
    end
end

end
