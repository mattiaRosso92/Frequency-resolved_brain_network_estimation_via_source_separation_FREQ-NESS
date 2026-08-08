function [decayCoeff, goodFit] = FREQNESS_ExponentialDK(FREQ, varargin)

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
%  This function estimates how rapidly the eigenvalues associated with a
%  given component decay as a function of frequency. It does so by fitting
%  an exponentially decaying function to the eigenspectrum across
%  frequencies, either for all frequencies or within a user-defined
%  frequency range. 
%  
%  The decay rate can be interpreted as a compact index of
%  how strongly the component is dominated by low- vs. high-frequency
%  contributions, either across the entire spectrum or within a defined 
%  frequency range. 
% 
%  The coefficients in the output and can be compared across experimental 
%  conditions or groups of participants.
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
% - Optional arguments:
%
%   - which_comp : component index to analyse (default: 1)
% 
%   - range2fit  : frequency range to fit, expressed in Hz, e.g. [8 12].
%                  When empty, the exponential fit is performed across all
%                  available frequencies.
%
%   - plot_all   : logical flag (true/false). When true, the function
%                  produces additional figures for each individual
%                  participant, plotting their eigenspectrum and fitted
%                  exponential decay. By default, plot_all = false and
%                  only the grand-average fit is visualized.
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - decayCoeff : Vector [nSubs x 1] containing the exponential decay
%                 coefficient (lambda) estimated for the selected
%                 component in each participant. Larger positive values
%                 correspond to a steeper decay of eigenvalues across
%                 frequencies (i.e., stronger dominance of low frequencies).
%
%  - goodFit    : Structure containing goodness-of-fit information for the
%                 exponential model. Currently:
%                   • goodFit.R2  -> [nSubs x 1] vector with the coefficient
%                                    of determination (R²) for each subject.
%                 These values quantify how well the exponential model
%                 captures the eigenspectrum shape within the fitted
%                 frequency range.
%
% ------------------------------------------------------------------------
%  REFERENCES:
% ------------------------------------------------------------------------
%
%  - Kudryavtsev et al., in preparation
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

%% Handle optional arguments (name-value pairs)

opts = struct('which_comp', [], 'range2fit', [], 'plot_all',false);
opts = parse_name_value_pairs(opts, varargin{:});

which_comp = opts.which_comp;
range2fit  = opts.range2fit;
plot_all   = opts.plot_all;   

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

% Display for the user
fprintf('\nFREQNESS ExponentialDK: modelling eigenvalue decay from %.1f to %.1f Hz for %d participants.\n', ...
    range2fit(1), range2fit(end), nsubs);

%% Fit exponentially decaying function

% -------------------------------------------------------------------------
% 1) Select component
% -------------------------------------------------------------------------
if isempty(which_comp)
    disp('Component not specified. Defaulting to analyzing the 1st component.');
    which_comp = 1;
elseif ~isscalar(which_comp) || which_comp < 1 || which_comp > ncomps
    warning(['Input variable "which_comp" must be a scalar between 1 and ', num2str(ncomps), ...
             '. Defaulting to analyzing the 1st component.']);
    which_comp = 1;
end

% Extract eigenspectrum for the selected component: [nFrex x nSubs]
eig_comp = squeeze(eigenspectrum(which_comp,:,:)); 
if isvector(eig_comp)
    eig_comp = eig_comp(:); % [nFrex x 1] if single subject
end

% -------------------------------------------------------------------------
% 2) Define x-axis (frequency or index)
% -------------------------------------------------------------------------
if exist('frex','var') && numel(frex)==nfrex
    x_all = frex(:);
    xlab  = 'Frequency (Hz)';
else
    x_all = (1:nfrex).';
    xlab  = 'Frequency index';
end

% -------------------------------------------------------------------------
% 3) Map range2fit (in Hz) to indices to fit
% -------------------------------------------------------------------------
if isempty(range2fit)
    % Default: fit across all frequencies
    idx_range2fit = 1:nfrex;
else
    if ~exist('frex','var')
        error('range2fit was provided in Hz, but FREQ.frex is missing.');
    end
    if numel(range2fit) ~= 2
        error(['The variable range2fit must be a 2-element vector containing the boundaries ' ...
               'of the frequency range to fit the exponential decay, e.g. [8 12].']);
    end
    
    % Find closest 1st and last frequencies
    [~, idx_first] = min(abs(frex - range2fit(1)));
    [~, idx_last]  = min(abs(frex - range2fit(2)));
    idx_range2fit  = idx_first:idx_last;
    
    % Warn if requested boundaries are not exact matches
    tol = 1e-3; % Hz
    if any(~ismembertol(range2fit(:), frex(:), tol))
        warning('Some requested frequencies in range2fit are not present in FREQ.frex. Using closest matches instead.');
    end
end

% Frequencies (or indices) actually used for the fit
x_fit = x_all(idx_range2fit);

% -------------------------------------------------------------------------
% 4) Fit exponential decay per subject: y = A * exp(-lambda * x)
% -------------------------------------------------------------------------
decayCoeff = nan(nsubs,1);
goodFit.R2 = nan(nsubs,1);
A_sub      = nan(nsubs,1);   % store amplitudes for optional individual plots

for subi = 1:nsubs
    
    y = eig_comp(idx_range2fit, subi);
    
    % Keep only positive, finite values
    mask = y > 0 & isfinite(y);
    if sum(mask) < 3
        % Not enough data points to fit reliably
        continue
    end
    
    x_sub = x_fit(mask);
    y_sub = y(mask);
    
    % Linearize: log(y) = log(A) - lambda * x
    logy = log(y_sub);
    X    = [x_sub(:) ones(numel(x_sub),1)];
    
    beta = X \ logy;          % [slope; intercept]
    lambda = -beta(1);        % decay coefficient
    A      = exp(beta(2));    % amplitude (not returned, but used for plotting)
    
    % Store decay coefficient
    decayCoeff(subi) = lambda;
    A_sub(subi)      = A;
    
    % Compute R² in the original (non-log) space
    y_hat = A * exp(-lambda * x_sub);
    sse   = sum((y_sub - y_hat).^2);
    sst   = sum((y_sub - mean(y_sub)).^2);
    if sst > 0
        goodFit.R2(subi) = 1 - sse/sst;
    else
        goodFit.R2(subi) = NaN;
    end
end

% For visualization, also compute a group-level fit on the average eigenspectrum
mean_eig = mean(eig_comp, 2, 'omitnan');
y_group  = mean_eig(idx_range2fit);
mask_g   = y_group > 0 & isfinite(y_group);
if sum(mask_g) >= 3
    xg    = x_fit(mask_g);
    yg    = y_group(mask_g);
    logyg = log(yg);
    Xg    = [xg(:) ones(numel(xg),1)];
    betag = Xg \ logyg;
    lambda_group = -betag(1);
    A_group      = exp(betag(2));
else
    lambda_group = NaN;
    A_group      = NaN;
end

%% Visualize fit

% X-axis: use 'frex' if available, otherwise frequency index
if exist('frex','var') && numel(frex)==nfrex
    x = frex(:);
    xlab  = 'Frequency (Hz)';
else
    x = (1:nfrex).';
    xlab  = 'Frequency index';
end

% Mean ± SEM across subjects
if nsubs > 1
    avg_eig = mean(eig_comp, 2, 'omitnan');
    sem_eig = std(eig_comp, 0, 2, 'omitnan') ./ sqrt(nsubs);
else
    avg_eig = eig_comp;
    sem_eig = [];
end

% Parula colors
cmap    = .8*parula(256);
col_eig = cmap(1,:);      % eigenspectrum (purple)
col_fit = cmap(200,:);    % fitted curve (yellow)

figure('Color','w','Units','normalized','Position',[0.3 0.3 0.35 0.5]);
hold on

% SEM band
if nsubs > 1 && ~all(isnan(sem_eig))
    h_sem = fill([x; flipud(x)], [avg_eig - sem_eig; flipud(avg_eig + sem_eig)], ...
        [0 0 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

% Mean eigenspectrum
h_eig = plot(x, avg_eig, 'Color', col_eig, 'LineWidth', 2);

% Fitted exponential decay
if ~isnan(A_group) && ~isnan(lambda_group)
    x_fit_line = x(idx_range2fit);
    y_fit_line = A_group * exp(-lambda_group * x_fit_line);
    h_fit = plot(x_fit_line, y_fit_line, 'Color', col_fit, 'LineWidth', 2.4);
end

% Range markers
x_min_fit = x(idx_range2fit(1));
x_max_fit = x(idx_range2fit(end));
h_min = xline(x_min_fit, 'k--', 'LineWidth', 1);
h_max = xline(x_max_fit, 'k--', 'LineWidth', 1);

xlabel(xlab, 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Eigenvalue','FontSize',13,'FontWeight','bold');
title(sprintf('Exponential decay fit - Component #%d', which_comp), ...
      'FontSize',14,'FontWeight','bold');
set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
grid on; grid minor; xlim([min(x) max(x)]);
legend( ...
    [h_eig, h_fit, h_min], ...   % include only one dashed line (they are identical)
    {'Mean eigenspectrum', 'Exponential fit', 'Fit range'}, ...
    'Location','northeast', ...
    'FontSize',12);

% -------------------------------------------------------------------------
% Optional individual-subject plots
% -------------------------------------------------------------------------
if plot_all
    for subi = 1:nsubs
        
        y_full = eig_comp(:,subi);
        
        figure('Color','w','Units','normalized','Position',[0.3 0.3 0.35 0.5]);
        hold on
        
        % Subject eigenspectrum
        h_eig_s = plot(x, y_full, 'Color', col_eig, 'LineWidth', 1.5);
      
        % Fit
        if ~isnan(A_sub(subi)) && ~isnan(decayCoeff(subi))
            x_fit_line = x(idx_range2fit);
            y_fit_line = A_sub(subi) * exp(-decayCoeff(subi) * x_fit_line);
            h_fit_s = plot(x_fit_line, y_fit_line, 'Color', col_fit, 'LineWidth', 2.4);
        end
        
        % Range
        h_min_s = xline(x_min_fit, 'k--', 'LineWidth', 1);
        xline(x_max_fit, 'k--', 'LineWidth', 1);
        
        xlabel(xlab,'FontSize',13,'FontWeight','bold');
        ylabel('Eigenvalue','FontSize',13,'FontWeight','bold');
        title(sprintf('Subject %d – exponential decay fit (component %d)', subi, which_comp), ...
              'FontSize',14,'FontWeight','bold');
        set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
        grid on; grid minor; xlim([min(x) max(x)]);

        % Legend (per subject)
        if exist('h_fit_s','var')
            legend([h_eig_s, h_fit_s, h_min_s], ...
                {'Eigenspectrum','Exponential fit','Fit range'}, ...
                'Location','northeast','FontSize',12);
        else
            legend([h_eig_s, h_min_s], ...
                {'Eigenspectrum','Fit range'}, ...
                'Location','northeast','FontSize',12);
        end

    end
end

end


%% Helper Function: Parse Name-Value Pairs
function opts = parse_name_value_pairs(opts, varargin)
% Simple name-value parser used to handle optional arguments.
% Example:
%   opts = struct('which_comp', [], 'range2fit', [], 'plot_all', false);
%   opts = parse_name_value_pairs(opts, 'which_comp', 2, 'range2fit', [8 12], 'plot_all', true);

if mod(length(varargin), 2) ~= 0
    error('Arguments must be given as name-value pairs.');
end

for i = 1:2:length(varargin)
    name = lower(varargin{i});
    if isfield(opts, name)
        opts.(name) = varargin{i+1};
    else
        error(['Unrecognized argument: ', name]);
    end
end
end
