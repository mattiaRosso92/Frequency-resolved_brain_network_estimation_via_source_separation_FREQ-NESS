function [FREQ] = FREQNESS_NetworkEstimation(data, frex, srate, varargin)

% ========================================================================

%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you use this toolbox, please cite:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P., Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%
%  This script takes a multivariate dataset (channels × time matrix)
%  and produces a frequency-resolved brain network landscape via
%  Generalized Eigendecomposition (GED).
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - data   : 2D matrix in channels-by-time format, or
%             3D matrix in channels-by-time-by-subjects format
%             NOTE: While this approach can be applied to scalp sensor data,
%             source-reconstructed brain voxel data is necessary to
%             interpret the output as brain networks.
%
%  - frex   : Vector of frequencies (ascending order).
%
%  - srate  : Sampling rate (Hz).
%
%  The following parameters can be specified using name-value pairs:
%
%  - 'duration'        : Duration of the recording to analyze (seconds).
%                        Default: full data duration.
%
%  - 'fwidth'          : Full-width at half-maximum (FWHM) of the frequency filter.
%                        Default: computed based on reference values.
%
%  - 'filter'          : Filter design method. Options:
%                        'logarithmic' (default) or 'linear'.
%
%  - 'regularisation'  : Shrinkage factor for covariance matrix regularization.
%                        Default: 0.01.
%
%  - 'ncomps'          : Number of components to retain per frequency.
%                        Default: 30.
%
%  - 'bad_segments'    : Vector of sample indices (e.g., 1-based timepoints)
%                        to be excluded from the analysis.
%                        These will be removed from covariance computations.
%
%  - 'rescale'         : Logical flag controlling automatic rescaling of
%                        very low-amplitude data. Default: false.
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - FREQ       : Structure containing eigenvalues, eigenvectors, spatial
%                 activation patterns, network time series, and additional 
%                 variables which will be used by second-order functions in
%                 the toolbox.
%                 Info on the specific fields is found in the comments where
%                 the respective variables are assigned to GED.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso & Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 12/06/2025

% ========================================================================


% NOTE: credit for the GED implementation goes to Mike X Cohen
% (Cohen, 2022 - 'A tutorial on generalized eigendecomposition for denoising,
% contrast enhancement, and dimension reduction in multichannel electrophysiology')


%% Check mandatory inputs

if nargin < 3
    error('data, frex, and srate are required inputs.');
end

if ~isnumeric(data) || ~isreal(data) || isempty(data) || ...
        ndims(data) > 3 || any(~isfinite(data(:)))
    error(['data must be a non-empty, real, finite numeric array in ' ...
        '[channels x time] or [channels x time x participants] format.']);
end

if ~isnumeric(srate) || ~isreal(srate) || ~isscalar(srate) || ...
        ~isfinite(srate) || srate <= 0
    error('srate must be one positive finite scalar.');
end

if ~isnumeric(frex) || ~isreal(frex) || ~isvector(frex) || isempty(frex) || ...
        any(~isfinite(frex)) || any(frex <= 0)
    error('frex must be a non-empty vector of positive finite frequencies.');
end

if any(frex >= srate/2)
    error('All frequencies in frex must be lower than the Nyquist frequency.');
end

%% Handle Optional Arguments

% Define default values
opts = struct('duration', [], 'fwidth', [], 'filter', 'logarithmic', ...
    'regularisation', 0.01, 'ncomps', 30, 'bad_segments',[], ...
    'rescale',false);

% Parse name-value pair arguments
opts = parse_name_value_pairs(opts, varargin{:});

% If duration is empty, use full data duration
if isempty(opts.duration)
    opts.duration = size(data, 2) / srate;
end

% Assign values
duration       = opts.duration;
fwidth         = opts.fwidth;
filter         = opts.filter;
regularisation = opts.regularisation;
ncomps         = opts.ncomps;
idx2remove     = opts.bad_segments; % in the current version, these apply to all subjects in the input
rescale        = opts.rescale;

% Validate optional inputs
if ~isempty(duration) && (~isnumeric(duration) || ~isreal(duration) || ...
        ~isscalar(duration) || ~isfinite(duration) || duration <= 0)
    error('duration must be one positive finite scalar expressed in seconds.');
end

if ~(ischar(filter) || (isstring(filter) && isscalar(filter))) || ...
        ~(strcmpi(filter,'logarithmic') || strcmpi(filter,'linear'))
    error('filter must be either ''logarithmic'' or ''linear''.');
end

if ~isnumeric(regularisation) || ~isreal(regularisation) || ...
        ~isscalar(regularisation) || ~isfinite(regularisation) || ...
        regularisation < 0 || regularisation >= 1
    error('regularisation must be one finite scalar in the interval [0, 1).');
end

if ~isnumeric(ncomps) || ~isreal(ncomps) || ~isscalar(ncomps) || ...
        ~isfinite(ncomps) || ncomps ~= round(ncomps) || ncomps < 1
    error('ncomps must be one positive integer.');
end


if ncomps > size(data,1)
    error('ncomps cannot exceed the number of input channels or voxels.');
end

if ~islogical(rescale) || ~isscalar(rescale)
    error('rescale must be one logical value (true or false).');
end

% Compute necessary parameters
pnts2keep = duration * srate;

if abs(pnts2keep - round(pnts2keep)) > eps(max(1,pnts2keep))
    error('duration multiplied by srate must define an integer number of samples.');
end
pnts2keep = round(pnts2keep);

% Check duration
if pnts2keep > size(data,2)
    error('Requested duration exceeds the length of the data.');
end

if pnts2keep < 2
    error('The requested duration must contain at least two samples.');
end

if ~isempty(idx2remove) && (~isnumeric(idx2remove) || ~isreal(idx2remove) || ...
        ~isvector(idx2remove) || any(~isfinite(idx2remove)) || ...
        any(idx2remove ~= round(idx2remove)))
    error('bad_segments must contain finite integer sample indices.');
end

idx2remove = unique(idx2remove(:)','stable');

% Check that the range of the bad segments falls within the requested duration
if any(idx2remove < 1) || any(idx2remove > pnts2keep)
    error('One or more indices in bad_segments fall outside the analyzed time range.');
end

% Define good indexes to keep for computing the covariance matrices
idx2cov = setdiff(1:pnts2keep,idx2remove);
if numel(idx2cov) < 2
    error('Fewer than two samples remain for covariance estimation.');
end

% Re-sort frequency vector in ascending order
frex = sort(frex,'ascend');
nfrex = length(frex);

% Handle missing filter width (based on Rosso et al., 2025 - Advanced Science)
if isempty(fwidth)

    ref_frex = 2.439; % set a reference frequency
    ref_fwhm = 0.35;  % set a reference filter width
    nfrex_above = 80; % how many frequencies above the reference
    nfrex_below = 6;  % how many frequencies below the reference
    % Compute the vector of frequencies
    frex_above = linspace(ref_frex, ref_frex * (nfrex_above / 2), nfrex_above);
    frex_below = linspace(ref_frex / 2, ref_frex / (2 * nfrex_below), nfrex_below);
    frex_all = [frex_below(end:-1:1), frex_above];
    % Compute the vector of filter width (log-spaced)
    fwhm_above = logspace(log10(ref_fwhm), log10(ref_fwhm * (nfrex_above)), nfrex_above);
    fwhm_below = logspace(log10(ref_fwhm), log10(ref_fwhm / (nfrex_below)), nfrex_below + 1);
    fwidth_all = [fwhm_below(end:-1:2), fwhm_above];
    % Match user-requested frequencies with the reference frequencies from Rosso et al. (2025)
    [~, idx] = min(abs(frex_all - frex(1))); % matching the lowest user-defined frequency
    fwidth = fwidth_all(idx);

    % Compute filter widths
    if strcmpi(filter, 'logarithmic')
        fwidth_all = logspace(log10(fwidth), log10(fwidth * nfrex), nfrex);
    else
        fwidth_all = linspace(fwidth, fwidth * nfrex, nfrex);
    end

else

    % Assign input filter width vector, if given by the user
    if isscalar(fwidth)
        fwidth_all = repmat(fwidth,1,nfrex);
    elseif isvector(fwidth) && numel(fwidth)==nfrex
        fwidth_all = fwidth(:)';
    else
        error('fwidth must be a scalar or contain one value per frequency.');
    end

end

if ~isnumeric(fwidth_all) || ~isreal(fwidth_all) || ...
        any(~isfinite(fwidth_all)) || any(fwidth_all <= 0)
    error('fwidth must contain positive finite values.');
end


%% Input data check

% Check input size to verify if more than one participant is being analyzed
nsubs = size(data,3);

% Occasionally, anatomical source reconstruction can return extremely low
% values. This can cause unreliable source separation due to the
% scale of the regularization factor of the covariance matrices.

fprintf('\nFREQNESS Network Estimation: estimating frequency-resolved networks for %d participants.\n', nsubs);
scale_factors = ones(1,nsubs);
for subi = 1:nsubs

    % Extract data for this participant
    this_data = data(:,1:pnts2keep,subi);

    % Robust estimate of the scale of the data
    scale_ref = median(abs(this_data(:)));

    if scale_ref == 0
        error(['The median absolute data value is zero for participant #' ...
            num2str(subi) '. Data scale cannot be estimated reliably.']);
    end

    % Detect current order of magnitude
    order_mag = floor(log10(scale_ref));

    if order_mag < -2
        warning(['The scale of your data is very low (order of magnitude: 10e' num2str(order_mag) ').' ...
            ' This can cause unreliable network estimation.']);

        if rescale
            % Define desired target order
            scale = 10^abs(order_mag - 1);  % re-scale to bring values into the hundreds range

            % Assign re-scaled data
            data(:,1:pnts2keep,subi) = scale * this_data;
            scale_factors(subi) = scale;
        end
    end

end

%% Generalized eigendecomposition (GED) 

% Initialize output structure
FREQ = [];

% Initialize GED outputs
GEDevals = zeros(ncomps, nfrex, nsubs);
[GEDevecs, GEDpats] = deal(zeros(size(data,1), ncomps, nfrex, nsubs));
GEDts = zeros(ncomps, pnts2keep, nfrex, nsubs);

% Loop over participants (also works with one single participant)
for subi = 1:nsubs
    disp(['FREQNESS Network Estimation: participant #' num2str(subi)]);

    % Assign broadband source data
    broadData = data(:, 1:pnts2keep, subi);

    % Remove bad segments, if any bad segments is provided as input
    if ~isempty(idx2remove)
        pct_removed = 100 * length(idx2remove) / pnts2keep;
        disp(['Removing ' num2str(length(idx2remove)) ' bad timepoints ' ...
            '(' num2str(pct_removed) '% of the input data).']);
    end


    % Compute and regularise the broadband covariance matrix once
    covR = cov(broadData(:,idx2cov)');
    evalsR = eig(covR);
    covR = (1-regularisation)*covR + ...
        regularisation * mean(evalsR) * eye(size(covR));

    % Perform frequency-resolved brain network separation via GED
    for frexi = 1:nfrex % loop over input frequencies
        disp(['Estimating network at ' num2str(frex(frexi)) ' Hz'])

        % Filter data
        narrowData = filterFGx(broadData,srate,frex(frexi),fwidth_all(frexi),0); % turn the last argument to 1 to visualize the filter

        % Compute narrowband covariance matrix
        covS = cov(narrowData(:,idx2cov)');

        % Regularisation (to bring the covariance matrix to full rank)
        covS = covS  + 1e-6*eye(size(covS)); % regularise the S covariance matrix by adding a small perturbation/noise

        % Eigendecomposition
        [evecs,evals] = eig(covS ,covR);
        [evals,sidx]  = sort( diag(evals),'descend' ); % the first output returns sorted evals extracted from diagonal
        evecs = evecs(:,sidx);          % sort eigenvectors
        if ~isfinite(sum(evals)) || abs(sum(evals)) < realmin
            error('Generalized eigenvalues cannot be normalized because their sum is zero or non-finite.');
        end
        evals = evals.*100./sum(evals); % normalize eigenvalues to percent variance explained
        % Assign temporary variables to output
        GEDevecs(:,:,frexi,subi) = evecs(:,1:ncomps);
        GEDevals(:,frexi,subi) = evals(1:ncomps);

        % Brain networks' spatial activation patterns and time series
        for compi = 1:ncomps

            % Compute spatial activation patterns and flip sign
            GEDpats(:,compi,frexi,subi) = covS * GEDevecs(:,compi,frexi,subi); % get component
            [~,idxmax] = max(abs(GEDpats(:,compi,frexi,subi)));     % find max magnitude
            GEDpats(:,compi,frexi,subi)  = GEDpats(:,compi,frexi,subi) * sign(GEDpats(idxmax,compi,frexi,subi)); % possible sign flip

            % Compute the network's activation time series
            GEDts(compi,:,frexi,subi) = GEDevecs(:,compi,frexi,subi)' * broadData;

        end

    end

end

% Assign outputs to structure
FREQ.evals = GEDevals;    % eigenvalues
FREQ.evecs = GEDevecs;    % eigenvectors
FREQ.pats  = GEDpats;     % spatial activation patterns
FREQ.ts    = GEDts;       % time series
FREQ.frex  = frex;        % to make use of the analyzed frequencies later on
FREQ.fwhm  = fwidth_all;  % to make use of the corresponding filter widths
FREQ.srate = srate;       % re-assign to output, so it can be used by secondary functions
FREQ.duration = pnts2keep/srate;        % analyzed duration in seconds
FREQ.bad_segments = idx2remove;         % excluded sample indices
FREQ.regularisation = regularisation;   % covariance shrinkage factor
FREQ.scale_factors = scale_factors;     % participant-wise internal scale factors

end


%% Helper Function: Parse Name-Value Pairs
function opts = parse_name_value_pairs(opts, varargin)
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

%%
