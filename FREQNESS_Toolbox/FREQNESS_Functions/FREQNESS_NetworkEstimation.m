function [FREQ] = FREQNESS_NetworkEstimation(data, frex, srate, varargin)

% ========================================================================

%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  Please cite the first FREQNESS paper:
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


%% Handle Optional Arguments

% Define default values
opts = struct('duration', [], 'fwidth', [], 'filter', 'logarithmic', 'regularisation', 0.01, 'ncomps', 30, 'bad_segments',[]);

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

% Compute necessary parameters
pnts2keep = duration * srate;

% Check duration
if pnts2keep > size(data,2)
    error('Requested duration exceeds the length of the data.');
end

% Check that the range of the bad segments falls within the requested duration
if any(idx2remove < 1) || any(idx2remove > pnts2keep)
    error('One or more indices in bad_segments fall outside the analyzed time range.');
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
    elseif strcmpi(filter, 'linear')
        fwidth_all = linspace(fwidth, fwidth * nfrex, nfrex);
    else
        disp("Invalid filter type. Defaulting to 'logarithmic'.");
        fwidth_all = logspace(log10(fwidth), log10(fwidth * nfrex), nfrex);
    end

else

    % Assign input filter width vector, if given by the user
    fwidth_all = fwidth;

end



%% Input data check

% Check input size to verify if more than one participant is being analyzed
if length(size(data)) < 3 % no extra dimension in the input
    nsubs = 1;
else
    nsubs = size(data,3);
end

% Occasionally, anatomical source reconstruction can return extremely low
% values. This can cause unreliable source separation due to the
% scale of the regularization factor of the covariance matrices.

fprintf('\nFREQNESS Network Estimation: estimating frequency-resolved networks for %d participants.\n', nsubs);
for subi = 1:nsubs

    % Extract data for this participants
    this_data = data(:,:,subi);

    % Robust estimate of the scale of the data
    scale_ref = median(abs(this_data(:)));

    % Detect current order of magnitude
    order_mag = floor(log10(scale_ref));

    if order_mag < -2
        warning(['The scale of your data is very low (order of magnitude: 10e' num2str(order_mag) ').' ...
            ' This can cause unreliable network estimation.']);

        flag_scale = input('Do you want to re-scale the data to bring values into the hundreds range? (1 = yes; 0 = no):');

        if flag_scale == 1
            % Define desired target order
            scale = 10^abs(order_mag - 1);  % re-scale to bring values into the hundreds range

            % Assign re-scaled data
            data(:,:,subi) = scale * this_data;
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

    % Define good indexes to keep for computing the covariance matrices
    idx2cov = setdiff(1:pnts2keep,idx2remove);

    % Remove bad segments, if any bad segments is provided as input
    if ~isempty(idx2remove)
        pct_removed = 100 * length(idx2remove) / pnts2keep;
        disp(['Removing ' num2str(length(idx2remove)) ' bad timepoints ' ...
            '(' num2str(pct_removed) '% of the input data).']);
    end


    % Perform frequency-resolved brain network separation via GED
    for frexi = 1:nfrex % loop over input frequencies
        disp(['Estimating network at ' num2str(frex(frexi)) ' Hz'])

        % Filter data
        narrowData = filterFGx(broadData,srate,frex(frexi),fwidth_all(frexi),0); % turn the last argument to 1 to visualize the filter

        % Compute covariance matrices
        covS = cov(narrowData(:,idx2cov)');
        covR = cov(broadData(:,idx2cov)');

        % Regularisation (to bring covariance matrices to full rank)
        covS = covS  + 1e-6*eye(size(covS)); % regularise the S covariance matrix by adding a small perturbation/noise
        evalsR = eig(covR ); % re-compute eigenvalues of R covariance matrix for its regularisation
        covR  = (1-regularisation)*covR  + regularisation * mean(evalsR) * eye(size(covR)); % regularise the R covariance matrix

        % Eigendecomposition
        [evecs,evals] = eig(covS ,covR);
        [evals,sidx]  = sort( diag(evals),'descend' ); % the first output returns sorted evals extracted from diagonal
        evecs = evecs(:,sidx);          % sort eigenvectors
        evals = evals.*100./sum(evals); % normalize eigenvalues to percent variance explained
        % Assign temporary variables to output
        GEDevecs(:,:,frexi,subi) = evecs(:,1:ncomps);
        GEDevals(:,frexi,subi) = evals(1:ncomps);

        % Brain networks' spatial activation patterns and time series
        for compi = 1:ncomps

            % Compute spatial activation patterns and flip sign
            GEDpats(:,compi,frexi,subi) = GEDevecs(:,compi,frexi,subi)' * covS; % get component
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


