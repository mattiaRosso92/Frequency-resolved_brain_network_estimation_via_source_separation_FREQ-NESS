function [HYPER] = HYPERNESS_NetworkEstimation(data, frex, srate, varargin)

% ========================================================================

%  HYPER-BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION TOOLBOX
%
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  PUBLICATION PLACEHOLDER - REPLACE BEFORE RELEASE
%
%  [AUTHORS]. ([YEAR]).
%  [HYPER-NESS PAPER TITLE].
%  [JOURNAL, VOLUME, PAGES].
%  [DOI].
%
%  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%
% ========================================================================
%
%  This function takes a multivariate hyperscanning dataset in
%  channels-by-time-by-participants format and produces a
%  frequency-resolved hyper-brain network landscape via Generalized
%  Eigendecomposition (GED).
%
%  Participant data are internally concatenated along the channel
%  dimension:
%
%      [channels x time x participants]
%
%                  becomes
%
%      [(channels*participants) x time]
%
%  GED is then performed once on the complete joint covariance matrices.
%  This preserves both the diagonal intra-participant covariance blocks
%  and the off-diagonal inter-participant cross-covariance blocks during
%  network estimation.
%
%  For every retained component, the joint eigenvector is partitioned back
%  into participant-specific subvectors. The narrowband component variance
%  is decomposed exactly as:
%
%      Total = sum(Intra_p) + sum(Inter_pq)
%
%  where:
%
%      Intra_p  = w_p' * S_pp * w_p
%
%      Inter_pq = 2 * w_p' * S_pq * w_q, for p < q
%
%  The factor of two in Inter_pq accounts for both symmetric off-diagonal
%  covariance blocks S_pq and S_qp.
%
%  NOTE: inter-participant contributions are signed. Their ratios can
%  therefore be negative or greater than one and must not be interpreted
%  as non-negative mixture proportions.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - data   : 3D numeric matrix in
%             channels-by-time-by-participants format.
%
%             All participants must have:
%               1) the same number of channels;
%               2) the same number of time samples;
%               3) sample-by-sample temporal synchronization;
%               4) data expressed in the same physical units.
%
%             The number of participants is inferred from size(data,3).
%             At least two participants are required.
%
%             NOTE: While this approach can be applied to scalp sensor
%             data, source-reconstructed brain voxel data are necessary to
%             interpret the spatial output as source-level brain networks.
%
%  - frex   : Vector of frequencies in Hz (ascending order).
%
%  - srate  : Sampling rate in Hz, common to all participants.
%
%  The following parameters can be specified using name-value pairs:
%
%  - 'duration'        : Duration of the recording to analyze (seconds).
%                        Default: full data duration.
%
%  - 'fwidth'          : Full-width at half-maximum (FWHM) of the
%                        frequency filter. A scalar applies the same width
%                        to every frequency. A vector must contain one
%                        value per frequency.
%                        Default: computed based on FREQ-NESS reference
%                        values.
%
%  - 'filter'          : Filter-width progression. Options:
%                        'logarithmic' (default) or 'linear'.
%                        This setting is used only when 'fwidth' is empty.
%
%  - 'regularisation'  : Shrinkage factor for covariance matrix R.
%                        Default: 0.01.
%
%  - 'ncomps'          : Number of components to retain per frequency.
%                        Default: 30.
%
%  - 'bad_segments'    : Vector of synchronized sample indices to exclude
%                        from covariance computations.
%
%                        The same indices are removed for every participant
%                        because inter-participant covariance requires
%                        aligned observations.
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENT:
% ------------------------------------------------------------------------
%
%  - HYPER : Structure containing eigenvalues, eigenvectors, spatial
%            activation patterns, network time series, participant-specific
%            projections, and the static inter/intra decomposition.
%
%    Fields corresponding directly to FREQNESS_NetworkEstimation:
%
%      HYPER.evals : Eigenvalues
%                    [nComponents x nFrequencies].
%
%      HYPER.evecs : Joint eigenvectors
%                    [nJointChannels x nComponents x nFrequencies].
%
%      HYPER.pats  : Joint spatial activation patterns
%                    [nJointChannels x nComponents x nFrequencies].
%
%      HYPER.ts    : Joint broadband network time series
%                    [nComponents x nTime x nFrequencies].
%
%      HYPER.frex  : Frequencies used in the analysis.
%
%      HYPER.fwhm  : Filter widths corresponding to HYPER.frex.
%
%      HYPER.srate : Sampling rate.
%
%    HYPER-NESS-specific fields:
%
%      HYPER.participant_evecs
%          Participant partitions of the joint eigenvectors
%          [nChannels x nComponents x nFrequencies x nParticipants].
%
%      HYPER.participant_pats
%          Participant partitions of the joint activation patterns
%          [nChannels x nComponents x nFrequencies x nParticipants].
%
%      HYPER.participant_ts
%          Participant-specific broadband component projections
%          [nComponents x nTime x nFrequencies x nParticipants].
%
%      HYPER.total
%          Total narrowband component variance
%          [nComponents x nFrequencies].
%
%      HYPER.intra
%          Participant-wise intra-brain contributions
%          [nComponents x nFrequencies x nParticipants].
%
%      HYPER.inter
%          Pairwise inter-brain contributions
%          [nComponents x nFrequencies x nPairs].
%
%      HYPER.intra_total
%      HYPER.inter_total
%          Contributions summed across participants or participant pairs
%          [nComponents x nFrequencies].
%
%      HYPER.intra_ratio
%      HYPER.inter_ratio
%          Individual contributions divided by HYPER.total.
%
%      HYPER.intra_total_ratio
%      HYPER.inter_total_ratio
%          Aggregate contributions divided by HYPER.total.
%
%      HYPER.pairs
%          Participant indices corresponding to the third dimension of
%          HYPER.inter [nPairs x 2].
%
%      HYPER.decomposition_error
%          Relative numerical error between HYPER.total and the sum of all
%          intra- and inter-participant contributions.
%
%      HYPER.nparticipants
%      HYPER.nchannels
%      HYPER.channel_indices
%          Metadata required to map the joint solution back to individual
%          participants.
%
% ------------------------------------------------------------------------
%  SIGN CONVENTION:
% ------------------------------------------------------------------------
%
%  GED eigenvectors are defined only up to a global sign. No sign
%  correction is applied to HYPER.evecs, HYPER.ts, participant-specific
%  filters, participant-specific time series, or variance contributions.
%
%  Following FREQNESS_NetworkEstimation, only each spatial activation
%  pattern is oriented so that its largest-magnitude coefficient is
%  positive. This convention affects visualization only.
%
%  Participant subvectors must never be sign-flipped independently:
%  independent flips would change the pairwise inter-brain contributions
%  and would no longer represent the original joint GED solution.
%
% ------------------------------------------------------------------------
%  AUTHORS:
% ------------------------------------------------------------------------
%
%  [TO BE COMPLETED]
%
% ========================================================================


% NOTE: credit for the GED implementation goes to Mike X Cohen
% (Cohen, 2022 - 'A tutorial on generalized eigendecomposition for
% denoising, contrast enhancement, and dimension reduction in multichannel
% electrophysiology').


%% Handle Optional Arguments

% Define default values
opts = struct('duration', [], 'fwidth', [], ...
    'filter', 'logarithmic', 'regularisation', 0.01, ...
    'ncomps', 30, 'bad_segments', []);

% Parse name-value pair arguments
opts = parse_name_value_pairs(opts, varargin{:});

% Assign values
duration       = opts.duration;
fwidth         = opts.fwidth;
filter         = opts.filter;
regularisation = opts.regularisation;
ncomps         = opts.ncomps;
idx2remove     = opts.bad_segments;


%% Input data check

if ~isnumeric(data) || isempty(data) || ~isreal(data)
    error('Data must be a non-empty real numeric matrix.');
end

if ndims(data) > 3
    error('Data must be provided as channels-by-time-by-participants.');
end

[nchannels, npnts, nparticipants] = size(data);

if nparticipants < 2
    error(['HYPER-NESS requires at least two simultaneously recorded ' ...
        'participants in the third dimension of the input data.']);
end

if any(~isfinite(data(:)))
    error('Data contain NaN or Inf values.');
end

if ~isnumeric(frex) || isempty(frex) || ~isvector(frex) || ...
        any(~isfinite(frex)) || any(frex <= 0)
    error('Frequencies must be provided as a vector of positive values.');
end

if ~isnumeric(srate) || ~isscalar(srate) || ...
        ~isfinite(srate) || srate <= 0
    error('Sampling rate must be a positive scalar.');
end

if any(frex >= srate/2)
    error('All frequencies must be lower than the Nyquist frequency.');
end

if ~isnumeric(regularisation) || ~isscalar(regularisation) || ...
        ~isfinite(regularisation) || ...
        regularisation < 0 || regularisation >= 1
    error('Regularisation must be a scalar in the interval [0,1).');
end

if ~isnumeric(ncomps) || ~isscalar(ncomps) || ...
        ~isfinite(ncomps) || ncomps < 1 || ncomps ~= round(ncomps)
    error('ncomps must be a positive integer.');
end


%% Define data duration

% If duration is empty, use full data duration
if isempty(duration)
    duration = npnts/srate;
end

% Compute necessary parameters
pnts2keep = round(duration*srate);

% Check duration
if pnts2keep > npnts
    error('Requested duration exceeds the length of the data.');
end

if pnts2keep < 2
    error('The requested duration contains fewer than two samples.');
end

% Check bad segments
if ~isempty(idx2remove)
    if ~isnumeric(idx2remove) || ~isvector(idx2remove) || ...
            any(~isfinite(idx2remove)) || ...
            any(idx2remove ~= round(idx2remove))
        error('bad_segments must contain integer sample indices.');
    end

    idx2remove = unique(idx2remove(:)');

    if any(idx2remove < 1) || any(idx2remove > pnts2keep)
        error(['One or more indices in bad_segments fall outside the ' ...
            'analyzed time range.']);
    end
end

% Define synchronized indexes used in all covariance blocks
idx2cov = setdiff(1:pnts2keep, idx2remove);

if length(idx2cov) < 2
    error('Fewer than two samples remain for covariance estimation.');
end


%% Define frequencies and filter widths

% Re-sort frequency vector in ascending order
frex = sort(frex, 'ascend');
frex = frex(:)';
nfrex = length(frex);

% Handle missing filter width (based on FREQ-NESS reference values)
if isempty(fwidth)

    ref_frex = 2.439; % set a reference frequency
    ref_fwhm = 0.35;  % set a reference filter width
    nfrex_above = 80; % how many frequencies above the reference
    nfrex_below = 6;  % how many frequencies below the reference

    % Compute the reference vector of frequencies
    frex_above = linspace( ...
        ref_frex, ref_frex*(nfrex_above/2), nfrex_above);
    frex_below = linspace( ...
        ref_frex/2, ref_frex/(2*nfrex_below), nfrex_below);
    frex_all = [frex_below(end:-1:1), frex_above];

    % Compute the reference vector of filter widths
    fwhm_above = logspace( ...
        log10(ref_fwhm), log10(ref_fwhm*nfrex_above), nfrex_above);
    fwhm_below = logspace( ...
        log10(ref_fwhm), log10(ref_fwhm/nfrex_below), ...
        nfrex_below+1);
    fwidth_ref = [fwhm_below(end:-1:2), fwhm_above];

    % Match the lowest user-defined frequency
    [~, idx] = min(abs(frex_all-frex(1)));
    fwidth_lowest = fwidth_ref(idx);

    % Compute filter widths
    if strcmpi(filter, 'logarithmic')
        fwidth_all = logspace( ...
            log10(fwidth_lowest), ...
            log10(fwidth_lowest*nfrex), nfrex);
    elseif strcmpi(filter, 'linear')
        fwidth_all = linspace( ...
            fwidth_lowest, fwidth_lowest*nfrex, nfrex);
    else
        disp("Invalid filter type. Defaulting to 'logarithmic'.");
        fwidth_all = logspace( ...
            log10(fwidth_lowest), ...
            log10(fwidth_lowest*nfrex), nfrex);
    end

else

    if ~isnumeric(fwidth) || any(~isfinite(fwidth)) || any(fwidth <= 0)
        error('fwidth must contain positive finite values.');
    end

    if isscalar(fwidth)
        fwidth_all = repmat(fwidth, 1, nfrex);
    elseif isvector(fwidth) && length(fwidth) == nfrex
        fwidth_all = fwidth(:)';
    else
        error(['fwidth must be a scalar or a vector containing one value ' ...
            'per frequency.']);
    end

end


%% Joint data construction

njoint_channels = nchannels*nparticipants;

% Define the joint-channel ranges associated with each participant
channel_indices = zeros(nparticipants, 2);

for participant_i = 1:nparticipants
    channel_indices(participant_i,1) = ...
        (participant_i-1)*nchannels + 1;
    channel_indices(participant_i,2) = ...
        participant_i*nchannels;
end

% Define the unique participant pairs
pairs = nchoosek(1:nparticipants, 2);
npairs = size(pairs,1);

% Retain the requested duration
data = data(:,1:pnts2keep,:);

% IMPORTANT:
% A direct reshape of channels x time x participants would incorrectly
% interleave time and participant samples. The participant dimension is
% first moved next to the channel dimension.
jointData = reshape( ...
    permute(data, [1 3 2]), ...
    njoint_channels, pnts2keep);

% Check that the stacking order is exactly:
% Participant 1 channels, Participant 2 channels, ..., Participant P.
for participant_i = 1:nparticipants
    this_idx = channel_indices(participant_i,1): ...
        channel_indices(participant_i,2);

    if ~isequal(jointData(this_idx,:), data(:,:,participant_i))
        error('Internal participant stacking failed.');
    end
end


%% Joint data scale check

% A common scaling factor can be applied to the entire hyperscanning
% system without changing the relative weighting between participants.
% Participants must never be rescaled independently here, because doing so
% would alter the cross-participant covariance blocks.

scale_ref = median(abs(jointData(:)));

if scale_ref == 0
    error('The median absolute data value is zero.');
end

order_mag = floor(log10(scale_ref));

if order_mag < -2
    warning(['The scale of the joint data is very low (order of magnitude: ' ...
        '10e' num2str(order_mag) '). This can cause unreliable network ' ...
        'estimation because of covariance regularization.']);

    flag_scale = input(['Do you want to apply one common scaling factor ' ...
        'to all participants? (1 = yes; 0 = no): ']);

    if flag_scale == 1
        scale = 10^abs(order_mag-1);
        jointData = scale*jointData;
        data = scale*data;
    end
end


%% Define number of retained components

if ncomps > njoint_channels
    warning(['Requested ncomps exceeds the number of joint channels. ' ...
        'ncomps was reduced to %d.'], njoint_channels);
    ncomps = njoint_channels;
end


%% Initialize GED outputs

HYPER = [];

GEDevals = zeros(ncomps, nfrex);
GEDevecs = zeros(njoint_channels, ncomps, nfrex);
GEDpats  = zeros(njoint_channels, ncomps, nfrex);
GEDts    = zeros(ncomps, pnts2keep, nfrex);

participant_evecs = zeros( ...
    nchannels, ncomps, nfrex, nparticipants);
participant_pats = zeros( ...
    nchannels, ncomps, nfrex, nparticipants);
participant_ts = zeros( ...
    ncomps, pnts2keep, nfrex, nparticipants);

total_contribution = zeros(ncomps, nfrex);
intra_contribution = zeros(ncomps, nfrex, nparticipants);
inter_contribution = zeros(ncomps, nfrex, npairs);
decomposition_error = zeros(ncomps, nfrex);


%% Generalized eigendecomposition (GED)

fprintf(['\nHYPERNESS Network Estimation: estimating frequency-resolved ' ...
    'joint networks for %d participants.\n'], nparticipants);

fprintf(['Joint data dimensions: %d channels ' ...
    '(%d channels x %d participants) x %d timepoints.\n'], ...
    njoint_channels, nchannels, nparticipants, pnts2keep);

broadData = jointData;

% Remove bad segments only from covariance computations
if ~isempty(idx2remove)
    pct_removed = 100*length(idx2remove)/pnts2keep;
    disp(['Removing ' num2str(length(idx2remove)) ...
        ' synchronized bad timepoints (' ...
        num2str(pct_removed) '% of the input data).']);
end

% Perform frequency-resolved hyper-brain network separation via GED
for frexi = 1:nfrex

    disp(['Estimating joint network at ' num2str(frex(frexi)) ' Hz']);

    % Filter complete joint data
    narrowData = filterFGx( ...
        broadData, srate, frex(frexi), fwidth_all(frexi), 0);

    % Compute complete joint covariance matrices
    covS = cov(narrowData(:,idx2cov)');
    covR = cov(broadData(:,idx2cov)');

    % Regularisation, following FREQNESS_NetworkEstimation
    covS = covS + 1e-6*eye(size(covS));
    evalsR = eig(covR);
    covR = (1-regularisation)*covR + ...
        regularisation*mean(evalsR)*eye(size(covR));

    % Eigendecomposition
    [evecs, evals] = eig(covS, covR);
    [evals, sidx] = sort(diag(evals), 'descend');
    evecs = evecs(:,sidx);
    evals = evals.*100./sum(evals);

    % Assign eigenvalues and joint eigenvectors
    GEDevals(:,frexi) = evals(1:ncomps);
    GEDevecs(:,:,frexi) = evecs(:,1:ncomps);

    % Spatial patterns, time series, and blockwise decomposition
    for compi = 1:ncomps

        this_evec = GEDevecs(:,compi,frexi);

        % Compute the joint spatial activation pattern
        this_pattern = covS*this_evec;

        % FREQ-NESS pattern-only sign convention:
        % orient the largest-magnitude pattern coefficient positively.
        % The eigenvector and time series are not sign-flipped.
        [~,idxmax] = max(abs(this_pattern));
        this_pattern = this_pattern*sign(this_pattern(idxmax));

        GEDpats(:,compi,frexi) = this_pattern;

        % Compute the joint network's broadband activation time series
        GEDts(compi,:,frexi) = this_evec'*broadData;

        % Total narrowband variance of this joint component
        total_contribution(compi,frexi) = ...
            real(this_evec'*covS*this_evec);

        % Participant-specific filter partitions, patterns, time series,
        % and diagonal intra-participant covariance contributions
        for participant_i = 1:nparticipants

            idx_p = channel_indices(participant_i,1): ...
                channel_indices(participant_i,2);

            evec_p = this_evec(idx_p);

            participant_evecs(:,compi,frexi,participant_i) = evec_p;
            participant_pats(:,compi,frexi,participant_i) = ...
                this_pattern(idx_p);
            participant_ts(compi,:,frexi,participant_i) = ...
                evec_p'*broadData(idx_p,:);

            covS_pp = covS(idx_p,idx_p);

            intra_contribution(compi,frexi,participant_i) = ...
                real(evec_p'*covS_pp*evec_p);
        end

        % Unique pairwise off-diagonal covariance contributions
        for pair_i = 1:npairs

            participant_p = pairs(pair_i,1);
            participant_q = pairs(pair_i,2);

            idx_p = channel_indices(participant_p,1): ...
                channel_indices(participant_p,2);
            idx_q = channel_indices(participant_q,1): ...
                channel_indices(participant_q,2);

            evec_p = this_evec(idx_p);
            evec_q = this_evec(idx_q);
            covS_pq = covS(idx_p,idx_q);

            % Factor two represents both S_pq and S_qp blocks.
            inter_contribution(compi,frexi,pair_i) = ...
                real(2*evec_p'*covS_pq*evec_q);
        end

        % Verify exact blockwise reconstruction of total variance
        reconstructed_total = ...
            sum(intra_contribution(compi,frexi,:),3) + ...
            sum(inter_contribution(compi,frexi,:),3);

        decomposition_error(compi,frexi) = ...
            abs(total_contribution(compi,frexi)-reconstructed_total) / ...
            max(abs(total_contribution(compi,frexi)),eps);

    end

end


%% Compute aggregate contributions and contribution ratios

intra_total = sum(intra_contribution,3);
inter_total = sum(inter_contribution,3);

intra_ratio = nan(size(intra_contribution));
inter_ratio = nan(size(inter_contribution));
intra_total_ratio = nan(size(intra_total));
inter_total_ratio = nan(size(inter_total));

valid_total = abs(total_contribution) > eps;

for participant_i = 1:nparticipants
    this_contribution = intra_contribution(:,:,participant_i);
    this_ratio = nan(size(this_contribution));
    this_ratio(valid_total) = ...
        this_contribution(valid_total)./ ...
        total_contribution(valid_total);
    intra_ratio(:,:,participant_i) = this_ratio;
end

for pair_i = 1:npairs
    this_contribution = inter_contribution(:,:,pair_i);
    this_ratio = nan(size(this_contribution));
    this_ratio(valid_total) = ...
        this_contribution(valid_total)./ ...
        total_contribution(valid_total);
    inter_ratio(:,:,pair_i) = this_ratio;
end

intra_total_ratio(valid_total) = ...
    intra_total(valid_total)./total_contribution(valid_total);
inter_total_ratio(valid_total) = ...
    inter_total(valid_total)./total_contribution(valid_total);


%% Assign outputs to structure

% Core outputs corresponding to FREQNESS_NetworkEstimation
HYPER.evals = GEDevals;
HYPER.evecs = GEDevecs;
HYPER.pats  = GEDpats;
HYPER.ts    = GEDts;
HYPER.frex  = frex;
HYPER.fwhm  = fwidth_all;
HYPER.srate = srate;

% Participant-specific partitions of the joint solution
HYPER.participant_evecs = participant_evecs;
HYPER.participant_pats  = participant_pats;
HYPER.participant_ts    = participant_ts;

% Static variance decomposition
HYPER.total = total_contribution;
HYPER.intra = intra_contribution;
HYPER.inter = inter_contribution;
HYPER.intra_total = intra_total;
HYPER.inter_total = inter_total;
HYPER.intra_ratio = intra_ratio;
HYPER.inter_ratio = inter_ratio;
HYPER.intra_total_ratio = intra_total_ratio;
HYPER.inter_total_ratio = inter_total_ratio;
HYPER.pairs = pairs;
HYPER.decomposition_error = decomposition_error;

% Participant and stacking metadata
HYPER.nparticipants = nparticipants;
HYPER.nchannels = nchannels;
HYPER.njoint_channels = njoint_channels;
HYPER.channel_indices = channel_indices;

% Analysis metadata used by second-order functions
HYPER.duration = pnts2keep/srate;
HYPER.bad_segments = idx2remove;
HYPER.regularisation = regularisation;

max_error = max(decomposition_error(:));

fprintf(['HYPERNESS Network Estimation complete. Maximum relative ' ...
    'decomposition error: %.3g.\n'], max_error);

end


%% Helper Function: Parse Name-Value Pairs

function opts = parse_name_value_pairs(opts, varargin)

if mod(length(varargin),2) ~= 0
    error('Arguments must be given as name-value pairs.');
end

for i = 1:2:length(varargin)
    name = lower(char(varargin{i}));

    if isfield(opts,name)
        opts.(name) = varargin{i+1};
    else
        error(['Unrecognized argument: ',name]);
    end
end

end

%%
