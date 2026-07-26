function IND = FREQNESS_InducedResponses(FREQ, events, epoch_window, baseline_window, varargin)

% ========================================================================
%
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you find this function useful, please cite the first FREQNESS paper:
%
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain
%  Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%  This function computes event-related induced oscillatory responses from
%  the broadband component time series produced by
%  FREQNESS_NetworkEstimation.
%
%  For each selected FREQNESS network frequency, the corresponding
%  component time series is filtered continuously at that same frequency
%  using FREQ.frex and FREQ.fwhm. Instantaneous power is obtained from the
%  squared magnitude of the Hilbert transform. Power is then segmented
%  around the supplied events, normalized trial-by-trial relative to the
%  pre-event baseline in decibels, and averaged across trials.
%
%  Filtering is performed before epoching to reduce trial-edge artifacts.
%  The resulting response retains the frequency-resolved organization of
%  the input FREQ structure.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - FREQ : structure with fields
%           • FREQ.ts    -> broadband component time series
%                           [nComp x nTime x nFrex x (nSubs)]
%           • FREQ.frex  -> analyzed frequency vector [nFrex x 1]
%           • FREQ.fwhm  -> filter FWHM vector [nFrex x 1]
%           • FREQ.srate -> sampling rate in Hz
%
%  - events : event onset samples in the continuous FREQ.ts time series.
%             For one participant, provide a numeric vector. For multiple
%             participants, provide a cell array with one numeric vector
%             per participant: events{subi}.
%
%  - epoch_window : two-element time window around each event, expressed
%                   in seconds, e.g. [-0.1 3.4]. Both endpoints are
%                   included in the output epoch.
%
%  - baseline_window : two-element baseline window in seconds relative to
%                      event onset, e.g. [-0.1 0]. The first endpoint is
%                      included and the second endpoint is excluded, so
%                      that a baseline ending at 0 does not include the
%                      event-onset sample.
%
% ------------------------------------------------------------------------
%  OPTIONAL NAME–VALUE PAIRS:
% ------------------------------------------------------------------------
%
%  - 'which_comp' : component index or vector of component indices to
%                   analyze. Default: 1.
%
%  - 'frex2model' : scalar frequency or two-element frequency range in Hz.
%                   The closest available frequencies in FREQ.frex are
%                   used. Default: all frequencies in FREQ.frex.
%
%  - 'keep_trials' : logical flag. When true, single-trial induced
%                    responses are retained in IND.power_trials.
%                    Default: false.
%
%  - 'scale_max' : logical flag. When true, IND.power_scaled contains an
%                  additional replication-oriented version of IND.power,
%                  divided by its maximum separately for each component
%                  and participant. IND.power always remains in dB.
%                  Default: false.
%
%  - 'plot_avg' : logical flag. When true, the function visualizes the
%                 grand-average time-frequency response for each selected
%                 component. Default: false.
%
%  - 'plot_all' : logical flag. When true, participant-level time-frequency
%                 maps are plotted.
%                 Default: false.
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%
%  - IND : structure with fields
%
%       • IND.power           -> trial-averaged induced power in dB
%                                [nComp x nEpochTime x nFrex x nSubs]
%       • IND.time            -> epoch time vector in seconds
%       • IND.frex            -> frequencies included in the analysis
%       • IND.fwhm            -> corresponding filter FWHM values
%       • IND.comps           -> analyzed component indices
%       • IND.events          -> valid event samples used per participant
%       • IND.ntrials         -> number of valid trials per participant
%       • IND.epoch_window    -> requested epoch window in seconds
%       • IND.baseline_window -> requested baseline window in seconds
%       • IND.srate           -> sampling rate in Hz
%
%    When keep_trials is true:
%
%       • IND.power_trials{subi}
%         [nComp x nEpochTime x nFrex x nTrials]
%
%    When scale_max is true:
%
%       • IND.power_scaled
%         [nComp x nEpochTime x nFrex x nSubs]
%
% ------------------------------------------------------------------------
%  REFERENCES:
% ------------------------------------------------------------------------
%
%  - Malvaso et al. (2025).
%    FREQ-NESS reveals age-related differences in frequency-resolved brain
%    networks during auditory recognition and resting state.
%    bioRxiv. https://doi.org/10.1101/2025.06.30.662094
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso, Chiara Malvaso, Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 26/07/2026
%
% ========================================================================

%% Check mandatory inputs

if nargin < 2 || isempty(events)
    error('events is required and must contain event-onset sample indices.');
end

if nargin < 3 || isempty(epoch_window)
    error('epoch_window is required and must be expressed in seconds.');
end

if nargin < 4 || isempty(baseline_window)
    error('baseline_window is required and must be expressed in seconds.');
end

if ~isstruct(FREQ)
    error('Input must be a FREQ structure.');
end

required_fields = {'ts','frex','fwhm','srate'};
for fieldi = 1:numel(required_fields)
    if ~isfield(FREQ,required_fields{fieldi}) || isempty(FREQ.(required_fields{fieldi}))
        error('FREQ.%s is missing or empty.',required_fields{fieldi});
    end
end

validate_time_window(epoch_window,'epoch_window');
validate_time_window(baseline_window,'baseline_window');

if baseline_window(1) < epoch_window(1) || baseline_window(2) > epoch_window(2)
    error('baseline_window must fall entirely within epoch_window.');
end

if ~isnumeric(FREQ.srate) || ~isscalar(FREQ.srate) || ...
        ~isfinite(FREQ.srate) || FREQ.srate <= 0
    error('FREQ.srate must be one positive finite scalar.');
end

%% Handle optional arguments (name-value pairs)

opts = struct( ...
    'which_comp', [], ...
    'frex2model', [], ...
    'keep_trials', false, ...
    'scale_max', false, ...
    'plot_avg', false, ...
    'plot_all', false);

opts = parse_name_value_pairs(opts,varargin{:});

which_comp = opts.which_comp;
frex2model = opts.frex2model;
keep_trials = opts.keep_trials;
scale_max = opts.scale_max;
plot_avg = opts.plot_avg;
plot_all = opts.plot_all;

if isempty(which_comp)
    which_comp = 1;
    fprintf('\nFREQNESS InducedResponses: no components specified. Defaulting to component #1.\n');
end

validate_logical_scalar(keep_trials,'keep_trials');
validate_logical_scalar(scale_max,'scale_max');
validate_logical_scalar(plot_avg,'plot_avg');
validate_logical_scalar(plot_all,'plot_all');

%% Map and validate inputs from FREQ

ts = FREQ.ts;
frex = FREQ.frex(:);
fwhm = FREQ.fwhm(:);
srate = FREQ.srate;

ncomps_all = size(ts,1);
ntime_cont = size(ts,2);
nfrex_all = size(ts,3);
nsubs = size(ts,4);

if numel(frex) ~= nfrex_all
    error('Length of FREQ.frex does not match the 3rd dimension of FREQ.ts.');
end

if numel(fwhm) ~= nfrex_all
    error('FREQ.fwhm must contain one value per frequency in FREQ.frex.');
end

if any(~isfinite(frex)) || any(frex <= 0)
    error('FREQ.frex must contain positive finite frequencies.');
end

if any(~isfinite(fwhm)) || any(fwhm <= 0)
    error('FREQ.fwhm must contain positive finite filter widths.');
end

if ~isnumeric(which_comp) || ~isvector(which_comp) || ...
        any(~isfinite(which_comp)) || any(which_comp ~= round(which_comp))
    error('which_comp must contain finite integer component indices.');
end

which_comp = unique(which_comp(:)','stable');

if any(which_comp < 1) || any(which_comp > ncomps_all)
    error('which_comp must contain indices between 1 and %d.',ncomps_all);
end

[idx_frex2model,frex2model_actual] = map_frequencies(frex,frex2model);
fwhm2model = fwhm(idx_frex2model);

ncomps = numel(which_comp);
nfrex = numel(idx_frex2model);

%% Validate and organize event samples

events_all = normalize_events(events,nsubs,ntime_cont);

%% Define epoch and baseline samples

epoch_offsets = round(epoch_window(1)*srate):round(epoch_window(2)*srate);
time = epoch_offsets/srate;
nTimeEpoch = numel(time);

% Use a half-open interval [start,end) to exclude event onset from a
% prestimulus baseline ending at zero.
baseline_idx = time >= baseline_window(1) & time < baseline_window(2);

if ~any(baseline_idx)
    error(['baseline_window does not contain any samples at the current ' ...
        'sampling rate. Increase its duration or check FREQ.srate.']);
end

%% Initialize outputs

power_avg = nan(ncomps,nTimeEpoch,nfrex,nsubs);
events_valid = cell(nsubs,1);
ntrials = zeros(nsubs,1);

if keep_trials
    power_trials = cell(nsubs,1);
end

fprintf(['\nFREQNESS InducedResponses: computing induced responses for ' ...
    '%d component(s), %d frequencies, and %d participants.\n'], ...
    ncomps,nfrex,nsubs);

%% Compute induced oscillatory responses

for subi = 1:nsubs

    % Retain only events whose complete epoch lies within FREQ.ts
    this_events = events_all{subi};
    valid_event_idx = this_events + epoch_offsets(1) >= 1 & ...
        this_events + epoch_offsets(end) <= ntime_cont;

    events_valid{subi} = this_events(valid_event_idx);
    ntrials(subi) = numel(events_valid{subi});

    nremoved = sum(~valid_event_idx);
    if nremoved > 0
        warning(['Participant #%d: excluding %d event(s) because the ' ...
            'requested epoch extends outside FREQ.ts.'],subi,nremoved);
    end

    if ntrials(subi) == 0
        error('Participant #%d has no valid events after epoch-boundary checks.',subi);
    end

    this_trials = nan(ncomps,nTimeEpoch,nfrex,ntrials(subi));

    fprintf('FREQNESS InducedResponses: participant #%d (%d trials).\n', ...
        subi,ntrials(subi));

    for compi = 1:ncomps

        for frexi = 1:nfrex

            idx_freq = idx_frex2model(frexi);

            % Broadband component time series for the GED network estimated
            % at the current frequency
            this_ts = reshape(ts(which_comp(compi),:,idx_freq,subi),1,ntime_cont);

            if any(~isfinite(this_ts))
                error(['FREQ.ts contains non-finite values for participant ' ...
                    '#%d, component #%d, frequency %.3f Hz.'], ...
                    subi,which_comp(compi),frex(idx_freq));
            end

            % Matched-frequency Gaussian filtering, followed by Hilbert power
            filt_ts = filterFGx(this_ts,srate,frex(idx_freq),fwhm(idx_freq),0);
            power_ts = abs(hilbert(filt_ts.')).^2;
            power_ts = reshape(power_ts,1,ntime_cont);

            for triali = 1:ntrials(subi)

                idx_epoch = events_valid{subi}(triali) + epoch_offsets;
                this_power = power_ts(idx_epoch);
                baseline_power = mean(this_power(baseline_idx),'omitnan');

                if ~isfinite(baseline_power) || baseline_power <= 0
                    warning(['Participant #%d, trial #%d, component #%d, ' ...
                        'frequency %.3f Hz: invalid baseline power. ' ...
                        'The trial is retained as NaN.'], ...
                        subi,triali,which_comp(compi),frex(idx_freq));
                    continue
                end

                % Trial-wise baseline normalization in decibels
                this_trials(compi,:,frexi,triali) = ...
                    10*log10(this_power./baseline_power);

            end

        end

    end

    % Average power after trial-wise power and baseline normalization:
    % this preserves induced activity that is not phase-locked across trials.
    power_avg(:,:,:,subi) = mean(this_trials,4,'omitnan');

    if keep_trials
        power_trials{subi} = this_trials;
    end

end

%% Pack output structure

IND = [];
IND.power = power_avg;
IND.time = time;
IND.frex = frex2model_actual;
IND.fwhm = fwhm2model;
IND.comps = which_comp;
IND.events = events_valid;
IND.ntrials = ntrials;
IND.epoch_window = epoch_window;
IND.baseline_window = baseline_window;
IND.srate = srate;

if keep_trials
    IND.power_trials = power_trials;
end

%% Optional maximum scaling for replication-oriented analyses

if scale_max

    power_scaled = nan(size(power_avg));

    for subi = 1:nsubs
        for compi = 1:ncomps

            this_power = reshape(power_avg(compi,:,:,subi),nTimeEpoch,nfrex);
            max_power = max(this_power(:),[],'omitnan');

            if isfinite(max_power) && max_power ~= 0
                power_scaled(compi,:,:,subi) = this_power./max_power;
            else
                warning(['Participant #%d, component #%d: maximum scaling ' ...
                    'was not applied because the maximum power was invalid or zero.'], ...
                    subi,which_comp(compi));
            end

        end
    end

    IND.power_scaled = power_scaled;

end

%% Visualize group-level induced responses

if plot_avg

    avg_power = mean(power_avg,4,'omitnan');

    for compi = 1:ncomps

        figure('Name','FREQNESS InducedResponses','Color','w'); clf

        imagesc(time,1:nfrex,squeeze(avg_power(compi,:,:))');
        axis xy
        colormap(parula(256))
        cbar = colorbar;
        ylabel(cbar,'Power change (dB)')
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        title(sprintf('Induced responses - component #%d',which_comp(compi)))
        set_frequency_ticks(frex2model_actual);

    end

end

%% Optional participant-level visualizations

if plot_all

    for subi = 1:nsubs
        for compi = 1:ncomps

            figure('Name','FREQNESS InducedResponses: participant','Color','w'); clf

            imagesc(time,1:nfrex,squeeze(power_avg(compi,:,:,subi))');
            axis xy
            colormap(parula(256))
            cbar = colorbar;
            ylabel(cbar,'Power change (dB)')
            xlabel('Time (s)')
            ylabel('Frequency (Hz)')
            title(sprintf('Participant #%d - component #%d', ...
                subi,which_comp(compi)))
            set_frequency_ticks(frex2model_actual);

        end
    end

end

end


%% Helper Function: Validate Time Windows
function validate_time_window(time_window,window_name)
if ~isnumeric(time_window) || numel(time_window) ~= 2 || ...
        any(~isfinite(time_window)) || time_window(1) >= time_window(2)
    error('%s must contain two finite increasing values in seconds.',window_name);
end
end


%% Helper Function: Validate Logical Scalars
function validate_logical_scalar(value,value_name)
if ~(islogical(value) || isnumeric(value)) || ~isscalar(value) || ...
        ~isfinite(value) || ~ismember(value,[0 1])
    error('%s must be a logical scalar.',value_name);
end
end


%% Helper Function: Normalize Events
function events_all = normalize_events(events,nsubs,ntime_cont)

if isnumeric(events)

    if nsubs ~= 1
        error(['For multiple participants, events must be a cell array ' ...
            'containing one vector of event samples per participant.']);
    end

    events_all = {events};

elseif iscell(events)

    if numel(events) ~= nsubs
        error('events must contain one cell per participant (%d cells).',nsubs);
    end

    events_all = reshape(events,nsubs,1);

else
    error('events must be a numeric vector or a cell array of numeric vectors.');
end

for subi = 1:nsubs

    this_events = events_all{subi};

    if ~isnumeric(this_events) || ~isvector(this_events) || isempty(this_events) || ...
            any(~isfinite(this_events)) || any(this_events ~= round(this_events))
        error('events{%d} must contain finite integer sample indices.',subi);
    end

    this_events = this_events(:)';

    if any(this_events < 1) || any(this_events > ntime_cont)
        error('events{%d} contains samples outside the FREQ.ts time range.',subi);
    end

    if numel(unique(this_events)) ~= numel(this_events)
        warning('events{%d} contains duplicate samples. Duplicates were removed.',subi);
        this_events = unique(this_events,'stable');
    end

    events_all{subi} = this_events;

end

end


%% Helper Function: Map Requested Frequencies
function [idx_frex2model,frex2model_actual] = map_frequencies(frex,frex2model)

nfrex = numel(frex);

if isempty(frex2model)

    idx_frex2model = 1:nfrex;

elseif ~isnumeric(frex2model) || ~isvector(frex2model) || ...
        any(~isfinite(frex2model))

    error('frex2model must be empty, one frequency, or a two-element range in Hz.');

elseif isscalar(frex2model)

    [~,idx_frex2model] = min(abs(frex-frex2model));

    if ~ismembertol(frex2model,frex,1e-10)
        warning(['Requested frex2model (%.3f Hz) is not present in FREQ.frex. ' ...
            'Using the closest available frequency: %.3f Hz.'], ...
            frex2model,frex(idx_frex2model));
    end

elseif numel(frex2model) == 2

    frex2model = sort(frex2model(:),'ascend');

    [~,idx_first] = min(abs(frex-frex2model(1)));
    [~,idx_last] = min(abs(frex-frex2model(2)));

    idx_frex2model = min(idx_first,idx_last):max(idx_first,idx_last);

    if ~all(ismembertol(frex2model,frex,1e-10))
        warning(['One or both frex2model limits are not present in FREQ.frex. ' ...
            'Using the closest available range: %.3f to %.3f Hz.'], ...
            frex(idx_frex2model(1)),frex(idx_frex2model(end)));
    end

else
    error('frex2model must contain either one frequency or a two-element range.');
end

frex2model_actual = frex(idx_frex2model);

end


%% Helper Function: Set Frequency-Axis Ticks
function set_frequency_ticks(frex)

nfrex = numel(frex);
nticks = min(nfrex,12);
idx_ticks = unique(round(linspace(1,nfrex,nticks)));

yticks(idx_ticks)
yticklabels(arrayfun(@(x) sprintf('%.2f',x), ...
    frex(idx_ticks),'UniformOutput',false))

end


%% Helper Function: Parse Name-Value Pairs
function opts = parse_name_value_pairs(opts,varargin)
if mod(length(varargin),2) ~= 0
    error('Arguments must be given as name-value pairs.');
end
for i = 1:2:length(varargin)
    name = lower(varargin{i});
    if isfield(opts,name)
        opts.(name) = varargin{i+1};
    else
        error(['Unrecognized argument: ',name]);
    end
end
end

%%
