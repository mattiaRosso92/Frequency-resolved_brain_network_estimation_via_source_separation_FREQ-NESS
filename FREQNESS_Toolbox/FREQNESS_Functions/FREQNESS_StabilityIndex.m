function STAB = FREQNESS_StabilityIndex(FREQ,varargin)

% ========================================================================
%
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you use this toolbox, please cite:
%
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain
%  Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%  This function computes the Stability Index of frequency-resolved network
%  time series. For each selected network, the corresponding broadband
%  component time series is filtered at its own FREQ.frex frequency using
%  the FREQ.fwhm value stored by FREQNESS_NetworkEstimation. Instantaneous
%  frequency is obtained from the first derivative of the unwrapped analytic
%  phase, smoothed with a moving median, and summarized by its standard
%  deviation over time.
%
%  The Stability Index is expressed in Hz. Lower values indicate more stable
%  oscillatory frequency over time; higher values indicate wider frequency
%  fluctuations. The calculation follows Rosso et al. (2021), without the
%  behavioral-correlation analyses of the original study.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENT:
% ------------------------------------------------------------------------
%
%  - FREQ : structure with fields
%           • FREQ.ts    -> broadband component time series
%                           [nComp x nTime x nFrex x (nSubs)]
%           • FREQ.frex  -> analyzed frequency vector [nFrex x 1]
%           • FREQ.fwhm  -> filter FWHM vector [nFrex x 1]
%           • FREQ.srate -> sampling rate in Hz
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
%  - 'fwhm' : Gaussian-filter FWHM in Hz. A scalar applies one width to all
%             selected frequencies; a vector supplies one width per
%             selected frequency. Default: the corresponding FREQ.fwhm
%             values stored by FREQNESS_NetworkEstimation.
%
%  - 'median_window' : moving-median duration in seconds. The window is
%                      converted to samples using FREQ.srate. Default: .4.
%
%  - 'time_window' : two-element [start end] interval in seconds over which
%                    the Stability Index and mean frequency are calculated.
%                    Filtering, analytic-signal estimation, phase
%                    differentiation, and moving-median smoothing are always
%                    performed over the complete FREQ.ts duration. Default:
%                    the complete available duration.
%
%  - 'keep_frequency_sliding' : logical flag. When true, the raw and
%                               median-smoothed instantaneous-frequency
%                               time series are retained. Default: false.
%
%  - 'plot_avg' : logical flag. When true, group-average Stability Index
%                 curves are plotted for every selected component.
%                 Default: false.
%
%  - 'plot_all' : logical flag. When true, participant-level Stability
%                 Index curves are plotted. Default: false.
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%
%  - STAB : structure with fields
%
%       • STAB.stability_index -> frequency SD in Hz
%                                  [nComp x nFrex x nSubs]
%       • STAB.mean_frequency -> mean instantaneous frequency in Hz
%                                  [nComp x nFrex x nSubs]
%       • STAB.time           -> instantaneous-frequency time vector
%       • STAB.frex           -> frequencies included in the analysis
%       • STAB.fwhm           -> corresponding filter FWHM values
%       • STAB.comps          -> analyzed component indices
%       • STAB.srate         -> sampling rate in Hz
%       • STAB.median_window -> effective median duration in seconds
%       • STAB.median_samples -> moving-median window in samples
%       • STAB.time_window    -> interval used for summary calculation
%       • STAB.analysis_mask  -> samples included in summary calculation
%
%    When keep_frequency_sliding is true:
%
%       • STAB.frequency_sliding
%       • STAB.frequency_sliding_raw
%         [nComp x (nTime-1) x nFrex x nSubs]
%
% ------------------------------------------------------------------------
%  REFERENCES:
% ------------------------------------------------------------------------
%
%  - Rosso, M., Leman, M., & Moumdjian, L. (2021). Neural Entrainment
%    Meets Behavior: The Stability Index as a Neural Outcome Measure of
%    Auditory-Motor Coupling. Frontiers in Human Neuroscience, 15, 668918.
%    https://doi.org/10.3389/fnhum.2021.668918
%
%  - Cohen, M. X. (2014). Fluctuations in Oscillation Frequency Control
%    Spike Timing and Coordinate Neural Networks. Journal of Neuroscience,
%    34(27), 8988-8998.
%    https://doi.org/10.1523/JNEUROSCI.0261-14.2014
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso, Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 12/10/2026
%
% ========================================================================

%% Check input structure

if ~isstruct(FREQ)
    error('Input must be a FREQ structure.');
end

required_fields = {'ts','frex','fwhm','srate'};
for fieldi = 1:numel(required_fields)
    if ~isfield(FREQ,required_fields{fieldi}) || isempty(FREQ.(required_fields{fieldi}))
        error('FREQ.%s is missing or empty.',required_fields{fieldi});
    end
end

ts = FREQ.ts;
frex = FREQ.frex(:);
fwhm = FREQ.fwhm(:);
srate = FREQ.srate;

if ~isnumeric(ts) || ~isreal(ts) || isempty(ts) || ndims(ts)>4 || ...
        any(~isfinite(ts(:)))
    error(['FREQ.ts must be a finite real numeric array with dimensions ' ...
        '[nComp x nTime x nFrex x (nSubs)].']);
end
if size(ts,2)<3
    error('FREQ.ts must contain at least three time samples.');
end
if ~isnumeric(srate) || ~isscalar(srate) || ~isfinite(srate) || srate<=0
    error('FREQ.srate must be one positive finite scalar.');
end

ncomps_all = size(ts,1);
ntime = size(ts,2);
nfrex_all = size(ts,3);
nsubs = size(ts,4);

if ~isnumeric(frex) || ~isreal(frex) || numel(frex)~=nfrex_all || ...
        any(~isfinite(frex)) || any(frex<=0) || any(frex>=srate/2)
    error(['FREQ.frex must contain one positive finite frequency below the ' ...
        'Nyquist frequency for every frequency in FREQ.ts.']);
end
if ~isnumeric(fwhm) || ~isreal(fwhm) || numel(fwhm)~=nfrex_all || ...
        any(~isfinite(fwhm)) || any(fwhm<=0)
    error('FREQ.fwhm must contain one positive finite value per frequency.');
end

%% Handle optional arguments

opts = struct( ...
    'which_comp', [], ...
    'frex2model', [], ...
    'fwhm', [], ...
    'median_window', .4, ...
    'time_window', [], ...
    'keep_frequency_sliding', false, ...
    'plot_avg', false, ...
    'plot_all', false);
opts = parse_name_value_pairs(opts,varargin{:});

which_comp = opts.which_comp;
frex2model = opts.frex2model;
fwhm_override = opts.fwhm;
median_window = opts.median_window;
time_window = opts.time_window;
keep_frequency_sliding = opts.keep_frequency_sliding;
plot_avg = opts.plot_avg;
plot_all = opts.plot_all;

if isempty(which_comp)
    which_comp = 1;
    fprintf('\nFREQNESS StabilityIndex: no components specified. Defaulting to component #1.\n');
end
if ~isnumeric(which_comp) || ~isreal(which_comp) || ~isvector(which_comp) || ...
        any(~isfinite(which_comp)) || any(which_comp~=round(which_comp))
    error('which_comp must contain finite integer component indices.');
end
which_comp = unique(which_comp(:)','stable');
if any(which_comp<1) || any(which_comp>ncomps_all)
    error('which_comp must contain indices between 1 and %d.',ncomps_all);
end

if ~isnumeric(median_window) || ~isreal(median_window) || ...
        ~isscalar(median_window) || ~isfinite(median_window) || median_window<=0
    error('median_window must be one positive finite duration in seconds.');
end
median_samples = max(1,round(median_window*srate));
if median_samples>ntime-1
    error(['median_window corresponds to %d samples, but instantaneous ' ...
        'frequency contains only %d samples.'],median_samples,ntime-1);
end
median_window_effective = median_samples/srate;

validate_logical_scalar(keep_frequency_sliding,'keep_frequency_sliding');
validate_logical_scalar(plot_avg,'plot_avg');
validate_logical_scalar(plot_all,'plot_all');

[idx_frex2model,frex2model_actual] = map_frequencies(frex,frex2model);
fwhm2model = resolve_fwhm(fwhm,idx_frex2model,fwhm_override);
ncomps = numel(which_comp);
nfrex = numel(idx_frex2model);

%% Initialize outputs

stability_index = nan(ncomps,nfrex,nsubs);
mean_frequency = nan(ncomps,nfrex,nsubs);
time = (1:ntime-1)/srate;
[time_window,analysis_mask] = resolve_time_window(time_window,time);

if nfrex>1
    warning(['Stability Index is expressed in absolute Hz. Comparisons across ' ...
        'network frequencies should be interpreted cautiously because the ' ...
        'corresponding Gaussian FWHM values and expected absolute frequency ' ...
        'variability may differ. Within-frequency comparisons across ' ...
        'participants or conditions are more directly interpretable.']);
end
if time_window(1)<=time(1) || time_window(2)>=time(end)
    warning(['The Stability Index window includes the beginning or end of ' ...
        'FREQ.ts. FFT-based filtering and analytic-signal estimation may be ' ...
        'less reliable near the recording boundaries. Consider specifying ' ...
        'an interior time_window when edge effects could materially influence ' ...
        'the analysis.']);
end

if keep_frequency_sliding
    frequency_sliding = nan(ncomps,ntime-1,nfrex,nsubs);
    frequency_sliding_raw = nan(ncomps,ntime-1,nfrex,nsubs);
end

fprintf(['\nFREQNESS StabilityIndex: computing frequency stability for ' ...
    '%d component(s), %d frequencies, and %d participants.\n'], ...
    ncomps,nfrex,nsubs);
fprintf('Moving-median window: %d samples (%.4g seconds).\n', ...
    median_samples,median_window_effective);

%% Compute Stability Index

for subi = 1:nsubs

    fprintf('FREQNESS StabilityIndex: participant #%d/%d.\n',subi,nsubs);

    for frexi = 1:nfrex

        idx_freq = idx_frex2model(frexi);
        this_ts = reshape(ts(which_comp,:,idx_freq,subi),ncomps,ntime);

        filt_ts = filterFGx( ...
            this_ts,srate,frex(idx_freq),fwhm2model(frexi),false);
        analytic_signal = FREQNESS_AnalyticSignal(filt_ts.').';
        phase = unwrap(angle(analytic_signal),[],2);
        this_frequency_raw = srate*diff(phase,1,2)/(2*pi);
        this_frequency = movmedian(this_frequency_raw,median_samples,2);

        stability_index(:,frexi,subi) = ...
            std(this_frequency(:,analysis_mask),0,2);
        mean_frequency(:,frexi,subi) = ...
            mean(this_frequency(:,analysis_mask),2);

        if keep_frequency_sliding
            frequency_sliding(:,:,frexi,subi) = this_frequency;
            frequency_sliding_raw(:,:,frexi,subi) = this_frequency_raw;
        end

    end

end

%% Pack output structure

STAB = [];
STAB.stability_index = stability_index;
STAB.mean_frequency = mean_frequency;
STAB.time = time;
STAB.frex = frex2model_actual;
STAB.fwhm = fwhm2model;
STAB.comps = which_comp;
STAB.srate = srate;
STAB.median_window = median_window_effective;
STAB.median_samples = median_samples;
STAB.time_window = time_window;
STAB.analysis_mask = analysis_mask;

if keep_frequency_sliding
    STAB.frequency_sliding = frequency_sliding;
    STAB.frequency_sliding_raw = frequency_sliding_raw;
end

%% Visualize group-level Stability Index

if plot_avg

    for compi = 1:ncomps

        this_index = reshape(stability_index(compi,:,:),nfrex,nsubs);
        group_mean = mean(this_index,2,'omitnan');
        if nsubs>1
            group_sem = std(this_index,0,2,'omitnan')/sqrt(nsubs);
        else
            group_sem = zeros(nfrex,1);
        end

        figure('Name','FREQNESS StabilityIndex','Color','w'); clf
        errorbar(frex2model_actual,group_mean,group_sem,'o-', ...
            'Color',[.15 .35 .65],'MarkerFaceColor',[.15 .35 .65], ...
            'LineWidth',1.5,'CapSize',8)
        xlabel('Network frequency (Hz)')
        ylabel('Stability Index (Hz)')
        title(sprintf('Network Stability Landscape - component #%d', ...
            which_comp(compi)))
        grid on

    end

end


%% Optional participant-level visualizations

if plot_all

    for subi = 1:nsubs
        for compi = 1:ncomps

            this_index = reshape(stability_index(compi,:,subi),nfrex,1);
            figure('Name','FREQNESS StabilityIndex: participant','Color','w'); clf
            plot(frex2model_actual,this_index,'o-', ...
                'Color',[.15 .35 .65],'MarkerFaceColor',[.15 .35 .65], ...
                'LineWidth',1.5)
            xlabel('Network frequency (Hz)')
            ylabel('Stability Index (Hz)')
            title(sprintf('Participant #%d - component #%d', ...
                subi,which_comp(compi)))
            grid on

        end
    end

end

end


%% Helper Function: Resolve Filter Widths
function fwhm2model = resolve_fwhm(fwhm,idx_frex2model,fwhm_override)

nfrex = numel(idx_frex2model);
if isempty(fwhm_override)
    fwhm2model = fwhm(idx_frex2model);
elseif ~isnumeric(fwhm_override) || ~isreal(fwhm_override) || ...
        ~isvector(fwhm_override) || any(~isfinite(fwhm_override)) || ...
        any(fwhm_override<=0)
    error('fwhm must contain positive finite values in Hz.');
elseif isscalar(fwhm_override)
    fwhm2model = repmat(fwhm_override,nfrex,1);
elseif numel(fwhm_override)==nfrex
    fwhm2model = fwhm_override(:);
else
    error(['fwhm must be empty, one scalar, or contain one value per ' ...
        'selected frequency (%d values).'],nfrex);
end
end


%% Helper Function: Resolve Analysis Time Window
function [time_window,analysis_mask] = resolve_time_window(time_window,time)

duration = time(end)+time(1);
if isempty(time_window)
    time_window = [0 duration];
elseif ~isnumeric(time_window) || ~isreal(time_window) || ...
        ~isvector(time_window) || numel(time_window)~=2 || ...
        any(~isfinite(time_window)) || time_window(1)<0 || ...
        time_window(1)>=time_window(2) || time_window(2)>duration
    error(['time_window must contain two finite increasing values in seconds ' ...
        'within the FREQ.ts duration [0, %.4g].'],duration);
else
    time_window = time_window(:)';
end

analysis_mask = time>=time_window(1) & time<=time_window(2);
if sum(analysis_mask)<2
    error('time_window must contain at least two instantaneous-frequency samples.');
end
end


%% Helper Function: Validate Logical Scalars
function validate_logical_scalar(value,value_name)
if ~(islogical(value) || isnumeric(value)) || ~isscalar(value) || ...
        ~isfinite(value) || ~ismember(value,[0 1])
    error('%s must be a logical scalar.',value_name);
end
end


%% Helper Function: Map Requested Frequencies
function [idx_frex2model,frex2model_actual] = map_frequencies(frex,frex2model)

nfrex = numel(frex);

if isempty(frex2model)

    idx_frex2model = 1:nfrex;

elseif ~isnumeric(frex2model) || ~isreal(frex2model) || ...
        ~isvector(frex2model) || any(~isfinite(frex2model))

    error('frex2model must be empty, one frequency, or a two-element range in Hz.');

elseif isscalar(frex2model)

    [~,idx_frex2model] = min(abs(frex-frex2model));

    if ~ismembertol(frex2model,frex,1e-10)
        warning(['Requested frex2model (%.3f Hz) is not present in FREQ.frex. ' ...
            'Using the closest available frequency: %.3f Hz.'], ...
            frex2model,frex(idx_frex2model));
    end

elseif numel(frex2model)==2

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


%% Helper Function: Parse Name-Value Pairs
function opts = parse_name_value_pairs(opts,varargin)
if mod(length(varargin),2)~=0
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

