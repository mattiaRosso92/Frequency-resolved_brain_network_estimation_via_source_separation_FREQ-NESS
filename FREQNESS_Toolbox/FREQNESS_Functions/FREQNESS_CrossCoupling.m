function CFC = FREQNESS_CrossCoupling(FREQ, lfo_freq, varargin)

% ========================================================================
%
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you use this toolbox, please cite:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved Brain
%  Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%  This function computes cross-frequency coupling (CFC) between a low-
%  frequency oscillator (LFO) and higher-frequency carrier networks, using
%  the FREQNESS component time series stored in FREQ.ts.
%
%  Specifically, it quantifies how the instantaneous phase of one
%  component at a low frequency (the modulator) modulates the power of the
%  *same* component across the eigenspectrum (the carriers).
%
%  CFC is estimated via:
%      Phase-amplitude coupling (PAC): mean carrier power in
%      phase bins across one LFO cycle, with a sinusoidal fit (sineFit.m)
%      to estimate modulation amplitude, phase shift, etc.
%
%  The function operates on ONE component only (which_comp), which acts
%  both as the LFO modulator and as the carrier component across
%  frequencies.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - FREQ : structure with fields
%
%           • FREQ.frex   -> frequency vector [nFrex x 1]
%           • FREQ.ts     -> component time series
%                            [nComp x nTime x nFrex x (nSubs)]
%                            or [nComp x nTime x nFrex] for single-subject
%           • FREQ.srate  -> sampling rate (Hz) of the component time
%                            series in FREQ.ts. Used here to re-filter
%                            components around each center frequency.
%           • FREQ.fwhm   -> vector of full-width at half-maximum (FWHM)
%                            values for the bandpass / wavelet used at
%                            each frequency, [nFrex x 1]. CFC re-filters
%                            the time series at each frex using
%                            filterFGx(FREQ.srate, FREQ.frex, FREQ.fwhm).
%
%           (optional but very relevant for visualization)
%           • FREQ.pats   -> activation patterns
%                            [nVoxels x nComp x nFrex x (nSubs)]
%
%  - lfo_freq : scalar with the frequency (in Hz) of the low-frequency
%               oscillator to be used as phase modulator. The closest
%               value in FREQ.frex will be used.
%
% ------------------------------------------------------------------------
%  OPTIONAL NAME–VALUE PAIRS:
% ------------------------------------------------------------------------
%
%  - 'MNI'        : MNI coordinates for voxels (N x 3), same order as
%                   FREQ.pats. Used only for visualization.
%
%  - 'frex2model' : [fmin fmax] frequency range (Hz) within which carrier
%                   frequencies are considered. The LFO frequency is
%                   automatically excluded from carriers.
%                   Default: [min(FREQ.frex) max(FREQ.frex)].
%
%  - 'which_comp' : scalar index of the component to use as both modulator
%                   and carrier across frequencies. Default: 1.
%
%  - 'plot_all'   : logical flag. When true, produces per-subject CFC
%                   plots in addition to the group-level visualization.
%                   Default: false.
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%
%  CFC : structure with fields
%
%     • CFC.PAC_all      -> [nCarriers x nSubs x nBins] PAC histograms
%     • CFC.PAC_avg      -> [nCarriers x nBins] group average PAC
%     • CFC.sAmpl        -> [nCarriers x nSubs] sine-fit amplitude
%     • CFC.pShift       -> [nCarriers x nSubs] sine-fit phase shift
%     • CFC.dcOff        -> [nCarriers x nSubs] sine-fit DC offset
%     • CFC.mFrex        -> [nCarriers x nSubs] sine-fit modulation freq
%     • CFC.goodFit      -> [nCarriers x nSubs] sine-fit error / goodness
%     • CFC.carrier_frex -> [nCarriers x 1] carrier frequencies (Hz)
%     • CFC.lfo_freq     -> scalar, actual LFO frequency used (Hz)
%     • CFC.comp         -> scalar, component index used
%     • CFC.phase_edges  -> [1 x (nBins+1)] phase bin edges (rad)
%
%  Additionally, the function produces group-level figures summarizing:
%   - PAC as a function of phase and frequency (surface plot)
%   - Optional 3D brain plots of LFO and peak carrier networks if MNI
%     coordinates and FREQ.pats are provided.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso & Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%
% ========================================================================

%% ----------------------- Parse and validate inputs ----------------------

% Handle optional arguments (name–value pairs) via helper
opts = struct( ...
    'mni',        [], ...
    'frex2model', [min(FREQ.frex) max(FREQ.frex)], ...
    'which_comp', 1, ...
    'plot_all',   false);

opts = parse_name_value_pairs(opts, varargin{:});

MNI        = opts.mni;
frex2model = opts.frex2model;
which_comp = opts.which_comp;
plot_all   = opts.plot_all;

% Check required fields
if ~isfield(FREQ,'frex')
    error('FREQNESS_CrossCoupling:MissingField', ...
        'FREQ.frex is required.');
end
if ~isfield(FREQ,'ts')
    error('FREQNESS_CrossCoupling:MissingField', ...
        'FREQ.ts is required (nComp x nTime x nFrex x (sub)).');
end
if ~isfield(FREQ,'srate')
    error('FREQNESS_CrossCoupling:MissingField', ...
        'FREQ.srate is required for re-filtering the network time series.');
end
if ~isfield(FREQ,'fwhm')
    error('FREQNESS_CrossCoupling:MissingField', ...
        'FREQ.fwhm is required for re-filtering the network time series and must contain one FWHM value per frequency in FREQ.frex.');
end

% Assignment to local variable
frex = FREQ.frex(:);
nfrex = numel(frex);

% Bring FREQ.ts to 4D: comp x time x frex x sub
ts = FREQ.ts;
switch ndims(ts)
    case 3
        [ncomps, npnts, nfrex_ts] = size(ts);
        nsubs = 1;
        ts    = reshape(ts, [ncomps, npnts, nfrex_ts, 1]);
    case 4
        [ncomps, npnts, nfrex_ts, nsubs] = size(ts);
    otherwise
        error('FREQ.ts must be a 3D or 4D array.');
end

if nfrex_ts ~= nfrex
    error('FREQ.ts and FREQ.frex have incompatible frequency dimensions.');
end

% Check FWHM dimensionality
fwhm = FREQ.fwhm(:);
if numel(fwhm) ~= nfrex
    error('FREQ.fwhm must have one entry per frequency in FREQ.frex.');
end

% Component sanity check
if which_comp < 1 || which_comp > ncomps
    error('which_comp (%d) is out of bounds. Available components: 1..%d.', ...
        which_comp, ncomps);
end

%% ---------------------- Define LFO and carriers -------------------------

% Find index of LFO in FREQ.frex
[~, idx_lfo] = min(abs(frex - lfo_freq));
lfo_freq_actual = frex(idx_lfo);

% Frequency mask for carriers
mask_frex = frex >= frex2model(1) & frex <= frex2model(2);
mask_frex(1:idx_lfo) = false; % exclude frequencies <= LFO itself

idx_carriers = find(mask_frex);
if isempty(idx_carriers)
    error('No carrier frequencies found in the specified frex2model range.');
end

% Assign carrier frequencies
carrier_frex = frex(idx_carriers);
ncars = numel(carrier_frex);

%% ------------------ Compute LFO phase for all subjects ------------------

% Re-filter the LFO component time series around lfo_freq_actual using the
% same filter width used during FREQNESS_NetworkEstimation, then compute 
% phase via Hilbert transform.

lfo_phase_all = zeros(npnts, nsubs);

for subi = 1:nsubs
    % Raw LFO component time series (1 x nTime for filterFGx)
    lfo_ts = squeeze(ts(which_comp,:,idx_lfo, subi));         % [1 x nTime]
    % Re-filter around the LFO frequency
    lfo_ts = filterFGx(lfo_ts, FREQ.srate, lfo_freq_actual, fwhm(idx_lfo), 0);
    % Hilbert operates along columns (time x chans), hence transpose back
    lfo_phase_all(:,subi) = angle( hilbert(lfo_ts.') );       % [nTime x 1], in [-pi, pi]
end

% Define phase binning
nbins        = 37; % as in the paper implementation
phase_edges  = linspace(-pi, pi, nbins+1);
phase_centers = phase_edges(1:end-1) + diff(phase_edges)/2;

%% --------------------- PAC across carriers & subjects -------------------

PAC_all = nan(ncars, nsubs, nbins);     % [nCarriers x nSubs x nBins]

for carri = 1:ncars
    
    idx_freq = idx_carriers(carri); % actual frequency index in FREQ.frex
    
    for subi = 1:nsubs
        
        % Carrier time series: re-filter around current carrier frequency
        carr_ts = squeeze(ts(which_comp, :, idx_freq, subi));      % [1 x nTime]
        carr_ts = filterFGx(carr_ts, FREQ.srate, frex(idx_freq), fwhm(idx_freq), 0);
        carr_pow = abs( hilbert(carr_ts.')).^2;                    % [nTime x 1] power
        
        % Phase of the LFO (same component, lfo_freq)
        phase = lfo_phase_all(:,subi);                             % [nTime x 1]
        
        % --- PAC: mean power per phase bin ---
        for bini = 1:nbins
            idx_phase = phase > phase_edges(bini) & phase <= phase_edges(bini+1);
            if any(idx_phase)
                PAC_all(carri,subi,bini) = mean(carr_pow(idx_phase), 'omitnan');
            end
        end
        
    end % subjects
    
end % carriers

% Group-level averages: [nCarriers x nBins]
PAC_avg = squeeze(mean(PAC_all, 2, 'omitnan'));

%% ------------------ Sinusoidal fit of PAC (sineFit.m) -------------------

smth = 1;  % smoothing factor for PAC curves before fitting

sAmpl   = nan(ncars, nsubs);  % amplitude
mFrex   = nan(ncars, nsubs);  % modulation frequency
pShift  = nan(ncars, nsubs);  % phase shift
dcOff   = nan(ncars, nsubs);  % DC offset
goodFit = nan(ncars, nsubs);  % MSE / goodness (depends on sineFit)

x2fit = 1:nbins;

for carri = 1:ncars
    for subi = 1:nsubs
        
        y = squeeze(PAC_all(carri,subi,:));   % [nBins x 1]
        if all(isnan(y))
            continue
        end
        
        % Smooth and fit sine: params = [dcOff, ampl, freq, phase, mse]
        y_sm = smooth(y, smth)';             % row vector
        
        params = sineFit(x2fit, y_sm, 0);        % no plot
        
        dcOff(carri,subi)   = params(1);
        sAmpl(carri,subi)   = params(2);
        mFrex(carri,subi)   = params(3);
        pShift(carri,subi)  = params(4);
        goodFit(carri,subi) = params(5);
        
    end
end

% Replace zeros with NaN (mirroring the paper script)
sAmpl(  sAmpl   == 0) = NaN;
pShift( pShift  == 0) = NaN;
dcOff(  dcOff   == 0) = NaN;

% Compute average amplitude modulation index, for plotting
sAmpl_avg = squeeze( mean(sAmpl,2,'omitnan') );

%% ---------------------------- Pack output -------------------------------

CFC.PAC_all       = PAC_all;
CFC.PAC_avg       = PAC_avg;

CFC.sAmpl         = sAmpl;
CFC.pShift        = pShift;
CFC.dcOff         = dcOff;
CFC.mFrex         = mFrex;
CFC.goodFit       = goodFit;

CFC.carrier_frex  = carrier_frex;
CFC.lfo_freq      = lfo_freq_actual;
CFC.comp          = which_comp;
CFC.phase_edges   = phase_edges;
CFC.phase_centers = phase_centers;

%% ----------------------- Group-level visualizations ---------------------

figure('Name','FREQNESS CrossCoupling: PAC','Color','w'); clf
cmap = parula(ncars);  % one color per carrier frequency

% PAC histograms (power over phase) for all carrier frequencies
subplot(2,1,1); hold on
for carri = 1:ncars
    plot(phase_centers, ...
         PAC_avg(carri,:), ...
         'Color', cmap(carri,:), ...
         'LineWidth', 1.6);
end
xlabel('Phase (rad)', 'FontSize', 13, 'FontWeight','bold');
ylabel('Modulation amplitude (a.u.)',       'FontSize', 13, 'FontWeight','bold');
title(sprintf('PAC across frequencies (Group Average): LFO %.3g Hz (comp %d)', ...
    lfo_freq_actual, which_comp), ...
    'FontSize', 14, 'FontWeight','bold');
set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
grid on; grid minor;
% Legend with carrier frequencies (on the second subplot)
leg_entries = cell(ncars,1);
for carri = 1:ncars
    leg_entries{carri} = sprintf('%.1f Hz', carrier_frex(carri));
end
legend(leg_entries, 'Location','eastoutside');

% Sine-fit amplitude (power modulation) across carrier frequencies
subplot(2,1,2); hold on
col_line = cmap(1,:);  % pick the first parula color for this line
plot(carrier_frex, sAmpl_avg, 'Color', col_line, 'LineWidth', 1.8);
xlabel('Carrier frequency (Hz)', 'FontSize', 13, 'FontWeight','bold');
ylabel('Power (a.u.)', 'FontSize', 13, 'FontWeight','bold');
title(sprintf('Power modulation: LFO %.3g Hz (comp %d)', ...
    lfo_freq_actual, which_comp), ...
    'FontSize', 14, 'FontWeight','bold');
set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
grid on; grid minor;


%% -------- Optional: brain plots if MNI & activation patterns provided ----

if ~isempty(MNI) && isfield(FREQ,'pats')
    
    pats = FREQ.pats;
    switch ndims(pats)
        case 3
            [nvoxs, ncomps_pat, nfrex_pat] = size(pats);
            nSubs_p = 1;
            pats    = reshape(pats,[nvoxs ncomps_pat nfrex_pat 1]);
        case 4
            [nvoxs, ncomps_pat, nfrex_pat, nSubs_p] = size(pats);
        otherwise
            warning('FREQ.pats has unexpected dimensionality. Skipping brain plots.');
            nSubs_p = 0;
    end
    
    if nSubs_p > 0 && ncomps_pat >= which_comp && nfrex_pat == nfrex
        
        % LFO pattern (group average if multi-subject)
        lfo_pat = squeeze(pats(:, which_comp, idx_lfo, :)); % [nVox x nSubs_p]
        if nSubs_p > 1
            lfo_pat = mean(lfo_pat, 2, 'omitnan');
        end
        
        % Carrier with maximal (mean-over-phase) PAC
        PAC_avg_scalar = mean(PAC_avg, 2, 'omitnan');   % [nCarriers x 1]
        [~, idx_peak_local] = max(PAC_avg_scalar);
        idx_peak_frex = idx_carriers(idx_peak_local);
        
        peak_pat = squeeze(pats(:, which_comp, idx_peak_frex, :));
        if nSubs_p > 1
            peak_pat = mean(peak_pat, 2, 'omitnan');
        end
        
        % Normalize 0–1 for visualization
        lfo_vis  = abs( lfo_pat  ./ max(abs(lfo_pat)) );
        peak_vis = abs( peak_pat ./ max(abs(peak_pat)) );
        
        figure('Name','FREQNESS CrossCoupling: Networks Spatial Patterns','Color','w'); clf
        subplot(1,2,1); hold on
        scatter3(MNI(:,1), MNI(:,2), MNI(:,3), ...
            20, lfo_vis, 'filled');
        title(sprintf('LFO network (%.3g Hz, comp %d)', lfo_freq_actual, which_comp));
        axis equal; grid on; view(135,30); colorbar;
        xlabel('X'); ylabel('Y'); zlabel('Z');
        
        subplot(1,2,2); hold on
        scatter3(MNI(:,1), MNI(:,2), MNI(:,3), ...
            20, peak_vis, 'filled');
        title(sprintf('Peak carrier (%.3g Hz, comp %d)', ...
            frex(idx_peak_frex), which_comp));
        axis equal; grid on; view(135,30); colorbar;
        xlabel('X'); ylabel('Y'); zlabel('Z');
        
    end
end

%% -------------------------- Per-subject plots ---------------------------

if plot_all
    cmap = parula(ncars);

    for subi = 1:nsubs
        
        figure('Name',sprintf('FREQNESS CFC: PAC Sub %d',subi), ...
               'Color','w','Units','normalized','Position',[0.32 0.32 0.35 0.55]); 
        clf
        
        % --- PAC histograms (power over phase) for this subject ---
        subplot(2,1,1); hold on
        for carri = 1:ncars
            plot(phase_centers, ...
                 squeeze(PAC_all(carri,subi,:)), ...
                 'Color', cmap(carri,:), ...
                 'LineWidth', 1.6);
        end
        xlabel('Phase (rad)', 'FontSize', 13, 'FontWeight','bold');
        ylabel('Modulation amplitude (a.u.)', 'FontSize', 13, 'FontWeight','bold');
        title(sprintf('PAC across frequencies (Sub %d): LFO %.3g Hz (comp %d)', ...
            subi, lfo_freq_actual, which_comp), ...
            'FontSize',14,'FontWeight','bold');
        set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
        grid on; grid minor;

        % Legend with carrier frequencies
        leg_entries = cell(ncars,1);
        for carri = 1:ncars
            leg_entries{carri} = sprintf('%.1f Hz', carrier_frex(carri));
        end
        legend(leg_entries,'Location','eastoutside');

        % --- Sine-fit amplitude (power modulation) for this subject ---
        subplot(2,1,2); hold on
        col_line   = cmap(1,:);              % first parula color
        sAmpl_sub  = sAmpl(:,subi);          % [nCarriers x 1]
        plot(carrier_frex, sAmpl_sub, ...
            'Color', col_line, 'LineWidth', 1.8);
        xlabel('Carrier frequency (Hz)', 'FontSize', 13, 'FontWeight','bold');
        ylabel('Power (a.u.)', 'FontSize', 13, 'FontWeight','bold');
        title(sprintf('Power modulation (Sub %d): LFO %.3g Hz (comp %d)', ...
            subi, lfo_freq_actual, which_comp), ...
            'FontSize',14,'FontWeight','bold');
        set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
        grid on; grid minor;

    end
end

end


%% Helper Function: Parse Name-Value Pairs
function opts = parse_name_value_pairs(opts, varargin)
% Simple name–value parser used to handle optional arguments.
%
% Example:
%   opts = struct('mni', [], 'frex2model', [], 'which_comp', 1, 'plot_all', false);
%   opts = parse_name_value_pairs(opts, 'mni', MNIcoords, 'plot_all', true);
%
% All field names are case-insensitive.

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
