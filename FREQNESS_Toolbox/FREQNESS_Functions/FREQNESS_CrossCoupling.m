function CFC = FREQNESS_CrossCoupling(FREQ,lfo_freq,varargin)

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
%  This function computes phase-amplitude cross-frequency coupling (PAC)
%  between a low-frequency oscillator (LFO) and higher-frequency carrier
%  networks using the FREQNESS component time series stored in FREQ.ts.
%
%  The selected component is re-filtered at the LFO and carrier frequencies.
%  The neural LFO phase and carrier power are obtained from their analytic
%  signals, and mean carrier power is computed in LFO phase bins.
%
%  PAC modulation is summarized through a deterministic first-harmonic
%  regression fitted over the phase-bin centres:
%
%      P(phi) = b0 + bc*cos(phi) + bs*sin(phi)
%
%  The function operates on one component only, which acts both as the LFO
%  modulator and as the carrier component across frequencies.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - FREQ : structure with fields
%           • FREQ.frex  -> frequency vector [nFrex x 1]
%           • FREQ.ts    -> component time series
%                          [nComp x nTime x nFrex x (nSubs)]
%           • FREQ.srate -> sampling rate (Hz)
%           • FREQ.fwhm  -> one filter width per frequency
%           • FREQ.pats  -> optional spatial activation patterns
%
%  - lfo_freq : positive scalar LFO frequency (Hz). The closest value in
%               FREQ.frex is used.
%
% ------------------------------------------------------------------------
%  OPTIONAL NAME-VALUE PAIRS:
% ------------------------------------------------------------------------
%  - 'mni'        : MNI coordinates [nVoxels x 3], used only for optional
%                   spatial visualization.
%
%  - 'frex2model' : ascending [fmin fmax] carrier-frequency range.
%                   Default: all frequencies above the selected LFO.
%
%  - 'which_comp' : component used as modulator and carrier. Default: 1.
%
%  - 'plot_all'   : plot participant-level PAC summaries. Default: false.
%
%  - 'nbins'      : number of LFO phase bins. Default: 37.
%
%  - 'min_valid_bins' : minimum populated phase bins required for harmonic
%                       regression. Default: max(6,ceil(nbins/2)).
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%  - CFC.PAC_all       : [nCarriers x nSubs x nBins] PAC histograms
%  - CFC.PAC_avg       : [nCarriers x nBins] group-average PAC
%  - CFC.fitted_PAC    : first-harmonic fits over phase bins
%  - CFC.coefficients  : [b0 bc bs] regression coefficients
%  - CFC.amplitude_raw : first-harmonic amplitude
%  - CFC.amplitude_normalized : amplitude as % of absolute DC offset
%  - CFC.preferred_phase : preferred LFO phase (rad)
%  - CFC.dc_offset     : first-harmonic DC offset
%  - CFC.mse / CFC.r2  : regression goodness-of-fit
%  - CFC.valid_bins    : populated bins available to each fit
%  - CFC.lfo_phase     : neural LFO phase time series
%
%  Legacy aliases are retained: sAmpl, pShift, dcOff, goodFit, and mFrex.
%  mFrex is fixed to 1/nbins cycles per bin, equivalent to one cycle per
%  LFO phase cycle.
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


%% Check mandatory inputs

if nargin < 2 || ~isstruct(FREQ)
    error('FREQ and lfo_freq are required inputs.');
end

required_fields = {'frex','ts','srate','fwhm'};
for fieldi = 1:numel(required_fields)
    if ~isfield(FREQ,required_fields{fieldi}) || isempty(FREQ.(required_fields{fieldi}))
        error('FREQ.%s is missing or empty.',required_fields{fieldi});
    end
end

frex = FREQ.frex(:);
if ~isnumeric(frex) || ~isreal(frex) || any(~isfinite(frex)) || ...
        any(frex <= 0) || any(diff(frex) <= 0)
    error('FREQ.frex must contain positive, finite, strictly increasing frequencies.');
end

if ~isnumeric(lfo_freq) || ~isreal(lfo_freq) || ~isscalar(lfo_freq) || ...
        ~isfinite(lfo_freq) || lfo_freq <= 0
    error('lfo_freq must be one positive finite scalar.');
end

if ~isnumeric(FREQ.srate) || ~isreal(FREQ.srate) || ...
        ~isscalar(FREQ.srate) || ~isfinite(FREQ.srate) || FREQ.srate <= 0
    error('FREQ.srate must be one positive finite scalar.');
end

if any(frex >= FREQ.srate/2)
    error('FREQ.frex must be lower than the Nyquist frequency.');
end

ts = FREQ.ts;
if ~isnumeric(ts) || ~isreal(ts) || any(~isfinite(ts(:))) || ...
        ndims(ts) < 3 || ndims(ts) > 4
    error(['FREQ.ts must be a real finite numeric array in ' ...
        '[nComp x nTime x nFrex x (nSubs)] format.']);
end

ncomps = size(ts,1);
npnts = size(ts,2);
nfrex = size(ts,3);
nsubs = size(ts,4);

if npnts < 2
    error('FREQ.ts must contain at least two timepoints.');
end

if nfrex ~= numel(frex)
    error('FREQ.ts and FREQ.frex have incompatible frequency dimensions.');
end

fwhm = FREQ.fwhm(:);
if ~isnumeric(fwhm) || ~isreal(fwhm) || numel(fwhm) ~= nfrex || ...
        any(~isfinite(fwhm)) || any(fwhm <= 0)
    error('FREQ.fwhm must contain one positive finite value per frequency.');
end


%% Handle optional arguments

opts = struct('mni',[],'frex2model',[],'which_comp',1, ...
    'plot_all',false,'nbins',37,'min_valid_bins',[]);
opts = parse_name_value_pairs(opts,varargin{:});

MNI = opts.mni;
frex2model = opts.frex2model;
which_comp = opts.which_comp;
plot_all = opts.plot_all;
nbins = opts.nbins;
min_valid_bins = opts.min_valid_bins;

if isempty(frex2model)
    frex2model = [frex(1) frex(end)];
elseif ~isnumeric(frex2model) || ~isreal(frex2model) || ...
        ~isvector(frex2model) || numel(frex2model) ~= 2 || ...
        any(~isfinite(frex2model)) || frex2model(1) > frex2model(2)
    error('frex2model must contain two finite ascending frequency boundaries.');
else
    frex2model = frex2model(:)';
end

if ~isnumeric(which_comp) || ~isreal(which_comp) || ~isscalar(which_comp) || ...
        ~isfinite(which_comp) || which_comp ~= round(which_comp) || ...
        which_comp < 1 || which_comp > ncomps
    error('which_comp must be an integer between 1 and %d.',ncomps);
end

if ~islogical(plot_all) || ~isscalar(plot_all)
    error('plot_all must be one logical value (true or false).');
end

if ~isnumeric(nbins) || ~isreal(nbins) || ~isscalar(nbins) || ...
        ~isfinite(nbins) || nbins ~= round(nbins) || nbins < 6
    error('nbins must be one integer greater than or equal to 6.');
end

if isempty(min_valid_bins)
    min_valid_bins = max(6,ceil(nbins/2));
elseif ~isnumeric(min_valid_bins) || ~isreal(min_valid_bins) || ...
        ~isscalar(min_valid_bins) || ~isfinite(min_valid_bins) || ...
        min_valid_bins ~= round(min_valid_bins) || min_valid_bins < 3 || ...
        min_valid_bins > nbins
    error('min_valid_bins must be one integer from 3 to nbins.');
end

if ~isempty(MNI) && (~isnumeric(MNI) || ~isreal(MNI) || ...
        ndims(MNI) ~= 2 || size(MNI,2) ~= 3 || any(isinf(MNI(:))))
    error('MNI must be a real [nVoxels x 3] numeric matrix without infinite values.');
end


%% Define LFO and carrier frequencies

[~,idx_lfo] = min(abs(frex-lfo_freq));
lfo_freq_actual = frex(idx_lfo);

if abs(lfo_freq_actual-lfo_freq) > 1e-3
    warning('lfo_freq %.3g Hz is unavailable. Using %.3g Hz.', ...
        lfo_freq,lfo_freq_actual);
end

mask_frex = frex >= frex2model(1) & frex <= frex2model(2);
mask_frex(1:idx_lfo) = false;
idx_carriers = find(mask_frex);

if isempty(idx_carriers)
    error('No carrier frequencies remain above the selected LFO.');
end

carrier_frex = frex(idx_carriers);
ncars = numel(carrier_frex);


%% Compute neural LFO phase

lfo_phase_all = zeros(npnts,nsubs);

for subi = 1:nsubs
    lfo_ts = squeeze(ts(which_comp,:,idx_lfo,subi));
    lfo_ts = filterFGx(lfo_ts,FREQ.srate,lfo_freq_actual,fwhm(idx_lfo),0);
    lfo_phase_all(:,subi) = angle(FREQNESS_AnalyticSignal(lfo_ts'));
end

phase_edges = linspace(-pi,pi,nbins+1);
phase_centers = phase_edges(1:end-1) + diff(phase_edges)/2;


%% Compute PAC across carriers and participants

PAC_all = nan(ncars,nsubs,nbins);

for carri = 1:ncars
    idx_freq = idx_carriers(carri);

    for subi = 1:nsubs
        carr_ts = squeeze(ts(which_comp,:,idx_freq,subi));
        carr_ts = filterFGx(carr_ts,FREQ.srate,frex(idx_freq),fwhm(idx_freq),0);
        carr_pow = abs(FREQNESS_AnalyticSignal(carr_ts')).^2;
        phase = lfo_phase_all(:,subi);

        for bini = 1:nbins
            idx_phase = phase > phase_edges(bini) & phase <= phase_edges(bini+1);
            if any(idx_phase)
                PAC_all(carri,subi,bini) = mean(carr_pow(idx_phase),'omitnan');
            end
        end
    end
end

PAC_avg = reshape(mean(PAC_all,2,'omitnan'),ncars,nbins);


%% Fit deterministic first-harmonic regression

coefficients = nan(ncars,nsubs,3);
fitted_PAC = nan(size(PAC_all));
amplitude_raw = nan(ncars,nsubs);
amplitude_normalized = nan(ncars,nsubs);
preferred_phase = nan(ncars,nsubs);
dc_offset = nan(ncars,nsubs);
mse = nan(ncars,nsubs);
r2 = nan(ncars,nsubs);
valid_bins = zeros(ncars,nsubs);

for carri = 1:ncars
    for subi = 1:nsubs
        y = squeeze(PAC_all(carri,subi,:));
        valid_fit = isfinite(phase_centers(:)) & isfinite(y);
        valid_bins(carri,subi) = sum(valid_fit);

        if valid_bins(carri,subi) < min_valid_bins
            continue
        end

        phi = phase_centers(valid_fit)';
        y_valid = y(valid_fit);
        design = [ones(numel(phi),1) cos(phi) sin(phi)];
        beta = design\y_valid;
        yfit = beta(1) + beta(2)*cos(phase_centers(:)) + ...
            beta(3)*sin(phase_centers(:));
        residuals = y_valid-design*beta;
        sse = sum(residuals.^2);
        sst = sum((y_valid-mean(y_valid)).^2);

        coefficients(carri,subi,:) = reshape(beta,1,1,3);
        fitted_PAC(carri,subi,:) = reshape(yfit,1,1,nbins);
        amplitude_raw(carri,subi) = hypot(beta(2),beta(3));
        preferred_phase(carri,subi) = atan2(beta(3),beta(2));
        dc_offset(carri,subi) = beta(1);
        mse(carri,subi) = mean(residuals.^2);

        if abs(beta(1)) > eps
            amplitude_normalized(carri,subi) = ...
                100*amplitude_raw(carri,subi)/abs(beta(1));
        end
        if sst > eps
            r2(carri,subi) = 1-sse/sst;
        end
    end
end

mFrex = nan(ncars,nsubs);
mFrex(isfinite(amplitude_raw)) = 1/nbins;
amplitude_avg = mean(amplitude_raw,2,'omitnan');


%% Pack output

CFC.PAC_all = PAC_all;
CFC.PAC_avg = PAC_avg;
CFC.carrier_frex = carrier_frex;
CFC.lfo_freq = lfo_freq_actual;
CFC.comp = which_comp;
CFC.phase_edges = phase_edges;
CFC.phase_centers = phase_centers;
CFC.fitted_PAC = fitted_PAC;
CFC.coefficients = coefficients;
CFC.amplitude_raw = amplitude_raw;
CFC.amplitude_normalized = amplitude_normalized;
CFC.preferred_phase = preferred_phase;
CFC.dc_offset = dc_offset;
CFC.mse = mse;
CFC.r2 = r2;
CFC.valid_bins = valid_bins;
CFC.lfo_phase = lfo_phase_all;

% Legacy output aliases
CFC.sAmpl = amplitude_raw;
CFC.pShift = preferred_phase;
CFC.dcOff = dc_offset;
CFC.mFrex = mFrex;
CFC.goodFit = mse;


%% Group-level PAC visualization

figure('Name','FREQNESS CrossCoupling: PAC','Color','w');
cmap = parula(ncars);

subplot(2,1,1); hold on
for carri = 1:ncars
    plot(phase_centers,PAC_avg(carri,:),'Color',cmap(carri,:), ...
        'LineWidth',1.6);
end
xlabel('LFO phase (rad)','FontSize',13,'FontWeight','bold');
ylabel('Carrier power (a.u.)','FontSize',13,'FontWeight','bold');
title(sprintf(['PAC across frequencies - Group average: LFO %.3g Hz ' ...
    '(component %d)'],lfo_freq_actual,which_comp), ...
    'FontSize',14,'FontWeight','bold');
set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
grid on
grid minor

leg_entries = arrayfun(@(x) sprintf('%.1f Hz',x),carrier_frex, ...
    'UniformOutput',false);
legend(leg_entries,'Location','eastoutside');

subplot(2,1,2); hold on
plot(carrier_frex,amplitude_avg,'Color',cmap(1,:),'LineWidth',1.8);
xlabel('Carrier frequency (Hz)','FontSize',13,'FontWeight','bold');
ylabel('First-harmonic amplitude (a.u.)','FontSize',13,'FontWeight','bold');
title(sprintf('Power modulation: LFO %.3g Hz (component %d)', ...
    lfo_freq_actual,which_comp),'FontSize',14,'FontWeight','bold');
set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
grid on
grid minor


%% Optional spatial visualization

if ~isempty(MNI)
    plot_spatial_patterns(FREQ,MNI,which_comp,idx_lfo,idx_carriers, ...
        PAC_avg,frex,nsubs)
end


%% Optional participant-level PAC visualizations

if plot_all
    for subi = 1:nsubs
        figure('Name',sprintf('FREQNESS CFC: PAC Participant %d',subi), ...
            'Color','w','Units','normalized','Position',[0.32 0.32 0.35 0.55]);

        subplot(2,1,1); hold on
        for carri = 1:ncars
            plot(phase_centers,squeeze(PAC_all(carri,subi,:)), ...
                'Color',cmap(carri,:),'LineWidth',1.6);
        end
        xlabel('LFO phase (rad)','FontSize',13,'FontWeight','bold');
        ylabel('Carrier power (a.u.)','FontSize',13,'FontWeight','bold');
        title(sprintf(['PAC - Participant #%d: LFO %.3g Hz ' ...
            '(component %d)'],subi,lfo_freq_actual,which_comp), ...
            'FontSize',14,'FontWeight','bold');
        set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
        grid on
        grid minor
        legend(leg_entries,'Location','eastoutside');

        subplot(2,1,2); hold on
        plot(carrier_frex,amplitude_raw(:,subi), ...
            'Color',cmap(1,:),'LineWidth',1.8);
        xlabel('Carrier frequency (Hz)','FontSize',13,'FontWeight','bold');
        ylabel('First-harmonic amplitude (a.u.)', ...
            'FontSize',13,'FontWeight','bold');
        title(sprintf('Power modulation - Participant #%d',subi), ...
            'FontSize',14,'FontWeight','bold');
        set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');
        grid on
        grid minor
    end
end

end


%% Helper function: Plot spatial patterns
function plot_spatial_patterns(FREQ,MNI,which_comp,idx_lfo,idx_carriers,PAC_avg,frex,nsubs)

if ~isfield(FREQ,'pats') || isempty(FREQ.pats)
    warning('MNI coordinates were provided but FREQ.pats is unavailable. Skipping spatial plots.');
    return
end

pats = FREQ.pats;
if ~isnumeric(pats) || ~isreal(pats) || any(isinf(pats(:))) || ...
        ndims(pats) < 3 || ndims(pats) > 4 || ...
        size(pats,2) < which_comp || size(pats,3) ~= numel(frex) || ...
        size(pats,4) ~= nsubs
    warning('FREQ.pats has incompatible dimensions. Skipping spatial plots.');
    return
end

nvoxs = size(pats,1);
if size(MNI,1) ~= nvoxs
    error('MNI must contain one row per voxel in FREQ.pats.');
end

valid_mni = all(isfinite(MNI),2);
if ~any(valid_mni)
    error('MNI contains no valid coordinates.');
end

lfo_pat = reshape(pats(:,which_comp,idx_lfo,:),nvoxs,nsubs);
lfo_pat = mean(lfo_pat,2,'omitnan');

carrier_power = mean(PAC_avg,2,'omitnan');
if ~any(isfinite(carrier_power))
    warning('PAC estimates are unavailable. Skipping spatial plots.');
    return
end
[~,idx_peak_local] = max(carrier_power);
idx_peak_frex = idx_carriers(idx_peak_local);
peak_pat = reshape(pats(:,which_comp,idx_peak_frex,:),nvoxs,nsubs);
peak_pat = mean(peak_pat,2,'omitnan');

lfo_vis = normalize_pattern(lfo_pat);
peak_vis = normalize_pattern(peak_pat);

figure('Name','FREQNESS CrossCoupling: Spatial Patterns','Color','w');
subplot(1,2,1); hold on
scatter3(MNI(valid_mni,1),MNI(valid_mni,2),MNI(valid_mni,3), ...
    20,lfo_vis(valid_mni),'filled');
title(sprintf('LFO network (%.3g Hz, component %d)', ...
    frex(idx_lfo),which_comp));
axis equal
grid on
view(135,30)
colorbar
xlabel('MNI X'); ylabel('MNI Y'); zlabel('MNI Z');

subplot(1,2,2); hold on
scatter3(MNI(valid_mni,1),MNI(valid_mni,2),MNI(valid_mni,3), ...
    20,peak_vis(valid_mni),'filled');
title(sprintf('Peak carrier (%.3g Hz, component %d)', ...
    frex(idx_peak_frex),which_comp));
axis equal
grid on
view(135,30)
colorbar
xlabel('MNI X'); ylabel('MNI Y'); zlabel('MNI Z');

end


%% Helper function: Normalize spatial pattern
function pattern = normalize_pattern(pattern)

pattern = abs(pattern);
valid_pattern = isfinite(pattern);

if any(valid_pattern)
    max_pattern = max(pattern(valid_pattern));
    if max_pattern > 0
        pattern(valid_pattern) = pattern(valid_pattern)/max_pattern;
    else
        pattern(valid_pattern) = 0;
    end
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
