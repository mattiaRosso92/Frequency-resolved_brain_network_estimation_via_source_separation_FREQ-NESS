function PEAK = FREQNESS_IndividualPeak(FREQ,frex_range,varargin)

% ========================================================================
%
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you use this toolbox, please cite:
%  Rosso, M., Fernández-Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ-NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved Brain
%  Networks During Auditory Stimulation. Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%  Find each participant's largest stored network-landscape value within a
%  frequency range. The returned frequency is an existing FREQ.frex bin.
%  Eigenvalues are used as supplied, without normalization, smoothing,
%  interpolation, or a requirement for an interior local maximum.
%
%  This is a network-landscape peak, not a power-spectrum peak. The method
%  follows the individual alpha-network peak analysis described in:
%  https://doi.org/10.64898/2026.09.07.749858
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - FREQ : structure containing
%           FREQ.evals -> [nComp x nFrex x (nSubs)]
%           FREQ.frex  -> positive, strictly increasing frequencies in Hz
%
%  - frex_range : mandatory [lower upper] limits in Hz, within FREQ.frex.
%                  Only bins inside this inclusive interval are searched.
%                  At least two stored frequency bins must be included.
%
% ------------------------------------------------------------------------
%  OPTIONAL NAME-VALUE PAIRS:
% ------------------------------------------------------------------------
%  - 'which_comp' : one-based component index or vector. Default: 1.
%  - 'plot_avg'   : plot the distribution of individual peak frequencies,
%                   with participant points, a boxplot, and the group mean.
%                   One figure per selected component. Default: false.
%  - 'plot_all'   : plot every participant's selected-component landscape,
%                   search interval, and selected peak. Default: false.
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%  - PEAK : structure containing
%           peak_frex   -> peak frequencies in Hz [nComp x nSubs]
%           peak_value  -> stored landscape values [nComp x nSubs]
%           idx_peak    -> one-based indices in the full FREQ.frex
%                          [nComp x nSubs]
%           which_comp  -> selected component indices
%           frex_range  -> requested frequency limits in Hz
%           frex        -> actual frequency bins searched
%           at_boundary -> maxima at the first or last searched bin
%                          [nComp x nSubs], logical
%           is_tied     -> multiple bins share the maximum
%                          [nComp x nSubs], logical
%
%  Ties select the lowest-frequency maximum. Non-finite values are ignored;
%  peak_frex, peak_value, and idx_peak are NaN when no valid bin remains.
%  Flags are false for missing results. Boundary and tied maxima are retained
%  for inspection; they do not establish a distinct oscillatory peak.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso, Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 19/10/2026
% ========================================================================

%% Check inputs and optional arguments

if nargin<2
    error('frex_range is required as [lower upper] limits in Hz.');
end
if ~isstruct(FREQ) || ~isscalar(FREQ) || ...
        ~isfield(FREQ,'evals') || ~isfield(FREQ,'frex')
    error('FREQ must be one structure containing evals and frex.');
end
evals = FREQ.evals;
frex = FREQ.frex;
if ~isnumeric(evals) || ~isreal(evals) || isempty(evals) || ndims(evals)>3
    error('FREQ.evals must be real numeric [nComp x nFrex x (nSubs)].');
end
if ~isnumeric(frex) || ~isreal(frex) || ~isvector(frex) || ...
        numel(frex)~=size(evals,2) || any(~isfinite(frex(:))) || ...
        any(frex(:)<=0) || any(diff(frex(:))<=0)
    error(['FREQ.frex must contain positive, strictly increasing frequencies ' ...
        'matching the frequency dimension of FREQ.evals.']);
end
frex = double(frex(:));
if ~isnumeric(frex_range) || ~isreal(frex_range) || ...
        ~isvector(frex_range) || numel(frex_range)~=2 || ...
        any(~isfinite(frex_range(:))) || frex_range(1)>=frex_range(2)
    error('frex_range must contain two finite ascending limits [lower upper].');
end
frex_range = double(frex_range(:)');
if frex_range(1)<frex(1) || frex_range(2)>frex(end)
    error('frex_range must lie within the analyzed FREQ.frex range.');
end
idx_frex = find(frex>=frex_range(1) & frex<=frex_range(2));
if numel(idx_frex)<2
    error('frex_range must include at least two stored frequency bins.');
end

opts = struct('which_comp',1,'plot_avg',false,'plot_all',false);
opts = parse_name_value_pairs(opts,varargin{:});
which_comp = opts.which_comp;
if ~isnumeric(which_comp) || ~isreal(which_comp) || ...
        ~isvector(which_comp) || isempty(which_comp) || ...
        any(~isfinite(which_comp(:))) || ...
        any(which_comp(:)~=round(which_comp(:))) || ...
        any(which_comp(:)<1 | which_comp(:)>size(evals,1))
    error('which_comp must contain valid one-based component indices.');
end
which_comp = unique(double(which_comp(:)'),'stable');
if ~islogical(opts.plot_avg) || ~isscalar(opts.plot_avg) || ...
        ~islogical(opts.plot_all) || ~isscalar(opts.plot_all)
    error('plot_avg and plot_all must be logical scalars.');
end

%% Find individual maxima within the requested range

ncomps = numel(which_comp);
nsubs = size(evals,3);
PEAK = struct;
PEAK.peak_frex = nan(ncomps,nsubs);
PEAK.peak_value = nan(ncomps,nsubs);
PEAK.idx_peak = nan(ncomps,nsubs);
PEAK.which_comp = which_comp;
PEAK.frex_range = frex_range;
PEAK.frex = frex(idx_frex);
PEAK.at_boundary = false(ncomps,nsubs);
PEAK.is_tied = false(ncomps,nsubs);

fprintf(['\nFREQNESS Individual Peak: searching %.4g-%.4g Hz ' ...
    '(%d bins) for %d component(s) and %d participant(s).\n'], ...
    frex_range(1),frex_range(2),numel(idx_frex),ncomps,nsubs);

for compi = 1:ncomps
    for subi = 1:nsubs
        this_landscape = double(reshape( ...
            evals(which_comp(compi),idx_frex,subi),[],1));
        mask_valid = isfinite(this_landscape);
        if ~any(mask_valid)
            continue
        end
        this_landscape(~mask_valid) = -Inf;
        [peak_value,idx_local] = max(this_landscape);
        idx_peak = idx_frex(idx_local);
        PEAK.peak_frex(compi,subi) = frex(idx_peak);
        PEAK.peak_value(compi,subi) = peak_value;
        PEAK.idx_peak(compi,subi) = idx_peak;
        PEAK.at_boundary(compi,subi) = ...
            idx_local==1 || idx_local==numel(idx_frex);
        PEAK.is_tied(compi,subi) = sum(this_landscape==peak_value)>1;
    end
end

%% Visualize group distributions and individual landscapes

if opts.plot_avg || opts.plot_all
    plot_peak_results(evals,frex,PEAK,opts.plot_avg,opts.plot_all);
end

end


%% Helper Function: Plot Peak Results
function plot_peak_results(evals,frex,PEAK,plot_avg,plot_all)

col_line = [0 .45 .70];
col_peak = [.84 .37 0];
nsubs = size(evals,3);
for compi = 1:numel(PEAK.which_comp)
    this_comp = PEAK.which_comp(compi);
    if plot_avg
        figure('Color','w','Position',[100 100 480 520]);
        hold on
        this_peaks = PEAK.peak_frex(compi,:);
        this_peaks = this_peaks(isfinite(this_peaks));
        if ~isempty(this_peaks)
            boxchart(ones(size(this_peaks)),this_peaks, ...
                'BoxFaceColor',col_line,'BoxFaceAlpha',.15, ...
                'MarkerStyle','none','BoxWidth',.35);
            offsets = .12*sin((1:numel(this_peaks))*2.399963);
            scatter(1+offsets,this_peaks,40,col_line,'filled', ...
                'MarkerFaceAlpha',.75,'MarkerEdgeColor','w');
            scatter(1,mean(this_peaks),90,col_peak,'d','filled', ...
                'MarkerEdgeColor','w');
        else
            text(1,mean(PEAK.frex_range),'No valid peaks', ...
                'HorizontalAlignment','center');
        end
        xlim([.6 1.4])
        padding = .08*diff(PEAK.frex_range);
        ylim(PEAK.frex_range+[-padding padding])
        set(gca,'XTick',1,'XTickLabel',{'Participants'},'Box','off', ...
            'FontSize',11,'LineWidth',1,'YGrid','on','GridAlpha',.15)
        ylabel('Individual peak frequency (Hz)','FontWeight','bold')
        title(sprintf('Individual Peak - Component %d',this_comp))
        subtitle(sprintf('n = %d | Diamond: group mean',numel(this_peaks)))
    end
    if plot_all
        for subi = 1:nsubs
            figure('Color','w','Position',[100 100 850 480]);
            hold on
            this_landscape = double(reshape(evals(this_comp,:,subi),[],1));
            this_landscape(~isfinite(this_landscape)) = NaN;
            plot(frex,this_landscape,'Color',col_line,'LineWidth',2);
            this_limits = ylim;
            padding = .06*diff(this_limits);
            this_limits = this_limits+[-padding padding];
            ylim(this_limits)
            patch(PEAK.frex_range([1 2 2 1]), ...
                this_limits([1 1 2 2]),col_line,'FaceAlpha',.07, ...
                'EdgeColor','none','HandleVisibility','off');
            xline(PEAK.frex_range(1),'--','Color',[.5 .5 .5]);
            xline(PEAK.frex_range(2),'--','Color',[.5 .5 .5]);
            if isfinite(PEAK.peak_frex(compi,subi))
                scatter(PEAK.peak_frex(compi,subi), ...
                    PEAK.peak_value(compi,subi),70,col_peak,'filled', ...
                    'MarkerEdgeColor','w');
                subtitle(sprintf('Selected maximum: %.4g Hz', ...
                    PEAK.peak_frex(compi,subi)))
            else
                subtitle('No valid values in the search interval')
            end
            xlim([frex(1) frex(end)])
            xlabel('Frequency (Hz)','FontWeight','bold')
            ylabel('Stored landscape value','FontWeight','bold')
            title(sprintf('Individual Peak - Participant %d - Component %d', ...
                subi,this_comp))
            set(gca,'Box','off','FontSize',11,'LineWidth',1, ...
                'YGrid','on','GridAlpha',.15)
        end
    end
end

end


%% Helper Function: Parse Name-Value Pairs
function opts = parse_name_value_pairs(opts,varargin)

if mod(numel(varargin),2)~=0
    error('Arguments must be given as name-value pairs.');
end
for argi = 1:2:numel(varargin)
    name = varargin{argi};
    if ~(ischar(name) && isrow(name)) && ~(isstring(name) && isscalar(name))
        error('Optional argument names must be text scalars.');
    end
    name = lower(char(name));
    if ~isfield(opts,name)
        error('Unrecognized argument: %s',name);
    end
    opts.(name) = varargin{argi+1};
end

end
