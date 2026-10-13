function STATS = FREQNESS_LandscapeComparison(LANDSCAPE_A,LANDSCAPE_B,varargin)

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
%  This function compares participant-level FREQ-NESS landscapes between
%  two samples. Paired within-subject and independent between-subject
%  designs are supported. Statistical tests are performed separately at
%  every requested frequency, with Benjamini-Hochberg FDR correction across
%  frequencies for each measure and component by default.
%
%  Recognized landscape fields are FREQ.evals, ENTR.H2, ENTR.ED,
%  CFC.amplitude_raw, CFC.amplitude_normalized, STAB.stability_index, and
%  STAB.mean_frequency. When 'measure' is omitted, every recognized field
%  shared by both inputs is analyzed separately and retains its original
%  field name in STATS. Raw numeric landscapes are also accepted when their
%  frequency axis is supplied through 'frex'.
% ========================================================================

% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - LANDSCAPE_A, LANDSCAPE_B : structures containing participant-level
%            landscapes, or numeric arrays. Supported numeric layouts are
%            [nFrex x nParticipants] and
%            [nFeatures x nFrex x nParticipants].
%
%  - measure : field name or cell array of field names. When omitted, every
%              recognized shared landscape is analyzed.
%
%  - design : 'within' (default) or 'between'.
%
%  - test : 'ttest' or 'signrank' for within-subject designs; 'ttest' or
%           'ranksum' for between-subject designs. Defaults are 'signrank'
%           for within-subject and 'ranksum' for between-subject designs.
%
%  - which_comp : one-based feature/component numbers. Default: 1.
%
%  - frex : frequency axis for raw numeric inputs. Structure inputs use
%           their stored frex or carrier_frex field.
%
%  - frex2compare : frequencies to compare. Values are mapped to the closest
%                   stored frequency. Default: every frequency.
%
%  - sample_labels : labels for the two conditions or groups.
%                    Default: {'A','B'}.
%
%  - tail : 'both' (default), 'left', or 'right'.
%
%  - alpha : significance threshold. Default: .05.
%
%  - correction : 'fdr' (default) or 'none'.
%
%  - frex_tolerance : maximum accepted frequency-matching distance in Hz.
%                     Default: half the median frequency spacing.
%
%  - plot_results : logical flag controlling result figures. Default: true.
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - STATS : structure containing analysis metadata and one field per
%            analyzed landscape. For example, STATS.evals contains selected
%            observations, descriptive statistics, test statistics, effect
%            sizes, uncorrected p values, FDR q values, and significance.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso, Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 13/10/2026
%
% ========================================================================

%% Map inputs

opts = struct('measure', [], 'design', 'within', 'test', [], ...
    'which_comp', 1, 'frex', [], 'frex2compare', [], ...
    'sample_labels', {{'A','B'}}, 'tail', 'both', 'alpha', .05, ...
    'correction', 'fdr', 'frex_tolerance', [], 'plot_results', true);
opts = parse_name_value_pairs(opts,varargin{:});

[design,test,tail] = validate_analysis_options(opts.design,opts.test,opts.tail);
sample_labels = validate_sample_labels(opts.sample_labels);
validate_probability(opts.alpha,'alpha');
correction = validatestring(opts.correction,{'fdr','none'});
if ~islogical(opts.plot_results) || ~isscalar(opts.plot_results)
    error('plot_results must be one logical value.');
end

measures = select_measures(LANDSCAPE_A,LANDSCAPE_B,opts.measure);

STATS = struct;
STATS.analysis = 'FREQ-NESS landscape comparison';
STATS.design = design;
STATS.test = test;
STATS.sample_labels = sample_labels;
STATS.alpha = opts.alpha;
STATS.correction = correction;
STATS.correction_scope = ...
    'frequencies separately for each measure and component';
STATS.measures = measures;

fprintf(['\nFREQNESS Landscape Comparison: comparing %s and %s using a ' ...
    '%s-subject %s analysis.\n'],sample_labels{1},sample_labels{2}, ...
    design,test);

%% Compare landscape measures

for measurei = 1:numel(measures)
    measure = measures{measurei};
    frex_a = frequency_axis(LANDSCAPE_A,measure,opts.frex);
    frex_b = frequency_axis(LANDSCAPE_B,measure,opts.frex);
    if numel(frex_a)~=numel(frex_b) || any(abs(frex_a-frex_b)>1e-10)
        error('Frequency grids do not match for measure ''%s''.',measure);
    end

    values_a = canonical_values(LANDSCAPE_A,measure,frex_a,'LANDSCAPE_A');
    values_b = canonical_values(LANDSCAPE_B,measure,frex_b,'LANDSCAPE_B');
    if size(values_a,1)~=size(values_b,1)
        error('Feature counts do not match for measure ''%s''.',measure);
    end
    if strcmp(design,'within') && size(values_a,3)~=size(values_b,3)
        error('Within-subject comparisons require matching participant counts.');
    end
    if size(values_a,3)<2 || size(values_b,3)<2
        error('Each sample must contain at least two participants.');
    end

    comps = select_components(opts.which_comp,size(values_a,1));
    [requested_frex,frex2compare,idx_frex2compare] = ...
        select_frequencies(frex_a,opts.frex2compare,opts.frex_tolerance);
    fprintf('  %s: %d frequencies, component(s) %s.\n',measure, ...
        numel(frex2compare),mat2str(comps));

    result = compare_measure(measure,values_a,values_b,requested_frex, ...
        frex2compare,idx_frex2compare,comps,design,test,tail, ...
        opts.alpha,correction);
    if opts.plot_results
        plot_landscape_results(result,sample_labels,correction,opts.alpha);
    end
    STATS.(measure) = result;
end

end

%% Helper Function: Compare One Measure
function result = compare_measure(measure,values_a,values_b,requested_frex, ...
    frex2compare,idx_frex2compare,comps,design,test,tail,alpha,correction)

selected_a = values_a(comps,idx_frex2compare,:);
selected_b = values_b(comps,idx_frex2compare,:);
ncomps = numel(comps);
nfrex = numel(frex2compare);

mean_a = nan(ncomps,nfrex);
mean_b = nan(ncomps,nfrex);
median_a = nan(ncomps,nfrex);
median_b = nan(ncomps,nfrex);
std_a = nan(ncomps,nfrex);
std_b = nan(ncomps,nfrex);
sem_a = nan(ncomps,nfrex);
sem_b = nan(ncomps,nfrex);
mean_difference = nan(ncomps,nfrex);
median_difference = nan(ncomps,nfrex);
statistic = nan(ncomps,nfrex);
effect_size = nan(ncomps,nfrex);
p_uncorrected = nan(ncomps,nfrex);
n_a = zeros(ncomps,nfrex);
n_b = zeros(ncomps,nfrex);
statistic_name = '';
effect_size_name = '';

for compi = 1:ncomps
    for frexi = 1:nfrex
        sample_a = squeeze(selected_a(compi,frexi,:));
        sample_b = squeeze(selected_b(compi,frexi,:));
        if strcmp(design,'within')
            valid = isfinite(sample_a) & isfinite(sample_b);
            sample_a = sample_a(valid);
            sample_b = sample_b(valid);
        else
            sample_a = sample_a(isfinite(sample_a));
            sample_b = sample_b(isfinite(sample_b));
        end
        n_a(compi,frexi) = numel(sample_a);
        n_b(compi,frexi) = numel(sample_b);
        if ~isempty(sample_a)
            mean_a(compi,frexi) = mean(sample_a);
            median_a(compi,frexi) = median(sample_a);
            if numel(sample_a)>1
                std_a(compi,frexi) = std(sample_a);
                sem_a(compi,frexi) = std_a(compi,frexi)/sqrt(numel(sample_a));
            end
        end
        if ~isempty(sample_b)
            mean_b(compi,frexi) = mean(sample_b);
            median_b(compi,frexi) = median(sample_b);
            if numel(sample_b)>1
                std_b(compi,frexi) = std(sample_b);
                sem_b(compi,frexi) = std_b(compi,frexi)/sqrt(numel(sample_b));
            end
        end
        if ~isempty(sample_a) && ~isempty(sample_b)
            mean_difference(compi,frexi) = mean(sample_a)-mean(sample_b);
            median_difference(compi,frexi) = median(sample_a)-median(sample_b);
        end
        if numel(sample_a)<2 || numel(sample_b)<2
            continue
        end
        if strcmp(design,'within')
            [statistic(compi,frexi),statistic_name, ...
                effect_size(compi,frexi),effect_size_name, ...
                p_uncorrected(compi,frexi)] = ...
                paired_test(sample_a,sample_b,test,tail);
        else
            [statistic(compi,frexi),statistic_name, ...
                effect_size(compi,frexi),effect_size_name, ...
                p_uncorrected(compi,frexi)] = ...
                independent_test(sample_a,sample_b,test,tail);
        end
    end
end

q_fdr = nan(size(p_uncorrected));
fdr_threshold = nan(ncomps,1);
if strcmp(correction,'fdr')
    for compi = 1:ncomps
        [q_fdr(compi,:),fdr_threshold(compi)] = ...
            bh_fdr(p_uncorrected(compi,:),alpha);
    end
else
    q_fdr = p_uncorrected;
    fdr_threshold(:) = alpha;
end

result = struct;
result.measure = measure;
result.requested_frex = requested_frex;
result.frex = frex2compare;
result.idx_frex = idx_frex2compare;
result.comps = comps;
result.values_A = selected_a;
result.values_B = selected_b;
result.mean_A = mean_a;
result.mean_B = mean_b;
result.median_A = median_a;
result.median_B = median_b;
result.std_A = std_a;
result.std_B = std_b;
result.sem_A = sem_a;
result.sem_B = sem_b;
result.mean_difference = mean_difference;
result.median_difference = median_difference;
result.statistic = statistic;
result.statistic_name = statistic_name;
result.effect_size = effect_size;
result.effect_size_name = effect_size_name;
result.p_uncorrected = p_uncorrected;
result.q_fdr = q_fdr;
result.significant = isfinite(q_fdr) & q_fdr<=alpha;
result.fdr_threshold = fdr_threshold;
result.n_A = n_a;
result.n_B = n_b;

end

%% Helper Function: Statistical Tests
function [statistic,statistic_name,effect,effect_name,p_value] = ...
    paired_test(sample_a,sample_b,test,tail)

differences = sample_a-sample_b;
if all(differences==0)
    statistic = 0;
    effect = 0;
    p_value = 1;
    if strcmp(test,'ttest')
        statistic_name = 't';
        effect_name = 'Cohen''s dz';
    else
        statistic_name = 'W+';
        effect_name = 'matched rank-biserial r';
    end
    return
end

if strcmp(test,'ttest')
    [~,p_value,~,stats] = ttest(sample_a,sample_b,'Tail',tail);
    statistic = stats.tstat;
    deviation = std(differences);
    if deviation>0
        effect = mean(differences)/deviation;
    else
        effect = sign(mean(differences))*Inf;
    end
    statistic_name = 't';
    effect_name = 'Cohen''s dz';
else
    p_value = signrank(sample_a,sample_b,'Tail',tail,'Method','approximate');
    differences = differences(differences~=0);
    ranks = tiedrank(abs(differences));
    w_positive = sum(ranks(differences>0));
    w_negative = sum(ranks(differences<0));
    statistic = w_positive;
    effect = (w_positive-w_negative)/(w_positive+w_negative);
    statistic_name = 'W+';
    effect_name = 'matched rank-biserial r';
end

end

function [statistic,statistic_name,effect,effect_name,p_value] = ...
    independent_test(sample_a,sample_b,test,tail)

if all(sample_a==sample_a(1)) && all(sample_b==sample_a(1))
    p_value = 1;
    effect = 0;
    if strcmp(test,'ttest')
        statistic = 0;
        statistic_name = 't';
        effect_name = 'Hedges'' g';
    else
        statistic = numel(sample_a)*numel(sample_b)/2;
        statistic_name = 'U';
        effect_name = 'rank-biserial r';
    end
    return
end

if strcmp(test,'ttest')
    [~,p_value,~,stats] = ttest2(sample_a,sample_b,'Tail',tail, ...
        'Vartype','unequal');
    statistic = stats.tstat;
    n_a = numel(sample_a);
    n_b = numel(sample_b);
    pooled = sqrt(((n_a-1)*var(sample_a)+(n_b-1)*var(sample_b))/ ...
        (n_a+n_b-2));
    if pooled>0
        cohens_d = (mean(sample_a)-mean(sample_b))/pooled;
        correction = 1-3/(4*(n_a+n_b)-9);
        effect = correction*cohens_d;
    else
        effect = sign(mean(sample_a)-mean(sample_b))*Inf;
    end
    statistic_name = 't';
    effect_name = 'Hedges'' g';
else
    p_value = ranksum(sample_a,sample_b,'Tail',tail,'Method','approximate');
    ranks = tiedrank([sample_a(:);sample_b(:)]);
    n_a = numel(sample_a);
    n_b = numel(sample_b);
    statistic = sum(ranks(1:n_a))-n_a*(n_a+1)/2;
    effect = 2*statistic/(n_a*n_b)-1;
    statistic_name = 'U';
    effect_name = 'rank-biserial r';
end

end

%% Helper Function: Benjamini-Hochberg FDR
function [q_values,threshold] = bh_fdr(p_values,alpha)

q_values = nan(size(p_values));
valid = isfinite(p_values);
values = p_values(valid);
if isempty(values)
    threshold = NaN;
    return
end

[sorted_p,order] = sort(values);
ranks = (1:numel(sorted_p))';
sorted_q = sorted_p(:).*numel(sorted_p)./ranks;
sorted_q = flipud(cummin(flipud(sorted_q)));
sorted_q = min(sorted_q,1);
restored = zeros(size(sorted_q));
restored(order) = sorted_q;
q_values(valid) = restored;
accepted = sorted_p(:)<=alpha*ranks/numel(sorted_p);
if any(accepted)
    threshold = sorted_p(find(accepted,1,'last'));
else
    threshold = NaN;
end

end

%% Helper Function: Input Mapping
function measures = select_measures(landscape_a,landscape_b,measure)

numeric_a = isnumeric(landscape_a);
numeric_b = isnumeric(landscape_b);
if numeric_a~=numeric_b
    error('LANDSCAPE_A and LANDSCAPE_B must use the same input format.');
end
if numeric_a
    if isempty(measure)
        measures = {'values'};
    else
        measures = normalize_measure_input(measure);
        if numel(measures)~=1 || ~strcmpi(measures{1},'values')
            error('Raw numeric inputs support only measure=''values''.');
        end
        measures = {'values'};
    end
    return
end
if ~isstruct(landscape_a) || ~isstruct(landscape_b)
    error('Landscape inputs must be structures or real numeric arrays.');
end

if isempty(measure)
    recognized = {'evals','H2','ED','amplitude_raw', ...
        'amplitude_normalized','stability_index','mean_frequency'};
    measures = cell(1,0);
    for fieldi = 1:numel(recognized)
        field_a = resolve_field(landscape_a,recognized{fieldi});
        field_b = resolve_field(landscape_b,recognized{fieldi});
        if ~isempty(field_a) && ~isempty(field_b)
            measures{end+1} = field_a; %#ok<AGROW>
        end
    end
    if isempty(measures)
        error(['No shared landscape fields were recognized; supply ' ...
            '''measure'' explicitly.']);
    end
    return
end

requested = normalize_measure_input(measure);
measures = cell(1,0);
for fieldi = 1:numel(requested)
    field_a = resolve_field(landscape_a,requested{fieldi});
    field_b = resolve_field(landscape_b,requested{fieldi});
    if isempty(field_a) || isempty(field_b)
        error('Both inputs must contain the measure ''%s''.',requested{fieldi});
    end
    if ~strcmpi(field_a,field_b)
        error('Measure names do not match for ''%s''.',requested{fieldi});
    end
    if ~any(strcmp(measures,field_a))
        measures{end+1} = field_a; %#ok<AGROW>
    end
end

end

function requested = normalize_measure_input(measure)

if ischar(measure)
    requested = {measure};
elseif isstring(measure)
    requested = cellstr(measure(:)');
elseif iscell(measure) && ~isempty(measure) && ...
        all(cellfun(@(value)ischar(value) && ~isempty(value),measure))
    requested = reshape(measure,1,[]);
else
    error('measure must be text or a non-empty cell array of field names.');
end

end

function actual = resolve_field(container,requested)

names = fieldnames(container);
idx = find(strcmpi(names,requested),1);
if isempty(idx)
    actual = '';
else
    actual = names{idx};
end

end

function values = field_value(container,requested)

actual = resolve_field(container,requested);
if isempty(actual)
    error('Input does not contain the requested field ''%s''.',requested);
end
values = container.(actual);

end

function frex = frequency_axis(container,measure,supplied_frex)

if ~isempty(supplied_frex)
    frex = validate_frequency_vector(supplied_frex,'frex');
    return
end
if isnumeric(container)
    error('frex is required for raw numeric landscape inputs.');
end
if any(strcmpi(measure,{'amplitude_raw','amplitude_normalized'}))
    axis_name = 'carrier_frex';
else
    axis_name = 'frex';
end
frex = validate_frequency_vector(field_value(container,axis_name), ...
    ['LANDSCAPE.' axis_name]);

end

function frex = validate_frequency_vector(frex,name)

if ~isnumeric(frex) || ~isreal(frex) || ~isvector(frex) || ...
        isempty(frex) || any(~isfinite(frex))
    error('%s must be a finite real vector.',name);
end
frex = reshape(double(frex),1,[]);
if any(diff(frex)<=0)
    error('%s must be strictly increasing.',name);
end

end

function values = canonical_values(container,measure,frex,label)

if strcmp(measure,'values') && isnumeric(container)
    values = container;
else
    values = field_value(container,measure);
end
if ~isnumeric(values) || ~isreal(values) || any(isinf(values(:)))
    error('%s.%s must contain finite-or-NaN real numeric values.',label,measure);
end
values = double(values);
if isvector(values) && numel(values)==numel(frex)
    values = reshape(values,[1 numel(frex) 1]);
elseif ismatrix(values) && size(values,1)==numel(frex)
    values = reshape(values,[1 numel(frex) size(values,2)]);
elseif ndims(values)==3 && size(values,2)==numel(frex)
    % Already arranged as [nFeatures x nFrex x nParticipants].
else
    error(['%s.%s must be [nFrex x nParticipants] or ' ...
        '[nFeatures x nFrex x nParticipants].'],label,measure);
end
if any(size(values)==0)
    error('%s.%s is empty.',label,measure);
end

end

%% Helper Function: Selection
function comps = select_components(which_comp,nfeatures)

if isempty(which_comp)
    which_comp = 1;
end
if ~isnumeric(which_comp) || ~isreal(which_comp) || ...
        ~isvector(which_comp) || any(~isfinite(which_comp)) || ...
        any(which_comp~=round(which_comp)) || any(which_comp<1) || ...
        any(which_comp>nfeatures)
    error('which_comp must contain integers between 1 and %d.',nfeatures);
end
comps = reshape(double(which_comp),1,[]);
if numel(unique(comps))~=numel(comps)
    error('which_comp cannot contain duplicates.');
end

end

function [requested,resolved,indices] = ...
    select_frequencies(frex,frex2compare,tolerance)

if isempty(frex2compare) || ...
        ((ischar(frex2compare) || isstring(frex2compare)) && ...
        strcmpi(frex2compare,'all'))
    requested = frex;
    resolved = frex;
    indices = 1:numel(frex);
    return
end
if ~isnumeric(frex2compare) || ~isreal(frex2compare) || ...
        ~isvector(frex2compare) || any(~isfinite(frex2compare))
    error('frex2compare must be a finite numeric vector, ''all'', or empty.');
end
requested = reshape(double(frex2compare),1,[]);
if isempty(tolerance)
    if numel(frex)==1
        tolerance = max(abs(frex(1))*1e-10,1e-10);
    else
        tolerance = median(diff(frex))/2;
    end
end
if ~isnumeric(tolerance) || ~isscalar(tolerance) || ...
        ~isfinite(tolerance) || tolerance<0
    error('frex_tolerance must be one non-negative finite scalar.');
end

indices = zeros(size(requested));
distances = zeros(size(requested));
for frexi = 1:numel(requested)
    [distances(frexi),indices(frexi)] = min(abs(frex-requested(frexi)));
end
if any(distances>tolerance)
    error('One or more requested frequencies exceed frex_tolerance.');
end
if numel(unique(indices))~=numel(indices)
    error('frex2compare maps more than one request to the same bin.');
end
resolved = frex(indices);

end

%% Helper Function: Plot Results
function plot_landscape_results(result,sample_labels,correction,alpha)

colors = parula(4);
colors = colors([1 4],:);
for compi = 1:numel(result.comps)
    figure('Color','w','Position',[100 100 900 700]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile, hold on
    plot_with_error(result.frex,result.mean_A(compi,:), ...
        result.sem_A(compi,:),colors(1,:),sample_labels{1});
    plot_with_error(result.frex,result.mean_B(compi,:), ...
        result.sem_B(compi,:),colors(2,:),sample_labels{2});
    ylabel(result.measure,'FontWeight','bold')
    title(sprintf('%s landscape comparison - component #%d', ...
        result.measure,result.comps(compi)),'Interpreter','none')
    legend('Location','best','Box','off')
    grid on
    set(gca,'FontSize',12,'Box','off')

    nexttile, hold on
    yline(0,'--k','LineWidth',1);
    plot(result.frex,result.mean_difference(compi,:),'k','LineWidth',1.8);
    significant = result.significant(compi,:);
    if any(significant)
        scatter(result.frex(significant), ...
            result.mean_difference(compi,significant),42,colors(2,:), ...
            'filled','DisplayName',significance_label(correction,alpha));
        legend('Location','best','Box','off')
    end
    xlabel('Frequency (Hz)','FontWeight','bold')
    ylabel(sprintf('%s - %s',sample_labels{1},sample_labels{2}), ...
        'FontWeight','bold')
    grid on
    set(gca,'FontSize',12,'Box','off')
end

end

function plot_with_error(frex,average,sem,color,label)

valid = isfinite(average) & isfinite(sem);
if any(valid)
    fill([frex(valid) fliplr(frex(valid))], ...
        [average(valid)-sem(valid) fliplr(average(valid)+sem(valid))], ...
        color,'FaceAlpha',.18,'EdgeColor','none','HandleVisibility','off');
end
plot(frex,average,'Color',color,'LineWidth',2,'DisplayName',label);

end

function label = significance_label(correction,alpha)

if strcmp(correction,'fdr')
    label = sprintf('FDR q <= %.3g',alpha);
else
    label = sprintf('p <= %.3g',alpha);
end

end

%% Helper Function: Validation
function [design,test,tail] = validate_analysis_options(design,test,tail)

if isstring(design) && isscalar(design)
    design = char(design);
end
if ~ischar(design)
    error('design must be ''within'' or ''between''.');
end
if any(strcmpi(design,{'within','paired'}))
    design = 'within';
elseif any(strcmpi(design,{'between','independent'}))
    design = 'between';
else
    error('design must be ''within'' or ''between''.');
end

if isempty(test)
    if strcmp(design,'within')
        test = 'signrank';
    else
        test = 'ranksum';
    end
elseif isstring(test) && isscalar(test)
    test = char(test);
end
if strcmp(design,'within')
    test = validatestring(test,{'ttest','signrank'});
else
    test = validatestring(test,{'ttest','ranksum'});
end
tail = validatestring(tail,{'both','left','right'});

end

function labels = validate_sample_labels(labels)

if isstring(labels)
    labels = cellstr(labels);
end
if ~iscell(labels) || numel(labels)~=2 || ...
        ~all(cellfun(@(value)ischar(value) && ~isempty(value),labels))
    error('sample_labels must contain two non-empty text labels.');
end
labels = reshape(labels,1,2);

end

function validate_probability(value,name)

if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value<=0 || value>=1
    error('%s must lie strictly between 0 and 1.',name);
end

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
