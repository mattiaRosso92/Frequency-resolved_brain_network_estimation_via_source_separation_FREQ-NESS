function STATS = FREQNESS_LandscapeComparison(LANDSCAPE_A,LANDSCAPE_B,varargin)

% Preserve the established two-landscape interface while accepting one
% multidimensional cell array for generalized factorial designs.
if iscell(LANDSCAPE_A)
    if nargin>=2
        factorial_options = [{LANDSCAPE_B},varargin];
    else
        factorial_options = {};
    end
    STATS = factorial_landscape_comparison(LANDSCAPE_A, ...
        factorial_options{:});
    return
elseif nargin<2
    error('LANDSCAPE_A and LANDSCAPE_B are required.');
end

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
%  This function compares participant-level FREQ-NESS landscapes using
%  pairwise or generalized factorial designs. Statistical tests are
%  performed separately at every requested frequency, with
%  Benjamini-Hochberg FDR correction across frequencies for each effect,
%  measure, and component by default.
%
%  Recognized landscape fields are FREQ.evals, ENTR.H2, ENTR.ED,
%  CFC.amplitude_raw, CFC.amplitude_normalized, STAB.stability_index, and
%  STAB.mean_frequency. When 'measure' is omitted, every recognized field
%  shared by both inputs is analyzed separately and retains its original
%  field name in STATS. Raw numeric landscapes are also accepted when their
%  frequency axis is supplied through 'frex'.
%
%  Generalized designs use one cell array whose dimensions represent the
%  experimental factors and whose entries contain the landscape structure
%  for every combination of factor levels. For example, a 2-by-3 cell array
%  represents two factors with two and three levels. Each factor can be
%  declared as within- or between-subject. Existing calls using
%  LANDSCAPE_A and LANDSCAPE_B remain unchanged.
% ========================================================================

% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - LANDSCAPE_A, LANDSCAPE_B : structures containing participant-level
%            landscapes, or numeric arrays. Supported numeric layouts are
%            [nFrex x nParticipants] and
%            [nFeatures x nFrex x nParticipants].
%        Alternatively, the first input can be an N-dimensional cell array
%        containing one landscape per factorial cell, followed directly by
%        the name-value inputs.
%
%  - measure : field name or cell array of field names. When omitted, every
%              recognized shared landscape is analyzed.
%
%  - design : 'within' (default) or 'between'. For factorial cell inputs,
%             supply one value per factor, for example
%             {'between','within'}. A single value applies to all factors.
%
%  - test : 'ttest' or 'signrank' for within-subject designs; 'ttest' or
%           'ranksum' for between-subject designs. Defaults are 'signrank'
%           for within-subject and 'ranksum' for between-subject designs.
%           An incompatible non-parametric choice is corrected to the
%           design-appropriate test and reported with a warning.
%           Factorial cell inputs additionally support 'auto', 'anova',
%           'friedman', 'kruskalwallis', and 'permutation_anova'. 'auto'
%           preserves the default pairwise rank test for two levels and
%           otherwise uses factorial ANOVA.
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
%                    Retained for established pairwise calls.
%
%  - factor_names : names of the factorial dimensions. Defaults to
%                   {'Factor1','Factor2',...}.
%
%  - level_labels : one cell array of level names per factor. Defaults to
%                   {'Level1','Level2',...} within each factor.
%
%  - subject_ids : optional cell array matching the factorial input. Each
%                  entry contains participant identifiers for that
%                  landscape. Repeated observations are reordered and
%                  matched by identifier when supplied. Otherwise
%                  participant order is assumed for within-subject factors
%                  with a warning.
%
%  - effects : 'all' (default) or a cell array selecting main effects and
%              interactions by factor name, for example
%              {'Task','Group_x_Task'}.
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
%  - n_permutations : number of Monte Carlo permutations used only when
%                     test='permutation_anova'. Default: 10000.
%
%  - random_seed : random seed for permutation_anova. Default: 0.
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - STATS : structure containing analysis metadata and one field per
%            analyzed landscape. Established pairwise calls retain fields
%            such as STATS.evals. Generalized calls store every result under
%            STATS.effects.<effect>.<measure>. Results include test
%            statistics, effect sizes, uncorrected p values, FDR q values,
%            and significance.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso, Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 16/10/2026
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

%% Factorial Landscape Comparison
function STATS = factorial_landscape_comparison(LANDSCAPE_GRID,varargin)

opts = struct('measure', [], 'design', 'within', 'test', 'auto', ...
    'which_comp', 1, 'frex', [], 'frex2compare', [], ...
    'sample_labels', [], 'factor_names', [], 'level_labels', [], ...
    'subject_ids', [], 'effects', 'all', 'tail', 'both', ...
    'alpha', .05, 'correction', 'fdr', 'frex_tolerance', [], ...
    'n_permutations', 10000, 'random_seed', 0, ...
    'plot_results', true);
opts = parse_name_value_pairs(opts,varargin{:});

[LANDSCAPE_GRID,factor_sizes,cell_levels] = ...
    normalize_factorial_grid(LANDSCAPE_GRID);
nFactors = numel(factor_sizes);
factor_names = validate_factor_names(opts.factor_names,nFactors);
level_labels = validate_level_labels(opts.level_labels, ...
    opts.sample_labels,factor_names,factor_sizes);
factor_design = validate_factor_design(opts.design,nFactors);
[test,test_requested] = validate_factorial_test(opts.test, ...
    factor_design,factor_sizes);
effect_definitions = factorial_effect_definitions(factor_names,opts.effects);
tail = validatestring(opts.tail,{'both','left','right'});
if ~strcmp(tail,'both') && ~ismember(test,{'ttest','signrank','ranksum'})
    warning(['Factorial omnibus tests are non-directional. Using ' ...
        'tail=''both'' instead.']);
    tail = 'both';
end
validate_probability(opts.alpha,'alpha');
correction = validatestring(opts.correction,{'fdr','none'});
validate_positive_integer(opts.n_permutations,'n_permutations');
if ~isnumeric(opts.random_seed) || ~isscalar(opts.random_seed) || ...
        ~isfinite(opts.random_seed) || opts.random_seed~=round(opts.random_seed)
    error('random_seed must be one finite integer.');
end
if ~islogical(opts.plot_results) || ~isscalar(opts.plot_results)
    error('plot_results must be one logical value.');
end

measures = select_factorial_measures(LANDSCAPE_GRID,opts.measure);
[values,frequency_axes,participant_counts] = ...
    validate_factorial_landscapes(LANDSCAPE_GRID,measures,opts.frex);
[orders,subject_ids] = match_factorial_subjects(participant_counts, ...
    cell_levels,factor_design,opts.subject_ids);
for measurei = 1:numel(measures)
    for celli = 1:numel(LANDSCAPE_GRID)
        values{measurei,celli} = ...
            values{measurei,celli}(:,:,orders{celli});
    end
end

% Preserve the complete established implementation for generalized calls
% that still represent one two-level pairwise comparison.
pairwise_tests = {'ttest','signrank','ranksum'};
if nFactors==1 && factor_sizes(1)==2 && ismember(test,pairwise_tests)
    pairwise_options = {'measure',opts.measure, ...
        'design',factor_design{1},'test',test, ...
        'which_comp',opts.which_comp,'frex',opts.frex, ...
        'frex2compare',opts.frex2compare, ...
        'sample_labels',level_labels{1},'tail',tail, ...
        'alpha',opts.alpha,'correction',correction, ...
        'frex_tolerance',opts.frex_tolerance, ...
        'plot_results',opts.plot_results};
    reordered = LANDSCAPE_GRID;
    for celli = 1:2
        reordered{celli} = reorder_landscape_participants( ...
            LANDSCAPE_GRID{celli},measures,frequency_axes, ...
            orders{celli},opts.frex);
    end
    STATS = FREQNESS_LandscapeComparison(reordered{1},reordered{2}, ...
        pairwise_options{:});
    STATS.analysis = 'pairwise';
    STATS.test_requested = test_requested;
    STATS.factor_names = factor_names;
    STATS.level_labels = level_labels;
    STATS.factor_design = factor_design;
    STATS.input_shape = factor_sizes;
    STATS.subject_ids = subject_ids;
    definition = effect_definitions(1);
    STATS.effects = struct;
    STATS.effects.(definition.key) = struct('name',definition.name, ...
        'key',definition.key,'factors',definition.factors);
    for measurei = 1:numel(measures)
        measure = measures{measurei};
        STATS.effects.(definition.key).(measure) = STATS.(measure);
    end
    return
end

fprintf(['\nFREQNESS Landscape Comparison: testing a %s design with %d ' ...
    'factors, %d factorial cells, and %d measure(s) using %s.\n'], ...
    strjoin(factor_design,'/'),nFactors,numel(LANDSCAPE_GRID), ...
    numel(measures),test);

STATS = struct;
STATS.analysis = 'factorial';
STATS.factor_names = factor_names;
STATS.level_labels = level_labels;
STATS.factor_design = factor_design;
STATS.design = factor_design;
STATS.input_shape = factor_sizes;
STATS.test = test;
STATS.test_requested = test_requested;
STATS.tail = tail;
STATS.alpha = opts.alpha;
STATS.correction = correction;
STATS.correction_scope = ...
    'frequencies separately for each effect, measure, and component';
STATS.measures = measures;
STATS.subject_ids = subject_ids;
STATS.cell_n_participants = reshape(participant_counts,size(LANDSCAPE_GRID));
STATS.effects = struct;
STATS.effect_names = {effect_definitions.key};
STATS.cell_means = struct;
STATS.n_permutations = 0;
STATS.random_seed = opts.random_seed;

for effecti = 1:numel(effect_definitions)
    definition = effect_definitions(effecti);
    STATS.effects.(definition.key) = struct('name',definition.name, ...
        'key',definition.key,'factors',definition.factors);
end

for measurei = 1:numel(measures)
    measure = measures{measurei};
    frex = frequency_axes{measurei};
    nfeatures = size(values{measurei,1},1);
    comps = select_components(opts.which_comp,nfeatures);
    [requested_frex,frex2compare,idx_frex2compare] = ...
        select_frequencies(frex,opts.frex2compare,opts.frex_tolerance);
    selected = cell(size(LANDSCAPE_GRID));
    for celli = 1:numel(LANDSCAPE_GRID)
        selected{celli} = values{measurei,celli}(comps, ...
            idx_frex2compare,:);
    end
    [Y,observation_levels,subject_index,cell_means] = ...
        factorial_landscape_long_data(selected,cell_levels,factor_design);
    STATS.cell_means.(measure) = reshape(cell_means,size(LANDSCAPE_GRID));

    fprintf('  %s: %d frequencies, component(s) %s.\n',measure, ...
        numel(frex2compare),mat2str(comps));

    if ismember(test,{'friedman','kruskalwallis'})
        [effect_results,models] = factorial_rank_statistics(Y, ...
            observation_levels,subject_index,test,effect_definitions, ...
            numel(comps),numel(frex2compare));
    else
        [effect_results,models] = factorial_anova_statistics(Y, ...
            observation_levels,subject_index,factor_sizes,factor_design, ...
            effect_definitions,numel(comps),numel(frex2compare));
    end

    if strcmp(test,'permutation_anova')
        fprintf('    Computing %d design-aware Monte Carlo permutations.\n', ...
            opts.n_permutations);
        rng_state = rng;
        rng(opts.random_seed,'twister');
        rng_cleanup = onCleanup(@()rng(rng_state));
        for effecti = 1:numel(effect_results)
            effect_results(effecti).p_uncorrected = ...
                permutation_anova_p_values(models(effecti), ...
                effect_results(effecti).statistic,opts.n_permutations);
        end
        clear rng_cleanup
        STATS.n_permutations = opts.n_permutations;
    end

    for effecti = 1:numel(effect_results)
        result = correct_factorial_frequencies(effect_results(effecti), ...
            correction,opts.alpha);
        result.measure = measure;
        result.requested_frex = requested_frex;
        result.frex = frex2compare;
        result.idx_frex = idx_frex2compare;
        result.comps = comps;
        result.values = reshape(selected,size(LANDSCAPE_GRID));
        key = result.key;
        STATS.effects.(key).(measure) = result;
        if opts.plot_results
            plot_factorial_landscape_results(result,correction,opts.alpha);
        end
    end
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

%% Helper Function: Normalize Factorial Grid
function [LANDSCAPE_GRID,factor_sizes,cell_levels] = ...
    normalize_factorial_grid(LANDSCAPE_GRID)

if ~iscell(LANDSCAPE_GRID) || isempty(LANDSCAPE_GRID) || ...
        any(cellfun(@isempty,LANDSCAPE_GRID(:)))
    error(['Factorial input must be a non-empty cell array of landscape ' ...
        'structures or numeric arrays.']);
end
grid_size = size(LANDSCAPE_GRID);
factor_dimensions = find(grid_size>1);
if isempty(factor_dimensions)
    error('Factorial input must contain at least two cells.');
end
factor_sizes = grid_size(factor_dimensions);
if isscalar(factor_sizes)
    LANDSCAPE_GRID = reshape(LANDSCAPE_GRID,[factor_sizes 1]);
else
    LANDSCAPE_GRID = reshape(LANDSCAPE_GRID,factor_sizes);
end
nFactors = numel(factor_sizes);
cell_levels = zeros(numel(LANDSCAPE_GRID),nFactors);
grid_shape = factor_sizes;
if nFactors==1
    grid_shape = [factor_sizes 1];
end
for celli = 1:numel(LANDSCAPE_GRID)
    indices = cell(1,nFactors);
    [indices{:}] = ind2sub(grid_shape,celli);
    cell_levels(celli,:) = cellfun(@double,indices);
end

end

%% Helper Function: Validate Factor Names
function factor_names = validate_factor_names(factor_names,nFactors)

if isempty(factor_names)
    factor_names = arrayfun(@(value)sprintf('Factor%d',value), ...
        1:nFactors,'UniformOutput',false);
elseif isstring(factor_names)
    factor_names = cellstr(factor_names);
end
if ~iscell(factor_names) || numel(factor_names)~=nFactors || ...
        ~all(cellfun(@(value)ischar(value) && ~isempty(strtrim(value)), ...
        factor_names))
    error('factor_names must contain one non-empty text name per factor.');
end
factor_names = reshape(factor_names,1,[]);
keys = cellfun(@matlab.lang.makeValidName,factor_names, ...
    'UniformOutput',false);
if numel(unique(lower(string(keys))))~=nFactors
    error('factor_names must remain unique when converted to MATLAB fields.');
end

end

%% Helper Function: Validate Level Labels
function level_labels = validate_level_labels(level_labels,sample_labels, ...
    factor_names,factor_sizes)

nFactors = numel(factor_sizes);
if isempty(level_labels)
    level_labels = cell(1,nFactors);
    for factori = 1:nFactors
        level_labels{factori} = arrayfun(@(value)sprintf('Level%d',value), ...
            1:factor_sizes(factori),'UniformOutput',false);
    end
    if ~isempty(sample_labels)
        if nFactors~=1 || factor_sizes(1)~=2
            error(['sample_labels can only be used with a one-factor, ' ...
                'two-level cell input. Use level_labels instead.']);
        end
        level_labels{1} = validate_sample_labels(sample_labels);
    end
elseif ~iscell(level_labels) || numel(level_labels)~=nFactors
    error('level_labels must contain one cell array per factor.');
end

for factori = 1:nFactors
    labels = level_labels{factori};
    if isstring(labels)
        labels = cellstr(labels);
    end
    if ~iscell(labels) || numel(labels)~=factor_sizes(factori) || ...
            ~all(cellfun(@(value)ischar(value) && ...
            ~isempty(strtrim(value)),labels))
        error(['level_labels{%d} must contain %d non-empty labels for ' ...
            '%s.'],factori,factor_sizes(factori),factor_names{factori});
    end
    level_labels{factori} = reshape(labels,1,[]);
end

end

%% Helper Function: Validate Factorial Design
function factor_design = validate_factor_design(design,nFactors)

if ischar(design) || (isstring(design) && isscalar(design))
    factor_design = repmat({char(design)},1,nFactors);
elseif isstring(design)
    factor_design = cellstr(design);
else
    factor_design = design;
end
if ~iscell(factor_design) || numel(factor_design)~=nFactors
    error('design must contain one ''within'' or ''between'' value per factor.');
end
for factori = 1:nFactors
    value = factor_design{factori};
    if isstring(value) && isscalar(value)
        value = char(value);
    end
    if ~ischar(value)
        error('design values must be ''within'' or ''between''.');
    elseif any(strcmpi(value,{'within','paired'}))
        factor_design{factori} = 'within';
    elseif any(strcmpi(value,{'between','independent'}))
        factor_design{factori} = 'between';
    else
        error('design values must be ''within'' or ''between''.');
    end
end
factor_design = reshape(factor_design,1,[]);

end

%% Helper Function: Validate Factorial Test
function [test,test_requested] = validate_factorial_test(test,design,factor_sizes)

if isempty(test)
    test = 'auto';
elseif isstring(test) && isscalar(test)
    test = char(test);
end
test = validatestring(test,{'auto','ttest','signrank','ranksum', ...
    'anova','friedman','kruskalwallis','permutation_anova'});
test_requested = test;
nFactors = numel(factor_sizes);

if strcmp(test,'auto')
    if nFactors==1 && factor_sizes(1)==2
        if strcmp(design{1},'within')
            test = 'signrank';
        else
            test = 'ranksum';
        end
    else
        test = 'anova';
    end
    return
end

if nFactors==1 && factor_sizes(1)==2
    if strcmp(design{1},'within') && strcmp(test,'ranksum')
        warning(['test=''ranksum'' is incompatible with a within-subject ' ...
            'factor. Using test=''signrank'' instead.']);
        test = 'signrank';
    elseif strcmp(design{1},'between') && strcmp(test,'signrank')
        warning(['test=''signrank'' is incompatible with a between-subject ' ...
            'factor. Using test=''ranksum'' instead.']);
        test = 'ranksum';
    elseif strcmp(test,'friedman')
        warning('Two-level Friedman analysis is equivalent to signrank here.');
        test = 'signrank';
    elseif strcmp(test,'kruskalwallis')
        warning('Two-level Kruskal-Wallis analysis is equivalent to ranksum here.');
        test = 'ranksum';
    end
    return
end

if nFactors==1
    if strcmp(test,'ttest')
        warning('ttest supports two levels only. Using anova instead.');
        test = 'anova';
    elseif strcmp(test,'signrank')
        warning('signrank supports two levels only. Using friedman instead.');
        test = 'friedman';
    elseif strcmp(test,'ranksum')
        warning('ranksum supports two levels only. Using kruskalwallis instead.');
        test = 'kruskalwallis';
    end
    if strcmp(design{1},'within') && strcmp(test,'kruskalwallis')
        warning(['kruskalwallis is an independent-samples test. Using ' ...
            'friedman for this within-subject factor instead.']);
        test = 'friedman';
    elseif strcmp(design{1},'between') && strcmp(test,'friedman')
        warning(['friedman is a repeated-measures test. Using ' ...
            'kruskalwallis for this between-subject factor instead.']);
        test = 'kruskalwallis';
    end
else
    if ismember(test,{'signrank','ranksum','friedman','kruskalwallis'})
        warning(['The selected rank test does not model factorial ' ...
            'interactions. Using permutation_anova instead.']);
        test = 'permutation_anova';
    elseif strcmp(test,'ttest')
        warning('ttest does not model factorial effects. Using anova instead.');
        test = 'anova';
    end
end

end

%% Helper Function: Define Factorial Effects
function effects = factorial_effect_definitions(factor_names,requested)

nFactors = numel(factor_names);
effects = repmat(struct('name','','key','','factors',[]), ...
    1,2^nFactors-1);
for mask = 1:(2^nFactors-1)
    factors = find(bitget(mask,1:nFactors));
    names = factor_names(factors);
    key_names = cellfun(@matlab.lang.makeValidName,names, ...
        'UniformOutput',false);
    effects(mask) = struct('name',strjoin(names,' x '), ...
        'key',strjoin(key_names,'_x_'),'factors',factors);
end

if ischar(requested) || (isstring(requested) && isscalar(requested))
    if strcmpi(requested,'all')
        return
    end
    requested = {char(requested)};
elseif isstring(requested)
    requested = cellstr(requested);
end
if ~iscell(requested) || isempty(requested) || ...
        ~all(cellfun(@(value)ischar(value) && ~isempty(value),requested))
    error('effects must be ''all'' or a non-empty cell array of effect names.');
end

keep = false(1,numel(effects));
for requesti = 1:numel(requested)
    value = requested{requesti};
    matches = strcmpi(value,{effects.key}) | strcmpi(value,{effects.name});
    if ~any(matches)
        value = regexprep(value,'\s+[xX]\s+|[:*]','_x_');
        value = matlab.lang.makeValidName(value);
        matches = strcmpi(value,{effects.key});
    end
    if ~any(matches)
        error('Unknown factorial effect: %s.',requested{requesti});
    end
    keep = keep | matches;
end
effects = effects(keep);

end

%% Helper Function: Select Factorial Measures
function measures = select_factorial_measures(LANDSCAPE_GRID,measure)

numeric_cells = cellfun(@isnumeric,LANDSCAPE_GRID(:));
if any(numeric_cells) && ~all(numeric_cells)
    error('All factorial cells must use the same input format.');
end
if all(numeric_cells)
    measures = select_measures(LANDSCAPE_GRID{1},LANDSCAPE_GRID{2},measure);
    return
end
if ~all(cellfun(@isstruct,LANDSCAPE_GRID(:)))
    error('Factorial cells must contain structures or real numeric arrays.');
end

if isempty(measure)
    recognized = {'evals','H2','ED','amplitude_raw', ...
        'amplitude_normalized','stability_index','mean_frequency'};
    measures = cell(1,0);
    for fieldi = 1:numel(recognized)
        actual = resolve_field(LANDSCAPE_GRID{1},recognized{fieldi});
        if isempty(actual)
            continue
        end
        shared = true;
        for celli = 2:numel(LANDSCAPE_GRID)
            shared = shared && ...
                ~isempty(resolve_field(LANDSCAPE_GRID{celli},actual));
        end
        if shared
            measures{end+1} = actual; %#ok<AGROW>
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
    actual = resolve_field(LANDSCAPE_GRID{1},requested{fieldi});
    if isempty(actual)
        error('Every input must contain the measure ''%s''.',requested{fieldi});
    end
    for celli = 2:numel(LANDSCAPE_GRID)
        current = resolve_field(LANDSCAPE_GRID{celli},requested{fieldi});
        if isempty(current) || ~strcmpi(current,actual)
            error('Every input must contain the measure ''%s''.',requested{fieldi});
        end
    end
    if ~any(strcmp(measures,actual))
        measures{end+1} = actual; %#ok<AGROW>
    end
end

end

%% Helper Function: Validate Factorial Landscapes
function [values,frequency_axes,participant_counts] = ...
    validate_factorial_landscapes(LANDSCAPE_GRID,measures,supplied_frex)

nMeasures = numel(measures);
nCells = numel(LANDSCAPE_GRID);
values = cell(nMeasures,nCells);
frequency_axes = cell(1,nMeasures);
participant_counts = zeros(nCells,1);

for measurei = 1:nMeasures
    measure = measures{measurei};
    reference_size = [];
    for celli = 1:nCells
        current_frex = frequency_axis(LANDSCAPE_GRID{celli},measure, ...
            supplied_frex);
        current = canonical_values(LANDSCAPE_GRID{celli},measure, ...
            current_frex,sprintf('LANDSCAPE_GRID{%d}',celli));
        if celli==1
            frequency_axes{measurei} = current_frex;
            reference_size = size(current);
            reference_size(end+1:3) = 1;
        elseif numel(current_frex)~=numel(frequency_axes{measurei}) || ...
                any(abs(current_frex-frequency_axes{measurei})>1e-10) || ...
                size(current,1)~=reference_size(1)
            error(['All factorial inputs must share the feature dimension ' ...
                'and frequency grid for measure ''%s''.'],measure);
        end
        current_n = size(current,3);
        if measurei==1
            participant_counts(celli) = current_n;
            if current_n<2
                error('Every factorial cell must contain at least two participants.');
            end
        elseif current_n~=participant_counts(celli)
            error(['Every selected measure in factorial cell %d must have ' ...
                'the same participant count.'],celli);
        end
        values{measurei,celli} = current;
    end
end

end

%% Helper Function: Match Repeated Participants
function [orders,subject_ids] = match_factorial_subjects(nsubs, ...
    cell_levels,design,subject_ids)

nCells = numel(nsubs);
orders = cell(nCells,1);
for celli = 1:nCells
    orders{celli} = 1:nsubs(celli);
end
within_factors = find(strcmp(design,'within'));
between_factors = find(strcmp(design,'between'));
if isempty(subject_ids)
    subject_ids = cell(nCells,1);
    for celli = 1:nCells
        subject_ids{celli} = arrayfun(@(value)sprintf('S%d',value), ...
            1:nsubs(celli),'UniformOutput',false);
    end
    if ~isempty(within_factors)
        warning(['subject_ids were not supplied. Repeated observations are ' ...
            'assumed to have identical participant order within each ' ...
            'between-subject stratum.']);
    end
else
    if ~iscell(subject_ids) || numel(subject_ids)~=nCells
        error('subject_ids must be a cell array matching the factorial input.');
    end
    subject_ids = reshape(subject_ids,[],1);
    for celli = 1:nCells
        subject_ids{celli} = normalize_subject_ids(subject_ids{celli}, ...
            nsubs(celli),celli);
    end
end

if isempty(within_factors)
    subject_ids = reshape(subject_ids,size(nsubs));
    orders = reshape(orders,size(nsubs));
    return
end
if isempty(between_factors)
    strata = ones(nCells,1);
else
    [~,~,strata] = unique(cell_levels(:,between_factors),'rows','stable');
end
for stratum = unique(strata(:))'
    cells = find(strata==stratum);
    reference_ids = subject_ids{cells(1)};
    for celli = cells(:)'
        ids = subject_ids{celli};
        if numel(ids)~=numel(reference_ids)
            error(['Within-subject cells in the same between-subject stratum ' ...
                'must contain matching participant counts.']);
        end
        [present,order] = ismember(reference_ids,ids);
        if ~all(present) || numel(unique(order))~=numel(order)
            error(['subject_ids must identify the same participants across ' ...
                'all repeated-measures cells in a stratum.']);
        end
        orders{celli} = order;
        subject_ids{celli} = reference_ids;
    end
end
subject_ids = reshape(subject_ids,size(nsubs));
orders = reshape(orders,size(nsubs));

end

function ids = normalize_subject_ids(ids,nsubs,celli)

if isnumeric(ids) || islogical(ids) || isstring(ids) || iscell(ids)
    ids = cellstr(string(ids(:)));
else
    error('subject_ids{%d} must contain numeric or text identifiers.',celli);
end
if numel(ids)~=nsubs || any(cellfun(@isempty,ids)) || ...
        numel(unique(ids))~=numel(ids)
    error(['subject_ids{%d} must contain one unique, non-empty identifier ' ...
        'per participant.'],celli);
end
ids = reshape(ids,1,[]);

end

function output = reorder_landscape_participants(input,measures, ...
    frequency_axes,order,supplied_frex)

if isnumeric(input)
    output = canonical_values(input,'values',frequency_axes{1},'LANDSCAPE');
    output = output(:,:,order);
    return
end
output = input;
for measurei = 1:numel(measures)
    measure = measures{measurei};
    actual = resolve_field(output,measure);
    current = canonical_values(output,actual,frequency_axes{measurei}, ...
        'LANDSCAPE');
    output.(actual) = current(:,:,order);
end
if ~isempty(supplied_frex) && isempty(resolve_field(output,'frex'))
    output.frex = supplied_frex;
end

end

%% Helper Function: Assemble Long Factorial Data
function [Y,levels,subject_index,cell_means] = ...
    factorial_landscape_long_data(values,cell_levels,design)

nCells = numel(values);
nFactors = size(cell_levels,2);
within_factors = find(strcmp(design,'within'));
between_factors = find(strcmp(design,'between'));
cell_n = cellfun(@(value)size(value,3),values(:));
cell_means = cell(size(values));
Y = nan(sum(cell_n),size(values{1},1)*size(values{1},2));
levels = zeros(sum(cell_n),nFactors);
subject_index = zeros(sum(cell_n),1);
row_start = 1;

if isempty(within_factors)
    strata = (1:nCells)';
elseif isempty(between_factors)
    strata = ones(nCells,1);
else
    [~,~,strata] = unique(cell_levels(:,between_factors),'rows','stable');
end
stratum_values = unique(strata(:))';
stratum_offsets = zeros(max(stratum_values),1);
next_subject = 0;
if ~isempty(within_factors)
    for stratum = stratum_values
        first_cell = find(strata==stratum,1);
        stratum_offsets(stratum) = next_subject;
        next_subject = next_subject+cell_n(first_cell);
    end
end

for celli = 1:nCells
    current = values{celli};
    cell_means{celli} = mean(current,3,'omitnan');
    nsubjects = size(current,3);
    rows = row_start:(row_start+nsubjects-1);
    Y(rows,:) = reshape(permute(current,[3 1 2]),nsubjects,[]);
    levels(rows,:) = repmat(cell_levels(celli,:),nsubjects,1);
    if isempty(within_factors)
        subject_index(rows) = next_subject+(1:nsubjects);
        next_subject = next_subject+nsubjects;
    else
        subject_index(rows) = ...
            stratum_offsets(strata(celli))+(1:nsubjects);
    end
    row_start = row_start+nsubjects;
end

end

%% Helper Function: Factorial ANOVA Statistics
function [results,models] = factorial_anova_statistics(Y,levels, ...
    subject_index,factor_sizes,design,effects,ncomps,nfrex)

nFactors = numel(factor_sizes);
factor_codings = cell(1,nFactors);
for factori = 1:nFactors
    factor_codings{factori} = effect_coding(levels(:,factori), ...
        factor_sizes(factori));
end
all_definitions = numeric_effect_definitions(nFactors);
all_matrices = cell(1,numel(all_definitions));
for effecti = 1:numel(all_definitions)
    all_matrices{effecti} = interaction_matrix( ...
        factor_codings(all_definitions{effecti}));
end

results = repmat(empty_factorial_effect,1,numel(effects));
models = repmat(empty_factorial_model,1,numel(effects));
for effecti = 1:numel(effects)
    target_factors = effects(effecti).factors;
    between_only = all(strcmp(design(target_factors),'between'));
    if between_only
        [current_Y,current_levels,current_subjects] = ...
            subject_average_data(Y,levels,subject_index);
        current_codings = cell(1,nFactors);
        for factori = 1:nFactors
            current_codings{factori} = effect_coding( ...
                current_levels(:,factori),factor_sizes(factori));
        end
        included = cellfun(@(factors) ...
            all(strcmp(design(factors),'between')),all_definitions);
        X_full = ones(size(current_Y,1),1);
        target_columns = [];
        for definitioni = find(included)
            columns = interaction_matrix( ...
                current_codings(all_definitions{definitioni}));
            new_columns = size(X_full,2)+(1:size(columns,2));
            X_full = [X_full columns]; %#ok<AGROW>
            if isequal(all_definitions{definitioni},target_factors)
                target_columns = new_columns;
            end
        end
        permutation_mode = 'between';
        block_index = current_subjects;
    else
        current_Y = Y;
        block_index = subject_index;
        subjects = unique(subject_index,'stable');
        X_full = zeros(size(Y,1),numel(subjects));
        for subi = 1:numel(subjects)
            X_full(subject_index==subjects(subi),subi) = 1;
        end
        target_columns = [];
        for definitioni = 1:numel(all_definitions)
            factors = all_definitions{definitioni};
            if ~any(strcmp(design(factors),'within'))
                continue
            end
            columns = all_matrices{definitioni};
            new_columns = size(X_full,2)+(1:size(columns,2));
            X_full = [X_full columns]; %#ok<AGROW>
            if isequal(factors,target_factors)
                target_columns = new_columns;
            end
        end
        permutation_mode = 'within';
    end
    if isempty(target_columns)
        error('The requested effect %s is not estimable.',effects(effecti).name);
    end
    X_reduced = X_full;
    X_reduced(:,target_columns) = [];
    [statistic,p_values,effect_size,df1,df2] = ...
        partial_f_columns(current_Y,X_full,X_reduced);
    results(effecti) = initialize_factorial_effect(effects(effecti), ...
        statistic,p_values,effect_size,'F','partial eta squared', ...
        df1,df2,ncomps,nfrex);
    models(effecti) = struct('Y',current_Y,'X_full',X_full, ...
        'X_reduced',X_reduced,'permutation_mode',permutation_mode, ...
        'block_index',block_index);
end

end

function coding = effect_coding(levels,nlevels)

coding = zeros(numel(levels),nlevels-1);
for leveli = 1:(nlevels-1)
    coding(levels==leveli,leveli) = 1;
end
coding(levels==nlevels,:) = -1;

end

function definitions = numeric_effect_definitions(nFactors)

definitions = cell(1,2^nFactors-1);
for mask = 1:(2^nFactors-1)
    definitions{mask} = find(bitget(mask,1:nFactors));
end

end

function output = interaction_matrix(matrices)

output = ones(size(matrices{1},1),1);
for matrixi = 1:numel(matrices)
    current = matrices{matrixi};
    combined = zeros(size(output,1),size(output,2)*size(current,2));
    column = 1;
    for lefti = 1:size(output,2)
        for righti = 1:size(current,2)
            combined(:,column) = output(:,lefti).*current(:,righti);
            column = column+1;
        end
    end
    output = combined;
end

end

function [Y_subject,levels_subject,subjects] = ...
    subject_average_data(Y,levels,subject_index)

subjects = unique(subject_index,'stable');
Y_subject = nan(numel(subjects),size(Y,2));
levels_subject = zeros(numel(subjects),size(levels,2));
for subi = 1:numel(subjects)
    rows = subject_index==subjects(subi);
    Y_subject(subi,:) = mean(Y(rows,:),1,'omitnan');
    levels_subject(subi,:) = levels(find(rows,1),:);
end

end

function [statistic,p_values,effect_size,df1,df2] = ...
    partial_f_columns(Y,X_full,X_reduced)

ncolumns = size(Y,2);
statistic = nan(1,ncolumns);
p_values = nan(1,ncolumns);
effect_size = nan(1,ncolumns);
df1 = nan(1,ncolumns);
df2 = nan(1,ncolumns);
for column = 1:ncolumns
    valid = isfinite(Y(:,column));
    [statistic(column),p_values(column),effect_size(column), ...
        df1(column),df2(column)] = partial_f_single(Y(valid,column), ...
        X_full(valid,:),X_reduced(valid,:));
end

end

function [statistic,p_value,effect_size,df1,df2] = ...
    partial_f_single(Y,X_full,X_reduced)

statistic = NaN;
p_value = NaN;
effect_size = NaN;
df1 = NaN;
df2 = NaN;
if numel(Y)<2
    return
end
tolerance = max(size(X_full))*eps(max(norm(X_full),1));
rank_full = rank(X_full,tolerance);
rank_reduced = rank(X_reduced,tolerance);
df1 = rank_full-rank_reduced;
df2 = size(X_full,1)-rank_full;
if df1<1 || df2<1
    return
end
residual_full = Y-X_full*(pinv(X_full)*Y);
residual_reduced = Y-X_reduced*(pinv(X_reduced)*Y);
sse_full = sum(residual_full.^2);
sse_reduced = sum(residual_reduced.^2);
numerator = max(sse_reduced-sse_full,0)/df1;
denominator = sse_full/df2;
if denominator>eps
    statistic = numerator/denominator;
elseif numerator>eps
    statistic = Inf;
else
    statistic = 0;
end
p_value = 1-fcdf(statistic,df1,df2);
denominator_eta = statistic*df1+df2;
if denominator_eta>0
    effect_size = statistic*df1/denominator_eta;
else
    effect_size = 0;
end

end

%% Helper Function: Factorial Rank Statistics
function [results,models] = factorial_rank_statistics(Y,levels, ...
    subject_index,test,effects,ncomps,nfrex)

if size(levels,2)~=1 || numel(effects)~=1
    error('friedman and kruskalwallis support one factorial effect only.');
end
nlevels = max(levels(:,1));
definition = effects(1);
model = empty_factorial_model;

if strcmp(test,'kruskalwallis')
    [statistic,p_values,effect_size,df1] = ...
        kruskalwallis_columns(Y,levels(:,1),nlevels);
    effect_name = 'epsilon squared';
    statistic_name = 'H';
else
    subjects = unique(subject_index,'stable');
    repeated = nan(numel(subjects),nlevels,size(Y,2));
    for subi = 1:numel(subjects)
        for leveli = 1:nlevels
            row = subject_index==subjects(subi) & levels(:,1)==leveli;
            if sum(row)~=1
                error(['Friedman analysis requires one observation per ' ...
                    'participant at every factor level.']);
            end
            repeated(subi,leveli,:) = Y(row,:);
        end
    end
    [statistic,p_values,effect_size,df1] = friedman_columns(repeated);
    effect_name = 'Kendall''s W';
    statistic_name = 'chi-square';
end

results = initialize_factorial_effect(definition,statistic,p_values, ...
    effect_size,statistic_name,effect_name,df1,inf(size(df1)), ...
    ncomps,nfrex);
models = model;

end

function [statistic,p_values,effect_size,df1] = ...
    kruskalwallis_columns(Y,groups,nlevels)

ncolumns = size(Y,2);
statistic = nan(1,ncolumns);
p_values = nan(1,ncolumns);
effect_size = nan(1,ncolumns);
df1 = repmat(nlevels-1,1,ncolumns);
for column = 1:ncolumns
    valid = isfinite(Y(:,column));
    values = Y(valid,column);
    current_groups = groups(valid);
    group_n = accumarray(current_groups,1,[nlevels 1]);
    if numel(values)<nlevels+1 || any(group_n==0)
        continue
    end
    [ranks,tie_sum] = rank_columns(values);
    rank_sums = accumarray(current_groups,ranks,[nlevels 1]);
    H = 12/(numel(values)*(numel(values)+1))* ...
        sum(rank_sums.^2./group_n)-3*(numel(values)+1);
    tie_correction = 1-tie_sum/(numel(values)^3-numel(values));
    if tie_correction>eps
        H = H/tie_correction;
    else
        H = 0;
    end
    statistic(column) = max(H,0);
    p_values(column) = 1-chi2cdf(statistic(column),nlevels-1);
    effect_size(column) = max((statistic(column)-nlevels+1)/ ...
        (numel(values)-nlevels),0);
end

end

function [statistic,p_values,effect_size,df1] = friedman_columns(Y)

[~,nlevels,ncolumns] = size(Y);
statistic = nan(1,ncolumns);
p_values = nan(1,ncolumns);
effect_size = nan(1,ncolumns);
df1 = repmat(nlevels-1,1,ncolumns);
for column = 1:ncolumns
    current = Y(:,:,column);
    current = current(all(isfinite(current),2),:);
    n = size(current,1);
    if n<2
        continue
    end
    rank_sums = zeros(nlevels,1);
    tie_sum = 0;
    for subi = 1:n
        [ranks,current_ties] = rank_columns(current(subi,:)');
        rank_sums = rank_sums+ranks;
        tie_sum = tie_sum+current_ties;
    end
    value = 12/(n*nlevels*(nlevels+1))*sum(rank_sums.^2)- ...
        3*n*(nlevels+1);
    tie_correction = 1-tie_sum/(n*nlevels*(nlevels^2-1));
    if tie_correction>eps
        value = value/tie_correction;
    else
        value = 0;
    end
    statistic(column) = max(value,0);
    p_values(column) = 1-chi2cdf(statistic(column),nlevels-1);
    effect_size(column) = statistic(column)/(n*(nlevels-1));
end

end

function [ranks,tie_sum] = rank_columns(values)

[nrows,ncolumns] = size(values);
ranks = zeros(nrows,ncolumns);
tie_sum = zeros(1,ncolumns);
for column = 1:ncolumns
    [sorted,order] = sort(values(:,column));
    current = 1;
    while current<=nrows
        last = current;
        while last<nrows && sorted(last+1)==sorted(current)
            last = last+1;
        end
        ranks(order(current:last),column) = (current+last)/2;
        count = last-current+1;
        if count>1
            tie_sum(column) = tie_sum(column)+count^3-count;
        end
        current = last+1;
    end
end

end

%% Helper Function: Initialize Factorial Result
function result = initialize_factorial_effect(definition,statistic,p_values, ...
    effect_size,statistic_name,effect_size_name,df1,df2,ncomps,nfrex)

result = empty_factorial_effect;
result.name = definition.name;
result.key = definition.key;
result.factors = definition.factors;
result.statistic = reshape(statistic,[ncomps nfrex]);
result.statistic_name = statistic_name;
result.effect_size = reshape(effect_size,[ncomps nfrex]);
result.effect_size_name = effect_size_name;
result.p_uncorrected = reshape(p_values,[ncomps nfrex]);
result.df1 = reshape(df1,[ncomps nfrex]);
result.df2 = reshape(df2,[ncomps nfrex]);

end

function result = empty_factorial_effect

result = struct('name','','key','','factors',[],'measure','', ...
    'requested_frex',[],'frex',[],'idx_frex',[],'comps',[], ...
    'values',[],'statistic',[],'statistic_name','', ...
    'effect_size',[],'effect_size_name','','p_uncorrected',[], ...
    'q_fdr',[],'significant',[],'fdr_threshold',[], ...
    'df1',[],'df2',[]);

end

function model = empty_factorial_model

model = struct('Y',[],'X_full',[],'X_reduced',[], ...
    'permutation_mode','','block_index',[]);

end

%% Helper Function: Permutation ANOVA
function p_values = permutation_anova_p_values(model,observed,n_permutations)

output_size = size(observed);
observed = observed(:)';
p_values = nan(size(observed));
for column = 1:numel(observed)
    if isnan(observed(column))
        continue
    end
    valid = isfinite(model.Y(:,column));
    Y = model.Y(valid,column);
    X_full = model.X_full(valid,:);
    X_reduced = model.X_reduced(valid,:);
    blocks = model.block_index(valid);
    [~,~,~,df1,df2] = partial_f_single(Y,X_full,X_reduced);
    if ~isfinite(df1) || ~isfinite(df2)
        continue
    end
    residual_reduced = Y-X_reduced*(pinv(X_reduced)*Y);
    fitted_reduced = Y-residual_reduced;
    count = 0;
    for permi = 1:n_permutations
        if strcmp(model.permutation_mode,'between')
            order = randperm(numel(Y));
        else
            order = (1:numel(Y))';
            subjects = unique(blocks,'stable');
            for subi = 1:numel(subjects)
                rows = find(blocks==subjects(subi));
                order(rows) = rows(randperm(numel(rows)));
            end
        end
        permuted = fitted_reduced+residual_reduced(order);
        current = partial_f_single(permuted,X_full,X_reduced);
        count = count+(current>=observed(column));
    end
    p_values(column) = (count+1)/(n_permutations+1);
end
p_values = reshape(p_values,output_size);

end

%% Helper Function: Correct Factorial Frequencies
function result = correct_factorial_frequencies(result,correction,alpha)

q_fdr = nan(size(result.p_uncorrected));
fdr_threshold = nan(size(result.p_uncorrected,1),1);
if strcmp(correction,'fdr')
    for compi = 1:size(result.p_uncorrected,1)
        [q_fdr(compi,:),fdr_threshold(compi)] = ...
            bh_fdr(result.p_uncorrected(compi,:),alpha);
    end
else
    q_fdr = result.p_uncorrected;
    fdr_threshold(:) = alpha;
end
result.q_fdr = q_fdr;
result.significant = isfinite(q_fdr) & q_fdr<=alpha;
result.fdr_threshold = fdr_threshold;

end

%% Helper Function: Plot Factorial Results
function plot_factorial_landscape_results(result,correction,alpha)

primary = [0 .45 .70];
accent = [.84 .37 0];
measure_label = landscape_measure_label(result.measure);
for compi = 1:numel(result.comps)
    figure('Color','w','Position',[100 100 920 530],'Renderer','painters');
    hold on
    values = result.statistic(compi,:);
    fill([result.frex fliplr(result.frex)], ...
        [zeros(size(values)) fliplr(values)],primary,'FaceAlpha',.09, ...
        'EdgeColor','none','HandleVisibility','off');
    plot(result.frex,values,'Color',primary,'LineWidth',2.4, ...
        'Marker','o','MarkerSize',5,'MarkerFaceColor','w', ...
        'MarkerEdgeColor',primary,'DisplayName',result.statistic_name);
    significant = result.significant(compi,:);
    if any(significant)
        scatter(result.frex(significant), ...
            values(significant),58,accent,'filled','MarkerEdgeColor','w', ...
            'LineWidth',.9,'DisplayName',significance_label(correction,alpha));
    end
    xlabel('Frequency (Hz)','FontWeight','bold','Color',[.09 .13 .17])
    ylabel(sprintf('%s statistic',result.statistic_name), ...
        'FontWeight','bold','Color',[.09 .13 .17])
    title({sprintf('%s comparison - %s',measure_label,result.name), ...
        sprintf('Component %d',result.comps(compi))}, ...
        'FontSize',14,'FontWeight','bold','Interpreter','none')
    xlim([result.frex(1) result.frex(end)])
    current_limits = ylim;
    ylim([0 max(current_limits(2),eps)])
    legend('Location','northeast','Box','off')
    style_publication_landscape_axis
end

end

function validate_positive_integer(value,name)

if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value<1 || value~=round(value)
    error('%s must be one positive integer.',name);
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
    if isscalar(frex)
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

colors = [0 .45 .70; .84 .37 0];
accent = [0 .62 .45];
measure_label = landscape_measure_label(result.measure);
for compi = 1:numel(result.comps)
    figure('Color','w','Position',[100 100 920 700],'Renderer','painters');
    layout = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile, hold on
    plot_with_error(result.frex,result.mean_A(compi,:), ...
        result.sem_A(compi,:),colors(1,:),sample_labels{1});
    plot_with_error(result.frex,result.mean_B(compi,:), ...
        result.sem_B(compi,:),colors(2,:),sample_labels{2});
    ylabel(measure_label,'FontWeight','bold','Color',[.09 .13 .17])
    legend('Location','best','Box','off','NumColumns',2)
    style_publication_landscape_axis

    nexttile, hold on
    difference = result.mean_difference(compi,:);
    yline(0,'--','Color',[.40 .44 .49],'LineWidth',.9);
    fill([result.frex fliplr(result.frex)], ...
        [zeros(size(difference)) fliplr(difference)],[.39 .45 .54], ...
        'FaceAlpha',.10,'EdgeColor','none','HandleVisibility','off');
    plot(result.frex,difference,'Color',[.20 .25 .31], ...
        'LineWidth',2,'Marker','o','MarkerSize',4, ...
        'MarkerFaceColor','w','MarkerEdgeColor',[.20 .25 .31]);
    significant = result.significant(compi,:);
    if any(significant)
        scatter(result.frex(significant), ...
            difference(significant),52,accent,'filled','MarkerEdgeColor','w', ...
            'LineWidth',.8,'DisplayName',significance_label(correction,alpha));
        legend('Location','best','Box','off')
    end
    xlabel('Frequency (Hz)','FontWeight','bold','Color',[.09 .13 .17])
    ylabel(sprintf('%s - %s',sample_labels{1},sample_labels{2}), ...
        'FontWeight','bold','Color',[.09 .13 .17])
    xlim([result.frex(1) result.frex(end)])
    style_publication_landscape_axis
    title(layout,{sprintf('%s comparison',measure_label), ...
        sprintf('Component %d',result.comps(compi))}, ...
        'FontSize',15,'FontWeight','bold','Interpreter','none')
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
plot(frex,average,'LineStyle','none','Marker','o','MarkerSize',4.5, ...
    'MarkerFaceColor','w','MarkerEdgeColor',color, ...
    'HandleVisibility','off');

end

function style_publication_landscape_axis

grid on
set(gca,'FontName','Helvetica','FontSize',10,'Box','off', ...
    'XGrid','off','YGrid','on','GridColor',[.85 .87 .90], ...
    'GridAlpha',.75,'LineWidth',.8,'Layer','top', ...
    'XColor',[.20 .25 .31],'YColor',[.20 .25 .31])

end

function label = landscape_measure_label(measure)

switch lower(measure)
    case 'evals'
        label = 'Network landscape';
    case 'h2'
        label = 'Renyi entropy landscape';
    case 'ed'
        label = 'Effective dimensionality landscape';
    case 'amplitude_raw'
        label = 'Cross-frequency coupling landscape';
    case 'amplitude_normalized'
        label = 'Normalized cross-frequency coupling landscape';
    case 'stability_index'
        label = 'Network stability landscape';
    case 'mean_frequency'
        label = 'Mean instantaneous frequency landscape';
    case 'values'
        label = 'Landscape';
    otherwise
        label = strrep(measure,'_',' ');
end

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
test = validatestring(test,{'ttest','signrank','ranksum'});
if strcmp(design,'within') && strcmp(test,'ranksum')
    warning(['test=''ranksum'' is an independent-samples test and is ' ...
        'incompatible with design=''within''. Using test=''signrank'' instead.']);
    test = 'signrank';
elseif strcmp(design,'between') && strcmp(test,'signrank')
    warning(['test=''signrank'' is a paired-samples test and is incompatible ' ...
        'with design=''between''. Using test=''ranksum'' instead.']);
    test = 'ranksum';
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
