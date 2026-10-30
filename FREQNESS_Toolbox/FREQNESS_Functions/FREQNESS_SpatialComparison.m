function STATS = FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI,varargin)

% Preserve the established two-structure interface while accepting one
% multidimensional cell array for generalized factorial designs.
if iscell(FREQ_A)
    if nargin<2
        error('A cell-array input must be followed by MNI coordinates.');
    end
    if nargin>=3
        factorial_options = [{MNI},varargin];
    else
        factorial_options = {};
    end
    STATS = factorial_spatial_comparison(FREQ_A,FREQ_B, ...
        factorial_options{:});
    return
elseif nargin<3
    error('FREQ_A, FREQ_B, and MNI coordinates are required.');
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
%  This function performs pairwise or factorial spatial comparisons of
%  FREQ-NESS activation patterns. Voxel-wise statistics are combined with a
%  non-parametric cluster-mass permutation test based on MNI adjacency.
%  The maximum cluster mass is calculated across all requested frequencies
%  and components in every permutation, controlling the family-wise error
%  rate across the complete analysis requested in one call.
%
%  For within-subject designs, participants must have the same order in
%  FREQ_A and FREQ_B. Between-group designs may have unequal sample sizes.
%  By default, absolute activation patterns are L2-normalized within
%  participant before comparison, making the analysis invariant to
%  arbitrary GED signs and global pattern scale.
%
%  Generalized designs use one cell array whose dimensions represent the
%  experimental factors and whose entries contain the FREQ structure for
%  every combination of factor levels. For example, a 2-by-3 cell array
%  represents two factors with two and three levels. Each factor can be
%  declared as within- or between-subject. Existing calls using FREQ_A and
%  FREQ_B remain unchanged.
%
%  Complete all-within ANOVA designs use orthonormal participant contrasts
%  for each effect, with its own participant-by-effect error term. If an
%  effect has q contrasts and n participants, nominal df are [q,q*(n-1)].
%  One sign per participant is shared across all contrasts, effects,
%  frequencies and components in a draw. Validity requires independent
%  participants and joint central symmetry of the null contrast fields
%  after map preparation. This is not arbitrary within-subject shuffling.
%  Greenhouse-Geisser p values form F clusters; epsilon is recalculated
%  in every draw. Cluster p values use the common maximum across the whole
%  requested family. Omnibus F clusters do not establish a direction.
%  Mixed designs retain legacy residual permutations and warn because
%  within-block exchangeability is a stronger covariance assumption.
% ========================================================================

% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - FREQ_A, FREQ_B : structures containing
%           • pats -> [nVoxels x nComp x nFrex x nSubs]
%           • frex -> analyzed frequencies in Hz
%        Alternatively, the first input can be an N-dimensional cell array
%        containing one FREQ structure per factorial cell, followed directly
%        by MNI as the second input.
%
%  - MNI : MNI coordinates [nVoxels x 3], ordered as FREQ.pats.
%
%  - frex2compare : mandatory name-value input. Numeric values are mapped to
%                   the closest analyzed bins. Use 'all' for every bin.
%
%  - comps2compare : one-based component numbers. Default: 1.
%
%  - design : 'within' (default) or 'between'.
%             For factorial cell inputs, supply one value per factor, for
%             example {'between','within'}. A single value is applied to all
%             factors.
%
%  - test : 'ttest' (default), 'signrank', or 'ranksum'. An incompatible
%           non-parametric choice is corrected with a warning: ranksum is
%           replaced by signrank for within-subject designs, and signrank is
%           replaced by ranksum for between-group designs.
%           Factorial cell inputs additionally support 'auto', 'anova',
%           'friedman', 'kruskalwallis', and 'permutation_anova'. 'auto'
%           preserves the pairwise t test for two levels and otherwise uses
%           factorial ANOVA with design-aware permutation inference.
%
%  - condition_labels : two labels for FREQ_A and FREQ_B. Default: {'A','B'}.
%                       Retained for established pairwise calls.
%
%  - factor_names : names of the factorial dimensions. Defaults to
%                   {'Factor1','Factor2',...}.
%
%  - level_labels : one cell array of level names per factor. Defaults to
%                   {'Level1','Level2',...} within each factor.
%
%  - subject_ids : optional cell array matching the factorial input. Each
%                  entry contains participant identifiers for that FREQ
%                  structure. Repeated observations are reordered and
%                  matched by identifier when supplied. Otherwise participant
%                  order is assumed for within-subject factors with a warning.
%
%  - effects : 'all' (default) or a cell array selecting main effects and
%              interactions by factor name, for example
%              {'Task','Group_x_Task'}.
%
%  - map_mode : 'magnitude' (default) or 'signed_aligned'.
%
%  - normalization : 'l2' (default) or 'none'.
%
%  - cluster_forming_alpha : voxel p threshold. Default: .001. Pairwise
%                            tests are two-sided; omnibus F is upper-tail.
%
%  - alpha : corrected cluster significance threshold. Default: .05.
%
%  - connectivity : MNI-grid connectivity (6, 18, or 26). Default: 6.
%
%  - frex_tolerance : maximum Hz distance accepted when matching a requested
%                     frequency. By default, half the median spacing of
%                     FREQ.frex is used.
%
%  - n_permutations : 'auto' (default), 'exact', or a positive integer.
%                     'auto' enumerates every sign configuration when at
%                     most 65536 configurations are possible, and otherwise
%                     uses 10000 Monte Carlo permutations.
%                     This also applies to complete all-within ANOVA grids:
%                     14 participants -> 16384 exact configurations. Mixed
%                     designs retain Monte Carlo residual permutations.
%
%  - random_seed : random seed for Monte Carlo permutations. Default: 0.
%
%  - plot_results : logical flag controlling result figures. Default: true.
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - STATS : structure containing resolved frequencies and components,
%            mean differences, test statistics, effect sizes, uncorrected
%            voxel-wise p values, cluster-corrected p values, significant
%            masks, cluster labels, and cluster summaries.
%            For all-within ANOVA, STATS.effects.(effect) additionally holds
%            ss_effect, ss_error, epsilon_gg, df1_gg, df2_gg and the cell-by-q
%            contrast_matrix. df1/df2 remain nominal. p_uncorrected is the
%            GG F approximation, NOT a permutation p value. All map fields
%            are [nVoxels x nSelectedFrex x nSelectedComps].
%            STATS.contrast_cell_levels lists one-based factor levels for
%            each contrast-matrix row; inference_method, permutation_unit,
%            permutation_assumption and n_participants describe inference.
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

%% Map inputs from FREQ

opts = struct('frex2compare', [], 'comps2compare', [], ...
    'design', 'within', 'test', [], ...
    'condition_labels', {{'A','B'}}, 'map_mode', 'magnitude', ...
    'normalization', 'l2', 'cluster_forming_alpha', .001, ...
    'alpha', .05, 'connectivity', 6, 'frex_tolerance', [], ...
    'correction_scope', 'all_maps', 'n_permutations', 'auto', ...
    'random_seed', 0, 'plot_results', true);
opts = parse_name_value_pairs(opts,varargin{:});
[design,test] = validate_analysis_options(opts.design,opts.test);

if isempty(opts.frex2compare)
    error(['frex2compare is a mandatory name-value input. Supply numeric ' ...
        'values in Hz or use ''all''.']);
end

[patterns_a,frex_a] = validate_freq(FREQ_A,'FREQ_A');
[patterns_b,frex_b] = validate_freq(FREQ_B,'FREQ_B');
if ~isequal(size(patterns_a,1),size(patterns_b,1)) || ...
        ~isequal(size(patterns_a,2),size(patterns_b,2)) || ...
        ~isequal(size(patterns_a,3),size(patterns_b,3))
    error(['FREQ_A.pats and FREQ_B.pats must have matching voxel, ' ...
        'component, and frequency dimensions.']);
end
if numel(frex_a)~=numel(frex_b) || any(abs(frex_a-frex_b)>1e-10)
    error('FREQ_A.frex and FREQ_B.frex must contain the same grid.');
end

[nvoxs,nComp,~,n_a] = size(patterns_a);
n_b = size(patterns_b,4);
if strcmp(design,'within') && n_a~=n_b
    error('Within-subject comparisons require matching participant counts.');
end
if n_a<2 || n_b<2
    error('Each sample must contain at least two participants.');
end
if ~isnumeric(MNI) || ~isreal(MNI) || any(~isfinite(MNI(:))) || ...
        ~isequal(size(MNI),[nvoxs 3])
    error('MNI must be [nVox x 3], finite, real, and match FREQ.pats.');
end

condition_labels = validate_condition_labels(opts.condition_labels);
map_mode = validatestring(opts.map_mode,{'magnitude','signed_aligned'});
normalization = validatestring(opts.normalization,{'l2','none'});
if ~strcmpi(opts.correction_scope,'all_maps')
    error('correction_scope currently supports only ''all_maps''.');
end
validate_probability(opts.cluster_forming_alpha,'cluster_forming_alpha');
validate_probability(opts.alpha,'alpha');
if ~isnumeric(opts.connectivity) || ~isscalar(opts.connectivity) || ...
        ~ismember(opts.connectivity,[6 18 26])
    error('connectivity must be 6, 18, or 26.');
end
if ~isnumeric(opts.random_seed) || ~isscalar(opts.random_seed) || ...
        ~isfinite(opts.random_seed) || opts.random_seed~=round(opts.random_seed)
    error('random_seed must be one finite integer.');
end
if ~islogical(opts.plot_results) || ~isscalar(opts.plot_results)
    error('plot_results must be one logical value.');
end

[requested_frex,frex2compare,idx_frex2compare] = ...
    select_frequencies(frex_a,opts.frex2compare,opts.frex_tolerance);
comps2compare = select_components(opts.comps2compare,nComp);

% Reorder selected maps as [nVox x nSelectedFrex x nSelectedComp x nSubs]
maps_a = permute(patterns_a(:,comps2compare,idx_frex2compare,:),[1 3 2 4]);
maps_b = permute(patterns_b(:,comps2compare,idx_frex2compare,:),[1 3 2 4]);
[maps_a,maps_b,sign_flips_a,sign_flips_b] = ...
    prepare_maps(maps_a,maps_b,map_mode,normalization,design);

fprintf(['\nFREQNESS Spatial Comparison: comparing %s and %s at %d frequencies ' ...
    'and %d components using a %s-subject %s analysis (n=%d, n=%d).\n'], ...
    condition_labels{1},condition_labels{2},numel(frex2compare), ...
    numel(comps2compare),design,test,n_a,n_b);

%% Compute observed statistics

if strcmp(design,'within')
    differences = maps_b-maps_a;
    [mean_difference,statistic,effect_size,p_uncorrected, ...
        statistic_name,effect_size_name] = ...
        within_statistics(differences,test);
else
    [mean_difference,statistic,effect_size,p_uncorrected, ...
        statistic_name,effect_size_name] = ...
        between_statistics(maps_a,maps_b,test);
end
[adjacency,grid_spacing] = spatial_adjacency(MNI,opts.connectivity);

%% Compute maximum-cluster permutation distribution

if strcmp(design,'within')
    [permutation_plan,exact_permutations] = ...
        permutation_plan_within(n_a,opts.n_permutations,opts.random_seed);
    permutation_description = 'paired sign-flip';
else
    [permutation_plan,exact_permutations] = ...
        permutation_plan_between(n_a,n_b,opts.n_permutations,opts.random_seed);
    permutation_description = 'independent label';
    pooled_maps = cat(4,maps_a,maps_b);
end
n_permutations = size(permutation_plan,1);
fprintf('Computing %d %s %s permutations...\n',n_permutations, ...
    ternary(exact_permutations,'exact','Monte Carlo'), ...
    permutation_description);
null_maxima = zeros(n_permutations,1);
progress_step = max(1,floor(n_permutations/10));
for permi = 1:n_permutations
    if strcmp(design,'within')
        signs = reshape(permutation_plan(permi,:),[1 1 1 n_a]);
        [~,permuted_statistic,~,permuted_p] = ...
            within_statistics(differences.*signs,test);
    else
        group_a = permutation_plan(permi,:);
        [~,permuted_statistic,~,permuted_p] = ...
            between_statistics(pooled_maps(:,:,:,group_a), ...
            pooled_maps(:,:,:,~group_a),test);
    end
    null_maxima(permi) = maximum_cluster_mass(permuted_statistic, ...
        permuted_p,opts.cluster_forming_alpha,adjacency);
    if mod(permi,progress_step)==0 || permi==n_permutations
        fprintf('  Permutation %d/%d\n',permi,n_permutations);
    end
end

%% Correct observed spatial clusters

p_cluster_corrected = ones(size(statistic));
significant_mask = false(size(statistic));
cluster_labels = zeros(size(statistic));
clusters = empty_cluster_structure;

for frexi = 1:numel(frex2compare)
    for compi = 1:numel(comps2compare)
        spatial_map = statistic(:,frexi,compi);
        supra_threshold = p_uncorrected(:,frexi,compi)<= ...
            opts.cluster_forming_alpha;
        positive_label = 0;
        negative_label = 0;
        for polarity = [1 -1]
            if polarity>0
                mask = supra_threshold & spatial_map>0;
            else
                mask = supra_threshold & spatial_map<0;
            end
            observed_clusters = connected_components(mask,adjacency);
            for clusteri = 1:numel(observed_clusters)
                members = observed_clusters{clusteri};
                mass = sum(abs(spatial_map(members)));
                if exact_permutations
                    p_corrected = mean(null_maxima>=mass);
                else
                    p_corrected = (1+sum(null_maxima>=mass))/(n_permutations+1);
                end
                p_cluster_corrected(members,frexi,compi) = p_corrected;
                significant_mask(members,frexi,compi) = p_corrected<=opts.alpha;
                if polarity>0
                    positive_label = positive_label+1;
                    cluster_label = positive_label;
                else
                    negative_label = negative_label+1;
                    cluster_label = -negative_label;
                end
                cluster_labels(members,frexi,compi) = cluster_label;
                [~,peak_local] = max(abs(spatial_map(members)));
                peak_voxel = members(peak_local);
                record = struct( ...
                    'frex',frex2compare(frexi), ...
                    'comp',comps2compare(compi), ...
                    'polarity',polarity, ...
                    'size',numel(members), ...
                    'mass',mass, ...
                    'p_corrected',p_corrected, ...
                    'peak_t',spatial_map(peak_voxel), ...
                    'peak_statistic',spatial_map(peak_voxel), ...
                    'peak_voxel',peak_voxel, ...
                    'peak_mni',MNI(peak_voxel,:));
                clusters(end+1) = record; %#ok<AGROW>
            end
        end
    end
end

%% Assemble output structure

STATS.requested_frex = requested_frex;
STATS.frex = frex2compare;
STATS.idx_frex = idx_frex2compare;
STATS.comps = comps2compare;
STATS.condition_labels = condition_labels;
STATS.design = design;
STATS.test = test;
if strcmp(design,'within')
    STATS.n_participants = n_a;
else
    STATS.n_participants = [n_a n_b];
end
STATS.n_A = n_a;
STATS.n_B = n_b;
STATS.mean_difference = mean_difference;
STATS.statistic = statistic;
STATS.statistic_name = statistic_name;
if strcmp(test,'ttest')
    STATS.t_values = statistic;
else
    STATS.t_values = nan(size(statistic));
end
STATS.effect_size = effect_size;
STATS.effect_size_name = effect_size_name;
if strcmp(design,'within') && strcmp(test,'ttest')
    STATS.effect_size_dz = effect_size;
else
    STATS.effect_size_dz = nan(size(effect_size));
end
STATS.p_uncorrected = p_uncorrected;
STATS.p_cluster_corrected = p_cluster_corrected;
STATS.significant_mask = significant_mask;
STATS.cluster_labels = cluster_labels;
STATS.clusters = clusters;
STATS.map_mode = map_mode;
STATS.normalization = normalization;
STATS.cluster_forming_alpha = opts.cluster_forming_alpha;
STATS.alpha = opts.alpha;
STATS.connectivity = opts.connectivity;
STATS.grid_spacing = grid_spacing;
STATS.correction_scope = 'all_maps';
STATS.n_permutations = n_permutations;
STATS.exact_permutations = exact_permutations;
STATS.random_seed = opts.random_seed;
STATS.sign_flips_a = sign_flips_a;
STATS.sign_flips_b = sign_flips_b;

if opts.plot_results
    plot_spatial_results(MNI,STATS);
end

end

%% Factorial Spatial Comparison
function STATS = factorial_spatial_comparison(FREQ_GRID,MNI,varargin)

opts = struct('frex2compare', [], 'comps2compare', [], ...
    'design', 'within', 'test', 'auto', 'condition_labels', [], ...
    'factor_names', [], 'level_labels', [], 'subject_ids', [], ...
    'effects', 'all', 'map_mode', 'magnitude', ...
    'normalization', 'l2', 'cluster_forming_alpha', .001, ...
    'alpha', .05, 'connectivity', 6, 'frex_tolerance', [], ...
    'correction_scope', 'all_maps', 'n_permutations', 'auto', ...
    'random_seed', 0, 'plot_results', true);
opts = parse_name_value_pairs(opts,varargin{:});

if isempty(opts.frex2compare)
    error(['frex2compare is a mandatory name-value input. Supply numeric ' ...
        'values in Hz or use ''all''.']);
end

[FREQ_GRID,factor_sizes,cell_levels] = normalize_factorial_grid(FREQ_GRID);
nFactors = numel(factor_sizes);
factor_names = validate_factor_names(opts.factor_names,nFactors);
level_labels = validate_level_labels(opts.level_labels, ...
    opts.condition_labels,factor_names,factor_sizes);
factor_design = validate_factor_design(opts.design,nFactors);
[test,test_requested] = validate_factorial_test(opts.test, ...
    factor_design,factor_sizes);
effect_definitions = factorial_effect_definitions(factor_names,opts.effects);

map_mode = validatestring(opts.map_mode,{'magnitude','signed_aligned'});
normalization = validatestring(opts.normalization,{'l2','none'});
if ~any(strcmpi(opts.correction_scope,{'all_maps','all_effects_and_maps'}))
    error(['correction_scope supports ''all_maps'' or ' ...
        '''all_effects_and_maps'' for factorial inputs.']);
end
validate_probability(opts.cluster_forming_alpha,'cluster_forming_alpha');
validate_probability(opts.alpha,'alpha');
if ~isnumeric(opts.connectivity) || ~isscalar(opts.connectivity) || ...
        ~ismember(opts.connectivity,[6 18 26])
    error('connectivity must be 6, 18, or 26.');
end
if ~isnumeric(opts.random_seed) || ~isscalar(opts.random_seed) || ...
        ~isfinite(opts.random_seed) || opts.random_seed~=round(opts.random_seed)
    error('random_seed must be one finite integer.');
end
if ~islogical(opts.plot_results) || ~isscalar(opts.plot_results)
    error('plot_results must be one logical value.');
end

[FREQ_GRID,subject_ids] = match_factorial_subjects(FREQ_GRID, ...
    cell_levels,factor_design,opts.subject_ids);

% Preserve the complete established implementation for generalized calls
% that still represent one two-level pairwise comparison.
pairwise_tests = {'ttest','signrank','ranksum'};
if nFactors==1 && factor_sizes(1)==2 && ismember(test,pairwise_tests)
    pairwise_options = {'frex2compare',opts.frex2compare, ...
        'comps2compare',opts.comps2compare,'design',factor_design{1}, ...
        'test',test,'condition_labels',level_labels{1}, ...
        'map_mode',map_mode,'normalization',normalization, ...
        'cluster_forming_alpha',opts.cluster_forming_alpha, ...
        'alpha',opts.alpha,'connectivity',opts.connectivity, ...
        'frex_tolerance',opts.frex_tolerance, ...
        'correction_scope',opts.correction_scope, ...
        'n_permutations',opts.n_permutations, ...
        'random_seed',opts.random_seed,'plot_results',opts.plot_results};
    STATS = FREQNESS_SpatialComparison(FREQ_GRID{1},FREQ_GRID{2},MNI, ...
        pairwise_options{:});
    STATS.analysis = 'pairwise';
    STATS.test_requested = test_requested;
    STATS.factor_names = factor_names;
    STATS.level_labels = level_labels;
    STATS.factor_design = factor_design;
    STATS.input_shape = factor_sizes;
    STATS.subject_ids = subject_ids;
    effect_key = effect_definitions(1).key;
    STATS.effects = struct;
    STATS.effects.(effect_key) = pairwise_effect_result(STATS, ...
        effect_definitions(1));
    return
end

[patterns,frex] = validate_factorial_freq(FREQ_GRID);
[nvoxs,nComp,~,~] = size(patterns{1});
if ~isnumeric(MNI) || ~isreal(MNI) || any(~isfinite(MNI(:))) || ...
        ~isequal(size(MNI),[nvoxs 3])
    error('MNI must be [nVox x 3], finite, real, and match FREQ.pats.');
end
[requested_frex,frex2compare,idx_frex2compare] = ...
    select_frequencies(frex,opts.frex2compare,opts.frex_tolerance);
comps2compare = select_components(opts.comps2compare,nComp);

[maps,sign_flips] = prepare_factorial_maps(patterns,comps2compare, ...
    idx_frex2compare,map_mode,normalization);
[Y,observation_levels,subject_index,cell_means,cell_n] = ...
    factorial_long_data(maps,cell_levels,factor_design);
[adjacency,grid_spacing] = spatial_adjacency(MNI,opts.connectivity);

fprintf(['\nFREQNESS Spatial Comparison: testing a %s design with %d factors, ' ...
    '%d factorial cells, %d frequencies, and %d components using %s.\n'], ...
    strjoin(factor_design,'/'),nFactors,numel(FREQ_GRID), ...
    numel(frex2compare),numel(comps2compare),test);

repeated_anova = all(strcmp(factor_design,'within')) && ...
    ismember(test,{'anova','permutation_anova'});
if repeated_anova
    [effect_results,permutation_models] = repeated_spatial_statistics(Y, ...
        observation_levels,subject_index,factor_sizes,effect_definitions, ...
        nvoxs,numel(frex2compare),numel(comps2compare));
elseif ismember(test,{'friedman','kruskalwallis'})
    [effect_results,permutation_models] = factorial_rank_statistics(Y, ...
        observation_levels,subject_index,test, ...
        effect_definitions,nvoxs,numel(frex2compare),numel(comps2compare));
else
    if any(strcmp(factor_design,'within'))
        warning('FREQNESS:SpatialComparison:LegacyMixedInference', ...
            ['Mixed-design spatial ANOVA retains legacy residual permutations ' ...
            'and requires within-block exchangeability. Participant-contrast ' ...
            'inference is currently available only for all-within designs.']);
    end
    [effect_results,permutation_models] = factorial_anova_statistics(Y, ...
        observation_levels,subject_index,factor_sizes,factor_design, ...
        effect_definitions,nvoxs,numel(frex2compare),numel(comps2compare));
end

if repeated_anova
    [sign_plan,exact_permutations] = permutation_plan_within( ...
        cell_n(1),opts.n_permutations,opts.random_seed);
    n_permutations = size(sign_plan,1);
    fprintf('Computing %d participant sign flips (exact: %d)...\n', ...
        n_permutations,exact_permutations);
else
    [n_permutations,exact_permutations] = factorial_permutation_count( ...
        opts.n_permutations);
    fprintf('Computing %d Monte Carlo design-aware permutations...\n', ...
        n_permutations);
end
rng_state = rng;
rng(opts.random_seed,'twister');
rng_cleanup = onCleanup(@()rng(rng_state));
null_maxima = zeros(n_permutations,1);
progress_step = max(1,floor(n_permutations/10));
for permi = 1:n_permutations
    maximum = 0;
    for effecti = 1:numel(permutation_models)
        model = permutation_models(effecti);
        if repeated_anova
            % One participant sign is shared across every effect and map.
            [permuted_statistic,permuted_p] = ...
                repeated_spatial_statistic(model,sign_plan(permi,:)');
        else
            [permuted_statistic,permuted_p] = permute_factorial_model(model);
        end
        permuted_statistic = reshape(permuted_statistic, ...
            [nvoxs numel(frex2compare) numel(comps2compare)]);
        permuted_p = reshape(permuted_p, ...
            [nvoxs numel(frex2compare) numel(comps2compare)]);
        maximum = max(maximum,maximum_positive_cluster_mass( ...
            permuted_statistic,permuted_p,opts.cluster_forming_alpha, ...
            adjacency));
    end
    null_maxima(permi) = maximum;
    if mod(permi,progress_step)==0 || permi==n_permutations
        fprintf('  Permutation %d/%d\n',permi,n_permutations);
    end
end
clear rng_cleanup

for effecti = 1:numel(effect_results)
    effect_results(effecti) = correct_factorial_clusters( ...
        effect_results(effecti),MNI,adjacency,null_maxima, ...
        opts.cluster_forming_alpha,opts.alpha,n_permutations, ...
        exact_permutations,frex2compare,comps2compare);
end

STATS.analysis = 'factorial';
STATS.requested_frex = requested_frex;
STATS.frex = frex2compare;
STATS.idx_frex = idx_frex2compare;
STATS.comps = comps2compare;
STATS.factor_names = factor_names;
STATS.level_labels = level_labels;
STATS.factor_design = factor_design;
STATS.design = factor_design;
STATS.input_shape = factor_sizes;
STATS.test = test;
STATS.test_requested = test_requested;
STATS.effects = struct;
STATS.effect_names = cell(1,numel(effect_results));
for effecti = 1:numel(effect_results)
    key = effect_results(effecti).key;
    STATS.effects.(key) = effect_results(effecti);
    STATS.effect_names{effecti} = key;
end
STATS.cell_means = reshape(cell_means,size(FREQ_GRID));
STATS.cell_n_participants = reshape(cell_n,size(FREQ_GRID));
STATS.subject_ids = subject_ids;
STATS.map_mode = map_mode;
STATS.normalization = normalization;
STATS.cluster_forming_alpha = opts.cluster_forming_alpha;
STATS.alpha = opts.alpha;
STATS.connectivity = opts.connectivity;
STATS.grid_spacing = grid_spacing;
STATS.correction_scope = 'all_effects_and_maps';
STATS.n_permutations = n_permutations;
STATS.exact_permutations = exact_permutations;
STATS.random_seed = opts.random_seed;
STATS.sign_flips = reshape(sign_flips,size(FREQ_GRID));
if repeated_anova
    STATS.inference_method = 'within_subject_contrast_signflip';
    STATS.permutation_unit = 'participant';
    STATS.permutation_assumption = ['independent participants; joint central ' ...
        'symmetry of null participant contrast fields'];
    STATS.n_participants = cell_n(1);
    STATS.contrast_cell_levels = cell_levels;
else
    STATS.inference_method = 'legacy_factorial';
end

if opts.plot_results
    plot_factorial_spatial_results(MNI,STATS);
end

end

%% Helper Function: Normalize Factorial Grid
function [FREQ_GRID,factor_sizes,cell_levels] = ...
    normalize_factorial_grid(FREQ_GRID)

if ~iscell(FREQ_GRID) || isempty(FREQ_GRID) || ...
        any(cellfun(@isempty,FREQ_GRID(:)))
    error('Factorial input must be a non-empty cell array of FREQ structures.');
end
grid_size = size(FREQ_GRID);
factor_dimensions = find(grid_size>1);
if isempty(factor_dimensions)
    error('Factorial input must contain at least two cells.');
end
factor_sizes = grid_size(factor_dimensions);
if isscalar(factor_sizes)
    FREQ_GRID = reshape(FREQ_GRID,[factor_sizes 1]);
else
    FREQ_GRID = reshape(FREQ_GRID,factor_sizes);
end
nFactors = numel(factor_sizes);
cell_levels = zeros(numel(FREQ_GRID),nFactors);
grid_shape = factor_sizes;
if nFactors==1
    grid_shape = [factor_sizes 1];
end
for celli = 1:numel(FREQ_GRID)
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
function level_labels = validate_level_labels(level_labels,condition_labels, ...
    factor_names,factor_sizes)

nFactors = numel(factor_sizes);
if isempty(level_labels)
    level_labels = cell(1,nFactors);
    for factori = 1:nFactors
        level_labels{factori} = arrayfun(@(value)sprintf('Level%d',value), ...
            1:factor_sizes(factori),'UniformOutput',false);
    end
    if ~isempty(condition_labels)
        if nFactors~=1 || factor_sizes(1)~=2
            error(['condition_labels can only be used with a one-factor, ' ...
                'two-level cell input. Use level_labels instead.']);
        end
        level_labels{1} = validate_condition_labels(condition_labels);
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
            ~all(cellfun(@(value)ischar(value) && ~isempty(strtrim(value)),labels))
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
        test = 'ttest';
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

%% Helper Function: Match Repeated Participants
function [FREQ_GRID,subject_ids] = match_factorial_subjects(FREQ_GRID, ...
    cell_levels,design,subject_ids)

nCells = numel(FREQ_GRID);
nsubs = zeros(nCells,1);
for celli = 1:nCells
    [patterns,~] = validate_freq(FREQ_GRID{celli},sprintf('FREQ_GRID{%d}',celli));
    nsubs(celli) = size(patterns,4);
    if nsubs(celli)<2
        error('Every factorial cell must contain at least two participants.');
    end
end

within_factors = find(strcmp(design,'within'));
between_factors = find(strcmp(design,'between'));
if isempty(subject_ids)
    subject_ids = cell(size(FREQ_GRID));
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
    subject_ids = reshape(subject_ids,size(FREQ_GRID));
    for celli = 1:nCells
        subject_ids{celli} = normalize_subject_ids(subject_ids{celli}, ...
            nsubs(celli),celli);
    end
end

if isempty(within_factors)
    return
end
if isempty(between_factors)
    strata = ones(nCells,1);
else
    [~,~,strata] = unique(cell_levels(:,between_factors),'rows','stable');
end
for stratum = unique(strata(:))'
    cells = find(strata==stratum);
    reference = cells(1);
    reference_ids = subject_ids{reference};
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
        FREQ_GRID{celli}.pats = FREQ_GRID{celli}.pats(:,:,:,order);
        subject_ids{celli} = reference_ids;
    end
end

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

%% Helper Function: Validate Factorial FREQ Inputs
function [patterns,frex] = validate_factorial_freq(FREQ_GRID)

patterns = cell(size(FREQ_GRID));
frex = [];
reference_size = [];
for celli = 1:numel(FREQ_GRID)
    [patterns{celli},current_frex] = validate_freq(FREQ_GRID{celli}, ...
        sprintf('FREQ_GRID{%d}',celli));
    current_size = size(patterns{celli});
    current_size(end+1:4) = 1;
    if celli==1
        frex = current_frex;
        reference_size = current_size(1:3);
    elseif any(current_size(1:3)~=reference_size) || ...
            numel(current_frex)~=numel(frex) || ...
            any(abs(current_frex-frex)>1e-10)
        error(['All factorial FREQ inputs must share voxel, component, and ' ...
            'frequency dimensions and the same FREQ.frex grid.']);
    end
end

end

%% Helper Function: Prepare Factorial Maps
function [maps,sign_flips] = prepare_factorial_maps(patterns,components, ...
    frequency_indices,map_mode,normalization)

maps = cell(size(patterns));
sign_flips = cell(size(patterns));
for celli = 1:numel(patterns)
    maps{celli} = permute(patterns{celli}(:,components, ...
        frequency_indices,:),[1 3 2 4]);
    sign_flips{celli} = ones(size(maps{celli},2),size(maps{celli},3), ...
        size(maps{celli},4));
end

if strcmp(map_mode,'magnitude')
    maps = cellfun(@abs,maps,'UniformOutput',false);
else
    reference = maps{1}(:,:,:,1);
    for frexi = 1:size(reference,2)
        for compi = 1:size(reference,3)
            if norm(reference(:,frexi,compi))<=eps
                reference(:,frexi,compi) = 1;
            end
            for celli = 1:numel(maps)
                for subi = 1:size(maps{celli},4)
                    if maps{celli}(:,frexi,compi,subi)'* ...
                            reference(:,frexi,compi)<0
                        maps{celli}(:,frexi,compi,subi) = ...
                            -maps{celli}(:,frexi,compi,subi);
                        sign_flips{celli}(frexi,compi,subi) = -1;
                    end
                end
            end
        end
    end
end

if strcmp(normalization,'l2')
    for celli = 1:numel(maps)
        scales = sqrt(sum(maps{celli}.^2,1));
        if any(scales(:)<=eps)
            error('A selected factorial pattern has zero L2 norm.');
        end
        maps{celli} = maps{celli}./scales;
    end
end

end

%% Helper Function: Assemble Long Factorial Data
function [Y,levels,subject_index,cell_means,cell_n] = ...
    factorial_long_data(maps,cell_levels,design)

nCells = numel(maps);
nFactors = size(cell_levels,2);
within_factors = find(strcmp(design,'within'));
between_factors = find(strcmp(design,'between'));
cell_n = cellfun(@(value)size(value,4),maps(:));
cell_means = cell(size(maps));
Y = [];
levels = zeros(sum(cell_n),nFactors);
subject_index = zeros(sum(cell_n),1);
row_start = 1;

if isempty(within_factors)
    strata = (1:nCells)';
else
    if isempty(between_factors)
        strata = ones(nCells,1);
    else
        [~,~,strata] = unique(cell_levels(:,between_factors),'rows','stable');
    end
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
else
    next_subject = 0;
end

for celli = 1:nCells
    current = maps{celli};
    cell_means{celli} = mean(current,4);
    nsubs = size(current,4);
    rows = row_start:(row_start+nsubs-1);
    current_Y = reshape(permute(current,[4 1 2 3]),nsubs,[]);
    Y = [Y;current_Y]; %#ok<AGROW>
    levels(rows,:) = repmat(cell_levels(celli,:),nsubs,1);
    if isempty(within_factors)
        subject_index(rows) = next_subject+(1:nsubs);
        next_subject = next_subject+nsubs;
    else
        subject_index(rows) = stratum_offsets(strata(celli))+(1:nsubs);
    end
    row_start = row_start+nsubs;
end

end

%% Helper Function: Convert Pairwise Effect
function effect = pairwise_effect_result(STATS,definition)

effect.name = definition.name;
effect.key = definition.key;
effect.factors = definition.factors;
effect.statistic = STATS.statistic;
effect.statistic_name = STATS.statistic_name;
effect.effect_size = STATS.effect_size;
effect.effect_size_name = STATS.effect_size_name;
effect.p_uncorrected = STATS.p_uncorrected;
effect.p_cluster_corrected = STATS.p_cluster_corrected;
effect.significant_mask = STATS.significant_mask;
effect.cluster_labels = STATS.cluster_labels;
effect.clusters = STATS.clusters;

end

%% Helper Function: Complete Within-Subject ANOVA
function [results,models] = repeated_spatial_statistics(Y,levels, ...
    subject_index,factor_sizes,effects,nvoxs,nfrex,ncomps)

% Y is [observations x maps]; maps flatten voxel, frequency, component.
% Original cell maps have already been rectified/aligned and normalized.
% Do not renormalize the averages or contrasts formed here.
subjects = unique(subject_index,'stable');
nsubs = numel(subjects);
ncells = prod(factor_sizes);
nmaps = size(Y,2);
cell_levels = levels(subject_index==subjects(1),:);
data = zeros(nsubs,ncells,nmaps);
for subi = 1:nsubs
    rows = subject_index==subjects(subi);
    if sum(rows)~=ncells || ~isequal(levels(rows,:),cell_levels)
        error('Within-subject ANOVA requires the same complete cells per participant.');
    end
    data(subi,:,:) = reshape(Y(rows,:),[1 ncells nmaps]);
end
% Center before projecting to avoid leakage from large subject intercepts.
data = data-mean(data,2);
% Suppress only projection roundoff, so an absent effect does not acquire
% evidence from division by its equally tiny estimated error variance.
projection_tolerance = 16*eps*ncells*reshape( ...
    max(max(abs(data),[],1),[],2),[1 nmaps]);
results = repmat(empty_factorial_effect,1,numel(effects));
models = repmat(struct('scores',[],'second_moment',[], ...
    'contrast_matrix',[]),1,numel(effects));
for effecti = 1:numel(effects)
    matrices = cell(1,numel(factor_sizes));
    for factori = 1:numel(factor_sizes)
        nlevels = factor_sizes(factori);
        if ismember(factori,effects(effecti).factors)
            coding = zeros(nlevels,nlevels-1);
            for leveli = 1:(nlevels-1)
                scale = sqrt(leveli*(leveli+1));
                coding(1:leveli,leveli) = -1/scale;
                coding(leveli+1,leveli) = leveli/scale;
            end
        else
            coding = ones(nlevels,1)/sqrt(nlevels);
        end
        matrices{factori} = coding(cell_levels(:,factori),:);
    end
    weights = interaction_matrix(matrices); % [cells x orthonormal contrasts]
    ncontrasts = size(weights,2);
    scores = zeros(nsubs,ncontrasts,nmaps);
    for contrasti = 1:ncontrasts
        values = reshape(sum(data.*reshape(weights(:,contrasti), ...
            [1 ncells 1]),2),[nsubs nmaps]);
        values(abs(values)<=projection_tolerance) = 0;
        scores(:,contrasti,:) = reshape(values,[nsubs 1 nmaps]);
    end
    second_moment = zeros(ncontrasts,ncontrasts,nmaps);
    for lefti = 1:ncontrasts
        for righti = 1:ncontrasts
            second_moment(lefti,righti,:) = sum( ...
                scores(:,lefti,:).*scores(:,righti,:),1);
        end
    end
    models(effecti) = struct('scores',scores, ...
        'second_moment',second_moment,'contrast_matrix',weights);
    [statistic,p_values,eta,epsilon_gg,ss_effect,ss_error] = ...
        repeated_spatial_statistic(models(effecti),ones(nsubs,1));
    df1 = ncontrasts;
    df2 = ncontrasts*(nsubs-1);
    result = initialize_factorial_effect(effects(effecti),statistic, ...
        p_values,eta,'F','partial eta squared',df1,df2,nvoxs,nfrex,ncomps);
    result.epsilon_gg = reshape(epsilon_gg,[nvoxs nfrex ncomps]);
    result.df1_gg = df1*result.epsilon_gg;
    result.df2_gg = df2*result.epsilon_gg;
    result.ss_effect = reshape(ss_effect,[nvoxs nfrex ncomps]);
    result.ss_error = reshape(ss_error,[nvoxs nfrex ncomps]);
    result.contrast_matrix = weights;
    result.p_uncorrected_method = 'Greenhouse-Geisser F approximation';
    results(effecti) = result;
end

end

function [statistic,p_values,eta,epsilon_gg,ss_effect,ss_error] = ...
    repeated_spatial_statistic(model,signs)

% Z: [participants x contrasts x maps]. Each effect has its own
% participant-by-effect error; F = (n-1)*SS_effect/SS_error, df=(q,q*(n-1)).
% A common sign across contrasts retains arbitrary within-effect covariance.
scores = model.scores;
[nsubs,ncontrasts,nmaps] = size(scores);
means = reshape(sum(scores.*reshape(signs,[nsubs 1 1]),1), ...
    [ncontrasts nmaps])/nsubs;
scatter = model.second_moment-nsubs* ...
    reshape(means,[ncontrasts 1 nmaps]).*reshape(means,[1 ncontrasts nmaps]);
diagonal = (1:ncontrasts+1:ncontrasts^2)';
scatter_flat = reshape(scatter,[ncontrasts^2 nmaps]);
moment_flat = reshape(model.second_moment,[ncontrasts^2 nmaps]);
ss_error = sum(scatter_flat(diagonal,:),1);
total = sum(moment_flat(diagonal,:),1);
% Recompute near-constant fields directly rather than subtracting moments.
for mapi = find(ss_error<=64*eps*total)
    values = reshape(scores(:,:,mapi),[nsubs ncontrasts]).*signs;
    centered = values-means(:,mapi)';
    if all(max(values,[],1)==min(values,[],1))
        centered(:) = 0;
    end
    current_scatter = centered'*centered;
    scatter_flat(:,mapi) = current_scatter(:);
end
ss_error = sum(scatter_flat(diagonal,:),1);
ss_effect = nsubs*sum(means.^2,1);
statistic = zeros(1,nmaps);
positive_error = ss_error>0;
statistic(positive_error) = (nsubs-1)*ss_effect(positive_error)./ ...
    ss_error(positive_error);
statistic(ss_error==0 & ss_effect>0) = Inf;
eta = zeros(1,nmaps);
positive_total = ss_effect+ss_error>0;
eta(positive_total) = ss_effect(positive_total)./ ...
    (ss_effect(positive_total)+ss_error(positive_total));
epsilon_gg = ones(1,nmaps);
if ncontrasts>1
    frobenius = sum(scatter_flat.^2,1);
    valid = frobenius>0;
    epsilon_gg(valid) = ss_error(valid).^2./(ncontrasts*frobenius(valid));
    epsilon_gg = max(1/ncontrasts,min(1,epsilon_gg));
end
% GG is an approximate voxel-wise reference, also recalculated in each
% draw. Corrected cluster p values come from the participant sign flips.
p_values = fcdf(statistic,ncontrasts*epsilon_gg, ...
    ncontrasts*(nsubs-1)*epsilon_gg,'upper');

end

%% Helper Function: Factorial ANOVA Statistics
function [results,models] = factorial_anova_statistics(Y,levels, ...
    subject_index,factor_sizes,design,effects,nvoxs,nfrex,ncomps)

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
models = repmat(empty_permutation_model,1,numel(effects));
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
        current_levels = levels; %#ok<NASGU>
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
    [statistic,p_values,effect_size,model] = partial_f_model( ...
        current_Y,X_full,X_reduced,permutation_mode,block_index);
    results(effecti) = initialize_factorial_effect(effects(effecti), ...
        statistic,p_values,effect_size,'F','partial eta squared', ...
        model.df1,model.df2,nvoxs,nfrex,ncomps);
    model.key = effects(effecti).key;
    models(effecti) = model;
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
Y_subject = zeros(numel(subjects),size(Y,2));
levels_subject = zeros(numel(subjects),size(levels,2));
for subi = 1:numel(subjects)
    rows = subject_index==subjects(subi);
    Y_subject(subi,:) = mean(Y(rows,:),1);
    levels_subject(subi,:) = levels(find(rows,1),:);
end

end

function [statistic,p_values,effect_size,model] = partial_f_model( ...
    Y,X_full,X_reduced,permutation_mode,block_index)

tolerance = max(size(X_full))*eps(norm(X_full));
rank_full = rank(X_full,tolerance);
rank_reduced = rank(X_reduced,tolerance);
df1 = rank_full-rank_reduced;
df2 = size(X_full,1)-rank_full;
if df1<1 || df2<1
    error(['The factorial design does not contain enough independent ' ...
        'observations to estimate the requested effect.']);
end
projection_full = eye(size(X_full,1))-X_full*pinv(X_full);
projection_reduced = eye(size(X_reduced,1))-X_reduced*pinv(X_reduced);
residual_full = projection_full*Y;
residual_reduced = projection_reduced*Y;
sse_full = sum(residual_full.^2,1);
sse_reduced = sum(residual_reduced.^2,1);
numerator = max(sse_reduced-sse_full,0)/df1;
denominator = sse_full/df2;
statistic = zeros(1,size(Y,2));
valid = denominator>eps;
statistic(valid) = numerator(valid)./denominator(valid);
statistic(~valid & numerator>eps) = inf;
p_values = 1-fcdf(statistic,df1,df2);
effect_size = zeros(size(statistic));
denominator_eta = statistic*df1+df2;
valid_eta = denominator_eta>0;
effect_size(valid_eta) = statistic(valid_eta)*df1./ ...
    denominator_eta(valid_eta);

model = empty_permutation_model;
model.kind = 'anova';
model.Y_fitted_reduced = Y-residual_reduced;
model.residuals_reduced = residual_reduced;
model.projection_full = projection_full;
model.projection_reduced = projection_reduced;
model.permutation_mode = permutation_mode;
model.block_index = block_index;
model.df1 = df1;
model.df2 = df2;

end

%% Helper Function: Factorial Rank Statistics
function [results,models] = factorial_rank_statistics(Y,levels, ...
    subject_index,test,effects,nvoxs,nfrex,ncomps)

if size(levels,2)~=1 || numel(effects)~=1
    error('friedman and kruskalwallis support one factorial effect only.');
end
nlevels = max(levels(:,1));
definition = effects(1);
model = empty_permutation_model;

if strcmp(test,'kruskalwallis')
    [statistic,p_values,effect_size] = kruskalwallis_columns( ...
        Y,levels(:,1),nlevels);
    model.kind = 'kruskalwallis';
    model.Y = Y;
    model.group_levels = levels(:,1);
    model.nlevels = nlevels;
    model.df1 = nlevels-1;
    model.df2 = inf;
    effect_name = 'epsilon squared';
else
    subjects = unique(subject_index,'stable');
    repeated = zeros(numel(subjects),nlevels,size(Y,2));
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
    [statistic,p_values,effect_size] = friedman_columns(repeated);
    model.kind = 'friedman';
    model.Y = repeated;
    model.nlevels = nlevels;
    model.df1 = nlevels-1;
    model.df2 = inf;
    effect_name = 'Kendall''s W';
end

results = initialize_factorial_effect(definition,statistic,p_values, ...
    effect_size,ternary(strcmp(test,'friedman'),'chi-square','H'), ...
    effect_name,model.df1,model.df2,nvoxs,nfrex,ncomps);
model.key = definition.key;
models = model;

end

function [statistic,p_values,effect_size] = ...
    kruskalwallis_columns(Y,groups,nlevels)

[ranks,tie_sum] = rank_columns(Y);
nobservations = size(Y,1);
rank_sums = zeros(nlevels,size(Y,2));
group_n = zeros(nlevels,1);
for leveli = 1:nlevels
    rows = groups==leveli;
    group_n(leveli) = sum(rows);
    rank_sums(leveli,:) = sum(ranks(rows,:),1);
end
statistic = 12/(nobservations*(nobservations+1))* ...
    sum(rank_sums.^2./group_n,1)-3*(nobservations+1);
tie_correction = 1-tie_sum/(nobservations^3-nobservations);
valid = tie_correction>eps;
statistic(valid) = statistic(valid)./tie_correction(valid);
statistic(~valid) = 0;
statistic = max(statistic,0);
p_values = 1-chi2cdf(statistic,nlevels-1);
effect_size = max((statistic-nlevels+1)/(nobservations-nlevels),0);

end

function [statistic,p_values,effect_size] = friedman_columns(Y)

[nsubjects,nlevels,nmaps] = size(Y);
rank_sums = zeros(nlevels,nmaps);
tie_sum = zeros(1,nmaps);
for subi = 1:nsubjects
    values = reshape(Y(subi,:,:),nlevels,nmaps);
    [ranks,current_ties] = rank_columns(values);
    rank_sums = rank_sums+ranks;
    tie_sum = tie_sum+current_ties;
end
statistic = 12/(nsubjects*nlevels*(nlevels+1))* ...
    sum(rank_sums.^2,1)-3*nsubjects*(nlevels+1);
tie_correction = 1-tie_sum/(nsubjects*nlevels*(nlevels^2-1));
valid = tie_correction>eps;
statistic(valid) = statistic(valid)./tie_correction(valid);
statistic(~valid) = 0;
statistic = max(statistic,0);
p_values = 1-chi2cdf(statistic,nlevels-1);
effect_size = statistic/(nsubjects*(nlevels-1));

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
    effect_size,statistic_name,effect_size_name,df1,df2,nvoxs,nfrex,ncomps)

result = empty_factorial_effect;
result.name = definition.name;
result.key = definition.key;
result.factors = definition.factors;
result.statistic = reshape(statistic,[nvoxs nfrex ncomps]);
result.statistic_name = statistic_name;
result.effect_size = reshape(effect_size,[nvoxs nfrex ncomps]);
result.effect_size_name = effect_size_name;
result.p_uncorrected = reshape(p_values,[nvoxs nfrex ncomps]);
result.df1 = df1;
result.df2 = df2;

end

function result = empty_factorial_effect

result = struct('name','','key','','factors',[],'statistic',[], ...
    'statistic_name','','effect_size',[],'effect_size_name','', ...
    'p_uncorrected',[],'p_cluster_corrected',[], ...
    'significant_mask',[],'cluster_labels',[],'clusters',[], ...
    'df1',[],'df2',[],'epsilon_gg',[],'df1_gg',[],'df2_gg',[], ...
    'ss_effect',[],'ss_error',[],'contrast_matrix',[], ...
    'p_uncorrected_method','parametric reference');

end

function model = empty_permutation_model

model = struct('key','','kind','','Y',[],'Y_fitted_reduced',[], ...
    'residuals_reduced',[],'projection_full',[], ...
    'projection_reduced',[],'permutation_mode','', ...
    'block_index',[],'group_levels',[],'nlevels',[], ...
    'df1',[],'df2',[]);

end

%% Helper Function: Factorial Permutations
function [count,exact] = factorial_permutation_count(requested)

exact = false;
if ischar(requested) || (isstring(requested) && isscalar(requested))
    requested = char(requested);
    validatestring(requested,{'auto','exact'});
    if strcmpi(requested,'exact')
        error(['Exact enumeration is not available for generalized ' ...
            'factorial designs. Use Monte Carlo permutations.']);
    end
    count = 10000;
elseif isnumeric(requested) && isscalar(requested) && ...
        isfinite(requested) && requested>0 && requested==round(requested)
    count = requested;
else
    error('n_permutations must be a positive integer or ''auto''.');
end

end

function [statistic,p_values] = permute_factorial_model(model)

switch model.kind
    case 'anova'
        if strcmp(model.permutation_mode,'between')
            order = randperm(size(model.residuals_reduced,1));
        else
            order = (1:size(model.residuals_reduced,1))';
            subjects = unique(model.block_index,'stable');
            for subi = 1:numel(subjects)
                rows = find(model.block_index==subjects(subi));
                order(rows) = rows(randperm(numel(rows)));
            end
        end
        permuted = model.Y_fitted_reduced+model.residuals_reduced(order,:);
        sse_full = sum((model.projection_full*permuted).^2,1);
        sse_reduced = sum((model.projection_reduced*permuted).^2,1);
        numerator = max(sse_reduced-sse_full,0)/model.df1;
        denominator = sse_full/model.df2;
        statistic = zeros(1,size(permuted,2));
        valid = denominator>eps;
        statistic(valid) = numerator(valid)./denominator(valid);
        statistic(~valid & numerator>eps) = inf;
        p_values = 1-fcdf(statistic,model.df1,model.df2);
    case 'kruskalwallis'
        groups = model.group_levels(randperm(numel(model.group_levels)));
        [statistic,p_values] = kruskalwallis_columns( ...
            model.Y,groups,model.nlevels);
    case 'friedman'
        permuted = model.Y;
        for subi = 1:size(permuted,1)
            permuted(subi,:,:) = permuted(subi,randperm(model.nlevels),:);
        end
        [statistic,p_values] = friedman_columns(permuted);
    otherwise
        error('Unknown factorial permutation model.');
end

end

function maximum = maximum_positive_cluster_mass(statistic,p_values, ...
    cluster_alpha,adjacency)

maximum = 0;
for frexi = 1:size(statistic,2)
    for compi = 1:size(statistic,3)
        spatial_map = statistic(:,frexi,compi);
        mask = p_values(:,frexi,compi)<=cluster_alpha & spatial_map>0;
        clusters = connected_components(mask,adjacency);
        for clusteri = 1:numel(clusters)
            maximum = max(maximum,sum(spatial_map(clusters{clusteri})));
        end
    end
end

end

%% Helper Function: Correct Factorial Clusters
function result = correct_factorial_clusters(result,MNI,adjacency, ...
    null_maxima,cluster_alpha,alpha,n_permutations,exact,frex,components)

result.p_cluster_corrected = ones(size(result.statistic));
result.significant_mask = false(size(result.statistic));
result.cluster_labels = zeros(size(result.statistic));
result.clusters = empty_factorial_cluster_structure;
for frexi = 1:numel(frex)
    for compi = 1:numel(components)
        spatial_map = result.statistic(:,frexi,compi);
        mask = result.p_uncorrected(:,frexi,compi)<=cluster_alpha & ...
            spatial_map>0;
        clusters = connected_components(mask,adjacency);
        for clusteri = 1:numel(clusters)
            members = clusters{clusteri};
            mass = sum(spatial_map(members));
            if exact
                p_corrected = mean(null_maxima>=mass);
            else
                p_corrected = (1+sum(null_maxima>=mass))/(n_permutations+1);
            end
            result.p_cluster_corrected(members,frexi,compi) = p_corrected;
            result.significant_mask(members,frexi,compi) = p_corrected<=alpha;
            result.cluster_labels(members,frexi,compi) = clusteri;
            [~,peak_local] = max(spatial_map(members));
            peak_voxel = members(peak_local);
            record = struct('effect',result.name,'frex',frex(frexi), ...
                'comp',components(compi),'size',numel(members), ...
                'mass',mass,'p_corrected',p_corrected, ...
                'peak_f',spatial_map(peak_voxel), ...
                'peak_statistic',spatial_map(peak_voxel), ...
                'peak_voxel',peak_voxel,'peak_mni',MNI(peak_voxel,:));
            result.clusters(end+1) = record;
        end
    end
end

end

function clusters = empty_factorial_cluster_structure

clusters = struct('effect',{},'frex',{},'comp',{},'size',{}, ...
    'mass',{},'p_corrected',{},'peak_f',{},'peak_statistic',{}, ...
    'peak_voxel',{},'peak_mni',{});

end

%% Helper Function: Plot Factorial Spatial Results
function plot_factorial_spatial_results(MNI,STATS)

effect_names = STATS.effect_names;
for effecti = 1:numel(effect_names)
    result = STATS.effects.(effect_names{effecti});
    for compi = 1:numel(STATS.comps)
        figure_width = max(760,440*numel(STATS.frex)+120);
        current_figure = figure('Color','w', ...
            'Name',sprintf('Spatial effect: %s',result.name), ...
            'Position',[100 100 figure_width 520],'Renderer','painters');
        layout = tiledlayout(1,numel(STATS.frex),'TileSpacing','compact', ...
            'Padding','compact');
        color_limit = max(result.statistic(:,:,compi),[],'all');
        if ~isfinite(color_limit) || color_limit<=0
            color_limit = 1;
        end
        for frexi = 1:numel(STATS.frex)
            nexttile
            values = result.statistic(:,frexi,compi);
            scatter3(MNI(:,1),MNI(:,2),MNI(:,3),12,[.80 .83 .86], ...
                'filled','MarkerFaceAlpha',.68,'MarkerEdgeColor','none');
            hold on
            significant = result.significant_mask(:,frexi,compi);
            if any(significant)
                scatter3(MNI(significant,1),MNI(significant,2), ...
                    MNI(significant,3),58,values(significant),'filled', ...
                    'MarkerEdgeColor',[.09 .13 .17],'LineWidth',.75);
            else
                text(.5,.04,'No corrected clusters','Units','normalized', ...
                    'HorizontalAlignment','center','Color',[.40 .44 .49], ...
                    'FontSize',9);
            end
            format_publication_spatial_axis(MNI)
            clim([0 color_limit])
            title(sprintf('%g Hz',STATS.frex(frexi)), ...
                'FontSize',12,'FontWeight','bold')
        end
        colormap(current_figure,publication_sequential_colormap(256))
        current_colorbar = colorbar;
        current_colorbar.Layout.Tile = 'east';
        current_colorbar.Label.String = sprintf('%s statistic', ...
            result.statistic_name);
        current_colorbar.FontSize = 10;
        title(layout,{sprintf('Spatial comparison - %s',result.name), ...
            sprintf(['Component %d  |  Outlined markers: ' ...
            'cluster-corrected p <= %.3g'],STATS.comps(compi),STATS.alpha)}, ...
            'FontSize',15,'FontWeight','bold','Interpreter','none');
    end
end

end

%% Helper Function: Validate FREQ Input
function [patterns,frex] = validate_freq(FREQ,label)

if ~isstruct(FREQ) || ~isfield(FREQ,'pats') || isempty(FREQ.pats)
    error('%s.pats is required.',label);
end
patterns = FREQ.pats;
if ~isnumeric(patterns) || ~isreal(patterns) || any(~isfinite(patterns(:))) || ...
        ndims(patterns)>4
    error(['%s.pats must be a finite real array with dimensions ' ...
        '[nVox x nComp x nFrex x (nSubs)].'],label);
end
if ~isfield(FREQ,'frex') || isempty(FREQ.frex)
    error('%s.frex is required.',label);
end
frex = FREQ.frex(:);
if ~isnumeric(frex) || ~isreal(frex) || any(~isfinite(frex)) || ...
        numel(frex)~=size(patterns,3) || any(diff(frex)<=0)
    error('%s.frex must be finite, strictly increasing, and match %s.pats.',label,label);
end

end

%% Helper Function: Select Frequencies
function [requested,resolved,indices] = select_frequencies(frex,requested_input,tolerance)

if ischar(requested_input) || (isstring(requested_input) && isscalar(requested_input))
    if ~strcmpi(requested_input,'all')
        error('frex2compare must be numeric or the string ''all''.');
    end
    requested = frex;
    resolved = frex;
    indices = (1:numel(frex))';
    return
end
if ~isnumeric(requested_input) || ~isreal(requested_input) || ...
        isempty(requested_input) || any(~isfinite(requested_input(:)))
    error('frex2compare must contain at least one finite numeric value.');
end
requested = requested_input(:);
if any(requested<frex(1)) || any(requested>frex(end))
    error('Requested frex2compare values fall outside the available range.');
end
if isempty(tolerance)
    if isscalar(frex)
        tolerance = inf;
    else
        tolerance = median(diff(frex))/2;
    end
elseif ~isnumeric(tolerance) || ~isscalar(tolerance) || ...
        ~isfinite(tolerance) || tolerance<0
    error('frex_tolerance must be one non-negative finite scalar.');
end

indices = zeros(numel(requested),1);
for frexi = 1:numel(requested)
    distances = abs(frex-requested(frexi));
    minimum = min(distances);
    nearest = find(abs(distances-minimum)<=eps(max(max(distances),1))*8);
    if numel(nearest)~=1
        error('Frequency %.4g Hz is equally close to two bins.',requested(frexi));
    end
    if minimum>tolerance+eps
        error('No analyzed frequency lies within %.4g Hz of %.4g Hz.', ...
            tolerance,requested(frexi));
    end
    indices(frexi) = nearest;
end
if numel(unique(indices))~=numel(indices)
    error('Two frex2compare values resolve to the same analyzed bin.');
end
resolved = frex(indices);
if any(abs(requested-resolved)>1e-12)
    warning(['Some requested frex2compare values are not present in FREQ.frex. ' ...
        'Using closest matches instead.']);
end

end

%% Helper Function: Select Components
function comps = select_components(requested,nComp)

if isempty(requested)
    disp('comps2compare not specified. Defaulting to the 1st component.');
    comps = 1;
    return
end
if ~isnumeric(requested) || ~isreal(requested) || ~isvector(requested) || ...
        any(~isfinite(requested)) || any(requested~=round(requested)) || ...
        any(requested<1) || any(requested>nComp) || ...
        numel(unique(requested))~=numel(requested)
    error('comps2compare must contain unique integers between 1 and %d.',nComp);
end
comps = requested(:)';

end

%% Helper Function: Prepare Maps
function [maps_a,maps_b,flips_a,flips_b] = ...
    prepare_maps(maps_a,maps_b,map_mode,normalization,design)

[~,nfrex,nComp,n_a] = size(maps_a);
n_b = size(maps_b,4);
flips_a = ones(nfrex,nComp,n_a);
flips_b = ones(nfrex,nComp,n_b);
if strcmp(map_mode,'magnitude')
    maps_a = abs(maps_a);
    maps_b = abs(maps_b);
else
    for frexi = 1:nfrex
        for compi = 1:nComp
            if strcmp(design,'within')
                for subi = 1:n_a
                    if maps_a(:,frexi,compi,subi)'* ...
                            maps_b(:,frexi,compi,subi)<0
                        maps_b(:,frexi,compi,subi) = ...
                            -maps_b(:,frexi,compi,subi);
                        flips_b(frexi,compi,subi) = -1;
                    end
                end
                reference_maps = .5*(squeeze(maps_a(:,frexi,compi,:))+ ...
                    squeeze(maps_b(:,frexi,compi,:)));
            else
                reference_maps = [squeeze(maps_a(:,frexi,compi,:)) ...
                    squeeze(maps_b(:,frexi,compi,:))];
            end
            [reference,~,~] = svd(reference_maps,'econ');
            reference = reference(:,1);
            [~,peak] = max(abs(reference));
            if reference(peak)<0
                reference = -reference;
            end
            for subi = 1:n_a
                if maps_a(:,frexi,compi,subi)'*reference<0
                    maps_a(:,frexi,compi,subi) = -maps_a(:,frexi,compi,subi);
                    flips_a(frexi,compi,subi) = -flips_a(frexi,compi,subi);
                end
            end
            for subi = 1:n_b
                if maps_b(:,frexi,compi,subi)'*reference<0
                    maps_b(:,frexi,compi,subi) = -maps_b(:,frexi,compi,subi);
                    flips_b(frexi,compi,subi) = -flips_b(frexi,compi,subi);
                end
            end
        end
    end
end

if strcmp(normalization,'l2')
    for frexi = 1:nfrex
        for compi = 1:nComp
            for subi = 1:n_a
                norm_a = norm(maps_a(:,frexi,compi,subi));
                if norm_a<=eps
                    error('A selected participant pattern has zero L2 norm.');
                end
                maps_a(:,frexi,compi,subi) = maps_a(:,frexi,compi,subi)/norm_a;
            end
            for subi = 1:n_b
                norm_b = norm(maps_b(:,frexi,compi,subi));
                if norm_b<=eps
                    error('A selected FREQ_B participant pattern has zero L2 norm.');
                end
                maps_b(:,frexi,compi,subi) = maps_b(:,frexi,compi,subi)/norm_b;
            end
        end
    end
end

end

%% Helper Function: Spatial Statistics
function [mean_difference,statistic,effect_size,p_values, ...
    statistic_name,effect_size_name] = within_statistics(differences,test)

if strcmp(test,'ttest')
    [mean_difference,statistic,effect_size,p_values] = ...
        paired_statistics(differences);
    statistic_name = 't';
    effect_size_name = 'Cohen''s dz';
    return
end

mean_difference = mean(differences,4);
statistic = zeros(size(mean_difference));
effect_size = zeros(size(mean_difference));
p_values = ones(size(mean_difference));
flat_differences = reshape(differences,[],size(differences,4));
for mapi = 1:size(flat_differences,1)
    this_difference = flat_differences(mapi,:);
    this_difference = this_difference(this_difference~=0);
    n_nonzero = numel(this_difference);
    if n_nonzero==0
        continue
    end
    ranks = tiedrank(abs(this_difference));
    signed_rank = sum(sign(this_difference).*ranks);
    variance = n_nonzero*(n_nonzero+1)*(2*n_nonzero+1)/6;
    statistic(mapi) = signed_rank/sqrt(variance);
    effect_size(mapi) = signed_rank/(n_nonzero*(n_nonzero+1)/2);
    p_values(mapi) = erfc(abs(statistic(mapi))/sqrt(2));
end
statistic_name = 'z';
effect_size_name = 'matched rank-biserial r';

end

function [mean_difference,statistic,effect_size,p_values, ...
    statistic_name,effect_size_name] = between_statistics(maps_a,maps_b,test)

n_a = size(maps_a,4);
n_b = size(maps_b,4);
mean_a = mean(maps_a,4);
mean_b = mean(maps_b,4);
mean_difference = mean_b-mean_a;

if strcmp(test,'ttest')
    variance_a = var(maps_a,0,4);
    variance_b = var(maps_b,0,4);
    term_a = variance_a/n_a;
    term_b = variance_b/n_b;
    denominator = sqrt(term_a+term_b);
    statistic = zeros(size(mean_difference));
    valid = denominator>eps;
    statistic(valid) = mean_difference(valid)./denominator(valid);
    zero_variance_effect = ~valid & abs(mean_difference)>eps;
    statistic(zero_variance_effect) = ...
        sign(mean_difference(zero_variance_effect))*Inf;
    degrees_freedom = inf(size(mean_difference));
    degrees_freedom(valid) = (term_a(valid)+term_b(valid)).^2./ ...
        (term_a(valid).^2/(n_a-1)+term_b(valid).^2/(n_b-1));
    x = degrees_freedom./(degrees_freedom+statistic.^2);
    p_values = betainc(x,degrees_freedom/2,.5);
    p_values(~valid & mean_difference==0) = 1;
    p_values(zero_variance_effect) = 0;

    pooled = sqrt(((n_a-1)*variance_a+(n_b-1)*variance_b)/(n_a+n_b-2));
    effect_size = zeros(size(mean_difference));
    valid_pooled = pooled>eps;
    effect_size(valid_pooled) = ...
        mean_difference(valid_pooled)./pooled(valid_pooled);
    effect_size(zero_variance_effect) = ...
        sign(mean_difference(zero_variance_effect))*Inf;
    correction = 1-3/(4*(n_a+n_b)-9);
    effect_size = correction*effect_size;
    statistic_name = 't';
    effect_size_name = 'Hedges'' g';
    return
end

statistic = zeros(size(mean_difference));
effect_size = zeros(size(mean_difference));
p_values = ones(size(mean_difference));
flat_a = reshape(maps_a,[],n_a);
flat_b = reshape(maps_b,[],n_b);
variance = n_a*n_b*(n_a+n_b+1)/12;
for mapi = 1:size(flat_a,1)
    ranks = tiedrank([flat_b(mapi,:) flat_a(mapi,:)]);
    u_b = sum(ranks(1:n_b))-n_b*(n_b+1)/2;
    statistic(mapi) = (u_b-n_a*n_b/2)/sqrt(variance);
    effect_size(mapi) = 2*u_b/(n_a*n_b)-1;
    p_values(mapi) = erfc(abs(statistic(mapi))/sqrt(2));
end
statistic_name = 'z';
effect_size_name = 'rank-biserial r';

end

function [mean_difference,t_values,effect_size,p_values] = ...
    paired_statistics(differences)

nsubs = size(differences,4);
mean_difference = mean(differences,4);
standard_deviation = std(differences,0,4);
denominator = standard_deviation/sqrt(nsubs);
t_values = zeros(size(mean_difference));
valid = denominator>eps;
t_values(valid) = mean_difference(valid)./denominator(valid);
zero_variance_effect = ~valid & abs(mean_difference)>eps;
t_values(zero_variance_effect) = sign(mean_difference(zero_variance_effect))*inf;
effect_size = zeros(size(mean_difference));
valid_sd = standard_deviation>eps;
effect_size(valid_sd) = mean_difference(valid_sd)./standard_deviation(valid_sd);
effect_size(zero_variance_effect) = sign(mean_difference(zero_variance_effect))*inf;
degrees_freedom = nsubs-1;
x = degrees_freedom./(degrees_freedom+t_values.^2);
p_values = betainc(x,degrees_freedom/2,.5);
p_values(t_values==0) = 1;

end

%% Helper Function: Spatial Adjacency
function [adjacency,spacing] = spatial_adjacency(MNI,connectivity)

nvoxs = size(MNI,1);
spacing = nan(1,3);
scaled = zeros(size(MNI));
for dim = 1:3
    values = unique(MNI(:,dim));
    gaps = diff(values);
    gaps = gaps(gaps>eps);
    if ~isempty(gaps)
        spacing(dim) = min(gaps);
        scaled(:,dim) = (MNI(:,dim)-values(1))/spacing(dim);
    end
end
adjacency = cell(nvoxs,1);
if connectivity==6
    maximum_axes = 1;
elseif connectivity==18
    maximum_axes = 2;
else
    maximum_axes = 3;
end
for voxi = 1:nvoxs
    delta = abs(scaled-scaled(voxi,:));
    on_grid = all(abs(delta-round(delta))<1e-6,2);
    immediate = all(delta<=1+1e-6,2);
    changed_axes = sum(delta>1e-6,2);
    adjacency{voxi} = find(on_grid & immediate & changed_axes>=1 & ...
        changed_axes<=maximum_axes);
end

end

%% Helper Function: Connected Components
function clusters = connected_components(mask,adjacency)

visited = false(size(mask));
clusters = {};
for seed = find(mask(:))'
    if visited(seed)
        continue
    end
    stack = seed;
    visited(seed) = true;
    members = zeros(0,1);
    while ~isempty(stack)
        current = stack(end);
        stack(end) = [];
        members(end+1,1) = current; %#ok<AGROW>
        neighbors = adjacency{current};
        new_neighbors = neighbors(mask(neighbors) & ~visited(neighbors));
        visited(new_neighbors) = true;
        stack = [stack; new_neighbors(:)]; %#ok<AGROW>
    end
    clusters{end+1} = members; %#ok<AGROW>
end

end

%% Helper Function: Maximum Cluster Mass
function maximum = maximum_cluster_mass(statistic,p_values, ...
    cluster_forming_alpha,adjacency)

maximum = 0;
for frexi = 1:size(statistic,2)
    for compi = 1:size(statistic,3)
        spatial_map = statistic(:,frexi,compi);
        significant = p_values(:,frexi,compi)<=cluster_forming_alpha;
        for polarity = [1 -1]
            if polarity>0
                mask = significant & spatial_map>0;
            else
                mask = significant & spatial_map<0;
            end
            clusters = connected_components(mask,adjacency);
            for clusteri = 1:numel(clusters)
                maximum = max(maximum,sum(abs(spatial_map(clusters{clusteri}))));
            end
        end
    end
end

end

%% Helper Function: Permutation Plans
function [signs,exact] = ...
    permutation_plan_within(nsubs,n_permutations,random_seed)

possible = 2^nsubs;
if ischar(n_permutations) || (isstring(n_permutations) && isscalar(n_permutations))
    mode = validatestring(n_permutations,{'auto','exact'});
    if strcmp(mode,'exact') && possible>1048576
        error('Exact enumeration is too large; use Monte Carlo permutations.');
    end
    exact = strcmp(mode,'exact') || possible<=65536;
    if exact
        count = possible;
    else
        count = 10000;
    end
elseif isnumeric(n_permutations) && isscalar(n_permutations) && ...
        isfinite(n_permutations) && n_permutations>0 && ...
        n_permutations==round(n_permutations)
    count = n_permutations;
    exact = count>=possible && possible<=65536;
    if exact
        count = possible;
    end
else
    error('n_permutations must be a positive integer, ''auto'', or ''exact''.');
end

if exact
    configuration = uint64((0:count-1)');
    signs = zeros(count,nsubs);
    for subi = 1:nsubs
        signs(:,subi) = 2*double(bitget(configuration,subi))-1;
    end
else
    previous_rng = rng;
    cleanup = onCleanup(@()rng(previous_rng));
    rng(random_seed,'twister');
    signs = 2*randi([0 1],count,nsubs)-1;
end

end

function [assignments,exact] = ...
    permutation_plan_between(n_a,n_b,n_permutations,random_seed)

total = n_a+n_b;
possible = nchoosek(total,n_a);
if ischar(n_permutations) || ...
        (isstring(n_permutations) && isscalar(n_permutations))
    mode = validatestring(n_permutations,{'auto','exact'});
    if strcmp(mode,'exact') && possible>1048576
        error('Exact enumeration is too large; use Monte Carlo permutations.');
    end
    exact = strcmp(mode,'exact') || possible<=65536;
    if exact
        count = possible;
    else
        count = 10000;
    end
elseif isnumeric(n_permutations) && isscalar(n_permutations) && ...
        isfinite(n_permutations) && n_permutations>0 && ...
        n_permutations==round(n_permutations)
    count = n_permutations;
    exact = count>=possible && possible<=65536;
    if exact
        count = possible;
    end
else
    error('n_permutations must be a positive integer, ''auto'', or ''exact''.');
end

assignments = false(count,total);
if exact
    combinations = nchoosek(1:total,n_a);
    for permi = 1:count
        assignments(permi,combinations(permi,:)) = true;
    end
else
    previous_rng = rng;
    cleanup = onCleanup(@()rng(previous_rng));
    rng(random_seed,'twister');
    for permi = 1:count
        assignments(permi,randperm(total,n_a)) = true;
    end
end

end

%% Helper Function: Plot Results
function plot_spatial_results(MNI,STATS)

cmap = diverging_colormap(256);
for frexi = 1:numel(STATS.frex)
    for compi = 1:numel(STATS.comps)
        difference = STATS.mean_difference(:,frexi,compi);
        significant = STATS.significant_mask(:,frexi,compi);
        color_limit = max(max(abs(difference)),eps);
        current_figure = figure('Color','w','Position',[100 100 1120 520], ...
            'Renderer','painters');
        layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

        nexttile, hold on
        scatter3(MNI(:,1),MNI(:,2),MNI(:,3),12,[.80 .83 .86], ...
            'filled','MarkerFaceAlpha',.68,'MarkerEdgeColor','none');
        scatter3(MNI(:,1),MNI(:,2),MNI(:,3),36,difference,'filled', ...
            'MarkerFaceAlpha',.95,'MarkerEdgeColor','none');
        format_spatial_axis('Mean difference',color_limit);

        nexttile, hold on
        scatter3(MNI(:,1),MNI(:,2),MNI(:,3),12,[.80 .83 .86], ...
            'filled','MarkerFaceAlpha',.68,'MarkerEdgeColor','none');
        if any(significant)
            scatter3(MNI(significant,1),MNI(significant,2),MNI(significant,3), ...
                54,difference(significant),'filled','MarkerFaceAlpha',.95, ...
                'MarkerEdgeColor',[.09 .13 .17],'LineWidth',.75);
        end
        format_spatial_axis('Cluster-corrected difference',color_limit);
        colormap(current_figure,cmap)
        current_colorbar = colorbar;
        current_colorbar.Layout.Tile = 'east';
        current_colorbar.Label.String = 'Mean difference';
        current_colorbar.FontSize = 10;
        title(layout,{sprintf('%s - %s  |  %.4g Hz  |  Component %d', ...
            STATS.condition_labels{2},STATS.condition_labels{1}, ...
            STATS.frex(frexi),STATS.comps(compi)), ...
            sprintf('Outlined markers: cluster-corrected p <= %.3g', ...
            STATS.alpha)},'FontSize',15,'FontWeight','bold');
    end
end

end

function format_spatial_axis(title_text,color_limit)

format_publication_spatial_axis([])
title(title_text,'FontSize',12,'FontWeight','bold')
clim([-color_limit color_limit])

end

function format_publication_spatial_axis(MNI)

axis equal off
axis vis3d
view(-62,18)
camproj orthographic
if ~isempty(MNI)
    xlim([min(MNI(:,1)) max(MNI(:,1))])
    ylim([min(MNI(:,2)) max(MNI(:,2))])
    zlim([min(MNI(:,3)) max(MNI(:,3))])
end
set(gca,'FontName','Helvetica','FontSize',10,'Clipping','off')

end

function cmap = publication_sequential_colormap(ncolors)

anchors = [0.2670 0.0049 0.3294; 0.2297 0.3224 0.5457; ...
    0.1276 0.5669 0.5506; 0.3692 0.7889 0.3829; ...
    0.9932 0.9062 0.1439];
source = linspace(0,1,size(anchors,1));
target = linspace(0,1,ncolors);
cmap = interp1(source,anchors,target,'pchip');
cmap = min(max(cmap,0),1);

end

function cmap = diverging_colormap(ncolors)

half = ceil(ncolors/2);
blue = [linspace(.15,1,half)' linspace(.3,1,half)' ones(half,1)];
red = [ones(ncolors-half,1) linspace(1,.2,ncolors-half)' ...
    linspace(1,.2,ncolors-half)'];
cmap = [blue;red];

end

%% Helper Function: Validation
function [design,test] = validate_analysis_options(design,test)

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
    test = 'ttest';
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

end

function labels = validate_condition_labels(labels)

if isstring(labels)
    labels = cellstr(labels);
end
if ~iscell(labels) || numel(labels)~=2 || ...
        ~all(cellfun(@(value)ischar(value) && ~isempty(value),labels))
    error('condition_labels must contain two non-empty text labels.');
end
labels = reshape(labels,1,2);

end

function validate_probability(value,name)

if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value<=0 || value>=1
    error('%s must lie strictly between 0 and 1.',name);
end

end

function clusters = empty_cluster_structure

clusters = struct('frex',{},'comp',{},'polarity',{},'size',{}, ...
    'mass',{},'p_corrected',{},'peak_t',{},'peak_statistic',{}, ...
    'peak_voxel',{},'peak_mni',{});

end

function output = ternary(condition,true_value,false_value)

if condition
    output = true_value;
else
    output = false_value;
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
