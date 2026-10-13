function STATS = FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI,varargin)

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
%  This function performs a within-subject or between-group spatial
%  comparison of FREQ-NESS activation patterns. Voxel-wise statistics are
%  combined with a non-parametric cluster-mass permutation test based on MNI
%  adjacency.
%  The maximum cluster mass is calculated across all requested frequencies
%  and components in every permutation, controlling the family-wise error
%  rate across the complete analysis requested in one call.
%
%  For within-subject designs, participants must have the same order in
%  FREQ_A and FREQ_B. Between-group designs may have unequal sample sizes.
%  By default, absolute activation patterns are L2-normalized within
%  participant before comparison, making the analysis invariant to
%  arbitrary GED signs and global pattern scale.
% ========================================================================

% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - FREQ_A, FREQ_B : structures containing
%           • pats -> [nVoxels x nComp x nFrex x nSubs]
%           • frex -> analyzed frequencies in Hz
%
%  - MNI : MNI coordinates [nVoxels x 3], ordered as FREQ.pats.
%
%  - frex2compare : mandatory name-value input. Numeric values are mapped to
%                   the closest analyzed bins. Use 'all' for every bin.
%
%  - comps2compare : one-based component numbers. Default: 1.
%
%  - design : 'within' (default) or 'between'.
%
%  - test : 'ttest' (default), 'signrank', or 'ranksum'. An incompatible
%           non-parametric choice is corrected with a warning: ranksum is
%           replaced by signrank for within-subject designs, and signrank is
%           replaced by ranksum for between-group designs.
%
%  - condition_labels : two labels for FREQ_A and FREQ_B. Default: {'A','B'}.
%
%  - map_mode : 'magnitude' (default) or 'signed_aligned'.
%
%  - normalization : 'l2' (default) or 'none'.
%
%  - cluster_forming_alpha : two-sided voxel threshold. Default: .001.
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
        figure('Color','w','Position',[100 100 1100 480]);

        subplot(1,2,1), hold on
        scatter3(MNI(:,1),MNI(:,2),MNI(:,3),8,'k','filled', ...
            'MarkerFaceAlpha',.18,'MarkerEdgeAlpha',.18);
        scatter3(MNI(:,1),MNI(:,2),MNI(:,3),28,difference,'filled', ...
            'MarkerFaceAlpha',.9,'MarkerEdgeAlpha',.25);
        format_spatial_axis('Mean difference',cmap,color_limit);

        subplot(1,2,2), hold on
        scatter3(MNI(:,1),MNI(:,2),MNI(:,3),8,'k','filled', ...
            'MarkerFaceAlpha',.18,'MarkerEdgeAlpha',.18);
        if any(significant)
            scatter3(MNI(significant,1),MNI(significant,2),MNI(significant,3), ...
                28,difference(significant),'filled','MarkerFaceAlpha',.9, ...
                'MarkerEdgeAlpha',.25);
        end
        format_spatial_axis('Cluster-corrected difference',cmap,color_limit);
        sgtitle(sprintf('%s - %s | %.4g Hz | Component #%d', ...
            STATS.condition_labels{2},STATS.condition_labels{1}, ...
            STATS.frex(frexi),STATS.comps(compi)),'FontSize',16);
    end
end

end

function format_spatial_axis(title_text,cmap,color_limit)

axis equal
grid on
view(42,24)
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
title(title_text)
colormap(gca,cmap)
clim([-color_limit color_limit])
colorbar
set(gca,'FontSize',12)

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
