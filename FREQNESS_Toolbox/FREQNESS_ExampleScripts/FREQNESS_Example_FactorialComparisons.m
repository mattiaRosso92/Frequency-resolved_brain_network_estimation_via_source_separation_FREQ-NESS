% ========================================================================
%  FREQNESS: EXAMPLE OF FACTORIAL STATISTICAL COMPARISONS
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
%  This template shows how to run FREQNESS_SpatialComparison and
%  FREQNESS_LandscapeComparison for a 2-by-4 mixed design:
%
%       Factor 1: Group     -> Control, Clinical        (between subjects)
%       Factor 2: Condition -> Rest, Task1, Task2, Recovery
%                                                           (within subjects)
%
%  Replace the eight example FREQ variable names and participant identifiers
%  with those used in your study. The statistical calls can then be copied
%  directly into an existing pipeline.
%
%  Required structure of every FREQ_GRID entry:
%       FREQ.pats  [nVoxels x nComponents x nFrequencies x nParticipants]
%       FREQ.evals [nComponents x nFrequencies x nParticipants]
%       FREQ.frex  [1 x nFrequencies]
%
%  All cells must share the same voxels, components, and frequency grid.
%  Different groups may contain different participant counts. Repeated
%  conditions within one group must contain the same participants.
%
%  MNI must have dimensions [nVoxels x 3] and follow the voxel order used
%  in FREQ.pats.
% ========================================================================

%% 1. ASSIGN THE STUDY DATA

% Load or create the following eight group-level FREQ structures before
% continuing. Replace these names with the variables in your pipeline:
%
%   FREQ_control_rest
%   FREQ_control_task1
%   FREQ_control_task2
%   FREQ_control_recovery
%   FREQ_clinical_rest
%   FREQ_clinical_task1
%   FREQ_clinical_task2
%   FREQ_clinical_recovery
%
% Also load the matching MNI coordinate matrix:
%
%   MNI [nVoxels x 3]

% The dimensions of this cell array define the experimental design.
% Rows represent Group; columns represent Condition.
FREQ_GRID = cell(2,4);
FREQ_GRID{1,1} = FREQ_control_rest;
FREQ_GRID{1,2} = FREQ_control_task1;
FREQ_GRID{1,3} = FREQ_control_task2;
FREQ_GRID{1,4} = FREQ_control_recovery;
FREQ_GRID{2,1} = FREQ_clinical_rest;
FREQ_GRID{2,2} = FREQ_clinical_task1;
FREQ_GRID{2,3} = FREQ_clinical_task2;
FREQ_GRID{2,4} = FREQ_clinical_recovery;


%% 2. DESCRIBE THE FACTORS AND PARTICIPANTS

% One name, design type, and set of level labels per cell-array dimension.
factor_names  = {'Group','Condition'};
factor_design = {'between','within'};
level_labels  = {{'Control','Clinical'}, ...
                 {'Rest','Task1','Task2','Recovery'}};

% Supply one identifier per participant in the last dimension of each
% FREQ structure. Within a group, the same identifiers must appear in every
% repeated condition. Their order may differ because the functions reorder
% them automatically. Replace these example identifiers with your own.
control_ids  = {'C01','C02','C03','C04'};
clinical_ids = {'P01','P02','P03','P04','P05'};

subject_ids = cell(size(FREQ_GRID));
subject_ids(1,:) = {control_ids,control_ids,control_ids,control_ids};
subject_ids(2,:) = {clinical_ids,clinical_ids,clinical_ids,clinical_ids};


%% 3. SPATIAL COMPARISON

% Tests FREQ.pats voxel by voxel. The observed omnibus statistics are
% calculated directly. Design-aware permutations are used to correct the
% maximum cluster mass across effects, frequencies, and components.

spatial_frex                 = [8 12]; % Frequencies in Hz, or 'all'
spatial_components           = 1;      % One-based component numbers
spatial_test                 = 'anova';
spatial_effects              = 'all';  % Or {'Group','Condition', ...
                                       %     'Group_x_Condition'}
spatial_map_mode             = 'magnitude'; % Or 'signed_aligned'
spatial_normalization        = 'l2';        % Or 'none'
spatial_cluster_alpha        = .001;   % Cluster-forming voxel threshold
spatial_corrected_alpha      = .05;    % Corrected cluster threshold
spatial_connectivity         = 6;      % 6, 18, or 26 neighbours
spatial_frequency_tolerance  = [];     % [] uses half the frequency spacing
spatial_n_permutations       = 10000;  % Use fewer only while debugging
spatial_random_seed          = 0;
spatial_plot_results         = true;

STATS_SPATIAL = FREQNESS_SpatialComparison(FREQ_GRID,MNI, ...
    'frex2compare',spatial_frex, ...
    'comps2compare',spatial_components, ...
    'design',factor_design, ...
    'test',spatial_test, ...
    'factor_names',factor_names, ...
    'level_labels',level_labels, ...
    'subject_ids',subject_ids, ...
    'effects',spatial_effects, ...
    'map_mode',spatial_map_mode, ...
    'normalization',spatial_normalization, ...
    'cluster_forming_alpha',spatial_cluster_alpha, ...
    'alpha',spatial_corrected_alpha, ...
    'connectivity',spatial_connectivity, ...
    'frex_tolerance',spatial_frequency_tolerance, ...
    'correction_scope','all_effects_and_maps', ...
    'n_permutations',spatial_n_permutations, ...
    'random_seed',spatial_random_seed, ...
    'plot_results',spatial_plot_results);

% Main effects and interactions retain the factor names:
%   STATS_SPATIAL.effects.Group
%   STATS_SPATIAL.effects.Condition
%   STATS_SPATIAL.effects.Group_x_Condition


%% 4. LANDSCAPE COMPARISON

% FREQ itself is a valid landscape input because it contains FREQ.evals.
% Other valid grids can contain outputs such as:
%
%   ENTR.H2 / ENTR.ED                    [nFrex x nParticipants]
%   CFC.amplitude_raw                    [nCarrierFrex x nParticipants]
%   CFC.amplitude_normalized             [nCarrierFrex x nParticipants]
%   STAB.stability_index                 [nComp x nFrex x nParticipants]
%   STAB.mean_frequency                  [nComp x nFrex x nParticipants]
%
% In every case, create a cell array with the same factorial organization.

LANDSCAPE_GRID = FREQ_GRID;

landscape_measure             = 'evals'; % Omit/[] to analyze all shared fields
landscape_components          = 1;       % One-based component numbers
landscape_frex                = [];      % [] tests every stored frequency
landscape_test                = 'anova';
landscape_effects             = 'all';   % Or selected effect names
landscape_tail                = 'both';  % Omnibus tests are non-directional
landscape_alpha               = .05;
landscape_correction          = 'fdr';   % BH-FDR across frequencies
landscape_frequency_tolerance = [];
landscape_n_permutations      = 10000;   % Used only by 'permutation_anova'
landscape_random_seed         = 0;
landscape_plot_results        = true;

STATS_LANDSCAPE = FREQNESS_LandscapeComparison(LANDSCAPE_GRID, ...
    'measure',landscape_measure, ...
    'design',factor_design, ...
    'test',landscape_test, ...
    'which_comp',landscape_components, ...
    'frex2compare',landscape_frex, ...
    'factor_names',factor_names, ...
    'level_labels',level_labels, ...
    'subject_ids',subject_ids, ...
    'effects',landscape_effects, ...
    'tail',landscape_tail, ...
    'alpha',landscape_alpha, ...
    'correction',landscape_correction, ...
    'frex_tolerance',landscape_frequency_tolerance, ...
    'n_permutations',landscape_n_permutations, ...
    'random_seed',landscape_random_seed, ...
    'plot_results',landscape_plot_results);

% The selected measure is nested under every tested effect:
%   STATS_LANDSCAPE.effects.Group.evals
%   STATS_LANDSCAPE.effects.Condition.evals
%   STATS_LANDSCAPE.effects.Group_x_Condition.evals


%% 5. SIMPLER DESIGNS

% One factor with four repeated conditions:
%   FREQ_GRID = {FREQ_rest,FREQ_task1,FREQ_task2,FREQ_recovery};
%   factor_names = {'Condition'};
%   factor_design = {'within'};
%   level_labels = {{'Rest','Task1','Task2','Recovery'}};

% One factor with two independent groups:
%   FREQ_GRID = {FREQ_control,FREQ_clinical};
%   factor_names = {'Group'};
%   factor_design = {'between'};
%   level_labels = {{'Control','Clinical'}};

% The two function calls remain otherwise unchanged.

