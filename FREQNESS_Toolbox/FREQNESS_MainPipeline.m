% ========================================================================
%  FREQNESS: Main analysis pipeline for brain network estimation
%
%  If you use this toolbox, please cite:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P., Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
%
% ========================================================================
%
%  This script constitutes the main analysis pipeline for
%  FREQNESS, designed to handle multi-participant, multi-group,
%  and multi-condition datasets in a fully automated workflow.
%
%  Simply place your data in the appropriate folders and press RUN.
%
%  FREQNESS will take care of everything else: loading all  participants,
%  estimating frequency-resolved brain networks, generating publication-ready
%  visualizations and group-level summaries for several secondary analyses.
%
%  The new pipeline operates on a simple folder structure:
%  each Group or Condition is represented by a folder, and each folder
%  contains one `.mat` file per participant, each storing a single
%  data matrix in voxels-by-time format.
%
%  Additionally, the user provides one `.mat` file containing the
%  voxel MNI coordinates (Nvoxels × 3).
%  Using this information, the script:
%    1) Imports all datasets across Groups/Conditions
%    2) Runs frequency-resolved GED (via FREQNESS_NetworkEstimation)
%    3) Aggregates outputs at group-level when applicable
%    4) Produces a complete set of visualizations, spatial maps, and
%       component summaries through FREQNESS_Visualizer
%
%  NOTE: while MNI coordinates are not necessary to run the core functions,
%  some of the secondary functions will not be able to process some
%  topographical aspects of the networks.
%
%  NOTE #2: after downloading the toolbox, please do NOT change the structure
%  of its folders and subfolders.
%
% ------------------------------------------------------------------------
%  FUNCTIONS OVERVIEW:
% ------------------------------------------------------------------------
%  The script loads the example datasets into the MATLAB workspace and
%  calls the following functions:
%
% - FREQNESS_Startup()
%    Initializes the environment, sets up the necessary directories, and
%    fetches data and MNI coordinates from the respective folders.
%
%  - FREQNESS_NetworkEstimation(Data, Frequencies, SamplingRate)
%    Performs Generalized Eigendecomposition (GED) over a user-defined
%    sample of frequencies, separating frequency-resolved networks from
%    source-reconstructed brain voxel data.
%
%  The output structure of the core function 'FREQNESS_NetworkEstimation'
%  will be processed by the following second-order functions to perform a
%  series of analyses.
%
%  - FREQNESS_Visualizer(FREQ, Landscape, Patterns, ...)
%    Takes as input selected outputs from `FREQNESS_NetworkEstimation`
%    to visualize both the network landscape and brain topographies.
%
%  - FREQNESS_EntropyLandscape(FREQ)
%    Computes quadratic Rényi entropy (H2) and effective dimensionality
%    (ED) across frequencies. Both measures are returned for numerical
%    analysis, while only the H2 entropy landscape is visualized.
%
%  - FREQNESS_ExponentialDK(FREQ, ...)
%    Fits an exponentially decaying function to the eigenvalues of a
%    selected component across frequencies. The resulting decay
%    coefficients provide a compact measure of how strongly that component
%    is dominated by low- vs high-frequency contributions, which can be
%    compared across groups or conditions.
%
%  - FREQNESS_FreqGradients(FREQ, MNI, ...)
%    Models spatial gradients of activation patterns across frequencies
%    for a chosen component. Using voxel-wise MNI coordinates, it fits
%    linear/quadratic polynomials to capture how dominant activations shift
%    along X, Y and Z as a function of frequency, returning subject-wise
%    gradient coefficients and goodness-of-fit metrics.
%
%  - FREQNESS_CompGradients(FREQ, MNI, ...)
%    Analogous to FREQNESS_FreqGradients, but parametrizes spatial
%    gradients across components at a selected frequency. It tests
%    whether network activation shifts systematically along X, Y and Z as a
%    function of component index, again providing subject-wise polynomial
%    coefficients and fit quality.
%
%  - FREQNESS_CrossCoupling(FREQ, lfo_freq, ...)
%    Computes phase–amplitude cross-frequency coupling (PAC) between a low-frequency
%    network and higher-frequency carrier networks for a selected component.
%    For each carrier frequency, carrier power is binned by LFO phase to obtain
%    PAC histograms and deterministic first-harmonic parameters, with
%    optional 3D visualizations of the LFO and peak carrier networks in MNI space.
%
%  - FREQNESS_BackProjection(FREQ, freq2project, ...)
%    Reconstructs selected frequency-resolved components in broadband voxel
%    space using the GED forward model and component time series.
%
%  - FREQNESS_NetworkRemoval(FREQ, data, freq2remove, ...)
%    Removes selected backprojected network activity from the original data.
%    The cleaned data can then be re-estimated to inspect the residual
%    frequency-resolved network landscape.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso & Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 23/11/2025
% ========================================================================

%% ========================================================================
% STARTUP
% ========================================================================

clear
close all
clc

% Setup directories, fetch data and MNI coordinates
[allData, MNI, path_home] = FREQNESS_Startup();

%% ========================================================================
% DEFINE ANALYSIS SETTINGS
%
%  Edit this section before running the pipeline
% ========================================================================

% ------------------------------------------------------------------------
% 1) FREQNESS_NetworkEstimation
% ------------------------------------------------------------------------
frex   = 1.2:1.2:20*1.2; % frequency vector (Hz) used for the GED analysis
srate  = 250;             % sampling rate (Hz) of your data

% Optional network-estimation settings
network_duration       = [];            % seconds; leave [] to use all data
network_fwidth         = [];            % scalar or one value per frequency
network_filter         = 'logarithmic'; % or 'linear'
network_regularisation = 0.01;          % covariance shrinkage factor
network_ncomps         = 30;            % components retained per frequency
network_bad_segments   = [];            % sample indices excluded from covariance
network_rescale        = false;         % rescale very low-amplitude data

% ------------------------------------------------------------------------
% 2) FREQNESS_Visualizer
% ------------------------------------------------------------------------
% Flag to visualize all participants (discouraged for large samples)
plot_all   = false;
save_nifti = true; % save 8-mm NIFTI activation maps
% Brain network landscape
Landscape = [];
Landscape.frex   = frex; % assign 'frex' to visualize all frequencies
Landscape.ncomps = 10;   % how many components in the network landscape
% Spatial activation patterns
Patterns = [];
Patterns.frex    = [2.4 8.4];      % select frequencies
Patterns.ncomps  = 1;              % set how many top components to visualize
Patterns.path_output =  path_home; % set output path to save nifti images
Patterns.MNI_coords = MNI; % assigning the MNI coordinates of your data
% NOTE: ensure that the order of the MNI coordinates matches YOUR OWN DATA!

% ------------------------------------------------------------------------
% 3) FREQNESS_EntropyLandscape
% ------------------------------------------------------------------------
% This function only requires the FREQ structure, so no additional settings
% are strictly necessary. Settings might be needed in future extensions.

% ------------------------------------------------------------------------
% 4) FREQNESS_ExponentialDK
% ------------------------------------------------------------------------
expdk_range2fit  = [9 20]; % frequency range (Hz) to fit; leave [] to use all frequencies
expdk_which_comp = 1;      % component index to analyse

% ------------------------------------------------------------------------
% 5) FREQNESS_FreqGradients
% ------------------------------------------------------------------------
freqgrad_frex2model = [8 12]; % frequency range (Hz) over which to model gradients
freqgrad_comp2model = 1;      % component index to analyse

% ------------------------------------------------------------------------
% 6) FREQNESS_CompGradients
% ------------------------------------------------------------------------
compgrad_freq2model  = 10;    % frequencies (Hz) at which to analyse gradients across components
compgrad_comps2model = [1 5]; % range of components to include in the gradient fit

% ------------------------------------------------------------------------
% 7) FREQNESS_CrossCoupling
% ------------------------------------------------------------------------
lfo_freq = 2; % low-frequency oscillator (Hz)

% ------------------------------------------------------------------------
% 8) FREQNESS_BackProjection
% ------------------------------------------------------------------------
backproj_freq2project  = 8.4; % network frequency to backproject (Hz)
backproj_comps2project = 1;   % component(s) to backproject

% ------------------------------------------------------------------------
% 9) FREQNESS_NetworkRemoval
% ------------------------------------------------------------------------
netrem_freq2remove      = 8.4; % network frequency to remove (Hz)
netrem_comps2remove     = 1;   % component(s) to remove
netrem_plot_landscape   = true; % re-estimate and plot the cleaned landscape
netrem_landscape_ncomps = 3;    % cleaned components to visualize

%% ========================================================================
% INITIALIZE OUTPUT VARIABLES (one per Group / Condition)
% ========================================================================

% Number of experimental Groups / Conditions
nconds = numel(allData);

% 1) Core FREQNESS output
FREQ        = cell(nconds,1);
% Entropy landscape
h2          = cell(nconds,1);
ed          = cell(nconds,1);
% Exponential decay of eigenvalues
decayCoeff  = cell(nconds,1);
goodFit_exp = cell(nconds,1);
% Spatial gradients across frequencies
gradCoeff_f = cell(nconds,1);
goodFit_f   = cell(nconds,1);
% Spatial gradients across components
gradCoeff_c = cell(nconds,1);
goodFit_c   = cell(nconds,1);
% Cross-frequency coupling structure
CFC         = cell(nconds,1);
% Backprojected broadband network activity
backProj    = cell(nconds,1);
% Network-removal outputs and cleaned-data re-estimation
dataClean       = cell(nconds,1);
removedActivity = cell(nconds,1);
FREQ_clean      = cell(nconds,1);

% Iterate the pipeline over Conditions or Groups
for condi = 1:nconds
    fprintf('\nANALYSING CONDITION/GROUP #%d\n',condi);

    % Skip empty conditions (no .mat files in folder)
    if isempty(allData{condi})
        warning('Skipping condition #%d: no data found in FREQNESS_Data/Dataset_%d.', condi,condi);
        continue
    end

    %% ========================================================================
    % 1) FREQNESS network estimation
    % ========================================================================

    % Core function with the settings defined above
    FREQ{condi} = FREQNESS_NetworkEstimation(allData{condi},frex,srate, ...
        'duration',network_duration, ...
        'fwidth',network_fwidth, ...
        'filter',network_filter, ...
        'regularisation',network_regularisation, ...
        'ncomps',network_ncomps, ...
        'bad_segments',network_bad_segments, ...
        'rescale',network_rescale);


    %% ========================================================================
    % 2) FREQNESS VISUALIZATION
    % ========================================================================

    % Plot network landscape and, when MNI coordinates are available,
    % spatial activation patterns and optional NIFTI images
    if isempty(MNI)
        warning(['MNI coordinates are unavailable or do not match the data. ' ...
            'Only the network landscape will be visualized.']);
        FREQNESS_Visualizer(FREQ{condi},Landscape,[], ...
            'plot_all',plot_all,'save_nifti',false)
    else
        FREQNESS_Visualizer(FREQ{condi},Landscape,Patterns, ...
            'plot_all',plot_all,'save_nifti',save_nifti)
    end

    % NOTE: computation and storage of NIFTI files supported only for 8mm MNI space

    %% ========================================================================
    % 3) FREQNESS ENTROPY LANDSCAPE
    % ========================================================================

    % Compute entropy-based indices of the eigenspectrum
    [h2{condi}, ed{condi}] = FREQNESS_EntropyLandscape(FREQ{condi});


    %% ========================================================================
    % 4) FREQNESS EXPONENTIAL DECAY (Eigenvalue decay across frequency)
    % ========================================================================

    [decayCoeff{condi}, goodFit_exp{condi}] = FREQNESS_ExponentialDK(FREQ{condi}, ...
        'which_comp', expdk_which_comp, ...
        'range2fit',  expdk_range2fit,...
        'plot_all',plot_all);

    %% ========================================================================
    % 5) FREQNESS FREQUENCY GRADIENTS (Spatial gradients across frequencies)
    % ========================================================================

    if ~isempty(MNI) % run only if coordinates are provided

        [gradCoeff_f{condi}, goodFit_f{condi}] = FREQNESS_FreqGradients(FREQ{condi}, MNI, ...
            'frex2model', freqgrad_frex2model, ...
            'comp2model', freqgrad_comp2model,...
            'plot_all',plot_all);
    end

    %% ========================================================================
    % 6) FREQNESS COMPONENT GRADIENTS (Spatial gradients across components)
    % ========================================================================

    if ~isempty(MNI) % run only if coordinates are provided

        [gradCoeff_c{condi}, goodFit_c{condi}] = FREQNESS_CompGradients(FREQ{condi}, MNI, ...
            'freq2model',  compgrad_freq2model, ...
            'comps2model', compgrad_comps2model,...
            'plot_all',plot_all);

    end

    %% ========================================================================
    % 7) FREQNESS CROSS FREQUENCY COUPLING (Phase-amplitude coupling)
    % ========================================================================

    CFC{condi} = FREQNESS_CrossCoupling(FREQ{condi}, lfo_freq, 'mni', MNI);

    %% ========================================================================
    % 8) FREQNESS BACKPROJECTION
    % ========================================================================

    backProj{condi} = FREQNESS_BackProjection(FREQ{condi}, ...
        backproj_freq2project, ...
        'comps2project',backproj_comps2project);

    %% ========================================================================
    % 9) FREQNESS NETWORK REMOVAL
    % ========================================================================

    % Match the exact data segment analyzed by FREQNESS_NetworkEstimation.
    data2remove = allData{condi}(:,1:size(FREQ{condi}.ts,2),:);
    [dataClean{condi}, removedActivity{condi}] = FREQNESS_NetworkRemoval( ...
        FREQ{condi}, data2remove, netrem_freq2remove, ...
        'comps2remove',netrem_comps2remove);

    % Re-estimation belongs to the pipeline rather than NetworkRemoval:
    % dataClean remains available for any subsequent analysis chosen by users.
    if netrem_plot_landscape
        FREQ_clean{condi} = FREQNESS_NetworkEstimation(dataClean{condi},frex,srate, ...
            'duration',network_duration, ...
            'fwidth',network_fwidth, ...
            'filter',network_filter, ...
            'regularisation',network_regularisation, ...
            'ncomps',network_ncomps, ...
            'bad_segments',network_bad_segments, ...
            'rescale',network_rescale);
        Landscape_clean = [];
        Landscape_clean.frex = frex;
        Landscape_clean.ncomps = netrem_landscape_ncomps;
        FREQNESS_Visualizer(FREQ_clean{condi},Landscape_clean,[], ...
            'plot_all',plot_all)
    end

end

%% ========================================================================
%  FINAL REMARKS
% ========================================================================
%
% This pipeline produced several FREQNESS outputs, enabling direct
% comparison of frequency-specific brain functional connectivity across
% experimental Conditions or Groups.
%
% For statistical analysis, decide which FREQNESS outputs
% you want to compare and analyze them across experimental
% Conditions, Groups, or any other relevant factor
% based on your study design.
%
% Please check the FREQNESS GitHub repository for new releases.
% Future updates will include statistical analysis tools,
% cross-frequency coupling between brain networks,
% brain network-induced responses, and more.
% https://github.com/mattiaRosso92/Frequency-resolved_brain_network_estimation_via_source_separation_FREQ-NESS.git
%
% Feel free to reach out to us if you need guidance or consultation.
% Mattia Rosso:     mattia.rosso@clin.au.dk
% Leonardo Bonetti: leonardo.bonetti@clin.au.dk
%                   leonardo.bonetti@psych.ox.ac.uk
%
%  If you use this toolbox, please cite:
%  M. Rosso, G. Fernández-Rubio, P. E. Keller, E. Brattico, P. Vuust, M. L. Kringelbach, L. Bonetti.
%  FREQ-NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved Brain Networks During Auditory Stimulation.
%  Adv. Sci. 2025, 2413195.
%  https://doi.org/10.1002/advs.202413195

%%
