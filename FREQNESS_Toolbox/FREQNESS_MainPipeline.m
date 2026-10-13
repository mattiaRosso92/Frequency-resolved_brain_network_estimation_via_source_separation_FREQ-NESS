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
%  FREQNESS loads one participant at a time, estimates frequency-resolved
%  brain networks, and saves the complete FREQ structure before releasing
%  the participant's data and continuing with the next file.
%
%  The new pipeline operates on a simple folder structure:
%  each Group or Condition is represented by a folder, and each folder
%  contains one `.mat` file per participant, each storing a single
%  data matrix in voxels-by-time format.
%
%  Additionally, the user provides one `.mat` file containing the
%  voxel MNI coordinates (Nvoxels × 3).
%  Using this information, the script:
%    1) Inspects participant files across Groups/Conditions
%    2) Runs frequency-resolved GED for each participant separately
%    3) Saves each FREQ in FREQNESS_Results/CORE/<condition>/<participant>_FREQ.mat
%    4) Optionally combines saved results for the existing group analyses
%
%  BLOCK 0 defines the settings used by both analysis blocks.
%  BLOCK 1 estimates and saves the participant core results.
%  BLOCK 2 runs the optional group analyses and is disabled by default.
%  Run Block 0 before running either analysis block separately.
%
%  Recommended workflow: first run with run_core = true and
%  run_group_analyses = false to keep CORE memory use low. Once CORE finishes,
%  set run_core = false and run_group_analyses = true, keep the same estimation
%  settings, and run again for group analyses.
%
%  NOTE: Block 2 loads the saved FREQ structures into group arrays. Its
%  memory requirements therefore increase with the number of participants.
%  Keep run_group_analyses false when only the participant core files are
%  required. The core functions still accept multi-participant inputs.
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
%  The script processes the example datasets using the following functions:
%
% - FREQNESS_Startup(false)
%    Initializes the environment, sets up the necessary directories, and
%    inspects data files and fetches MNI coordinates from the respective
%    folders. Calling FREQNESS_Startup() without an argument still loads
%    all participant data, preserving the original startup workflow.
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
%  - ENTR = FREQNESS_EntropyLandscape(FREQ)
%    Computes quadratic Rényi entropy (H2) and effective dimensionality
%    (ED) across frequencies and returns them with the frequency axis in
%    one ENTR structure. Only the H2 entropy landscape is visualized.
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
%  - FREQNESS_InducedResponses(FREQ, events, epoch_window, baseline_window, ...)
%    Computes event-related induced oscillatory responses from the broadband
%    component time series. Because event samples and epoch definitions are
%    experiment-specific, this analysis is called separately after the core
%    pipeline; a template is provided at the end of this script.
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
% BLOCK 0 - STARTUP AND ANALYSIS SETTINGS
% ========================================================================

clear
close all
clc

% Setup directories and inspect data files without loading participant arrays
[~, MNI, path_home, folder_data] = FREQNESS_Startup(false);

% ========================================================================
% DEFINE ANALYSIS SETTINGS
%
%  Edit this section before running the pipeline
% ========================================================================

% Select the blocks to run. Group analyses require saved participant FREQ files.
run_core              = true;
run_group_analyses     = false;

% Core results are saved outside FREQNESS_Data, in one folder per condition.
path_results = fullfile(path_home,'FREQNESS_Results');
overwrite_core = false; % set true to replace previously saved core results

% Select the optional analyses for Block 2 (used only when run_group_analyses
% is true). These sections preserve the existing group-level calculations.
run_visualizer        = true;
run_entropy_landscape = true;
run_exponential_dk    = true;
run_freq_gradients    = true;
run_comp_gradients    = true;
run_cross_coupling    = true;
run_backprojection    = true;
run_network_removal   = true;

% ------------------------------------------------------------------------
% 1) FREQNESS_NetworkEstimation
% ------------------------------------------------------------------------
frex   = 1.2:1.2:20*1.2; % frequency vector (Hz) used for the GED analysis
srate  = 250;             % sampling rate (Hz) of your data

% Optional network-estimation settings
network_duration       = [];            % seconds; leave [] to use all data
network_fwidth         = [];            % [] for automatic widths; scalar or one value per frequency otherwise
network_filter         = 'logarithmic'; % automatic mode: Q = 7 to 3.5; 'linear' uses constant Q = 5
network_regularisation = 0.01;          % covariance shrinkage factor
network_ncomps         = 5;             % components retained in evecs, pats and ts
network_bad_segments   = [];            % sample indices excluded from covariance
network_rescale        = false;         % rescale very low-amplitude data
% NOTE: FREQ.evals always retains the complete eigenspectrum, independently
% of network_ncomps. The five-component default applies only to this pipeline.

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
Patterns.path_output = path_results; % condition folders are assigned in Block 2
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

% Store network settings with each participant result. Existing files are
% reused only when these settings and their source-file metadata match.
NETWORK_SETTINGS = struct('frex',frex,'srate',srate, ...
    'duration',network_duration,'fwidth',network_fwidth, ...
    'filter',network_filter,'regularisation',network_regularisation, ...
    'ncomps',network_ncomps,'bad_segments',network_bad_segments, ...
    'rescale',network_rescale);

% Number of experimental Groups / Conditions
nconds = numel(folder_data);

%% ========================================================================
% BLOCK 1 - NETWORK ESTIMATION AND CORE EXPORT (one participant at a time)
% ========================================================================

% Run Block 0 first. The full FREQ structure, including component time series,
% is saved for subsequent analyses. No participant arrays are collected here.
if run_core

    for condi = 1:nconds
        fprintf('\nESTIMATING CONDITION/GROUP: %s\n',folder_data(condi).name);

        files_mat = folder_data(condi).files_mat;
        nsubs = numel(files_mat);
        if nsubs == 0
            continue
        end

        core_path = fullfile(path_results,'CORE',folder_data(condi).name);
        if ~exist(core_path,'dir')
            mkdir(core_path)
        end

        for subi = 1:nsubs
            fprintf('\nPARTICIPANT #%d/%d: %s\n',subi,nsubs,files_mat(subi).name);

            % Identify the source file and its dedicated core output
            source_file = fullfile(files_mat(subi).folder,files_mat(subi).name);
            [~,participant_name] = fileparts(files_mat(subi).name);
            output_file = fullfile(core_path,[participant_name '_FREQ.mat']);
            source_info = dir(source_file);
            if isempty(source_info)
                error('Participant file no longer exists: %s',source_file);
            end
            PARTICIPANT = struct('condition',folder_data(condi).name, ...
                'name',participant_name,'source_file',source_file, ...
                'source_size',folder_data(condi).data_size, ...
                'source_bytes',source_info.bytes,'source_datenum',source_info.datenum);

            % Retain existing results only when their inputs and settings match
            if exist(output_file,'file') && ~overwrite_core
                S = load(output_file,'PARTICIPANT','NETWORK_SETTINGS');
                saved_info = whos('-file',output_file);
                if ~isfield(S,'PARTICIPANT') || ~isfield(S,'NETWORK_SETTINGS') || ...
                        ~ismember('FREQ',{saved_info.name}) || ...
                        ~isequaln(S.PARTICIPANT,PARTICIPANT) || ...
                        ~isequaln(S.NETWORK_SETTINGS,NETWORK_SETTINGS)
                    error(['Existing core file does not match the current inputs or settings: %s. ' ...
                        'Set overwrite_core to true, or choose another path_results.'],output_file);
                end
                fprintf('Core output already exists; skipping participant.\n');
                clear S PARTICIPANT
                continue
            end

            % Load and validate one voxels-by-time data matrix
            S = load(source_file);
            fn = fieldnames(S);
            if numel(fn) ~= 1
                error('File %s must contain exactly one data matrix.',files_mat(subi).name);
            end
            this_data = S.(fn{1});
            clear S
            if ~isnumeric(this_data) || ~isreal(this_data) || isempty(this_data) || ...
                    ~ismatrix(this_data) || any(~isfinite(this_data(:))) || ...
                    ~isequal(size(this_data),PARTICIPANT.source_size)
                error(['File %s must contain one real finite voxels-by-time matrix ' ...
                    'with the dimensions inspected during startup.'],files_mat(subi).name);
            end

            % Match the original startup's double precision. The participant
            % dimension remains a singleton: size(this_data,3) is equal to 1.
            this_data = double(this_data);
            FREQ = FREQNESS_NetworkEstimation(this_data,frex,srate, ...
                'duration',network_duration, ...
                'fwidth',network_fwidth, ...
                'filter',network_filter, ...
                'regularisation',network_regularisation, ...
                'ncomps',network_ncomps, ...
                'bad_segments',network_bad_segments, ...
                'rescale',network_rescale);
            clear this_data

            % Complete the save before replacing the final participant file.
            % Version 7.3 supports large component time-series arrays.
            temporary_file = [tempname(core_path) '.mat'];
            try
                save(temporary_file,'FREQ','PARTICIPANT','NETWORK_SETTINGS','-v7.3');
                [saved,message] = movefile(temporary_file,output_file,'f');
                if ~saved
                    error('Could not finalize core file %s: %s',output_file,message);
                end
            catch exception
                if exist(temporary_file,'file')
                    delete(temporary_file)
                end
                rethrow(exception)
            end
            fprintf('Saved core output: %s\n',output_file);

            % Release this participant before importing the next data file
            clear FREQ PARTICIPANT
        end
    end

end

%% ========================================================================
% BLOCK 2 - OPTIONAL GROUP ANALYSES (combine saved participant FREQ outputs)
% ========================================================================

% Run Block 0 first. Set run_group_analyses to true to enable this block.
% The following sections retain their group outputs in the workspace, as in
% the original pipeline. Unlike Block 1, they require memory for group arrays.
if run_group_analyses

    % Initialize output variables (one per Group / Condition)

    % 1) Core FREQNESS output
    FREQ        = cell(nconds,1);
    % Entropy landscape
    ENTR        = cell(nconds,1);
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
        files_mat = folder_data(condi).files_mat;
        nsubs = numel(files_mat);
        if nsubs == 0
            continue
        end

        % ========================================================================
        % 1) LOAD SAVED FREQNESS CORE OUTPUTS
        % ========================================================================

        core_path = fullfile(path_results,'CORE',folder_data(condi).name);
        for subi = 1:nsubs
            [~,participant_name] = fileparts(files_mat(subi).name);
            input_file = fullfile(core_path,[participant_name '_FREQ.mat']);
            if ~exist(input_file,'file')
                error('Core output missing: %s. Run Block 1 first.',input_file);
            end
            source_file = fullfile(files_mat(subi).folder,files_mat(subi).name);
            source_info = dir(source_file);
            if isempty(source_info)
                error('Participant file no longer exists: %s',source_file);
            end
            PARTICIPANT = struct('condition',folder_data(condi).name, ...
                'name',participant_name,'source_file',source_file, ...
                'source_size',folder_data(condi).data_size, ...
                'source_bytes',source_info.bytes,'source_datenum',source_info.datenum);
            S = load(input_file,'FREQ','PARTICIPANT','NETWORK_SETTINGS');
            if ~isfield(S,'FREQ') || ~isfield(S,'PARTICIPANT') || ...
                    ~isfield(S,'NETWORK_SETTINGS') || ...
                    ~isequaln(S.PARTICIPANT,PARTICIPANT) || ...
                    ~isequaln(S.NETWORK_SETTINGS,NETWORK_SETTINGS)
                error(['Core output does not match the current inputs or settings: %s. ' ...
                    'Run Block 1 with overwrite_core set to true.'],input_file);
            end

            % Allocate once, then restore each participant along the last dimension
            if subi == 1
                FREQ{condi} = S.FREQ;
                nvoxs = size(S.FREQ.evals,1);
                nfrex = numel(S.FREQ.frex);
                ncomps = size(S.FREQ.ts,1);
                ntime = size(S.FREQ.ts,2);
                FREQ{condi}.evals = zeros(nvoxs,nfrex,nsubs);
                FREQ{condi}.evecs = zeros(nvoxs,ncomps,nfrex,nsubs);
                FREQ{condi}.pats  = zeros(nvoxs,ncomps,nfrex,nsubs);
                FREQ{condi}.ts   = zeros(ncomps,ntime,nfrex,nsubs);
                FREQ{condi}.scale_factors = ones(1,nsubs);
            end
            FREQ{condi}.evals(:,:,subi) = S.FREQ.evals;
            FREQ{condi}.evecs(:,:,:,subi) = S.FREQ.evecs;
            FREQ{condi}.pats(:,:,:,subi) = S.FREQ.pats;
            FREQ{condi}.ts(:,:,:,subi) = S.FREQ.ts;
            FREQ{condi}.scale_factors(subi) = S.FREQ.scale_factors;
            clear S PARTICIPANT
        end

        % Keep NIFTI exports from different conditions in separate folders
        Patterns.path_output = fullfile(path_results,'GROUP',folder_data(condi).name);


        % ========================================================================
        % 2) FREQNESS VISUALIZATION
        % ========================================================================

        % Plot network landscape and, when MNI coordinates are available,
        % spatial activation patterns and optional NIFTI images
        if run_visualizer
            if isempty(MNI)
                warning(['MNI coordinates are unavailable or do not match the data. ' ...
                    'Only the network landscape will be visualized.']);
                FREQNESS_Visualizer(FREQ{condi},Landscape,[], ...
                    'plot_all',plot_all,'save_nifti',false)
            else
                FREQNESS_Visualizer(FREQ{condi},Landscape,Patterns, ...
                    'plot_all',plot_all,'save_nifti',save_nifti)
            end
        end

        % NOTE: computation and storage of NIFTI files supported only for 8mm MNI space

        % ========================================================================
        % 3) FREQNESS ENTROPY LANDSCAPE
        % ========================================================================

        % Compute entropy-based indices of the eigenspectrum
        if run_entropy_landscape
            ENTR{condi} = FREQNESS_EntropyLandscape(FREQ{condi});
        end


        % ========================================================================
        % 4) FREQNESS EXPONENTIAL DECAY (Eigenvalue decay across frequency)
        % ========================================================================

        if run_exponential_dk
            [decayCoeff{condi}, goodFit_exp{condi}] = FREQNESS_ExponentialDK(FREQ{condi}, ...
                'which_comp', expdk_which_comp, ...
                'range2fit',  expdk_range2fit,...
                'plot_all',plot_all);
        end

        % ========================================================================
        % 5) FREQNESS FREQUENCY GRADIENTS (Spatial gradients across frequencies)
        % ========================================================================

        if run_freq_gradients && ~isempty(MNI) % run only if requested and coordinates are provided

            [gradCoeff_f{condi}, goodFit_f{condi}] = FREQNESS_FreqGradients(FREQ{condi}, MNI, ...
                'frex2model', freqgrad_frex2model, ...
                'comp2model', freqgrad_comp2model,...
                'plot_all',plot_all);
        end

        % ========================================================================
        % 6) FREQNESS COMPONENT GRADIENTS (Spatial gradients across components)
        % ========================================================================

        if run_comp_gradients && ~isempty(MNI) % run only if requested and coordinates are provided

            [gradCoeff_c{condi}, goodFit_c{condi}] = FREQNESS_CompGradients(FREQ{condi}, MNI, ...
                'freq2model',  compgrad_freq2model, ...
                'comps2model', compgrad_comps2model,...
                'plot_all',plot_all);

        end

        % ========================================================================
        % 7) FREQNESS CROSS FREQUENCY COUPLING (Phase-amplitude coupling)
        % ========================================================================

        if run_cross_coupling
            CFC{condi} = FREQNESS_CrossCoupling(FREQ{condi}, lfo_freq, 'mni', MNI);
        end

        % ========================================================================
        % 8) FREQNESS BACKPROJECTION
        % ========================================================================

        if run_backprojection
            backProj{condi} = FREQNESS_BackProjection(FREQ{condi}, ...
                backproj_freq2project, ...
                'comps2project',backproj_comps2project);
        end

        % ========================================================================
        % 9) FREQNESS NETWORK REMOVAL
        % ========================================================================

        if run_network_removal
            % Match the exact data segment analyzed by FREQNESS_NetworkEstimation.
            % Raw data are imported only when this optional analysis needs them.
            ntime = size(FREQ{condi}.ts,2);
            data2remove = zeros(size(FREQ{condi}.evecs,1),ntime,nsubs);
            for subi = 1:nsubs
                S = load(fullfile(files_mat(subi).folder,files_mat(subi).name));
                fn = fieldnames(S);
                if numel(fn) ~= 1
                    error('File %s must contain exactly one data matrix.',files_mat(subi).name);
                end
                this_data = S.(fn{1});
                clear S
                if ~isnumeric(this_data) || ~isreal(this_data) || ...
                        ~isequal(size(this_data),folder_data(condi).data_size) || ...
                        any(~isfinite(this_data(:)))
                    error('Participant data changed after startup: %s',files_mat(subi).name);
                end
                data2remove(:,:,subi) = this_data(:,1:ntime);
                clear this_data
            end
            [dataClean{condi}, removedActivity{condi}] = FREQNESS_NetworkRemoval( ...
                FREQ{condi}, data2remove, netrem_freq2remove, ...
                'comps2remove',netrem_comps2remove);
            clear data2remove

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

    end

end

%% ========================================================================
%  FINAL REMARKS
% ========================================================================
%
% Block 1 saves one complete FREQ structure per participant under
% FREQNESS_Results/CORE/<condition>. Load an individual file to inspect its
% networks or run subsequent analyses without repeating network estimation.
% Block 2 optionally combines these files and retains group outputs in the
% workspace for comparisons across experimental Conditions or Groups.
%
% For statistical analysis, decide which FREQNESS outputs
% you want to compare and analyze them across experimental
% Conditions, Groups, or any other relevant factor
% based on your study design.
%
% Please check the FREQNESS GitHub repository for new releases.
% Event-related induced responses are available through
% FREQNESS_InducedResponses. They are not run automatically because event
% samples, epoch windows, and baseline windows depend on the experiment.
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

%% ========================================================================
%  OPTIONAL TEMPLATE: EVENT-RELATED INDUCED RESPONSES
% ========================================================================
%
% Load one participant core file and provide that participant's event samples:
%
% load(fullfile(path_results,'CORE','Dataset_1','participant_01_FREQ.mat'),'FREQ')
% events = [...];           % event samples for this participant (1-based)
% epoch_window = [-0.1 3.4]; % seconds relative to event onset
% baseline_window = [-0.1 0];
% IND = FREQNESS_InducedResponses(FREQ,events,epoch_window, ...
%     baseline_window,'which_comp',1,'plot_avg',true);

%%
