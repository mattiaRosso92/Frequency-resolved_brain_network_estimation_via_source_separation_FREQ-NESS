% ========================================================================
%  FREQNESS: EXAMPLE SCRIPT FOR BRAIN NETWORK ESTIMATION
%
%  Please cite the first FREQNESS paper:
%  M. Rosso, G. Fernandez-Rubio, P. E. Keller, E. Brattico, P. Vuust, M. L. Kringelbach, L. Bonetti.
%  FREQ-NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved Brain Networks During Auditory Stimulation.
%  Adv. Sci. 2025, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
%
% ========================================================================
%
%  This script demonstrates the application of FREQNESS 
%  (Frequency-Resolved Brain Network Estimation via Source Separation).
%
%  The example dataset consists of a 3-minute recording with 
%  two conditions:
%  - Resting state ('rest')
%  - Passive listening to an isochronous metronome at 2.4 Hz ('beat')
%  This data is available at the following link (Zenodo repository):
%  https://zenodo.org/records/14922536
%  We recommend to place the data in the folder 'Data_Example' which can be
%  found in the FREQNESS_Toolbox main folder.  
%
%  The script loads the example datasets into the MATLAB workspace and 
%  calls three core functions:
%
% ------------------------------------------------------------------------
%  FUNCTIONS OVERVIEW:
% ------------------------------------------------------------------------
%  - FREQNESS_Startup(path_home) 
%    Initializes the environment. Given `path_home` as input, it sets up 
%    the necessary directories.
%
%  - FREQNESS_NetworkEstimation(...) 
%    Performs Generalized Eigendecomposition (GED) over a user-defined 
%    sample of frequencies, separating frequency-resolved networks.
%
%  - FREQNESS_Visualizer(...) 
%    Takes as input selected outputs from `FREQNESS_NetworkEstimation` 
%    to visualize both the network landscape and brain topographies.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso & Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 24/02/2025
% ========================================================================


%% STARTUP

clear 
close all
clc

% Setup directories
path_home = '/Users/mattiaipem/Desktop/FREQNESS_push/FREQNESS_Toolbox';
addpath(path_home)
FREQNESS_Startup(path_home);


%% PERFORM FREQNESS (ONLY ESSENTIAL INPUTS)

%%% ------------------- USER SETTINGS ------------------- %%%

% Select dataset: (1 = Resting state, 2 = Passive beat listening)
import_data_label = 2;  

% Define frequency vector and sampling rate
testFrex  = [1.2:1.2:1.2*20];  % Frequencies (Hz, ascending order)
testSrate = 250;        % Sampling rate (Hz)


%%% ------------------ COMPUTATION --------------------- %%%

% Import selected data (channel-by-time matrix)
if import_data_label == 1
    load([path_home '/Data_Example/dataExample_rest.mat'])  % Resting-state data
elseif import_data_label == 2
    load([path_home '/Data_Example/dataExample_beat.mat'])  % Beat-listening data
else
    error('Invalid data selection (choose 1 or 2).')
end


% Run FREQNESS network estimation (default parameters)
GED = FREQNESS_NetworkEstimation(testData, testFrex, testSrate);


%% PERFORM FREQNESS (ALTERNATIVE SCENARIO WITH OPTIONAL INPUTS)

% This section demonstrates the same function as above,  
% but with optional settings provided. Any missing arguments  
% will automatically use their default values.  
% NOTE: YOU ONLY NEED TO RUN ONE SECTION: EITHER THIS ONE OR THE PREVIOUS ONE.  


% %%% ------------------- USER SETTINGS ------------------- %%%
% 
% % Select dataset: (1 = Resting state, 2 = Passive beat listening)
% import_data_label = 1;  
% 
% % Define frequency vector and sampling rate
% testFrex  = [7:.3:13];  % Frequencies (Hz, ascending order)
% testSrate = 250;        % Sampling rate (Hz)
% testTime   = 20;     % duration of the signal to analyze, in seconds (equal or less than the input data matrix)
% testFwhm   = .1;     % full width at half the maximum of a Gaussian filter, only set for the lowest frequency 
% testFilter = 'logarithmic'; % defines the filter width for all the remaining frequencies: set the spacing as 'logarithmic' or 'linear'
% testRegularization    = .01;     % shrinkage factor for regularization of covariance matrices (0.01 is recommended)
% testNcomps = 30;   % number of signal components (and associated brain networks) to retain (equal or less than the input data matrix)
% 
% 
% %%% ------------------ COMPUTATION --------------------- %%%
% 
% % Import selected data (channel-by-time matrix)
% if import_data_label == 1
%     load([path_home 'Data_Example/dataExample_rest.mat'])  % Resting-state data
% elseif import_data_label == 2
%     load([path_home 'Data_Example/dataExample_beat.mat'])  % Beat-listening data
% else
%     error('Invalid data selection (choose 1 or 2).')
% end
% 
% 
% GED = FREQNESS_NetworkEstimation(testData, testFrex, testSrate, ...
%                                 'duration', testTime, 'filter', testFilter, 'regularisation', testRegularization, 'ncomps', testNcomps); %% Call with optional parameters


%% FREQNESS VISUALIZATION

% NOTE:
% Outputs #1 (Landscape) and #2 (3D brain plot): supported for any MNI brain space (e.g. 1,2,8mm, etc.)
% Output #3: supported only for 8mm MNI space


%%% ------------------- USER SETTINGS ------------------- %%%

% Settings (visualization of network landscape)
Landscape = [];
Landscape.frex   = testFrex; 
Landscape.ncomps = 10;

% Settings (visualization of network activation patterns)
Patterns = [];
Patterns.frex    = [12 2.4];
Patterns.ncomps  = 1; % set how many top components to visualize
Patterns.path_output =  '/Users/mattiaipem/Desktop'; % set output path to save nifti images

% Load MNI coordinates related to the brain voxels of the current dataset
% NOTE: provide the MNI coordinates in the order matching YOUR OWN DATA!
load('MNI152_8mm_coord_dyi.mat'); %all voxels MNI coordinates
Patterns.MNI_coords = MNI8; %assigning the MNI coordinates of your data for visualization purposes (Outputs #2 and #3)

%%%   if you wish to remove the cerebellum voxels (not included in template of Output #2), please uncomment the following lines; NOTE: this removal works only for 8mm brain %%% 
%     load('cerebellum_coords.mat'); %only cerebellar voxels
%     % Remove cerebellar voxels since they are not covered by the template in Output #2
%     % Please, note that you DO NOT NECESSARILY have to remove them
%     [~, idx_cerebellum] = ismember(MNI8, cerebellum_coords, 'rows');  % find cerebellum indexes in MNI coordinates matrix (all voxels)
%     MNI8(idx_cerebellum~=0,:) = nan; %assigning nans to MNI coordinates matrix
%     Patterns.MNI_coords = MNI8; %assigning the MNI coordinates of your data for visualization purposes (Outputs #2 and #3)


%%% ------------------ COMPUTATION --------------------- %%%

% Plot network landscape and save nifti images
FREQNESS_Visualizer(GED,Landscape,Patterns)


%%


%% FINAL REMARKS

% This demonstration uses data from a single subject.  
% If you have multiple subjects, you can average the outputs  
% across subjects and use the same FREQNESS_Visualizer function  
% to plot the group-level results.  
%
% For statistical analysis, decide which FREQNESS outputs  
% you want to compare and analyze them across experimental  
% conditions, subject groups, or any other relevant factors  
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
%  Please cite the first FREQNESS paper:
%  M. Rosso, G. Fernández-Rubio, P. E. Keller, E. Brattico, P. Vuust, M. L. Kringelbach, L. Bonetti.
%  FREQ-NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved Brain Networks During Auditory Stimulation.
%  Adv. Sci. 2025, 2413195.
%  https://doi.org/10.1002/advs.202413195

%%

