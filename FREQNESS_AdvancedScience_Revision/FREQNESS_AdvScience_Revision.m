%% LEONARDO - REVISION I - ADVANCED SCIENCES

%%% OBS!! I RAN THE GRADIOMETER SOURCE RECONSTRUCTION FOR BOTH RESTING
%%% STATE AND BEAT LISTENING FOR THE DATASET WE ALREADY USED. DONE IT WITH
%%% THE CODE ABOVE, AS IN THE FIRST PLACE

%%% HERE, I'M RUNNING SOURCE RECONSTRUCTION (WITH MAGNETOMETERS) FOR THE
%%% RESTING STATE FOR ANOTHER OF MY DATASETS (LEARNINGBACH2017)
%%% HOWEVER, THIS DATASET WAS PREPROCESSED IN A SLIGHTLY DIFFERENT WAY AND
%%% HAS A SAMPLING RATE OF 150HZ

%%

%% BEAMFORMING for anatomical source reconstruction

% Setting up cluster (parallel computing)
% clusterconfig('scheduler', 'none');  % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 2); % slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') % the function being run on the cluster must be contained in this folder

%%

%%% OBS!! I DECIDED TO RECONSTRUCT ALL TIME-POINTS, THEN WE NEED TO DECIDE HOW MANY DO WE WANT TO TAKE

% User settings
clust_l = 1; % 1 = using cluster of computers (CFIN-MIB, Aarhus University); 0 = running locally
freqq = [];  % frequency range (empty [] for broad band)
% freqq = [0.1 1]; %frequency range (empty [] for broad band)
% freqq = [2 8]; %frequency range (empty [] for broad band)
sensl = 1; %1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers (SUGGESTED 1!)
workingdir_home = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD'; %high-order working directory (a subfolder for each analysis with information about frequency, time and absolute value will be created)
invers = 1; %1-4 = different ways (e.g. mean, t-values, etc.) to aggregate trials and then source reconstruct only one trial; 5 for single trial independent source reconstruction
absl = 0; % 1 = absolute value of sources; 0 = not

% Actual computation
list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/ers*rest*mat');
condslabels = {{'Undefined'},{'Beat_listening'}};
% removing SUBJ0027 for technical reasons in Beat Listening

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');
workingdir = [workingdir_home '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_broadband_invers_' num2str(invers) '_LB2017_SR150Hz'];
mkdir(workingdir)

for subi = 1:length(list) %over subjects
    S = [];
%     if ~isempty(freqq) %if you want to apply the bandpass filter, you need to provide continuous data
%         %             disp(['copying continuous data for subj ' num2str(ii)])
%         %thus pasting it here
%         %             copyfile([list_c(ii).folder '/' list_c(ii).name],[workingdir '/' list_c(ii).name]); %.mat file
%         %             copyfile([list_c(ii).folder '/' list_c(ii).name(1:end-3) 'dat'],[workingdir '/' list_c(ii).name(1:end-3) 'dat']); %.dat file
%         %and assigning the path to the structure S
%         S.norm_megsensors.MEGdata_c = [list(subi).folder '/' list(subi).name(10:end)];
%     end
    %copy-pasting epoched files
    %         disp(['copying epoched data for subj ' num2str(ii)])
    %         copyfile([list(ii).folder '/' list(ii).name],[workingdir '/' list(ii).name]); %.mat file
    %         copyfile([list(ii).folder '/' list(ii).name(1:end-3) 'dat'],[workingdir '/' list(ii).name(1:end-3) 'dat']); %.dat file
    
    S.Aarhus_cluster = clust_l; %1 for parallel computing; 0 for local computation
    
    S.norm_megsensors.zscorel_cov = 1; % 1 for zscore normalization; 0 otherwise
    S.norm_megsensors.workdir = workingdir;
    S.norm_megsensors.MEGdata_e = [list(subi).folder '/' list(subi).name];
    S.norm_megsensors.freq = freqq; %frequency range
    S.norm_megsensors.forward = 'Single Shell'; %forward solution (for now better to stick to 'Single Shell')
    
    S.beamfilters.sensl = sensl; %1 = magnetometers; 2 = gradiometers; 3 = both MEG sensors (mag and grad) (SUGGESTED 3!)
    S.beamfilters.maskfname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz'; % path to brain mask: (e.g. 8mm MNI152-T1: '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz')
    
    %%% CHECK THIS ONE ESPECIALLY!!!!!!! %%%
    S.inversion.znorml = 0; % 1 for inverting MEG data using the zscored normalized one; (SUGGESTED 0 IN BOTH CASES!)
    %                                 0 to normalize the original data with respect to maximum and minimum of the experimental conditions if you have both magnetometers and gradiometers.
    %                                 0 to use original data in the inversion if you have only mag or grad (while e.g. you may have used zscored-data for covariance matrix)
    %
    S.inversion.timef = []; % vector of data-points to be extracted; leave it empty [] for working on the full length of the epoch
    S.inversion.conditions = condslabels{1}; %cell with characters for the labels of the experimental conditions (e.g. {'Old_Correct','New_Correct'})
    S.inversion.bc = []; %extreme time-samples for baseline correction (leave empty [] if you do not want to apply it)
    S.inversion.abs = absl; %1 for absolute values of sources time-series (recommendnded 1!)
    S.inversion.effects = invers;
    
    S.smoothing.spatsmootl = 0; %1 for spatial smoothing; 0 otherwise
    S.smoothing.spat_fwhm = 100; %spatial smoothing fwhm (suggested = 100)
    S.smoothing.tempsmootl = 0; %1 for temporal smoothing; 0 otherwise
    S.smoothing.temp_param = 0.01; %temporal smoothing parameter (suggested = 0.01)
    S.smoothing.tempplot = []; %vector with sources indices to be plotted (original vs temporally smoothed timeseries; e.g. [1 2030 3269]). Leave empty [] for not having any plot.
    
    S.nifti = 0; %1 for plotting nifti images of the reconstructed sources of the experimental conditions
    S.out_name = ['SUBJ_' list(subi).name(18:21)]; %name (character) for output nifti images (conditions name is automatically detected and added)
    
    if clust_l ~= 1 %useful  mainly for begugging purposes
        MEG_SR_Beam_LBPD(S);
    else
        jobid = job2cluster(@MEG_SR_Beam_LBPD,S); %running with parallel computing
    end
end

%%

%%% HERE, I'M RUNNING SOURCE RECONSTRUCTION (WITH MAGNETOMETERS) FOR THE
%%% RESTING STATE FOR ANOTHER OF MY DATASETS (TSA2021)
%%% HERE, SAME PREPROCESSING AS IN THE ORIGINAL PIPELINE THAT WE USED, BUT
%%% IN TSA2021 WE HAVE BOTH YOUNG AND OLDER ADULTS (A BIT LESS THAN 40 IN
%%% BOTH GROUPS)

%% (FAKE) EPOCHS FOR RESTING STATE (USEFUL FOR COMPUTER-SCIENCE PURPOSES)

prefix_tobeadded = 'e_long'; %adds this prefix to epoched files
spm_list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuos/Resting/spmeeg*rest*.mat');
for ii = 1:length(spm_list) %over .mat files
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]); %load spm_list .mat files
    D = D.montage('switch',0);
    dummy = D.fname; %OBS! D.fname does not work, so we need to use a 'dummy' variable instead
    
    end_time_sec = D.time(end)-10; %total time minus 10 seconds (just to avoid any potential issues with boundaries, etc.)
    start_time_sec = 30; %30 seconds

    % matrix with the information about the beginning of the epoch in seconds
    trl_sam = zeros(1,3); %prepare the samples matrix with 0's in all its cells
    % matrix with the information about the beginning of the epoch in points (use sample frequency for the conversion
    trl_sec = zeros(1,3); %prepare the seconds matrix with 0's in all its cells
    
    %deftrig(k,1) = 0.012 + trigcor(k,1); %adding a 0.012 seconds delay to the triggers sent during the experiment (this delay was due to technical reasons related to the stimuli)
    trl_sec(1,1) = start_time_sec; % extracting a random trigger time
    %remove 1000ms of baseline
    trl_sec(1,2) = end_time_sec; %end time-window epoch in seconds
    trl_sec(1,3) = trl_sec(1,2) - trl_sec(1,1); %range time-windows in seconds
    trl_sam(1,1) = round(trl_sec(1,1) * 250) + 1; %beginning time-window epoch in samples %250Hz per second
    trl_sam(1,2) = round(trl_sec(1,2) * 250) + 1; %end time-window epoch in samples
    trl_sam(1,3) = -25; %sample before the onset of the stimulus (corresponds to 0.100ms)
    
    %creates the epochinfo structure that is required for the source reconstruction later
    epochinfo.trl = trl_sam;
    epochinfo.time_continuous = D.time;
    %switch the montage to 0 because for some reason OSL people prefer to do the epoching with the not denoised data
    D = D.montage('switch',0);
    %build structure for spm_eeg_epochs
    S = [];
    S.D = D;
    S.trl = trl_sam;
    S.prefix = prefix_tobeadded;
    D = spm_eeg_epochs(S);
    %store the epochinfo structure inside the D object
    D.epochinfo = epochinfo;
    D = D.montage('switch',1);
    D.epochinfo.conditionlabels = D.conditions; %to add for later use in the source reconstruction
    D.save();
    disp(['Subject ' num2str(ii) ' is done'])
end


%% ACTUAL SOURCE RECONSTRUCTION

%%% OBS!! I DECIDED TO RECONSTRUCT ALL TIME-POINTS, THEN WE NEED TO DECIDE HOW MANY DO WE WANT TO TAKE

% User settings
clust_l = 1; % 1 = using cluster of computers (CFIN-MIB, Aarhus University); 0 = running locally
freqq = [];  % frequency range (empty [] for broad band)
% freqq = [0.1 1]; %frequency range (empty [] for broad band)
% freqq = [2 8]; %frequency range (empty [] for broad band)
sensl = 1; %1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers (SUGGESTED 1!)
workingdir_home = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD'; %high-order working directory (a subfolder for each analysis with information about frequency, time and absolute value will be created)
invers = 1; %1-4 = different ways (e.g. mean, t-values, etc.) to aggregate trials and then source reconstruct only one trial; 5 for single trial independent source reconstruction
absl = 0; % 1 = absolute value of sources; 0 = not

% Actual computation
list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/GED_TSA2021/Continuos/Resting/e_longspmeeg*rest*.mat');
condslabels = {{'Undefined'}};
% removing SUBJ0027 for technical reasons in Beat Listening

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');
workingdir = [workingdir_home '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_broadband_invers_' num2str(invers) '_TSA2021_SR250Hz'];
mkdir(workingdir)

for subi = 1:length(list) %over subjects
    S = [];
%     if ~isempty(freqq) %if you want to apply the bandpass filter, you need to provide continuous data
%         %             disp(['copying continuous data for subj ' num2str(ii)])
%         %thus pasting it here
%         %             copyfile([list_c(ii).folder '/' list_c(ii).name],[workingdir '/' list_c(ii).name]); %.mat file
%         %             copyfile([list_c(ii).folder '/' list_c(ii).name(1:end-3) 'dat'],[workingdir '/' list_c(ii).name(1:end-3) 'dat']); %.dat file
%         %and assigning the path to the structure S
%         S.norm_megsensors.MEGdata_c = [list(subi).folder '/' list(subi).name(10:end)];
%     end
    %copy-pasting epoched files
    %         disp(['copying epoched data for subj ' num2str(ii)])
    %         copyfile([list(ii).folder '/' list(ii).name],[workingdir '/' list(ii).name]); %.mat file
    %         copyfile([list(ii).folder '/' list(ii).name(1:end-3) 'dat'],[workingdir '/' list(ii).name(1:end-3) 'dat']); %.dat file
    
    S.Aarhus_cluster = clust_l; %1 for parallel computing; 0 for local computation
    
    S.norm_megsensors.zscorel_cov = 1; % 1 for zscore normalization; 0 otherwise
    S.norm_megsensors.workdir = workingdir;
    S.norm_megsensors.MEGdata_e = [list(subi).folder '/' list(subi).name];
    S.norm_megsensors.freq = freqq; %frequency range
    S.norm_megsensors.forward = 'Single Shell'; %forward solution (for now better to stick to 'Single Shell')
    
    S.beamfilters.sensl = sensl; %1 = magnetometers; 2 = gradiometers; 3 = both MEG sensors (mag and grad) (SUGGESTED 3!)
    S.beamfilters.maskfname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz'; % path to brain mask: (e.g. 8mm MNI152-T1: '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz')
    
    %%% CHECK THIS ONE ESPECIALLY!!!!!!! %%%
    S.inversion.znorml = 0; % 1 for inverting MEG data using the zscored normalized one; (SUGGESTED 0 IN BOTH CASES!)
    %                                 0 to normalize the original data with respect to maximum and minimum of the experimental conditions if you have both magnetometers and gradiometers.
    %                                 0 to use original data in the inversion if you have only mag or grad (while e.g. you may have used zscored-data for covariance matrix)
    %
    S.inversion.timef = []; % vector of data-points to be extracted; leave it empty [] for working on the full length of the epoch
    S.inversion.conditions = condslabels{1}; %cell with characters for the labels of the experimental conditions (e.g. {'Old_Correct','New_Correct'})
    S.inversion.bc = []; %extreme time-samples for baseline correction (leave empty [] if you do not want to apply it)
    S.inversion.abs = absl; %1 for absolute values of sources time-series (recommendnded 1!)
    S.inversion.effects = invers;
    
    S.smoothing.spatsmootl = 0; %1 for spatial smoothing; 0 otherwise
    S.smoothing.spat_fwhm = 100; %spatial smoothing fwhm (suggested = 100)
    S.smoothing.tempsmootl = 0; %1 for temporal smoothing; 0 otherwise
    S.smoothing.temp_param = 0.01; %temporal smoothing parameter (suggested = 0.01)
    S.smoothing.tempplot = []; %vector with sources indices to be plotted (original vs temporally smoothed timeseries; e.g. [1 2030 3269]). Leave empty [] for not having any plot.
    
    S.nifti = 0; %1 for plotting nifti images of the reconstructed sources of the experimental conditions
    S.out_name = ['SUBJ_' list(subi).name(18:21)]; %name (character) for output nifti images (conditions name is automatically detected and added)
    
    if clust_l ~= 1 %useful  mainly for begugging purposes
        MEG_SR_Beam_LBPD(S);
    else
        jobid = job2cluster(@MEG_SR_Beam_LBPD,S); %running with parallel computing
    end
end

%%


%%

%% MATTIA - REVISION I - ADVANCED SCIENCES



%% Replicate GENERALISED EIGENVECTOR DECOMPOSITION (GED) on LB2017 and TSA2021 datasets (lines 560-...)


% Setting up cluster (parallel computing)
% clusterconfig('scheduler', 'none');  % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 8); % slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia')


% Choose whether to test whole sample of frequencies (set to 1), or just
% the target stimulation frequency (set to 0)
allfrex_flag = 1; 

% NOTE: with the current approach, we anchor the width of the filter to the
% first value based on Rosso et al., (2021a,2023a), and then separately
% compute the other widths as logarithmic functions of the frequency bins 

% Settings
S = struct();     % initialize input structure to submit to cluster
S.time = 304;     % in seconds (to match the original FREQNESS dataset)
S.targetfrex = 2.439; % stimulation frequency: hard-coded 1 / 0.410 s (range is 0.400 and 0.420 but 0.410 is the mode) 
S.fwhm  = .35;        % filter width
S.shr   = 0.01;       % shrinking proportion for regularization; value between 0 and 1 (until 0.001 it works for restoring the full rank of the matrix)
S.ncomps2keep = 100;  % how many components to keep (given the huge amount of components = N voxels)
S.comps2see   = [1 2 3]; % list components of interest, of which you want to save nifti image
S.path_script    = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia'; % path of the present script
S.path_source    = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1_LB2017_SR150Hz/'; % '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1_TSA2021_SR250Hz/';%
S.path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/LB2017/'; % ; '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/TSA2021/'; %
S.name_folder = ['GEDoutput_frex_' num2str(S.targetfrex) 'Hz/']; 
if contains(S.path_GEDoutput, 'LB2017') % NOTE: select the right srate for the right dataset
    S.srate = 150;
else
    S.srate = 250; % 250 Hz is the standard after downsampling
end
% The following defines whether and how we are going to shuffle the sourde data 
% in permutation testing
% 0 = no permutation applied 
% 1 = randomize sources idxs (rows)
% 2 = randomize source idx indebendently per timepoint (rows, for every column)
% 3 = randomize timepoints (columns)
% 4 = randomize all together
S.permu_flag = 0;
S.permu_name = {'randLabel', 'randLabelPointWise', 'randTime', 'randAll'}; % to append to files/folder names for saving the data
% Double - check your permutation strategy
if S.permu_flag > 0
	warning('You have selected a permutation strategy. Double-check that you have selected the correct one, as submitting the job to the cluster might result in a very long computation time. Press SPACE BAR to continue')
    pause()
end

% Set all frequencies to test, and associated filter widths
if allfrex_flag == 1 

    % Set with respect to stimulation frequency
    nfrex_above = 80; % tested for even numbers
    nfrex_below = 6;   % tested for even numbers
    % Compute all frequencies
    frex_above = linspace(S.targetfrex, S.targetfrex * (nfrex_above/2), nfrex_above);
    frex_below = linspace(S.targetfrex / 2, S.targetfrex / (2 * nfrex_below), nfrex_below);
    frex_all = [frex_below(end:-1:1), frex_above];
    % compute all filter widths
    fwhm_above = logspace(log10(S.fwhm), log10(S.fwhm * (nfrex_above)), nfrex_above);
    fwhm_below = logspace(log10(S.fwhm), log10(S.fwhm / (nfrex_below)), nfrex_below+1);
    fwhm_all   = [fwhm_below(end:-1:2) fwhm_above];

    % Visualize frewquencies and filters
    figure
    plot(frex_all, fwhm_all, 's-')
    xlabel('Center frequencies (Hz)')
    ylabel('Filter FWHM (Hz)')
    title('Filter width as a function of center frequency')


    % Run GED with parallel computing for a set of frequencies
    nfrex = length(frex_all);
    for frexi = 1:nfrex
        S.targetfrex  = frex_all(frexi); % overwrite fields of input structure
        S.fwhm        = fwhm_all(frexi);
        S.name_folder = ['GEDoutput_frex_' num2str(S.targetfrex) 'Hz/'];
        % Running with parallel computing
        jobid = job2cluster(@FREQNESS_AdvScience_GED_SingleSubjects,S); % always give a structure as input
        
        disp(['Submitting GED for frequency = ' num2str(S.targetfrex) ' Hz and width = ' num2str(S.fwhm) ' Hz on cluster'])
        
    end

else
    
    % Run GED with parallel computing for single stimulation frequency
    jobid = job2cluster(@FREQNESS_AdvScience_GED_SingleSubjects,S); % always give a structure as input

end


%% GED assessment (TSA2021 and LB2017)

% This blocks enables the visual assessment of GED, for individual
% participants and grand-averages, and produces NIFTI images
% - Power spectrum
% - Eigenspectrum
% - Activations maps (cleaned and saved as NIFTI for visualization in the next block)

clear all
close all

path_dataset = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/LB2017/';

% Subject IDs of TSA2021 dataset, divided by age group
idx_young = {'002.mat','005.mat','009.mat','013.mat','014.mat','015.mat','017.mat','019.mat','021.mat','022.mat','024.mat','025.mat','037.mat','039.mat','040.mat','041.mat','044.mat','046.mat','050.mat','053.mat','054.mat','055.mat','057.mat','060.mat','063.mat','064.mat','065.mat','067.mat','069.mat','070.mat','071.mat','072.mat','073.mat','074.mat','075.mat','077.mat','078.mat'}; %young
%idx_old   = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,35,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76]; %elderly
% Subjects IDs to be removed from LB2017 dataset
idx_2remove = {'002.mat','003.mat'};

%GED list
list_up = dir([path_dataset 'GED*']);

% Frequency to inspect
% frex2see = 2.439;   %  stimfrex;
freq_vect_idx = 1:length(list_up);
flag_figures = 0;

% Select 'UniversitA della strada' approach
streetuni_flag = 1; % when set to 1, processes activation patterns for interpretation

% General settings
if contains(path_dataset,'LB2017') % this dataset was downsampled to 150 Hz
    srate = 150;
else
    srate = 250;
end
ncomps = 30;         % components to carry from previous block, to produce scree plots
comps2see = [1 2 3]; % expected top components, to assess in more detail


for frexi = 1:length(freq_vect_idx)
    
    path_frex = [list_up(freq_vect_idx(frexi)).folder '/' list_up(freq_vect_idx(frexi)).name];
    path_GEDoutput = path_frex;
    list = dir([path_frex '/SUBJ*.mat']);
    
    % Filter the young ones for TSA2021 dataset and the bad ones for LB2017
    if contains(path_dataset,'TSA2021')
        young_finder = {list(:).name};
        young_filter = contains(young_finder, idx_young);
        list = list(young_filter);
    elseif contains(path_dataset,'LB2017')
        bad_finder = {list(:).name};
        bad_filter = ~contains(bad_finder, idx_2remove);  
        list = list(bad_filter);
    end
    
    nsubs = length(list);
    for subi = 1:nsubs % over subjects
        
        load([list(subi).folder '/' list(subi).name],'evals','GEDts','GEDmap')
        npnts  = size(GEDts,2);       % number of time points in the signal
        nevals = length(evals);
        
        % Initialize matrix for all participants, once upon first iteration
        if subi == 1
            GEDpow_all = zeros(ncomps,npnts,nsubs);
            evals_all  = zeros(nevals,nsubs);
            GEDts_all  = zeros(ncomps,npnts,nsubs);
            GEDmap_all = zeros(ncomps,nevals,nsubs);
            
            % Settings for power spectrum
            tmax = npnts/srate; % duration, in seconds
            frexres = 1/tmax;   % Rayleigh frequency 1 / T(in seconds)
            % FFT parameters
            nfft = ceil( srate/frexres );  % length sufficiently high to guarantee resolution at denominator (= pnts)
            hz   = linspace(0,srate,nfft); % vector of frequencies
        end
        evals_all(:,subi)   = (evals.*100)./sum(evals); % assign and normalize, for plotting
        GEDts_all(:,:,subi) = GEDts(1:ncomps,:);        % assign for storing
        
        % Compute power spectrum
        for compi = comps2see
            % Demeaning
            GEDts_all(compi,:,subi) = GEDts_all(compi,:,subi) - mean(GEDts_all(compi,:,subi),2); % if not performed, some subjects showed prominent 0-Hz frequency component
            % Power spectrum
            GEDpow_all(compi,:,subi) = abs(fft(GEDts_all(compi,:,subi)',nfft,1)/npnts).^2;
        end
        
        % Assign activation patterns
        for compi = comps2see
            GEDmap_all(compi,:,subi) = GEDmap(compi,:);
        end
        
        disp(['GED assessment of freq idx #' num2str(frexi) ' - participant #' num2str(subi)])
    end
    
    % Compute grand-averages over subjects
    GEDts_avg =  squeeze( mean(GEDts_all,3) );
    GEDpow_avg = squeeze( mean(GEDpow_all,3) );
    evals_avg  = squeeze( mean(evals_all,2) );
    % Process activation pattern, if flag is enabled
    if streetuni_flag ==1
        GEDmap_avg = squeeze( mean(abs(GEDmap_all),3) );
    else
        GEDmap_avg = squeeze( mean(GEDmap_all,3) );
    end
    GEDmap_avg = GEDmap_avg' ./ (max(GEDmap_avg'));  % normalize to 1
    GEDmap_avg = GEDmap_avg';
    
    
    % Calculate distance of every participant's activation pattern from average, per component
    % (check the assumption that network A = COMP#1 and netowrk B = COMP#2)
    z_thresh = 2.3;  %~.01  % 1.96;
    idx_toofar = zeros(length(comps2see),nsubs);
    %Remove 'bad' matrices based on distance
    [GEDmap_dist,GEDmap_z] = deal ( zeros(ncomps,nsubs) );
    for compi = comps2see
        for subi = 1:nsubs
            %Compute distance
            GEDmap_dist(compi,subi) = abs( sqrt(trace(GEDmap_all(compi,:,subi)'*GEDmap_avg(compi,:))) );  % i take abs becasue I get complex numbers with one component = 0; it does not change the result
            %Normalize distance (temporary variable)
            GEDmap_z(compi,subi) = ( GEDmap_dist(compi,subi) - mean(GEDmap_dist(compi,:)) ) / std(GEDmap_dist(compi,:));
            
            
            
            %         %remove distances just for visualization, using logical indexing
            %         GEDmap_z(idx_toofar) = NaN;
        end
        % find indexes of outliers
        idx_toofar(compi,:) = abs(GEDmap_z(compi,:)) > z_thresh;
        disp( [ num2str(sum(idx_toofar(compi,:))) ' outliers identified for component #' num2str(compi) ] )
    end
    
    
    outliers_flag = 0;  % when set to 1, excludes activation patterns to deviate from average within component
    
    % plotting weights in the brain (nifti images)
    maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); % getting the mask for creating the figure
    for compi = comps2see % over components of interest
        
        % Start with assignment to temporary variable, to avoid overwriting later on
        if outliers_flag == 1
            GEDmap2see = squeeze( mean(abs(GEDmap_all(compi,:,~idx_toofar(compi,:))),3) ); % exclude outliers, if enabled
            GEDmap2see = GEDmap2see'./(max(GEDmap2see')); % and normalize
            GEDmap2see = GEDmap2see';
        else
            GEDmap2see = GEDmap_avg(compi, :);
        end
        
        %thresholding activation patterns according to mean + 1 standard deviation
        thresh = mean(GEDmap2see) + std(GEDmap2see);
        GEDmap2see(GEDmap2see<thresh) = 0;
        
        % using nifti toolbox
        SO = GEDmap2see';
        % building nifti image
        SS = size(maskk.img);
        dumimg = zeros(SS(1),SS(2),SS(3),1); % no time here so size of 4D = 1
        for sourci = 1:size(SO,1) % over brain sources
            dumm = find(maskk.img == sourci); % finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
            [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); % getting subscript in 3D from index
            dumimg(i1,i2,i3,:) = SO(sourci,:); % storing values for all time-points in the image matrix
        end
        nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
        nii.img = dumimg; %storing matrix within image structure
        nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
        
        % Save output
        disp(['saving nifti image - component ' num2str(compi) ' frequency ' num2str(list_up(freq_vect_idx(frexi)).name(16:end))])
        %     if adj_compentr == 1
        %         save_nii(nii,['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1/Test_GED/CompEntrAdj_SubjsAverage_Comp_' num2str(frexi) '_Var_' num2str(round(evals_avg(compi))) '.nii']); %printing image
        %     else
        if outliers_flag ==1
            save_nii(nii,[path_GEDoutput '/SubjsAverage_Comp_' num2str(compi) '_Var_' num2str(round(evals_avg(compi))) 'frex_' num2str(list_up(freq_vect_idx(frexi)).name(16:end)) '_NO_outliers.nii']); % printing image
        else
            save_nii(nii,[path_GEDoutput '/SubjsAverage_Comp_' num2str(compi) '_Var_' num2str(round(evals_avg(compi))) 'frex_' num2str(list_up(freq_vect_idx(frexi)).name(16:end)) '.nii']); % printing image
        end
    end
end




%Visualize distances
if flag_figures == 1
    figure(1000), clf
    for compi = comps2see

        subplot(3,1,compi), hold on
        plot(GEDmap_z(compi,:),'bs-','linew',2,'markerfacecolor','w')
        legend('GED maps z-scores')
        plot(xlim, [z_thresh, z_thresh],'--');
        plot(xlim, -1*[z_thresh, z_thresh],'--');
        title('Distances before removal')
        %     %visualize difference after removal
        %     subplot(212)
        %     plot(GEDmap_z,'bs-','linew',2,'markerfacecolor','w')
        %     legend('GED maps z-scores' , 'CovR z-scores')
        %     yline(z_thresh,'--');
        %     yline(-z_thresh,'--');
        %     title('Distances after removal')


    end
end

% %Remove actual outlier matrices, using logical indexing
% GEDmap_all(compi,:,idx_toofar) = NaN;
% %Re-compute grand-averages (excluding NaNs)
% GEDmap_avg = squeeze(mean(covS, 3 , 'omitnan')); 

% Visualization
if flag_figures == 1

    maxfrex2see = 30;
    targetfrex = 2.439;
    [peakfrex, idx_peak] = deal( nan(length(comps2see),nsubs) );
    idx_stim = dsearchn(hz',targetfrex);
    % Individual participants
    for subi = 1:nsubs
        figure(subi),clf
        for compi = comps2see
            % Store peaks and their indexes, for check and eventual further analyses
            [peakfrex(compi,subi), idx_peak(compi,subi)] = max(smooth(GEDpow_all(compi,:,subi),1));
            % Spectrum
            subplot(3,2,(compi-1)*2+1), hold on
            plot(hz,smooth(GEDpow_all(compi,:,subi),1),'k','linew',1.7)
            set(gca,'xlim',[0 maxfrex2see])
            plot(hz(idx_peak(compi,subi)), peakfrex(compi,subi),'*','markerfacecolor','r','linew',5) % empirical peak
            plot(hz(idx_stim), peakfrex(compi,subi),'*','markerfacecolor','b','linew',5)         % expected stimulation frequency
            xlabel('Frequency (Hz)'), ylabel('Power (muV^2)') %or 'Signal-to-noise ratio (%)'
            title(['Power spectrum - C #' num2str(compi)])
            axis square
            % Eigenspectrum
            subplot(3,2,(compi-1)*2+2);     hold on
            plot(evals_all(:,subi),'s-','markerfacecolor','k','linew',2)       % eigenspectrum
            plot(compi,evals_all(compi,subi),'*','markerfacecolor','r','linew',5) % current component
            xlim([1 ncomps])
            xlabel('Component'), ylabel('Eigenvalue')
            title(['Eigenspectrum - C #' num2str(compi)])
            axis square
        end
    end
    % Grand-average
    figure(100),clf
    [peakfrex_avg, idx_peak_avg] = deal( nan(length(comps2see)) );
    for compi = comps2see
        
        [peakfrex_avg(compi), idx_peak_avg(compi)] = max(smooth(GEDpow_avg(compi,:),1));
        % Spectrum
        subplot(3,2,(compi-1)*2+1), hold on
        plot(hz,smooth(GEDpow_avg(compi,:),1),'k','linew',1.7)
        set(gca,'xlim',[0 maxfrex2see])
        plot(hz(idx_peak_avg(compi)), peakfrex_avg(compi),'*','markerfacecolor','r','linew',5) % empirical peak
        plot(hz(idx_stim), peakfrex_avg(compi),'*','markerfacecolor','b','linew',5)         % expected stimulation frequency
        xlabel('Frequency (Hz)'), ylabel('Power (muV^2)') %or 'Signal-to-noise ratio (%)'
        title(['Power spectrum - C #' num2str(compi)])
        axis square
        % Eigenspectrum
        subplot(3,2,(compi-1)*2+2);     hold on
        plot(evals_avg,'s-','markerfacecolor','k','linew',2)       % eigenspectrum
        plot(compi,evals_avg(compi),'*','markerfacecolor','r','linew',5) % current component
        xlim([1 ncomps])
        xlabel('Component'), ylabel('Eigenvalue')
        title(['Eigenspectrum - C #' num2str(compi)])
        axis square
        
    end
end

%%


%% Assess frequency-specificity of GED

clear all
close all
clc


% Subjects IDs of TSA2021 dataset, divided by age group
idx_young = {'002.mat','005.mat','009.mat','013.mat','014.mat','015.mat','017.mat','019.mat','021.mat','022.mat','024.mat','025.mat','037.mat','039.mat','040.mat','041.mat','044.mat','046.mat','050.mat','053.mat','054.mat','055.mat','057.mat','060.mat','063.mat','064.mat','065.mat','067.mat','069.mat','070.mat','071.mat','072.mat','073.mat','074.mat','075.mat','077.mat','078.mat'}; %young
%idx_old   = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,35,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76]; %elderly
% Subjects IDs to be removed from LB2017 dataset
idx_2remove = {'002.mat','003.mat'};


% Select components of interest
comps2see = [1:10];
stimfrex = 2.439;

% Set home path
path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/LB2017/';

% Get all_frex from home path
list_GED = dir(path_GEDoutput);
list_GED = list_GED( [list_GED.isdir] );
list_GED = list_GED(startsWith({list_GED.name}, 'GED')); % added since GED and PCA coexist in the same folder


% Initialize vectors for storing frequencies
nfrex = length(list_GED);
all_frex = zeros(nfrex,1);
% Initialize topEvals
topEvals_avg = zeros(length(comps2see),nfrex);% avg of top components across subjects, per frequency

% Import and process top eigenvalues, per frequency
for frexi = 1:nfrex
    
    % Get folder name
    this_folder = [path_GEDoutput, list_GED(frexi).name '/'];
    % Retrieve frequency
    this_frex = regexp(this_folder, '([\d.]+)Hz', 'tokens');
    % Assign frequency to vector
    all_frex(frexi) = str2double(this_frex{1,1});
    
    % All subjects within frequency frexi
    fulllist = dir([this_folder '/SUBJ*.mat']);
    
    % Sublist: Filter the young ones for TSA2021 dataset and the bad ones for LB2017    
    if contains(path_GEDoutput, 'TSA2021')             
        young_finder = {fulllist(:).name};
        young_filter = contains(young_finder, idx_young);
        list_subs = fulllist(young_filter);
 
    elseif contains(path_GEDoutput,'LB2017')
        bad_finder = {fulllist(:).name};
        bad_filter = ~contains(bad_finder, idx_2remove);  
        list_subs = fulllist(bad_filter);
    end

    nsubs = length(list_subs);
    % Process evals from subjects
    for subi = 1:nsubs % over subjects
        
        % Get eigenvalues (already sorted and normalized)
        load([this_folder, list_subs(subi).name],'evals')
        nevals = length(evals);
        
        if subi == 1 && frexi == 1 
            % Initialize matrix for all participants, once upon first iteration
            evals_all = zeros(nevals,nsubs,nfrex); % for Rest and Listen conditions
        end
        
        evals_all(:,subi,frexi) = evals; % assign for computing average and for exporting
      
        
    end
    
    % Compute grand-averages over subjects (just for visualization)
    topEvals_avg(:,frexi) = squeeze( mean(evals_all(comps2see,:,frexi),2) );
    
    disp(['Processing evals of frex ' this_frex{1,1} 'Hz'])
    
end
% Re-sort vector of frequencies (strings follow alphabetic order)
[all_frex, idx_sort] = sort(all_frex, 'ascend');
topEvals_avg = topEvals_avg(:,idx_sort);

% Visualize Rest 
figure(300), clf
subplot(311), hold on
plot(all_frex,topEvals_avg,'s-','markerfacecolor','k','linew',2)
line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
xticks(round(all_frex,1)), xtickangle(45)
%xlim([0 nfrex+1])
xlabel('Frequencies (Hz)')
ylabel('Explained variance')
legendEntries = cell(length(comps2see),1);
for legi = 1:length(legendEntries)
    legendEntries{legi} = ['Component #' num2str(comps2see(legi))];
end
legend(legendEntries, 'Location', 'Best');
title(['Top ' num2str(length(comps2see)) ' components'])
subplot(312), hold on
plot(all_frex,sum(topEvals_avg,1),'s-','color','k','markerfacecolor','k','linew',2)
line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
xticks(all_frex), xtickangle(45)
%xlim([0 nfrex+1])
xlabel('Frequencies (Hz)')
ylabel('Explained variance')
title(['Sum of top ' num2str(length(comps2see)) ' components'])
subplot(313), hold on
plot(all_frex,std(topEvals_avg,0,1),'s-','color','r','markerfacecolor','k','linew',2)
line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
xticks(all_frex), xtickangle(45)
%xlim([0 nfrex+1])
% ylim([0 50])
xlabel('Frequencies (Hz)')
ylabel('Explained variance')
title(['STD of top ' num2str(length(comps2see)) ' components'])




%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Replicate GENERALISED EIGENVECTOR DECOMPOSITION (GED) with GRADIOMETERS on FREQNESS dataset (lines 560-...)

% Setting up cluster (parallel computing)
% clusterconfig('scheduler', 'none');  % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 8); % slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia')

%%

which_cond = 1;

% This will determin which directories you import data from, and to which
% directories you save your output
cond_suffix = {'_rest'; '_beat'};

% Choose whether to test whole sample of frequencies (set to 1), or just
% the target stimulation frequency (set to 0)
allfrex_flag = 1; 

% NOTE: with the current approach, we anchor the width of the filter to the
% first value based on Rosso et al., (2021a,2023a), and then separately
% compute the other widths as logarithmic functions of the frequency bins 

% Settings
S = struct();     % initialize input structure to submit to cluster
S.targetfrex = 2.439; % stimulation frequency: 
S.fwhm  = .35;        % filter width
S.shr   = 0.01;       % shrinking proportion for regularization; value between 0 and 1 (until 0.001 it works for restoring the full rank of the matrix)
S.srate = 250;        % recorded at 1000 Hz, downsampled by 4
S.ncomps2keep = 100;  % how many components to keep (given the huge amount of components = N voxels)
S.comps2see   = [1 2 3]; % list components of interest, of which you want to save nifti image
S.path_script    = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia'; % path of the present script
S.path_source    = ['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_2_freq_broadband_invers_1' cond_suffix{which_cond} '/'];
S.path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Gradiometers/';
S.name_folder = ['GEDoutput_frex_' num2str(S.targetfrex) 'Hz' cond_suffix{which_cond} '/']; 
% The following defines whether and how we are going to shuffle the sourde data 
% in permutation testing
% 0 = no permutation applied 
% 1 = randomize sources idxs (rows)
% 2 = randomize source idx indebendently per timepoint (rows, for every column)
% 3 = randomize timepoints (columns)
% 4 = randomize all together
S.permu_flag = 0;
S.permu_name = {'randLabel', 'randLabelPointWise', 'randTime', 'randAll'}; % to append to files/folder names for saving the data
% Double - check your permutation strategy
if S.permu_flag > 0
	warning('You have selected a permutation strategy. Double-check that you have selected the correct one, as submitting the job to the cluster might result in a very long computation time. Press SPACE BAR to continue')
    pause()
end

% Set all frequencies to test, and associated filter widths
if allfrex_flag == 1 

    % Set with respect to stimulation frequency
    nfrex_above = 80; % tested for even numbers
    nfrex_below = 6;   % tested for even numbers
    % Compute all frequencies
    frex_above = linspace(S.targetfrex, S.targetfrex * (nfrex_above/2), nfrex_above);
    frex_below = linspace(S.targetfrex / 2, S.targetfrex / (2 * nfrex_below), nfrex_below);
    frex_all = [frex_below(end:-1:1), frex_above];
    % compute all filter widths
    fwhm_above = logspace(log10(S.fwhm), log10(S.fwhm * (nfrex_above)), nfrex_above);
    fwhm_below = logspace(log10(S.fwhm), log10(S.fwhm / (nfrex_below)), nfrex_below+1);
    fwhm_all   = [fwhm_below(end:-1:2) fwhm_above];

    % Visualize frewquencies and filters
    figure
    plot(frex_all, fwhm_all, 's-')
    xlabel('Center frequencies (Hz)')
    ylabel('Filter FWHM (Hz)')
    title('Filter width as a function of center frequency')


    % Run GED with parallel computing for a set of frequencies
    nfrex = length(frex_all);
    for frexi = 1:nfrex
        S.targetfrex  = frex_all(frexi); % overwrite fields of input structure
        S.fwhm        = fwhm_all(frexi);
        S.name_folder = ['GEDoutput_frex_' num2str(S.targetfrex) 'Hz' cond_suffix{which_cond} '/'];
        % Running with parallel computing
        jobid = job2cluster(@GED_SingleSubjects_v3,S); % always give a structure as input
        
        disp(['Submitting GED for frequency = ' num2str(S.targetfrex) ' Hz and width = ' num2str(S.fwhm) ' Hz on cluster'])
        
    end

else
    
    % Run GED with parallel computing for single stimulation frequency
    jobid = job2cluster(@GED_SingleSubjects_v3,S); % always give a structure as input

end


%%

%% Replicate GENERALISED EIGENVECTOR DECOMPOSITION (GED) with MAGNETOMETERS AND SOURCE LEAKAGE CORRECTION on FREQNESS dataset (lines 560-...)

% Setting up cluster (parallel computing)
% clusterconfig('scheduler', 'none');  % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 8); % slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia')
%%
which_cond = 2;
% This will determin which directories you import data from, and to which
% directories you save your output
cond_suffix = {'_rest'; '_beat'};

% Choose whether to test whole sample of frequencies (set to 1), or just
% the target stimulation frequency (set to 0)
allfrex_flag = 1; 

% NOTE: with the current approach, we anchor the width of the filter to the
% first value based on Rosso et al., (2021a,2023a), and then separately
% compute the other widths as logarithmic functions of the frequency bins 

% Settings
S = struct();     % initialize input structure to submit to cluster
S.targetfrex = 2.439; % stimulation frequency: hard-coded 1 / 0.410 s (range is 0.400 and 0.420 but 0.410 is the mode) 
S.fwhm  = .35;        % filter width
S.shr   = 0.01;       % shrinking proportion for regularization; value between 0 and 1 (until 0.001 it works for restoring the full rank of the matrix)
S.srate = 250;        % recorded at 1000 Hz, downsampled by 4
S.ncomps2keep = 100;  % how many components to keep (given the huge amount of components = N voxels)
S.comps2see   = [1 2 3]; % list components of interest, of which you want to save nifti image
S.path_script    = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia'; % path of the present script
S.path_source    = ['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1' cond_suffix{which_cond} '/'];
S.path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Mag_SourceLeakCorr/';
mkdir(S.path_GEDoutput)
S.name_folder = ['GEDoutput_frex_' num2str(S.targetfrex) 'Hz' cond_suffix{which_cond} '/']; 
% The following defines whether and how we are going to shuffle the sourde data 
% in permutation testing
% 0 = no permutation applied 
% 1 = randomize sources idxs (rows)
% 2 = randomize source idx indebendently per timepoint (rows, for every column)
% 3 = randomize timepoints (columns)
% 4 = randomize all together
S.permu_flag = 0;
S.permu_name = {'randLabel', 'randLabelPointWise', 'randTime', 'randAll'}; % to append to files/folder names for saving the data
% Double - check your permutation strategy
if S.permu_flag > 0
	warning('You have selected a permutation strategy. Double-check that you have selected the correct one, as submitting the job to the cluster might result in a very long computation time. Press SPACE BAR to continue')
    pause()
end

% Set all frequencies to test, and associated filter widths
if allfrex_flag == 1 

    % Set with respect to stimulation frequency
    nfrex_above = 80; % tested for even numbers
    nfrex_below = 6;   % tested for even numbers
    % Compute all frequencies
    frex_above = linspace(S.targetfrex, S.targetfrex * (nfrex_above/2), nfrex_above);
    frex_below = linspace(S.targetfrex / 2, S.targetfrex / (2 * nfrex_below), nfrex_below);
    frex_all = [frex_below(end:-1:1), frex_above];
    % compute all filter widths
    fwhm_above = logspace(log10(S.fwhm), log10(S.fwhm * (nfrex_above)), nfrex_above);
    fwhm_below = logspace(log10(S.fwhm), log10(S.fwhm / (nfrex_below)), nfrex_below+1);
    fwhm_all   = [fwhm_below(end:-1:2) fwhm_above];

    % Visualize frewquencies and filters
    figure
    plot(frex_all, fwhm_all, 's-')
    xlabel('Center frequencies (Hz)')
    ylabel('Filter FWHM (Hz)')
    title('Filter width as a function of center frequency')


    % Run GED with parallel computing for a set of frequencies
    nfrex = length(frex_all);
    for frexi = 1:nfrex
        S.targetfrex  = frex_all(frexi); % overwrite fields of input structure
        S.fwhm        = fwhm_all(frexi);
        S.name_folder = ['GEDoutput_frex_' num2str(S.targetfrex) 'Hz' cond_suffix{which_cond} '/'];
        % Running with parallel computing
        jobid = job2cluster(@GED_SingleSubjects_v3,S); % always give a structure as input
        
        disp(['Submitting GED for frequency = ' num2str(S.targetfrex) ' Hz and width = ' num2str(S.fwhm) ' Hz on cluster'])
        
    end

else
    
    % Run GED with parallel computing for single stimulation frequency
    jobid = job2cluster(@GED_SingleSubjects_v3,S); % always give a structure as input

end


%%



%% ASSESSMENT ONLY - Replicate original analysis with shorter durations

%% GED assessment (Rest and Listen)

% This blocks enables the visual assessment of GED, for individual
% participants and grand-averages, and produces NIFTI images
% - Power spectrum
% - Eigenspectrum
% - Activations maps (cleaned and saved as NIFTI for visualization in the next block)

clear all
close all

% Select duration and condition
duration_suffix = '_240_sec';



%GED list
list_up = dir(['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Duration/GED*' duration_suffix]);

% Frequency to inspect
% frex2see = 2.439;   %  stimfrex;
freq_vect_idx = 1:length(list_up);
flag_figures = 0;

% Select 'UniversitA della strada' approach
streetuni_flag = 1; % when set to 1, processes activation patterns for interpretation

%load('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/time.mat');
% Set home path
% path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/';
% path_suffix = {'Hz_rest', 'Hz_beat'};
% path_frex = [path_GEDoutput, 'GEDoutput_frex_', num2str(frex2see), path_suffix{which_cond}];


for frexi = 1:length(freq_vect_idx)
    
    path_frex = [list_up(freq_vect_idx(frexi)).folder '/' list_up(freq_vect_idx(frexi)).name];
    path_GEDoutput = path_frex;
    list = dir([path_frex '/SUBJ*.mat']);
    nsubs = length(list);
    ncomps = 30;         % components to carry from previous block, to produce scree plots
    comps2see = [1 2 3]; % expected top components, to assess in more detail
    srate = 250;
    
    for subi = 1:nsubs % over subjects
        
        load([list(subi).folder '/' list(subi).name],'evals','GEDts','GEDmap')
        npnts  = size(GEDts,2);       % number of time points in the signal
        nevals = length(evals);
        
        % Initialize matrix for all participants, once upon first iteration
        if subi == 1
            GEDpow_all = zeros(ncomps,npnts,nsubs);
            evals_all  = zeros(nevals,nsubs);
            GEDts_all  = zeros(ncomps,npnts,nsubs);
            GEDmap_all = zeros(ncomps,nevals,nsubs);
            
            % Settings for power spectrum
            tmax = npnts/srate; % duration, in seconds
            frexres = 1/tmax;   % Rayleigh frequency 1 / T(in seconds)
            % FFT parameters
            nfft = ceil( srate/frexres );  % length sufficiently high to guarantee resolution at denominator (= pnts)
            hz   = linspace(0,srate,nfft); % vector of frequencies
        end
        evals_all(:,subi)   = (evals.*100)./sum(evals); % assign and normalize, for plotting
        GEDts_all(:,:,subi) = GEDts(1:ncomps,:);        % assign for storing
        
        % Compute power spectrum
        for compi = comps2see
            % Demeaning
            GEDts_all(compi,:,subi) = GEDts_all(compi,:,subi) - mean(GEDts_all(compi,:,subi),2); % if not performed, some subjects showed prominent 0-Hz frequency component
            % Power spectrum
            GEDpow_all(compi,:,subi) = abs(fft(GEDts_all(compi,:,subi)',nfft,1)/npnts).^2;
        end
        
        % Assign activation patterns
        for compi = comps2see
            GEDmap_all(compi,:,subi) = GEDmap(compi,:);
        end
        
        disp(['GED assessment of freq idx #' num2str(frexi) ' - participant #' num2str(subi)])
    end
    
    % Compute grand-averages over subjects
    GEDts_avg =  squeeze( mean(GEDts_all,3) );
    GEDpow_avg = squeeze( mean(GEDpow_all,3) );
    evals_avg  = squeeze( mean(evals_all,2) );
    % Process activation pattern, if flag is enabled
    if streetuni_flag ==1
        GEDmap_avg = squeeze( mean(abs(GEDmap_all),3) );
    else
        GEDmap_avg = squeeze( mean(GEDmap_all,3) );
    end
    GEDmap_avg = GEDmap_avg' ./ (max(GEDmap_avg'));  % normalize to 1
    GEDmap_avg = GEDmap_avg';
    
    
    % Calculate distance of every participant's activation pattern from average, per component
    % (check the assumption that network A = COMP#1 and netowrk B = COMP#2)
    z_thresh = 2.3;  %~.01  % 1.96;
    idx_toofar = zeros(length(comps2see),nsubs);
    %Remove 'bad' matrices based on distance
    [GEDmap_dist,GEDmap_z] = deal ( zeros(ncomps,nsubs) );
    for compi = comps2see
        for subi = 1:nsubs
            %Compute distance
            GEDmap_dist(compi,subi) = abs( sqrt(trace(GEDmap_all(compi,:,subi)'*GEDmap_avg(compi,:))) );  % i take abs becasue I get complex numbers with one component = 0; it does not change the result
            %Normalize distance (temporary variable)
            GEDmap_z(compi,subi) = ( GEDmap_dist(compi,subi) - mean(GEDmap_dist(compi,:)) ) / std(GEDmap_dist(compi,:));
            
            
            
            %         %remove distances just for visualization, using logical indexing
            %         GEDmap_z(idx_toofar) = NaN;
        end
        % find indexes of outliers
        idx_toofar(compi,:) = abs(GEDmap_z(compi,:)) > z_thresh;
        disp( [ num2str(sum(idx_toofar(compi,:))) ' outliers identified for component #' num2str(compi) ] )
    end
    
    
    outliers_flag = 0;  % when set to 1, excludes activation patterns to deviate from average within component
    
    % plotting weights in the brain (nifti images)
    maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); % getting the mask for creating the figure
    for compi = comps2see % over components of interest
        
        % Start with assignment to temporary variable, to avoid overwriting later on
        if outliers_flag == 1
            GEDmap2see = squeeze( mean(abs(GEDmap_all(compi,:,~idx_toofar(compi,:))),3) ); % exclude outliers, if enabled
            GEDmap2see = GEDmap2see'./(max(GEDmap2see')); % and normalize
            GEDmap2see = GEDmap2see';
        else
            GEDmap2see = GEDmap_avg(compi, :);
        end
        
        %thresholding activation patterns according to mean + 1 standard deviation 
        thresh = mean(GEDmap2see) + std(GEDmap2see);
        GEDmap2see(GEDmap2see<thresh) = 0;
        
        % using nifti toolbox
        SO = GEDmap2see';
        % building nifti image
        SS = size(maskk.img);
        dumimg = zeros(SS(1),SS(2),SS(3),1); % no time here so size of 4D = 1
        for sourci = 1:size(SO,1) % over brain sources
            dumm = find(maskk.img == sourci); % finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
            [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); % getting subscript in 3D from index
            dumimg(i1,i2,i3,:) = SO(sourci,:); % storing values for all time-points in the image matrix
        end
        nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
        nii.img = dumimg; %storing matrix within image structure
        nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
        
        % Save output
        disp(['saving nifti image - component ' num2str(compi) ' frequency ' num2str(list_up(freq_vect_idx(frexi)).name(16:end))])
        %     if adj_compentr == 1
        %         save_nii(nii,['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1/Test_GED/CompEntrAdj_SubjsAverage_Comp_' num2str(frexi) '_Var_' num2str(round(evals_avg(compi))) '.nii']); %printing image
        %     else
        if outliers_flag ==1
            save_nii(nii,[path_GEDoutput '/SubjsAverage_Comp_' num2str(compi) '_Var_' num2str(round(evals_avg(compi))) 'frex_' num2str(list_up(freq_vect_idx(frexi)).name(16:end)) '_NO_outliers.nii']); % printing image
        else
            save_nii(nii,[path_GEDoutput '/SubjsAverage_Comp_' num2str(compi) '_Var_' num2str(round(evals_avg(compi))) 'frex_' num2str(list_up(freq_vect_idx(frexi)).name(16:end)) '.nii']); % printing image
        end
    end
end




%Visualize distances
if flag_figures == 1
    figure(1000), clf
    for compi = comps2see

        subplot(3,1,compi), hold on
        plot(GEDmap_z(compi,:),'bs-','linew',2,'markerfacecolor','w')
        legend('GED maps z-scores')
        plot(xlim, [z_thresh, z_thresh],'--');
        plot(xlim, -1*[z_thresh, z_thresh],'--');
        title('Distances before removal')
        %     %visualize difference after removal
        %     subplot(212)
        %     plot(GEDmap_z,'bs-','linew',2,'markerfacecolor','w')
        %     legend('GED maps z-scores' , 'CovR z-scores')
        %     yline(z_thresh,'--');
        %     yline(-z_thresh,'--');
        %     title('Distances after removal')


    end
end

% %Remove actual outlier matrices, using logical indexing
% GEDmap_all(compi,:,idx_toofar) = NaN;
% %Re-compute grand-averages (excluding NaNs)
% GEDmap_avg = squeeze(mean(covS, 3 , 'omitnan')); 

% Visualization
if flag_figures == 1

    maxfrex2see = 30;
    targetfrex = 2.439;
    [peakfrex, idx_peak] = deal( nan(length(comps2see),nsubs) );
    idx_stim = dsearchn(hz',targetfrex);
    % Individual participants
    for subi = 1:nsubs
        figure(subi),clf
        for compi = comps2see
            % Store peaks and their indexes, for check and eventual further analyses
            [peakfrex(compi,subi), idx_peak(compi,subi)] = max(smooth(GEDpow_all(compi,:,subi),1));
            % Spectrum
            subplot(3,2,(compi-1)*2+1), hold on
            plot(hz,smooth(GEDpow_all(compi,:,subi),1),'k','linew',1.7)
            set(gca,'xlim',[0 maxfrex2see])
            plot(hz(idx_peak(compi,subi)), peakfrex(compi,subi),'*','markerfacecolor','r','linew',5) % empirical peak
            plot(hz(idx_stim), peakfrex(compi,subi),'*','markerfacecolor','b','linew',5)         % expected stimulation frequency
            xlabel('Frequency (Hz)'), ylabel('Power (muV^2)') %or 'Signal-to-noise ratio (%)'
            title(['Power spectrum - C #' num2str(compi)])
            axis square
            % Eigenspectrum
            subplot(3,2,(compi-1)*2+2);     hold on
            plot(evals_all(:,subi),'s-','markerfacecolor','k','linew',2)       % eigenspectrum
            plot(compi,evals_all(compi,subi),'*','markerfacecolor','r','linew',5) % current component
            xlim([1 ncomps])
            xlabel('Component'), ylabel('Eigenvalue')
            title(['Eigenspectrum - C #' num2str(compi)])
            axis square
        end
    end
    % Grand-average
    figure(100),clf
    [peakfrex_avg, idx_peak_avg] = deal( nan(length(comps2see)) );
    for compi = comps2see
        
        [peakfrex_avg(compi), idx_peak_avg(compi)] = max(smooth(GEDpow_avg(compi,:),1));
        % Spectrum
        subplot(3,2,(compi-1)*2+1), hold on
        plot(hz,smooth(GEDpow_avg(compi,:),1),'k','linew',1.7)
        set(gca,'xlim',[0 maxfrex2see])
        plot(hz(idx_peak_avg(compi)), peakfrex_avg(compi),'*','markerfacecolor','r','linew',5) % empirical peak
        plot(hz(idx_stim), peakfrex_avg(compi),'*','markerfacecolor','b','linew',5)         % expected stimulation frequency
        xlabel('Frequency (Hz)'), ylabel('Power (muV^2)') %or 'Signal-to-noise ratio (%)'
        title(['Power spectrum - C #' num2str(compi)])
        axis square
        % Eigenspectrum
        subplot(3,2,(compi-1)*2+2);     hold on
        plot(evals_avg,'s-','markerfacecolor','k','linew',2)       % eigenspectrum
        plot(compi,evals_avg(compi),'*','markerfacecolor','r','linew',5) % current component
        xlim([1 ncomps])
        xlabel('Component'), ylabel('Eigenvalue')
        title(['Eigenspectrum - C #' num2str(compi)])
        axis square
        
    end
end



%% Assess frequency-specificity of GED

clear all
close all
clc

duration_suffix = '_240_sec';

% Select components of interest
comps2see = [1:10];
stimfrex = 2.439;

% Set home path
path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Duration/'; % '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/';
path_suffix = {'Hz_rest', 'Hz_beat'};


% Load beat listening data
list_beat = dir(['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Duration/GEDoutput_frex_' num2str(stimfrex) 'Hz_beat_30_sec/SUBJ*.mat']); % dir(['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/GEDoutput_frex_' num2str(stimfrex) 'Hz_beat/SUBJ*.mat']); % do it for stimfrex, but it's the same for any frequency
% Get indexes of actual subjects (to match the other list)
subs_beat = zeros(length(list_beat),1);
for subi = 1:length(list_beat)
    % Retrieve subject
    this_sub = regexp(list_beat(subi).name, '([\d.]+).mat', 'tokens');
    % Assign frequency to vector
    subs_beat(subi) = str2double(this_sub{1,1});
end

% Load resting state data
list_rest = dir(['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Duration/GEDoutput_frex_' num2str(stimfrex) 'Hz_rest_30_sec/SUBJ*.mat']); % dir(['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/GEDoutput_frex_' num2str(stimfrex) 'Hz_rest/SUBJ*.mat']); % do it for stimfrex, but it's the same for any frequency
subs_rest = zeros(length(list_rest),1);
% Get indexes of actual subjects (to match the other list)
for subi = 1:length(list_rest)
    % Retrieve subject
    this_sub = regexp(list_rest(subi).name, '([\d.]+).mat', 'tokens');
    % Assign frequency to vector
    subs_rest(subi) = str2double(this_sub{1,1});
end    

% Align participants for pair-wise comparisons
[label_pairs,idx_beat,idx_rest] = intersect(subs_beat',subs_rest'); 
% Double-check that all participants are paired across conditions
if sum(subs_beat(idx_beat) - subs_rest(idx_rest)) == 0
    disp('All participants are matched across conditions!')
else
    error('WARNING: not all participants are matched across conditions!')
end

% Indices for subjects within condition to be used in the loop
idx4loop{1} = idx_rest; idx4loop{2} = idx_beat;
save('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Subjs_idx_match.mat','idx4loop')

% Loop over Rest and Listen
for condi = 1:2

    % Get all_frex from home path
    list_GED = dir(path_GEDoutput);
    list_GED = list_GED( [list_GED.isdir] );
    list_GED = list_GED( contains({list_GED.name}, path_suffix{condi}) & contains({list_GED.name}, duration_suffix));
    list_GED = list_GED(startsWith({list_GED.name}, 'GED')); % added since GED and PCA coexist in the same folder


    % Initialize vectors for storing frequencies
    nfrex = length(list_GED);
    all_frex = zeros(nfrex,1);
    % Initialize topEvals
    topEvals_avg = zeros(length(comps2see),nfrex);% avg of top components across subjects, per frequency

    % Import and process top eigenvalues, per frequency
    for frexi = 1:nfrex

        % Get folder name
        this_folder = [path_GEDoutput, list_GED(frexi).name '/'];   
        % Retrieve frequency
        this_frex = regexp(this_folder, ['([\d.]+)' path_suffix{condi}], 'tokens');
        % Assign frequency to vector
        all_frex(frexi) = str2double(this_frex{1,1});

        % All subjects within frequency frexi
        fulllist = dir([this_folder '/SUBJ*.mat']);
        
        % Sublist of only subjects which are matched between the two conditions
        list_subs = fulllist(idx4loop{condi});
        nsubs = length(list_subs);
        % Process evals from subjects
        for subi = 1:nsubs % over subjects

            % Get eigenvalues (already sorted and normalized)
            load([this_folder, list_subs(subi).name],'evals')
            nevals = length(evals);
            
            if subi == 1 && frexi == 1 && condi == 1
                % Initialize matrix for all participants, once upon first iteration
                evals_all = zeros(nevals,nsubs,nfrex,2); % for Rest and Listen conditions
            end
            
            evals_all(:,subi,frexi,condi) = evals; % assign for computing average and for exporting
            

        end

        % Compute grand-averages over subjects (just for visualization)   
        topEvals_avg(:,frexi) = squeeze( mean(evals_all(comps2see,:,frexi,condi),2) );

        disp(['Processing evals of frex ' this_frex{1,1} 'Hz'])

    end
    % Re-sort vector of frequencies (strings follow alphabetic order)
    [all_frex, idx_sort] = sort(all_frex, 'ascend');
    topEvals_avg = topEvals_avg(:,idx_sort);

    % Visualize Rest and Listen
    figure(300 + condi), clf
    subplot(311), hold on
    plot(all_frex,topEvals_avg,'s-','markerfacecolor','k','linew',2)
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    xticks(round(all_frex,1)), xtickangle(45)
    %xlim([0 nfrex+1])
    xlabel('Frequencies (Hz)')
    ylabel('Explained variance')
    legendEntries = cell(length(comps2see),1);
    for legi = 1:length(legendEntries)
       legendEntries{legi} = ['Component #' num2str(comps2see(legi))];
    end
    legend(legendEntries, 'Location', 'Best');
    title(['Top ' num2str(length(comps2see)) ' components'])
    subplot(312), hold on
    plot(all_frex,sum(topEvals_avg,1),'s-','color','k','markerfacecolor','k','linew',2)
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    xticks(all_frex), xtickangle(45)
    %xlim([0 nfrex+1])
    xlabel('Frequencies (Hz)')
    ylabel('Explained variance')
    title(['Sum of top ' num2str(length(comps2see)) ' components'])
    subplot(313), hold on
    plot(all_frex,std(topEvals_avg,0,1),'s-','color','r','markerfacecolor','k','linew',2)
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    xticks(all_frex), xtickangle(45)
    %xlim([0 nfrex+1])
    % ylim([0 50])
    xlabel('Frequencies (Hz)')
    ylabel('Explained variance')
    title(['STD of top ' num2str(length(comps2see)) ' components'])

end





%% ASSESSMENT ONLY - GRADIOMETERES

%% GED assessment (Rest and Listen)

% This blocks enables the visual assessment of GED, for individual
% participants and grand-averages, and produces NIFTI images
% - Power spectrum
% - Eigenspectrum
% - Activations maps (cleaned and saved as NIFTI for visualization in the next block)

clear all
close all


%GED list
list_up = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Gradiometers/GED*');

% Frequency to inspect
% frex2see = 2.439;   %  stimfrex;
freq_vect_idx = 1:length(list_up);
flag_figures = 0;

% Select 'UniversitA della strada' approach
streetuni_flag = 1; % when set to 1, processes activation patterns for interpretation

%load('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/time.mat');
% Set home path
% path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/';
% path_suffix = {'Hz_rest', 'Hz_beat'};
% path_frex = [path_GEDoutput, 'GEDoutput_frex_', num2str(frex2see), path_suffix{which_cond}];


for frexi = 1:length(freq_vect_idx)
    
    path_frex = [list_up(freq_vect_idx(frexi)).folder '/' list_up(freq_vect_idx(frexi)).name];
    path_GEDoutput = path_frex;
    list = dir([path_frex '/SUBJ*.mat']);
    nsubs = length(list);
    ncomps = 30;         % components to carry from previous block, to produce scree plots
    comps2see = [1 2 3]; % expected top components, to assess in more detail
    srate = 250;
    
    for subi = 1:nsubs % over subjects
        
        load([list(subi).folder '/' list(subi).name],'evals','GEDts','GEDmap')
        npnts  = size(GEDts,2);       % number of time points in the signal
        nevals = length(evals);
        
        % Initialize matrix for all participants, once upon first iteration
        if subi == 1
            GEDpow_all = zeros(ncomps,npnts,nsubs);
            evals_all  = zeros(nevals,nsubs);
            GEDts_all  = zeros(ncomps,npnts,nsubs);
            GEDmap_all = zeros(ncomps,nevals,nsubs);
            
            % Settings for power spectrum
            tmax = npnts/srate; % duration, in seconds
            frexres = 1/tmax;   % Rayleigh frequency 1 / T(in seconds)
            % FFT parameters
            nfft = ceil( srate/frexres );  % length sufficiently high to guarantee resolution at denominator (= pnts)
            hz   = linspace(0,srate,nfft); % vector of frequencies
        end
        evals_all(:,subi)   = (evals.*100)./sum(evals); % assign and normalize, for plotting
        GEDts_all(:,:,subi) = GEDts(1:ncomps,:);        % assign for storing
        
        % Compute power spectrum
        for compi = comps2see
            % Demeaning
            GEDts_all(compi,:,subi) = GEDts_all(compi,:,subi) - mean(GEDts_all(compi,:,subi),2); % if not performed, some subjects showed prominent 0-Hz frequency component
            % Power spectrum
            GEDpow_all(compi,:,subi) = abs(fft(GEDts_all(compi,:,subi)',nfft,1)/npnts).^2;
        end
        
        % Assign activation patterns
        for compi = comps2see
            GEDmap_all(compi,:,subi) = GEDmap(compi,:);
        end
        
        disp(['GED assessment of freq idx #' num2str(frexi) ' - participant #' num2str(subi)])
    end
    
    % Compute grand-averages over subjects
    GEDts_avg =  squeeze( mean(GEDts_all,3) );
    GEDpow_avg = squeeze( mean(GEDpow_all,3) );
    evals_avg  = squeeze( mean(evals_all,2) );
    % Process activation pattern, if flag is enabled
    if streetuni_flag ==1
        GEDmap_avg = squeeze( mean(abs(GEDmap_all),3) );
    else
        GEDmap_avg = squeeze( mean(GEDmap_all,3) );
    end
    GEDmap_avg = GEDmap_avg' ./ (max(GEDmap_avg'));  % normalize to 1
    GEDmap_avg = GEDmap_avg';
    
    
    % Calculate distance of every participant's activation pattern from average, per component
    % (check the assumption that network A = COMP#1 and netowrk B = COMP#2)
    z_thresh = 2.3;  %~.01  % 1.96;
    idx_toofar = zeros(length(comps2see),nsubs);
    %Remove 'bad' matrices based on distance
    [GEDmap_dist,GEDmap_z] = deal ( zeros(ncomps,nsubs) );
    for compi = comps2see
        for subi = 1:nsubs
            %Compute distance
            GEDmap_dist(compi,subi) = abs( sqrt(trace(GEDmap_all(compi,:,subi)'*GEDmap_avg(compi,:))) );  % i take abs becasue I get complex numbers with one component = 0; it does not change the result
            %Normalize distance (temporary variable)
            GEDmap_z(compi,subi) = ( GEDmap_dist(compi,subi) - mean(GEDmap_dist(compi,:)) ) / std(GEDmap_dist(compi,:));
            
            
            
            %         %remove distances just for visualization, using logical indexing
            %         GEDmap_z(idx_toofar) = NaN;
        end
        % find indexes of outliers
        idx_toofar(compi,:) = abs(GEDmap_z(compi,:)) > z_thresh;
        disp( [ num2str(sum(idx_toofar(compi,:))) ' outliers identified for component #' num2str(compi) ] )
    end
    
    
    outliers_flag = 0;  % when set to 1, excludes activation patterns to deviate from average within component
    
    % plotting weights in the brain (nifti images)
    maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); % getting the mask for creating the figure
    for compi = comps2see % over components of interest
        
        % Start with assignment to temporary variable, to avoid overwriting later on
        if outliers_flag == 1
            GEDmap2see = squeeze( mean(abs(GEDmap_all(compi,:,~idx_toofar(compi,:))),3) ); % exclude outliers, if enabled
            GEDmap2see = GEDmap2see'./(max(GEDmap2see')); % and normalize
            GEDmap2see = GEDmap2see';
        else
            GEDmap2see = GEDmap_avg(compi, :);
        end
        
        %thresholding activation patterns according to mean + 1 standard deviation 
        thresh = mean(GEDmap2see) + std(GEDmap2see);
        GEDmap2see(GEDmap2see<thresh) = 0;
        
        % using nifti toolbox
        SO = GEDmap2see';
        % building nifti image
        SS = size(maskk.img);
        dumimg = zeros(SS(1),SS(2),SS(3),1); % no time here so size of 4D = 1
        for sourci = 1:size(SO,1) % over brain sources
            dumm = find(maskk.img == sourci); % finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
            [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); % getting subscript in 3D from index
            dumimg(i1,i2,i3,:) = SO(sourci,:); % storing values for all time-points in the image matrix
        end
        nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
        nii.img = dumimg; %storing matrix within image structure
        nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
        
        % Save output
        disp(['saving nifti image - component ' num2str(compi) ' frequency ' num2str(list_up(freq_vect_idx(frexi)).name(16:end))])
        %     if adj_compentr == 1
        %         save_nii(nii,['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1/Test_GED/CompEntrAdj_SubjsAverage_Comp_' num2str(frexi) '_Var_' num2str(round(evals_avg(compi))) '.nii']); %printing image
        %     else
        if outliers_flag ==1
            save_nii(nii,[path_GEDoutput '/SubjsAverage_Comp_' num2str(compi) '_Var_' num2str(round(evals_avg(compi))) 'frex_' num2str(list_up(freq_vect_idx(frexi)).name(16:end)) '_NO_outliers.nii']); % printing image
        else
            save_nii(nii,[path_GEDoutput '/SubjsAverage_Comp_' num2str(compi) '_Var_' num2str(round(evals_avg(compi))) 'frex_' num2str(list_up(freq_vect_idx(frexi)).name(16:end)) '.nii']); % printing image
        end
    end
end

%% Assess frequency-specificity of GED

clear all
close all
clc

% Select components of interest
comps2see = [1:10];
stimfrex = 2.439;

% Set home path
path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Gradiometers/';
path_suffix = {'Hz_rest', 'Hz_beat'};


% Select the index of 
%which_lfo  = 1; % the component to use as low-frequency modulator
which_comp = 1; % the index of the component for the high-frequency carrier

% Load beat listening data
list_beat = dir(['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Gradiometers/GEDoutput_frex_' num2str(stimfrex) 'Hz_beat/SUBJ*.mat']); % do it for stimfrex, but it's the same for any frequency
% Get indexes of actual subjects (to match the other list)
subs_beat = zeros(length(list_beat),1);
for subi = 1:length(list_beat)
    % Retrieve subject
    this_sub = regexp(list_beat(subi).name, '([\d.]+).mat', 'tokens');
    % Assign frequency to vector
    subs_beat(subi) = str2double(this_sub{1,1});
end

% Load resting state data
list_rest = dir(['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Gradiometers/GEDoutput_frex_' num2str(stimfrex) 'Hz_rest/SUBJ*.mat']); % do it for stimfrex, but it's the same for any frequency
subs_rest = zeros(length(list_rest),1);
% Get indexes of actual subjects (to match the other list)
for subi = 1:length(list_rest)
    % Retrieve subject
    this_sub = regexp(list_rest(subi).name, '([\d.]+).mat', 'tokens');
    % Assign frequency to vector
    subs_rest(subi) = str2double(this_sub{1,1});
end    

% Align participants for pair-wise comparisons
[label_pairs,idx_beat,idx_rest] = intersect(subs_beat',subs_rest'); 
% Double-check that all participants are paired across conditions
if sum(subs_beat(idx_beat) - subs_rest(idx_rest)) == 0
    disp('All participants are matched across conditions!')
else
    error('WARNING: not all participants are matched across conditions!')
end

% Indices for subjects within condition to be used in the loop
idx4loop{1} = idx_rest; idx4loop{2} = idx_beat;
save('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Subjs_idx_match.mat','idx4loop')

% Loop over Rest and Listen
for condi = 1:2

    % Get all_frex from home path
    list_GED = dir(path_GEDoutput);
    list_GED = list_GED( [list_GED.isdir] );
    list_GED = list_GED(endsWith({list_GED.name}, path_suffix{condi}));
    list_GED = list_GED(startsWith({list_GED.name}, 'GED')); % added since GED and PCA coexist in the same folder


    % Initialize vectors for storing frequencies
    nfrex = length(list_GED);
    all_frex = zeros(nfrex,1);
    % Initialize topEvals
    topEvals_avg = zeros(length(comps2see),nfrex);% avg of top components across subjects, per frequency

    % Import and process top eigenvalues, per frequency
    for frexi = 1:nfrex

        % Get folder name
        this_folder = [path_GEDoutput, list_GED(frexi).name '/'];   
        % Retrieve frequency
        this_frex = regexp(this_folder, ['([\d.]+)' path_suffix{condi}], 'tokens');
        % Assign frequency to vector
        all_frex(frexi) = str2double(this_frex{1,1});

        % All subjects within frequency frexi
        fulllist = dir([this_folder '/SUBJ*.mat']);
        
        % Sublist of only subjects which are matched between the two conditions
        list_subs = fulllist(idx4loop{condi}); %:
        nsubs = length(list_subs);
        % Process evals from subjects
        for subi = 1:nsubs % over subjects

            % Get eigenvalues (already sorted and normalized)
            load([this_folder, list_subs(subi).name],'evals')
            nevals = length(evals);
            
            if subi == 1 && frexi == 1 && condi == 1
                % Initialize matrix for all participants, once upon first iteration
                evals_all = zeros(nevals,nsubs,nfrex,2); % for Rest and Listen conditions
            end
            
            evals_all(:,subi,frexi,condi) = evals; % assign for computing average and for exporting
            
        end

        % Compute grand-averages over subjects (just for visualization)   
        topEvals_avg(:,frexi) = squeeze( mean(evals_all(comps2see,:,frexi,condi),2) );

        disp(['Processing evals of frex ' this_frex{1,1} 'Hz'])

    end
    % Re-sort vector of frequencies (strings follow alphabetic order)
    [all_frex, idx_sort] = sort(all_frex, 'ascend');
    topEvals_avg = topEvals_avg(:,idx_sort);

    % Visualize Rest and Listen
    figure(300 + condi), clf
    subplot(311), hold on
    plot(all_frex,topEvals_avg,'s-','markerfacecolor','k','linew',2)
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    xticks(round(all_frex,1)), xtickangle(45)
    %xlim([0 nfrex+1])
    xlabel('Frequencies (Hz)')
    ylabel('Explained variance')
    legendEntries = cell(length(comps2see),1);
    for legi = 1:length(legendEntries)
       legendEntries{legi} = ['Component #' num2str(comps2see(legi))];
    end
    legend(legendEntries, 'Location', 'Best');
    title(['Top ' num2str(length(comps2see)) ' components'])
    subplot(312), hold on
    plot(all_frex,sum(topEvals_avg,1),'s-','color','k','markerfacecolor','k','linew',2)
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    xticks(all_frex), xtickangle(45)
    %xlim([0 nfrex+1])
    xlabel('Frequencies (Hz)')
    ylabel('Explained variance')
    title(['Sum of top ' num2str(length(comps2see)) ' components'])
    subplot(313), hold on
    plot(all_frex,std(topEvals_avg,0,1),'s-','color','r','markerfacecolor','k','linew',2)
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    xticks(all_frex), xtickangle(45)
    %xlim([0 nfrex+1])
    % ylim([0 50])
    xlabel('Frequencies (Hz)')
    ylabel('Explained variance')
    title(['STD of top ' num2str(length(comps2see)) ' components'])

end







%%

%% STATISTICAL TESTING

%% Statistical testing of randomization strategies

clear all
close all

comps2test = [1:3]; % [1:20]
stimfrex = 2.439;
load('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Stats/TestRands.mat'); %loading GED components eigenvalues

% Sorting frequencies
evals_all = evals_all(:,:,idx_sort,:); % only sort these: the randomizations where correctly sorted already
%evals_all_rands = evals_all_rands(:,:,idx_sort,:);

% Assign to panel of Figure 4
figA = evals_all(:,:,:,1);       % Resting state, no randomization
figB = evals_all_rands(:,:,:,1); % Resting state, label randomization
figC = evals_all_rands(:,:,:,2); % Resting state, label-pointwise randomization

% Select components of interest
ncomps = length(comps2test);
nfrex = size(evals_all,3);
% Test eigenvalues
[PevalsAB,TevalsAB,PevalsAC,TevalsAC] = deal(zeros(ncomps,nfrex));
for compi = comps2test
    for frexi = 1:nfrex
        % AB comparison
        [p,~,stats] = signrank( squeeze(figA(compi,:,frexi)) , squeeze(figB(compi,:,frexi)));
        PevalsAB(compi,frexi) = p;
        TevalsAB(compi,frexi) = stats.zval;
        % AC comparison
        [p,~,stats] = signrank( squeeze(figA(compi,:,frexi)) , squeeze(figC(compi,:,frexi)));
        PevalsAC(compi,frexi) = p;
        TevalsAC(compi,frexi) = stats.zval;
    end
    
    % Visualize eigenvalues
    figure(1000 + compi), clf, title('AB comparison')
    subplot(311)
    plot(all_frex, squeeze( figA(compi,:,:) ) , 'k'), hold on
    plot(all_frex, squeeze( figB(compi,:,:) ) , 'r')
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    subplot(312)
    plot(all_frex, TevalsAB(compi,:)), hold on
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    subplot(313)
    plot(all_frex, PevalsAB(compi,:)), hold on
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    ylim([0 1])
     % Visualize eigenvalues
    figure(2000 + compi), clf, title('AC comparison')
    subplot(311)
    plot(all_frex, squeeze( figA(compi,:,:) ) , 'k'), hold on
    plot(all_frex, squeeze( figC(compi,:,:) ) , 'r')
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    subplot(312)
    plot(all_frex, TevalsAC(compi,:)), hold on
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    subplot(313)
    plot(all_frex, PevalsAC(compi,:)), hold on
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    ylim([0 1])
end


%FDR correction
% AB comparison
PfdrAB = zeros(size(PevalsAB,1),1); %FDR thresholds (independently for each GED component)
P_threshAB = zeros(size(PevalsAB)); %frequencies (variance explained) which are significant after FDR correction
for ii = 1:size(PevalsAB,1) %over GED components
    PfdrAB(ii,1) = fdr(PevalsAB(ii,:)); %FDR correction (getting the threshold)
    P_threshAB(ii,PevalsAB(ii,:)<=PfdrAB(ii,1)) = 1; %getting the significant frequencies (independently for each component)
end
%AC comparison
PfdrAC = zeros(size(PevalsAC,1),1); %FDR thresholds (independently for each GED component)
P_threshAC = zeros(size(PevalsAC)); %frequencies (variance explained) which are significant after FDR correction
for ii = 1:size(PevalsAC,1) %over GED components
    PfdrAC(ii,1) = fdr(PevalsAC(ii,:)); %FDR correction (getting the threshold)
    P_threshAC(ii,PevalsAC(ii,:)<=PfdrAC(ii,1)) = 1; %getting the significant frequencies (independently for each component)
end




%% Statistical testing of GED on FIGURES material (all supplementary info)

clear all
close all

% Set figure to analyze
figName = 'S6_PCA.mat';


% Set home path
path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Images_revisions/'; % '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/';

comps2test = [1:20];
stimfrex = 2.439;
load([path_GEDoutput, figName]);%loading GED components eigenvalues
evals_all = evals_all(:,:,idx_sort,:); %sorting frequencies
% Select components of interest
ncomps = length(comps2test);
nfrex = size(evals_all,3);
% Test eigenvalues
[Pevals,Tevals] = deal(zeros(ncomps,nfrex));
for compi = comps2test
    for frexi = 1:nfrex
%         [~,p,~,stats] = ttest( squeeze(evals_rest(compi,:,frexi)) , squeeze(evals_beat(compi,:,frexi)), 'Tail','right');
%         P(compi,frexi) = p;
%         T(compi,frexi) = stats.tstat;
        [p,~,stats] = signrank( squeeze(evals_all(compi,:,frexi,2)) , squeeze(evals_all(compi,:,frexi,1)));
        Pevals(compi,frexi) = p;
        Tevals(compi,frexi) = stats.zval;
    end
    
%     % Visualize eigenvalues
%     figure(1000 + compi), clf
%     subplot(311)
%     plot(all_frex, squeeze( evals_all(compi,:,:,1) ) , 'k'), hold on
%     plot(all_frex, squeeze( evals_all(compi,:,:,2) ) , 'r')
%     line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
%     subplot(312)
%     plot(all_frex, Tevals(compi,:)), hold on
%     line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
%     subplot(313)
%     plot(all_frex, Pevals(compi,:)), hold on
%     line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
%     ylim([0 1])
end


%FDR correction
PFDR = zeros(size(Pevals,1),1); %FDR thresholds (independently for each GED component)
P_thresh = zeros(size(Pevals)); %frequencies (variance explained) which are significant after FDR correction
for ii = 1:size(Pevals,1) %over GED components
    PFDR(ii,1) = fdr(Pevals(ii,:)); %FDR correction (getting the threshold)
    P_thresh(ii,Pevals(ii,:)<PFDR(ii,1)) = 1; %getting the significant frequencies (independently for each component)
end

% Display significant frequencies
for compi = 1:3
   disp(['Significant frequencies for component # ' num2str(compi)])
   disp( all_frex(find(P_thresh(compi,:))) );
    
end






%% SUPPLEMENTARY TABLES

% Generate Table S3, containing cumulative varaince for top components over
% frequencies

load('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Images/GEDallsubjs.mat')
path_output = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Revisions/Tables/';

% Settings
nconds = 2;
ncomps2see = [1:10];
condlabels = {'cumvariance_rest'; 'cumvariance_beat'};

% Re-sort frequencies
evals_all = evals_all(:,:,idx_sort,:);

% Compute average across subjects
avg_evals = squeeze(mean(evals_all,2));

% Just to inspect
for ni = ncomps2see
    for fi = idx2pick
        disp(['N = ' num2str(ni) ' ; F = ' num2str(all_frex(fi)) 'Hz'])
        
        sum( avg_evals(1:ni,fi,which_cond) );
        
        out_sum = avg_evals(1:ni,fi,which_cond);
        
    end
end

% Store output for all frequencies
out_sum = zeros(length(ncomps2see),length(all_frex),nconds);
for condi = 1:nconds
    this_output = [];
    for ni = ncomps2see
        for frexi = 1:length(all_frex)
            % compute
            out_sum(ni,frexi,condi) = sum( avg_evals(1:ni,frexi,condi) );

            
        end
    end
    % Append labels for frex and numbers
   
    this_output = [all_frex' ;  out_sum(:,:,condi)];
    this_output = [[nan;ncomps2see'], this_output];
    save([path_output condlabels{condi} '.mat'], 'this_output')
end

