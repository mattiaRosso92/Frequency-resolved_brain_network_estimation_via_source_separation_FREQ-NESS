
%% Pipeline for separation of frequency-specific brain networks on brain voxel data 

% This script performs the following macro steps on MEG data

% - Preprocessing
% - Source reconstruction with Beamforming
% - Generalized eigendecomposition on source data

% The experimental design includes a condition of resting state and
% a condition of passive listening to a 2.4Hz beat


% Mattia Rosso & Leonardo Bonetti
% Aarhus, 18/03/2024

%% SELECT CONDITION

% Choose whether you want to process the condition of Resting State (1)
% or Beat Listening (2)

clear
which_cond = 2;


%% Maxfilter

% To set-up maxfilter the first time, please follow the following steps (source:  wiki.pet.auh.dk/wiki/Utility_functions_for_CFIN_servers)
% 1. open terminal window and type: 
%       gedit .bashrc
% 2. copy the following line to the beginning of the file, save and close
%    the window: 
%       source /usr/local/common/meeg-cfin/configurations/setup_environment.sh
% 3. in same terminal window, type: 
%       gedit .bash_profile
% 4. copy following 3 lines into the file, save and close window:
%       if [ -f $HOME/.bashrc ]; then
%           source $HOME/.bashrc
%       fi

% Before running maxfilter you need to close matlab, open the terminal and write: 'use anaconda', then open matlab and run maxfilter script

maxfilter_path = '/neuro/bin/util/maxfilter'; % insert your path here
project = 'MINDLAB2018_MEG-LearningBach-MemoryInformation';
maxDir = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc'; %output path
movement_comp = 0; %1 = yes; 0 = no

% Depending on which_cond, we are going to select files containing the
% following chunks of strings in their name
chunk_cond = {'est'; 'eat'};

path = '/raw/sorted/MINDLAB2018_MEG-LearningBach-MemoryInformation'; % path with all the subjects' folders
foldi = dir([path '/0*']); % list all the folders starting with '0' in order to avoid hidden files
for subi = 2:length(foldi) % over subjects (NOTE: we skip SUBJ0001 becasue of wrong beat stimuli)
    cart = [foldi(subi).folder '/' foldi(subi).name]; % create a path which combines the folder name and the file name
    pnana = dir([cart '/2*']); % search for folders starting with '2'
    for pp = 1:length(pnana) % loop to explore ad analyze all the folders inside the path above
        cart2 = [pnana(1).folder '/' pnana(pp).name];
        pr = dir([cart2 '/ME*']); % look for MEG folder
        if ~isempty(pr) % if pr is not empty, proceed with subfolders inside the meg path
            pnunu = dir([pr(1).folder '/' pr(1).name '/0*']);
            for dd = 1:length(pnunu) % over experimental blocks
                if contains(pnunu(dd).name,chunk_cond{which_cond})
                    fpath = dir([pnunu(1).folder '/' pnunu(dd).name '/files/*.fif']); % looks for .fif file
                    rawName = ([fpath.folder '/' fpath.name]); % assigns the final path of the .fif file to the rawName path used in the maxfilter command
                    maxfName = ['SUBJ' foldi(subi).name '_' fpath.name(1:end-4)]; % define the output name of the maxfilter processing
                    if movement_comp == 1
                        % movement compensation
                        % NOTE: took -ds 4 out for now
                        cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    else %no movement compensation (to be used if HPI coils did not work properly)
                        cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssds.fif'] ' -st 4 -corr 0.98 -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    end
                    system(cmd);
                end
            end
        end
    end
end


%% LBPD_startup_D

% NOTE: Do not change any of the functions that you can open from the paths below

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; % path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') % add the path where the function is, to submit the jobs to the server


%% Converting the .fif files into SPM objects

% NOTE: This section never runs in parallel; it puts the output in the current directory
% Make sure you are in the right path before running this section 

% SPM objects consist of a .dat file, containing the actual data, and a .mat
% file, being their mask

% Path for dir() to either Resting state {1} or Beat listening {2} condition
dir_cond = {'/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/*rest10_tss*.fif';... 
            '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/*eat*.fif'};

% Dir to .fif files, dir to the pre-selected condition
fif_list = dir(dir_cond{which_cond}); 

% Skip the first participant for this dataset (wrong audio stimulation)
for subi = 2:length(fif_list) % loop over the .fif files for each subject
    S = []; % initialize structure 'S'                   
    S.dataset = [fif_list(subi).folder '/' fif_list(subi).name];
    D = spm_eeg_convert(S); 
%     D = job2cluster(@cluster_spmobject, S); %actual function for conversion
end

% NOTE on D structures: after typing D. in the command window, you can pres TAB to
% visalize its fields in a tiny window hovering on the command line


%% Removing bad segments using OSLVIEW

% Check data for potential bad segments in the recordings (periods).
% Marking is done by right-clicking in the proximity of the event and click
% on 'mark event'.
% A first click (green dashed label) marks the beginning of a bad period
% A second click indicates the end of a bad period (red)
% Push the disk button to save to disk (no prefix will be added, same name is kept)

% NOTE: remember to check for bad segments of the signal both at 'megplanar' and 'megmag' channels (you can change the channels in the OSLVIEW interface)

% NOTE: remember to mark the trial within the bad segments as 'badtrials' and use the label for removing them from the Averaging (after Epoching) 

% Dir to either Resting state {1} or Beat listening {2} condition
dir_cond = {'/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/spmeeg*rest10*.mat';... 
            '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/spmeeg*eat*.mat'};

% Dir to SPM objects, dir to the pre-selected condition
spm_list = dir(dir_cond{which_cond}); 

for subi = 1:length(spm_list) 
    D = spm_eeg_load([spm_list(subi).folder '/' spm_list(subi).name]);
    D = oslview(D);
    D.save(); % save the selected bad segments and/or channels in OSLVIEW
    disp(subi)
end


%% AFRICA denoising (part I)

% If this is the first time running AFRICA denoising, type 'select_repository_version' 
% in your command window and select the latest option

% Setting up the cluster
clusterconfig('scheduler', 'cluster'); % 'none' for local testing
clusterconfig('long_running', 1); %for now keep 1; %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path where is the function for submit the jobs to the server


%% ICA calculation on the cluster

% Dir to either Resting state {1} or Beat listening {2} condition
dir_cond = {'/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/spmeeg*rest10*.mat';... 
            '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/spmeeg*eat*.mat'};

% Dir to SPM objects, dir to the pre-selected condition
spm_list = dir(dir_cond{which_cond}); 

for subi = 1:length(spm_list)
    S = [];
    D = spm_eeg_load([spm_list(subi).folder '/' spm_list(subi).name]);
    S.D = D;
    % Submit the job to the cluster and assign id
    jobid = job2cluster(@cluster_africa,S);
%   D = osl_africa(D,'do_ica',true,'do_ident',false,'do_remove',false,'used_maxfilter',true); 
%   D.save();
end

%% AFRICA denoising (part II)

% Visual inspection and removal of artifacted components
% Look into EOG and ECG channels: components are sorted by correlation with 
% these sensors in these visualization modes

% Dir to either Resting state {1} or Beat listening {2} condition
dir_cond = {'/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/spmeeg*rest10*.mat';... 
            '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/spmeeg*eat*.mat'};

% Dir to SPM objects, dir to the pre-selected condition
spm_list = dir(dir_cond{which_cond}); 

for subi = 1:length(spm_list) 
    D = spm_eeg_load([spm_list(subi).folder '/' spm_list(subi).name]);
    D = osl_africa(D,'do_ident','manual','do_remove',false,'artefact_channels',{'EOG','ECG'}); % manual identification; don't remove identified componentss yet; sort by EOG and then by ECG
    % hacking the function to manage to get around the OUT OF MEMORY problem
    S = [];
    S.D = D;
    jobid = job2cluster(@cluster_rembadcomp,S); %here we send the identified components to the cluster for removal
%   D.save();
    disp(subi)
end

%% OPTIONAL: Visual check for events timing
% 
% for subi = 1:length(spm_list)
% 
% D = spm_eeg_load([spm_list(subi).folder '/' spm_list(subi).name]);
% events2see = D.events;
% figure(100),clf
% plot(diff([events2see(:).time]))
% title(['Events timing particiant #' num2str(subi)])
% %plot(diff([event(1:200).time]))
% pause;
% 
% end

%% Epoching

% Produce one single epoch, taking the shortest trial duration across
% subjects and conditions as the standard epoch length

% Adds this prefix to epoched files
prefix2add = 'e'; 

% Dir to either Resting state {1} or Beat listening {2} condition
dir_cond = {'/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/spmeeg*rest10*.mat';... 
            '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/spmeeg*eat*.mat'};
        
% Dir to SPM objects, dir to the pre-selected condition
spm_list = dir(dir_cond{which_cond}); 

% more condition-specific settings
condlabels = {'Resting_state','Beat_listening'};
eventvalue = [1,160];
idx_wrongrest = [5 8]; % indexes of participants with wrong first trigger, in RESTING STATE condition; exclude from the loop
list4loop = {[1:idx_wrongrest(1)-1, idx_wrongrest(1)+1:idx_wrongrest(2)-1  idx_wrongrest(2)+1:length(spm_list)] ; 1:length(spm_list)}; % define list for loop, excluding participants with wrong trigger from resting state
out_epochfilename = {'/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/ActualLengthEpoch_rest.mat'; ...
                     '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/ActualLengthEpoch_beat.mat'}; % name for saving output file containing the actual epochlength

% Process event information
timeinfo = cell(length(spm_list),3);
for subi = list4loop{which_cond} % over .mat files
    D = spm_eeg_load([spm_list(subi).folder '/' spm_list(subi).name]); %load spm_list .mat files
    D = D.montage('switch',0);
    %if strcmp(dummy(22:26), 'speed') %checks whether characters 22 to 26 are equal to 'speed'; the loop continues if this is true (1) and it stops if this is false (0)
    events = D.events; %look for triggers
    
    pretrig = -100; % epoch start in ms (prestimulus/trigger)
    %codes for a later (end of this section) detection of (possible) bad trials
    posttrig = (D.time(end) - events(1).time)*1000; % epoch end in ms (posttrigger/stimulus)
    
    if which_cond == 1
        timeinfo{subi,1} = events(end).time - events(1).time - 0.1; %this makes sense only for beat since there are multiple events (the sounds), and not just one (the beginning of the resting state) as in resting state
    end
    timeinfo{subi,2} = D.time(end) - events(1).time - 0.1; %time between end of recording and first event (e.g. first beat)
    timeinfo{subi,3} = spm_list(subi);
    
    % Settings for the epoching
    S2 = [];
    S2.timewin = [pretrig posttrig]; %creating the timewindow of interest
    S2.D = D;
    % event definitions
    
    S2.trialdef(1).conditionlabel = condlabels{which_cond};     % name of the category
    S2.trialdef(1).eventtype = 'STI101_up';         % the event type based on D for you trial/experiment
    S2.trialdef(1).eventvalue = eventvalue(which_cond);                  % event value based on D for your trial/experiment
    
    % Other settings
    S2.prefix = [prefix2add];
    S2.reviewtrials = 0;
    S2.save = 0;
    S2.epochinfo.padding = 0;
    S2.event = D.events;
    S2.fsample = D.fsample;
    S2.timeonset = D.timeonset;
    % function to define trials specifications
    [epochinfo.trl, epochinfo.conditionlabels, MT] = spm_eeg_definetrial(S2);
    % function for the actual epoching
    S3 = [];
    S3 = epochinfo;
    S3.prefix = S2.prefix;
    S3.D = D;
    D = spm_eeg_epochs(S3);    %     D = spm_eeg_load([spm_list(ii).folder '/e' spm_list(ii).name]); %loading epoched data
    D = D.montage('switch',1);
    % from here you go at the end of the if state
    
    D.save();
    disp([num2str(subi)])
end


save(out_epochfilename{which_cond},'timeinfo')

%%

%%

%% *** SOURCE RECONSTRUCTION ***

%% CREATING 8mm PARCELLATION FOR EASIER INSPECTION IN FSLEYES

% NOTE: This section is only executed for the better handling of some visualization purposes, but it does not affect any of the beamforming algorithm;
% In order not to mix up the MNI coordinates, we recommend running the following lines

% Define the number of voxels for the current parcellation
nvoxels = 3559;

% 1) Use load_nii to load a previous NIFTI image
imag_8mm = load_nii('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_T1_8mm_brain.nii.gz');
Minfo = size(imag_8mm.img); % get info about the size of the original image
M8 = zeros(Minfo(1), Minfo(2), Minfo(3)); % initialize an empty matrix with the same dimensions as the original .nii image
counter = 0; % set a counter
M1 = imag_8mm.img;
for xi = 1:Minfo(1) % loop across each voxel for every dimension
    for yi = 1:Minfo(2)
        for zi = 1:Minfo(3)
            if M1(xi,yi,zi) ~= 0 %if we have an actual brain voxel
                counter = counter+1;
                M8(xi,yi,zi) = counter;
            end
        end
    end
end
% 2) Put your matrix in the field ".img"
imag_8mm.img = M8; % assign values to new matrix 
% 3) Save NIFTI image using save_nii
save_nii(imag_8mm,'/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.gz');
% 4) Use FSLEYES to inspect the figure
% Create parcellation on the 8mm template
for voxi = 1:nvoxels %for each 8mm voxel
    cmd = ['fslmaths /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.nii.gz -thr ' num2str(ii) ' -uthr ' num2str(ii) ' -bin /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/AAL_80mm_3559ROIs/' num2str(ii) '.nii.gz'];
    system(cmd)
    disp(voxi)
end
% 5) GET MNI COORDINATES OF THE NEW FIGURE AND SAVE THEM ON DISK
MNI8 = zeros(nvoxels,3); % 3-D matrix
for mm = 1:nvoxels %over brain voxel
    path_8mm = ['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/' num2str(mm) '.nii.gz']; %path for each of the 3559 parcels
    [mni_coord,pkfo] = osl_mnimask2mnicoords(path_8mm);  %getting MNI coordinates
    MNI8(mm,:) = mni_coord; %storing MNI coordinates
end
%  saving on disk
save('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_coord_dyi.mat', 'MNI8');


%% CONVERSION T1 - DICOM TO NIFTI

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/dicm2nii'); % adds path to the dcm2nii folder in osl
MRIsubj = dir('/raw/sorted/MINDLAB2018_MEG-LearningBach-MemoryInformation/0*'); % path with all the subjects folders (shared for both conditions)
MRIoutput = {'/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/MRI_nifti_rest'; ...
             '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/MRI_nifti_beat'}; % this depends on the condition
mkdir(MRIoutput{which_cond});

SUBJ = zeros(length(MRIsubj),1);
for subi = 1:length(MRIsubj) % over subjects
    if strcmp('0022',MRIsubj(subi).name) % SUBJ0022 had corrupted MRI DICOM files, so we provided a template
        niiFolder = [MRIoutput{which_cond} '/' MRIsubj(subi).name];
        copyfile('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_2mm.nii',[niiFolder '/MNI152_T1_2mm.nii'])
    else
        cnt = 0;
        asd = [MRIoutput{which_cond} '/' MRIsubj(subi).name];
        mkdir(asd); %if not, creating it
        if isempty(dir([asd '/*.nii'])) %if there are no nifti images.. I need to convert them
            flagg = 0;
            MRIMEGdate = dir([MRIsubj(subi).folder '/' MRIsubj(subi).name '/20*']);
            niiFolder = [MRIoutput{which_cond} '/' MRIsubj(subi).name];
            for dati = 1:length(MRIMEGdate) % over dates of recording
                if ~isempty(dir([MRIMEGdate(dati).folder '/' MRIMEGdate(dati).name '/MR*'])) %if we get an MRI recording
                    MRI2 = dir([MRIMEGdate(dati).folder '/' MRIMEGdate(dati).name '/MR/*INV2']); %looking for T1
                    if ~isempty(MRI2) %if we have it
                        cnt = cnt + 1;
                        flagg = 1; %determining that I could convert MRI T1
                        dcmSource = [MRI2(1).folder '/' MRI2(1).name '/files/'];
                        %                         if ii ~= 68 || jj ~= 3 %this is because subject 0068 got two MRIs stored.. but the second one (indexed by jj = 3) is of another subject (0086); in this moment, subject 0086 is perfectly fine, but in subject 0068 there are still the two MRIs (for 0068 (jj = 2) and for 0086 (jj = 3))
                        dicm2nii(dcmSource, niiFolder, '.nii');
                        %                         end
                    end
                end
            end
            if flagg == 0
                warning(['subject ' MRIsubj(subi).name ' has no MRI T1 - copying template']);
                copyfile('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_2mm.nii',[niiFolder '/MNI152_T1_2mm.nii'])
            end
            SUBJ(subi) = cnt;
        end
    end
    disp(subi)
end

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% RHINO coregistration

% The algorithm strips the brain from the skull and co-registers it with
% MEG data; ensures that the MRI coordinates are aigned with the 'digital
% pen'. Taken from Oxford (OSL software) by Leonardo Bonetti

% Dir to either Resting state {1} or Beat listening {2} condition
dir_cond = {'/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/es*rest10*.mat';... 
            '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/es*eat*.mat'};        
% Dir to epoched files (encoding), dir to the pre-selected condition
list = dir(dir_cond{which_cond});
% Select as input the previously outputted NIFTI folders 
MRIinput = {'/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/MRI_nifti_rest'; ...
            '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/MRI_nifti_beat'};
        
for subi = 1:length(list) %over subjects
    S = [];
    S.ii = subi;
    S.D = [list(subi).folder '/' list(subi).name]; %path to major files
    D = spm_eeg_load(S.D);
    dummyname = D.fname;
    dummymri = dir([MRIinput{which_cond} '/' dummyname(13:16) '/*.nii']); %path to nifti files
    S.mri = [dummymri(1).folder '/' dummymri(1).name];
    %standard parameters
    S.useheadshape = 1;
    S.use_rhino = 1; %set 1 for rhino, 0 for no rhino
    %         S.forward_meg = 'MEG Local Spheres';
    S.forward_meg = 'Single Shell'; %CHECK WHY IT SEEMS TO WORK ONLY WITH SINGLE SHELL!!
    S.fid.label.nasion = 'Nasion';
    S.fid.label.lpa = 'LPA';
    S.fid.label.rpa = 'RPA';
    jobid = job2cluster(@coregfunc,S); %running with parallel computing
end


%% OPTIONAL: checking (or copying) RHINO

% Dir to either Resting state {1} or Beat listening {2} condition
dir_cond = {'/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/es*rest10*.mat';... 
            '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/es*eat*.mat'};        
% Dir to epoched files (encoding), dir to 
list = dir(dir_cond{which_cond});

for subi = 1:length(list) %over subjects
    D = spm_eeg_load([list(subi).folder '/' list(subi).name]);
    if isfield(D,'inv')
        rhino_display(D)
    end
    disp(['Subject ' num2str(subi)])
end

% NOTE: SUBJ0027 is misplaced after co-registration in Beat Listening condition, therefore we will discard it in the next section

%%

%% BEAMFORMING for anatomical source reconstruction

% Setting up cluster (parallel computing)
% clusterconfig('scheduler', 'none');  % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 2); % slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') % the function being run on the cluster must be contained in this folder

%%

% Loading timeinfo only for Beat Listening condition, regardless of the condition; we are going to take the shortest duration from this condition and apply it to all data
load('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/ActualLengthEpoch_beat.mat'); 

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

% Dir to either Resting state {1} or Beat listening {2} condition for epoched files
dir_cond = {'/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/es*rest10*.mat';... 
            '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/es*eat*.mat'};
% List of subjects with coregistration (RHINO - OSL/FSL) - epoched
list = dir(dir_cond{which_cond}); %dir to epoched files (encoding)
condslabels = {{'Resting_state'},{'Beat_listening'}};
% removing SUBJ0027 for technical reasons in Beat Listening
if which_cond == 2 
    list(24) = []; 
end

if isempty(freqq)
    workingdir = {[workingdir_home '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_broadband_invers_' num2str(invers) '_rest'] ;...
                  [workingdir_home '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_broadband_invers_' num2str(invers) '_beat']};
else
    workingdir = {[workingdir_home '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_' num2str(freqq(1)) '_' num2str(freqq(2)) '_invers_' num2str(invers) '_rest'];...
                  [workingdir_home '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_' num2str(freqq(1)) '_' num2str(freqq(2)) '_invers_' num2str(invers) '_beat']};
end
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');
if ~exist(workingdir{which_cond},'dir') %creating working folder if it does not exist
    mkdir(workingdir{which_cond})
end
for subi = 1:length(list) %over subjects
    S = [];
    if ~isempty(freqq) %if you want to apply the bandpass filter, you need to provide continuous data
        %             disp(['copying continuous data for subj ' num2str(ii)])
        %thus pasting it here
        %             copyfile([list_c(ii).folder '/' list_c(ii).name],[workingdir '/' list_c(ii).name]); %.mat file
        %             copyfile([list_c(ii).folder '/' list_c(ii).name(1:end-3) 'dat'],[workingdir '/' list_c(ii).name(1:end-3) 'dat']); %.dat file
        %and assigning the path to the structure S
        S.norm_megsensors.MEGdata_c = [list(subi).folder '/' list(subi).name(10:end)];
    end
    %copy-pasting epoched files
    %         disp(['copying epoched data for subj ' num2str(ii)])
    %         copyfile([list(ii).folder '/' list(ii).name],[workingdir '/' list(ii).name]); %.mat file
    %         copyfile([list(ii).folder '/' list(ii).name(1:end-3) 'dat'],[workingdir '/' list(ii).name(1:end-3) 'dat']); %.dat file
    
    S.Aarhus_cluster = clust_l; %1 for parallel computing; 0 for local computation
    
    S.norm_megsensors.zscorel_cov = 1; % 1 for zscore normalization; 0 otherwise
    S.norm_megsensors.workdir = workingdir{which_cond};
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
    timepnts = 1:min(cell2mat(timeinfo(:,1)))*250; % vector of time-points, expressed as indexes of samples (based on minimal length over subjects, so that now all subjects have the same length)
    S.inversion.timef = timepnts; % vector of data-points to be extracted; leave it empty [] for working on the full length of the epoch
    S.inversion.conditions = condslabels{which_cond}; %cell with characters for the labels of the experimental conditions (e.g. {'Old_Correct','New_Correct'})
    S.inversion.bc = []; %extreme time-samples for baseline correction (leave empty [] if you do not want to apply it)
    S.inversion.abs = absl; %1 for absolute values of sources time-series (recommendnded 1!)
    S.inversion.effects = invers;
    
    S.smoothing.spatsmootl = 0; %1 for spatial smoothing; 0 otherwise
    S.smoothing.spat_fwhm = 100; %spatial smoothing fwhm (suggested = 100)
    S.smoothing.tempsmootl = 0; %1 for temporal smoothing; 0 otherwise
    S.smoothing.temp_param = 0.01; %temporal smoothing parameter (suggested = 0.01)
    S.smoothing.tempplot = []; %vector with sources indices to be plotted (original vs temporally smoothed timeseries; e.g. [1 2030 3269]). Leave empty [] for not having any plot.
    
    S.nifti = 0; %1 for plotting nifti images of the reconstructed sources of the experimental conditions
    S.out_name = ['SUBJ_' list(subi).name(13:16)]; %name (character) for output nifti images (conditions name is automatically detected and added)
    
    if clust_l ~= 1 %useful  mainly for begugging purposes
        MEG_SR_Beam_LBPD(S);
    else
        jobid = job2cluster(@MEG_SR_Beam_LBPD,S); %running with parallel computing
    end
end


%%

%%

%% GENERALISED EIGENVECTOR DECOMPOSITION (GED) 

% Setting up cluster (parallel computing)
% clusterconfig('scheduler', 'none');  % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 8); % slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia')

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
S.path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/';
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


%% GED assessment (Rest and Listen)

% This blocks enables the visual assessment of GED, for individual
% participants and grand-averages, and produces NIFTI images
% - Power spectrum
% - Eigenspectrum
% - Activations maps (cleaned and saved as NIFTI for visualization in the next block)

clear all
close all

%GED list
list_up = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/GED*');

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


%% GED Assessment (Randomizations of Rest)

clear all
close all

% Randomization settings
rand_names = {'_randLabel', '_randLabelPointWise', '_randTime', '_randAll'}; % to access folders inside the Rest conditon
nrands = length(rand_names);

%GED list
list_rest = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/GED*rest');

% Frequency to inspect
% frex2see = 2.439;   %  stimfrex;
freq_vect_idx = 1:length(list_rest);
flag_figures = 0;

% Select 'UniversitA della strada' approach
streetuni_flag = 1; % when set to 1, processes activation patterns for interpretation

% Loop over folders, by frequency
for frexi = 1:length(freq_vect_idx)
    % Loop over subfolders, by randomization strategy
    for randi = 1:nrands
        % Define path
        path_frex = [list_rest(freq_vect_idx(frexi)).folder '/' list_rest(freq_vect_idx(frexi)).name '/' rand_names{randi}];
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

            disp(['Randomization: ' rand_names{randi} ' - GED assessment of freq idx #' num2str(frexi) ' - participant #' num2str(subi)])
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
            disp(['Randomization' rand_names{randi} ': saving nifti image - component ' num2str(compi) ' frequency ' num2str(list_rest(freq_vect_idx(frexi)).name(16:end))])
            %     if adj_compentr == 1
            %         save_nii(nii,['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1/Test_GED/CompEntrAdj_SubjsAverage_Comp_' num2str(frexi) '_Var_' num2str(round(evals_avg(compi))) '.nii']); %printing image
            %     else
            if outliers_flag ==1
                save_nii(nii,[path_GEDoutput '/SubjsAverage_Comp_' num2str(compi) '_Var_' num2str(round(evals_avg(compi))) 'frex_' num2str(list_rest(freq_vect_idx(frexi)).name(16:end)) '_NO_outliers.nii']); % printing image
            else
                save_nii(nii,[path_GEDoutput '/SubjsAverage_Comp_' num2str(compi) '_Var_' num2str(round(evals_avg(compi))) 'frex_' num2str(list_rest(freq_vect_idx(frexi)).name(16:end)) '.nii']); % printing image
            end
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


%% GED assessment (removing and exporting NIFTI files)

% %delete average GED over subjects in case you need it
% %GED list
% list_up = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/GED*');
% freq_vect_idx = 1:172;
% for ii = 1:length(freq_vect_idx)
%     listss = dir([list_up(freq_vect_idx(ii)).folder '/' list_up(freq_vect_idx(ii)).name '/SubjsAverage*']);
%     for pp = 1:length(listss)
%         delete([listss(pp).folder '/' listss(pp).name])
%     end
% end

% Copy-pasting figures to '/aux' to export them from the server

% For Rest and Listen
% GED List
list_up = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/GED*');
freq_vect_idx = 1:length(list_up);
for ii = 1:length(freq_vect_idx)
    listss = dir([list_up(freq_vect_idx(ii)).folder '/' list_up(freq_vect_idx(ii)).name '/SubjsAverage*']);
    for pp = 1:length(listss)
        copyfile([listss(pp).folder '/' listss(pp).name],['/aux/MINDLAB2023_MEG-AuditMemDement/Mattia/NiftiImages/' listss(pp).name])
    end
end

% For randomizations of Rest
% Randomization settings
rand_names = {'_randLabel', '_randLabelPointWise', '_randTime', '_randAll'}; % to access folders inside the Rest conditon
nrands = length(rand_names);
%GED list
list_rest = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/GED*rest');
freq_vect_idx = 1:length(list_rest);
for ii = 1:length(freq_vect_idx)
    for randi = 1:nrands
        listss = dir([list_rest(freq_vect_idx(ii)).folder '/' list_rest(freq_vect_idx(ii)).name '/' rand_names{randi} '/SubjsAverage*']);
        for pp = 1:length(listss)
            copyfile([listss(pp).folder '/' listss(pp).name],['/aux/MINDLAB2023_MEG-AuditMemDement/Mattia/NiftiImages/' listss(pp).name(1:end-4) rand_names{randi} listss(pp).name(end-3:end) ])
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
path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/';
path_suffix = {'Hz_rest', 'Hz_beat'};
rand_names = {'_randLabel', '_randLabelPointWise', '_randTime', '_randAll'}; % to access folders inside the Rest conditon
nrands = length(rand_names);

% Select the index of 
%which_lfo  = 1; % the component to use as low-frequency modulator
which_comp = 1; % the index of the component for the high-frequency carrier

% Load beat listening data
list_beat = dir(['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/GEDoutput_frex_' num2str(stimfrex) 'Hz_beat/SUBJ*.mat']); % do it for stimfrex, but it's the same for any frequency
% Get indexes of actual subjects (to match the other list)
subs_beat = zeros(length(list_beat),1);
for subi = 1:length(list_beat)
    % Retrieve subject
    this_sub = regexp(list_beat(subi).name, '([\d.]+).mat', 'tokens');
    % Assign frequency to vector
    subs_beat(subi) = str2double(this_sub{1,1});
end

% Load resting state data
list_rest = dir(['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/GEDoutput_frex_' num2str(stimfrex) 'Hz_rest/SUBJ*.mat']); % do it for stimfrex, but it's the same for any frequency
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
                evals_all_rands = zeros(nevals,nsubs,nfrex,nrands); % for randomizations of Rest condition
            end
            
            evals_all(:,subi,frexi,condi) = evals; % assign for computing average and for exporting
            
            % If inside Rest condition, assign randomizations
            if condi == 1
                for randi = 1:nrands
                    % Get eigenvalues from relevant folder (already sorted and normalized)
                    load([this_folder, rand_names{randi}, '/', list_subs(subi).name],'evals')
                    evals_all_rands(:,subi,frexi,randi) = evals;
                end
            end

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

% % Re-sort frequencies
evals_all_rands = evals_all_rands(:,:,idx_sort,:);
% Visualize randomizations of Rest
figure(1300), clf
for randi = 1:nrands
    % compute average
    topEvals_avg_rand = squeeze( mean(evals_all_rands(:,:,:,randi),2) );
    subplot(nrands,1,randi), hold on
    plot(all_frex,topEvals_avg_rand(comps2see,:),'s-','markerfacecolor','k','linew',2)
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    xticks(round(all_frex,1)), xtickangle(45)
    %xlim([0 nfrex+1])
    ylim([0 30])
    xlabel('Frequencies (Hz)')
    ylabel('Explained variance')
%     legendEntries = cell(length(comps2see),1);
%     for legi = 1:length(legendEntries)
%        legendEntries{legi} = ['Component #' num2str(comps2see(legi))];
%     end
%     legend(legendEntries, 'Location', 'Best');
    title([rand_names{randi}, ': top ' num2str(length(comps2see)) ' components'])
%     subplot(312), hold on
%     plot(all_frex,sum(topEvals_avg,1),'s-','color','k','markerfacecolor','k','linew',2)
%     line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
%     xticks(all_frex), xtickangle(45)
%     xlim([0 nfrex+1])
%     xlabel('Frequencies (Hz)')
%     ylabel('Explained variance')
%     title(['Sum of top ' num2str(length(comps2see)) ' components'])
%     subplot(313), hold on
%     plot(all_frex,std(topEvals_avg,0,1),'s-','color','r','markerfacecolor','k','linew',2)
%     line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
%     xticks(all_frex), xtickangle(45)
%     xlim([0 nfrex+1])
%     % ylim([0 50])
%     xlabel('Frequencies (Hz)')
%     ylabel('Explained variance')
%     title(['STD of top ' num2str(length(comps2see)) ' components'])
end
    
%% Statistical testing of GED components eigenvalues (explained variance)

clear all
close all

comps2test = [1:20];
stimfrex = 2.439;
load('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Images/GEDallsubjs.mat'); %loading GED components eigenvalues
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
    
    % Visualize eigenvalues
    figure(1000 + compi), clf
    subplot(311)
    plot(all_frex, squeeze( evals_all(compi,:,:,1) ) , 'k'), hold on
    plot(all_frex, squeeze( evals_all(compi,:,:,2) ) , 'r')
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    subplot(312)
    plot(all_frex, Tevals(compi,:)), hold on
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    subplot(313)
    plot(all_frex, Pevals(compi,:)), hold on
    line([stimfrex stimfrex], get(gca, 'Ylim'), 'color', 'r', 'LineStyle', '--')
    ylim([0 1])
end


%FDR correction
PFDR = zeros(size(Pevals,1),1); %FDR thresholds (independently for each GED component)
P_thresh = zeros(size(Pevals)); %frequencies (variance explained) which are significant after FDR correction
for ii = 1:size(Pevals,1) %over GED components
    PFDR(ii,1) = fdr(Pevals(ii,:)); %FDR correction (getting the threshold)
    P_thresh(ii,Pevals(ii,:)<PFDR(ii,1)) = 1; %getting the significant frequencies (independently for each component)
end


%% Cross-frequency coupling (PAC - Phase to Amplitude Coupling)


clear all
close all
clc

% Select the index of the component to use as low-frequency modulator
which_lfo = 2;

% Select components of interest
comps2see = [1 2 3];


% Settings for low-frequency modulator
lfo_frex  = 2.439; % modulator's frequency
lfo_fwhm  = .35;

% Settings for filter
srate = 250;
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia'); % where the function is

suffix{1} = '_rest'; suffix{2} = '_beat';
% Set home path
path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED';
% Get all_frex from home path
list_GED = dir([path_GEDoutput '/GED*']); % keep only folders

load('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/Subjs_idx_match.mat');
for condi = 1:2
    list_GED_cond = list_GED(endsWith({list_GED.name}, ['Hz', suffix{condi}])); % only ending with 'Hz'
    % Initialize vectors for storing frequencies
    nfrex = length(list_GED_cond);
    all_frex = zeros(nfrex,1);
    % Retrieve and store the frequency
    for frexi = 1:numel(list_GED_cond)
        % Check if the folder name contains 'Hz'\
        this_frex = regexp(list_GED_cond(frexi).name, '([\d.]+)Hz', 'tokens');
        all_frex(frexi) = str2double (this_frex{1,1});
    end
    % Re-sort vector of frequencies (strings follow alphabetic order)
    [all_frex, idx_sort] = sort(all_frex, 'ascend');
    % Drop all frequencies below the modulator
    all_frex(all_frex < lfo_frex) = [];
    nfrex = length(all_frex); % Re-compute N
    % Re-compute fwhm for all frequencies (NOTE: should be the same used to perform GED)
    % all_fwhm = logspace(log10(lfo_fwhm), log10(lfo_fwhm * (nfrex)), nfrex); % Based on Cohen ANT course
    all_fwhm = linspace(lfo_fwhm, lfo_fwhm * nfrex, nfrex); % based on Brainstorm, fwhm for Morlet wavelet is linear function of frequency (~17 is max value for 60 hz) (Morillon uses this)
    
    
    % Initialize CFC
    PAC_avg = zeros(nfrex,length(comps2see));% avg of modulation index across subjects, per frequency
    nbins = 37; % set bins for MPA
    if condi == 1
        MPA_avg = zeros([size(PAC_avg),nbins,2]);
    end
    
    % Import and process low-frequency modulator (only once, per subjects)
    this_folder =  [path_GEDoutput, '/GEDoutput_frex_' num2str(lfo_frex) 'Hz' suffix{condi} '/']; % Get folder name
    % Compute phase timeseries from subjects
    list_subs_all = dir(fullfile(this_folder, '*.mat')); %all subjects
    list_subs = list_subs_all(idx4loop{condi}); %only subjects who match between conditions
    nsubs = length(list_subs);
    for subi = 1:nsubs % over subjects
        
        % Get components timeseries (already sorted and normalized)
        load([this_folder, list_subs(subi).name],'GEDts')
        lfo_ts  = GEDts(comps2see,:); % only keep relevant components
        % Filter timeseries
        lfo_ts = filterFGx (lfo_ts, srate, lfo_frex, lfo_fwhm, 0);
        % Initialize matrices for all participants, once upon first iteration
        if subi == 1
            lfo_phase_all = zeros([size(lfo_ts), nsubs]);
        end
        % Compute phase timeseries
        lfo_phase_all(:,:,subi) = angle( hilbert(lfo_ts') )'; % Hilbert operates along the columns ('vertical' input, time x chans)
        
    end
    % Compute edges to define phase bins (for MPA) - keep soft-coded, in case of re-wrapping
    phase_edges = linspace(min(lfo_phase_all(:)),max(lfo_phase_all(:)),nbins+1);
    
    % Import and process carriers, PAC, and periodogram, per frequency
    for frexi = 1:nfrex
        
        % Get folder name
        this_folder =  [path_GEDoutput,'/GEDoutput_frex_' num2str(all_frex(frexi)) 'Hz' suffix{condi} '/'];
        
        % Assign frequency to this carrier
        carr_frex = all_frex(frexi);
        % Assigned optimized filter width to this carrier
        carr_fwhm = all_fwhm(frexi);  % lfo_fwhm;
        
        % Compute CFC from subjects
        list_subs_all = dir(fullfile(this_folder, '*.mat')); %all subjects
        list_subs = list_subs_all(idx4loop{condi}); %only subjects who match between conditions
        nsubs = length(list_subs);
        for subi = 1:nsubs % over subjects
            
            % Get components timeseries (already sorted and normalized)
            load([this_folder, list_subs(subi).name],'GEDts')
            carr_ts  = GEDts(comps2see,:); % only keep relevant components
            % Filter timeseries
            carr_ts = filterFGx (carr_ts, srate, carr_frex, carr_fwhm, 0);
            % Initialize matrices for all participants, once upon first iteration
            if subi == 1 && frexi == 1 && condi == 1
                carr_pow_all = zeros([size(carr_ts), nsubs]);
                PAC_all = zeros(nfrex,length(comps2see),nsubs);
                MPA_all = zeros([size(PAC_all),nbins,2]);
            end
            
            % Compute power timeseries
            carr_pow_all(:,:,subi)  = abs( hilbert(carr_ts') )' .^2; % Hilbert operates along the columns ('vertical' input, time x chans)
            
            % Copute PAC and MPA
            phase = wrapToPi( lfo_phase_all(which_lfo,:,subi) ); % the modulator is alwas the same
            for compi = comps2see
                
                % Assignment to temporary variables, for better readability
                pow   = carr_pow_all(compi,:,subi);
                %             pow   = 10*log10(pow); % in Morillon's pipeline... but watch
                %             %out if you end up computing sqrt of negative numbers...
                % Compute phase-amplitue coupling (PAC) with Morillon's formula
                %PAC_all(frexi,compi,subi) = abs( sum( sqrt(pow) * exp(1i*phase') ) ) ./ sqrt( sum(pow) * length(pow) ); % normalized [0 1]
                %PAC_all(frexi,compi,subi) = abs( sqrt(pow) * exp(1i*phase') )  ./ sum( sqrt(pow) ) ; % normalized [0 1]
                % from Zalta et al.(2024)
                PAC_all(frexi,compi,subi) = abs(sum( sqrt(pow).*exp(1i*phase) )) / (length(pow)*sqrt(sum(pow)));
                
                % Compute movement phase-related amplitudes (MPA) (from Seeber
                % et al., 2016; Rosso et al., 2022)
                for bini = 1:nbins
                    
                    % Compute modulation over all run (no issues related to moving mean or median)
                    MPA_all(frexi,compi,subi,bini,condi) = mean(pow(phase>phase_edges(bini) & phase<phase_edges(bini+1)) , 'omitnan');
                    
                end
            end            
        end
        
        % Compute average CFC measures subjects, per frequency
        PAC_avg(frexi,:) = squeeze( mean(PAC_all(frexi,:,:),3) );
        MPA_avg(frexi,:,:,condi) = squeeze( mean(MPA_all(frexi,:,:,:,condi),3) );
        
        disp(['Computing power modulation of frex ' num2str(carr_frex) 'Hz by ' num2str(lfo_frex) 'Hz'])
    end
    
    
    % Sinusoudal fit for MPA
    % Initialize output
    % Sine amplitude, phase shift, DC offset (aka intercept, aka sine mean value)
    if condi == 1
        [sAmpl, mFrex, pShift , dcOff, goodFit] = deal(zeros(nfrex,nsubs,length(comps2see),2)); % fit to grand-average
    end
    
    for compi = comps2see
        smth = 1; % smoothing factor
        for frexi = 1:nfrex
            for subi = 1:nsubs
                
                disp(['SineFit for MPA - Frex = ' num2str(all_frex(frexi)) ' Hz - Sub#' num2str(subi)])
                
                %Initialize temporal variable (reset within loop)
%                 tempParam = [];
                % Compute sine parameters (offset, amplitude, frequency,
                % phaseshift, MSE)
                x = 1:nbins;
                y = smooth( squeeze(MPA_all(frexi,compi,subi,:,condi)), smth)';
                [tempParam] = sineFit(x,y, 0); % last argument suppressess output figure
                % Assign amplitude to output matrix
                sAmpl(frexi,subi,compi,condi)   = tempParam(2);
                mFrex(frexi,subi,compi,condi)   = tempParam(3);
                pShift(frexi,subi,compi,condi)  = tempParam(4);
                dcOff(frexi,subi,compi,condi)   = tempParam(1);
                goodFit(frexi,subi,compi,condi) = tempParam(5);
                
                
            end
        end
    end
end
% Set Zeros To ?NaN?
sAmpl(sAmpl   == 0) = NaN;
pShift(pShift == 0) = NaN;
dcOff(dcOff   == 0) = NaN;
% Save outputs
save([path_GEDoutput '/CFC_by_' num2str(lfo_frex) '_comp#' num2str(which_lfo) '.mat'], 'PAC_all','PAC_avg','MPA_all','MPA_avg','sAmpl','pShift','dcOff','mFrex','goodFit','all_frex','phase_edges');

%% Visualize CFC

clear all
close all
clc

% Select the index of the component to use as low-frequency modulator
which_lfo = 1;
lfo_frex  = 2.439; % modulator's frequency
condition = 1; %1 = resting state; 0 = beat


% Set home path
path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/';

% Load CFC data
load([path_GEDoutput '/CFC_by_' num2str(lfo_frex) '_comp#' num2str(which_lfo) '.mat']);


% Visualize PAC by stimulation frequency
frexbounds = [all_frex(2) all_frex(end)+1]; % start from 2 to skip the modulator
comps2plot = [1 2]; %in case you want a subset
nfrex = length(all_frex);

figure(500), clf
subplot(311), hold on
plot(all_frex, PAC_avg(:,comps2plot)','s-','markerfacecolor','k','linew',2)
xticks(round(all_frex,1)), xtickangle(45)
xlim(frexbounds)
legendEntries = cell(length(comps2plot),1);
for legi = 1:length(legendEntries)
   legendEntries{legi} = ['Component #' num2str(comps2plot(legi))];
end
legend(legendEntries, 'Location', 'Best');
title(['PAC by stimulation frequency (component #' num2str(which_lfo) ')'])
xlabel('Frequencies (Hz)')
ylabel('R')
subplot(312), hold on
plot(all_frex,sum(PAC_avg(:,comps2plot),2),'s-','color','k','markerfacecolor','k','linew',2)
xticks(all_frex), xtickangle(45)
xlim(frexbounds) 
xlabel('Frequencies (Hz)')
ylabel('R')
title(['Sum of top ' num2str(length(comps2plot)) ' components'])
subplot(313), hold on
plot(all_frex,std(PAC_avg(:,comps2plot),0,2),'s-','color','r','markerfacecolor','k','linew',2)
xticks(all_frex), xtickangle(45)
xlim(frexbounds)
xlabel('Frequencies (Hz)')
ylabel('R')
title(['STD of top ' num2str(length(comps2plot)) ' components'])


% Visualize MPA by stimulation frequency
phasebounds = [phase_edges(1) phase_edges(end)]; 
comps2plot = [1]; %in case you want a subset
frex2plot = dsearchn(all_frex,10.5); % target the one expected to be modulated
colvec = zeros(nfrex,3); % rgb vector
colvec(:,1) = linspace(0,1,nfrex);
colvec(:,2) = linspace(1,0,nfrex);
%offset_frex = sort( repmat([1:nfrex], length(phasebounds)) );
smth = 1;
nsubs = size(MPA_all,3);

% Individual participants
x = phase_edges(1:end-1);
for subi = 1:nsubs
    figure(600+subi),clf,hold on
    for frexi = 1:nfrex   
        y = smooth( squeeze(MPA_all(frexi,comps2plot,subi,:,condition))' - squeeze(mean(MPA_all(frexi,comps2plot,subi,:,condition),4,'omitnan')), smth )';
        %y = y(end:-1:1);
        z = all_frex(frexi) * ones(size(x));
        %z = (nfrex - (frexi-1)) * ones(size(x));
        if frexi == frex2plot
            plot3(x,y,z, '*-','color',colvec(frexi,:),'linew',2.6)
        else
            plot3(x,y,z, '-','color',colvec(frexi,:),'linew',1.6)
       end
    end
    legendEntries = cell(nfrex,1);
    for legi = 1:length(legendEntries)
       legendEntries{legi} = [num2str(all_frex(legi)) ' Hz'];
    end
    legend(legendEntries, 'FontSize',6, 'Location', 'Best');
    grid on, rotate3d on
    view(0,50+180)
    %xticklabels(phase_edges)
    xtickangle(45)
    xlabel('Phase (rad)')
    ylabel('Power (/mu^2)')
    title(['MPA'])

end

% Grand-average
figure(600), clf, hold on
x = phase_edges(1:end-1);
for frexi = 1:nfrex   
    y = smooth( squeeze(MPA_avg(frexi,comps2plot,:,condition))' - squeeze(mean(MPA_avg(frexi,comps2plot,:,condition),3,'omitnan')), smth )';
    %y = y(end:-1:1);
    z = all_frex(frexi) * ones(size(x));
    %z = (nfrex - (frexi-1)) * ones(size(x));
    if frexi == frex2plot
        plot3(x,y,z, '*-','color',colvec(frexi,:),'linew',2.6)
    else
        plot3(x,y,z, '-','color',colvec(frexi,:),'linew',1.6)
   end
end
legendEntries = cell(nfrex,1);
for legi = 1:length(legendEntries)
   legendEntries{legi} = [num2str(all_frex(legi)) ' Hz'];
end
legend(legendEntries, 'FontSize',6, 'Location', 'Best');
grid on, rotate3d on
view(0,50+180)
%xticklabels(phase_edges), 
xtickangle(45)
xlabel('Phase (rad)')
ylabel('Power (/mu^2)')
title(['MPA'])

% MPA as surface
figure(700), clf
x = phase_edges(1:end-1)';
y = all_frex(1:end);
Z = squeeze(MPA_avg(1:length(y),comps2plot,:,condition))' - squeeze(mean(MPA_avg(1:length(y),comps2plot,:,condition),3))';
% for frexi = 1:length(y)
%     Z(:,comps2plot,frexi) = smooth(Z(frexi,comps2plot,:),smth);
% end
surf(x,y,Z')
zlim(max(Z(:))*[-1 1])
% shading interp
rotate3d on, axis square
xlabel('Phase'), ylabel('Frequency (Hz)'), zlabel('Power')

% Sinusoidal fit 
for subi = 1:nsubs
figure(800+subi), clf, hold on
subplot(411)
plot(all_frex,sAmpl(:,subi,comps2plot,condition),'s-')
ylim([0 0.1])
title('Modulation amplitude')
subplot(412)
plot(all_frex,pShift(:,subi,comps2plot,condition),'s-')
title('Modulation phase')
subplot(413)
plot(all_frex,mFrex(:,subi,comps2plot,condition),'s-')
title('Modulation frequency') % 1 cycle is expected
subplot(414)
plot(all_frex,goodFit(:,subi,comps2plot,condition),'s-')
title('Goodness of fit') 
end

figure(800+nsubs+1), clf, hold on
subplot(411)
plot(all_frex,mean(sAmpl(:,:,comps2plot,condition),2),'s-')
ylim([0 0.1])
title('Modulation amplitude')
subplot(412)
plot(all_frex,mean(pShift(:,:,comps2plot,condition),2),'s-')
title('Modulation phase')
subplot(413)
plot(all_frex,mean(mFrex(:,:,comps2plot,condition),2),'s-')
title('Modulation frequency') % 1 cycle is expected
subplot(414)
plot(all_frex,mean(goodFit(:,:,comps2plot,condition),2),'s-')
title('Goodness of fit') 


%% Statistical tests fro CFC
clear all
close all

% Select the index of the component to use as low-frequency modulator
which_lfo = 2;
lfo_frex  = 2.439; % modulator's frequency

% Set home path
path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Paper_GED/';
% Import data
load([path_GEDoutput '/CFC_by_' num2str(lfo_frex) '_comp#' num2str(which_lfo) '.mat'], 'PAC_all','PAC_avg','MPA_all','MPA_avg','sAmpl','pShift','dcOff','mFrex','goodFit','all_frex','phase_edges');

% Tests
[P,T] = deal(zeros(size(sAmpl,1),size(sAmpl,3)));
for frexi = 1:size(sAmpl,1)
    for compi = 1:size(sAmpl,3)
        %     [~,p,~,stats] = ttest( squeeze(sAmpl_beat(frexi,:)) , squeeze(sAmpl_rest(frexi,:)), 'Tail','right');
        %     P(frexi) = p;
        %     T(frexi) = stats.tstat;
        
        %     [p,~,stats] = signrank( squeeze(sAmpl(frexi,:,compi,2)) , squeeze(sAmpl_rest(frexi,:,compi,1)), 'Tail','right');
        [p,~,stats] = signrank( squeeze(sAmpl(frexi,:,compi,2)) , squeeze(sAmpl(frexi,:,compi,1)),'Tail','right');
        P(frexi,compi) = p;
        T(frexi,compi) = stats.zval;
    end
end


figure
plot(squeeze(mean(sAmpl(:,:,:,2),2,'omitnan'))), 
figure
plot(squeeze(median(sAmpl(:,:,:,1),2,'omitnan')))
figure
plot(squeeze(median(sAmpl(:,:,:,2),2,'omitnan'))-squeeze(median(sAmpl(:,:,:,1),2,'omitnan')))
figure
T(1:10,:) = 0;
plot(T)
figure
plot(P)

%%

%%
