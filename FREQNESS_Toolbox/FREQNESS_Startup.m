function [allData, MNI, path_home, folder_data] = FREQNESS_Startup(load_data)
% ========================================================================
%  FREQNESS STARTUP FUNCTION
% ========================================================================
%
%  If you use this toolbox, please cite:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ-NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved
%  Brain Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%
%  This function initializes the FREQNESS toolbox by automatically
%  determining its installation directory and adding all required
%  subfolders to the MATLAB path.
%
%  By default, all participant matrices are loaded into allData, with one
%  cell per Group or Condition. Each matrix has dimensions:
%
%      voxels x time x participants
%
%  Set load_data to false to inspect participant files without loading
%  their data matrices. In this case, allData contains empty cells and the
%  fourth output, folder_data, lists the folders and their files_mat entries.
%  The data_size field contains the common [nVoxels nTime] size per folder.
%  MNI coordinates and toolbox paths are initialized in both modes.
%
%  Examples:
%      [allData,MNI,path_home] = FREQNESS_Startup();
%      [~,MNI,path_home,folder_data] = FREQNESS_Startup(false);
%
%  AUTHORS:
%  Mattia Rosso & Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%
% ========================================================================

%% Select whether to load participant data

if nargin < 1
    load_data = true;
end
if ~islogical(load_data) || ~isscalar(load_data)
    error('load_data must be one logical value (true or false).');
end

%% Define main path based on the directory of this function
path_home = fileparts(mfilename('fullpath'));

%% Add toolbox directories
addpath(path_home);
addpath(fullfile(path_home, 'FREQNESS_Functions'));
addpath(fullfile(path_home, 'FREQNESS_ExternalFunctions'));
addpath(fullfile(path_home, 'FREQNESS_ExternalFunctions', 'nifti_tools'));
addpath(fullfile(path_home, 'FREQNESS_Data'));
addpath(fullfile(path_home, 'FREQNESS_MNI_Coordinates'));

if load_data
    fprintf('\nFREQNESS Startup: loading data and MNI coordinates.\n');
else
    fprintf('\nFREQNESS Startup: inspecting data files and loading MNI coordinates.\n');
end
fprintf('\nFREQNESS successfully initialized.\n');
fprintf('Base directory: %s\n\n', path_home);

%% Fetch data from experimental Groups or Conditions

% Find folders inside FREQNESS_Data
folder_data = dir(fullfile(path_home, 'FREQNESS_Data'));
folder_data = folder_data([folder_data.isdir] & ~ismember({folder_data.name},{'.','..'}));
if ~isempty(folder_data)
    [~,folder_order] = sort(lower({folder_data.name}));
    folder_data = folder_data(folder_order);
end

% Prepare container for group / condition data
nfolds  = numel(folder_data);
allData = cell(nfolds,1);
nsubs   = zeros(nfolds,1);   % number of subjects per folder

for foldi = 1:nfolds

    % Current group/condition folder
    this_folder = fullfile(folder_data(foldi).folder, folder_data(foldi).name);

    % List .mat files (one per participant)
    files_mat = dir(fullfile(this_folder, '*.mat'));
    if ~isempty(files_mat)
        [~,file_order] = sort(lower({files_mat.name}));
        files_mat = files_mat(file_order);
    end
    folder_data(foldi).files_mat = files_mat;
    folder_data(foldi).data_size = [];
    nsubs(foldi) = numel(files_mat);
    % Check whether the files are missing
    if nsubs(foldi) == 0
        warning('No .mat files found in folder: %s', this_folder);
        allData{foldi} = [];
    else

        % Inspect or load each participant in this folder
        for subi = 1:nsubs(foldi)

            if load_data
                % Load .mat file (must contain exactly one data matrix)
                S  = load(fullfile(files_mat(subi).folder, files_mat(subi).name));
                fn = fieldnames(S);
                if numel(fn) ~= 1
                    error('File %s must contain exactly one data matrix.', files_mat(subi).name);
                end
                this_data = S.(fn{1});   % assign voxels-by-time data matrix

                if ~isnumeric(this_data) || ~isreal(this_data) || isempty(this_data) || ...
                        ~ismatrix(this_data) || any(~isfinite(this_data(:)))
                    error(['File %s must contain one non-empty real numeric ' ...
                        'voxels-by-time matrix with finite values.'],files_mat(subi).name);
                end
                data_size = size(this_data);
            else
                % Read dimensions and class without importing the data matrix
                file_info = whos('-file',fullfile(files_mat(subi).folder,files_mat(subi).name));
                if numel(file_info) ~= 1
                    error('File %s must contain exactly one data matrix.',files_mat(subi).name);
                end
                numeric_classes = {'double','single','int8','uint8','int16', ...
                    'uint16','int32','uint32','int64','uint64'};
                if ~ismember(file_info.class,numeric_classes) || file_info.complex || ...
                        numel(file_info.size) ~= 2 || any(file_info.size == 0)
                    error(['File %s must contain one non-empty real numeric ' ...
                        'voxels-by-time matrix.'],files_mat(subi).name);
                end
                data_size = file_info.size;
            end

            % On first subject, store dimensions and optionally allocate data
            if subi == 1
                nvoxs = data_size(1);
                ntime = data_size(2);
                folder_data(foldi).data_size = data_size;
                if load_data
                    groupData = nan(nvoxs, ntime, nsubs(foldi));
                end
            else
                % Basic consistency check
                if ~isequal(data_size, [nvoxs ntime])
                    error('Size mismatch in file %s relative to first subject in folder %s.', ...
                        files_mat(subi).name, this_folder);
                end
            end

            % Assign into group matrix along 3rd dimension
            if load_data
                groupData(:,:,subi) = this_data;
            end
        end

        % Store in cell array (one entry per group / condition)
        if load_data
            allData{foldi} = groupData;
        end

    end

end



%% Fetch MNI coordinates (if available)

% Initialize default: empty
MNI = [];

% Find .mat file containing the MNI coordinates
folder_mni = fullfile(path_home, 'FREQNESS_MNI_Coordinates');
files_mni  = dir(fullfile(folder_mni, '*.mat'));
if ~isempty(files_mni)
    [~,mni_file_order] = sort(lower({files_mni.name}));
    files_mni = files_mni(mni_file_order);
end

% Proceed only if exactly one file is found
if numel(files_mni) == 1

    S  = load(fullfile(files_mni(1).folder, files_mni(1).name));
    fn = fieldnames(S);

    % Must contain exactly one variable
    if numel(fn) == 1
        temp_mni = S.(fn{1});

        % Must contain a real numeric XYZ coordinate matrix
        if ~isnumeric(temp_mni) || ~isreal(temp_mni) || isempty(temp_mni) || ...
                any(isinf(temp_mni(:)))
            warning(['MNI file %s ignored: coordinates must be provided as ' ...
                'a real numeric matrix without infinite values.'],files_mni(1).name);
        elseif size(temp_mni,2) == 3
            MNI = temp_mni;
        elseif size(temp_mni,1) == 3 && size(temp_mni,2) ~= 3
            MNI = temp_mni'; % transpose if has coordinates along the rows
        else
            warning('MNI file %s ignored: coordinate matrix does not have 3 columns.', ...
                files_mni(1).name);
        end
    else
        warning('MNI file %s ignored: file must contain exactly one variable, consisting of a 3-column coordinate matrix.', ...
            files_mni(1).name);
    end

elseif numel(files_mni) > 1
    warning('Multiple MNI .mat files found. None loaded.');

end

% Ensure that the coordinate grid matches every non-empty dataset
if ~isempty(MNI)
    data_nvoxs = [];
    for foldi = 1:nfolds
        if ~isempty(folder_data(foldi).data_size)
            data_nvoxs(end+1) = folder_data(foldi).data_size(1); %#ok<AGROW>
        end
    end

    if ~isempty(data_nvoxs) && any(data_nvoxs ~= size(MNI,1))
        warning(['MNI coordinates ignored: the coordinate file contains %d rows, ' ...
            'but the datasets contain %s voxels. Replace the coordinate ' ...
            'file with one matching the voxel order and source space of your data.'], ...
            size(MNI,1),mat2str(unique(data_nvoxs)));
        MNI = [];
    end
end


end
