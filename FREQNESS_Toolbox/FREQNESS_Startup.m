function [allData, MNI, path_home] = FREQNESS_Startup()
% ========================================================================
%  FREQNESS STARTUP FUNCTION
% ========================================================================
%
%  If you find this function useful, please cite the first FREQNESS paper:
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
% ========================================================================

%% Define main path based on the directory of this function
path_home = fileparts(mfilename('fullpath'));

%% Add toolbox directories
addpath(path_home);
addpath(fullfile(path_home, 'FREQNESS_Functions'));
addpath(fullfile(path_home, 'FREQNESS_ExternalFunctions'));
addpath(fullfile(path_home, 'FREQNESS_ExternalFunctions', 'nifti_tools'));
addpath(fullfile(path_home, 'FREQNESS_Data'));
addpath(fullfile(path_home, 'FREQNESS_MNI_Coordinates'));

fprintf('\nFREQNESS Startup: loading data and MNI coordinates.\n');
fprintf('\nFREQNESS successfully initialized.\n');
fprintf('Base directory: %s\n\n', path_home);

%% Fetch data from experimental Groups or Conditions

% Find folders inside FREQNESS_Data
folder_data = dir(fullfile(path_home, 'FREQNESS_Data'));
folder_data = folder_data([folder_data.isdir] & ~ismember({folder_data.name},{'.','..'}));

% Prepare container for group / condition data
nfolds  = numel(folder_data);
allData = cell(nfolds,1);
nsubs   = zeros(nfolds,1);   % number of subjects per folder

for foldi = 1:nfolds

    % Current group/condition folder
    this_folder = fullfile(folder_data(foldi).folder, folder_data(foldi).name);

    % List .mat files (one per participant)
    files_mat = dir(fullfile(this_folder, '*.mat'));
    nsubs(foldi) = numel(files_mat);
    % Check whether the files are missing
    if nsubs(foldi) == 0
        warning('No .mat files found in folder: %s', this_folder);
        allData{foldi} = [];
    else

        % Store all participants in this folder
        for subi = 1:nsubs(foldi)

            % Load .mat file (must contain exactly one data matrix)
            S  = load(fullfile(files_mat(subi).folder, files_mat(subi).name));
            fn = fieldnames(S);
            if numel(fn) ~= 1
                error('File %s must contain exactly one data matrix.', files_mat(subi).name);
            end
            this_data = S.(fn{1});   % assign voxels-by-time data matrix

            % On first subject, initialize 3D array
            if subi == 1
                [nvoxs, ntime] = size(this_data);
                groupData = nan(nvoxs, ntime, nsubs(foldi));
            else
                % Basic consistency check
                if ~isequal(size(this_data), [nvoxs ntime])
                    error('Size mismatch in file %s relative to first subject in folder %s.', ...
                        files_mat(subi).name, this_folder);
                end
            end

            % Assign into group matrix along 3rd dimension
            groupData(:,:,subi) = this_data;
        end

        % Store in cell array (one entry per group / condition)
        allData{foldi} = groupData;

    end

end



%% Fetch MNI coordinates (if available)

% Initialize default: empty
MNI = [];

% Find .mat file containing the MNI coordinates
folder_mni = fullfile(path_home, 'FREQNESS_MNI_Coordinates');
files_mni  = dir(fullfile(folder_mni, '*.mat'));

% Proceed only if exactly one file is found
if numel(files_mni) == 1

    S  = load(fullfile(files_mni(1).folder, files_mni(1).name));
    fn = fieldnames(S);

    % Must contain exactly one variable
    if numel(fn) == 1
        temp_mni = S.(fn{1});

        % Must have 3 columns (XYZ coordinates)
        if size(temp_mni,2) == 3
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


end
