function app = FREQNESS_GUI(varargin)
%FREQNESS_GUI Launch the graphical interface for the FREQ-NESS Toolbox.
%
%   FREQNESS_GUI opens the interface.
%   APP = FREQNESS_GUI returns the application handle.
%   FREQNESS_GUI(DATASET_FOLDER) opens the interface and imports a Dataset*
%   folder immediately.

guiRoot = fileparts(mfilename('fullpath'));
toolboxRoot = fileparts(guiRoot);

requiredPaths = {
    guiRoot
    fullfile(guiRoot,'thirdparty','uiFileDnD')
    fullfile(toolboxRoot,'FREQNESS_Functions')
    fullfile(toolboxRoot,'FREQNESS_ExternalFunctions')
    fullfile(toolboxRoot,'FREQNESS_ExternalFunctions','nifti_tools')
    };

for pathi = 1:numel(requiredPaths)
    if isfolder(requiredPaths{pathi})
        addpath(requiredPaths{pathi});
    end
end

% uiFileDnD requires independently hosted figures in newer MATLAB releases.
% Defining the property here keeps the GUI self-contained and avoids edits to
% a user or toolbox startup file.
try
    addprop(groot,'ForceIndependentlyHostedFigures');
catch
    % The property already exists, or this MATLAB release does not need it.
end

app = FREQNESSApp(varargin{:});

end
