function mniFile = findDefaultMNI(toolboxFolder)
%FINDDEFAULTMNI Locate the bundled MNI coordinate MAT-file, if present.

if ~(ischar(toolboxFolder) || ...
        (isstring(toolboxFolder) && isscalar(toolboxFolder)))
    error('FREQNESS:GUI:InvalidToolboxFolder', ...
        'The toolbox folder must be one text scalar.');
end

mniFolder = fullfile(char(toolboxFolder),'FREQNESS_MNI_Coordinates');
mniFile = '';
if ~isfolder(mniFolder)
    return
end

files = dir(fullfile(mniFolder,'*.mat'));
if isempty(files)
    return
end

preferredName = 'MNI152_8mm_coord_dyi.mat';
preferredIndex = find(strcmpi({files.name},preferredName),1);
if ~isempty(preferredIndex)
    mniFile = fullfile(files(preferredIndex).folder,files(preferredIndex).name);
elseif isscalar(files)
    mniFile = fullfile(files(1).folder,files(1).name);
end

end
