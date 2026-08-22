function networkSet = inspectNetworkFolder(networkFolder)
%INSPECTNETWORKFOLDER Validate a FREQ_Networks* folder for secondary use.

if ~(ischar(networkFolder) || ...
        (isstring(networkFolder) && isscalar(networkFolder))) || ...
        ~isfolder(networkFolder)
    error('FREQNESS:GUI:InvalidNetworkFolder', ...
        'Select an existing FREQ_Networks* folder.');
end

networkFolder = char(networkFolder);
while numel(networkFolder) > 1 && any(networkFolder(end) == ['/' '\'])
    networkFolder(end) = [];
end
[parentFolder,folderName] = fileparts(networkFolder);
prefix = 'FREQ_Networks';
if numel(folderName) < numel(prefix) || ...
        ~strcmpi(folderName(1:numel(prefix)),prefix)
    error('FREQNESS:GUI:InvalidNetworkFolder', ...
        'The folder must be named FREQ_Networks* (received "%s").',folderName);
end

resultFiles = dir(fullfile(networkFolder,'*_FREQ.mat'));
if ~isempty(resultFiles)
    [~,order] = sort(lower({resultFiles.name}));
    resultFiles = resultFiles(order);
end
if isempty(resultFiles)
    error('FREQNESS:GUI:EmptyNetworkFolder', ...
        'No participant *_FREQ.mat files were found in %s.',networkFolder);
end

valid = false(numel(resultFiles),1);
participantIds = cell(numel(resultFiles),1);
for filei = 1:numel(resultFiles)
    filePath = fullfile(resultFiles(filei).folder,resultFiles(filei).name);
    variables = whos('-file',filePath);
    valid(filei) = any(strcmp({variables.name},'FREQ'));
    participantIds{filei} = regexprep(resultFiles(filei).name,'_FREQ\.mat$','');
end
if any(~valid)
    error('FREQNESS:GUI:InvalidNetworkResult', ...
        'Every participant result file must contain a FREQ structure.');
end

suffix = folderName(numel(prefix)+1:end);
networkSet = struct();
networkSet.schemaVersion = 1;
networkSet.name = folderName;
networkSet.folder = networkFolder;
networkSet.analysisFolder = fullfile(parentFolder,['FREQ_Analyses' suffix]);
networkSet.files = resultFiles;
networkSet.participantIds = participantIds;
networkSet.nParticipants = numel(resultFiles);

manifestFile = fullfile(networkFolder,'FREQNESS_Manifest.mat');
networkSet.manifestFile = manifestFile;
networkSet.hasManifest = isfile(manifestFile);

end
