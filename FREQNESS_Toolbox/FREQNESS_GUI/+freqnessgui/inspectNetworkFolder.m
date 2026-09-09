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
networkSet.frequencies = [];
networkSet.nComponents = [];
networkSet.samplingRate = [];
networkSet.nVoxels = [];
networkSet.datasetFolder = '';
networkSet.mniFile = '';

if networkSet.hasManifest
    loadedManifest = load(manifestFile,'manifest');
    if isfield(loadedManifest,'manifest') && isstruct(loadedManifest.manifest)
        manifest = loadedManifest.manifest;
        if isfield(manifest,'datasetFolder')
            networkSet.datasetFolder = char(manifest.datasetFolder);
        end
        if isfield(manifest,'config') && isstruct(manifest.config)
            config = manifest.config;
            if isfield(config,'mniFile')
                networkSet.mniFile = char(config.mniFile);
            end
            if isfield(config,'network') && isstruct(config.network)
                network = config.network;
                if isfield(network,'frequencies')
                    networkSet.frequencies = network.frequencies(:)';
                end
                if isfield(network,'ncomps')
                    networkSet.nComponents = network.ncomps;
                end
                if isfield(network,'samplingRate')
                    networkSet.samplingRate = network.samplingRate;
                end
            end
        end
    end
end

firstResultFile = fullfile(resultFiles(1).folder,resultFiles(1).name);
firstResultVariables = whos('-file',firstResultFile);
if any(strcmp({firstResultVariables.name},'participant'))
    participantMetadata = load(firstResultFile,'participant');
    if isstruct(participantMetadata.participant) && ...
            isfield(participantMetadata.participant,'sourceSize') && ...
            ~isempty(participantMetadata.participant.sourceSize)
        networkSet.nVoxels = participantMetadata.participant.sourceSize(1);
    end
end

firstResult = load(firstResultFile,'FREQ');
FREQ = firstResult.FREQ;
if ~isstruct(FREQ) || ~isfield(FREQ,'frex') || isempty(FREQ.frex)
    error('FREQNESS:GUI:InvalidNetworkResult', ...
        'Participant FREQ results must contain a nonempty FREQ.frex vector.');
end
networkSet.frequencies = FREQ.frex(:)';
if isempty(networkSet.nComponents)
    if isfield(FREQ,'evecs') && ~isempty(FREQ.evecs)
        networkSet.nComponents = size(FREQ.evecs,2);
    elseif isfield(FREQ,'ts') && ~isempty(FREQ.ts)
        networkSet.nComponents = size(FREQ.ts,1);
    end
end
if isempty(networkSet.samplingRate) && isfield(FREQ,'srate')
    networkSet.samplingRate = FREQ.srate;
end
if isempty(networkSet.nVoxels) && isfield(FREQ,'evecs') && ...
        ~isempty(FREQ.evecs)
    networkSet.nVoxels = size(FREQ.evecs,1);
end

networkSet.sourceDataAvailable = ~isempty(networkSet.datasetFolder) && ...
    isfolder(networkSet.datasetFolder);

end
