function manifest = markNetworkEstimationCancelled( ...
        dataset,configuration,message)
%MARKNETWORKESTIMATIONCANCELLED Finalize a stopped core-analysis manifest.

if nargin < 3 || isempty(message)
    message = 'Network estimation cancelled by the user.';
end
outputFolder = configuration.outputFolder;
if isempty(outputFolder)
    outputFolder = dataset.networkFolder;
end
if ~isfolder(outputFolder)
    [created,folderMessage] = mkdir(outputFolder);
    if ~created
        error('FREQNESS:GUI:OutputCreationFailed', ...
            'Could not create %s: %s',outputFolder,folderMessage);
    end
end
manifestFile = fullfile(outputFolder,'FREQNESS_Manifest.mat');
if isfile(manifestFile)
    loaded = load(manifestFile,'manifest');
else
    loaded = struct();
end
if isfield(loaded,'manifest') && isstruct(loaded.manifest)
    manifest = loaded.manifest;
    if isfield(manifest,'status') && strcmp(manifest.status,'completed')
        return
    end
else
    manifest = newManifest(dataset,configuration,outputFolder);
end

if ~isfield(manifest,'participants') || isempty(manifest.participants)
    manifest.participants = participantRecords(dataset,outputFolder);
end
terminalStatuses = {'completed','skipped','failed'};
for participanti = 1:numel(manifest.participants)
    status = manifest.participants(participanti).status;
    if ~any(strcmp(status,terminalStatuses))
        manifest.participants(participanti).status = 'cancelled';
        manifest.participants(participanti).message = char(message);
    end
end
manifest.status = 'cancelled';
manifest.message = char(message);
manifest.updatedAt = timestamp();
saveAtomic(manifestFile,manifest);

end

function manifest = newManifest(dataset,configuration,outputFolder)
manifest = struct();
manifest.schemaVersion = 1;
manifest.kind = 'FREQNESS participant network results';
manifest.datasetName = dataset.name;
manifest.datasetFolder = dataset.folder;
manifest.outputFolder = outputFolder;
manifest.createdAt = timestamp();
manifest.updatedAt = manifest.createdAt;
manifest.status = 'cancelled';
manifest.message = '';
manifest.config = configuration;
manifest.participants = participantRecords(dataset,outputFolder);
end

function records = participantRecords(dataset,outputFolder)
template = struct('id','','sourceFile','','outputFile','', ...
    'status','cancelled','message','');
records = repmat(template,numel(dataset.participants),1);
for participanti = 1:numel(dataset.participants)
    participant = dataset.participants(participanti);
    records(participanti).id = participant.id;
    records(participanti).sourceFile = participant.path;
    records(participanti).outputFile = fullfile( ...
        outputFolder,[participant.id '_FREQ.mat']);
end
end

function saveAtomic(manifestFile,manifest)
temporaryFile = [tempname(fileparts(manifestFile)) '.mat'];
temporaryCleanup = onCleanup(@()deleteIfPresent(temporaryFile));
save(temporaryFile,'manifest','-v7');
[moved,message] = movefile(temporaryFile,manifestFile,'f');
if ~moved
    error('FREQNESS:GUI:ManifestMoveFailed', ...
        'Could not update %s: %s',manifestFile,message);
end
clear temporaryCleanup
end

function deleteIfPresent(filePath)
if isfile(filePath)
    delete(filePath);
end
end

function value = timestamp()
value = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z'));
end
