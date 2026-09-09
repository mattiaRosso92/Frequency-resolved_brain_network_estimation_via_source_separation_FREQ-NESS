function manifest = markSecondaryAnalysisCancelled( ...
        configuration,message)
%MARKSECONDARYANALYSISCANCELLED Finalize a stopped background-run manifest.

if nargin < 2 || isempty(message)
    message = 'Cancelled by the user from the FREQ-NESS GUI.';
end
requiredFields = {'moduleId','outputFolder','participants'};
if ~isstruct(configuration) || ...
        ~all(isfield(configuration,requiredFields))
    error('FREQNESS:GUI:InvalidSecondaryConfig', ...
        'The cancelled analysis configuration is incomplete.');
end

modules = freqnessgui.moduleRegistry();
moduleIndex = find(strcmp({modules.id},configuration.moduleId),1);
if isempty(moduleIndex)
    error('FREQNESS:GUI:UnknownSecondaryModule', ...
        'Unknown secondary-analysis module: %s.',configuration.moduleId);
end
module = modules(moduleIndex);
outputFolder = char(configuration.outputFolder);
if ~isfolder(outputFolder)
    [created,folderMessage] = mkdir(outputFolder);
    if ~created
        error('FREQNESS:GUI:OutputCreationFailed', ...
            'Could not create %s: %s',outputFolder,folderMessage);
    end
end
manifestFile = fullfile(outputFolder,'FREQNESS_AnalysisManifest.mat');
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
    manifest = newManifest(configuration,module,outputFolder);
end

if ~isfield(manifest,'participants') || isempty(manifest.participants)
    manifest.participants = participantRecords( ...
        configuration.participants,module,outputFolder);
end
for participanti = 1:numel(manifest.participants)
    status = manifest.participants(participanti).status;
    if ~strcmp(status,'completed')
        manifest.participants(participanti).status = 'cancelled';
        manifest.participants(participanti).message = message;
    end
end
manifest.status = 'cancelled';
manifest.message = char(message);
manifest.updatedAt = timestamp();
saveAtomic(manifestFile,struct('manifest',manifest));

end

function manifest = newManifest(configuration,module,outputFolder)
manifest = struct();
manifest.schemaVersion = 1;
manifest.kind = 'FREQNESS secondary-analysis results';
manifest.moduleId = module.id;
manifest.moduleName = module.name;
manifest.functionName = module.functionName;
manifest.networkFolder = configuration.networkFolder;
manifest.outputFolder = outputFolder;
manifest.groupOutputFile = fullfile(outputFolder, ...
    ['Group_' module.folderName '.mat']);
manifest.createdAt = timestamp();
manifest.updatedAt = manifest.createdAt;
manifest.status = 'cancelled';
manifest.message = '';
manifest.configuration = configuration;
manifest.participants = participantRecords( ...
    configuration.participants,module,outputFolder);
manifest.figureFiles = repmat( ...
    struct('name','','figFile','','pngFile',''),0,1);
end

function records = participantRecords(participantIds,module,outputFolder)
participantIds = cellstr(participantIds);
template = struct('id','','inputFile','','outputFile','', ...
    'status','cancelled','message','');
records = repmat(template,numel(participantIds),1);
for participanti = 1:numel(participantIds)
    records(participanti).id = participantIds{participanti};
    records(participanti).outputFile = fullfile(outputFolder, ...
        sprintf('%s_%s.mat',participantIds{participanti},module.folderName));
end
end

function saveAtomic(outputFile,payload)
temporaryFile = [tempname(fileparts(outputFile)) '.mat'];
temporaryCleanup = onCleanup(@()deleteIfPresent(temporaryFile));
save(temporaryFile,'-struct','payload','-v7.3');
[moved,message] = movefile(temporaryFile,outputFile,'f');
if ~moved
    error('FREQNESS:GUI:OutputMoveFailed', ...
        'Could not finalize %s: %s',outputFile,message);
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
