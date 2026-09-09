function report = runSecondaryAnalysis( ...
        networkSet,configuration,MNI,progressFcn,cancellationFcn)
%RUNSECONDARYANALYSIS Execute and persist one configured secondary analysis.

if nargin < 3
    MNI = [];
end
if nargin < 4 || isempty(progressFcn)
    progressFcn = @(~,~)[];
end
if nargin < 5 || isempty(cancellationFcn)
    cancellationFcn = @()false;
end
if ~isa(progressFcn,'function_handle') || ...
        ~isa(cancellationFcn,'function_handle')
    error('FREQNESS:GUI:InvalidSecondaryConfig', ...
        'Progress and cancellation callbacks must be function handles.');
end
if ~isstruct(networkSet) || ~isfield(networkSet,'folder') || ...
        ~isfolder(networkSet.folder)
    error('FREQNESS:GUI:MissingNetworkResults', ...
        'Import a valid FREQ_Networks* folder before execution.');
end
if ~isstruct(configuration) || ~isfield(configuration,'moduleId') || ...
        ~isfield(configuration,'participants') || ...
        ~isfield(configuration,'outputFolder')
    error('FREQNESS:GUI:InvalidSecondaryConfig', ...
        'The secondary-analysis configuration is incomplete.');
end

modules = freqnessgui.moduleRegistry();
moduleIndex = find(strcmp({modules.id},configuration.moduleId),1);
if isempty(moduleIndex)
    error('FREQNESS:GUI:UnknownSecondaryModule', ...
        'Unknown secondary-analysis module: %s.',configuration.moduleId);
end
module = modules(moduleIndex);
outputFolder = char(configuration.outputFolder);
ensureFolder(outputFolder);
figureFolder = fullfile(outputFolder,'Figures');
ensureFolder(figureFolder);

selectedIds = cellstr(configuration.participants);
selectedIds = selectedIds(:)';
if isempty(selectedIds) || numel(unique(selectedIds)) ~= numel(selectedIds)
    error('FREQNESS:GUI:InvalidSecondaryConfig', ...
        'Select one or more unique participant IDs.');
end
nParticipants = numel(selectedIds);

resultTemplate = struct('id','','inputFile','','outputFile','', ...
    'status','pending','message','');
participantResults = repmat(resultTemplate,nParticipants,1);
for participanti = 1:nParticipants
    participantResults(participanti).id = selectedIds{participanti};
    participantResults(participanti).outputFile = fullfile(outputFolder, ...
        sprintf('%s_%s.mat',selectedIds{participanti},module.folderName));
end

manifest = struct();
manifest.schemaVersion = 1;
manifest.kind = 'FREQNESS secondary-analysis results';
manifest.moduleId = module.id;
manifest.moduleName = module.name;
manifest.functionName = module.functionName;
manifest.networkFolder = networkSet.folder;
manifest.outputFolder = outputFolder;
manifest.groupOutputFile = fullfile(outputFolder, ...
    ['Group_' module.folderName '.mat']);
manifest.createdAt = timestamp();
manifest.updatedAt = manifest.createdAt;
manifest.status = 'running';
manifest.message = 'Preparing secondary analysis.';
manifest.configuration = configuration;
manifest.participants = participantResults;
manifest.figureFiles = emptyFigureRecords();
manifestFile = fullfile(outputFolder,'FREQNESS_AnalysisManifest.mat');
saveManifest(manifestFile,manifest);

progressFcn(0,sprintf('%s: preparing %d participant(s).', ...
    module.name,nParticipants));
figuresBefore = findall(groot,'Type','figure');
figureCleanup = onCleanup(@()closeNewFigures(figuresBefore));

try
    checkCancellation(cancellationFcn,module.name);
    [participantFREQ,participantMetadata,inputFiles] = ...
        loadParticipantResults(networkSet,selectedIds,progressFcn, ...
        cancellationFcn,module.name);
    for participanti = 1:nParticipants
        participantResults(participanti).inputFile = inputFiles{participanti};
    end
    manifest.participants = participantResults;
    manifest.updatedAt = timestamp();
    manifest.message = 'Participant FREQ results loaded.';
    saveManifest(manifestFile,manifest);

    groupFREQ = freqnessgui.combineParticipantFREQ(participantFREQ);
    checkCancellation(cancellationFcn,module.name);
    events = [];
    if strcmp(module.id,'induced_responses')
        events = freqnessgui.loadSecondaryEvents( ...
            configuration.values.eventsFile,selectedIds, ...
            networkSet.participantIds);
    end
    sourceData = [];
    if strcmp(module.id,'network_removal')
        sourceData = freqnessgui.loadSecondarySourceData( ...
            networkSet,participantMetadata,groupFREQ);
    end
    checkCancellation(cancellationFcn,module.name);

    progressFcn(0.22,sprintf( ...
        'Running %s for %d participant(s); frequencies: %s Hz.', ...
        module.functionName,nParticipants, ...
        formatSelectedFrequencies(configuration,networkSet.frequencies)));
    call = freqnessgui.buildSecondaryFunctionCall( ...
        configuration,groupFREQ,MNI,events,sourceData);
    outputValues = cell(1,numel(call.outputNames));
    if isempty(outputValues)
        feval(call.functionName,call.arguments{:});
        groupOutput = struct();
    else
        [outputValues{:}] = feval(call.functionName,call.arguments{:});
        groupOutput = cell2struct(outputValues,call.outputNames,2);
    end
    checkCancellation(cancellationFcn,module.name);

    progressFcn(0.72,'Numerical analysis completed; saving figures.');
    figureFiles = saveNewFigures(figuresBefore,figureFolder);
    checkCancellation(cancellationFcn,module.name);
    completedAt = timestamp();
    groupPayload = struct( ...
        'groupOutput',groupOutput, ...
        'configuration',configuration, ...
        'participants',{participantMetadata}, ...
        'figureFiles',figureFiles, ...
        'completedAt',completedAt);
    saveAtomic(manifest.groupOutputFile,groupPayload);

    for participanti = 1:nParticipants
        checkCancellation(cancellationFcn,module.name);
        participantOutput = freqnessgui.extractParticipantOutput( ...
            module.id,groupOutput,participanti,nParticipants);
        participant = participantMetadata{participanti};
        participant.completedAt = completedAt;
        participant.sourceResultFile = inputFiles{participanti};
        participantPayload = struct( ...
            'participantOutput',participantOutput, ...
            'participant',participant, ...
            'configuration',configuration);
        saveAtomic(participantResults(participanti).outputFile, ...
            participantPayload);
        participantResults(participanti).status = 'completed';
        participantResults(participanti).message = ...
            'Participant output saved from the group analysis.';
        progressFcn(0.78+0.20*participanti/nParticipants,sprintf( ...
            'Participant #%d/%d saved: %s.',participanti,nParticipants, ...
            selectedIds{participanti}));
        checkCancellation(cancellationFcn,module.name);
    end

    manifest.status = 'completed';
    manifest.message = sprintf('%s completed successfully.',module.name);
    manifest.updatedAt = timestamp();
    manifest.participants = participantResults;
    manifest.figureFiles = figureFiles;
    saveManifest(manifestFile,manifest);
    progressFcn(1,manifest.message);
catch exception
    isCancellation = strcmp(exception.identifier, ...
        'FREQNESS:GUI:AnalysisCancelled');
    if isCancellation
        terminalStatus = 'cancelled';
        terminalMessage = sprintf('%s cancelled by the user.',module.name);
    else
        terminalStatus = 'failed';
        terminalMessage = exception.message;
    end
    for participanti = 1:nParticipants
        if strcmp(participantResults(participanti).status,'pending')
            participantResults(participanti).status = terminalStatus;
            participantResults(participanti).message = terminalMessage;
        end
    end
    manifest.status = terminalStatus;
    manifest.message = terminalMessage;
    manifest.updatedAt = timestamp();
    manifest.participants = participantResults;
    try
        saveManifest(manifestFile,manifest);
    catch
        % Preserve the original analysis exception.
    end
    progressFcn(1,terminalMessage);
    if ~isCancellation
        rethrow(exception)
    end
end

clear figureCleanup
statuses = {participantResults.status};
report = struct();
report.moduleId = module.id;
report.outputFolder = outputFolder;
report.groupOutputFile = manifest.groupOutputFile;
report.manifestFile = manifestFile;
report.figureFiles = manifest.figureFiles;
report.participants = participantResults;
report.nCompleted = sum(strcmp(statuses,'completed'));
report.nFailed = sum(strcmp(statuses,'failed'));
report.nCancelled = sum(strcmp(statuses,'cancelled'));
report.status = manifest.status;

end

function [participantFREQ,metadata,inputFiles] = loadParticipantResults( ...
        networkSet,selectedIds,progressFcn,cancellationFcn,moduleName)
nParticipants = numel(selectedIds);
participantFREQ = cell(1,nParticipants);
metadata = cell(1,nParticipants);
inputFiles = cell(1,nParticipants);
for participanti = 1:nParticipants
    checkCancellation(cancellationFcn,moduleName);
    participantId = selectedIds{participanti};
    resultIndex = find(strcmp(networkSet.participantIds,participantId),1);
    if isempty(resultIndex)
        error('FREQNESS:GUI:MissingParticipantResult', ...
            'No FREQ result was found for participant %s.',participantId);
    end
    resultFile = fullfile(networkSet.files(resultIndex).folder, ...
        networkSet.files(resultIndex).name);
    resultVariables = whos('-file',resultFile);
    if any(strcmp({resultVariables.name},'participant'))
        loaded = load(resultFile,'FREQ','participant');
    else
        loaded = load(resultFile,'FREQ');
    end
    if ~isfield(loaded,'FREQ') || ~isstruct(loaded.FREQ)
        error('FREQNESS:GUI:InvalidNetworkResult', ...
            '%s does not contain a valid FREQ structure.',resultFile);
    end
    participantFREQ{participanti} = loaded.FREQ;
    if isfield(loaded,'participant') && isstruct(loaded.participant)
        participant = loaded.participant;
    else
        participant = struct();
    end
    participant.id = participantId;
    participant.resultFile = resultFile;
    metadata{participanti} = participant;
    inputFiles{participanti} = resultFile;
    progressFcn(0.03+0.15*participanti/nParticipants,sprintf( ...
        'Loaded participant #%d/%d: %s.',participanti,nParticipants, ...
        participantId));
    checkCancellation(cancellationFcn,moduleName);
end
end

function checkCancellation(cancellationFcn,moduleName)
if cancellationFcn()
    error('FREQNESS:GUI:AnalysisCancelled', ...
        '%s was cancelled by the user.',moduleName);
end
end

function figureFiles = saveNewFigures(figuresBefore,figureFolder)
drawnow
allFigures = findall(groot,'Type','figure');
newFigures = allFigures;
for figurei = 1:numel(figuresBefore)
    newFigures(newFigures == figuresBefore(figurei)) = [];
end
newFigures = flipud(newFigures(:));
figureFiles = repmat(struct('name','','figFile','','pngFile',''), ...
    numel(newFigures),1);
for figurei = 1:numel(newFigures)
    figureHandle = newFigures(figurei);
    figureName = char(string(figureHandle.Name));
    if isempty(strtrim(figureName))
        figureName = sprintf('Figure %02d',figurei);
    end
    safeName = regexprep(figureName,'[^A-Za-z0-9_-]+','_');
    safeName = regexprep(safeName,'^_+|_+$','');
    if isempty(safeName)
        safeName = sprintf('Figure_%02d',figurei);
    end
    baseName = sprintf('%02d_%s',figurei,safeName);
    figFile = fullfile(figureFolder,[baseName '.fig']);
    pngFile = fullfile(figureFolder,[baseName '.png']);
    try
        savefig(figureHandle,figFile);
    catch exception
        warning('FREQNESS:GUI:FigureSaveFailed', ...
            'Could not save %s: %s',figFile,exception.message);
        figFile = '';
    end
    try
        exportgraphics(figureHandle,pngFile,'Resolution',150);
    catch exception
        warning('FREQNESS:GUI:FigureExportFailed', ...
            'Could not export %s: %s',pngFile,exception.message);
        pngFile = '';
    end
    figureFiles(figurei).name = figureName;
    figureFiles(figurei).figFile = figFile;
    figureFiles(figurei).pngFile = pngFile;
end
closeNewFigures(figuresBefore)
end

function closeNewFigures(figuresBefore)
allFigures = findall(groot,'Type','figure');
for figurei = 1:numel(allFigures)
    if ~any(allFigures(figurei) == figuresBefore)
        try
            close(allFigures(figurei));
        catch
            % Figure cleanup must not mask analysis results.
        end
    end
end
end

function ensureFolder(folderPath)
if isfolder(folderPath)
    return
end
[created,message] = mkdir(folderPath);
if ~created
    error('FREQNESS:GUI:OutputCreationFailed', ...
        'Could not create %s: %s',folderPath,message);
end
end

function saveManifest(manifestFile,manifest)
saveAtomic(manifestFile,struct('manifest',manifest));
end

function saveAtomic(outputFile,payload)
outputFolder = fileparts(outputFile);
temporaryFile = [tempname(outputFolder) '.mat'];
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

function textValue = formatSelectedFrequencies(configuration,defaultFrequencies)
frequencyValues = [];
fieldNames = fieldnames(configuration.selectedFrequencies);
for fieldi = 1:numel(fieldNames)
    value = configuration.selectedFrequencies.(fieldNames{fieldi});
    frequencyValues = [frequencyValues reshape(value,1,[])]; %#ok<AGROW>
end
if isempty(frequencyValues)
    frequencyValues = defaultFrequencies;
end
frequencyValues = unique(frequencyValues,'stable');
if isscalar(frequencyValues)
    textValue = sprintf('%.4g',frequencyValues);
elseif numel(frequencyValues) <= 6
    pieces = arrayfun(@(value)sprintf('%.4g',value),frequencyValues, ...
        'UniformOutput',false);
    textValue = strjoin(pieces,', ');
else
    textValue = sprintf('%.4g–%.4g (%d selected)', ...
        min(frequencyValues),max(frequencyValues),numel(frequencyValues));
end
end

function records = emptyFigureRecords()
records = repmat(struct('name','','figFile','','pngFile',''),0,1);
end
