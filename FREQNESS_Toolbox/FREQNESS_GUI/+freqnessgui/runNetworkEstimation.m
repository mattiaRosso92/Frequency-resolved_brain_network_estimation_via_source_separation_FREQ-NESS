function report = runNetworkEstimation( ...
        dataset,config,progressFcn,cancellationFcn)
%RUNNETWORKESTIMATION Estimate and persist one FREQ result per participant.

if nargin < 3 || isempty(progressFcn)
    progressFcn = @(~,~)[];
end
if nargin < 4 || isempty(cancellationFcn)
    cancellationFcn = @()false;
end
if ~isa(progressFcn,'function_handle') || ...
        ~isa(cancellationFcn,'function_handle')
    error('FREQNESS:GUI:InvalidNetworkConfig', ...
        'Progress and cancellation callbacks must be function handles.');
end

freqnessgui.validateNetworkConfig(config,dataset);
outputFolder = config.outputFolder;
if isempty(outputFolder)
    outputFolder = dataset.networkFolder;
end
if ~isfolder(outputFolder)
    [created,message] = mkdir(outputFolder);
    if ~created
        error('FREQNESS:GUI:OutputCreationFailed', ...
            'Could not create %s: %s',outputFolder,message);
    end
end

participants = dataset.participants;
nParticipants = numel(participants);
resultTemplate = struct('id','','sourceFile','','outputFile','', ...
    'status','pending','message','');
results = repmat(resultTemplate,nParticipants,1);
for participanti = 1:nParticipants
    participantInfo = participants(participanti);
    results(participanti).id = participantInfo.id;
    results(participanti).sourceFile = participantInfo.path;
    results(participanti).outputFile = fullfile( ...
        outputFolder,[participantInfo.id '_FREQ.mat']);
end

manifest = struct();
manifest.schemaVersion = 1;
manifest.kind = 'FREQNESS participant network results';
manifest.datasetName = dataset.name;
manifest.datasetFolder = dataset.folder;
manifest.outputFolder = outputFolder;
manifest.createdAt = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z'));
manifest.updatedAt = manifest.createdAt;
manifest.status = 'running';
manifest.message = 'Preparing participant network estimation.';
manifest.config = config;
manifest.participants = results;
saveManifest(outputFolder,manifest);

frequencyText = formatFrequencyList(config.network.frequencies);
progressFcn(0,sprintf( ...
    'Network estimation: %d participant(s), %d frequencies (%s Hz).', ...
    nParticipants,numel(config.network.frequencies),frequencyText));

wasCancelled = false;
for participanti = 1:nParticipants
    if cancellationFcn()
        wasCancelled = true;
        break
    end
    participantInfo = participants(participanti);
    outputFile = results(participanti).outputFile;

    progressFcn((participanti-1)/nParticipants,sprintf( ...
        'Participant #%d/%d: %s.',participanti,nParticipants,participantInfo.id));

    if isfile(outputFile) && ~config.execution.recomputeExisting
        results(participanti).status = 'skipped';
        results(participanti).message = 'Existing participant result retained.';
        manifest.participants = results;
        manifest.updatedAt = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z'));
        saveManifest(outputFolder,manifest);
        progressFcn(participanti/nParticipants,sprintf( ...
            'Participant #%d/%d retained existing output.', ...
            participanti,nParticipants));
        continue
    end

    try
        checkCancellation(cancellationFcn);
        loaded = load(participantInfo.path);
        variableNames = fieldnames(loaded);
        if numel(variableNames) ~= 1
            error('FREQNESS:GUI:InvalidParticipantFile', ...
                'Participant file must contain exactly one variable.');
        end
        data = loaded.(variableNames{1});
        if ~isnumeric(data) || ~isreal(data) || isempty(data) || ...
                ~ismatrix(data) || any(~isfinite(data(:)))
            error('FREQNESS:GUI:InvalidParticipantData', ...
                'Participant data must be one non-empty real finite numeric matrix.');
        end

        network = config.network;
        progressFcn((participanti-1)/nParticipants,sprintf( ...
            'Participant #%d/%d — analysing frequencies: %s Hz.', ...
            participanti,nParticipants,frequencyText));
        checkCancellation(cancellationFcn);
        FREQ = FREQNESS_NetworkEstimation(data,network.frequencies, ...
            network.samplingRate, ...
            'duration',network.duration, ...
            'fwidth',network.fwidth, ...
            'filter',network.filter, ...
            'regularisation',network.regularisation, ...
            'ncomps',network.ncomps, ...
            'bad_segments',network.badSegments, ...
            'rescale',network.rescale);
        checkCancellation(cancellationFcn);

        participant = struct();
        participant.schemaVersion = 1;
        participant.id = participantInfo.id;
        participant.sourceFile = participantInfo.path;
        participant.sourceVariable = variableNames{1};
        participant.sourceSize = size(data);
        participant.completedAt = char(datetime( ...
            'now','Format','yyyy-MM-dd HH:mm:ss Z'));
        configuration = config;

        temporaryFile = [tempname(outputFolder) '.mat'];
        temporaryCleanup = onCleanup(@()deleteIfPresent(temporaryFile));
        save(temporaryFile,'FREQ','participant','configuration','-v7.3');
        [moved,message] = movefile(temporaryFile,outputFile,'f');
        if ~moved
            error('FREQNESS:GUI:OutputMoveFailed', ...
                'Could not finalize %s: %s',outputFile,message);
        end
        clear temporaryCleanup

        results(participanti).status = 'completed';
        results(participanti).message = 'Network estimation completed.';
        progressFcn(participanti/nParticipants,sprintf( ...
            'Participant #%d/%d completed.',participanti,nParticipants));
    catch exception
        if strcmp(exception.identifier,'FREQNESS:GUI:AnalysisCancelled')
            wasCancelled = true;
        else
            results(participanti).status = 'failed';
            results(participanti).message = exception.message;
            progressFcn(participanti/nParticipants,sprintf( ...
                'Participant #%d/%d failed: %s', ...
                participanti,nParticipants,exception.message));
        end
    end

    manifest.participants = results;
    manifest.updatedAt = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z'));
    saveManifest(outputFolder,manifest);
    if wasCancelled
        break
    end
end

if wasCancelled
    cancellationMessage = 'Network estimation cancelled by the user.';
    for participanti = 1:nParticipants
        if strcmp(results(participanti).status,'pending')
            results(participanti).status = 'cancelled';
            results(participanti).message = cancellationMessage;
        end
    end
    manifest.status = 'cancelled';
    manifest.message = cancellationMessage;
    progressFcn(1,cancellationMessage);
else
    manifest.status = 'completed';
    manifest.message = 'Network estimation finished.';
    progressFcn(1,manifest.message);
end
manifest.participants = results;
manifest.updatedAt = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z'));
saveManifest(outputFolder,manifest);

statuses = {results.status};
report = struct();
report.outputFolder = outputFolder;
report.manifestFile = fullfile(outputFolder,'FREQNESS_Manifest.mat');
report.participants = results;
report.nCompleted = sum(strcmp(statuses,'completed'));
report.nSkipped = sum(strcmp(statuses,'skipped'));
report.nFailed = sum(strcmp(statuses,'failed'));
report.nCancelled = sum(strcmp(statuses,'cancelled'));
report.status = manifest.status;

end

function checkCancellation(cancellationFcn)
if cancellationFcn()
    error('FREQNESS:GUI:AnalysisCancelled', ...
        'Network estimation was cancelled by the user.');
end
end

function saveManifest(outputFolder,manifest)
manifestFile = fullfile(outputFolder,'FREQNESS_Manifest.mat');
temporaryFile = [tempname(outputFolder) '.mat'];
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

function text = formatFrequencyList(frequencies)
frequencies = frequencies(:)';
if numel(frequencies) <= 10
    pieces = arrayfun(@(value)sprintf('%.4g',value),frequencies, ...
        'UniformOutput',false);
    text = strjoin(pieces,', ');
else
    firstPieces = arrayfun(@(value)sprintf('%.4g',value),frequencies(1:4), ...
        'UniformOutput',false);
    lastPieces = arrayfun(@(value)sprintf('%.4g',value),frequencies(end-2:end), ...
        'UniformOutput',false);
    text = sprintf('%s, ..., %s',strjoin(firstPieces,', '), ...
        strjoin(lastPieces,', '));
end
end
