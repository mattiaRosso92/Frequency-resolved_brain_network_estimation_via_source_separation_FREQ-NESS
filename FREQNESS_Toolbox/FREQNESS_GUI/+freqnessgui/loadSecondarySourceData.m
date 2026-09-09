function sourceData = loadSecondarySourceData( ...
        networkSet,participantMetadata,FREQ)
%LOADSECONDARYSOURCEDATA Load analyzed broadband data for NetworkRemoval.

if ~isstruct(networkSet) || ~isfield(networkSet,'datasetFolder') || ...
        ~isfolder(networkSet.datasetFolder)
    error('FREQNESS:GUI:MissingSourceData', ...
        'The original Dataset* folder is unavailable.');
end
if ~iscell(participantMetadata) || isempty(participantMetadata)
    error('FREQNESS:GUI:MissingSourceData', ...
        'Participant source metadata are unavailable.');
end
if ~isstruct(FREQ) || ~isfield(FREQ,'ts') || isempty(FREQ.ts)
    error('FREQNESS:GUI:InvalidParticipantFREQ', ...
        'FREQ.ts is required to align the source data.');
end

nParticipants = numel(participantMetadata);
nVoxels = size(FREQ.evecs,1);
nTime = size(FREQ.ts,2);
sourceData = zeros(nVoxels,nTime,nParticipants);
for participanti = 1:nParticipants
    participant = participantMetadata{participanti};
    if ~isstruct(participant) || ~isfield(participant,'id')
        error('FREQNESS:GUI:MissingSourceData', ...
            'Participant source metadata are incomplete.');
    end
    sourceFile = '';
    if isfield(participant,'sourceFile') && isfile(participant.sourceFile)
        sourceFile = participant.sourceFile;
    else
        candidate = fullfile(networkSet.datasetFolder, ...
            [char(participant.id) '.mat']);
        if isfile(candidate)
            sourceFile = candidate;
        end
    end
    if isempty(sourceFile)
        error('FREQNESS:GUI:MissingSourceData', ...
            'The source MAT-file for participant %s was not found.', ...
            char(participant.id));
    end

    loaded = load(sourceFile);
    variableNames = fieldnames(loaded);
    sourceVariable = '';
    if isfield(participant,'sourceVariable') && ...
            isfield(loaded,participant.sourceVariable)
        sourceVariable = participant.sourceVariable;
    elseif isscalar(variableNames)
        sourceVariable = variableNames{1};
    end
    if isempty(sourceVariable)
        error('FREQNESS:GUI:MissingSourceData', ...
            'The source variable for participant %s is ambiguous.', ...
            char(participant.id));
    end
    value = loaded.(sourceVariable);
    if ~isnumeric(value) || ~isreal(value) || ~ismatrix(value) || ...
            any(~isfinite(value(:))) || size(value,1) ~= nVoxels || ...
            size(value,2) < nTime
        error('FREQNESS:GUI:IncompatibleSourceData', ...
            ['Source data for participant %s must contain at least the ' ...
            'same voxels and analyzed samples as FREQ.ts.'],char(participant.id));
    end
    sourceData(:,:,participanti) = value(:,1:nTime);
end

if nParticipants == 1
    sourceData = sourceData(:,:,1);
end

end
