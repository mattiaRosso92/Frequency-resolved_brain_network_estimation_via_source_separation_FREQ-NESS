function groupFREQ = combineParticipantFREQ(participantFREQ)
%COMBINEPARTICIPANTFREQ Stack participant FREQ structures for group functions.

if isstruct(participantFREQ) && isscalar(participantFREQ)
    participantFREQ = {participantFREQ};
end
if ~iscell(participantFREQ) || isempty(participantFREQ) || ...
        any(~cellfun(@(value)isstruct(value) && isscalar(value),participantFREQ))
    error('FREQNESS:GUI:InvalidParticipantFREQ', ...
        'Participant FREQ inputs must be a nonempty cell array of structures.');
end

participantFREQ = participantFREQ(:)';
nParticipants = numel(participantFREQ);
groupFREQ = participantFREQ{1};

arrayFields = {
    'evals',3
    'evecs',4
    'pats',4
    'ts',4
    };
for fieldi = 1:size(arrayFields,1)
    fieldName = arrayFields{fieldi,1};
    participantDimension = arrayFields{fieldi,2};
    present = cellfun(@(value)isfield(value,fieldName) && ...
        ~isempty(value.(fieldName)),participantFREQ);
    if any(present) && ~all(present)
        error('FREQNESS:GUI:IncompatibleParticipantFREQ', ...
            'FREQ.%s is missing for one or more selected participants.',fieldName);
    end
    if ~all(present)
        continue
    end

    participantArrays = cell(1,nParticipants);
    referenceSize = [];
    for participanti = 1:nParticipants
        value = participantFREQ{participanti}.(fieldName);
        if ~isnumeric(value) || ~isreal(value) || any(~isfinite(value(:)))
            error('FREQNESS:GUI:InvalidParticipantFREQ', ...
                'FREQ.%s must be a real finite numeric array.',fieldName);
        end
        valueSize = size(value);
        if numel(valueSize) > participantDimension && ...
                any(valueSize(participantDimension+1:end) ~= 1)
            error('FREQNESS:GUI:IncompatibleParticipantFREQ', ...
                'FREQ.%s contains unexpected dimensions.',fieldName);
        end
        paddedSize = ones(1,participantDimension);
        paddedSize(1:min(numel(valueSize),participantDimension)) = ...
            valueSize(1:min(numel(valueSize),participantDimension));
        if paddedSize(participantDimension) ~= 1
            error('FREQNESS:GUI:IncompatibleParticipantFREQ', ...
                ['Each participant file must contain one participant in ' ...
                'FREQ.%s.'],fieldName);
        end
        nonParticipantSize = paddedSize;
        nonParticipantSize(participantDimension) = 1;
        if isempty(referenceSize)
            referenceSize = nonParticipantSize;
        elseif ~isequal(nonParticipantSize,referenceSize)
            error('FREQNESS:GUI:IncompatibleParticipantFREQ', ...
                'Selected participants have incompatible FREQ.%s sizes.',fieldName);
        end
        participantArrays{participanti} = reshape(value,paddedSize);
    end
    groupFREQ.(fieldName) = cat(participantDimension,participantArrays{:});
end

sharedFields = {'frex','fwhm','srate','duration','bad_segments', ...
    'regularisation'};
for fieldi = 1:numel(sharedFields)
    fieldName = sharedFields{fieldi};
    present = cellfun(@(value)isfield(value,fieldName),participantFREQ);
    if any(present) && ~all(present)
        error('FREQNESS:GUI:IncompatibleParticipantFREQ', ...
            'FREQ.%s is missing for one or more selected participants.',fieldName);
    end
    if ~all(present)
        continue
    end
    referenceValue = participantFREQ{1}.(fieldName);
    for participanti = 2:nParticipants
        comparisonValue = participantFREQ{participanti}.(fieldName);
        if ismember(fieldName,{'frex','fwhm'})
            isCompatible = isequaln(referenceValue(:),comparisonValue(:));
        else
            isCompatible = isequaln(referenceValue,comparisonValue);
        end
        if ~isCompatible
            error('FREQNESS:GUI:IncompatibleParticipantFREQ', ...
                'Selected participants have different FREQ.%s values.',fieldName);
        end
    end
end

scaleFactors = ones(1,nParticipants);
hasScaleFactors = false;
for participanti = 1:nParticipants
    thisFREQ = participantFREQ{participanti};
    if isfield(thisFREQ,'scale_factors') && ~isempty(thisFREQ.scale_factors)
        value = thisFREQ.scale_factors;
        if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
                ~isfinite(value) || value <= 0
            error('FREQNESS:GUI:InvalidParticipantFREQ', ...
                'Each participant FREQ.scale_factors value must be positive.');
        end
        scaleFactors(participanti) = value;
        hasScaleFactors = true;
    end
end
if hasScaleFactors
    groupFREQ.scale_factors = scaleFactors;
elseif isfield(groupFREQ,'scale_factors')
    groupFREQ = rmfield(groupFREQ,'scale_factors');
end

end
