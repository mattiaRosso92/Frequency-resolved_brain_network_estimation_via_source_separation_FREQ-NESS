function participantOutput = extractParticipantOutput( ...
        moduleId,groupOutput,participantIndex,nParticipants)
%EXTRACTPARTICIPANTOUTPUT Extract one participant from a group backend result.

if ~isstruct(groupOutput) || ~isscalar(groupOutput)
    error('FREQNESS:GUI:InvalidSecondaryOutput', ...
        'The group secondary-analysis output must be one structure.');
end
if ~isnumeric(participantIndex) || ~isscalar(participantIndex) || ...
        participantIndex ~= round(participantIndex) || participantIndex < 1 || ...
        participantIndex > nParticipants
    error('FREQNESS:GUI:InvalidParticipantIndex', ...
        'Participant index is outside the group output.');
end

moduleId = lower(char(moduleId));
switch moduleId
    case 'entropy'
        participantOutput = struct( ...
            'H2',selectColumn(groupOutput.H2,participantIndex,nParticipants), ...
            'ED',selectColumn(groupOutput.ED,participantIndex,nParticipants));

    case 'exponential_dk'
        participantOutput = struct( ...
            'decayCoeff',groupOutput.decayCoeff(participantIndex), ...
            'goodFit',selectGoodFit(groupOutput.goodFit, ...
            participantIndex,nParticipants));

    case {'freq_gradients','comp_gradients'}
        participantOutput = struct( ...
            'gradCoeff',groupOutput.gradCoeff(:,:,participantIndex), ...
            'goodFit',selectGoodFit(groupOutput.goodFit, ...
            participantIndex,nParticipants));

    case 'cross_coupling'
        participantOutput = struct('CFC',selectCFCParticipant( ...
            groupOutput.CFC,participantIndex));

    case 'induced_responses'
        participantOutput = struct('IND',selectINDParticipant( ...
            groupOutput.IND,participantIndex));

    case 'backprojection'
        participantOutput = struct('backProj',selectThirdDimension( ...
            groupOutput.backProj,participantIndex,nParticipants));

    case 'network_removal'
        participantOutput = struct( ...
            'dataClean',selectThirdDimension(groupOutput.dataClean, ...
            participantIndex,nParticipants), ...
            'removedActivity',selectThirdDimension( ...
            groupOutput.removedActivity,participantIndex,nParticipants));

    case 'visualizer'
        participantOutput = struct( ...
            'includedInGroupVisualization',true, ...
            'participantIndex',participantIndex);

    otherwise
        error('FREQNESS:GUI:UnknownSecondaryModule', ...
            'Unknown secondary-analysis module: %s.',moduleId);
end

end

function value = selectColumn(value,participantIndex,nParticipants)
if nParticipants == 1
    return
end
if size(value,2) ~= nParticipants
    error('FREQNESS:GUI:InvalidSecondaryOutput', ...
        'A participant-indexed output has an unexpected size.');
end
value = value(:,participantIndex);
end

function value = selectThirdDimension(value,participantIndex,nParticipants)
if nParticipants == 1
    return
end
if size(value,3) ~= nParticipants
    error('FREQNESS:GUI:InvalidSecondaryOutput', ...
        'A participant-indexed output has an unexpected size.');
end
value = value(:,:,participantIndex);
end

function goodFit = selectGoodFit(goodFit,participantIndex,nParticipants)
fieldNames = fieldnames(goodFit);
for fieldi = 1:numel(fieldNames)
    fieldName = fieldNames{fieldi};
    value = goodFit.(fieldName);
    if nParticipants == 1
        continue
    end
    if isvector(value) && numel(value) == nParticipants
        goodFit.(fieldName) = value(participantIndex);
    elseif size(value,2) == nParticipants
        goodFit.(fieldName) = value(:,participantIndex);
    else
        error('FREQNESS:GUI:InvalidSecondaryOutput', ...
            'goodFit.%s has an unexpected participant dimension.',fieldName);
    end
end
end

function CFC = selectCFCParticipant(CFC,participantIndex)
threeDimensionalFields = {'PAC_all','fitted_PAC','coefficients'};
for fieldi = 1:numel(threeDimensionalFields)
    fieldName = threeDimensionalFields{fieldi};
    CFC.(fieldName) = CFC.(fieldName)(:,participantIndex,:);
end
columnFields = {'amplitude_raw','amplitude_normalized','preferred_phase', ...
    'dc_offset','mse','r2','valid_bins','sAmpl','pShift','dcOff','mFrex', ...
    'goodFit'};
for fieldi = 1:numel(columnFields)
    fieldName = columnFields{fieldi};
    CFC.(fieldName) = CFC.(fieldName)(:,participantIndex);
end
CFC.lfo_phase = CFC.lfo_phase(:,participantIndex);
CFC.PAC_avg = reshape(CFC.PAC_all(:,1,:), ...
    size(CFC.PAC_all,1),size(CFC.PAC_all,3));
end

function IND = selectINDParticipant(IND,participantIndex)
IND.power = IND.power(:,:,:,participantIndex);
IND.events = IND.events(participantIndex);
IND.ntrials = IND.ntrials(participantIndex);
if isfield(IND,'power_trials')
    IND.power_trials = IND.power_trials(participantIndex);
end
if isfield(IND,'power_scaled')
    IND.power_scaled = IND.power_scaled(:,:,:,participantIndex);
end
end
