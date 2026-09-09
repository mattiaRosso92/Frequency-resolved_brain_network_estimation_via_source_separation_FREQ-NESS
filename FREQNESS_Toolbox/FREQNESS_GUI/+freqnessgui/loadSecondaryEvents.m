function events = loadSecondaryEvents(eventsFile,selectedIds,allIds)
%LOADSECONDARYEVENTS Load and align event samples with selected participants.

if ~(ischar(eventsFile) || (isstring(eventsFile) && isscalar(eventsFile))) || ...
        ~isfile(eventsFile)
    error('FREQNESS:GUI:InvalidEventsFile', ...
        'Select an existing events MAT-file.');
end
selectedIds = cellstr(selectedIds);
allIds = cellstr(allIds);
loaded = load(char(eventsFile));
variableNames = fieldnames(loaded);
preferredNames = {'events','eventOnsets','event_onsets'};
selectedVariable = '';
for namei = 1:numel(preferredNames)
    if isfield(loaded,preferredNames{namei})
        selectedVariable = preferredNames{namei};
        break
    end
end
if isempty(selectedVariable)
    if numel(variableNames) ~= 1
        error('FREQNESS:GUI:InvalidEventsFile', ...
            ['Events files with multiple variables must contain one named ' ...
            'events or eventOnsets.']);
    end
    selectedVariable = variableNames{1};
end
rawEvents = loaded.(selectedVariable);

if isnumeric(rawEvents)
    if numel(selectedIds) ~= 1
        error('FREQNESS:GUI:InvalidEventsFile', ...
            'A numeric events vector can be used only with one participant.');
    end
    events = {rawEvents};
elseif iscell(rawEvents)
    rawEvents = rawEvents(:);
    if numel(rawEvents) == numel(allIds)
        [found,indices] = ismember(selectedIds,allIds);
        if ~all(found)
            error('FREQNESS:GUI:InvalidEventsFile', ...
                'Selected participant IDs could not be aligned with events.');
        end
        events = rawEvents(indices);
    elseif numel(rawEvents) == numel(selectedIds)
        events = rawEvents;
    else
        error('FREQNESS:GUI:InvalidEventsFile', ...
            ['The events cell array must contain one entry per selected ' ...
            'participant or per participant in FREQ_Networks*.']);
    end
elseif isstruct(rawEvents) && isscalar(rawEvents)
    events = cell(numel(selectedIds),1);
    for participanti = 1:numel(selectedIds)
        participantId = selectedIds{participanti};
        validFieldName = matlab.lang.makeValidName(participantId);
        if isvarname(participantId) && isfield(rawEvents,participantId)
            events{participanti} = rawEvents.(participantId);
        elseif isfield(rawEvents,validFieldName)
            events{participanti} = rawEvents.(validFieldName);
        else
            error('FREQNESS:GUI:InvalidEventsFile', ...
                'No events were found for participant %s.',participantId);
        end
    end
else
    error('FREQNESS:GUI:InvalidEventsFile', ...
        'Events must be a numeric vector, cell array, or participant structure.');
end

events = events(:);
for participanti = 1:numel(events)
    value = events{participanti};
    if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
            isempty(value) || any(~isfinite(value)) || ...
            any(value ~= round(value)) || any(value < 1)
        error('FREQNESS:GUI:InvalidEventsFile', ...
            ['Events for participant %s must be a nonempty vector of ' ...
            'positive finite sample indices.'],selectedIds{participanti});
    end
    events{participanti} = unique(value(:)','stable');
end

end
