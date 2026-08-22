function dataset = inspectDataset(datasetFolder)
%INSPECTDATASET Inspect participant MAT files without loading their arrays.

paths = freqnessgui.deriveOutputPaths(datasetFolder);
files = dir(fullfile(datasetFolder,'*.mat'));
if ~isempty(files)
    [~,order] = sort(lower({files.name}));
    files = files(order);
end

if isempty(files)
    error('FREQNESS:GUI:EmptyDataset', ...
        'No .mat participant files were found directly inside %s.',datasetFolder);
end

emptyParticipant = struct( ...
    'id','', ...
    'fileName','', ...
    'path','', ...
    'variableName','', ...
    'dataClass','', ...
    'dataSize',[], ...
    'valid',false, ...
    'message','');
participants = repmat(emptyParticipant,numel(files),1);

numericClasses = {'double','single','int8','uint8','int16','uint16', ...
    'int32','uint32','int64','uint64'};

for filei = 1:numel(files)
    filePath = fullfile(files(filei).folder,files(filei).name);
    [~,participantId] = fileparts(files(filei).name);

    participants(filei).id = participantId;
    participants(filei).fileName = files(filei).name;
    participants(filei).path = filePath;

    try
        variables = whos('-file',filePath);
        if numel(variables) ~= 1
            participants(filei).message = ...
                'File must contain exactly one participant data variable.';
            continue
        end

        variable = variables(1);
        participants(filei).variableName = variable.name;
        participants(filei).dataClass = variable.class;
        participants(filei).dataSize = variable.size;

        isTwoDimensional = numel(variable.size) <= 2 || ...
            all(variable.size(3:end) == 1);
        isNumeric = ismember(variable.class,numericClasses);
        isReal = ~isfield(variable,'complex') || ~variable.complex;
        isNonempty = prod(double(variable.size)) > 0;

        if ~isNumeric || ~isReal || ~isTwoDimensional || ~isNonempty
            participants(filei).message = ...
                'Data must be one non-empty real numeric voxels-by-time matrix.';
            continue
        end

        participants(filei).valid = true;
        participants(filei).message = 'Ready';
    catch exception
        participants(filei).message = exception.message;
    end
end

dataset = struct();
dataset.schemaVersion = 1;
dataset.name = paths.datasetName;
dataset.folder = char(datasetFolder);
dataset.networkFolder = paths.networkFolder;
dataset.analysisFolder = paths.analysisFolder;
dataset.participants = participants;
dataset.nParticipants = numel(participants);
dataset.nValid = sum([participants.valid]);
dataset.isValid = dataset.nValid == dataset.nParticipants;

end
