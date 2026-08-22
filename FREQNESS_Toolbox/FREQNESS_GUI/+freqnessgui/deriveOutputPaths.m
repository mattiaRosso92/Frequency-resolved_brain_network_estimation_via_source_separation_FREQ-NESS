function paths = deriveOutputPaths(inputFolder)
%DERIVEOUTPUTPATHS Derive core and secondary output folders from Dataset*.
%
% Dataset_1 becomes the sibling folders FREQ_Networks_1 and
% FREQ_Analyses_1. The text following the Dataset prefix is preserved.

if ~(ischar(inputFolder) || (isstring(inputFolder) && isscalar(inputFolder)))
    error('FREQNESS:GUI:InvalidFolder','Input folder must be one text scalar.');
end

inputFolder = char(inputFolder);
while numel(inputFolder) > 1 && any(inputFolder(end) == ['/' '\'])
    inputFolder(end) = [];
end

if ~isfolder(inputFolder)
    error('FREQNESS:GUI:InvalidFolder','Folder does not exist: %s',inputFolder);
end

[parentFolder,datasetName] = fileparts(inputFolder);
prefix = 'Dataset';
if numel(datasetName) < numel(prefix) || ...
        ~strcmpi(datasetName(1:numel(prefix)),prefix)
    error('FREQNESS:GUI:InvalidDatasetName', ...
        'The imported folder must be named Dataset* (received "%s").', ...
        datasetName);
end

suffix = datasetName(numel(prefix)+1:end);
paths = struct();
paths.datasetName = datasetName;
paths.suffix = suffix;
paths.parentFolder = parentFolder;
paths.networkFolder = fullfile(parentFolder,['FREQ_Networks' suffix]);
paths.analysisFolder = fullfile(parentFolder,['FREQ_Analyses' suffix]);

end
