function requestAnalysisCancellation(cancellationFile)
%REQUESTANALYSISCANCELLATION Create a marker visible to an analysis process.

if ~(ischar(cancellationFile) || ...
        (isstring(cancellationFile) && isscalar(cancellationFile))) || ...
        isempty(cancellationFile)
    error('FREQNESS:GUI:InvalidCancellationFile', ...
        'A cancellation-marker path is required.');
end
cancellationFile = char(cancellationFile);
parentFolder = fileparts(cancellationFile);
if ~isfolder(parentFolder)
    error('FREQNESS:GUI:InvalidCancellationFile', ...
        'The cancellation-marker folder does not exist: %s',parentFolder);
end
[fileId,message] = fopen(cancellationFile,'w');
if fileId < 0
    error('FREQNESS:GUI:CancellationWriteFailed', ...
        'Could not request cancellation: %s',message);
end
fileCleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,'Cancellation requested at %s\n', ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z')));
clear fileCleanup

end
