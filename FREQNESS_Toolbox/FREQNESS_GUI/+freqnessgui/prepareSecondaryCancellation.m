function cancellationFile = prepareSecondaryCancellation(outputFolder)
%PREPARESECONDARYCANCELLATION Reset and return a run cancellation marker.

if ~(ischar(outputFolder) || ...
        (isstring(outputFolder) && isscalar(outputFolder)))
    error('FREQNESS:GUI:InvalidSecondaryConfig', ...
        'The secondary output folder must be a text scalar.');
end
outputFolder = char(outputFolder);
if ~isfolder(outputFolder)
    [created,message] = mkdir(outputFolder);
    if ~created
        error('FREQNESS:GUI:OutputCreationFailed', ...
            'Could not create %s: %s',outputFolder,message);
    end
end
cancellationFile = fullfile(outputFolder, ...
    '.FREQNESS_CancelRequested');
if isfile(cancellationFile)
    delete(cancellationFile);
end

end
