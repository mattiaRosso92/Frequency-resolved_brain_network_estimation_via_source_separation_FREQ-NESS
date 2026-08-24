function [process,control] = launchSecondaryAnalysisProcess( ...
        networkSet,configuration,MNI,cancellationFile)
%LAUNCHSECONDARYANALYSISPROCESS Start one isolated MATLAB analysis process.

if ~freqnessgui.secondaryBackgroundAvailable()
    error('FREQNESS:GUI:BackgroundUnavailable', ...
        'A separate background MATLAB process is unavailable.');
end
outputFolder = char(configuration.outputFolder);
if ~isfolder(outputFolder)
    [created,message] = mkdir(outputFolder);
    if ~created
        error('FREQNESS:GUI:OutputCreationFailed', ...
            'Could not create %s: %s',outputFolder,message);
    end
end

controlFolder = tempname(outputFolder);
[created,message] = mkdir(controlFolder);
if ~created
    error('FREQNESS:GUI:OutputCreationFailed', ...
        'Could not create background control folder: %s',message);
end
control = struct();
control.controlFolder = controlFolder;
control.jobFile = fullfile(controlFolder,'job.mat');
control.progressFile = fullfile(controlFolder,'progress.mat');
control.resultFile = fullfile(controlFolder,'result.mat');
control.logFile = fullfile(outputFolder,'FREQNESS_Background.log');

try
    job = struct();
    job.networkSet = networkSet;
    job.configuration = configuration;
    job.MNI = MNI;
    job.cancellationFile = cancellationFile;
    job.progressFile = control.progressFile;
    job.resultFile = control.resultFile;
    save(control.jobFile,'job','-v7.3');

    guiRoot = fileparts(fileparts(mfilename('fullpath')));
    matlabExecutable = fullfile(matlabroot,'bin','matlab');
    if ispc && ~isfile(matlabExecutable)
        matlabExecutable = [matlabExecutable '.exe'];
    end
    jobLiteral = strrep(control.jobFile,'''','''''');
    batchCommand = sprintf( ...
        'freqnessgui.runSecondaryAnalysisProcess(''%s'')',jobLiteral);
    arguments = {matlabExecutable,'-sd',guiRoot,'-batch',batchCommand};
    if ismac && strcmp(computer('arch'),'maci64')
        arguments = [{'/usr/bin/arch','-x86_64'} arguments];
    end

    argumentList = javaObject('java.util.ArrayList');
    for argumenti = 1:numel(arguments)
        argumentList.add(javaObject('java.lang.String',arguments{argumenti}));
    end
    processBuilder = javaObject('java.lang.ProcessBuilder',argumentList);
    processBuilder.redirectErrorStream(true);
    processBuilder.redirectOutput(javaObject('java.io.File',control.logFile));
    process = processBuilder.start();
catch exception
    removeControlFolder(controlFolder);
    rethrow(exception)
end

end

function removeControlFolder(controlFolder)
if isfolder(controlFolder)
    rmdir(controlFolder,'s');
end
end
