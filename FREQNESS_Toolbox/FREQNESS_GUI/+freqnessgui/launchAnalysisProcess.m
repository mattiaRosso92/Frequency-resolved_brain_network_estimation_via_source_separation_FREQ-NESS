function [process,control] = launchAnalysisProcess( ...
        job,outputFolder,workerFunction,logFileName)
%LAUNCHANALYSISPROCESS Start one isolated MATLAB analysis process.

if ~freqnessgui.backgroundProcessAvailable()
    error('FREQNESS:GUI:BackgroundUnavailable', ...
        'A separate background MATLAB process is unavailable.');
end
if nargin < 4 || isempty(logFileName)
    logFileName = 'FREQNESS_Background.log';
end
outputFolder = char(outputFolder);
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
control.logFile = fullfile(outputFolder,logFileName);

try
    job.progressFile = control.progressFile;
    job.resultFile = control.resultFile;
    save(control.jobFile,'job','-v7.3');

    guiRoot = fileparts(fileparts(mfilename('fullpath')));
    matlabExecutable = fullfile(matlabroot,'bin','matlab');
    if ispc && ~isfile(matlabExecutable)
        matlabExecutable = [matlabExecutable '.exe'];
    end
    jobLiteral = strrep(control.jobFile,'''','''''');
    batchCommand = sprintf('%s(''%s'')',workerFunction,jobLiteral);
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
