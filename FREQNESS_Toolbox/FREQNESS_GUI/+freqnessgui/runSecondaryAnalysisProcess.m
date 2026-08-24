function runSecondaryAnalysisProcess(jobFile)
%RUNSECONDARYANALYSISPROCESS Execute a serialized GUI job in another MATLAB.

configureProcessPath();
loaded = load(jobFile,'job');
job = loaded.job;
progressSequence = 0;
try
    report = freqnessgui.runSecondaryAnalysis( ...
        job.networkSet,job.configuration,job.MNI,@publishProgress, ...
        @()isfile(job.cancellationFile));
    workerResult = struct('status',report.status,'report',report, ...
        'error',struct());
catch exception
    errorInfo = struct( ...
        'identifier',exception.identifier, ...
        'message',exception.message, ...
        'report',getReport(exception,'extended','hyperlinks','off'));
    workerResult = struct('status','failed','report',struct(), ...
        'error',errorInfo);
end
saveAtomic(job.resultFile,struct('workerResult',workerResult),'-v7.3');

    function publishProgress(fraction,message)
        progressSequence = progressSequence+1;
        progress = struct( ...
            'sequence',progressSequence, ...
            'fraction',fraction, ...
            'message',char(message), ...
            'updatedAt',char(datetime('now', ...
            'Format','yyyy-MM-dd HH:mm:ss Z')));
        saveAtomic(job.progressFile,struct('progress',progress),'-v7');
    end

end

function configureProcessPath()
packageFolder = fileparts(mfilename('fullpath'));
guiRoot = fileparts(packageFolder);
toolboxRoot = fileparts(guiRoot);
requiredPaths = {
    guiRoot
    fullfile(toolboxRoot,'FREQNESS_Functions')
    fullfile(toolboxRoot,'FREQNESS_ExternalFunctions')
    fullfile(toolboxRoot,'FREQNESS_ExternalFunctions','nifti_tools')
    };
for pathi = 1:numel(requiredPaths)
    if isfolder(requiredPaths{pathi})
        addpath(requiredPaths{pathi});
    end
end
end

function saveAtomic(outputFile,payload,version)
temporaryFile = [tempname(fileparts(outputFile)) '.mat'];
temporaryCleanup = onCleanup(@()deleteIfPresent(temporaryFile));
save(temporaryFile,'-struct','payload',version);
[moved,message] = movefile(temporaryFile,outputFile,'f');
if ~moved
    error('FREQNESS:GUI:OutputMoveFailed', ...
        'Could not finalize %s: %s',outputFile,message);
end
clear temporaryCleanup
end

function deleteIfPresent(filePath)
if isfile(filePath)
    delete(filePath);
end
end
