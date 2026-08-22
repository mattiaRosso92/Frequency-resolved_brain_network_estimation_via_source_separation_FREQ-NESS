function tests = testGuiServices
%TESTGUISERVICES Unit tests for GUI service functions.
tests = functiontests(localfunctions);
end

function testOutputNaming(testCase)
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,'s')); %#ok<NASGU>
datasetFolder = fullfile(temporaryRoot,'Dataset_Visual');
mkdir(datasetFolder);

paths = freqnessgui.deriveOutputPaths(datasetFolder);
verifyEqual(testCase,paths.networkFolder, ...
    fullfile(temporaryRoot,'FREQ_Networks_Visual'));
verifyEqual(testCase,paths.analysisFolder, ...
    fullfile(temporaryRoot,'FREQ_Analyses_Visual'));
end

function testInvalidDatasetPrefix(testCase)
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,'s')); %#ok<NASGU>
invalidFolder = fullfile(temporaryRoot,'Participants_1');
mkdir(invalidFolder);

verifyError(testCase,@()freqnessgui.deriveOutputPaths(invalidFolder), ...
    'FREQNESS:GUI:InvalidDatasetName');
end

function testDatasetInspection(testCase)
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,'s')); %#ok<NASGU>
datasetFolder = fullfile(temporaryRoot,'Dataset_1');
mkdir(datasetFolder);

participantData = randn(5,100); %#ok<NASGU>
save(fullfile(datasetFolder,'sub-001.mat'),'participantData');
dataset = freqnessgui.inspectDataset(datasetFolder);

verifyTrue(testCase,dataset.isValid);
verifyEqual(testCase,dataset.nParticipants,1);
verifyEqual(testCase,dataset.participants(1).id,'sub-001');
verifyEqual(testCase,dataset.participants(1).dataSize,[5 100]);
end

function testMultipleVariablesAreRejected(testCase)
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,'s')); %#ok<NASGU>
datasetFolder = fullfile(temporaryRoot,'Dataset_1');
mkdir(datasetFolder);

first = randn(5,100); %#ok<NASGU>
second = randn(5,100); %#ok<NASGU>
save(fullfile(datasetFolder,'sub-001.mat'),'first','second');
dataset = freqnessgui.inspectDataset(datasetFolder);

verifyFalse(testCase,dataset.isValid);
verifyFalse(testCase,dataset.participants(1).valid);
end

function testNumericVectorParser(testCase)
verifyEqual(testCase,freqnessgui.parseNumericVector( ...
    '1:0.5:2, 4 6',false,'Test'),[1 1.5 2 4 6]);
verifyEmpty(testCase,freqnessgui.parseNumericVector( ...
    '[]',true,'Test'));
verifyError(testCase,@()freqnessgui.parseNumericVector( ...
    '1:0:2',false,'Test'),'FREQNESS:GUI:InvalidNumericText');
end

function testModuleFoldersAreUnique(testCase)
modules = freqnessgui.moduleRegistry();
verifyEqual(testCase,numel(unique({modules.folderName})),numel(modules));
verifyTrue(testCase,all(~cellfun('isempty',{modules.functionName})));
end
