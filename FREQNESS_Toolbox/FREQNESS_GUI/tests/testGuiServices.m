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
verifyEqual(testCase,numel(unique({modules.id})),numel(modules));
verifyGreaterThanOrEqual(testCase,numel(unique({modules.category})),5);
end

function testEverySecondaryModuleHasAValidSchema(testCase)
modules = freqnessgui.moduleRegistry();
context = struct('frequencies',2:2:20,'nComponents',8);
for modulei = 1:numel(modules)
    schema = freqnessgui.secondaryModuleSchema(modules(modulei).id,context);
    fields = [schema.mandatory schema.optional];
    verifyEqual(testCase,numel(unique({fields.key})),numel(fields));
    verifyTrue(testCase,all(~cellfun('isempty',{fields.type})));
end

induced = modules(strcmp({modules.id},'induced_responses'));
removal = modules(strcmp({modules.id},'network_removal'));
verifyTrue(testCase,induced.requiresEvents);
verifyTrue(testCase,removal.requiresSourceData);
end

function testNetworkInspectionReadsManifestContext(testCase)
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,'s')); %#ok<NASGU>
datasetFolder = fullfile(temporaryRoot,'Dataset_Context');
networkFolder = fullfile(temporaryRoot,'FREQ_Networks_Context');
mkdir(datasetFolder);
mkdir(networkFolder);

FREQ = struct('frex',2:2:10,'evecs',zeros(5,4,5), ...
    'ts',zeros(4,20,5),'srate',100); %#ok<NASGU>
save(fullfile(networkFolder,'sub-001_FREQ.mat'),'FREQ');
config = freqnessgui.defaultConfig();
config.network.frequencies = [99 100];
config.network.ncomps = 4;
config.network.samplingRate = 100;
manifest = struct('datasetFolder',datasetFolder,'config',config); %#ok<NASGU>
save(fullfile(networkFolder,'FREQNESS_Manifest.mat'),'manifest');

networkSet = freqnessgui.inspectNetworkFolder(networkFolder);
verifyEqual(testCase,networkSet.frequencies,2:2:10);
verifyEqual(testCase,networkSet.nComponents,4);
verifyEqual(testCase,networkSet.samplingRate,100);
verifyEqual(testCase,networkSet.nVoxels,5);
verifyTrue(testCase,networkSet.sourceDataAvailable);
end

function testFrequencyRangeUsesExactStepGrid(testCase)
[frequencies,snappedRange] = freqnessgui.buildFrequencyVector( ...
    [1.1 4.9],1.2);
verifyEqual(testCase,snappedRange,[1.2 4.8],'AbsTol',1e-12);
verifyEqual(testCase,frequencies,[1.2 2.4 3.6 4.8],'AbsTol',1e-12);
end

function testSecondaryFrequencySelectionUsesImportedGrid(testCase)
available = [1.2 2.4 4.8 9.6 19.2];
[rangeIndices,rangeValues] = freqnessgui.mapFrequencySelection( ...
    available,[2 10]);
verifyEqual(testCase,rangeIndices,[2 4]);
verifyEqual(testCase,rangeValues,[2.4 4.8 9.6]);
[singleIndex,singleValue] = freqnessgui.mapFrequencySelection(available,5);
verifyEqual(testCase,singleIndex,3);
verifyEqual(testCase,singleValue,4.8);
end

function testFWHMPreviewMatchesSupportedForms(testCase)
frequencies = 1.2:1.2:6;
automatic = freqnessgui.computeFWHM(frequencies,[],'logarithmic');
verifySize(testCase,automatic,size(frequencies));
verifyGreaterThan(testCase,automatic,zeros(size(automatic)));
verifyEqual(testCase,freqnessgui.computeFWHM( ...
    frequencies,0.5,'linear'),0.5*ones(size(frequencies)));
verifyEqual(testCase,freqnessgui.computeFWHM( ...
    frequencies,1:numel(frequencies),'linear'),1:numel(frequencies));
end

function testMNICoordinatesLoadAndTranspose(testCase)
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,'s')); %#ok<NASGU>
mniFile = fullfile(temporaryRoot,'custom_mni.mat');
coordinates = reshape(1:12,3,4); %#ok<NASGU>
description = 'ignored metadata'; %#ok<NASGU>
save(mniFile,'coordinates','description');

mni = freqnessgui.loadMNICoordinates(mniFile);
verifyEqual(testCase,mni.coordinates,coordinates.');
verifyEqual(testCase,mni.nPoints,4);
verifyEqual(testCase,mni.variableName,'coordinates');
end

function testDefaultMNIDiscovery(testCase)
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,'s')); %#ok<NASGU>
mniFolder = fullfile(temporaryRoot,'FREQNESS_MNI_Coordinates');
mkdir(mniFolder);
MNI8 = randn(8,3); %#ok<NASGU>
expectedFile = fullfile(mniFolder,'MNI152_8mm_coord_dyi.mat');
save(expectedFile,'MNI8');

verifyEqual(testCase,freqnessgui.findDefaultMNI(temporaryRoot),expectedFile);
end
