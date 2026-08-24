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

function testFilterResponsesMatchFilterFGxGaussian(testCase)
frequencies = [5 10];
fwhm = [1 2];
[frequencyAxis,responses] = freqnessgui.computeFilterResponses( ...
    frequencies,fwhm,[0 20],4001);

verifySize(testCase,responses,[2 numel(frequencyAxis)]);
verifyEqual(testCase,max(responses,[],2),ones(2,1),'AbsTol',1e-12);
for frequencyi = 1:numel(frequencies)
    [~,peakIndex] = max(responses(frequencyi,:));
    verifyEqual(testCase,frequencyAxis(peakIndex),frequencies(frequencyi), ...
        'AbsTol',0.005);
end

[~,peakIndex] = max(responses(2,:));
[~,lowerRelative] = min(abs(responses(2,1:peakIndex)-0.5));
[~,upperRelative] = min(abs(responses(2,peakIndex:end)-0.5));
upperIndex = peakIndex-1+upperRelative;
empiricalFWHM = frequencyAxis(upperIndex)-frequencyAxis(lowerRelative);
verifyEqual(testCase,empiricalFWHM,fwhm(2),'AbsTol',0.03);
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

function testParticipantFREQCombination(testCase)
first = syntheticFREQ(1);
second = syntheticFREQ(2);
group = freqnessgui.combineParticipantFREQ({first,second});

verifySize(testCase,group.evals,[4 3 2]);
verifySize(testCase,group.evecs,[4 2 3 2]);
verifySize(testCase,group.pats,[4 2 3 2]);
verifySize(testCase,group.ts,[2 20 3 2]);
verifyEqual(testCase,group.evals(:,:,1),first.evals);
verifyEqual(testCase,group.evals(:,:,2),second.evals);
verifyEqual(testCase,group.scale_factors,[1 2]);
end

function testEverySecondaryModuleBuildsBackendCall(testCase)
modules = freqnessgui.moduleRegistry();
context = struct('frequencies',[2 4 8 16],'nComponents',5);
FREQ = syntheticFREQ(1);
MNI = zeros(4,3);
for modulei = 1:numel(modules)
    configuration = defaultSecondaryConfiguration( ...
        modules(modulei).id,context,tempdir,{'sub-001'});
    call = freqnessgui.buildSecondaryFunctionCall( ...
        configuration,FREQ,MNI,{10},zeros(4,20));
    verifyEqual(testCase,call.functionName,modules(modulei).functionName);
    verifyTrue(testCase,iscell(call.arguments));
    verifyTrue(testCase,iscell(call.outputNames));
end

configuration = defaultSecondaryConfiguration( ...
    'comp_gradients',context,tempdir,{'sub-001'});
call = freqnessgui.buildSecondaryFunctionCall( ...
    configuration,FREQ,MNI,[],[]);
rangeIndex = find(strcmp(call.arguments,'comps2model'),1)+1;
verifyEqual(testCase,call.arguments{rangeIndex},[1 3]);

configuration = defaultSecondaryConfiguration( ...
    'visualizer',context,tempdir,{'sub-001'});
call = freqnessgui.buildSecondaryFunctionCall( ...
    configuration,FREQ,MNI,[],[]);
verifyEqual(testCase,call.arguments{2}.ncomps,3);
verifyEqual(testCase,call.arguments{3}.ncomps,1);
end

function testEventsAlignToSelectedParticipants(testCase)
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,'s')); %#ok<NASGU>
eventsFile = fullfile(temporaryRoot,'events.mat');
eventOnsets = {[10 20],[30 40],[50 60]}; %#ok<NASGU>
save(eventsFile,'eventOnsets');

events = freqnessgui.loadSecondaryEvents(eventsFile, ...
    {'sub-003','sub-001'},{'sub-001','sub-002','sub-003'});
verifyEqual(testCase,events,{[50 60];[10 20]});
end

function testParticipantOutputExtraction(testCase)
nParticipants = 2;
entropy = struct('H2',reshape(1:6,3,2), ...
    'ED',reshape(11:16,3,2));
participant = freqnessgui.extractParticipantOutput( ...
    'entropy',entropy,2,nParticipants);
verifyEqual(testCase,participant.H2,entropy.H2(:,2));

gradient = struct('gradCoeff',reshape(1:18,3,3,2), ...
    'goodFit',struct('R2_best',reshape(1:6,3,2), ...
    'bestOrder',ones(3,2)));
participant = freqnessgui.extractParticipantOutput( ...
    'freq_gradients',gradient,2,nParticipants);
verifySize(testCase,participant.gradCoeff,[3 3]);
verifyEqual(testCase,participant.goodFit.R2_best, ...
    gradient.goodFit.R2_best(:,2));

CFC = mockCFC();
participant = freqnessgui.extractParticipantOutput( ...
    'cross_coupling',struct('CFC',CFC),2,nParticipants);
verifySize(testCase,participant.CFC.PAC_all,[3 1 4]);
verifySize(testCase,participant.CFC.PAC_avg,[3 4]);
verifySize(testCase,participant.CFC.lfo_phase,[10 1]);

IND = struct('power',zeros(2,5,3,2),'time',1:5,'frex',[2 4 8], ...
    'fwhm',[0.5 0.8 1.2],'comps',[1 2], ...
    'events',{{[10 20],[12 22]}},'ntrials',[2;2], ...
    'epoch_window',[-0.1 0.3],'baseline_window',[-0.1 0], ...
    'srate',10,'power_trials',{{zeros(2,5,3,2),zeros(2,5,3,2)}}, ...
    'power_scaled',zeros(2,5,3,2));
participant = freqnessgui.extractParticipantOutput( ...
    'induced_responses',struct('IND',IND),2,nParticipants);
verifySize(testCase,participant.IND.power,[2 5 3]);
verifyEqual(testCase,participant.IND.events,{[12 22]});

backprojection = struct('backProj',zeros(4,20,2));
participant = freqnessgui.extractParticipantOutput( ...
    'backprojection',backprojection,2,nParticipants);
verifySize(testCase,participant.backProj,[4 20]);

visualizer = freqnessgui.extractParticipantOutput( ...
    'visualizer',struct(),2,nParticipants);
verifyTrue(testCase,visualizer.includedInGroupVisualization);
end

function testSourceDataAlignmentForNetworkRemoval(testCase)
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,'s')); %#ok<NASGU>
datasetFolder = fullfile(temporaryRoot,'Dataset_Source');
mkdir(datasetFolder);
metadata = cell(1,2);
for participanti = 1:2
    sourceData = participanti*reshape(1:100,4,25); %#ok<NASGU>
    sourceFile = fullfile(datasetFolder,sprintf('sub-%03d.mat',participanti));
    save(sourceFile,'sourceData');
    metadata{participanti} = struct( ...
        'id',sprintf('sub-%03d',participanti), ...
        'sourceFile',sourceFile, ...
        'sourceVariable','sourceData');
end
networkSet = struct('datasetFolder',datasetFolder);
groupFREQ = freqnessgui.combineParticipantFREQ( ...
    {syntheticFREQ(1),syntheticFREQ(2)});
loadedData = freqnessgui.loadSecondarySourceData( ...
    networkSet,metadata,groupFREQ);
verifySize(testCase,loadedData,[4 20 2]);
expectedFirst = reshape(1:100,4,25);
verifyEqual(testCase,loadedData(:,:,1),expectedFirst(:,1:20));
end

function testEntropyExecutionPersistsGroupAndParticipantOutputs(testCase)
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,'s')); %#ok<NASGU>
networkFolder = fullfile(temporaryRoot,'FREQ_Networks_Adapter');
mkdir(networkFolder);
participantIds = {'sub-001','sub-002'};
for participanti = 1:numel(participantIds)
    FREQ = syntheticFREQ(participanti); %#ok<NASGU>
    participant = struct('id',participantIds{participanti}); %#ok<NASGU>
    save(fullfile(networkFolder, ...
        [participantIds{participanti} '_FREQ.mat']),'FREQ','participant');
end
networkSet = freqnessgui.inspectNetworkFolder(networkFolder);
moduleFolder = fullfile(networkSet.analysisFolder,'EntropyLandscape');
context = struct('frequencies',networkSet.frequencies, ...
    'nComponents',networkSet.nComponents);
configuration = defaultSecondaryConfiguration( ...
    'entropy',context,moduleFolder,participantIds);

guiRoot = fileparts(fileparts(mfilename('fullpath')));
functionFolder = fullfile(fileparts(guiRoot),'FREQNESS_Functions');
addpath(functionFolder);
pathCleanup = onCleanup(@()rmpath(functionFolder)); %#ok<NASGU>
previousVisibility = get(groot,'DefaultFigureVisible');
set(groot,'DefaultFigureVisible','off');
visibilityCleanup = onCleanup( ...
    @()set(groot,'DefaultFigureVisible',previousVisibility)); %#ok<NASGU>

report = freqnessgui.runSecondaryAnalysis( ...
    networkSet,configuration,[],[]);
verifyEqual(testCase,report.status,'completed');
verifyEqual(testCase,report.nCompleted,2);
verifyEqual(testCase,report.nFailed,0);
verifyEqual(testCase,report.nCancelled,0);
verifyTrue(testCase,isfile(report.groupOutputFile));
verifyTrue(testCase,isfile(report.manifestFile));

groupResult = load(report.groupOutputFile,'groupOutput');
verifySize(testCase,groupResult.groupOutput.H2,[3 2]);
for participanti = 1:numel(participantIds)
    participantResult = load(report.participants(participanti).outputFile, ...
        'participantOutput');
    verifyEqual(testCase,participantResult.participantOutput.H2, ...
        groupResult.groupOutput.H2(:,participanti));
end
loadedManifest = load(report.manifestFile,'manifest');
verifyEqual(testCase,loadedManifest.manifest.status,'completed');

unchangedManifest = freqnessgui.markSecondaryAnalysisCancelled( ...
    configuration,'Late cancellation must not replace completion.');
verifyEqual(testCase,unchangedManifest.status,'completed');
end

function testCancellationPersistsCancelledManifest(testCase)
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,'s')); %#ok<NASGU>
networkFolder = fullfile(temporaryRoot,'FREQ_Networks_Cancel');
mkdir(networkFolder);
participantIds = {'sub-001','sub-002'};
for participanti = 1:numel(participantIds)
    FREQ = syntheticFREQ(participanti); %#ok<NASGU>
    participant = struct('id',participantIds{participanti}); %#ok<NASGU>
    save(fullfile(networkFolder, ...
        [participantIds{participanti} '_FREQ.mat']),'FREQ','participant');
end
networkSet = freqnessgui.inspectNetworkFolder(networkFolder);
moduleFolder = fullfile(networkSet.analysisFolder,'EntropyLandscape');
context = struct('frequencies',networkSet.frequencies, ...
    'nComponents',networkSet.nComponents);
configuration = defaultSecondaryConfiguration( ...
    'entropy',context,moduleFolder,participantIds);
cancellationFile = freqnessgui.prepareSecondaryCancellation(moduleFolder);
progressFcn = @(~,message)cancelAfterFirstParticipant( ...
    message,cancellationFile);

report = freqnessgui.runSecondaryAnalysis( ...
    networkSet,configuration,[],progressFcn,@()isfile(cancellationFile));
verifyEqual(testCase,report.status,'cancelled');
verifyEqual(testCase,report.nCompleted,0);
verifyEqual(testCase,report.nFailed,0);
verifyEqual(testCase,report.nCancelled,2);
verifyFalse(testCase,isfile(report.groupOutputFile));

loadedManifest = load(report.manifestFile,'manifest');
verifyEqual(testCase,loadedManifest.manifest.status,'cancelled');
verifyTrue(testCase,all(strcmp( ...
    {loadedManifest.manifest.participants.status},'cancelled')));
end

function FREQ = syntheticFREQ(scaleFactor)
FREQ = struct();
FREQ.evals = scaleFactor*(reshape(1:12,4,3)+1);
FREQ.evecs = scaleFactor*reshape(1:24,4,2,3);
FREQ.pats = scaleFactor*reshape(1:24,4,2,3);
FREQ.ts = scaleFactor*reshape(1:120,2,20,3);
FREQ.frex = [2 4 8];
FREQ.fwhm = [0.5 0.8 1.2];
FREQ.srate = 100;
FREQ.duration = 0.2;
FREQ.bad_segments = [];
FREQ.regularisation = 0.01;
FREQ.scale_factors = scaleFactor;
end

function CFC = mockCFC()
CFC = struct();
CFC.PAC_all = zeros(3,2,4);
CFC.PAC_avg = zeros(3,4);
CFC.carrier_frex = [4 8 16];
CFC.lfo_freq = 2;
CFC.comp = 1;
CFC.phase_edges = 1:5;
CFC.phase_centers = 1:4;
CFC.fitted_PAC = zeros(3,2,4);
CFC.coefficients = zeros(3,2,3);
columnFields = {'amplitude_raw','amplitude_normalized','preferred_phase', ...
    'dc_offset','mse','r2','valid_bins','sAmpl','pShift','dcOff','mFrex', ...
    'goodFit'};
for fieldi = 1:numel(columnFields)
    CFC.(columnFields{fieldi}) = zeros(3,2);
end
CFC.lfo_phase = zeros(10,2);
end

function cancelAfterFirstParticipant(message,cancellationFile)
if contains(message,'Loaded participant #1/') && ~isfile(cancellationFile)
    freqnessgui.requestSecondaryCancellation(cancellationFile);
end
end

function configuration = defaultSecondaryConfiguration( ...
        moduleId,context,outputFolder,participants)
modules = freqnessgui.moduleRegistry();
module = modules(strcmp({modules.id},moduleId));
schema = freqnessgui.secondaryModuleSchema(moduleId,context);
fields = [schema.mandatory schema.optional];
values = struct();
selectedFrequencies = struct();
for fieldi = 1:numel(fields)
    field = fields(fieldi);
    if strcmp(field.type,'componentVector')
        values.(field.key) = freqnessgui.parseNumericVector( ...
            field.default,false,field.label);
    elseif strcmp(field.type,'numberOrEmpty')
        values.(field.key) = freqnessgui.parseNumericVector( ...
            field.default,true,field.label);
    else
        values.(field.key) = field.default;
    end
    if strcmp(field.type,'frequency')
        selectedFrequencies.(field.key) = field.default;
    elseif strcmp(field.type,'frequencyRange')
        [~,selectedFrequencies.(field.key)] = ...
            freqnessgui.mapFrequencySelection( ...
            context.frequencies,field.default);
    end
end
configuration = struct( ...
    'schemaVersion',1, ...
    'moduleId',moduleId, ...
    'functionName',module.functionName, ...
    'networkFolder','', ...
    'outputFolder',outputFolder, ...
    'values',values, ...
    'selectedFrequencies',selectedFrequencies, ...
    'participants',{participants});
end
