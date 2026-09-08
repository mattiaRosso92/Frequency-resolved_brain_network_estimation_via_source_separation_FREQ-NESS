function tests = test_FREQNESS_FilterWidths
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
test_directory = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(test_directory));
toolbox_root = fullfile(repo_root,'FREQNESS_Toolbox');

addpath(fullfile(toolbox_root,'FREQNESS_Functions'));
addpath(fullfile(toolbox_root,'FREQNESS_ExternalFunctions'));

testCase.TestData.toolbox_root = toolbox_root;
testCase.TestData.figure_visibility = get(groot,'defaultFigureVisible');
set(groot,'defaultFigureVisible','off');
end

function teardownOnce(testCase)
close all;
set(groot,'defaultFigureVisible',testCase.TestData.figure_visibility);
end

function testLogarithmicReferenceCalibration(testCase)
frex = [1 10 100];
expected = [0.142857142857143 2.020305089104422 28.571428571428573];

actual = FREQNESS_ComputeFilterWidths(frex,'logarithmic');

verifyEqual(testCase,actual,expected,'AbsTol',1e-12);
end

function testAutomaticWidthsAreFrequencyDensityInvariant(testCase)
frex_sparse = 1:1:100;
frex_dense = 1:0.25:100;

widths_sparse = FREQNESS_ComputeFilterWidths(frex_sparse,'logarithmic');
widths_dense = FREQNESS_ComputeFilterWidths(frex_dense,'logarithmic');

verifyEqual(testCase,widths_sparse,widths_dense(1:4:end),'AbsTol',0);
end

function testLinearModeUsesConstantQ(testCase)
frex = [1 10 100];
actual = FREQNESS_ComputeFilterWidths(frex,'linear');

verifyEqual(testCase,actual,frex/5,'AbsTol',0);
end

function testSingleFrequencyUsesLogarithmicMidpointQ(testCase)
actual = FREQNESS_ComputeFilterWidths(10,'logarithmic');
expected = 10/sqrt(7*3.5);

verifyEqual(testCase,actual,expected,'AbsTol',1e-12);
end

function testNetworkEstimationIsDensityInvariant(testCase)
data = deterministic_data();
frex_sparse = [2 4 6];
frex_dense = [2 3 4 5 6];

FREQ_sparse = FREQNESS_NetworkEstimation(data,frex_sparse,64, ...
    'ncomps',3);
FREQ_dense = FREQNESS_NetworkEstimation(data,frex_dense,64, ...
    'ncomps',3);
idx_shared = 1:2:numel(frex_dense);

verifyEqual(testCase,FREQ_sparse.fwhm,FREQ_dense.fwhm(idx_shared), ...
    'AbsTol',0);
verifyEqual(testCase,FREQ_sparse.evals,FREQ_dense.evals(:,idx_shared,:), ...
    'AbsTol',1e-12);
verifyEqual(testCase,FREQ_sparse.evecs, ...
    FREQ_dense.evecs(:,:,idx_shared,:),'AbsTol',1e-12);
verifyEqual(testCase,FREQ_sparse.pats,FREQ_dense.pats(:,:,idx_shared,:), ...
    'AbsTol',1e-12);
verifyEqual(testCase,FREQ_sparse.ts,FREQ_dense.ts(:,:,idx_shared,:), ...
    'AbsTol',1e-12);
end

function testExplicitWidthsPreserveReferenceOutputs(testCase)
fixture_file = fullfile(testCase.TestData.toolbox_root,'python','tests', ...
    'fixtures','matlab_core_reference.mat');
reference = load(fixture_file);

actual = FREQNESS_NetworkEstimation(reference.data,reference.frex, ...
    reference.srate,'fwidth',reference.fwidth,'regularisation',0.01, ...
    'ncomps',3,'bad_segments',reference.bad_segments);

verifyEqual(testCase,actual.frex,reference.FREQ.frex,'AbsTol',0);
verifyEqual(testCase,actual.fwhm,reference.FREQ.fwhm,'AbsTol',0);
verifyEqual(testCase,actual.evals,reference.FREQ.evals, ...
    'RelTol',1e-10,'AbsTol',1e-10);
verifyEqual(testCase,actual.pats,reference.FREQ.pats, ...
    'RelTol',1e-9,'AbsTol',1e-10);

[aligned_evecs,aligned_ts] = align_component_signs( ...
    actual.evecs,actual.ts,reference.FREQ.evecs);
verifyEqual(testCase,aligned_evecs,reference.FREQ.evecs, ...
    'RelTol',1e-9,'AbsTol',1e-10);
verifyEqual(testCase,aligned_ts,reference.FREQ.ts, ...
    'RelTol',1e-9,'AbsTol',1e-10);
end

function testExplicitWidthsOverrideAutomaticMode(testCase)
data = deterministic_data();
frex = [2 4 6];
fwidth = [0.4 0.8 1.2];

FREQ_log = FREQNESS_NetworkEstimation(data,frex,64, ...
    'fwidth',fwidth,'filter','logarithmic','ncomps',3);
FREQ_linear = FREQNESS_NetworkEstimation(data,frex,64, ...
    'fwidth',fwidth,'filter','linear','ncomps',3);
FREQ_scalar = FREQNESS_NetworkEstimation(data,frex,64, ...
    'fwidth',0.8,'ncomps',1);

verifyEqual(testCase,FREQ_log.fwhm,fwidth,'AbsTol',0);
verifyEqual(testCase,FREQ_linear.fwhm,fwidth,'AbsTol',0);
verifyEqual(testCase,FREQ_scalar.fwhm,[0.8 0.8 0.8],'AbsTol',0);
verifyEqual(testCase,FREQ_log.evals,FREQ_linear.evals,'AbsTol',0);
verifyEqual(testCase,FREQ_log.evecs,FREQ_linear.evecs,'AbsTol',0);
verifyEqual(testCase,FREQ_log.pats,FREQ_linear.pats,'AbsTol',0);
verifyEqual(testCase,FREQ_log.ts,FREQ_linear.ts,'AbsTol',0);
end

function testDownstreamFunctionsReusePhysicalWidths(testCase)
data = deterministic_data();
FREQ = FREQNESS_NetworkEstimation(data,[2 4 6],64,'ncomps',2);

CFC = FREQNESS_CrossCoupling(FREQ,2,'frex2model',[4 6], ...
    'which_comp',1,'plot_all',false,'nbins',12,'min_valid_bins',6);
verifyEqual(testCase,CFC.lfo_freq,2,'AbsTol',0);
verifyEqual(testCase,CFC.carrier_frex,[4; 6],'AbsTol',0);

IND = FREQNESS_InducedResponses(FREQ,[160 320],[-0.5 0.5], ...
    [-0.5 0],'which_comp',1,'frex2model',[2 6], ...
    'plot_avg',false,'plot_all',false);
verifyEqual(testCase,IND.frex,[2; 4; 6],'AbsTol',0);
verifyEqual(testCase,IND.fwhm,FREQ.fwhm(:),'AbsTol',0);
close all;
end

function data = deterministic_data
srate = 64;
time = (0:511)/srate;
data = [ ...
    sin(2*pi*2*time) + 0.20*cos(2*pi*11*time); ...
    cos(2*pi*3*time) + 0.10*sin(2*pi*13*time); ...
    0.70*sin(2*pi*4*time + 0.3) + 0.30*cos(2*pi*7*time); ...
    0.40*cos(2*pi*6*time - 0.2) + 0.25*sin(2*pi*9*time)];
end

function [aligned_evecs,aligned_ts] = align_component_signs( ...
        actual_evecs,actual_ts,reference_evecs)
aligned_evecs = actual_evecs;
aligned_ts = actual_ts;

for subi = 1:size(actual_evecs,4)
    for frexi = 1:size(actual_evecs,3)
        for compi = 1:size(actual_evecs,2)
            actual = actual_evecs(:,compi,frexi,subi);
            reference = reference_evecs(:,compi,frexi,subi);
            if dot(actual,reference) < 0
                aligned_evecs(:,compi,frexi,subi) = -actual;
                aligned_ts(compi,:,frexi,subi) = ...
                    -aligned_ts(compi,:,frexi,subi);
            end
        end
    end
end
end
