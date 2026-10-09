function tests = test_FREQNESS_CenterOfMassGradients

tests = functiontests(localfunctions);

end
function setupOnce(testCase)

testDirectory = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(fileparts(testDirectory));
addpath(fullfile(repositoryRoot,'FREQNESS_Toolbox','FREQNESS_Functions'));
testCase.TestData.figureVisibility = get(groot,'defaultFigureVisible');
set(groot,'defaultFigureVisible','off');

end

function teardownOnce(testCase)

set(groot,'defaultFigureVisible',testCase.TestData.figureVisibility);
close all

end

function testFrequencyGradientUsesPhysicalFrequency(testCase)

[FREQ,MNI] = synthetic_gradient_data;
[coefficients,goodFit] = FREQNESS_FreqGradients(FREQ,MNI,'comp2model',1);

verifySize(testCase,coefficients,[3 3 2]);
verifyEqual(testCase,squeeze(coefficients(1,2,:)),[1.5;1.5],'AbsTol',0.02);
verifyEqual(testCase,goodFit.bestOrder(1,:),[1 1]);
verifyGreaterThan(testCase,goodFit.R2_best(1,:),0.99);
close all

end

function testComponentGradientUsesCenterTrajectory(testCase)

[FREQ,MNI] = synthetic_gradient_data;
[coefficients,goodFit] = FREQNESS_CompGradients(FREQ,MNI,'freq2model',8);

verifySize(testCase,coefficients,[3 3 2]);
verifyEqual(testCase,squeeze(coefficients(1,2,:)),[8;8],'AbsTol',0.05);
verifyEqual(testCase,goodFit.bestOrder(1,:),[1 1]);
verifyGreaterThan(testCase,goodFit.R2_best(1,:),0.99);
close all

end

function testCenterTrajectoryIsSignAndScaleInvariant(testCase)

[FREQ,MNI] = synthetic_gradient_data;
baseline = FREQNESS_FreqGradients(FREQ,MNI,'comp2model',1);
FREQ.pats = -7.5*FREQ.pats;
transformed = FREQNESS_FreqGradients(FREQ,MNI,'comp2model',1);

verifyEqual(testCase,transformed,baseline,'AbsTol',1e-10);
close all

end

function testZeroPatternsReturnNaN(testCase)

FREQ.pats = zeros(20,3,4);
FREQ.frex = [2 4 6 8];
MNI = [(0:19)' (0:19)' (0:19)'];
[coefficients,goodFit] = FREQNESS_FreqGradients(FREQ,MNI);

verifyTrue(testCase,all(isnan(coefficients(:))));
verifyTrue(testCase,all(isnan(goodFit.R2_best(:))));
close all

end

function [FREQ,MNI] = synthetic_gradient_data

x = linspace(-45,45,181)';
MNI = [x 0.5*x+2 -0.25*x+5];
frex = [2 4 8 12 20];
nComp = 4;
nSubs = 2;
patterns = zeros(numel(x),nComp,numel(frex),nSubs);

for subi = 1:nSubs
    for compi = 1:nComp
        for frexi = 1:numel(frex)
            center = -25+1.5*frex(frexi)+8*(compi-1)+0.5*(subi-1);
            patterns(:,compi,frexi,subi) = ...
                exp(-0.5*((x-center)/2.5).^2)+0.01;
        end
    end
end

FREQ.pats = patterns;
FREQ.frex = frex;
FREQ.evals = ones(nComp,numel(frex),nSubs);
FREQ.evals(1,2,:) = 8;

end

