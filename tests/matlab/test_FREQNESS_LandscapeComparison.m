function tests = test_FREQNESS_LandscapeComparison

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

function testWithinSubjectTTestAndFDR(testCase)

[FREQ_A,FREQ_B] = paired_freq;
STATS = FREQNESS_LandscapeComparison(FREQ_A,FREQ_B, ...
    'measure','evals','test','ttest','which_comp',[1 2], ...
    'plot_results',false);

verifyEqual(testCase,STATS.design,'within');
verifyEqual(testCase,STATS.test,'ttest');
verifyEqual(testCase,STATS.correction,'fdr');
verifyEqual(testCase,STATS.measures,{'evals'});
verifyEqual(testCase,STATS.evals.comps,[1 2]);
verifySize(testCase,STATS.evals.statistic,[2 4]);
verifyTrue(testCase,all(STATS.evals.significant(1,2:3)));
verifyFalse(testCase,any(STATS.evals.significant(2,:)));

end

function testOmittedMeasurePreservesEntropyFields(testCase)

rng(8,'twister');
ENTR_A.frex = [2 4 6];
ENTR_A.H2 = 1+.03*randn(3,10);
ENTR_A.ED = 3+.08*randn(3,10);
ENTR_B.frex = ENTR_A.frex;
ENTR_B.H2 = ENTR_A.H2+.2;
ENTR_B.ED = ENTR_A.ED-.3;

STATS = FREQNESS_LandscapeComparison(ENTR_A,ENTR_B, ...
    'plot_results',false);

verifyEqual(testCase,STATS.measures,{'H2','ED'});
verifyEqual(testCase,STATS.H2.measure,'H2');
verifyEqual(testCase,STATS.ED.measure,'ED');
verifyEqual(testCase,STATS.H2.comps,1);
verifyEqual(testCase,STATS.H2.frex,[2 4 6]);

end

function testBetweenSubjectDesignAllowsDifferentGroupSizes(testCase)

rng(19,'twister');
GROUP_A.frex = [5 10 15];
GROUP_A.stability_index = .4+.04*randn(2,3,9);
GROUP_B.frex = GROUP_A.frex;
GROUP_B.stability_index = .65+.04*randn(2,3,13);

nonparametric = FREQNESS_LandscapeComparison(GROUP_A,GROUP_B, ...
    'measure','stability_index','design','between','which_comp',2, ...
    'plot_results',false);
parametric = FREQNESS_LandscapeComparison(GROUP_A,GROUP_B, ...
    'measure','stability_index','design','between','test','ttest', ...
    'which_comp',2,'plot_results',false);

verifyEqual(testCase,nonparametric.test,'ranksum');
verifyEqual(testCase,nonparametric.stability_index.statistic_name,'U');
verifyEqual(testCase,parametric.stability_index.effect_size_name,'Hedges'' g');
verifyEqual(testCase,nonparametric.stability_index.n_A,[9 9 9]);
verifyEqual(testCase,nonparametric.stability_index.n_B,[13 13 13]);

end

function testCrossCouplingMeasuresUseCarrierFrequencyAxis(testCase)

rng(22,'twister');
CFC_A.carrier_frex = [12 16 20];
CFC_A.amplitude_raw = .1+.01*randn(3,8);
CFC_A.amplitude_normalized = 10+randn(3,8);
CFC_B.carrier_frex = CFC_A.carrier_frex;
CFC_B.amplitude_raw = CFC_A.amplitude_raw+.03;
CFC_B.amplitude_normalized = CFC_A.amplitude_normalized+2;

STATS = FREQNESS_LandscapeComparison(CFC_A,CFC_B, ...
    'plot_results',false);

verifyEqual(testCase,STATS.measures, ...
    {'amplitude_raw','amplitude_normalized'});
verifyEqual(testCase,STATS.amplitude_raw.frex,[12 16 20]);
verifyEqual(testCase,STATS.amplitude_normalized.frex,[12 16 20]);

end

function testRawNumericLandscapeAndFigures(testCase)

values_a = [1 1.1 .9;2 2.1 1.9];
values_b = values_a+.3;
STATS = FREQNESS_LandscapeComparison(values_a,values_b, ...
    'frex',[6 12],'plot_results',true);

verifyEqual(testCase,STATS.measures,{'values'});
verifyEqual(testCase,STATS.values.frex,[6 12]);
verifyEqual(testCase,numel(findall(groot,'Type','figure')),1);
close all

end

function testIncompatibleNonparametricTestsWarnAndAreCorrected(testCase)

[FREQ_A,FREQ_B] = paired_freq;
lastwarn('');
within = FREQNESS_LandscapeComparison(FREQ_A,FREQ_B, ...
    'design','within','test','ranksum','plot_results',false);
[warning_message,~] = lastwarn;
verifyTrue(testCase,contains(warning_message,'Using test=''signrank'''));
verifyEqual(testCase,within.test,'signrank');

FREQ_B.evals = FREQ_B.evals(:,:,1:10);
lastwarn('');
between = FREQNESS_LandscapeComparison(FREQ_A,FREQ_B, ...
    'design','between','test','signrank','plot_results',false);
[warning_message,~] = lastwarn;
verifyTrue(testCase,contains(warning_message,'Using test=''ranksum'''));
verifyEqual(testCase,between.test,'ranksum');

end

function [FREQ_A,FREQ_B] = paired_freq

rng(31,'twister');
FREQ_A.frex = [4 8 12 16];
FREQ_A.evals = .3+.04*randn(2,4,14);
FREQ_B.frex = FREQ_A.frex;
FREQ_B.evals = FREQ_A.evals+.01*randn(size(FREQ_A.evals));
FREQ_B.evals(1,2:3,:) = FREQ_B.evals(1,2:3,:)+.15;

end
