function tests = test_FREQNESS_SpatialComparison

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

function testReturnsClusterCorrectedStructure(testCase)

[FREQ_A,FREQ_B,MNI,central_cluster] = synthetic_paired_data;
STATS = FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI, ...
    'frex2compare',8,'normalization','none', ...
    'cluster_forming_alpha',.05,'n_permutations','exact', ...
    'plot_results',false);

verifyEqual(testCase,STATS.design,'within');
verifyEqual(testCase,STATS.test,'ttest');
verifyEqual(testCase,STATS.n_participants,6);
verifyEqual(testCase,STATS.n_A,6);
verifyEqual(testCase,STATS.n_B,6);
verifyTrue(testCase,STATS.exact_permutations);
verifyEqual(testCase,STATS.n_permutations,64);
verifySize(testCase,STATS.mean_difference,[27 1]);
verifyEqual(testCase,STATS.comps,1);
verifyEqual(testCase,STATS.idx_frex,2);
verifyTrue(testCase,any(STATS.significant_mask(central_cluster)));
verifyTrue(testCase,any([STATS.clusters.p_corrected]<=.05));

end

function testClosestFrequencyAndMultipleMaps(testCase)

[FREQ_A,FREQ_B,MNI] = synthetic_paired_data;
STATS = FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI, ...
    'frex2compare',[7.8 12],'comps2compare',[1 2], ...
    'n_permutations',8,'plot_results',false);

verifyEqual(testCase,STATS.requested_frex,[7.8;12]);
verifyEqual(testCase,STATS.frex,[8;12]);
verifyEqual(testCase,STATS.idx_frex,[2;3]);
verifyEqual(testCase,STATS.comps,[1 2]);
verifySize(testCase,STATS.mean_difference,[27 2 2]);
verifyEqual(testCase,STATS.correction_scope,'all_maps');

end

function testMagnitudeL2ComparisonIsSignAndScaleInvariant(testCase)

[FREQ_A,FREQ_B,MNI] = synthetic_paired_data;
baseline = FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI, ...
    'frex2compare',8,'n_permutations',8,'plot_results',false);
FREQ_A.pats = -3*FREQ_A.pats;
FREQ_B.pats = 7*FREQ_B.pats;
transformed = FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI, ...
    'frex2compare',8,'n_permutations',8,'plot_results',false);

verifyEqual(testCase,transformed.mean_difference,baseline.mean_difference, ...
    'AbsTol',1e-12);
verifyEqual(testCase,transformed.t_values,baseline.t_values,'AbsTol',1e-10);
verifyEqual(testCase,transformed.p_cluster_corrected, ...
    baseline.p_cluster_corrected,'AbsTol',1e-12);

end

function testSignedAlignmentRemovesArbitrarySign(testCase)

[FREQ_A,~,MNI] = synthetic_paired_data;
FREQ_B = FREQ_A;
FREQ_B.pats = -FREQ_A.pats;
STATS = FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI, ...
    'frex2compare',4,'map_mode','signed_aligned', ...
    'n_permutations',8,'plot_results',false);

verifyEqual(testCase,STATS.mean_difference,zeros(27,1),'AbsTol',1e-14);
verifyTrue(testCase,all(STATS.sign_flips_b(:)==-1));
verifyFalse(testCase,any(STATS.significant_mask(:)));

end

function testAllFrequenciesProduceFigures(testCase)

[FREQ_A,FREQ_B,MNI] = synthetic_paired_data;
STATS = FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI, ...
    'frex2compare','all','comps2compare',[1 2], ...
    'n_permutations',8,'plot_results',true);

verifyEqual(testCase,STATS.frex,FREQ_A.frex(:));
verifyEqual(testCase,numel(findall(groot,'Type','figure')),6);
close all

end

function testMandatoryFrequencyAndRangeSafeguards(testCase)

[FREQ_A,FREQ_B,MNI] = synthetic_paired_data;
did_error = false;
try
    FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI,'plot_results',false);
catch ME
    did_error = contains(ME.message,'mandatory');
end
verifyTrue(testCase,did_error);

did_error = false;
try
    FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI, ...
        'frex2compare',30,'plot_results',false);
catch ME
    did_error = contains(ME.message,'outside');
end
verifyTrue(testCase,did_error);

end

function testBetweenSubjectTTestSupportsUnequalGroups(testCase)

[FREQ_A,FREQ_B,MNI] = synthetic_between_data;
STATS = FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI, ...
    'frex2compare',8,'design','between','test','ttest', ...
    'normalization','none','cluster_forming_alpha',.05, ...
    'n_permutations','exact','plot_results',false);

verifyEqual(testCase,STATS.design,'between');
verifyEqual(testCase,STATS.test,'ttest');
verifyEqual(testCase,STATS.n_participants,[4 5]);
verifyEqual(testCase,STATS.n_A,4);
verifyEqual(testCase,STATS.n_B,5);
verifyEqual(testCase,STATS.statistic_name,'t');
verifyEqual(testCase,STATS.effect_size_name,'Hedges'' g');
verifyTrue(testCase,STATS.exact_permutations);
verifyEqual(testCase,STATS.n_permutations,126);
verifyEqual(testCase,size(STATS.sign_flips_a,3),4);
verifyEqual(testCase,size(STATS.sign_flips_b,3),5);

end

function testIncompatibleNonparametricTestsWarnAndAreCorrected(testCase)

[FREQ_A,FREQ_B,MNI] = synthetic_paired_data;
lastwarn('');
within = FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI, ...
    'frex2compare',8,'design','within','test','ranksum', ...
    'cluster_forming_alpha',.05,'n_permutations',8, ...
    'plot_results',false);
[warning_message,~] = lastwarn;
verifyTrue(testCase,contains(warning_message,'Using test=''signrank'''));
verifyEqual(testCase,within.test,'signrank');
verifyEqual(testCase,within.statistic_name,'z');
verifyEqual(testCase,within.effect_size_name,'matched rank-biserial r');
verifyTrue(testCase,all(isnan(within.t_values(:))));

FREQ_B.pats = FREQ_B.pats(:,:,:,1:5);
lastwarn('');
between = FREQNESS_SpatialComparison(FREQ_A,FREQ_B,MNI, ...
    'frex2compare',8,'design','between','test','signrank', ...
    'cluster_forming_alpha',.05,'n_permutations',8, ...
    'plot_results',false);
[warning_message,~] = lastwarn;
verifyTrue(testCase,contains(warning_message,'Using test=''ranksum'''));
verifyEqual(testCase,between.test,'ranksum');
verifyEqual(testCase,between.statistic_name,'z');
verifyEqual(testCase,between.effect_size_name,'rank-biserial r');

end

function [FREQ_A,FREQ_B,MNI,central_cluster] = synthetic_paired_data

[x,y,z] = ndgrid(0:2,0:2,0:2);
MNI = 8*[x(:) y(:) z(:)];
frex = [4 8 12];
rng(12,'twister');
patterns_a = .2+.03*randn(27,2,3,6);
patterns_b = patterns_a+.015*randn(size(patterns_a));
central_cluster = [5 6 14 15];
patterns_b(central_cluster,1,2,:) = ...
    patterns_b(central_cluster,1,2,:)+.25;
FREQ_A.pats = patterns_a;
FREQ_A.frex = frex;
FREQ_B.pats = patterns_b;
FREQ_B.frex = frex;

end

function [FREQ_A,FREQ_B,MNI] = synthetic_between_data

[x,y,z] = ndgrid(0:2,0:2,0:2);
MNI = 8*[x(:) y(:) z(:)];
rng(44,'twister');
FREQ_A.pats = .2+.03*randn(27,1,1,4);
FREQ_A.frex = 8;
FREQ_B.pats = .2+.03*randn(27,1,1,5);
FREQ_B.pats([5 6 14 15],1,1,:) = ...
    FREQ_B.pats([5 6 14 15],1,1,:)+.3;
FREQ_B.frex = 8;

end
