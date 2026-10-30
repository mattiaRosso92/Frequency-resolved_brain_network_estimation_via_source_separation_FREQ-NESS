function tests = test_FREQNESS_RepeatedSpatial
% Independent Statistics Toolbox references for the public spatial API.
tests = functiontests(localfunctions);
end

function setupOnce(~)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'FREQNESS_Toolbox','FREQNESS_Functions'));
end

function testMatchesRanovaAcrossFactorLayouts(testCase)
layouts = {[2 4],3,[2 3 2]};
for layouti = 1:numel(layouts)
    shape = layouts{layouti};
    [grid,ids,MNI] = dataset(shape,14);
    stats = run_comparison(grid,ids,MNI,'normalization','none');
    names = stats.factor_names;
    within = array2table(stats.contrast_cell_levels,'VariableNames',names);
    for factori = 1:numel(names)
        within.(names{factori}) = categorical(within.(names{factori}));
    end
    % Two voxels and both frequency/component axes check flattening too.
    for mapi = [1 2 9 17]
        Y = zeros(14,numel(grid));
        for celli = 1:numel(grid)
            maps = permute(abs(grid{celli}.pats),[1 3 2 4]);
            flat = reshape(maps,[],14);
            Y(:,celli) = flat(mapi,:)';
        end
        variables = arrayfun(@(i)sprintf('Y%d',i),1:numel(grid),'UniformOutput',false);
        measurements = array2table(Y,'VariableNames',variables);
        model = fitrm(measurements,sprintf('Y1-Y%d~1',numel(grid)), ...
            'WithinDesign',within);
        reference = ranova(model,'WithinModel',strjoin(names,'*'));
        for effecti = 1:numel(stats.effect_names)
            effect = stats.effects.(stats.effect_names{effecti});
            term = strjoin(names(effect.factors),':');
            row = reference(['(Intercept):' term],:);
            error_row = reference(['Error(' term ')'],:);
            adjustment = epsilon(model,effect.contrast_matrix);
            verifyEqual(testCase,effect.statistic(mapi),row.F,'AbsTol',1e-8);
            verifyEqual(testCase,effect.ss_effect(mapi),row.SumSq,'AbsTol',1e-9);
            verifyEqual(testCase,effect.ss_error(mapi),error_row.SumSq,'AbsTol',1e-9);
            verifyEqual(testCase,effect.df1,row.DF);
            verifyEqual(testCase,effect.df2,error_row.DF);
            verifyEqual(testCase,effect.effect_size(mapi), ...
                row.SumSq/(row.SumSq+error_row.SumSq),'AbsTol',1e-11);
            verifyEqual(testCase,effect.epsilon_gg(mapi), ...
                adjustment.GreenhouseGeisser,'AbsTol',1e-10);
            verifyEqual(testCase,effect.p_uncorrected(mapi),row.pValueGG,'AbsTol',1e-10);
        end
    end
end
end

function testDrugMatchesPairedTestAfterCellNormalization(testCase)
[grid,ids,MNI] = dataset([2 4],14);
stats = run_comparison(grid,ids,MNI,'factor_names',{'Drug','Condition'});
averages = cell(1,2);
for drugi = 1:2
    averages{drugi} = zeros(8,2,2,14);
    for condi = 1:4
        maps = permute(abs(grid{drugi,condi}.pats),[1 3 2 4]);
        averages{drugi} = averages{drugi}+maps./sqrt(sum(maps.^2,1))/4;
    end
end
[~,~,~,reference] = ttest(averages{2}-averages{1},0,'Dim',4);
verifyEqual(testCase,stats.effects.Drug.statistic,reference.tstat.^2,'AbsTol',1e-9);
verifyEqual(testCase,stats.effects.Drug.df2,13);
verifyEqual(testCase,stats.effects.Condition.df2,39);
verifyEqual(testCase,stats.effects.Drug_x_Condition.df2,39);
verifyEqual(testCase,stats.inference_method,'within_subject_contrast_signflip');
verifyEqual(testCase,stats.permutation_unit,'participant');
verifySize(testCase,stats.effects.Drug.statistic,[8 2 2]);
end

function testRepeatedConditionsDoNotInflateDrugEvidence(testCase)
[grid,ids,MNI] = dataset([2 4],14);
difference = reshape(linspace(-1,2,14),[1 1 1 14]);
for drugi = 1:2
    for condi = 1:4
        grid{drugi,condi}.pats(:,:,:,:) = repmat(20+(drugi-1)*difference,[8 2 2 1]);
    end
end
stats = run_comparison(grid,ids,MNI,'normalization','none');
[~,~,~,reference] = ttest(difference,0,'Dim',4);
verifyEqual(testCase,stats.effects.Factor1.statistic, ...
    repmat(reference.tstat.^2,[8 2 2]),'AbsTol',1e-10);
verifyEqual(testCase,stats.effects.Factor2.statistic,zeros(8,2,2));
verifyEqual(testCase,stats.effects.Factor1_x_Factor2.statistic,zeros(8,2,2));
end

function testPureConditionAndInteractionStaySeparate(testCase)
[grid,ids,MNI] = dataset([2 4],6);
amplitude = reshape(linspace(.1,1.1,6),[1 1 1 6]);
condition = [-3 -1 1 3];
for interaction_flag = [false true]
    for drugi = 1:2
        for condi = 1:4
            drug_weight = 1;
            if interaction_flag
                drug_weight = 2*drugi-3;
            end
            grid{drugi,condi}.pats = repmat( ...
                20+drug_weight*condition(condi)*amplitude,[8 2 2 1]);
        end
    end
    stats = run_comparison(grid,ids,MNI,'normalization','none');
    verifyEqual(testCase,stats.effects.Factor1.statistic,zeros(8,2,2));
    if interaction_flag
        verifyEqual(testCase,stats.effects.Factor2.statistic,zeros(8,2,2));
        verifyGreaterThan(testCase,stats.effects.Factor1_x_Factor2.statistic,0);
    else
        verifyEqual(testCase,stats.effects.Factor1_x_Factor2.statistic,zeros(8,2,2));
        verifyGreaterThan(testCase,stats.effects.Factor2.statistic,0);
    end
end
end

function testConstantAndZeroContrasts(testCase)
[grid,ids,MNI] = dataset([2 4],6);
for celli = 1:8
    grid{celli}.pats(:) = 20;
end
zero = run_comparison(grid,ids,MNI,'normalization','none');
verifyEqual(testCase,zero.effects.Factor1.statistic,zeros(8,2,2));
verifyEqual(testCase,zero.effects.Factor1.p_uncorrected,ones(8,2,2));
for condi = 1:4
    grid{2,condi}.pats(:) = 22;
end
constant = run_comparison(grid,ids,MNI,'normalization','none');
verifyTrue(testCase,all(isinf(constant.effects.Factor1.statistic(:))));
verifyEqual(testCase,constant.effects.Factor1.p_uncorrected,zeros(8,2,2));
verifyEqual(testCase,constant.effects.Factor1.effect_size,ones(8,2,2));
end

function testExactReorderingSignScaleAndSeedContracts(testCase)
[grid,ids,MNI] = dataset([2 3],6);
baseline = run_comparison(grid,ids,MNI,'n_permutations','exact', ...
    'cluster_forming_alpha',.2);
verifyTrue(testCase,baseline.exact_permutations);
verifyEqual(testCase,baseline.n_permutations,64);
for celli = 1:numel(grid)
    order = randperm(6);
    grid{celli}.pats = -3*grid{celli}.pats(:,:,:,order);
    ids{celli} = ids{celli}(order);
end
changed = run_comparison(grid,ids,MNI,'n_permutations','exact', ...
    'cluster_forming_alpha',.2);
for effecti = 1:numel(baseline.effect_names)
    key = baseline.effect_names{effecti};
    verifyEqual(testCase,changed.effects.(key).statistic, ...
        baseline.effects.(key).statistic,'AbsTol',1e-9);
    verifyEqual(testCase,changed.effects.(key).p_cluster_corrected, ...
        baseline.effects.(key).p_cluster_corrected);
end
first = run_comparison(grid,ids,MNI,'random_seed',37);
second = run_comparison(grid,ids,MNI,'random_seed',37);
verifyEqual(testCase,first,second);
end

function testZeroPatternsAndMismatchedParticipantsFail(testCase)
[grid,ids,MNI] = dataset([2 4],6);
ids{2,4}(6) = 99;
did_error = false;
try
    run_comparison(grid,ids,MNI);
catch exception
    did_error = contains(exception.message,'same participant');
end
verifyTrue(testCase,did_error);
ids{2,4}(6) = 6;
grid{2,4}.pats(:) = 0;
did_error = false;
try
    run_comparison(grid,ids,MNI);
catch exception
    did_error = contains(exception.message,'zero L2 norm');
end
verifyTrue(testCase,did_error);
end

function stats = run_comparison(grid,ids,MNI,varargin)
stats = FREQNESS_SpatialComparison(grid,MNI,'frex2compare','all', ...
    'design','within','test','anova','subject_ids',ids, ...
    'comps2compare',[1 2],'n_permutations',7,'plot_results',false,varargin{:});
end

function [grid,ids,MNI] = dataset(shape,nsubs)
rng(17,'twister');
if isscalar(shape)
    shape = [shape 1];
end
grid = cell(shape);
ids = cell(shape);
[x,y,z] = ndgrid(0:1,0:1,0:1);
MNI = 8*[x(:) y(:) z(:)];
baseline = 20+randn(8,2,2,nsubs);
for celli = 1:numel(grid)
    grid{celli}.pats = baseline+.4*randn(size(baseline))+.05*celli;
    grid{celli}.frex = [8 9];
    ids{celli} = 1:nsubs;
end
end
