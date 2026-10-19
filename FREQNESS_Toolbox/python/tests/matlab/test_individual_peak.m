function test_individual_peak(output_dir)
% Deterministic MATLAB regression checks and optional cross-language fixture.

path_toolbox = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(fullfile(path_toolbox,'FREQNESS_Functions'))
frex = [4 7 8 9.5 11 13 16];
values = [999 1 2 9 4 2 999; ...
          999 9 1 9 1 0 999; ...
          999 NaN Inf -Inf NaN NaN 999; ...
          999 0 1 2 3 10 999; ...
          999 -5 -2 -1 -2 -4 999]';
FREQ = struct;
FREQ.frex = frex;
FREQ.evals = permute(cat(3,values,flipud(values)),[3 1 2]);
original = FREQ.evals;

PEAK = FREQNESS_IndividualPeak(FREQ,[7 13]);
assert(isequaln(PEAK.peak_frex,[9.5 7 NaN 13 9.5]))
assert(isequaln(PEAK.peak_value,[9 9 NaN 10 -1]))
assert(isequaln(PEAK.idx_peak,[4 2 NaN 6 4]))
assert(isequal(PEAK.at_boundary,logical([0 1 0 1 0])))
assert(isequal(PEAK.is_tied,logical([0 1 0 0 0])))
assert(isequal(PEAK.frex,[7 8 9.5 11 13]'))
assert(isequaln(FREQ.evals,original))

inside = FREQNESS_IndividualPeak(FREQ,[7.1 12.9]);
assert(isequal(inside.frex,[8 9.5 11]'))
assert(isequaln(inside.peak_frex,[9.5 9.5 NaN 11 9.5]))

single = FREQ;
single.evals = FREQ.evals(:,:,1);
selected = FREQNESS_IndividualPeak(single,[7 13],'which_comp',[2 1 2]);
assert(isequal(size(selected.peak_frex),[2 1]))
assert(isequal(selected.which_comp,[2 1]))

flat = struct('evals',reshape([Inf 5 NaN; 2 2 2]',[1 3 2]), ...
    'frex',[7 9 13]);
flat_peak = FREQNESS_IndividualPeak(flat,[7 13]);
assert(isequal(flat_peak.peak_frex,[9 7]))
assert(isequal(flat_peak.is_tied,logical([0 1])))

expect_error(@()FREQNESS_IndividualPeak(FREQ,[13 7]));
expect_error(@()FREQNESS_IndividualPeak(FREQ,[3 13]));
expect_error(@()FREQNESS_IndividualPeak(FREQ,[7 17]));
expect_error(@()FREQNESS_IndividualPeak(FREQ,[8.1 9.6]));
expect_error(@()FREQNESS_IndividualPeak(FREQ,[7 13],'which_comp',0));
expect_error(@()FREQNESS_IndividualPeak(FREQ,[7 13],'plot_avg',1));
expect_error(@()FREQNESS_IndividualPeak(FREQ,[7 13],'unknown',true));

% Exercise all plot modes, including an entirely missing participant.
visibility = get(groot,'defaultFigureVisible');
cleanup = onCleanup(@()set(groot,'defaultFigureVisible',visibility));
set(groot,'defaultFigureVisible','off')
plot_freq = FREQ;
plot_freq.evals = FREQ.evals(:,:,1:3);
for flags = [true false; false true; true true]'
    before = findall(groot,'Type','figure');
    state_before = rng;
    plotted = FREQNESS_IndividualPeak(plot_freq,[7 13], ...
        'which_comp',[1 2],'plot_avg',flags(1),'plot_all',flags(2));
    assert(isequal(rng,state_before))
    plain = FREQNESS_IndividualPeak(plot_freq,[7 13],'which_comp',[1 2]);
    assert(isequaln(plotted.peak_frex,plain.peak_frex))
    created = setdiff(findall(groot,'Type','figure'),before);
    assert(numel(created)==2*flags(1)+6*flags(2))
    close(created)
end

if nargin>0
    save(fullfile(output_dir,'individual_peak_reference.mat'),'FREQ','PEAK')
end
fprintf('MATLAB IndividualPeak regression checks passed.\n')

end

function expect_error(callback)

failed = false;
try
    callback();
catch
    failed = true;
end
assert(failed,'Expected invalid input to raise an error.')

end
