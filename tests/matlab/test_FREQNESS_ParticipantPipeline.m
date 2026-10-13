function tests = test_FREQNESS_ParticipantPipeline
tests = functiontests(localfunctions);
end

function setup(testCase)
test_directory = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(test_directory));
toolbox_root = fullfile(repo_root,'FREQNESS_Toolbox');
test_root = tempname;
mkdir(test_root);

testCase.TestData.original_path = path;
testCase.TestData.figure_visibility = get(groot,'defaultFigureVisible');
testCase.TestData.test_root = test_root;
set(groot,'defaultFigureVisible','off');

% Use an isolated toolbox so tests never read or overwrite study results.
copyfile(fullfile(toolbox_root,'FREQNESS_Startup.m'),test_root);
copyfile(fullfile(toolbox_root,'FREQNESS_Functions'), ...
    fullfile(test_root,'FREQNESS_Functions'));
mkdir(fullfile(test_root,'FREQNESS_ExternalFunctions','nifti_tools'));
copyfile(fullfile(toolbox_root,'FREQNESS_ExternalFunctions','filterFGx.m'), ...
    fullfile(test_root,'FREQNESS_ExternalFunctions'));
mkdir(fullfile(test_root,'FREQNESS_MNI_Coordinates'));
mkdir(fullfile(test_root,'FREQNESS_Data','Dataset_1'));
mkdir(fullfile(test_root,'FREQNESS_Data','Dataset_2'));
mkdir(fullfile(test_root,'FREQNESS_Data','Empty'));

rng(15);
time = (0:511)/64;
data = sin(2*pi*(1:6)'*time) + 0.2*cos(2*pi*(7:12)'*time) + ...
    0.01*randn(6,numel(time));
testData = single(1e-5*data);
save(fullfile(test_root,'FREQNESS_Data','Dataset_1','Participant_a.mat'),'testData');
testData = data + 0.2*randn(size(data));
save(fullfile(test_root,'FREQNESS_Data','Dataset_1','participant_B.mat'),'testData');
testData = 2*data;
save(fullfile(test_root,'FREQNESS_Data','Dataset_2','Participant_a.mat'),'testData');
MNI = [(1:6)' (2:7)' (3:8)'];
save(fullfile(test_root,'FREQNESS_MNI_Coordinates','MNI.mat'),'MNI');

source = fileread(fullfile(toolbox_root,'FREQNESS_MainPipeline.m'));
settings = {
    'frex','[4 8 12]'
    'srate','64'
    'network_duration','4'
    'network_fwidth','[0.8 1.6 2.4]'
    'network_bad_segments','[1 25]'
    'network_rescale','true'
    'save_nifti','false'
    'Landscape.ncomps','5'
    'Patterns.frex','[4 8]'
    'expdk_range2fit','[4 12]'
    'freqgrad_frex2model','[4 12]'
    'compgrad_freq2model','8'
    'lfo_freq','4'
    'backproj_freq2project','8'
    'netrem_freq2remove','8'
    };
for settingi = 1:size(settings,1)
    expression = ['(?m)^' regexptranslate('escape',settings{settingi,1}) '[ \t]*=[^\r\n]*'];
    source = regexprep(source,expression, ...
        [settings{settingi,1} ' = ' settings{settingi,2} ';']);
end

% Check array lifetime at the actual load boundary in the copied pipeline.
source = strrep(source,'            S = load(source_file);',sprintf([ ...
    '            assert(~exist(''FREQ'',''var'') && ~exist(''this_data'',''var'') && ~exist(''S'',''var''));\n' ...
    '            S = load(source_file);']));
script_file = fullfile(test_root,'FREQNESS_MainPipeline.m');
fid = fopen(script_file,'w');
fwrite(fid,source);
fclose(fid);
testCase.TestData.script_file = script_file;

addpath(test_root,'-begin');
clear FREQNESS_Startup
end

function teardown(testCase)
evalin('base','clear');
close all;
set(groot,'defaultFigureVisible',testCase.TestData.figure_visibility);
path(testCase.TestData.original_path);
clear FREQNESS_Startup
rmdir(testCase.TestData.test_root,'s');
end

function testDiscoveryPreservesOrderingAndLegacyStartup(testCase)
[allData,MNI,path_home,folder_data] = FREQNESS_Startup(false);
verifyEqual(testCase,path_home,testCase.TestData.test_root);
verifyTrue(testCase,all(cellfun(@isempty,allData)));
verifyEqual(testCase,{folder_data.name},{'Dataset_1','Dataset_2','Empty'});
verifyEqual(testCase,{folder_data(1).files_mat.name}, ...
    {'Participant_a.mat','participant_B.mat'});
verifyEqual(testCase,folder_data(1).data_size,[6 512]);
verifyEmpty(testCase,folder_data(3).files_mat);
verifySize(testCase,MNI,[6 3]);

[legacyData,legacyMNI] = FREQNESS_Startup();
verifySize(testCase,legacyData{1},[6 512 2]);
verifyClass(testCase,legacyData{1},'double');
verifyEqual(testCase,MNI,legacyMNI);
verifyEmpty(testCase,legacyData{3});
end

function testSavedParticipantsMatchBatchAndRetainFiveComponents(testCase)
run_pipeline(testCase.TestData.script_file);
verifyEqual(testCase,evalin('base','exist(''FREQ'',''var'')'),0);
verifyEqual(testCase,evalin('base','exist(''this_data'',''var'')'),0);
verifyEqual(testCase,evalin('base','exist(''S'',''var'')'),0);
verifyEqual(testCase,evalin('base','exist(''allData'',''var'')'),0);

[allData,~,~,folder_data] = FREQNESS_Startup();
for condi = 1:2
    expected = estimate_reference(allData{condi});
    core_path = fullfile(testCase.TestData.test_root,'FREQNESS_Results', ...
        'CORE',folder_data(condi).name);
    files_core = dir(fullfile(core_path,'*.mat'));
    verifyEqual(testCase,numel(files_core),numel(folder_data(condi).files_mat));
    for subi = 1:numel(files_core)
        saved = load(fullfile(core_path,files_core(subi).name));
        verifyEqual(testCase,size(saved.FREQ.evals,1),6);
        verifyEqual(testCase,size(saved.FREQ.evals,3),1);
        verifyEqual(testCase,size(saved.FREQ.evecs,2),5);
        verifyEqual(testCase,size(saved.FREQ.pats,2),5);
        verifyEqual(testCase,size(saved.FREQ.ts,1),5);
        verifyEqual(testCase,size(saved.FREQ.ts,4),1);
        verifyEqual(testCase,saved.FREQ.duration,4);
        verifyEqual(testCase,saved.PARTICIPANT.condition,folder_data(condi).name);
        verifyEqual(testCase,saved.NETWORK_SETTINGS.ncomps,5);
        verifyEqual(testCase,saved.FREQ_evals,saved.FREQ.evals);
        verifyEqual(testCase,saved.FREQ_pats,saved.FREQ.pats);
        verifyEqual(testCase,saved.FREQ_frex,saved.FREQ.frex);
        verify_participant(testCase,saved.FREQ,expected,subi);
    end
end
end

function testMatchingResultsAreReusedAndChangedSettingsRequireOverwrite(testCase)
run_pipeline(testCase.TestData.script_file);
core_file = fullfile(testCase.TestData.test_root,'FREQNESS_Results', ...
    'CORE','Dataset_1','Participant_a_FREQ.mat');
before = dir(core_file);
run_pipeline(testCase.TestData.script_file);
after = dir(core_file);
verifyEqual(testCase,after.datenum,before.datenum);
verifyEqual(testCase,after.bytes,before.bytes);

set_pipeline_setting(testCase.TestData.script_file,'network_ncomps','2');
try
    run_pipeline(testCase.TestData.script_file);
    verifyFail(testCase,'Changed network settings must not reuse an existing result.');
catch exception
    verifySubstring(testCase,exception.message,'Set overwrite_core to true');
end
saved = load(core_file,'FREQ');
verifyEqual(testCase,size(saved.FREQ.ts,1),5);

set_pipeline_setting(testCase.TestData.script_file,'overwrite_core','true');
run_pipeline(testCase.TestData.script_file);
saved = load(core_file,'FREQ');
verifyEqual(testCase,size(saved.FREQ.ts,1),2);
verifyEqual(testCase,size(saved.FREQ.evals,1),6);
end

function testChangedSourceIsNotSilentlyReused(testCase)
run_pipeline(testCase.TestData.script_file);
source_file = fullfile(testCase.TestData.test_root,'FREQNESS_Data', ...
    'Dataset_1','Participant_a.mat');
S = load(source_file);
testData = double(S.testData)*2;
save(source_file,'testData','-v7.3');
try
    run_pipeline(testCase.TestData.script_file);
    verifyFail(testCase,'Changed participant data must not reuse an existing result.');
catch exception
    verifySubstring(testCase,exception.message,'Set overwrite_core to true');
end
end

function testDiscoveryDefersFiniteValueCheckUntilLoading(testCase)
testData = nan(6,512);
save(fullfile(testCase.TestData.test_root,'FREQNESS_Data', ...
    'Dataset_1','Participant_a.mat'),'testData');
[allData,~,~,folder_data] = FREQNESS_Startup(false);
verifyTrue(testCase,all(cellfun(@isempty,allData)));
verifyEqual(testCase,folder_data(1).data_size,[6 512]);
try
    run_pipeline(testCase.TestData.script_file);
    verifyFail(testCase,'Non-finite participant data must fail when loaded.');
catch exception
    verifySubstring(testCase,exception.message,'real finite voxels-by-time matrix');
end
core_path = fullfile(testCase.TestData.test_root,'FREQNESS_Results','CORE','Dataset_1');
verifyEmpty(testCase,dir(fullfile(core_path,'*.mat')));
end

function testOptionalGroupAnalysesRestoreOriginalResults(testCase)
run_pipeline(testCase.TestData.script_file);
set_pipeline_setting(testCase.TestData.script_file,'run_core','false');
set_pipeline_setting(testCase.TestData.script_file,'run_group_analyses','true');
run_pipeline(testCase.TestData.script_file);

actual = evalin('base','FREQ');
ENTR = evalin('base','ENTR');
backProj = evalin('base','backProj');
dataClean = evalin('base','dataClean');
removedActivity = evalin('base','removedActivity');
FREQ_clean = evalin('base','FREQ_clean');
CFC = evalin('base','CFC');
[allData,~,~,folder_data] = FREQNESS_Startup();
for condi = 1:2
    expected = estimate_reference(allData{condi});
    for subi = 1:numel(folder_data(condi).files_mat)
        one = actual{condi};
        one.evals = one.evals(:,:,subi);
        one.evecs = one.evecs(:,:,:,subi);
        one.pats = one.pats(:,:,:,subi);
        one.ts = one.ts(:,:,:,subi);
        one.scale_factors = one.scale_factors(subi);
        verify_participant(testCase,one,expected,subi);
    end
    expected_ENTR = FREQNESS_EntropyLandscape(expected);
    verifyEqual(testCase,ENTR{condi}.H2,expected_ENTR.H2,'AbsTol',1e-12);
    verifyEqual(testCase,ENTR{condi}.ED,expected_ENTR.ED,'AbsTol',1e-12);
    verifyEqual(testCase,ENTR{condi}.frex,expected_ENTR.frex,'AbsTol',1e-12);
    expected_back = FREQNESS_BackProjection(expected,8,'comps2project',1);
    verifyEqual(testCase,backProj{condi},expected_back,'AbsTol',1e-10);
    verifyEqual(testCase,dataClean{condi}+removedActivity{condi}, ...
        allData{condi}(:,1:256,:),'AbsTol',1e-12);
    expected_clean = estimate_reference(dataClean{condi});
    verifyEqual(testCase,FREQ_clean{condi}.evals,expected_clean.evals,'AbsTol',1e-10);
    expected_cfc = FREQNESS_CrossCoupling(expected,4);
    verifyEqual(testCase,CFC{condi}.PAC_all,expected_cfc.PAC_all,'AbsTol',1e-10);
    verifyEqual(testCase,CFC{condi}.PAC_avg,expected_cfc.PAC_avg,'AbsTol',1e-10);
end
verifyEmpty(testCase,actual{3});
verifyEqual(testCase,evalin('base','exist(''data2remove'',''var'')'),0);
end

function testSelectiveComparisonBlockUsesThinStructures(testCase)

source = fullfile(testCase.TestData.test_root,'FREQNESS_Data', ...
    'Dataset_1','participant_B.mat');
destination = fullfile(testCase.TestData.test_root,'FREQNESS_Data', ...
    'Dataset_2','participant_B.mat');
copyfile(source,destination);
run_pipeline(testCase.TestData.script_file);

set_pipeline_setting(testCase.TestData.script_file,'run_core','false');
set_pipeline_setting(testCase.TestData.script_file,'run_spatial_comparison','true');
set_pipeline_setting(testCase.TestData.script_file,'run_landscape_comparison','true');
set_pipeline_setting(testCase.TestData.script_file,'comparison_conditions', ...
    '{''Dataset_1'',''Dataset_2''}');
set_pipeline_setting(testCase.TestData.script_file,'comparison_design','''within''');
set_pipeline_setting(testCase.TestData.script_file,'comparison_plot_results','false');
set_pipeline_setting(testCase.TestData.script_file,'spatialcomparison_frex','8');
set_pipeline_setting(testCase.TestData.script_file,'spatialcomparison_comps','1');
set_pipeline_setting(testCase.TestData.script_file, ...
    'spatialcomparison_n_permutations','8');
set_pipeline_setting(testCase.TestData.script_file,'landscapecomparison_frex','8');
set_pipeline_setting(testCase.TestData.script_file,'landscapecomparison_comps','1');
run_pipeline(testCase.TestData.script_file);

spatial = evalin('base','STATS_spatial');
landscape = evalin('base','STATS_landscape');
verifyEqual(testCase,spatial.design,'within');
verifyEqual(testCase,spatial.n_A,2);
verifyEqual(testCase,spatial.n_B,2);
verifyEqual(testCase,spatial.frex,8);
verifyEqual(testCase,landscape.design,'within');
verifyEqual(testCase,landscape.evals.frex,8);
verifyEqual(testCase,evalin('base','exist(''FREQ'',''var'')'),0);
verifyEqual(testCase,evalin('base','exist(''comparison_a'',''var'')'),0);
verifyEqual(testCase,evalin('base','exist(''comparison_b'',''var'')'),0);

set_pipeline_setting(testCase.TestData.script_file,'run_group_analyses','true');
set_pipeline_setting(testCase.TestData.script_file,'run_spatial_comparison','false');
set_pipeline_setting(testCase.TestData.script_file,'run_visualizer','false');
set_pipeline_setting(testCase.TestData.script_file,'run_entropy_landscape','true');
set_pipeline_setting(testCase.TestData.script_file,'run_exponential_dk','false');
set_pipeline_setting(testCase.TestData.script_file,'run_freq_gradients','false');
set_pipeline_setting(testCase.TestData.script_file,'run_comp_gradients','false');
set_pipeline_setting(testCase.TestData.script_file,'run_cross_coupling','false');
set_pipeline_setting(testCase.TestData.script_file,'run_backprojection','false');
set_pipeline_setting(testCase.TestData.script_file,'run_network_removal','false');
set_pipeline_setting(testCase.TestData.script_file, ...
    'landscapecomparison_measure','''H2''');
run_pipeline(testCase.TestData.script_file);

landscape = evalin('base','STATS_landscape');
verifyEqual(testCase,landscape.measures,{'H2'});
verifyEqual(testCase,landscape.H2.frex,8);
verifySize(testCase,landscape.H2.values_A,[1 1 2]);
verifySize(testCase,landscape.H2.values_B,[1 1 2]);

end

function run_pipeline(script_file)
command = sprintf('run(''%s'')',strrep(script_file,'''',''''''));
evalin('base',command);
end

function set_pipeline_setting(script_file,name,value)
source = fileread(script_file);
expression = ['(?m)^' regexptranslate('escape',name) '[ \t]*=[^\r\n]*'];
source = regexprep(source,expression,[name ' = ' value ';']);
fid = fopen(script_file,'w');
fwrite(fid,source);
fclose(fid);
end

function FREQ = estimate_reference(data)
FREQ = FREQNESS_NetworkEstimation(data,[4 8 12],64, ...
    'duration',4,'fwidth',[0.8 1.6 2.4],'ncomps',5, ...
    'bad_segments',[1 25],'rescale',true);
end

function verify_participant(testCase,actual,expected,subi)
verifyEqual(testCase,actual.evals,expected.evals(:,:,subi),'AbsTol',1e-10);
verifyEqual(testCase,actual.pats,expected.pats(:,:,:,subi),'AbsTol',1e-10);
verifyEqual(testCase,actual.scale_factors,expected.scale_factors(subi));
verifyEqual(testCase,actual.frex,expected.frex);
verifyEqual(testCase,actual.fwhm,expected.fwhm);
verifyEqual(testCase,actual.bad_segments,expected.bad_segments);
verifyEqual(testCase,actual.regularisation,expected.regularisation);
for frexi = 1:numel(expected.frex)
    for compi = 1:size(actual.evecs,2)
        if dot(actual.evecs(:,compi,frexi),expected.evecs(:,compi,frexi,subi)) < 0
            actual.evecs(:,compi,frexi) = -actual.evecs(:,compi,frexi);
            actual.ts(compi,:,frexi) = -actual.ts(compi,:,frexi);
        end
    end
end
verifyEqual(testCase,actual.evecs,expected.evecs(:,:,:,subi),'AbsTol',1e-10);
verifyEqual(testCase,actual.ts,expected.ts(:,:,:,subi),'AbsTol',1e-10);
end
