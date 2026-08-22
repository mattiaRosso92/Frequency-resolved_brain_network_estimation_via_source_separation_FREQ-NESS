clc, clear, close all

% Import structures, dipole forward model and simulation dependencies
path_simulations = fileparts(mfilename('fullpath'));
path_toolbox = fileparts(path_simulations);
path_dependencies = fullfile(path_simulations,'Simulation_dependencies');

if ~isfolder(path_dependencies)
    error(['Simulation_dependencies was not found. Place it inside: ' ...
           path_simulations]);
end

addpath(path_dependencies);
addpath(fullfile(path_toolbox,'FREQNESS_Functions'));
addpath(fullfile(path_toolbox,'FREQNESS_ExternalFunctions'));
load(fullfile(path_dependencies,'emptyEEG.mat'),'EEG','lf');

% Reproducible simulation
rng(1)


%% Simulation settings

% Target frequency and sampled frequency range
targetfrex = 10;
frexrange = 5;
frex = targetfrex-frexrange:1:targetfrex+frexrange;

% Time vector
srate = EEG.srate;
Tsim = 10; % duration of the simulation, in seconds
npnts = Tsim * srate;
tvec = 0:1/srate:Tsim-1/srate;

% Stationarity manipulation
switchtime = Tsim/2;
silenceduration = 0; % duration of the interval where both sources are silent
silencestart = switchtime-silenceduration/2;
silenceend = switchtime+silenceduration/2;
transitionguard = .5; % seconds excluded around activation boundaries
transitionplotwindow = .25; % seconds displayed around the phase transition

% Source parameters
nsources = 2;
naway = 1; % neighbor-rank spacing from source A; 0 uses the same dipole
sourceamp = 10*ones(nsources,1);
phasedifference = pi/2; % source B phase offset relative to source A
sourcephase = [0; phasedifference];
phaseoffsets = 0:pi/6:2*pi-pi/6; % source-B offsets from 0 to 330 degrees

% Network-estimation parameters
ncomps = 10;
fwhm = 2;
regularisation = .01;

% Noise
snr = 5;
noiseamp = mean(sourceamp)/snr;

% Initialize sensor space
nchans = length(EEG.chanlocs);
eegData = zeros(nchans,npnts);
chanlabels = {EEG.chanlocs.labels};

% Initialize visualization settings
offvis = 3*max(sourceamp) * (1:nchans)';


%% Generate signals

% Initialize source envelopes
sourceenv = zeros(nsources,npnts);
sourceenv(1,tvec < silencestart) = 1; % source A: before silence
sourceenv(2,tvec >= silenceend) = 1;  % source B: after silence

% Generate two equal-frequency sources with complementary activity
sources = zeros(nsources,npnts);
for sourcei = 1:nsources

    sourceoscillation = sin(2*pi*tvec*targetfrex + sourcephase(sourcei));
    sources(sourcei,:) = sourceamp(sourcei) * sourceenv(sourcei,:) .* ...
                         sourceoscillation;

end

% White noise added independently to each channel
wnoise = noiseamp * 2*(rand(size(eegData))-.5);


%% Dipole injection (ground-truth)

ndips = size(lf.Gain,3);

% Rank dipoles by Euclidean distance from the reference source
referencedip = ceil(ndips*rand);
dipdistance = sqrt(sum((lf.GridLoc-lf.GridLoc(referencedip,:)).^2,2));
[~,diporder] = sort(dipdistance);
neighboridx = 1 + (0:nsources-1)*naway;

if neighboridx(end) > ndips
    error('naway is too large for the requested number of sources.')
end

mydips = diporder(neighboridx);
selecteddistance = 1000*dipdistance(mydips); % distance from A, in mm
flag_overlapping = length(unique(mydips)) == 1;

% Assign source signals to dipoles
dipsData = sources;


%% Sensor data

% Normalize gain separately for each source projection
fwd_weights = zeros(nchans,nsources);
for sourcei = 1:nsources

    this_gain = squeeze(-lf.Gain(:,1,mydips(sourcei)));
    fwd_weights(:,sourcei) = this_gain ./ max(abs(this_gain));

end

% Project source data to scalp electrodes and add sensor noise
eegData = wnoise + fwd_weights * dipsData;


%% Visualize ground-truth simulation

figemptycolor = [245,245,220] / 255;
figemptysize = 10;
figemptyidx = ~ismember(1:ndips,mydips);
figcolors = ['r','b'];
figlinesize = [1.7,1.7];
figmarksize = 2*figemptysize*ones(1,nsources);
if flag_overlapping
    figmarksize = [2.8 1.6]*figemptysize;
end

figure(1), clf

% Dipoles map
subplot(2,4,[1 2])
hold on
for dipi = 1:length(mydips)

    plot3(lf.GridLoc(mydips(dipi),1),lf.GridLoc(mydips(dipi),2), ...
          lf.GridLoc(mydips(dipi),3),'o','MarkerSize',figmarksize(dipi), ...
          'MarkerFaceColor',figcolors(dipi))

end
plot3(lf.GridLoc(figemptyidx,1),lf.GridLoc(figemptyidx,2), ...
      lf.GridLoc(figemptyidx,3),'ok','MarkerFacecolor',figemptycolor, ...
      'MarkerSize',figemptysize)
xlabel('X'), ylabel('Y'), zlabel('Z')
grid on
axis square
legend('Source A','Source B')
title('Dipoles 3D map')

% Dipole signals
subplot(2,4,[3 4])
hold on
for sourcei = 1:nsources

    plot(tvec,dipsData(sourcei,:) + 3*sourceamp(sourcei)*(sourcei-1), ...
         'color',figcolors(sourcei),'LineWidth',figlinesize(sourcei))

end
if silenceduration > 0
    xline(silencestart,'k:','Silence starts')
    xline(silenceend,'k--','Source B starts')
else
    xline(switchtime,'k--','Switch')
end
yticks([])
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
legend('Source A','Source B')
title('Ground-truth source signals')

% Dipole projection maps
subplot(2,4,5)
topoplotIndie(fwd_weights(:,1),EEG.chanlocs,'numcontour',0, ...
              'electrodes','off','shading','interp');
title('Source A projection')
subplot(2,4,6)
topoplotIndie(fwd_weights(:,2),EEG.chanlocs,'numcontour',0, ...
              'electrodes','off','shading','interp');
title('Source B projection')
colormap jet

% Data in EEG sensor space
subplot(2,4,[7 8])
plot(tvec,repmat(offvis,1,npnts)' + eegData')
if silenceduration > 0
    xline(silencestart,'k:')
    xline(silenceend,'k--')
else
    xline(switchtime,'k--')
end
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
yticks(offvis(1:4:nchans))
yticklabels(chanlabels(1:4:nchans))
title('Data in sensor space')

sgtitle('Ground-truth stationarity manipulation')


%% -------------------------- %%
%          FREQ-NESS          %
%   ------------------------   %%

% Estimate the frequency-resolved network landscape
FREQ = FREQNESS_NetworkEstimation(eegData,frex,srate, ...
                                  'duration',Tsim, ...
                                  'fwidth',fwhm, ...
                                  'filter','linear', ...
                                  'regularisation',regularisation, ...
                                  'ncomps',ncomps);

% Identify the target-frequency output
[~,targetfrexi] = min(abs(FREQ.frex-targetfrex));

% Compute effective dimensionality of the frequency-resolved eigenspectrum
[~,ED] = FREQNESS_EntropyLandscape(FREQ);
close(gcf)
targetED = ED(targetfrexi,1);
targetevals = squeeze(FREQ.evals(1:2,targetfrexi,1));

% Visualize the top ten components in the network landscape
Landscape.frex = FREQ.frex;
Landscape.ncomps = ncomps;
FREQNESS_Visualizer(FREQ,Landscape,[],'save_nifti',false);


%% --------------------------- %%
%          ASSESSMENT          %
%   -------------------------   %%

% Extract the first two components at the target frequency from FREQ
targetPats = squeeze(FREQ.pats(:,1:2,targetfrexi,1));
targetTs = squeeze(FREQ.ts(1:2,:,targetfrexi,1));
combinedSource = sum(sources,1);

if flag_overlapping

    % The shared dipole provides one spatial ground truth across both epochs
    comporder = 1:2;
    matchedPats = targetPats;
    matchedTs = targetTs;
    FREQcorr_ts_signed = corr(matchedTs',combinedSource');
    FREQcorr_ts = abs(FREQcorr_ts_signed);

else

    % Match two spatially distinct components to their respective sources
    corrmatrix = abs(corr(targetTs',sources'));
    comporder = 1:2;
    if sum(diag(corrmatrix)) < sum(diag(fliplr(corrmatrix)))
        comporder = [2 1];
    end

    matchedPats = targetPats(:,comporder);
    matchedTs = targetTs(comporder,:);
    FREQcorr_ts_signed = diag(corr(matchedTs',sources'));
    FREQcorr_ts = abs(FREQcorr_ts_signed);

end

% Resolve arbitrary GED signs and normalize raw time series for visualization
FREQsign = sign(FREQcorr_ts_signed);
FREQsign(FREQsign == 0) = 1;
sourceTS = sources ./ max(abs(sources),[],2);
combinedSourceTS = combinedSource ./ max(abs(combinedSource));
FREQts = FREQsign .* matchedTs;
FREQts = FREQts ./ max(abs(FREQts),[],2);

% Assess component activity away from activation boundaries
firsthalfidx = tvec >= transitionguard & ...
               tvec < silencestart-transitionguard;
secondhalfidx = tvec >= silenceend+transitionguard & ...
                tvec < Tsim-transitionguard;
silenceidx = tvec >= silencestart & tvec < silenceend;

FREQactivity = [ sqrt(mean(FREQts(1,firsthalfidx).^2)) ...
                 sqrt(mean(FREQts(1,secondhalfidx).^2));
                 sqrt(mean(FREQts(2,firsthalfidx).^2)) ...
                 sqrt(mean(FREQts(2,secondhalfidx).^2)) ];

FREQselectivity = [ FREQactivity(1,1)/sum(FREQactivity(1,:));
                    FREQactivity(2,2)/sum(FREQactivity(2,:)) ];

% Raw component RMS energy during the silent interval, when present
if any(silenceidx)
    FREQsilence = sqrt(mean(FREQts(:,silenceidx).^2,2));
else
    FREQsilence = nan(nsources,1);
end


%% Visualize target-frequency FREQ-NESS results

figure(3), clf

% FREQ-NESS activation patterns
subplot(4,2,1)
topoplotIndie(matchedPats(:,1),EEG.chanlocs,'numcontour',0,'electrodes','off');
if flag_overlapping
    title({'FREQ-NESS component #1', ...
           ['Shared-pattern r = ' ...
            num2str(abs(corr(matchedPats(:,1),fwd_weights(:,1))),2)]})
else
    title({'FREQ-NESS component matched to source A', ...
           ['Pattern r = ' ...
            num2str(abs(corr(matchedPats(:,1),fwd_weights(:,1))),2)]})
end
subplot(4,2,2)
topoplotIndie(matchedPats(:,2),EEG.chanlocs,'numcontour',0,'electrodes','off');
if flag_overlapping
    title({'FREQ-NESS component #2', ...
           ['Shared-pattern r = ' ...
            num2str(abs(corr(matchedPats(:,2),fwd_weights(:,1))),2)]})
else
    title({'FREQ-NESS component matched to source B', ...
           ['Pattern r = ' ...
            num2str(abs(corr(matchedPats(:,2),fwd_weights(:,2))),2)]})
end
colormap jet

% Ground-truth source activation gates
subplot(4,2,[3 4])
stairs(tvec,sourceenv(1,:),'r','LineWidth',1.5)
hold on
stairs(tvec,sourceenv(2,:),'b','LineWidth',1.5)
if silenceduration > 0
    xline(silencestart,'k:','Silence starts')
    xline(silenceend,'k--','Source B starts')
else
    xline(switchtime,'k--','Switch')
end
ylim([0 1.1])
ylabel('Gate')
title({'Ground-truth source activation gates', ...
       'Source A: red; source B: blue'})

% Source A injected oscillation and raw FREQ-NESS component time series
subplot(4,2,[5 6])
if flag_overlapping
    plot(tvec,combinedSourceTS,'k--','LineWidth',1.2)
else
    plot(tvec,sourceTS(1,:),'k--','LineWidth',1.2)
end
hold on
plot(tvec,FREQts(1,:),'r','LineWidth',1)
if silenceduration > 0
    xline(silencestart,'k:','Silence starts')
    xline(silenceend,'k--','Source B starts')
else
    xline(switchtime,'k--','Switch')
end
ylim([-1.1 1.1])
ylabel('Normalized amplitude')
if flag_overlapping
    title({'Combined phase-reset signal (black) and raw FREQ.ts (red)', ...
           ['Component #1: r = ' num2str(FREQcorr_ts(1),2) ...
            ', ED = ' num2str(targetED,3)]})
else
    title({'Source A injected oscillation (black) and raw FREQ.ts (red)', ...
           ['Component #' num2str(comporder(1)) ': r = ' ...
            num2str(FREQcorr_ts(1),2) ', selectivity = ' ...
            num2str(FREQselectivity(1),2)]})
end

% Source B injected oscillation and raw FREQ-NESS component time series
subplot(4,2,[7 8])
if flag_overlapping
    plot(tvec,combinedSourceTS,'k--','LineWidth',1.2)
else
    plot(tvec,sourceTS(2,:),'k--','LineWidth',1.2)
end
hold on
plot(tvec,FREQts(2,:),'b','LineWidth',1)
if silenceduration > 0
    xline(silencestart,'k:','Silence starts')
    xline(silenceend,'k--','Source B starts')
else
    xline(switchtime,'k--','Switch')
end
ylim([-1.1 1.1])
xlabel('Time (s)')
ylabel('Normalized amplitude')
if flag_overlapping
    title({'Combined phase-reset signal (black) and component #2 (blue)', ...
           ['Residual component correlation: r = ' ...
            num2str(FREQcorr_ts(2),2)]})
else
    title({'Source B injected oscillation (black) and raw FREQ.ts (blue)', ...
           ['Component #' num2str(comporder(2)) ': r = ' ...
            num2str(FREQcorr_ts(2),2) ', selectivity = ' ...
            num2str(FREQselectivity(2),2)]})
end

if silenceduration > 0
    transitionlabel = [num2str(silenceduration) '-s silence'];
else
    transitionlabel = 'immediate transition';
end
sgtitle(['FREQ-NESS stationarity test with ' transitionlabel ' and ' ...
         num2str(rad2deg(phasedifference)) '-deg phase offset at ' ...
         num2str(FREQ.frex(targetfrexi)) ' Hz; naway = ' num2str(naway)])


%% Report assessment

fprintf('\nFREQ-NESS stationarity assessment at %.1f Hz\n',FREQ.frex(targetfrexi));
fprintf('Source B phase offset: %.1f degrees\n',rad2deg(phasedifference));
fprintf('Neighbor-rank spacing: %d\n',naway);
fprintf('Source separation: %.3f mm\n',selecteddistance(2));
fprintf('Source projection correlation: %.3f\n', ...
        abs(corr(fwd_weights(:,1),fwd_weights(:,2))));
fprintf('Effective dimensionality: %.3f\n',targetED);
fprintf('Top two eigenvalues: %.3f%% and %.3f%%\n', ...
        targetevals(1),targetevals(2));
if flag_overlapping
    fprintf('Component #1 / combined signal correlation: %.3f\n', ...
            FREQcorr_ts(1));
    fprintf('Component #2 / combined signal correlation: %.3f\n', ...
            FREQcorr_ts(2));
else
    fprintf('Component #%d / Source A time-series correlation: %.3f\n', ...
            comporder(1),FREQcorr_ts(1));
    fprintf('Component #%d / Source B time-series correlation: %.3f\n', ...
            comporder(2),FREQcorr_ts(2));
    fprintf('Source A component first-half selectivity: %.3f\n', ...
            FREQselectivity(1));
    fprintf('Source B component second-half selectivity: %.3f\n', ...
            FREQselectivity(2));
end
if silenceduration > 0
    fprintf('Component #1 silent-interval RMS: %.3f\n',FREQsilence(1));
    fprintf('Component #2 silent-interval RMS: %.3f\n\n',FREQsilence(2));
else
    fprintf('Transition interval: immediate (no silent samples)\n\n');
end


%% Parametric phase-offset experiment

nphaseoffsets = length(phaseoffsets);
[paramTimeCorrelation,paramSelectivity,paramPatternCorrelation] = ...
    deal(nan(nphaseoffsets,nsources));
paramEvals = nan(nphaseoffsets,nsources);
paramED = nan(nphaseoffsets,1);
paramSubspaceAngles = nan(nphaseoffsets,nsources);
paramCompOrder = nan(nphaseoffsets,nsources);
paramFREQts = nan(nphaseoffsets,nsources,npnts);
paramSourceTs = nan(nphaseoffsets,nsources,npnts);

for phasei = 1:nphaseoffsets

    % Change only the phase of source B
    this_sourcephase = [0; phaseoffsets(phasei)];
    this_sources = zeros(nsources,npnts);
    for sourcei = 1:nsources

        this_sourceoscillation = ...
            sin(2*pi*tvec*targetfrex + this_sourcephase(sourcei));
        this_sources(sourcei,:) = sourceamp(sourcei) * ...
            sourceenv(sourcei,:) .* this_sourceoscillation;

    end

    % Reuse the same dipoles and sensor noise at every phase offset
    this_eegData = wnoise + fwd_weights * this_sources;

    this_FREQ = FREQNESS_NetworkEstimation(this_eegData,frex,srate, ...
                                           'duration',Tsim, ...
                                           'fwidth',fwhm, ...
                                           'filter','linear', ...
                                           'regularisation',regularisation, ...
                                           'ncomps',ncomps);

    [~,this_ED] = FREQNESS_EntropyLandscape(this_FREQ);
    close(gcf)
    paramED(phasei) = this_ED(targetfrexi,1);

    this_targetPats = squeeze(this_FREQ.pats(:,1:2,targetfrexi,1));
    this_targetTs = squeeze(this_FREQ.ts(1:2,:,targetfrexi,1));

    if flag_overlapping

        % Rank-one ground truth: assess both eigenvalue-ordered components
        this_comporder = 1:2;
        this_matchedPats = this_targetPats;
        this_matchedTs = this_targetTs;
        this_combinedSource = sum(this_sources,1);
        this_corrsigned = corr(this_matchedTs',this_combinedSource');
        paramTimeCorrelation(phasei,:) = abs(this_corrsigned)';

    else

        % Match components to two spatially distinct source signals
        this_corrmatrix = abs(corr(this_targetTs',this_sources'));
        this_comporder = 1:2;
        if sum(diag(this_corrmatrix)) < ...
                sum(diag(fliplr(this_corrmatrix)))
            this_comporder = [2 1];
        end
        this_matchedPats = this_targetPats(:,this_comporder);
        this_matchedTs = this_targetTs(this_comporder,:);
        this_corrsigned = diag(corr(this_matchedTs',this_sources'));
        paramTimeCorrelation(phasei,:) = abs(this_corrsigned)';

    end
    paramCompOrder(phasei,:) = this_comporder;

    % Resolve component signs and compute epoch selectivity
    this_sign = sign(this_corrsigned);
    this_sign(this_sign == 0) = 1;
    this_FREQts = this_sign .* this_matchedTs;
    this_FREQtsmax = max(abs(this_FREQts),[],2);
    this_FREQtsmax(this_FREQtsmax == 0) = 1;
    this_FREQts = this_FREQts ./ this_FREQtsmax;
    paramFREQts(phasei,:,:) = this_FREQts;

    if flag_overlapping
        this_sourceTSmax = max(abs(this_combinedSource));
        paramSourceTs(phasei,1,:) = ...
            this_combinedSource ./ this_sourceTSmax;
        paramSourceTs(phasei,2,:) = ...
            this_combinedSource ./ this_sourceTSmax;
    else
        this_sourceTSmax = max(abs(this_sources),[],2);
        this_sourceTSmax(this_sourceTSmax == 0) = 1;
        paramSourceTs(phasei,:,:) = ...
            this_sources ./ this_sourceTSmax;
    end

    this_activity = [sqrt(mean(this_FREQts(1,firsthalfidx).^2)) ...
                     sqrt(mean(this_FREQts(1,secondhalfidx).^2));
                     sqrt(mean(this_FREQts(2,firsthalfidx).^2)) ...
                     sqrt(mean(this_FREQts(2,secondhalfidx).^2))];
    paramSelectivity(phasei,:) = ...
        [this_activity(1,1)/sum(this_activity(1,:)) ...
         this_activity(2,2)/sum(this_activity(2,:))];

    % Spatial-pattern and ground-truth subspace recovery
    if flag_overlapping
        paramPatternCorrelation(phasei,:) = ...
            abs(corr(this_matchedPats,fwd_weights(:,1)))';
    else
        paramPatternCorrelation(phasei,:) = ...
            abs(diag(corr(this_matchedPats,fwd_weights)))';
    end
    groundtruthSpatialBasis = orth(fwd_weights);
    estimatedSpatialBasis = ...
        orth(this_targetPats(:,1:size(groundtruthSpatialBasis,2)));
    this_subspacesingular = ...
        svd(groundtruthSpatialBasis'*estimatedSpatialBasis);
    this_subspacesingular = min(max(this_subspacesingular,-1),1);
    this_subspaceangles = acosd(this_subspacesingular)';
    paramSubspaceAngles(phasei,1:length(this_subspaceangles)) = ...
        this_subspaceangles;

    paramEvals(phasei,:) = ...
        squeeze(this_FREQ.evals(1:2,targetfrexi,1))';

end


%% Visualize phase-dependent component separation

phaseoffsetdegrees = rad2deg(phaseoffsets);

figure(4), clf
set(gcf,'Position',[100 100 1200 800])

if flag_overlapping

subplot(2,2,1)
plot(phaseoffsetdegrees,paramTimeCorrelation(:,1),'-or', ...
     'LineWidth',1.5,'MarkerFaceColor','r')
hold on
plot(phaseoffsetdegrees,paramTimeCorrelation(:,2),'-ob', ...
     'LineWidth',1.5,'MarkerFaceColor','b')
ylabel('Absolute time-series correlation')
title('Combined phase-reset signal recovery')
ylim([0 1.05])
legend('Component #1','Component #2','Location','best')
grid on
xticks(0:30:330)
xlim([0 330])

subplot(2,2,2)
plot(phaseoffsetdegrees,paramED,'-ok','LineWidth',1.5, ...
     'MarkerFaceColor','k')
ylabel('Effective dimensionality (ED)')
title('Dimensionality of the shared-dipole signal')
grid on
xticks(0:30:330)
xlim([0 330])

subplot(2,2,3)
plot(phaseoffsetdegrees,paramPatternCorrelation(:,1),'-or', ...
     'LineWidth',1.5,'MarkerFaceColor','r')
hold on
plot(phaseoffsetdegrees,paramPatternCorrelation(:,2),'-ob', ...
     'LineWidth',1.5,'MarkerFaceColor','b')
xlabel('Source B phase offset (degrees)')
ylabel('Absolute pattern correlation')
title('Recovery of the shared spatial pattern')
ylim([0 1.05])
grid on
xticks(0:30:330)
xlim([0 330])

subplot(2,2,4)
plot(phaseoffsetdegrees,paramEvals(:,1),'-or', ...
     'LineWidth',1.5,'MarkerFaceColor','r')
hold on
plot(phaseoffsetdegrees,paramEvals(:,2),'-ob', ...
     'LineWidth',1.5,'MarkerFaceColor','b')
xlabel('Source B phase offset (degrees)')
ylabel('Explained variance (%)')
title('Leading eigenvalues')
legend('Component #1','Component #2','Location','best')
grid on
xticks(0:30:330)
xlim([0 330])

else

subplot(2,2,1)
plot(phaseoffsetdegrees,paramTimeCorrelation(:,1),'-or', ...
     'LineWidth',1.5,'MarkerFaceColor','r')
hold on
plot(phaseoffsetdegrees,paramTimeCorrelation(:,2),'-ob', ...
     'LineWidth',1.5,'MarkerFaceColor','b')
ylabel('Absolute time-series correlation')
title('Temporal source recovery')
ylim([0 1.05])
legend('Source A','Source B','Location','best')
grid on
xticks(0:30:330)
xlim([0 330])

subplot(2,2,2)
plot(phaseoffsetdegrees,paramSelectivity(:,1),'-or', ...
     'LineWidth',1.5,'MarkerFaceColor','r')
hold on
plot(phaseoffsetdegrees,paramSelectivity(:,2),'-ob', ...
     'LineWidth',1.5,'MarkerFaceColor','b')
ylabel('Epoch selectivity')
title('Temporal separation')
ylim([0 1.05])
grid on
xticks(0:30:330)
xlim([0 330])

subplot(2,2,3)
plot(phaseoffsetdegrees,paramPatternCorrelation(:,1),'-or', ...
     'LineWidth',1.5,'MarkerFaceColor','r')
hold on
plot(phaseoffsetdegrees,paramPatternCorrelation(:,2),'-ob', ...
     'LineWidth',1.5,'MarkerFaceColor','b')
xlabel('Source B phase offset (degrees)')
ylabel('Absolute pattern correlation')
title('Spatial source recovery')
ylim([0 1.05])
grid on
xticks(0:30:330)
xlim([0 330])

subplot(2,2,4)
paramEvalDifference = paramEvals(:,1)-paramEvals(:,2);
evalhandle = plot(phaseoffsetdegrees,paramEvalDifference,'-om', ...
                  'LineWidth',1.5,'MarkerFaceColor','m');
hold on
subspaceangle = max(paramSubspaceAngles,[],2);
subspacehandle = plot(phaseoffsetdegrees,subspaceangle,'-ok', ...
                      'LineWidth',1.5,'MarkerFaceColor','k');
ylabel('Metric value')
ylim([0 1.1*max(paramEvalDifference)])
xlabel('Source B phase offset (degrees)')
title('Eigenvalue split and joint-subspace recovery')
legend([evalhandle subspacehandle], ...
       {'Eigenvalue difference (percentage points)', ...
        'Maximum subspace angle (degrees)'}, ...
       'Location','best')
grid on
xticks(0:30:330)
xlim([0 330])

end

sgtitle(['Parametric FREQ-NESS stationarity test at ' ...
         num2str(FREQ.frex(targetfrexi)) ' Hz; naway = ' num2str(naway)])


%% Visualize raw FREQ.ts for every tested phase offset

figure(5), clf
set(gcf,'Position',[50 50 1600 1200])
phaseplots = tiledlayout(4,3,'TileSpacing','compact','Padding','compact');
if flag_overlapping
    componentoffset = [0; 2.5; 5];
else
    componentoffset = [0; 2.5];
end

for phasei = 1:nphaseoffsets

    nexttile
    hold on
    this_sourceTS = squeeze(paramSourceTs(phasei,:,:));
    this_FREQts = squeeze(paramFREQts(phasei,:,:));

    if flag_overlapping
        sourcehandle = plot(tvec,this_sourceTS(1,:)+componentoffset(1), ...
                            'k--','LineWidth',.8);
        freqAhandle = plot(tvec,this_FREQts(1,:)+componentoffset(2), ...
                           'r','LineWidth',.8);
        freqBhandle = plot(tvec,this_FREQts(2,:)+componentoffset(3), ...
                           'b','LineWidth',.8);
    else
        sourcehandle = plot(tvec,this_sourceTS(1,:)+componentoffset(1), ...
                            'k--','LineWidth',.8);
        plot(tvec,this_sourceTS(2,:)+componentoffset(2), ...
             'k--','LineWidth',.8)
        freqAhandle = plot(tvec,this_FREQts(1,:)+componentoffset(1), ...
                           'r','LineWidth',.8);
        freqBhandle = plot(tvec,this_FREQts(2,:)+componentoffset(2), ...
                           'b','LineWidth',.8);
    end
    xline(switchtime,'k:')

    if flag_overlapping
        xlim([switchtime-transitionplotwindow ...
              switchtime+transitionplotwindow])
        ylim([-1.2 6.2])
        yticks(componentoffset)
        yticklabels({'Injected';'Component 1';'Component 2'})
        title({[num2str(phaseoffsetdegrees(phasei)) '-deg offset'], ...
               ['r = ' num2str(paramTimeCorrelation(phasei,1),'%.2f') ...
                '/' num2str(paramTimeCorrelation(phasei,2),'%.2f') ...
                ', ED = ' num2str(paramED(phasei),'%.2f')]})
    else
        xlim([0 Tsim])
        ylim([-1.2 3.7])
        yticks(componentoffset)
        yticklabels({'Source A';'Source B'})
        title({[num2str(phaseoffsetdegrees(phasei)) '-deg offset'], ...
               ['FREQ components ' num2str(paramCompOrder(phasei,1)) ...
                '/' num2str(paramCompOrder(phasei,2)) ...
                ', selectivity ' ...
                num2str(paramSelectivity(phasei,1),'%.2f') ...
                '/' num2str(paramSelectivity(phasei,2),'%.2f')]})
    end

    if phasei > nphaseoffsets-3
        xlabel('Time (s)')
    end

    if phasei == 1
        if flag_overlapping
            legend([sourcehandle freqAhandle freqBhandle], ...
                   {'Combined injected signal','FREQ.ts component #1', ...
                    'FREQ.ts component #2'},'Location','best')
        else
            legend([sourcehandle freqAhandle freqBhandle], ...
                   {'Injected sources','Matched FREQ.ts A', ...
                    'Matched FREQ.ts B'},'Location','best')
        end
    end

end


title(phaseplots,{['Raw FREQ.ts across source-B phase offsets at ' ...
                   num2str(FREQ.frex(targetfrexi)) ' Hz'], ...
                  ['Black dashed: injected oscillation; vertical line: ' ...
                   'phase transition']})


%% Report parametric phase-offset experiment

fprintf('Parametric source-B phase-offset experiment\n');
if flag_overlapping
    fprintf('ED range across phases: %.3f-%.3f\n', ...
            min(paramED),max(paramED));
    fprintf('Component #1 combined-signal correlation range: %.3f-%.3f\n', ...
            min(paramTimeCorrelation(:,1)), ...
            max(paramTimeCorrelation(:,1)));
    fprintf('Component #2 combined-signal correlation range: %.3f-%.3f\n', ...
            min(paramTimeCorrelation(:,2)), ...
            max(paramTimeCorrelation(:,2)));
    fprintf('Component #1 shared-pattern correlation range: %.3f-%.3f\n\n', ...
            min(paramPatternCorrelation(:,1)), ...
            max(paramPatternCorrelation(:,1)));
else
    [~,bestphasei] = max(mean(paramSelectivity,2));
    [~,worstphasei] = min(mean(paramSelectivity,2));
    fprintf('Best mean selectivity: %.3f at %.1f degrees\n', ...
            mean(paramSelectivity(bestphasei,:)), ...
            phaseoffsetdegrees(bestphasei));
    fprintf('Lowest mean selectivity: %.3f at %.1f degrees\n\n', ...
            mean(paramSelectivity(worstphasei,:)), ...
            phaseoffsetdegrees(worstphasei));
end
