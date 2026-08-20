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
transitionguard = .5; % seconds excluded around the switch for assessment

% Source parameters
nsources = 2;
sourceamp = 10*ones(nsources,1);
sourcephase = zeros(nsources,1);

% Network-estimation parameters
ncomps = 10;
fwhm = .5;
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
sourceenv(1,tvec < switchtime) = 1;  % source A: first half
sourceenv(2,tvec >= switchtime) = 1; % source B: second half

% Generate two equal-frequency sources with complementary activity
sources = zeros(nsources,npnts);
for sourcei = 1:nsources

    sourcecarrier = sin(2*pi*tvec*targetfrex + sourcephase(sourcei));
    sources(sourcei,:) = sourceamp(sourcei) * sourceenv(sourcei,:) .* sourcecarrier;

end

% White noise added independently to each channel
wnoise = noiseamp * 2*(rand(size(eegData))-.5);


%% Dipole injection (ground-truth)

ndips = size(lf.Gain,3);

% Calculate normalized projection maps for all candidate dipoles
allfwd = squeeze(-lf.Gain(:,1,:));
allfwdmax = max(abs(allfwd),[],1);
allfwdmax(allfwdmax == 0) = 1;
allfwd = allfwd ./ allfwdmax;

% Select one random dipole and the dipole with the least-correlated map
mydips = zeros(nsources,1);
mydips(1) = ceil(ndips*rand);
spatcorr = abs(corr(allfwd,allfwd(:,mydips(1))));
spatcorr(~isfinite(spatcorr)) = inf;
spatcorr(mydips(1)) = inf;
[~,mydips(2)] = min(spatcorr);

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
xline(switchtime,'k--','Switch')
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
xline(switchtime,'k--')
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

% Match the two target-frequency components to ground truth for assessment
corrmatrix = abs(corr(targetTs',sources'));
comporder = 1:2;
if sum(diag(corrmatrix)) < sum(diag(fliplr(corrmatrix)))
    comporder = [2 1];
end

matchedPats = targetPats(:,comporder);
matchedTs = targetTs(comporder,:);

% Correlation with ground truth: source time series and activation envelope
FREQcorr_ts = diag(abs(corr(matchedTs',sources')));

filtSources = filterFGx(sources,FREQ.srate,targetfrex,FREQ.fwhm(targetfrexi),0);
filtFREQts = filterFGx(matchedTs,FREQ.srate,targetfrex,FREQ.fwhm(targetfrexi),0);
[sourceActivation,FREQactivation] = deal(zeros(nsources,npnts));

for sourcei = 1:nsources

    sourceActivation(sourcei,:) = abs(hilbert(filtSources(sourcei,:)));
    FREQactivation(sourcei,:) = abs(hilbert(filtFREQts(sourcei,:)));

end

FREQcorr_activation = diag(abs(corr(FREQactivation',sourceActivation')));

% Normalize activation time series for visualization
sourceActivation = sourceActivation ./ max(sourceActivation,[],2);
FREQactivation = FREQactivation ./ max(FREQactivation,[],2);

% Assess component activity away from filtering transients
firsthalfidx = tvec >= transitionguard & tvec < switchtime-transitionguard;
secondhalfidx = tvec >= switchtime+transitionguard & tvec < Tsim-transitionguard;

FREQactivity = [ mean(FREQactivation(1,firsthalfidx)) mean(FREQactivation(1,secondhalfidx));
                 mean(FREQactivation(2,firsthalfidx)) mean(FREQactivation(2,secondhalfidx)) ];

FREQselectivity = [ FREQactivity(1,1)/sum(FREQactivity(1,:));
                    FREQactivity(2,2)/sum(FREQactivity(2,:)) ];


%% Visualize target-frequency FREQ-NESS results

figure(3), clf

% FREQ-NESS activation patterns
subplot(3,2,1)
topoplotIndie(matchedPats(:,1),EEG.chanlocs,'numcontour',0,'electrodes','off');
title({'FREQ-NESS component matched to source A', ...
       ['Pattern r = ' num2str(abs(corr(matchedPats(:,1),fwd_weights(:,1))),2)]})
subplot(3,2,2)
topoplotIndie(matchedPats(:,2),EEG.chanlocs,'numcontour',0,'electrodes','off');
title({'FREQ-NESS component matched to source B', ...
       ['Pattern r = ' num2str(abs(corr(matchedPats(:,2),fwd_weights(:,2))),2)]})
colormap jet

% Ground-truth source activation time series
subplot(3,2,[3 4])
plot(tvec,sourceActivation(1,:),'r','LineWidth',1.7)
hold on
plot(tvec,sourceActivation(2,:),'b','LineWidth',1.7)
xline(switchtime,'k--','Switch')
ylim([0 1.1])
xlabel('Time (s)')
ylabel('Normalized amplitude')
legend('Source A','Source B')
title('Ground-truth source activation time series')

% Estimated FREQ-NESS component activation time series
subplot(3,2,[5 6])
plot(tvec,FREQactivation(1,:),'r','LineWidth',1.7)
hold on
plot(tvec,FREQactivation(2,:),'b','LineWidth',1.7)
xline(switchtime,'k--','Switch')
ylim([0 1.1])
xlabel('Time (s)')
ylabel('Normalized amplitude')
legend(['Component #' num2str(comporder(1)) ': r = ' ...
        num2str(FREQcorr_activation(1),2) ', selectivity = ' ...
        num2str(FREQselectivity(1),2)], ...
       ['Component #' num2str(comporder(2)) ': r = ' ...
        num2str(FREQcorr_activation(2),2) ', selectivity = ' ...
        num2str(FREQselectivity(2),2)])
title('Estimated FREQ-NESS component activation time series')

sgtitle(['FREQ-NESS stationarity test at ' num2str(FREQ.frex(targetfrexi)) ' Hz'])


%% Report assessment

fprintf('\nFREQ-NESS stationarity assessment at %.1f Hz\n',FREQ.frex(targetfrexi));
fprintf('Source projection correlation: %.3f\n', ...
        abs(corr(fwd_weights(:,1),fwd_weights(:,2))));
fprintf('Component #%d / Source A time-series correlation: %.3f\n', ...
        comporder(1),FREQcorr_ts(1));
fprintf('Component #%d / Source B time-series correlation: %.3f\n', ...
        comporder(2),FREQcorr_ts(2));
fprintf('Source A component first-half selectivity: %.3f\n',FREQselectivity(1));
fprintf('Source B component second-half selectivity: %.3f\n\n',FREQselectivity(2));
