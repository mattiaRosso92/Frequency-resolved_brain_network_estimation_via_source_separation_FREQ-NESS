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

% Intermittency manipulation
nactivations = 5;
activationduration = 1.5; % duration of each activation, in seconds
silenceduration = 1;      % duration of each intervening silence, in seconds

% Time vector
srate = EEG.srate;
Tsim = nactivations*activationduration + ...
       (nactivations-1)*silenceduration;
npnts = Tsim * srate;
tvec = 0:1/srate:Tsim-1/srate;

% Source parameters
nsources = 1;
sourceamp = 10;
sourcephase = 0;

% Network-estimation parameters
ncomps = 10;
fwhm = 2;
regularisation = .01;

% Noise
snr = 20;
noiseamp = sourceamp/snr;

% Initialize sensor space
nchans = length(EEG.chanlocs);
eegData = zeros(nchans,npnts);
chanlabels = {EEG.chanlocs.labels};

% Initialize visualization settings
offvis = 3*sourceamp * (1:nchans)';


%% Generate signals

% Initialize source envelope and activation boundaries
sourceenv = zeros(nsources,npnts);
activationbounds = zeros(nactivations,2);

% Generate equal-duration activations separated by shorter silent periods
for activi = 1:nactivations

    activationstart = (activi-1)*(activationduration+silenceduration);
    activationend = activationstart + activationduration;
    activationbounds(activi,:) = [activationstart activationend];

    activationidx = tvec >= activationstart & tvec < activationend;
    sourceenv(activationidx) = 1;

end

% Generate one continuous carrier gated by the intermittent envelope
sourcecarrier = sin(2*pi*tvec*targetfrex + sourcephase);
sources = sourceamp * sourceenv .* sourcecarrier;

%% Dipole injection (ground-truth)

ndips = size(lf.Gain,3);
mydips = ceil(ndips*rand(nsources,1));

% Assign source signal to dipole
dipsData = sources;

% White noise added independently to each channel after dipole selection
wnoise = noiseamp * 2*(rand(size(eegData))-.5);


%% Sensor data

% Normalize gain for source projection
this_gain = squeeze(-lf.Gain(:,1,mydips));
fwd_weights = this_gain ./ max(abs(this_gain));

% Project source data to scalp electrodes and add sensor noise
eegData = wnoise + fwd_weights * dipsData;


%% Visualize ground-truth simulation

figemptycolor = [245,245,220] / 255;
figemptysize = 10;
figemptyidx = ~ismember(1:ndips,mydips);
figmarksize = 2*figemptysize;

figure(1), clf

% Dipoles map
subplot(2,4,[1 2])
hold on
plot3(lf.GridLoc(mydips,1),lf.GridLoc(mydips,2),lf.GridLoc(mydips,3), ...
      'or','MarkerSize',figmarksize,'MarkerFaceColor','r')
plot3(lf.GridLoc(figemptyidx,1),lf.GridLoc(figemptyidx,2), ...
      lf.GridLoc(figemptyidx,3),'ok','MarkerFacecolor',figemptycolor, ...
      'MarkerSize',figemptysize)
xlabel('X'), ylabel('Y'), zlabel('Z')
grid on
axis square
legend('Intermittent source')
title('Dipoles 3D map')

% Dipole signal
subplot(2,4,[3 4])
plot(tvec,dipsData,'r','LineWidth',1.2)
hold on
plot(tvec,sourceamp*sourceenv,'k--','LineWidth',1.2)
plot(tvec,-sourceamp*sourceenv,'k--','LineWidth',1.2)
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
legend('Source signal','Activation envelope')
title(['Ground truth: ' num2str(nactivations) ' intermittent activations'])

% Dipole projection map
subplot(2,4,[5 6])
topoplotIndie(fwd_weights,EEG.chanlocs,'numcontour',0, ...
              'electrodes','off','shading','interp');
title('Source projection')
colormap jet

% Data in EEG sensor space
subplot(2,4,[7 8])
plot(tvec,repmat(offvis,1,npnts)' + eegData')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
yticks(offvis(1:4:nchans))
yticklabels(chanlabels(1:4:nchans))
title('Data in sensor space')

sgtitle('Ground-truth intermittency manipulation')


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

% Compute effective dimensionality from the complete eigenspectrum
[~,ED] = FREQNESS_EntropyLandscape(FREQ);

% Replace the default H2 visualization with the ED landscape
close(gcf)
figure(3), clf
plot(FREQ.frex,ED(:,1),'Color',[.2 .1 .6],'LineWidth',2)
hold on
yline(nactivations,'r--',['Hypothesized N = ' num2str(nactivations)])
xline(targetfrex,'k--','Target frequency')
xlabel('Frequency (Hz)')
ylabel('Effective dimensionality')
ylim([0 1.1*max([ED(:);nactivations])])
title('Effective-dimensionality landscape')
grid on
grid minor


%% --------------------------- %%
%          ASSESSMENT          %
%   -------------------------   %%

% Extract target-frequency FREQ-NESS outputs
targetevals = squeeze(FREQ.evals(:,targetfrexi,1));
targetpat = squeeze(FREQ.pats(:,1,targetfrexi,1));
targetts = squeeze(FREQ.ts(1,:,targetfrexi,1));

% Extract effective dimensionality at the target frequency
targetED = ED(targetfrexi,1);

% Correlation with ground-truth spatial pattern and source time series
FREQcorr_pattern = abs(corr(targetpat,fwd_weights));
FREQcorr_ts_signed = corr(targetts',sources');
FREQcorr_ts = abs(FREQcorr_ts_signed);

% Resolve the arbitrary GED sign and normalize raw time series for visualization
FREQsign = sign(FREQcorr_ts_signed);
if FREQsign == 0
    FREQsign = 1;
end
sourceTS = sources ./ max(abs(sources));
FREQts = FREQsign * targetts ./ max(abs(targetts));


%% Visualize target-frequency FREQ-NESS results

figure(4), clf

% Ground-truth activation pattern
subplot(3,3,1)
topoplotIndie(fwd_weights,EEG.chanlocs,'numcontour',0,'electrodes','off');
title('Ground-truth source pattern')

% First FREQ-NESS activation pattern
subplot(3,3,2)
topoplotIndie(targetpat,EEG.chanlocs,'numcontour',0,'electrodes','off');
title({'FREQ-NESS component #1 pattern', ...
       ['Ground-truth r = ' num2str(FREQcorr_pattern,2)]})
colormap jet

% Target-frequency eigenspectrum
subplot(3,3,3)
bar(1:ncomps,targetevals(1:ncomps),'FaceColor',[.2 .4 .8])
xlabel('Component')
ylabel('Explained variance (%)')
title({['Effective dimensionality = ' num2str(targetED,3)], ...
       ['Hypothesized N = ' num2str(nactivations)]})
axis square

% Ground-truth activation gate
subplot(3,3,[4 5 6])
stairs(tvec,sourceenv,'k','LineWidth',1.5)
ylabel('Gate')
ylim([0 1.1])
title('Ground-truth activation gate')

% Injected dipole oscillation and raw FREQ-NESS component time series
subplot(3,3,[7 8 9])
plot(tvec,sourceTS,'r','LineWidth',1.2)
hold on
plot(tvec,FREQts,'b','LineWidth',1)
xlabel('Time (s)')
ylabel('Normalized amplitude')
ylim([-1.1 1.1])
title({'Injected dipole oscillation (red) and raw FREQ.ts (blue)', ...
       ['Component #1 time-series correlation: r = ' num2str(FREQcorr_ts,2)]})

sgtitle(['FREQ-NESS intermittency test at ' num2str(FREQ.frex(targetfrexi)) ' Hz'])


%% Report assessment

fprintf('\nFREQ-NESS intermittency assessment at %.1f Hz\n',FREQ.frex(targetfrexi));
fprintf('Number of ground-truth activations: %d\n',nactivations);
fprintf('Target-frequency effective dimensionality: %.3f\n',targetED);
fprintf('Component #1 pattern correlation: %.3f\n',FREQcorr_pattern);
fprintf('Raw FREQ.ts source time-series correlation: %.3f\n\n',FREQcorr_ts);
