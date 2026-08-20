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
naway = 1;                % neighbor-rank spacing from the reference source
nparamsteps = 15;         % number of distances in the parametric experiment

% Time vector
srate = EEG.srate;
Tsim = nactivations*activationduration + ...
       (nactivations-1)*silenceduration;
npnts = Tsim * srate;
tvec = 0:1/srate:Tsim-1/srate;

% Source parameters
nsources = nactivations;
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

% Initialize source envelopes and activation boundaries
sourceenv = zeros(nsources,npnts);
activationbounds = zeros(nactivations,2);

% Generate equal-duration activations separated by shorter silent periods
for activi = 1:nactivations

    activationstart = (activi-1)*(activationduration+silenceduration);
    activationend = activationstart + activationduration;
    activationbounds(activi,:) = [activationstart activationend];

    activationidx = tvec >= activationstart & tvec < activationend;
    sourceenv(activi,activationidx) = 1;

end

% Inject the same oscillation into a different dipole at each activation
sourceoscillation = sin(2*pi*tvec*targetfrex + sourcephase);
sources = sourceamp * sourceenv .* sourceoscillation;
combinedEnvelope = sum(sourceenv,1);
combinedSource = sum(sources,1);

%% Dipole injection (ground-truth)

ndips = size(lf.Gain,3);
referencedip = ceil(ndips*rand);

% Rank all dipoles by Euclidean distance from the reference source
dipdistance = sqrt(sum((lf.GridLoc-lf.GridLoc(referencedip,:)).^2,2));
[~,diporder] = sort(dipdistance);
neighboridx = 1 + (0:nactivations-1)*naway;

if neighboridx(end) > ndips
    error('naway is too large for the requested number of activations.')
end

mydips = diporder(neighboridx);
selecteddistance = 1000*dipdistance(mydips); % distance from reference, in mm

% Assign source signal to dipole
dipsData = sources;

% White noise added independently to each channel after dipole selection
wnoise = noiseamp * 2*(rand(size(eegData))-.5);


%% Sensor data

% Normalize gain for source projection
this_gain = squeeze(-lf.Gain(:,1,mydips));
fwd_weights = this_gain ./ max(abs(this_gain),[],1);

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
plot3(lf.GridLoc(mydips(1),1),lf.GridLoc(mydips(1),2), ...
      lf.GridLoc(mydips(1),3), ...
      'or','MarkerSize',figmarksize,'MarkerFaceColor','r')
plot3(lf.GridLoc(mydips(2:end),1),lf.GridLoc(mydips(2:end),2), ...
      lf.GridLoc(mydips(2:end),3), ...
      'ob','MarkerSize',figmarksize,'MarkerFaceColor','b')
plot3(lf.GridLoc(figemptyidx,1),lf.GridLoc(figemptyidx,2), ...
      lf.GridLoc(figemptyidx,3),'ok','MarkerFacecolor',figemptycolor, ...
      'MarkerSize',figemptysize)
xlabel('X'), ylabel('Y'), zlabel('Z')
grid on
axis square
legend('Reference source','Subsequent sources')
title('Dipoles 3D map')

% Dipole signal
subplot(2,4,[3 4])
plot(tvec,combinedSource,'r','LineWidth',1.2)
hold on
plot(tvec,sourceamp*combinedEnvelope,'k--','LineWidth',1.2)
plot(tvec,-sourceamp*combinedEnvelope,'k--','LineWidth',1.2)
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
legend('Source signal','Activation envelope')
title(['Ground truth: ' num2str(nactivations) ' intermittent activations'])

% Dipole projection map
subplot(2,4,[5 6])
topoplotIndie(fwd_weights(:,1),EEG.chanlocs,'numcontour',0, ...
              'electrodes','off','shading','interp');
title('Reference-source projection')
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

% Correlation with ground-truth spatial patterns and complete source sequence
FREQcorr_patterns = abs(corr(targetpat,fwd_weights));
[FREQcorr_pattern,matchedsource] = max(FREQcorr_patterns);
FREQcorr_ts_signed = corr(targetts',combinedSource');
FREQcorr_ts = abs(FREQcorr_ts_signed);

% Resolve the arbitrary GED sign and normalize raw time series for visualization
FREQsign = sign(FREQcorr_ts_signed);
if FREQsign == 0
    FREQsign = 1;
end
sourceTS = combinedSource ./ max(abs(combinedSource));
FREQts = FREQsign * targetts ./ max(abs(targetts));


%% Visualize target-frequency FREQ-NESS results

figure(4), clf

% Ground-truth activation pattern
subplot(3,3,1)
topoplotIndie(fwd_weights(:,matchedsource),EEG.chanlocs, ...
              'numcontour',0,'electrodes','off');
title(['Matched source #' num2str(matchedsource) ' pattern'])

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
stairs(tvec,combinedEnvelope,'k','LineWidth',1.5)
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


%% Parametric source-distance experiment

% Sample the complete feasible neighbor-count range logarithmically
maxnaway = floor((ndips-1)/(nactivations-1));
nawayrange = [0 unique(round(logspace(0,log10(maxnaway),nparamsteps)))];

paramED = zeros(size(nawayrange));
paramdistance = zeros(size(nawayrange));
parammaxdistance = zeros(size(nawayrange));
parammydips = zeros(nactivations,length(nawayrange));
paramsourcepats = zeros(nchans,nactivations,length(nawayrange));
paramFREQpats = zeros(nchans,ncomps,length(nawayrange));

for nawayi = 1:length(nawayrange)

    this_naway = nawayrange(nawayi);
    this_neighboridx = 1 + (0:nactivations-1)*this_naway;
    this_mydips = diporder(this_neighboridx);

    % Project the same intermittent source sequence through the new dipoles
    this_gain = squeeze(-lf.Gain(:,1,this_mydips));
    this_fwd_weights = this_gain ./ max(abs(this_gain),[],1);
    this_eegData = wnoise + this_fwd_weights * dipsData;

    % Estimate FREQ-NESS and compute effective dimensionality
    this_FREQ = FREQNESS_NetworkEstimation(this_eegData,frex,srate, ...
                                           'duration',Tsim, ...
                                           'fwidth',fwhm, ...
                                           'filter','linear', ...
                                           'regularisation',regularisation, ...
                                           'ncomps',ncomps);
    [~,this_ED] = FREQNESS_EntropyLandscape(this_FREQ);
    close(gcf)

    % Store target-frequency ED and physical source distances
    this_distance = 1000*dipdistance(this_mydips); % mm from reference
    paramED(nawayi) = this_ED(targetfrexi,1);
    paramdistance(nawayi) = mean(this_distance(2:end));
    parammaxdistance(nawayi) = max(this_distance);
    parammydips(:,nawayi) = this_mydips;
    paramsourcepats(:,:,nawayi) = this_fwd_weights;

    % Normalize target-frequency FREQ.pats for topographic comparison
    this_FREQpats = squeeze(this_FREQ.pats(:,:,targetfrexi,1));
    paramFREQpats(:,:,nawayi) = this_FREQpats ./ ...
        max(abs(this_FREQpats),[],1);

end

figure(5), clf

subplot(1,2,1)
plot(1:length(nawayrange),paramED,'o-','Color',[.2 .1 .6], ...
     'MarkerFaceColor',[.2 .1 .6],'LineWidth',1.7)
hold on
yline(nactivations,'r--',['N = ' num2str(nactivations)])
xlim([1 length(nawayrange)])
xticks(1:length(nawayrange))
xticklabels(string(nawayrange))
xtickangle(45)
xlabel('Neighbor-rank spacing from reference source')
ylabel('Effective dimensionality at target frequency')
title('ED by neighbor-count distance')
grid on
grid minor

subplot(1,2,2)
plot(paramdistance,paramED,'o-','Color',[.2 .1 .6], ...
     'MarkerFaceColor',[.2 .1 .6],'LineWidth',1.7)
hold on
yline(nactivations,'r--',['N = ' num2str(nactivations)])
xlabel('Mean distance from reference source (mm)')
ylabel('Effective dimensionality at target frequency')
title('ED by physical source distance')
grid on
grid minor

sgtitle(['Parametric source-distance experiment at ' ...
         num2str(FREQ.frex(targetfrexi)) ' Hz'])


%% Visualize spatial patterns across source distances

% Ground-truth source projection patterns
figure(6), clf
set(gcf,'Position',[50 50 1000 1800])

for nawayi = 1:length(nawayrange)
    for sourcei = 1:nactivations

        subplot(length(nawayrange),nactivations, ...
                (nawayi-1)*nactivations+sourcei)
        topoplotIndie(paramsourcepats(:,sourcei,nawayi),EEG.chanlocs, ...
                      'numcontour',0,'electrodes','off');
        caxis([-1 1])

        if nawayi == 1
            title(['Source #' num2str(sourcei)])
        end
        if sourcei == 1
            text(-.7,0,{['naway = ' num2str(nawayrange(nawayi))], ...
                        ['ED = ' num2str(paramED(nawayi),3)]}, ...
                 'HorizontalAlignment','right','FontSize',7, ...
                 'Clipping','off')
        end

    end
end

colormap jet
sgtitle('Ground-truth spatial patterns across source distances')

% Target-frequency FREQ-NESS spatial activation patterns
figure(7), clf
npatstoplot = min(nsources,ncomps);
set(gcf,'Position',[50 50 1000 1800])

for nawayi = 1:length(nawayrange)
    for compi = 1:npatstoplot

        subplot(length(nawayrange),npatstoplot, ...
                (nawayi-1)*npatstoplot+compi)
        topoplotIndie(paramFREQpats(:,compi,nawayi),EEG.chanlocs, ...
                      'numcontour',0,'electrodes','off');
        caxis([-1 1])

        if nawayi == 1
            title(['FREQ.pats #' num2str(compi)])
        end
        if compi == 1
            text(-.7,0,{['naway = ' num2str(nawayrange(nawayi))], ...
                        ['ED = ' num2str(paramED(nawayi),3)]}, ...
                 'HorizontalAlignment','right','FontSize',7, ...
                 'Clipping','off')
        end

    end
end


colormap jet
sgtitle(['Target-frequency FREQ.pats across source distances at ' ...
         num2str(FREQ.frex(targetfrexi)) ' Hz'])


%% Report assessment

fprintf('\nFREQ-NESS intermittency assessment at %.1f Hz\n',FREQ.frex(targetfrexi));
fprintf('Number of ground-truth activations: %d\n',nactivations);
fprintf('Neighbor-rank spacing: %d\n',naway);
fprintf('Mean source distance from reference: %.2f mm\n', ...
        mean(selecteddistance(2:end)));
fprintf('Target-frequency effective dimensionality: %.3f\n',targetED);
fprintf('Component #1 pattern correlation: %.3f\n',FREQcorr_pattern);
fprintf('Raw FREQ.ts source time-series correlation: %.3f\n\n',FREQcorr_ts);
