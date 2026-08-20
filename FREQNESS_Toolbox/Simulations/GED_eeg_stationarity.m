clc, clear, close all

% Import structures, dipole forward model and simulation dependencies
path_simulations = fileparts(mfilename('fullpath'));
path_dependencies = fullfile(path_simulations,'Simulation_dependencies');

if ~isfolder(path_dependencies)
    error(['Simulation_dependencies was not found. Place it inside: ' ...
           path_simulations]);
end

addpath(path_dependencies);
load(fullfile(path_dependencies,'emptyEEG.mat'),'EEG','lf');

% Reproducible simulation
rng(1)


%% Simulation settings

% Target frequency
stimfrex = 10;

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

    sourcecarrier = sin(2*pi*tvec*stimfrex + sourcephase(sourcei));
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


%% ------------------- %%
%          GED          %
%   -----------------   %%

% Initialization
ncomps = 2;
[GEDmap,GEDts] = deal(zeros(ncomps,nchans),zeros(ncomps,npnts));

% Assign and filter data
fwhm = .5;
broadData = eegData;
narrowData = filterFGx(eegData,srate,stimfrex,fwhm,0);

% Compute covariance matrices over the entire simulation
covS = cov(narrowData');
covR = cov(broadData');

% Visualize covariance matrices
figure(2), clf

subplot(2,3,[1 2])
plot(tvec,repmat(offvis,1,npnts)' + broadData')
xline(switchtime,'k--')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
title('Broadband sensor data')
subplot(2,3,[4 5])
plot(tvec,repmat(offvis,1,npnts)' + narrowData')
xline(switchtime,'k--')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
title('Narrowband sensor data')
subplot(2,3,3)
imagesc(covR), axis square
title('Covariance R (Broadband)')
subplot(2,3,6)
imagesc(covS), axis square
title('Covariance S (Narrowband)')
sgtitle('GED covariance separation')

% GED
[evecs,evals] = eig(covS,covR);
[evals,sidx] = sort(diag(evals),'descend');
evecs = evecs(:,sidx);

% Compute filter forward models and component time series
for compi = 1:ncomps

    GEDmap(compi,:) = evecs(:,compi)'*covS;
    [~,idxmax] = max(abs(GEDmap(compi,:)));
    compsign = sign(GEDmap(compi,idxmax));
    if compsign == 0
        compsign = 1;
    end
    GEDmap(compi,:) = GEDmap(compi,:)*compsign;
    GEDts(compi,:) = compsign*evecs(:,compi)'*eegData;

end


%% --------------------------- %%
%          ASSESSMENT          %
%   -------------------------   %%

% Match GED components to ground truth for assessment and visualization
corrmatrix = abs(corr(GEDts',sources'));
if sum(diag(corrmatrix)) < sum(diag(fliplr(corrmatrix)))

    GEDmap = GEDmap([2 1],:);
    GEDts = GEDts([2 1],:);
    evecs(:,1:2) = evecs(:,[2 1]);

end

% Match signs to corresponding ground-truth source signals
for compi = 1:ncomps

    compsign = sign(corr(GEDts(compi,:)',sources(compi,:)'));
    if compsign == 0
        compsign = 1;
    end
    GEDmap(compi,:) = GEDmap(compi,:)*compsign;
    GEDts(compi,:) = GEDts(compi,:)*compsign;

end

% Correlation with ground truth: source time series and activation envelope
GEDcorr_ts = diag(abs(corr(GEDts',sources')));

filtSources = filterFGx(sources,srate,stimfrex,fwhm,0);
filtGEDts = filterFGx(GEDts,srate,stimfrex,fwhm,0);
[sourceActivation,GEDactivation] = deal(zeros(nsources,npnts));

for sourcei = 1:nsources

    sourceActivation(sourcei,:) = abs(hilbert(filtSources(sourcei,:)));
    GEDactivation(sourcei,:) = abs(hilbert(filtGEDts(sourcei,:)));

end

GEDcorr_activation = diag(abs(corr(GEDactivation',sourceActivation')));

% Normalize activation time series for visualization
sourceActivation = sourceActivation ./ max(sourceActivation,[],2);
GEDactivation = GEDactivation ./ max(GEDactivation,[],2);

% Assess component activity away from filtering transients
firsthalfidx = tvec >= transitionguard & tvec < switchtime-transitionguard;
secondhalfidx = tvec >= switchtime+transitionguard & tvec < Tsim-transitionguard;

GEDactivity = [ mean(GEDactivation(1,firsthalfidx)) mean(GEDactivation(1,secondhalfidx));
                mean(GEDactivation(2,firsthalfidx)) mean(GEDactivation(2,secondhalfidx)) ];

GEDselectivity = [ GEDactivity(1,1)/sum(GEDactivity(1,:));
                   GEDactivity(2,2)/sum(GEDactivity(2,:)) ];


%% Visualize GED results

figure(3), clf

% GED activation patterns
subplot(3,2,1)
topoplotIndie(GEDmap(1,:),EEG.chanlocs,'numcontour',0,'electrodes','off');
title({'GED component #1 map', ...
       ['Source A map r = ' num2str(abs(corr(GEDmap(1,:)',fwd_weights(:,1))),2)]})
subplot(3,2,2)
topoplotIndie(GEDmap(2,:),EEG.chanlocs,'numcontour',0,'electrodes','off');
title({'GED component #2 map', ...
       ['Source B map r = ' num2str(abs(corr(GEDmap(2,:)',fwd_weights(:,2))),2)]})
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

% Estimated GED component activation time series
subplot(3,2,[5 6])
plot(tvec,GEDactivation(1,:),'r','LineWidth',1.7)
hold on
plot(tvec,GEDactivation(2,:),'b','LineWidth',1.7)
xline(switchtime,'k--','Switch')
ylim([0 1.1])
xlabel('Time (s)')
ylabel('Normalized amplitude')
legend(['GED #1: r = ' num2str(GEDcorr_activation(1),2) ...
        ', selectivity = ' num2str(GEDselectivity(1),2)], ...
       ['GED #2: r = ' num2str(GEDcorr_activation(2),2) ...
        ', selectivity = ' num2str(GEDselectivity(2),2)])
title('Estimated GED component activation time series')

sgtitle('GED test of spatial stationarity')


%% Report assessment

fprintf('\nGED stationarity assessment\n');
fprintf('Source projection correlation: %.3f\n', ...
        abs(corr(fwd_weights(:,1),fwd_weights(:,2))));
fprintf('Component #1 / Source A time-series correlation: %.3f\n',GEDcorr_ts(1));
fprintf('Component #2 / Source B time-series correlation: %.3f\n',GEDcorr_ts(2));
fprintf('Component #1 first-half selectivity: %.3f\n',GEDselectivity(1));
fprintf('Component #2 second-half selectivity: %.3f\n\n',GEDselectivity(2));
