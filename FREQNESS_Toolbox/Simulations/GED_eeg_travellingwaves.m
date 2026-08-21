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
plotduration = .5; % duration displayed in the time-series plots

% Travelling-wave parameters
nsources = 12;
sourceamp = 10;
phaseincrement = pi/6; % phase delay between contiguous sources
sourcephase = (0:nsources-1)' * phaseincrement;
nneighborcandidates = 6; % local neighbors considered at each path step
npathattempts = 500;

% Network-estimation parameters
ncomps = 10;
fwhm = 2;
regularisation = .01;

% Component-gradient parameters
gradientcomps = [1 5];
gradientthreshold = 2;

% Noise
snr = 20;
noiseamp = sourceamp/snr;

% Initialize sensor space
nchans = length(EEG.chanlocs);
eegData = zeros(nchans,npnts);
chanlabels = {EEG.chanlocs.labels};


%% Select a contiguous path through 3D source space

ndips = size(lf.Gain,3);
dipcoords = 1000*lf.GridLoc; % convert source coordinates from metres to mm
pathfound = false;

for attempti = 1:npathattempts

    % Sample a random propagation direction in 3D space
    thisdirection = randn(1,3);
    thisdirection = thisdirection/norm(thisdirection);

    % Begin within the lowest 20% of projections along that direction
    directionprojection = dipcoords*thisdirection';
    [~,projectionorder] = sort(directionprojection);
    startpool = projectionorder(1:max(1,round(.2*ndips)));

    thisdips = zeros(nsources,1);
    thisdips(1) = startpool(randi(length(startpool)));

    % Move through contiguous dipoles in the sampled direction
    for sourcei = 2:nsources

        dipvectors = dipcoords-dipcoords(thisdips(sourcei-1),:);
        dipdistance = sqrt(sum(dipvectors.^2,2));
        dipdistance(thisdips(1:sourcei-1)) = inf;

        [~,neighbororder] = sort(dipdistance);
        neighboridx = neighbororder(1:min(nneighborcandidates, ...
                                           ndips-sourcei+1));

        neighborprogress = dipvectors(neighboridx,:)*thisdirection';
        neighboralignment = neighborprogress ./ dipdistance(neighboridx);
        validneighbors = neighborprogress > 0;

        if ~any(validneighbors)
            break
        end

        neighboridx = neighboridx(validneighbors);
        neighboralignment = neighboralignment(validneighbors);
        [~,bestneighbori] = max(neighboralignment);
        thisdips(sourcei) = neighboridx(bestneighbori);

    end

    if all(thisdips > 0)
        mydips = thisdips;
        wavedirection = thisdirection;
        pathfound = true;
        break
    end

end

if ~pathfound
    error(['A contiguous dipole path could not be generated. Increase ' ...
           'npathattempts or nneighborcandidates.'])
end

% Quantify the selected path geometry
pathvectors = diff(dipcoords(mydips,:),1,1);
pathsteps = sqrt(sum(pathvectors.^2,2));
pathdistance = [0; cumsum(pathsteps)];


%% Generate travelling-wave signals

% Simultaneous oscillations with progressive phase delays along the path
sources = zeros(nsources,npnts);
for sourcei = 1:nsources

    sources(sourcei,:) = sourceamp * ...
        sin(2*pi*tvec*targetfrex-sourcephase(sourcei));

end

% Assign source signals to dipoles
dipsData = sources;

% White noise added independently to each channel
wnoise = noiseamp * 2*(rand(size(eegData))-.5);


%% Sensor data

% Normalize gain separately for each source projection
fwd_weights = zeros(nchans,nsources);
for sourcei = 1:nsources

    this_gain = squeeze(-lf.Gain(:,1,mydips(sourcei)));
    fwd_weights(:,sourcei) = this_gain ./ max(abs(this_gain));

end

% Project the travelling wave to scalp electrodes and add sensor noise
eegData = wnoise + fwd_weights * dipsData;


%% Visualize ground-truth travelling wave and sensor data

figemptycolor = [245,245,220] / 255;
figemptyidx = ~ismember(1:ndips,mydips);
wavecolors = turbo(nsources);
plotidx = tvec <= plotduration;

figure(1), clf
set(gcf,'Position',[100 100 950 800])

% Dipole path in 3D MNI space
hold on
scatter3(dipcoords(figemptyidx,1),dipcoords(figemptyidx,2), ...
         dipcoords(figemptyidx,3),8,figemptycolor,'filled', ...
         'MarkerEdgeColor',[.75 .75 .75])
plot3(dipcoords(mydips,1),dipcoords(mydips,2), ...
      dipcoords(mydips,3),'k-','LineWidth',1.2)
wavehandle = scatter3(dipcoords(mydips,1),dipcoords(mydips,2), ...
                      dipcoords(mydips,3),90,sourcephase,'filled', ...
                      'MarkerEdgeColor','k');

directionlength = max(pathdistance);
directionhandle = quiver3(dipcoords(mydips(1),1), ...
                          dipcoords(mydips(1),2), ...
                          dipcoords(mydips(1),3), ...
                          wavedirection(1)*directionlength, ...
                          wavedirection(2)*directionlength, ...
                          wavedirection(3)*directionlength,0, ...
                          'k','LineWidth',2,'MaxHeadSize',.5);

xlabel('X (mm)'), ylabel('Y (mm)'), zlabel('Z (mm)')
grid on
axis equal
view(45,20)
colormap(turbo)
clim([sourcephase(1) sourcephase(end)])
phasebar = colorbar;
phasebar.Label.String = 'Phase delay (rad)';
legend([wavehandle directionhandle], ...
       {'Contiguous wave sources','Propagation direction'}, ...
       'Location','best')
title('Travelling-wave path through 3D MNI space')

figure(2), clf
set(gcf,'Position',[100 100 1300 850])

% Ground-truth source signals
subplot(2,1,1)
hold on
sourceoffvis = 2.5*sourceamp*(0:nsources-1)';
for sourcei = 1:nsources

    plot(tvec(plotidx),sources(sourcei,plotidx)+sourceoffvis(sourcei), ...
         'Color',wavecolors(sourcei,:),'LineWidth',1.2)

end
xlabel('Time (s)')
ylabel('Source and propagation order')
yticks(sourceoffvis)
yticklabels(compose('Source %d',1:nsources))
title(['Ground-truth oscillations: ' ...
       num2str(rad2deg(phaseincrement)) '-deg phase steps'])

% Data in EEG sensor space
subplot(2,1,2)
channelspacing = 1.2*max(abs(eegData(:)));
offvis = channelspacing*(1:nchans)';
plot(tvec(plotidx),repmat(offvis,1,sum(plotidx))' + ...
     eegData(:,plotidx)')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
yticks(offvis(1:4:nchans))
yticklabels(chanlabels(1:4:nchans))
title('Travelling wave projected into EEG sensor space')

sgtitle(['Ground-truth EEG travelling wave at ' ...
         num2str(targetfrex) ' Hz'])


%% Report simulation geometry

fprintf('\nGround-truth EEG travelling wave at %.1f Hz\n',targetfrex);
fprintf('Number of contiguous sources: %d\n',nsources);
fprintf('Phase increment: %.1f degrees per source\n', ...
        rad2deg(phaseincrement));
fprintf('Random direction: [%.3f %.3f %.3f]\n',wavedirection);
fprintf('Path length: %.1f mm\n',pathdistance(end));
fprintf('Dipole step distance: median %.1f mm, range %.1f-%.1f mm\n\n', ...
        median(pathsteps),min(pathsteps),max(pathsteps));


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

% Compute quadratic entropy and effective dimensionality
[H2,ED] = FREQNESS_EntropyLandscape(FREQ);
fig_entropy = gcf;

% Identify the target-frequency output
[~,targetfrexi] = min(abs(FREQ.frex-targetfrex));
targetevals = squeeze(FREQ.evals(:,targetfrexi,1));
targetPats = squeeze(FREQ.pats(:,1:ncomps,targetfrexi,1));
targetTs = squeeze(FREQ.ts(1:ncomps,:,targetfrexi,1));
targetED = ED(targetfrexi,1);

% Relate each component pattern to the forward model at every wave phase
componentSourceCorrelationSigned = corr(targetPats,fwd_weights);
componentSourceCorrelation = abs(componentSourceCorrelationSigned);
[component2CorrelationSorted,component2SourceOrder] = ...
    sort(componentSourceCorrelation(2,:),'descend');
component2PhaseSorted = sourcephase(component2SourceOrder);
[component2PeakCorrelation,component2PeakSource] = ...
    max(componentSourceCorrelation(2,:));
component2PeakPhase = sourcephase(component2PeakSource);

% Quantify phase coupling within the leading GED component pair
targetTs_filtered = filterFGx(targetTs(1:2,:),srate,targetfrex,fwhm,0);
targetTs_analytic = FREQNESS_AnalyticSignal(targetTs_filtered');
componentphasedifference = angle(targetTs_analytic(:,2) .* ...
                                 conj(targetTs_analytic(:,1)));
componentmeanphase = angle(mean(exp(1i*componentphasedifference)));
componentphasePLV = abs(mean(exp(1i*componentphasedifference)));

% Compare the leading GED pair with the true sensor-space wave subspace
groundtruthWaveBasis = [fwd_weights*cos(sourcephase), ...
                        fwd_weights*sin(sourcephase)];
subspacesingularvalues = svd(orth(groundtruthWaveBasis)' * ...
                             orth(targetPats(:,1:2)));
subspacesingularvalues = min(max(subspacesingularvalues,-1),1);
wavesubspaceangles = acosd(subspacesingularvalues);


%% Visualize FREQ-NESS network landscape

figure(4), clf
set(gcf,'Position',[100 100 1200 500])

subplot(1,2,1)
hold on
componentcolors = turbo(ncomps);
for compi = 1:ncomps

    plot(FREQ.frex,squeeze(FREQ.evals(compi,:,1)), ...
         'Color',componentcolors(compi,:),'LineWidth',1.2)

end
xline(targetfrex,'k--','Target frequency')
xlabel('Frequency (Hz)')
ylabel('Explained variance (%)')
title(['Top ' num2str(ncomps) ' FREQ-NESS eigenvalues'])
grid on

subplot(1,2,2)
plot(FREQ.frex,ED(:,1),'k','LineWidth',2)
hold on
plot(FREQ.frex(targetfrexi),targetED,'or','MarkerFaceColor','r', ...
     'MarkerSize',8)
xline(targetfrex,'k--','Target frequency')
xlabel('Frequency (Hz)')
ylabel('Effective dimensionality (ED)')
title(['ED at target frequency = ' num2str(targetED,'%.2f')])
grid on

sgtitle('FREQ-NESS travelling-wave network landscape')


%% Visualize target-frequency component patterns

figure(5), clf
set(gcf,'Position',[100 100 1300 650])

for compi = 1:ncomps

    subplot(2,ceil(ncomps/2),compi)
    topoplotIndie(targetPats(:,compi),EEG.chanlocs, ...
                  'numcontour',0,'electrodes','off','shading','interp');
    title({['Component #' num2str(compi)], ...
           [num2str(targetevals(compi),'%.2f') '% variance']})

end
colormap jet
sgtitle(['FREQ-NESS patterns at ' num2str(FREQ.frex(targetfrexi)) ' Hz'])


%% Visualize raw target-frequency component time series

figure(6), clf
set(gcf,'Position',[100 100 1300 700])

% Normalize only for stacked visualization; FREQ.ts remains unchanged
targetTsmax = max(abs(targetTs),[],2);
targetTsmax(targetTsmax == 0) = 1;
targetTs_plot = targetTs ./ targetTsmax;
componentoffvis = 2.5*(0:ncomps-1)';
hold on
for compi = 1:ncomps

    plot(tvec(plotidx),targetTs_plot(compi,plotidx) + ...
         componentoffvis(compi),'Color',componentcolors(compi,:), ...
         'LineWidth',1.1)

end
xlabel('Time (s)')
ylabel('Component and eigenvalue order')
yticks(componentoffvis)
yticklabels(compose('Component %d',1:ncomps))
title({['Raw FREQ.ts at ' num2str(FREQ.frex(targetfrexi)) ' Hz'], ...
       ['Components 1-2: phase difference = ' ...
        num2str(rad2deg(componentmeanphase),'%.1f') ...
        ' deg, PLV = ' num2str(componentphasePLV,'%.4f')]})


%% Match component patterns to the dipole forward model

% Use the known forward model to obtain an exploratory source-space match
allfwd = squeeze(-lf.Gain(:,1,:));
allfwdmax = max(abs(allfwd),[],1);
allfwdmax(allfwdmax == 0) = 1;
allfwd = allfwd ./ allfwdmax;

componentdipcorr = abs(corr(allfwd,targetPats));
[componentdipcorr,componentdips] = max(componentdipcorr,[],1);
componentprojection = ...
    (dipcoords(componentdips,:)-dipcoords(mydips(1),:))*wavedirection';
componentordercorr = corr((1:ncomps)',componentprojection,'Type','Spearman');

figure(7), clf
set(gcf,'Position',[100 100 950 800])
hold on
scatter3(dipcoords(figemptyidx,1),dipcoords(figemptyidx,2), ...
         dipcoords(figemptyidx,3),8,figemptycolor,'filled', ...
         'MarkerEdgeColor',[.75 .75 .75])
plot3(dipcoords(mydips,1),dipcoords(mydips,2), ...
      dipcoords(mydips,3),'k-','LineWidth',1.2)
sourcehandle = scatter3(dipcoords(mydips,1),dipcoords(mydips,2), ...
                        dipcoords(mydips,3),65,'ok','LineWidth',1.2);
componenthandle = scatter3(dipcoords(componentdips,1), ...
                           dipcoords(componentdips,2), ...
                           dipcoords(componentdips,3),130,1:ncomps, ...
                           'd','filled','MarkerEdgeColor','k');
xlabel('X (mm)'), ylabel('Y (mm)'), zlabel('Z (mm)')
grid on
axis equal
view(45,20)
colormap(turbo)
componentbar = colorbar;
componentbar.Label.String = 'GED component index';
legend([sourcehandle componenthandle], ...
       {'Ground-truth wave path','Best-matching dipoles'}, ...
       'Location','best')
title({'Exploratory source match of FREQ-NESS patterns', ...
       ['Component-order/direction rho = ' ...
        num2str(componentordercorr,'%.2f')]})


%% Relate the leading component patterns to travelling-wave phase

sourcephasedeg = rad2deg(sourcephase);
firsthalfidx = sourcephase <= pi;

figure(8), clf
set(gcf,'Position',[100 100 1500 500])

% Inspect both leading components over the complete simulated wave
subplot(1,3,1)
hold on
plot(sourcephasedeg,componentSourceCorrelation(1,:),'-o', ...
     'Color',componentcolors(1,:),'LineWidth',1.5, ...
     'MarkerFaceColor',componentcolors(1,:))
plot(sourcephasedeg,componentSourceCorrelation(2,:),'-o', ...
     'Color',componentcolors(2,:),'LineWidth',1.5, ...
     'MarkerFaceColor',componentcolors(2,:))
xline(90,'k--','90 deg')
xline(180,'k:','180 deg')
xlabel('Ground-truth source phase (deg)')
ylabel('|Spatial correlation|')
xticks(sourcephasedeg)
xtickangle(45)
ylim([0 1])
grid on
legend({'Component #1','Component #2'},'Location','best')
title('Forward-pattern match across the full wave')

% Focus on the hypothesized critical phase range from zero to pi
subplot(1,3,2)
hold on
plot(sourcephasedeg(firsthalfidx), ...
     componentSourceCorrelation(2,firsthalfidx),'-o', ...
     'Color',componentcolors(2,:),'LineWidth',2, ...
     'MarkerFaceColor',componentcolors(2,:))
xline(90,'k--','90 deg')
plot(rad2deg(component2PeakPhase),component2PeakCorrelation,'pk', ...
     'MarkerFaceColor','y','MarkerSize',12)
xlabel('Ground-truth source phase (deg)')
ylabel('|Spatial correlation|')
xticks(sourcephasedeg(firsthalfidx))
xlim([sourcephasedeg(1) 180])
ylim([0 1])
grid on
title({['Component #2 over 0-' char(960)], ...
       ['Full-wave peak: ' ...
        num2str(rad2deg(component2PeakPhase),'%.0f') ...
        ' deg, |r| = ' num2str(component2PeakCorrelation,'%.2f')]})

% Rank the phase-associated forward patterns without imposing phase order
subplot(1,3,3)
rankhandle = bar(component2CorrelationSorted,'FaceColor','flat');
rankhandle.CData = wavecolors(component2SourceOrder,:);
xlabel('Descending correlation rank')
ylabel('|Spatial correlation|')
xticks(1:nsources)
xticklabels(compose('%g deg',rad2deg(component2PhaseSorted)))
xtickangle(45)
ylim([0 1])
grid on
title('Component #2 phase ranking')

sgtitle(['Phase association of FREQ-NESS patterns at ' ...
         num2str(FREQ.frex(targetfrexi)) ' Hz'])


%% FREQ-NESS component gradients in EEG sensor space

% Reconstruct approximate 3D sensor coordinates from EEGLAB polar values
channeltheta = deg2rad([EEG.chanlocs.theta]');
channelcolatitude = pi*[EEG.chanlocs.radius]';
channelcoords = [sin(channelcolatitude).*cos(channeltheta), ...
                 sin(channelcolatitude).*sin(channeltheta), ...
                 cos(channelcolatitude)];

[gradCoeff_components,goodFit_components] = ...
    FREQNESS_CompGradients(FREQ,channelcoords, ...
                           'freq2model',targetfrex, ...
                           'comps2model',gradientcomps, ...
                           'threshold_sd',gradientthreshold, ...
                           'plot_all',false);
fig_componentsgradient = gcf;
sgtitle(['Sensor-space FREQ-NESS component gradients at ' ...
         num2str(FREQ.frex(targetfrexi)) ' Hz'])


%% Report FREQ-NESS assessment

fprintf('\nFREQ-NESS travelling-wave assessment at %.1f Hz\n', ...
        FREQ.frex(targetfrexi));
fprintf('Effective dimensionality: %.3f\n',targetED);
fprintf('Top two eigenvalues: %.3f%% and %.3f%%\n', ...
        targetevals(1),targetevals(2));
fprintf('Component #2-#1 phase difference: %.3f degrees (PLV %.4f)\n', ...
        rad2deg(componentmeanphase),componentphasePLV);
fprintf('True/estimated two-pattern subspace angles: %.3f and %.3f degrees\n', ...
        wavesubspaceangles(1),wavesubspaceangles(2));
fprintf('Component #2 strongest source-phase match: %.1f degrees (|r| %.3f)\n', ...
        rad2deg(component2PeakPhase),component2PeakCorrelation);
fprintf('Component #2 descending phase matches (deg): %s\n', ...
        mat2str(rad2deg(component2PhaseSorted)',3));
fprintf('Component #2 descending correlations: %s\n', ...
        mat2str(component2CorrelationSorted,3));
fprintf('Component-order/direction Spearman rho: %.3f\n', ...
        componentordercorr);
fprintf('Sensor-gradient best R^2 [X Y Z]: [%.3f %.3f %.3f]\n\n', ...
        goodFit_components.R2_best(:,1));
