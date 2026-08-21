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

% Target frequency
targetfrex = 10;

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
