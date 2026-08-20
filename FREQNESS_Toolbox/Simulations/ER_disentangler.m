clc, clear, close all

% Import structures and dipole forward model
cd('/Users/au658766/Documents/MIB/learning_Doing/RSS_2022_(Linear Algebra)/all_datafiles');
load('emptyEEG.mat');


%% Simulation settings

% Stimulation frequency
stimfrex = 1.67; 

% Time vector
srate = 1000;
Tsim = 10; % duration of the simulation, in seconds;
npnts = Tsim * srate; 
tvec = 0:1/srate:Tsim-1/srate; % time vector, in seconds

% Define oscillators
noscs = 40;
nevos = 1;
oscshape = {'sin'; 'square'; 'ramp1'; 'ramp2'}; % ... ; test different waveforms

% Signal parameters
snr = 10; % signal-to-noise ratio, vector related to white noise; 5 is excellent, 4 is acceptable

% Oscillations
ratemod = ones(noscs,npnts); % ones for fixed oscillators
oscsamp  = ceil(10*rand(noscs,1)); % you may randomize here 
frexcutoff = 100;
oscsfrex   =  1 + (frexcutoff-1) * rand(noscs,1) ; % for the time being, fixed; flexible for becoming dynamic oscillator   ratemod *..
oscsfrex(1) = stimfrex; % force first oscillation to be the target entrained oscillation
oscsphase = pi*rand(size(oscsfrex)); % randomize phase offsets


% Evoked responses (for the time being, use Gaussian parameters)
evo_fwhm = .1;
evo_s  = evo_fwhm*(2*pi-1)/(4*pi); % normalized width
evo_amp = 2*oscsamp(1); % relataive to target entrained oscillation
nstims = 1.5*round( Tsim/stimfrex ); % how many stimulus repetitions

% Noise
noiseamp = mean( mean(oscsamp) ) / snr;

% Initialize sensor space
nchans = length(EEG.chanlocs);
eegData = zeros(nchans, npnts);
chanlabels = {EEG.chanlocs.labels};

% Initialize visualization settings
offvis = max(oscsamp) * [1:1:nchans]'; %  visualization offset, for channels; scale by amplitude


%% Generate signals

% Initialize signal matrices
oscs = zeros(noscs,npnts); 
evo = zeros(1,npnts); % only one evoked response

% Generate components

% Oscillators
for osci = 1:noscs
    
    oscs(osci,:) = oscsamp(osci) * sin( 2*pi*tvec*oscsfrex(osci) + oscsphase(osci) );
    
end

% Evoked-response
for stimi = 1:nstims % for each event
    
    evo_center = (stimi - 1) / stimfrex; % calculate the center of the Gaussian for the current stimi   
    ercomp1 = evo_amp * ( -exp(-.5*((tvec - evo_center)/evo_s).^2) );  % add the first component of the Gaussian with the current center to the signal
    ercomp2 = evo_amp * ( exp(-.5*((tvec - evo_center+evo_fwhm)/evo_s).^2) ); % shift the second component to get the 'heart-beat' shape
    evo = evo + ercomp1 + ercomp2;
    
end

% Broadband signal 
ed = 5; %exponential decay parameter (the higher it its, the close to white noise
broadgain = .1*max(oscsamp); % set with respect to largest oscillation
pinkSpctr = rand(1,floor(npnts/2)-1) .* exp(-(1:floor(npnts/2)-1)/ed); % white noise + 1/f component
% figure
% subplot(311)
% plot(pinkSpctr)
pinkSpctr = [pinkSpctr(1) pinkSpctr 0 0 pinkSpctr(:,end:-1:1)]'; % 'mirror' spectrum: dc - frequencies - nyquist - negative frequencies
% subplot(312)
% plot(pinkSpctr)
fc = pinkSpctr .* exp(1i*2*pi*rand(size(pinkSpctr)));% fourier coefficients
pinkNoise = broadgain * real(ifft(fc)) * npnts;
% subplot(313)
% plot(pinkNoise)

% white noise (to be added to channels; uncorrelated noise across channels)
wnoise   = noiseamp * 2*(rand(size(eegData))-.5); % 0-1, every component has its own random noise - uncorrelated



%% Dipole injection (ground-truth)
% only inject signals; noise will be added later on

ndips = size(lf.Gain,3); %  number of dipoles
mydips = ceil(ndips*rand(noscs+nevos,1)); % select dipole indexes
dummydip = datasample(setdiff([1:ndips]',mydips), 1); % dummy dipole to project aperiodic 1/f component

% assign signals to dipoles 
dipsData = [evo ; oscs]; 


%% Sensor data
 

% project from any dipole other than mydips; trick to impose spatial correlation assuming a single source for aperiodic component

% normalize gain for projection (not to create abnormal values) - my idea
% to include, might be a bad one
fwd_weights_broad = squeeze( -lf.Gain(:,1,dummydip) ./ max(abs(lf.Gain(:,1,dummydip))) ); % why negative? It was in a solution by Cohen...
% NOTE: it looks like t
broadsig = fwd_weights_broad * pinkNoise'; % project broadband component in spensor space as projected pink noise
% TODO: add alpha/beta bumps


% normalize gain for projection (not to create abnormal values) - my idea
% to include, might be a bad one
fwd_weights_narrow = squeeze( -lf.Gain(:,1,mydips) ./ max(abs(lf.Gain(:,1,mydips))) ); % why negative? It was in a solution by Cohen...

% simulate data
eegData = wnoise + broadsig(:,end-1); 
% now project dipole data to scalp electrodes
eegData = eegData + fwd_weights_narrow * dipsData;


% Visualize
figemptycolor = [245, 245, 220] / 255;
figemptysize = 10;
figemptyidx = ~ismember(1:ndips,(mydips));
figcolors = [ 'r','b',repmat(['k'],1,noscs-1) ];
figlinesize = [ 1.7,1.7,repmat(.8,1,noscs-1) ];
figmarksize = [ 2*figemptysize,2*figemptysize,repmat(1.3*figemptysize,1,noscs-1) ];


figure(1), clf
% Dipoles map
subplot(2,4,[1 2])
hold on % plot the dipoles with signal first
for dipi = 1:length(mydips)
    plot3(lf.GridLoc(mydips(dipi),1),lf.GridLoc(mydips(dipi),2),lf.GridLoc(mydips(dipi),3),'o','MarkerSize',figmarksize(dipi),'MarkerFaceColor',figcolors(dipi))
end
% then the rest
plot3(lf.GridLoc(figemptyidx,1),lf.GridLoc(figemptyidx,2),lf.GridLoc(figemptyidx,3),'ok','MarkerFacecolor',figemptycolor,'MarkerSize',figemptysize)
xlabel('X'), ylabel('Y'), zlabel('Z')
grid on
axis square
legend('Evoked response','Entrained oscillation')
title('Dipoles 3D map')

% Dipole signals (ground-truth)
subplot(2,4,[3 4])
hold on
for dipi = 1:length(mydips)
plot(repmat(offvis(dipi),1,npnts)' + dipsData(dipi,:)', 'color', figcolors(dipi), 'LineWidth', figlinesize(dipi))
end
yticks([])
xlabel('Time (ms)')
ylabel('Amplitude (a.u.)')
legend('Evoked response','Entrained oscillation')
title('Signals') 

% Dipole projection maps (forward model, ground-truth)
subplot(2,4,5)
clim = [-45 45];
topoplotIndie(-lf.Gain(:,1,mydips(2)), EEG.chanlocs,'maplimits',clim,'numcontour',0,'electrodes','off','shading','interp');
title('1st Oscillator projection')
subplot(2,4,6)
topoplotIndie(-lf.Gain(:,1,mydips(1)), EEG.chanlocs,'maplimits',clim,'numcontour',0,'electrodes','off','shading','interp');
title('1st Evoked response projection')
colormap jet

% Data in EEG sensor space
subplot(2,4,[7:8])
plot(repmat(offvis,1,npnts)' + eegData')
xlabel('Time (ms)')
ylabel('Amplitude (a.u.)')
yticks(offvis(1:4:nchans))
yticklabels(chanlabels(nchans:-4:1))
title('Data in sensor space')






%% SEE LAN_FRIDAY_PRACT_1_SOL FOR STRUCTURE OF PLOTTING FOR COMPATING SOURCE SEPARATIONS

% Initialization
ncomps = 2; % number of target  components
[GEDmap,PCAmap,ICAmap] = deal( zeros(ncomps,nchans) );
[GEDts,PCAts,ICAts] = deal( zeros(ncomps,npnts) );


%% ------------------- %%
%          GED          %
%   -----------------   %

% Assing and filter data
fwhm = .3; % filter width
broadData = eegData;
narrowData = filterFGx(eegData,srate,stimfrex,fwhm,0);

% Compute covariance matrices
covS = cov(narrowData');
covR = cov(broadData');

% Visualize
figure(2), clf

subplot(2,3,[1:2])
plot(repmat(offvis,1,npnts)' + broadData')
ylim([min(broadData(:)) max(broadData(:))+max(offvis(:))])
xlabel('Time (ms)')
ylabel('Amplitude (a.u.)')
title('Broadband sensor data')
subplot(2,3,[4:5])
plot(repmat(offvis,1,npnts)' + narrowData')
ylim([min(broadData(:)) max(broadData(:))+max(offvis(:))])
xlabel('Time (ms)')
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
[evals,sidx]  = sort( diag(evals),'descend' );
evecs = evecs(:,sidx);

% compute filter forward model and flip sign

for compi = 1:ncomps
GEDmap(compi,:) = evecs(:,compi)'*covS; % get component
[~,idxmax] = max(abs(GEDmap(:,1)));  % find max magnitude
GEDmap(compi,:)  = GEDmap(compi,:)*sign(GEDmap(compi,idxmax)); % possible sign flip

% compute first 2 components time series (projections)
GEDts(compi,:) = evecs(:,compi)'*eegData;
end


% %% ------------------- %%
% %          PCA          %
% %   -----------------   %
% 
% Compute covariance matrices
covS = cov(broadData');

% PCA
[evecs,evals] = eig(covS);
[evals,sidx]  = sort( diag(evals),'descend' );
evecs = evecs(:,sidx);

% compute filter forward model and flip sign
for compi = 1:ncomps
PCAmap(compi,:) = evecs(:,compi)'*covS; % get component
[~,idxmax] = max(abs(PCAmap(:,1)));  % find max magnitude
% idx=31; ??
PCAmap(compi,:)  = PCAmap(compi,:)*sign(PCAmap(compi,idxmax)); % possible sign flip

% compute first 2 components time series (projections)
PCAts(compi,:) = evecs(:,compi)'*eegData;
end



% ------------------- %%
%          ICA          %
%   -----------------   %


% reshape to 2D and push through jade
ivecs  = jader(eegData,30); % second argument limits the N of sources to extract

% compute ICA map and time series
allmaps = pinv(ivecs'); %compute in one step
for compi = 1:30
    ICAmap(compi,:) = allmaps(compi,:);
    ICAts(compi,:)  = ivecs(compi,:)*eegData;
end




%% TF analysis
      
% frequencies in Hz
frex = linspace(2,20,20);

% convenient to have component time series data as 2D
allts = [ oscs(1,:) % groundtruth (dipole oscillator)
          evo  
          GEDts(1,:)
          GEDts(2,:)
          PCAts(1,:)
          PCAts(2,:)
          ICAts(1,:)
          ICAts(2,:)];
         
nts = size(allts,1);

% initialize time-frequency matrix
ctf = zeros(nts,length(frex),npnts);


% loop over frequencies
for fi=1:length(frex)
    
    % filter data for both components at this frequency
    filtcomp = filterFGx(allts,EEG.srate,frex(fi),4);
    
    % loop over components
    for compi=1:size(ctf,1)
        
        % compute power time series as envelope of Hilbert transform
        as = hilbert(filtcomp(compi,:));
        
        % TF power is trial-average power
        ctf(compi,fi,:) = mean( abs(as).^2 ,2);
    end
end



%% --------------------------- %%
%          ASSESSMENT        %
%   -------------------------   %


% Correlation with groundtruth: activation timeseries 
% 'at best': pick maximum, and absolute value 
GEDcorr_osc = max( abs(corr(GEDts(1,:)',oscs(1,:)')) , abs(corr(GEDts(1,:)',evo'))  );
GEDcorr_evo = max( abs(corr(GEDts(2,:)',oscs(1,:)')) , abs(corr(GEDts(2,:)',evo'))  );
PCAcorr_osc = max( abs(corr(PCAts(1,:)',oscs(1,:)')) , abs(corr(PCAts(1,:)',evo'))  );
PCAcorr_evo = max( abs(corr(PCAts(2,:)',oscs(1,:)')) , abs(corr(PCAts(2,:)',evo'))  );
ICAcorr_osc = max( abs(corr(ICAts(1,:)',oscs(1,:)')) , abs(corr(ICAts(1,:)',evo'))  );
ICAcorr_evo = max( abs(corr(ICAts(2,:)',oscs(1,:)')) , abs(corr(ICAts(2,:)',evo'))  );
% conctenate (for plots)
corrs_ts = [ nan;
          nan;
          GEDcorr_osc;
          GEDcorr_evo;
          PCAcorr_osc;
          PCAcorr_evo;
          ICAcorr_osc;
          ICAcorr_evo; ];
      
% Correlation with groundtruth: activation pattern 
GEDcorr_osc_map = max( abs(corr(GEDmap(1,:)',fwd_weights_narrow(:,1))) , abs(corr(GEDmap(1,:)',fwd_weights_narrow(:,2)))  );
GEDcorr_evo_map = max( abs(corr(GEDmap(2,:)',fwd_weights_narrow(:,1))) , abs(corr(GEDmap(2,:)',fwd_weights_narrow(:,2)))  );
PCAcorr_osc_map = max( abs(corr(PCAmap(1,:)',fwd_weights_narrow(:,1))) , abs(corr(PCAmap(1,:)',fwd_weights_narrow(:,2)))  );
PCAcorr_evo_map = max( abs(corr(PCAmap(2,:)',fwd_weights_narrow(:,1))) , abs(corr(PCAmap(2,:)',fwd_weights_narrow(:,2)))  );
ICAcorr_osc_map = max( abs(corr(ICAmap(1,:)',fwd_weights_narrow(:,1))) , abs(corr(ICAmap(1,:)',fwd_weights_narrow(:,2)))  );
ICAcorr_evo_map = max( abs(corr(ICAmap(2,:)',fwd_weights_narrow(:,1))) , abs(corr(ICAmap(2,:)',fwd_weights_narrow(:,2)))  );
% conctenate (for comparison plots)
corrs_map = [ nan;
          nan;
          GEDcorr_osc_map;
          GEDcorr_evo_map;
          PCAcorr_osc_map;
          PCAcorr_evo_map;
          ICAcorr_osc_map;
          ICAcorr_evo_map; ];



%% Visualize


% GED: show activation timeseries, topographical maps and eigenspectrum
figure(3)
subplot(3,3,[1:3])
plot(GEDts(1,:))
title('Timeseries component #1')
hold on
plot(GEDts(2,:))
title('Timeseries component #2')
legend(['Component #1 - r = ' num2str(GEDcorr_osc)] , ['Component #2 - r = ' num2str(GEDcorr_evo)])
subplot(3,3,4)
topoplotIndie(GEDmap(1,:),EEG.chanlocs,'numcontour',0,'electrodes','off');
title( { 'Map component #1' , ['r = ' num2str(GEDcorr_osc_map)] } )
subplot(3,3,5)
topoplotIndie(GEDmap(2,:),EEG.chanlocs,'numcontour',0,'electrodes','off');
title( { 'Map component #2' , ['r = ' num2str(GEDcorr_evo_map)] } )
subplot(3,3,6)
plot(evals,'s-','linew',2,'markersize',10,'markerfacecolor','k')
set(gca,'xlim',[0 15])
ylabel('\lambda'), title('GED eigenspectra'), axis square
sgtitle('GED results')


% PCA: show activation timeseries, topographical maps and eigenspectrum
figure(4)
subplot(3,3,[1:3])
plot(PCAts(1,:))
title('Timeseries component #1')
hold on
plot(PCAts(2,:))
title('Timeseries component #2')
legend(['Component #1 - r = ' num2str(PCAcorr_osc)] , ['Component #2 - r = ' num2str(PCAcorr_evo)])
subplot(3,3,4)
topoplotIndie(PCAmap(1,:),EEG.chanlocs,'numcontour',0,'electrodes','off');
title( { 'Map component #1' , ['r = ' num2str(PCAcorr_osc_map)] } )
subplot(3,3,5)
topoplotIndie(PCAmap(2,:),EEG.chanlocs,'numcontour',0,'electrodes','off');
title( { 'Map component #2' , ['r = ' num2str(PCAcorr_evo_map)] } )
subplot(3,3,6)
plot(evals,'s-','linew',2,'markersize',10,'markerfacecolor','k')
set(gca,'xlim',[0 15])
ylabel('\lambda'), title('PCA eigenspectra'), axis square
sgtitle('PCA results')

% ICA: show activation timeseries, topographical maps and eigenspectrum
figure(5)
subplot(3,3,[1:3])
plot(ICAts(1,:))
title('Timeseries component #1')
hold on
plot(ICAts(2,:))
legend(['Component #1 - r = ' num2str(ICAcorr_osc)] , ['Component #2 - r = ' num2str(ICAcorr_evo)])
title('Timeseries component #2')
subplot(3,3,4)
topoplotIndie(ICAmap(1,:),EEG.chanlocs,'numcontour',0,'electrodes','off');
title( { 'Map component #1' , ['r = ' num2str(ICAcorr_osc_map)] } )
subplot(3,3,5)
topoplotIndie(ICAmap(2,:),EEG.chanlocs,'numcontour',0,'electrodes','off');
title( { 'Map component #2' , ['r = ' num2str(ICAcorr_evo_map)] } )
subplot(3,3,6)
% NOTE: ICA doesn't have an eigenspectrum, so we take 
% the variance of the maps as a proxy measure.
plot(var(ICAmap,[],2),'s-','linew',2,'markersize',10,'markerfacecolor','k')
set(gca,'xlim',[0 15])
ylabel('ICA energy'), title('ICA spectrum'), axis square
sgtitle('ICA results')


% TF comparison
figure(6), clf
decompName = { 'Ground-truth oscillation','Groundtruth evoked response','GED #1','GED #2','PCA #1','PCA #2','ICA #1','ICA #2' };
for tsi = nts:-1:1
    % Calculate the actual subplot position, skipping 5 and 10
    subplotPosition = tsi + floor((tsi-1)/4);

    % Skip positions 5 and 10
    if subplotPosition ~= 5 && subplotPosition ~= 10
        subplot(2,5,subplotPosition)
        contourf(tvec, frex, squeeze(ctf(tsi,:,:)), 40, 'linecolor', 'none')
        axis square
        title({decompName{tsi}, ['r = ' num2str(corrs_ts(tsi))]})
        xlabel('Time (s)'), ylabel('Frequency (Hz)')
    end
end






