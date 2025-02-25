
%% METHOD FIGURE - EXAMPLE OF TIME-SERIES FROM THE 3559 BRAIN VOXELS (10 SECONDS IN SECOND LINE)

clear
close all
load('/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Mattia/timeseries.mat')


bigT = 10; %seconds to be modulated
% Define parameters
fs = 250;            % Sampling frequency (Hz)
t = 0:1/fs:bigT-1/fs;    % Time vector (1 second)
sbam = 600*[1:6];
bim = transpose(bamba(1:6,1:2500));
sbam(1) = 0;
% Plot the results
figure;
% for i = 1:length(noise_levels)
%subplot(length(noise_levels), 1, i);
plot(t, sbam + bim);
%     title(['Sine wave with noise level ', num2str(noise_levels(i))]);
xlabel('Time (s)');
ylabel('Brain voxels');
set(gcf,'color','w');
set(0, 'DefaultAxesFontName', 'Helvetica Neue')
set(0, 'DefaultTextFontName', 'Helvetica Neue')
%     ylim([-1.5 1.5] * (1 + noise_levels(end))); % Adjust y-limits based on the highest noise level

exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Images/TimeSeries_3559BrainVoxels.pdf'],'Resolution',300)

%% METHOD FIGURE - SIMULATING GED COMPONENTS (10 SECONDS IN SECOND LINE)

clear
close all
bigT = 10; %seconds to be modulated

% Define parameters
fs = 1000;            % Sampling frequency (Hz)
t = 0:1/fs:bigT-1/fs;    % Time vector (1 second)
f = 2.4;              % Frequency of the sine wave (Hz)

% Define noise levels (you can adjust these values to control the noise gradient)
noise_levels = [0.5, 1, 2, 2.5, 5, 10]./8;
signal_level = [5,2,0.5,0.3,0.1,0.05]; %5*noise_levels(end:-1:1);

% Preallocate a matrix to hold the noisy signals
noisy_signals = zeros(length(noise_levels), length(t));

% Generate the noisy signals
for i = 1:length(noise_levels)
    % Generate the clean sine wave
    sine_wave = signal_level(i)*sin(2*pi*f*t);
    % Generate pink noise
    N = length(t);
    white_noise = randn(1, N);  % Generate white noise
    fft_white = fft(white_noise);  % Perform FFT
    frequencies = (1:N) / N;  % Normalized frequency
    pink_filter = 1 ./ sqrt(frequencies);  % 1/f filter
    pink_filter(1) = pink_filter(2);  % Avoid division by zero
    fft_pink = fft_white .* pink_filter;  % Apply the filter
    pink_noise_signal = real(ifft(fft_pink));  % Inverse FFT to get pink noise
    % Normalize pink noise to the desired level
    pink_noise_signal = pink_noise_signal / max(abs(pink_noise_signal));
    noisy_signals(i, :) = sine_wave + noise_levels(i) * randn(size(t)) + pink_noise_signal;
end

sbam = 2*max(signal_level)*[1:6]';
bim = noisy_signals(end:-1:1, :);
sbam(1) = 0;
% Plot the results
figure;
% for i = 1:length(noise_levels)
%subplot(length(noise_levels), 1, i);
plot(t, sbam + bim,'color','k','LineWidth',1.2);
%     title(['Sine wave with noise level ', num2str(noise_levels(i))]);
xlabel('Time (s)');
ylabel('GED components');
set(gcf,'color','w');
set(0, 'DefaultAxesFontName', 'Helvetica Neue')
set(0, 'DefaultTextFontName', 'Helvetica Neue')
%     ylim([-1.5 1.5] * (1 + noise_levels(end))); % Adjust y-limits based on the highest noise level
% end

exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Images/FakeGEDComponentsTimeSeries_10secs.pdf'],'Resolution',300)


%% METHOD FIGURE - SIMULATING GED COMPONENTS (2 SECONDS IN THIRD LINE)

clear
close all
bigT = 2; %seconds to be modulated

% Define parameters
fs = 1000;            % Sampling frequency (Hz)
t = 0:1/fs:bigT-1/fs;    % Time vector (1 second)
f = 2.4;              % Frequency of the sine wave (Hz)

% Define noise levels (you can adjust these values to control the noise gradient)
noise_levels = [0.5, 1, 2, 2.5, 5, 10]./8;
signal_level = [5,2,0.5,0.3,0.1,0.05]; %5*noise_levels(end:-1:1);

addpath('/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED')
% Preallocate a matrix to hold the noisy signals
noisy_signals = zeros(length(noise_levels), length(t));

% Generate the noisy signals
for i = 1:length(noise_levels)
    % Generate the clean sine wave
    sine_wave = signal_level(i)*sin(2*pi*f*t);
%     if i == 1
%         lfo = sine_wave;
%     end
    % Generate pink noise
    N = length(t);
    white_noise = randn(1, N);  % Generate white noise
    fft_white = fft(white_noise);  % Perform FFT
    frequencies = (1:N) / N;  % Normalized frequency
    pink_filter = 1 ./ sqrt(frequencies);  % 1/f filter
    pink_filter(1) = pink_filter(2);  % Avoid division by zero
    fft_pink = fft_white .* pink_filter;  % Apply the filter
    pink_noise_signal = real(ifft(fft_pink));  % Inverse FFT to get pink noise
    % Normalize pink noise to the desired level
    pink_noise_signal = pink_noise_signal / max(abs(pink_noise_signal));
    noisy_signals(i, :) = sine_wave + noise_levels(i) * randn(size(t)) + pink_noise_signal;
end

lfo = noisy_signals(1,:)./signal_level(1);

% instantaneous phase
lfo = filterFGx(lfo,fs,f,0.3,0);

% Step 1: Compute the FFT of the signal
X = fft(lfo);
% Step 2: Create a frequency domain filter
N = length(lfo);
H = zeros(1, N);
H(1) = 1;               % DC component
if mod(N, 2) == 0       % Even length
    H(2:N/2) = 2;
    H(N/2+1) = 1;       % Nyquist frequency
else                    % Odd length
    H(2:(N+1)/2) = 2;
end
% Step 3: Apply the filter to the FFT of the signal
X_H = X .* H;
% Step 4: Compute the inverse FFT to get the analytic signal
x_analytic = ifft(X_H);
% Step 5: Extract the phase of the analytic signal
phase = angle(x_analytic);

% Plot the original signal and its phase
figure;
plot(t, phase,'color','k','LineWidth',1.5);
title('Phase of the target GED component');
xlabel('Time (s)');
ylabel('Phase (radians)');
set(gcf,'color','w');
exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Images/PhaseHilbert.pdf'],'Resolution',300)


sbam = 2*max(signal_level)*[1:6]';
bim = noisy_signals(end:-1:1, :);
sbam(1) = 0;
% Plot the results
figure;
% for i = 1:length(noise_levels)
%subplot(length(noise_levels), 1, i);
plot(t, sbam + bim,'color','k','LineWidth',1.2);
%     title(['Sine wave with noise level ', num2str(noise_levels(i))]);
xlabel('Time (s)');
ylabel('GED components');
set(gcf,'color','w');
set(0, 'DefaultAxesFontName', 'Helvetica Neue')
set(0, 'DefaultTextFontName', 'Helvetica Neue')
%     ylim([-1.5 1.5] * (1 + noise_levels(end))); % Adjust y-limits based on the highest noise level
% end

exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Images/FakeGEDComponentsTimeSeries_2secs.pdf'],'Resolution',300)

%%

%% METHOD FIGURE - CROSS-FREQUENCY COUPLING (THIRD LINE)

%%% OBS!! This section requires the previous section to be run first %%%

close all
% bigT = 1; %seconds to be modulated
idx_carrier = 5;

% Define parameters
bigT = 2; %seconds to be modulated
fs = 1000;            % Sampling frequency (Hz)
t = 0:1/fs:bigT-1/fs;    % Time vector (1 second)
f = 2.4;              % Frequency of the sine wave (Hz)

% Define parameters
signal_freqs = [8 12 22 35 60 97]; %frequencies to be plotted
signal_level = ones(length(signal_freqs),1); %5*noise_levels(end:-1:1);
noise_levels = ones(length(signal_freqs),1)*0.1;


% Preallocate a matrix to hold the noisy signals
noisy_signals = zeros(length(noise_levels), length(t));

% Generate the noisy signals
for i = 1:length(noise_levels)
    % Generate the clean sine wave
    sine_wave = signal_level(i)*sin(2*pi*signal_freqs(i)*t);
    % Generate pink noise
    N = length(t);
    white_noise = randn(1, N);  % Generate white noise
    fft_white = fft(white_noise);  % Perform FFT
    frequencies = (1:N) / N;  % Normalized frequency
    pink_filter = 1 ./ sqrt(frequencies);  % 1/f filter
    pink_filter(1) = pink_filter(2);  % Avoid division by zero
    fft_pink = fft_white .* pink_filter;  % Apply the filter
    pink_noise_signal = real(ifft(fft_pink));  % Inverse FFT to get pink noise
    % Normalize pink noise to the desired level
    pink_noise_signal = pink_noise_signal / (2*max(abs(pink_noise_signal)));
    noisy_signals(i, :) = sine_wave + noise_levels(i) * randn(size(t)) + pink_noise_signal;
end

%ad-hoc modulation
noisy_signals(idx_carrier,:) = noisy_signals(idx_carrier,:) .* (lfo + 1);

sbam = 4*max(signal_level)*[1:6]';
bim = noisy_signals(end:-1:1, :);
sbam(1) = 0;
% Plot the results
figure;
% for i = 1:length(noise_levels)
%subplot(length(noise_levels), 1, i);
plot(t, sbam + bim, 'LineWidth', 1);
%     title(['Sine wave with noise level ', num2str(noise_levels(i))]);
xlabel('Time (s)');
ylabel('Frequencies (Hz)');
ylim([-2 26])
set(gcf,'color','w');
set(0, 'DefaultAxesFontName', 'Helvetica Neue')
set(0, 'DefaultTextFontName', 'Helvetica Neue')
title('1st GED component across frequencies')
exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Images/GEDComponentsWithCrossFrequencyCoupling.pdf'],'Resolution',300)

% Compute power-by-phase modulation
nbins = 37;
clear sbarba
phase_edges = linspace(min(min(phase)),max(max(phase)),nbins+1); %keep soft-coded, in case of re-wrapping
for oo = 1:length(signal_freqs)
    % Step 1: Compute the FFT of the signal
    X = fft(noisy_signals(oo, :));
    % Step 2: Create a frequency domain filter
    N = length(noisy_signals(oo, :));
    H = zeros(1, N);
    H(1) = 1;               % DC component
    if mod(N, 2) == 0       % Even length
        H(2:N/2) = 2;
        H(N/2+1) = 1;       % Nyquist frequency
    else                    % Odd length
        H(2:(N+1)/2) = 2;
    end
    % Step 3: Apply the filter to the FFT of the signal
    X_H = X .* H;
    % Step 4: Compute the inverse FFT to get the analytic signal
    x_analytic = ifft(X_H);
    % Step 5: Extract the phase of the analytic signal
    POWER = abs(x_analytic);
    for i=1:nbins-1
        % Compute modulation over all run (no issues related to moving median)
        sbarba(oo,i) = mean(POWER(phase>phase_edges(i) & phase<phase_edges(i+1)) , 'omitnan');
    end
end

sbam = 1.5*[1:6]';
bim = sbarba(end:-1:1,:);
sbam(1) = -1;
% Plot the results
figure;
% for i = 1:length(noise_levels)
%subplot(length(noise_levels), 1, i);
plot(1:nbins-1, sbam + bim, 'LineWidth', 1.4);
%     title(['Sine wave with noise level ', num2str(noise_levels(i))]);
xlabel('Phase');
ylabel('Frequencies (Hz)');
set(gcf,'color','w');
ylim([-1 11])
xlim([1 nbins-1])
set(0, 'DefaultAxesFontName', 'Helvetica Neue')
set(0, 'DefaultTextFontName', 'Helvetica Neue')
title('Amplitude distribution over phase')
exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Images/SineFittingExample.pdf'],'Resolution',300)


% Parameters
fs = 1000;                 % Sampling frequency
N = 100;                  % Number of samples
f = linspace(1, 100, N);   % Frequency vector from 1 to 100 Hz

% Generate Gaussian signal
center_freq = 60;          % Center frequency of the Gaussian
sigma = 2;                 % Standard deviation (controls the width)

gaussian_signal = 0.4*(exp(-(f - center_freq).^2 / (2 * sigma^2)));

% Add random fluctuations to the Gaussian signal with reduced noise amplitude
fluctuations = 0.05 * randn(size(f));  % Reduced amplitude of fluctuations
noisy_gaussian_signal = gaussian_signal + fluctuations;

% Ensure the signal values are non-negative
noisy_gaussian_signal = max(noisy_gaussian_signal, 0);

% Plot the Gaussian signal in the frequency domain
figure;
plot(f, noisy_gaussian_signal, 'LineWidth', 1.2,'color','k');
xlabel('Frequency (Hz)');
ylabel('Modulation strength');
xlim([1 100]);  % Limit x-axis to the frequency range
set(gcf,'color','w');
set(0, 'DefaultAxesFontName', 'Helvetica Neue')
set(0, 'DefaultTextFontName', 'Helvetica Neue')
title('Phase-amplitude coupling')
exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Images/Phase_AmplitudeModulationOverFrequencies.pdf'],'Resolution',300)

%% METHOD FIGURE - BROADBAND AND NARROW BAND COVARIANCE(S)

clear
close all
load('/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Mattia/cov_broad.mat')

% Define the number of colors
n = 256;
bu = 0.9;
% Create the custom colormap
% First half: dark blue to white
blue_to_white = [linspace(0, 1, n/2)', linspace(0, 1, n/2)', linspace(bu, 1, n/2)'];
% Second half: white to dark red
white_to_red = [linspace(1, bu, n/2)', linspace(1, 0, n/2)', linspace(1, 0, n/2)'];
% Combine the two halves
cmap = [blue_to_white; white_to_red];

% Plot the data
figure; imagesc(1:3559,1:3559,bomba); xlabel('Brain voxels'); ylabel('Brain voxels'); set(gca,'YDir','normal');
colormap(cmap);
% colorbar
set(gcf,'Color','w')
caxis([-20000 20000])
exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Images/BroadBand_Matrix.pdf'],'Resolution',300)

load('/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Mattia/cov_narrow.mat')
figure; imagesc(1:3559,1:3559,bomba); xlabel('Brain voxels'); ylabel('Brain voxels'); set(gca,'YDir','normal');
colormap(cmap);
% colorbar
set(gcf,'Color','w')
caxis([-400 400])
exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Images/NarrowBand_Matrix.pdf'],'Resolution',300)



%%

addpath('/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED')

%% FIGURE 2 - GED EIGENVALUES

col_l = 1; %1 for significant time-windows with different colors; 0 = only grey color
ylimm = [0 0.2]; %the following is the limit for the main figure of the actual data for the paper [2.5 21] %amplitude limits; leave empty [] for automatic adjustment
sbrim = -0.01; %number to adjust the placement of the lines/markers for significance
Bonf_thresh = 0; %1 for Bonferroni correction, 0 for FDR correction
lineplotdum = [20 2];
inter_l = 1; %1 = interpolation of frequencies around power line artifacts; 0 = no interpolation
load_realdata = 3; %1 = real data; 2 = permutation number 1; 3 = permutation number 2

frex = 1:39;
task_l = 1; %1 = resting state; 2 = beat
comps = [1:3]; %1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
workingdir = '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Mattia';

close all

if load_realdata == 1
    load([workingdir '/GEDallsubjs'])
else
    load([workingdir '/evals_randomizations'])
end
load([workingdir '/GED_stats'],'all_frex','stimfrex')
load([workingdir '/FDR_GEDVariance'])
if load_realdata > 1
    load([workingdir '/GEDallsubjs'],'idx_sort')
    Pevals(:,:) = 0.5; %simplistic way to code no indication of significance
end

BUM{1} = 'RS'; BUM{2} = 'beat';
S = [];
S.sbrim = sbrim;
S.ii = 1;
S.conds = {'1','2','3','4','5','6','7','8','9','10'};
% ROIN{1} = lab(ROIs_to_AAL{ii,1}(pp),:);
% load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/ROIs_6.mat'); %loading data
%structure for the function
if load_realdata == 1
    evals_all = evals_all(:,:,idx_sort,:);
elseif load_realdata == 2
    evals_all = evals_all_randLabel(:,:,:,:);
elseif load_realdata == 3
    evals_all = evals_all_randLabelPointWise(:,:,:,:);
end
data = permute(evals_all,[4 3 2 1]);


if inter_l == 1
    %interpolation to remove power line disturbance
    % Set to zero indexes to interpolate (frequencies affected by power line)
    data(:,45:48,:,:) = 0;
    % Linearly interpolate
    xc = 1:size(data,2);
    for ii = 1:size(data,1)
        for jj = 1:length(comps)
            for ss = 1:size(data,3)
                sbarba = squeeze(data(ii,:,ss,jj));
                sbirba = interp1(xc([1:44, 49:end]), sbarba([1:44, 49:end]), 45:48, 'linear');
                data(ii,45:48,ss,jj) = sbirba;
            end
        end
    end
end

if length(frex) > 80
    namee = ['GEDOverFrequencies_Supplementary_Data' num2str(load_realdata)];
else
    namee = ['GEDOverFrequencies_Main_Data' num2str(load_realdata)];
end

if load_realdata > 1 && task_l == 2
    error('only resting state available for randomisation, so set task_l = 1')
end
S.data = data(task_l,frex,:,comps);
%         S.data = data2(:,1:1026,:,:);
S.STE = 2; %1 = dot lines for standard error; 2 = shadows
S.transp = 0.3; %transparency for standard errors shadow
S.time_real = all_frex(frex)';
S.colorline = [
    0.8 0 0;
    0.3686 0.6314 0.7412;
    0.1882 0.4902 0.8118
    0.0784, 0.1569, 0.5255;
    0.0653, 0.1308, 0.4379;
    0.0523, 0.1046, 0.3504;
    0.0392, 0.0785, 0.2627;
    0.0261, 0.0523, 0.1751;
    0.0131, 0.0262, 0.0876;
    0.0000, 0.0000, 0.0000
    ];
S.legendl = 0;
S.x_lim = []; % Set x limits
S.y_lim = ylimm; %Set y limits
S.ROI_n = 1;
S.condition_n = comps;
S.ROIs_labels = BUM(task_l);
S.lineplot = lineplotdum; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    
%         S.subplot = [];


%extracting relevant p-values and doing Bonferroni correction
P_sel = Pevals(comps,frex);
if Bonf_thresh == 1
    Thresh = repmat(5.8140e-04,length(comps),1); %Bonf_thresh
else
    Thresh = PFDR(comps); %FDR thresh
end
vect = zeros(size(P_sel));
vect(P_sel<Thresh) = 1;
vect(:,45:48) = 0;

clear bum2 signtp_col
cnt = 0;
for cc = comps
    % Find the indices where '1's start and end
    d = diff([0, vect(cc,:), 0]); % Take the difference, with padding to detect edges
    starts = find(d == 1); % Find where '1' clusters start
    ends = find(d == -1) - 1; % Find where '1' clusters end
    % Combine the start and end indices into a matrix
    clusters = [starts; ends]';
%     % Display the clusters
%     disp('Clusters of 1s:');
%     disp(clusters);
    FREQ = all_frex(clusters);
    if size(FREQ,2) == 1
        FREQ = FREQ';
    end
    for cl = 1:size(clusters,1)
        cnt = cnt + 1;
        bum2{1,cnt} = FREQ(cl,:)
        signtp_col(cnt) = cc;
    end
end
if ~exist('bum2','var')
    bum2 = []; signtp_col = [];
end
S.signtp = bum2;
S.colorsign = {'*';'*';'*';'*';'*';'*';'*';'*';'*';'*'}; %colors/different symbols for significant clusters (different colors for different effects of the ANOVA)
S.x_lim = [S.time_real(1) S.time_real(end)];
if col_l == 1
    S.signtp_col = signtp_col;
else
    S.signtp_col = [];
end

waveform_plotting_local_v3(S) %actual function

exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Images/Cond' num2str(task_l) '_' namee '.pdf'],'Resolution',300)


%%

addpath('/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED')

%% FIGURE 3 - CROSS-FREQUENCY COUPLING

lfo_l = 1; %low frequency (stimulation frequency 2.48 Hz) modulator (component 1 or 2)


col_l = 1; %1 for significant time-windows with different colors; 0 = only grey color
ylimm = [0 0.12]; %amplitude limits; leave empty [] for automatic adjustment
sbrim = -0.15; %number to adjust the placement of the lines/markers for significance
lineplotdum = [20 2];

conds = [1:2];
export_l = 0; %1 = export images; 0 = not
close all

workingdir = '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Mattia';
load([workingdir '/CFC_by_2.439_comp#' num2str(lfo_l) '.mat'],'sAmpl')
load([workingdir '/lfo' num2str(lfo_l) '.mat'])
load([workingdir '/GEDallsubjs.mat'],'all_frex')
frex = 1:length(all_frex)-6;

S = [];
S.sbrim = sbrim;
S.ii = 1;
S.conds = {'RS','beat'};
% ROIN{1} = lab(ROIs_to_AAL{ii,1}(pp),:);
% load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/ROIs_6.mat'); %loading data
%structure for the function

data = permute(sAmpl,[3 1 2 4]);
data2 = data; data2(:,:,:,2) = data(:,:,:,1); data2(:,:,:,1) = data(:,:,:,2); 
S.data = data2(1,frex,:,:);
%         S.data = data2(:,1:1026,:,:);
S.STE = 2; %1 = dot lines for standard error; 2 = shadows
S.transp = 0.3; %transparency for standard errors shadow
S.time_real = all_frex(frex+6)';
% S.colorline = [0.8 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];
S.colorline = [0.8 0 0; 0.3686 0.6314 0.7412];
if export_l == 1
    S.legendl = 0;
else
    S.legendl = 0;
end
S.x_lim = []; % Set x limits
S.y_lim = ylimm; %Set y limits
S.ROI_n = 1;
S.condition_n = conds;
S.ROIs_labels = {'dum'};
S.lineplot = lineplotdum; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    
%         S.subplot = [];


%extracting relevant p-values and doing Bonferroni correction
P_sel = zeros(1,length(all_frex));
P_sel(7:end) = P(frex,1)';
if lfo_l == 1
    FDR = 0.0145; %comp 2 = 3.65-e04; comp 3 = 0.003
else
    FDR = 0.0211; %comp 2 = 8.53-e05; comp 3 = 0
end
vect = zeros(size(P_sel));
vect(P_sel<FDR) = 1;
vect(1:10) = 0; %non-meaningful comparisons


clear bum2 signtp_col
% Find the indices where '1's start and end
d = diff([0, vect, 0]); % Take the difference, with padding to detect edges
starts = find(d == 1); % Find where '1' clusters start
ends = find(d == -1) - 1; % Find where '1' clusters end
% Combine the start and end indices into a matrix
clusters = [starts; ends]';
%     % Display the clusters
%     disp('Clusters of 1s:');
%     disp(clusters);
FREQ = all_frex(clusters);
if size(FREQ,2) == 1
    FREQ = FREQ';
end
for cl = 1:size(clusters,1)
    bum2{1,cl} = FREQ(cl,:);
    signtp_col(cl) = 1;
end

    
S.signtp = bum2;
S.colorsign = {'*'}; %colors/different symbols for significant clusters (different colors for different effects of the ANOVA)
S.x_lim = [S.time_real(1) S.time_real(end)];
if col_l == 1
    S.signtp_col = signtp_col;
else
    S.signtp_col = [];
end

waveform_plotting_local_v3(S) %actual function


exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_FREQNESS_GED/Images/Cross_Freq_Modulation_GED_2.4Hz_Comp_' num2str(lfo_l) '.pdf'],'Resolution',300)

%%

%%
