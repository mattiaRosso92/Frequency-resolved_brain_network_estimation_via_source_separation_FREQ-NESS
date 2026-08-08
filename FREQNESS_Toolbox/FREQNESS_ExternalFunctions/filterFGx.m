function [filtdat,empVals,fx] = filterFGx(data,srate,f,fwhm,showplot)
% filterFGx   Narrow-band filter via frequency-domain Gaussian
%  [filtdat,empVals] = filterFGx(data,srate,f,fwhm,showplot)
% 
% 
%    INPUTS
%       data : 1 X time or chans X time
%      srate : sampling rate in Hz
%          f : peak frequency of filter
%       fhwm : standard deviation of filter, 
%              defined as full-width at half-maximum in Hz
%   showplot : set to true to show the frequency-domain filter shape
% 
%    OUTPUTS
%    filtdat : filtered data
%    empVals : the empirical frequency and FWHM (in Hz and in ms)
% 
% Empirical frequency and FWHM depend on the sampling rate and the
% number of time points, and may thus be slightly different from
% the requested values.
% 
% mikexcohen@gmail.com

%% input check

if nargin < 4
    help filterFGx
    error('Not enough inputs')
end

if ~isnumeric(data) || ~isreal(data) || isempty(data) || ndims(data) > 2 || ...
        any(~isfinite(data(:)))
    error('data must be a non-empty real numeric matrix containing finite values.')
end

if size(data,1) > size(data,2)
    error('data must be arranged as channels-by-time, with time along the second dimension.')
end

if size(data,2) < 2
    error('data must contain at least two time samples.')
end

if ~isnumeric(srate) || ~isscalar(srate) || ~isfinite(srate) || srate <= 0
    error('srate must be one positive finite scalar.')
end

if ~isnumeric(f) || ~isscalar(f) || ~isfinite(f) || f <= 0 || f >= srate/2
    error('f must be one positive finite scalar below the Nyquist frequency.')
end

if ~isnumeric(fwhm) || ~isscalar(fwhm) || ~isfinite(fwhm) || fwhm <= 0
    error('fwhm must be one positive finite scalar.')
end

if nargin<5
    showplot=false;
end

%% compute filter

% frequencies
ntime = size(data,2);
hz = linspace(0,srate,ntime);

% create Gaussian
s  = fwhm*(2*pi-1)/(4*pi); % normalized width
x  = hz-f;                 % shifted frequencies
fx = exp(-.5*(x/s).^2);    % gaussian
fx = fx./max(fx);          % gain-normalized

%% filter

filtdat = 2*real( ifft( bsxfun(@times,fft(data,[],2),fx) ,[],2) );

%% compute empirical frequency and standard deviation

[~,idx] = min(abs(hz-f));
empVals(1) = hz(idx);

% find values closest to .5 after MINUS before the peak
[~,idx_half_low] = min(abs(fx(1:idx)-.5));
[~,idx_half_high_rel] = min(abs(fx(idx:end)-.5));
idx_half_high = idx-1+idx_half_high_rel;
empVals(2) = hz(idx_half_high)-hz(idx_half_low);

% also temporal FWHM
tmp = abs(FREQNESS_AnalyticSignal(real(fftshift(ifft(fx)))));
tmp = tmp./max(tmp);
tx = (0:ntime-1)/srate;
[~,idxt] = max(tmp);
[~,idxt_half_low] = min(abs(tmp(1:idxt)-.5));
[~,idxt_half_high_rel] = min(abs(tmp(idxt:end)-.5));
idxt_half_high = idxt-1+idxt_half_high_rel;
empVals(3) = (tx(idxt_half_high)-tx(idxt_half_low))*1000;

%% inspect the Gaussian (turned off by default)

if showplot
    figure(10001+showplot),clf
    %subplot(211)
    plot(hz,fx,'k', 'LineWidth' , 1.3)
    hold on
    set(gca,'xlim',[max(f-10,0) f+10]);
    legend({'Filter kernel'})
    title(['Empirical filter: ' num2str(empVals(1)) ', ' num2str(empVals(2)) ' Hz' ])
    xlabel('Frequency (Hz)'), ylabel('Filter gain')
%     
%     subplot(212)
%     tmp1 = real(fftshift(ifft(fx))); tmp1 = tmp1./max(tmp1);
%     tmp2 = abs(FREQNESS_AnalyticSignal(tmp1));
%     plot(tx,tmp1, tx,tmp2), zoom on
%     xlabel('Time (s)'), ylabel('Amplitude gain')
end

%% done.
