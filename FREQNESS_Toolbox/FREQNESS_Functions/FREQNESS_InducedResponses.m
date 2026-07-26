function P = FREQNESS_InducedResponses(GED,evtimes,fwhm,srate,baseline)


% This function takes as input GED as outputted from
% FREQNESS_NetworkEstimation().

% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - GED                                : Structure from FREQNESS_NetworkEstimation
%      - .ts                            : 
%      - .frex
%
%
%  - frex2embed                         : vector containing the frequencies to define the networks to embed 
%                                        
%
%
%  - Optional arguments (name-value pairs):
%      - 'frex'                         : Vector of frequencies 
%      - 'idx_comp'                     : Scalar of component index to embed in phase space, for given frequencies (default: 1)
%      - 'timeinterval'                 : [start_time end_time] in seconds for analysis (default: full range)
%      - 'threshold'                    : Fraction of max distance to define recurrences (default: 0.1)
%      - 'video'                        : 'on' or 'off' to show animated phase space plot (default: 'on')
%      - 'figure'                       : 'on' or 'off' to show the figures (default: 'on')
%


% INPUTS:
%   ts: time series for each trial (time x trial)
%   fc: central frequency of the wavelet
%   fwhm: full width half maximum of the wavelet 
%   srate : sampling rate of the data
%   baseline: array that contains the time indices for normalization
% OUTPUT:
%   P: (ntimeseries) power time serie averaged over trials
%
    P_db = zeros(size(ts));
    for tt = 1:(size(ts,2)) %over trials
        clear filtdat
        % applying a gaussian filter
        [filtdat,~] = filterFGx(ts(:,tt)',srate,fc,fwhm);
        P_raw = abs(hilbert(filtdat)).^2;
        reference = mean(P_raw(baseline));
        P_db(:,tt) = 10*log10(P_raw./reference);
        
    end
    P = mean(P_db(:,:),2);
    

end