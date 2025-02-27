function [GED] = FREQNESS_NetworkEstimation(data, frex, srate, varargin)
% ========================================================================
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  Please cite the first FREQNESS paper:
%  M. Rosso, G. Fernández-Rubio, P. E. Keller, E. Brattico, P. Vuust, M. L. Kringelbach, L. Bonetti.
%  FREQ-NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved Brain Networks During Auditory Stimulation.
%  Adv. Sci. 2025, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%
%  This script takes a multivariate dataset (channels × time matrix) 
%  and produces a frequency-resolved brain network landscape.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - data   : 2D matrix in channels-by-time format.
%             NOTE: While this approach can be applied to scalp sensor data,
%             source-reconstructed brain voxel data is necessary to 
%             interpret the output as brain networks.
% 
%  - frex   : Vector of frequencies (ascending order).
%
%  - srate  : Sampling rate (Hz).
%
%  The following parameters can be specified using name-value pairs:
%
%  - 'duration'        : Duration of the recording to analyze (seconds). 
%                        Default: full data duration.
%
%  - 'fwidth'          : Full-width at half-maximum (FWHM) of the frequency filter.
%                        Default: computed based on reference values.
%
%  - 'filter'          : Filter design method. Options:
%                        'logarithmic' (default) or 'linear'.
%
%  - 'regularisation'  : Shrinkage factor for covariance matrix regularization.
%                        Default: 0.01.
%
%  - 'ncomps'          : Number of components to retain per frequency.
%                        Default: 30.
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - GED       : Structure containing eigenvalues, eigenvectors, spatial 
%                activation patterns, and network time series.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso & Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 24/02/2025
% ========================================================================


% NOTE: credit for the GED implementation goes to Mike X Cohen 
% (Cohen, 2022 - 'A tutorial on generalized eigendecomposition for denoising, 
% contrast enhancement, and dimension reduction in multichannel electrophysiology')


%% Handle Optional Arguments

% Define default values
opts = struct('duration', [], 'fwidth', [], 'filter', 'logarithmic', 'regularisation', 0.01, 'ncomps', 30);

% Parse name-value pair arguments
opts = parse_name_value_pairs(opts, varargin{:});

% If duration is empty, use full data duration
if isempty(opts.duration)
    disp("Duration not provided. Using full data duration.");
    opts.duration = size(data, 2) / srate;
end

% Assign values
duration       = opts.duration;
fwidth         = opts.fwidth;
filter         = opts.filter;
regularisation = opts.regularisation;
ncomps         = opts.ncomps;

% Compute necessary parameters
pnts2keep = duration * srate;
nfrex = length(frex);
% Re-sort frequency vector in ascending order
frex = sort(frex,'ascend');

% Handle missing filter width (based on Rosso et al., 2025 - Advanced Science)
if isempty(fwidth)
    disp('Filter width not provided. Using default reference values.');
    ref_frex = 2.439; % set a reference frequency
    ref_fwhm = 0.35;  % set a reference filter width
    nfrex_above = 80; % how many frequencies above the reference
    nfrex_below = 6;  % how many frequencies below the reference
    % Compute the vector of frequencies 
    frex_above = linspace(ref_frex, ref_frex * (nfrex_above / 2), nfrex_above);
    frex_below = linspace(ref_frex / 2, ref_frex / (2 * nfrex_below), nfrex_below);
    frex_all = [frex_below(end:-1:1), frex_above];
    % Compute the vector of filter width (log-spaced)
    fwhm_above = logspace(log10(ref_fwhm), log10(ref_fwhm * (nfrex_above)), nfrex_above);
    fwhm_below = logspace(log10(ref_fwhm), log10(ref_fwhm / (nfrex_below)), nfrex_below + 1);
    fwidth_all = [fwhm_below(end:-1:2), fwhm_above];
    % Match user-requested frequencies with the reference frequencies from Rosso et al. (2025)
    [~, idx] = min(abs(frex_all - frex(1))); % matching the lowest user-defined frequency
    fwidth = fwidth_all(idx);
end

% Compute filter widths
if strcmpi(filter, 'logarithmic')
    fwidth_all = logspace(log10(fwidth), log10(fwidth * nfrex), nfrex);
elseif strcmpi(filter, 'linear')
    fwidth_all = linspace(fwidth, fwidth * nfrex, nfrex);
else
    disp("Invalid filter type. Defaulting to 'logarithmic'.");
    fwidth_all = logspace(log10(fwidth), log10(fwidth * nfrex), nfrex);
end

%% GED Computation

% Initialize output structure
GED = [];

% Assign broadband source data
broadData = data(:, 1:pnts2keep);

% Initialize GED outputs
GEDevals = zeros(ncomps, nfrex);
[GEDevecs, GEDpats] = deal(zeros(size(broadData,1), ncomps, nfrex));
GEDts = zeros(ncomps, size(broadData,2), nfrex);

% Perform frequency-resolved brain network separation via GED
for frexi = 1:nfrex % loop over input frequencies
    display(['Processing frequency #' num2str(frexi)])
    
    % Filter data
    narrowData = filterFGx(broadData,srate,frex(frexi),fwidth_all(frexi),0); % turn the last argument of the function filterFGx() to 1 to get a visualization of the Gaussian filter in frequency domain
   
    % Compute covariance matrices
    covS = cov(narrowData');  
    covR = cov(broadData');
    
    % Regularisation (to bring covariance matrices to full rank)
    covS = covS  + 1e-6*eye(size(covS)); % regularise the S covariance matrix by adding a small perturbation/noise 
    evalsR = eig(covR ); % re-compute eigenvalues of R covariance matrix for its regularisation 
    covR  = (1-regularisation)*covR  + regularisation * mean(evalsR) * eye(size(covR)); % regularise the R covariance matrix
    
    % Eigendecomposition
    [evecs,evals] = eig(covS ,covR);
    [evals,sidx]  = sort( diag(evals),'descend' ); % the first output returns sorted evals extracted from diagonal
    evecs = evecs(:,sidx);          % sort eigenvectors
    evals = evals.*100./sum(evals); % normalize eigenvalues to percent variance explained
    % Assign temporary variables to output
    GEDevecs(:,:,frexi) = evecs(:,1:ncomps);
    GEDevals(:,frexi) = evals(1:ncomps);
    
    % Brain networks' spatial activation patterns and time series
    for compi = 1:ncomps
        
        % Compute spatial activation patterns and flip sign
        GEDpats(:,compi,frexi) = GEDevecs(:,compi,frexi)' * covS; % get component
        [~,idxmax] = max(abs(GEDpats(:,compi,frexi)));     % find max magnitude
        GEDpats(:,compi,frexi)  = GEDpats(:,compi,frexi) * sign(GEDpats(idxmax,compi,frexi)); % possible sign flip
        % Take the absolute value and normalize 0-to-1 (this is necessary due to sign ambiguity from source reconstruction)
        GEDpats(:,compi,frexi) = abs( GEDpats(:,compi,frexi) / max(GEDpats(:,compi,frexi)) );
        
        % Compute the network's activation time series 
        GEDts(compi,:,frexi) = GEDevecs(:,compi,frexi)' * broadData;
        
    end
    
end


% Assign outputs to structure
GED.evals = GEDevals; % eigenvalues
GED.evecs = GEDevecs; % eigenvectors
GED.pats  = GEDpats;  % spatial activation patterns
GED.ts    = GEDts;    % time series
GED.frex  = frex;     % to make use of the analyzed frequencies later on  

end

%% Helper Function: Parse Name-Value Pairs
function opts = parse_name_value_pairs(opts, varargin)
    if mod(length(varargin), 2) ~= 0
        error('Arguments must be given as name-value pairs.');
    end
    for i = 1:2:length(varargin)
        name = lower(varargin{i});
        if isfield(opts, name)
            opts.(name) = varargin{i+1};
        else
            error(['Unrecognized argument: ', name]);
        end
    end
end

%%

