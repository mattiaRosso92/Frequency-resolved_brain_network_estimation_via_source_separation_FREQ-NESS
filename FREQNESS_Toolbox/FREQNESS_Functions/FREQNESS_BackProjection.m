function backProj = FREQNESS_BackProjection(FREQ, freq2project, varargin)

% ========================================================================
%
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you find this function useful, please cite the first FREQNESS paper:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain
%  Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%  This function backprojects selected FREQNESS components from component
%  space into voxel space. The output represents the broadband voxel-space
%  contribution of networks estimated at one selected frequency.
%
%  FREQNESS_NetworkEstimation computes the component time series as
%
%      Y = W' * X
%
%  where W contains the GED eigenvectors and X is the broadband voxel
%  data. Here, a reconstruction-consistent forward model is obtained from
%  the complete set of retained eigenvectors:
%
%      A = pinv(W')
%
%  Selected network contributions are then reconstructed as:
%
%      X_network = A(:,components) * Y(components,:)
%
%  The forward model is computed before selecting components. This
%  preserves the dual relationship between each retained spatial filter
%  and its corresponding reconstruction vector.
%
%  NOTE: FREQ.pats contains target-frequency activation patterns intended
%  for spatial interpretation. These patterns are not used here because
%  their scaling is not the inverse mapping paired with FREQ.ts.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - FREQ : structure with fields
%           • FREQ.evecs -> GED eigenvectors
%                           [nVoxels x nComp x nFrex x (nSubs)]
%           • FREQ.ts    -> broadband component time series
%                           [nComp x nTime x nFrex x (nSubs)]
%           • FREQ.frex  -> frequency vector [nFrex x 1]
%
%  - freq2project : scalar frequency in Hz identifying the GED solution
%                   to backproject. If the requested value is not present
%                   in FREQ.frex, the closest available frequency is used.
%
% ------------------------------------------------------------------------
%  OPTIONAL NAME–VALUE PAIRS:
% ------------------------------------------------------------------------
%
%  - 'comps2project' : component index or vector of component indices to
%                      backproject. Default: 1. When multiple components
%                      are selected, their voxel-space contributions are
%                      summed.
%
% ------------------------------------------------------------------------
%  OUTPUT:
% ------------------------------------------------------------------------
%
%  - backProj : reconstructed broadband network activity in voxel space.
%               Its dimensions are [nVoxels x nTime] for one participant
%               and [nVoxels x nTime x nSubs] for multiple participants.
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso & Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK)
%
% ========================================================================

%% Check mandatory inputs

if nargin < 2 || isempty(freq2project)
    error(['freq2project is required. Specify one frequency in Hz from ' ...
        'the range covered by FREQ.frex.']);
end

if ~isstruct(FREQ)
    error('Input must be a FREQ structure.');
end

required_fields = {'evecs','ts','frex'};
for fieldi = 1:numel(required_fields)
    if ~isfield(FREQ,required_fields{fieldi}) || isempty(FREQ.(required_fields{fieldi}))
        error('FREQ.%s is missing or empty.',required_fields{fieldi});
    end
end

if ~isnumeric(freq2project) || ~isscalar(freq2project) || ~isfinite(freq2project)
    error('freq2project must be one finite scalar frequency expressed in Hz.');
end

%% Handle optional arguments (name-value pairs)

opts = struct('comps2project', []);
opts = parse_name_value_pairs(opts, varargin{:});

comps2project = opts.comps2project;

if isempty(comps2project)
    comps2project = 1;
    fprintf('\nFREQNESS BackProjection: no components specified. Defaulting to component #1.\n');
end

%% Map and validate inputs from FREQ

frex  = FREQ.frex(:);
evecs = FREQ.evecs;
ts    = FREQ.ts;

nvoxs  = size(evecs,1);
ncomps = size(evecs,2);
nfrex  = size(evecs,3);
nsubs  = size(evecs,4);

if numel(frex) ~= nfrex
    error('Length of FREQ.frex does not match the 3rd dimension of FREQ.evecs.');
end

if size(ts,1) ~= ncomps || size(ts,3) ~= nfrex || size(ts,4) ~= nsubs
    error('FREQ.evecs and FREQ.ts have incompatible dimensions.');
end

if ~isnumeric(comps2project) || ~isvector(comps2project) || ...
        any(~isfinite(comps2project)) || any(comps2project ~= round(comps2project))
    error('comps2project must contain finite integer component indices.');
end

comps2project = unique(comps2project(:)','stable');

if any(comps2project < 1) || any(comps2project > ncomps)
    error('comps2project must contain indices between 1 and %d.',ncomps);
end

% Map requested frequency to the closest frequency analyzed by FREQNESS
[~, which_frex] = min(abs(frex - freq2project));
freq2project_actual = frex(which_frex);

tol = 1e-10;
if ~ismembertol(freq2project,frex,tol)
    warning(['Requested freq2project (%.3f Hz) is not present in FREQ.frex. ' ...
        'Using the closest available frequency: %.3f Hz.'], ...
        freq2project,freq2project_actual);
end

fprintf(['\nFREQNESS BackProjection: backprojecting component(s) %s at ' ...
    '%.3f Hz for %d participants.\n'], ...
    mat2str(comps2project),freq2project_actual,nsubs);

%% Backproject selected components into voxel space

ntime = size(ts,2);
backProj = zeros(nvoxs,ntime,nsubs);

for subi = 1:nsubs

    % Complete retained GED filter set for this frequency and participant
    W = reshape(evecs(:,:,which_frex,subi),nvoxs,ncomps);

    % Reconstruction-consistent forward model paired with Y = W' * X
    A = pinv(W');

    % Stored broadband component time series
    Y = reshape(ts(:,:,which_frex,subi),ncomps,ntime);

    % Summed voxel-space contribution of the selected networks
    backProj(:,:,subi) = A(:,comps2project) * Y(comps2project,:);

end

% Preserve a 2D output for single-participant inputs
if nsubs == 1
    backProj = backProj(:,:,1);
end

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
