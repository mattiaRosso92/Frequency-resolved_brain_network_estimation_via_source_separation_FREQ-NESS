function [dataClean, removedActivity] = FREQNESS_NetworkRemoval(FREQ, data, freq2remove, varargin)

% ========================================================================
%
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you use this toolbox, please cite:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain
%  Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%  This function removes the broadband voxel-space contribution of selected
%  FREQNESS networks from the original data.
%
%  Network activity is reconstructed through FREQNESS_BackProjection and
%  subtracted from the broadband input:
%
%      dataClean = data - removedActivity
%
%  The supplied data must be the same data, with the same participant and
%  voxel ordering, used to compute the input FREQ structure.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - FREQ : output structure produced by FREQNESS_NetworkEstimation.
%
%  - data : original broadband data in [nVoxels x nTime] format for one
%           participant or [nVoxels x nTime x nSubs] format for multiple
%           participants.
%
%  - freq2remove : scalar frequency in Hz identifying the GED solution
%                  whose network contribution will be removed. If the
%                  requested value is not present in FREQ.frex, the closest
%                  available frequency is used.
%
% ------------------------------------------------------------------------
%  OPTIONAL NAME–VALUE PAIRS:
% ------------------------------------------------------------------------
%
%  - 'comps2remove' : component index or vector of component indices to
%                     remove. Default: 1.
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - dataClean       : broadband data after subtraction of the selected
%                      network contribution. Dimensions match data.
%
%  - removedActivity : backprojected voxel-space contribution that was
%                      subtracted from data. Dimensions match data.
%
%  The outputs satisfy, within numerical precision:
%
%      data = dataClean + removedActivity
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

if nargin < 3 || isempty(freq2remove)
    error(['freq2remove is required. Specify one frequency in Hz from ' ...
        'the range covered by FREQ.frex.']);
end

if ~isnumeric(data) || isempty(data)
    error('data must be a non-empty numeric matrix.');
end

if ndims(data) > 3
    error('data must be [nVoxels x nTime] or [nVoxels x nTime x nSubs].');
end

%% Handle optional arguments (name-value pairs)

opts = struct('comps2remove', []);
opts = parse_name_value_pairs(opts, varargin{:});

comps2remove = opts.comps2remove;

if isempty(comps2remove)
    comps2remove = 1;
    fprintf('\nFREQNESS NetworkRemoval: no components specified. Defaulting to component #1.\n');
end

%% Validate correspondence between FREQ and data

if ~isstruct(FREQ)
    error('Input must be a FREQ structure.');
end

if ~isfield(FREQ,'evecs') || isempty(FREQ.evecs)
    error('FREQ.evecs is missing or empty.');
end

if ~isfield(FREQ,'ts') || isempty(FREQ.ts)
    error('FREQ.ts is missing or empty.');
end

nvoxs = size(FREQ.evecs,1);
ntime = size(FREQ.ts,2);
nsubs = size(FREQ.evecs,4);

if size(data,1) ~= nvoxs
    error('The number of voxels in data does not match FREQ.evecs.');
end

if size(data,2) ~= ntime
    error(['The number of timepoints in data does not match FREQ.ts. ' ...
        'Provide the same analyzed data segment used by FREQNESS_NetworkEstimation.']);
end

if size(data,3) ~= nsubs
    error('The number of participants in data does not match the FREQ structure.');
end

%% Backproject and remove selected network activity

removedActivity = FREQNESS_BackProjection(FREQ,freq2remove, ...
    'comps2project',comps2remove);

dataClean = data - removedActivity;

% Verify the subtraction identity within floating-point precision
reconstructionError = max(abs(data(:) - dataClean(:) - removedActivity(:)));
tolerance = 10 * eps(max(1,max(abs(data(:)))));

if reconstructionError > tolerance
    warning(['The subtraction identity exceeded the expected numerical ' ...
        'tolerance (maximum error: %.3g).'],reconstructionError);
end

fprintf(['FREQNESS NetworkRemoval: removed component(s) %s. ' ...
    'Maximum subtraction error: %.3g.\n'], ...
    mat2str(comps2remove),reconstructionError);

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
