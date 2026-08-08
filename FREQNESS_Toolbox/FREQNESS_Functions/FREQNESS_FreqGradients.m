function [gradCoeff, goodFit] = FREQNESS_FreqGradients(FREQ, MNI, varargin)

% ========================================================================
%
%  FREQUENCY-RESOLVED NETWORK ESTIMATION TOOLBOX
%
%  If you use this toolbox, please cite:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved Brain
%  Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%  This function visualizes and models spatial gradients of the spatial
%  activation patterns produced by FREQNESS. This is useful to test whether 
%  the activation coefficients in the networks are distributed across the X,Y,Z
%  dimensions as a function of their frequency, following spatial gradients. 
%
%  These operations are carried out on the FREQ.pats output produced by 
%  FREQNESS_NetworkEstimation(). The input can consist of either a 2D or
%  3D matrix, depending on whether FREQNESS was run on individual
%  participants or at the group level. 
%
%  If FREQ.pats contains participants as a 3rd dimension, [plotting for
%  all participants aggregated; modelling for individual participants
%  producing parameters for all of them; fit visualized as grand-average in
%  frequency range of interest
% ========================================================================

% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - FREQ : structure with fields
%           • FREQ.pats  -> activation patterns matrix [nVoxels x nComp x nFrex x (nSubs)]
%           • FREQ.frex  -> frequency vector [nFrex x 1] (optional but
%                           very relevant for visualization and interpretation)
%   
%  - MNI  : MNI coordinates provided in the same order as your data
%                               (N x 3, where N is the brain voxel number)
%
%  - frex2model : optional argument to select a range of frequencies of which 
%                 you intend to model the gradient. E.g., frex2model = [8
%                 12] will fit the best polynomial model to all frequencies
%                 between 8 and 12 Hz, as a function the X,Y,Z coordinates
%                 of their dominant activation coefficients.
%                 When no input is given, the function will fit to all 
%                 frequencies in the FREQ.frex field.  
%
%  - comp2model : optional argument to select one component to for all frequencies
%                 you intend to model. E.g., comp2model = 1 will operate on
%                 the 1st component for all frequencies.
%                 When no input is given, the function will default to the
%                 1st component.
%
%  - plot_all   : logical flag (true/false). When true, the function
%                 produces additional figures for each individual
%                 participant, plotting their spatial activation patterns
%                 across frequencies and the corresponding polynomial fits
%                 (in addition to the grand-average visualization).
%                 By default, plot_all = false and only the aggregated
%                 group-level scatterplots and fits are shown.
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - gradCoeff : Polynomial coefficients modelling the gradient for each
%                individual participant. Size: [3 x 3 x nSubs], where:
%                dim = 1,2,3 -> X,Y,Z; coefficients = [b0 b1 b2].
%                When the best fitting model is linear, b2 = 0.
%
%  - goodFit   : Structure containing goodness-of-fit metrics per
%                participant and dimension. Fields:
%                  • goodFit.R2_linear    [3 x nSubs]
%                  • goodFit.R2_quadratic [3 x nSubs]
%                  • goodFit.R2_best      [3 x nSubs]
%                  • goodFit.bestOrder    [3 x nSubs] (1 or 2)
%
% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso, Nikita Kudryavtsev, Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), 16/11/2025
%
% ========================================================================

%% Map inputs from FREQ 

% Handle optional arguments (name-value pairs)
opts = struct('frex2model', [], 'comp2model', [], 'plot_all', false); % default values (plot_all added)
opts = parse_name_value_pairs(opts, varargin{:});

frex2model = opts.frex2model;
comp2model = opts.comp2model;
plot_all   = opts.plot_all;   % local flag to control individual-subject plotting

% Check main structure
if ~isstruct(FREQ)
    error('Input must be a struct with fields FREQ.pats (and optionally FREQ.frex).');
end

% Assign spatial patterns
if ~isfield(FREQ,'pats') || isempty(FREQ.pats)
    error('FREQ.pats is missing or empty.');
end
[nvoxs, nComp, nfrex, nsubs] = size(FREQ.pats);

% Check MNI coordinates
if size(MNI,1) ~= nvoxs || size(MNI,2) ~= 3
    error('The MNI matrix must be [nVox x 3] and match the first dimension of FREQ.pats.');
end

% Assign frequencies
if isfield(FREQ,'frex') && ~isempty(FREQ.frex)
    frex = FREQ.frex(:); % all frequencies in Hz
    if numel(frex) ~= nfrex
        error('Length of FREQ.frex does not match the 3rd dimension of FREQ.pats.');
    end
else
    if ~isempty(frex2model)
        error(['FREQ.frex is missing or empty, but frex2model was provided in Hz. ' ...
               'Either supply FREQ.frex or omit frex2model.']);
    end
    warning(['FREQ.frex is missing or empty. Indices will be used instead of actual frequencies expressed in Hz. ' ...
             'We recommend explicitly including FREQ.frex in the input structure.']);
    frex = (1:nfrex)'; % fallback to indices
end

% Pick component
if isempty(comp2model)
    disp('Component not specified. Defaulting to analyzing the 1st component.');
    which_comp = 1;
elseif ~isscalar(comp2model) || comp2model < 1 || comp2model > nComp
    warning(['Input variable "comp2model" must be a scalar between 1 and ', num2str(nComp), ...
             '. Defaulting to analyzing the 1st component.']);
    which_comp = 1;
else
    which_comp = comp2model;
end

% Extract patterns for the selected component
patterns = squeeze( FREQ.pats(:,which_comp,:,:) ); % [nVox x nFrex x nSubs]

% Map requested freqs to nearest available in FREQ.frex
if isempty(frex2model)
    % default: use all frequencies
    idx_range2model = 1:nfrex;
else
    if numel(frex2model) ~= 2
        error(['The variable frex2model must be a 2-element vector containing the boundaries ' ...
               'of the frequency range to model as a gradient. E.g., [8 12] will model all ' ...
               'frequencies ranging from 8 to 12 Hz in FREQ.frex.']);
    end

    % Find closest 1st and last frequencies
    [~, idx_first] = min(abs(frex - frex2model(1)));
    [~, idx_last]  = min(abs(frex - frex2model(2)));
    idx_range2model  = idx_first:idx_last;

    % Warn if requested boundaries are not exact matches
    tol = 1e-3; % Hz (or "index units" if frex are indices)
    if any( ~ismembertol(frex2model(:), frex(:), tol) )
        warning('Some requested frequencies in frex2model are not present in FREQ.frex. Using closest matches instead.');
    end
end

 
% Display for the user
fprintf('\nFREQNESS Frequency Gradients: modelling spatial gradients from %.1f to %.1f Hz for component %d across %d participants.\n', ...
    frex2model(1), frex2model(end), which_comp, nsubs);

%% Process spatial activation patters

% Convert to absolute values
patterns = abs(patterns);

% Threshold within subjects
for subi = 1:nsubs
    for frexi = 1:nfrex
        % Assign
        this_pat = patterns(:,frexi,subi);
        % Normalize 0-1
        this_pat = this_pat / max(this_pat);
        % Threshold
        thresh = mean(this_pat) + std(this_pat);
        this_pat(this_pat<thresh) = nan;
        % Re-assign
        patterns(:,frexi,subi) = this_pat;
    end
end

%% Create group-level variables

% Initialize for concatenation
[cat_x,cat_y,cat_z] = deal([]); % coordinates
cat_pats = []; % patterns

% Concatenate participants
for subi = 1:nsubs

    % coordinates
    cat_x = [cat_x; MNI(:,1)];
    cat_y = [cat_y; MNI(:,2)];
    cat_z = [cat_z; MNI(:,3)];
    % patterns
    cat_pats = [cat_pats; patterns(:,:,subi)];

end

% Concatenate only coordinates along frequency dimension
cat_x = repmat(cat_x,[1,nfrex]);
cat_y = repmat(cat_y,[1,nfrex]);
cat_z = repmat(cat_z,[1,nfrex]);


%% Visualize spatial distribution across XYZ spatial dimentions

% Produce color map for the frequencies
Cmap = parula(nfrex);

% Labels
coordlabs = {'X-coordinates';'Y-coordinates';'Z-coordinates'};

% The following operations are for visualization only; modelling will be
% carried out on the original variables
shuffler = .15;   % jitter for spatial coordinates (horizontal)
[pats2plot,x2plot,y2plot,z2plot,size2plot] = deal(nan(size(cat_pats)));
for frexi = 1:nfrex

    % Adding frequency offsets for visual readability
    pats2plot(:,frexi) = frexi + (cat_pats(:,frexi)./cat_pats(:,frexi) - 1) ;  % with frequency offsets
    x2plot(:,frexi) = cat_x(:,frexi) + shuffler * randn(size(cat_x,1),1);
    y2plot(:,frexi) = cat_y(:,frexi) + shuffler * randn(size(cat_y,1),1);
    z2plot(:,frexi) = cat_z(:,frexi) + shuffler * randn(size(cat_z,1),1);

    % Marker size: voxel amplitudes from cat_pats (after thresholding)
    size2plot(:,frexi) = 50*cat_pats(:,frexi);
    
end


% Produce plot (group-level, all subjects aggregated)
figure
% X
subplot(131), hold on
scatter(x2plot, pats2plot, size2plot, Cmap, 'filled','MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2)
xlim([min(x2plot(:)) max(x2plot(:))])
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)
set(gca, 'YTick', 1:nfrex, 'YTickLabel', compose('%.4g Hz', frex));
xlabel(coordlabs{1}, 'FontSize', 18)
ylabel('Frequency', 'FontSize', 18)
title('X gradient', 'FontSize', 14)
% Y
subplot(132), hold on 
scatter(y2plot, pats2plot, size2plot, Cmap, 'filled','MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2)
xlim([min(y2plot(:)) max(y2plot(:))])
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)
set(gca, 'YTick', 1:nfrex, 'YTickLabel', compose('%.4g Hz', frex));
xlabel(coordlabs{2}, 'FontSize', 18)
ylabel('Frequency', 'FontSize', 18)
title('Y gradient', 'FontSize', 14)
% Z
subplot(133), hold on
scatter(z2plot, pats2plot, size2plot, Cmap, 'filled','MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2)
xlim([min(z2plot(:)) max(z2plot(:))])
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)
set(gca, 'YTick', 1:nfrex, 'YTickLabel', compose('%.4g Hz', frex));
xlabel(coordlabs{3}, 'FontSize', 18)
ylabel('Frequency', 'FontSize', 18)
title('Z gradient', 'FontSize', 14)
sgtitle(['Spatial gradient across frequencies - Component #' num2str(comp2model) ' - Group level'],'FontSize', 18)



%% Model spatial gradients and overlay fits on scatterplots

% -------------------------------------------------------------------------
% Here we parametrize spatial gradients by fitting polynomial models
% to the same variables used in the plots:
%   y = frequency index (one row per frex)
%   x = spatial coordinates (X, Y, Z)
%
% We:
%   1) Select the frequency range specified by frex2model (if given),
%      using FREQ.frex (in Hz) to find the closest indices.
%   2) For each dimension (X, Y, Z), collect all suprathreshold voxels
%      across the selected frequencies.
%   3) Fit 1st- and 2nd-order polynomials:
%           freqIdx ~ poly(coord)
%      and compute R^2 for both.
%   4) Use BIC only to select the best model; we then store:
%         - gradCoeff(dim,coeff,subi) for each participant
%         - goodFit.R2_* and goodFit.bestOrder per participant
%   5) Overlay a group-level polynomial (based on concatenated data) on the
%      existing scatterplots, and add dashed y-lines marking the modelled
%      frequency range.
% -------------------------------------------------------------------------

% Make sure we have a frequency vector in Hz for mapping frex2model
if exist('frex','var') ~= 1 || isempty(frex)
    % If FREQ.frex is missing, we cannot interpret frex2model in Hz
    if nargin >= 3 && ~isempty(frex2model)
        error('frex2model was provided in Hz, but FREQ.frex is missing.');
    else
        frex = (1:nfrex)'; % fall back to indices
    end
end

% Determine frequency indices to model
if nargin < 3 || isempty(frex2model)
    idx_range2model = 1:nfrex;
else
    if numel(frex2model) ~= 2
        error('frex2model must be a 1x2 vector [fmin fmax].');
    end
    [~, idx_low]  = min(abs(frex - frex2model(1)));
    [~, idx_high] = min(abs(frex - frex2model(2)));
    idx_range2model = idx_low:idx_high;
end

% -------------------------------------------------------------------------
% Prepare outputs (subject-wise)
% -------------------------------------------------------------------------
% gradCoeff(dim,coeff,subi) = [b0 b1 b2] for X,Y,Z; quadratic term may be 0
gradCoeff = nan(3,3,nsubs);

goodFit.R2_linear     = nan(3,nsubs);
goodFit.R2_quadratic  = nan(3,nsubs);
goodFit.R2_best       = nan(3,nsubs);
goodFit.bestOrder     = nan(3,nsubs);   % 1 or 2

% -------------------------------------------------------------------------
% Group-level fit for visualization only (not stored in outputs)
% -------------------------------------------------------------------------
gradCoeff_group = nan(3,3);

goodFit_group.R2_linear     = nan(3,1);
goodFit_group.R2_quadratic  = nan(3,1);
goodFit_group.R2_best       = nan(3,1);
goodFit_group.bestOrder     = nan(3,1);

% Coordinate matrices in the same format as cat_pats (group-level)
coordMats_group = {cat_x, cat_y, cat_z};

for dim = 1:3  % 1 = X, 2 = Y, 3 = Z (group-level modelling)
    
    coordVec   = [];
    freqIdxVec = [];
    
    % Collect coord/frequency index pairs across selected frequencies
    for frexi = idx_range2model
        this_pat = cat_pats(:,frexi);
        mask     = ~isnan(this_pat);  % suprathreshold voxels only
        
        if any(mask)
            coordVec   = [coordVec;   coordMats_group{dim}(mask,frexi)];
            freqIdxVec = [freqIdxVec; frexi*ones(sum(mask),1)];
        end
    end
    
    nDat = numel(freqIdxVec);
    if nDat < 3
        % Not enough data to fit a quadratic model reliably
        continue
    end
    
    % Total variance for R^2
    sst = sum( (freqIdxVec - mean(freqIdxVec)).^2 );
    
    % ----- Linear model: freqIdx = b1*coord + b0 -----
    p_lin = polyfit(coordVec, freqIdxVec, 1);  % [b1 b0]
    y_lin = polyval(p_lin, coordVec);
    sse_lin = sum((freqIdxVec - y_lin).^2);
    r2_lin  = 1 - sse_lin/sst;
    
    % ----- Quadratic model: freqIdx = b2*coord^2 + b1*coord + b0 -----
    p_quad = polyfit(coordVec, freqIdxVec, 2); % [b2 b1 b0]
    y_quad = polyval(p_quad, coordVec);
    sse_quad = sum((freqIdxVec - y_quad).^2);
    r2_quad  = 1 - sse_quad/sst;
    
    % Store R^2 (group-level)
    goodFit_group.R2_linear(dim)    = r2_lin;
    goodFit_group.R2_quadratic(dim) = r2_quad;
    
    % ----- Model selection via BIC (group-level, for visualization only) -----
    k_lin   = 2;  % parameters: slope + intercept
    k_quad  = 3;  % parameters: quad + slope + intercept
    bic_lin  = nDat*log(sse_lin/nDat)  + k_lin*log(nDat);
    bic_quad = nDat*log(sse_quad/nDat) + k_quad*log(nDat);
    
    if bic_quad < bic_lin
        % Quadratic wins
        gradCoeff_group(dim,:)      = p_quad;    % [b2 b1 b0]
        goodFit_group.R2_best(dim)  = r2_quad;
        goodFit_group.bestOrder(dim)= 2;
    else
        % Linear wins (pad quadratic term with 0)
        gradCoeff_group(dim,:)      = [0 p_lin]; % [0 b1 b0]
        goodFit_group.R2_best(dim)  = r2_lin;
        goodFit_group.bestOrder(dim)= 1;
    end
    
    % ----- Overlay selected group-level model on the corresponding subplot -----
    
    % Coordinate range for plotting the model
    coordMin = min(coordVec);
    coordMax = max(coordVec);
    coordLine = linspace(coordMin, coordMax, 200);
    
    % Evaluate model in frequency-index space
    if goodFit_group.bestOrder(dim) == 2
        modelIdx = polyval(gradCoeff_group(dim,:), coordLine);          % [b2 b1 b0]
    else
        modelIdx = polyval(gradCoeff_group(dim,2:3), coordLine);        % [b1 b0]
    end
    
    % Bound model within the selected frequency index range
    minIdx = min(idx_range2model);
    maxIdx = max(idx_range2model);
    modelIdx(modelIdx < minIdx | modelIdx > maxIdx) = nan;
    
    % Select subplot (same as your scatter layout)
    subplot(1,3,dim); hold on
    
    % Plot fitted model as a black line
    plot(coordLine, modelIdx, 'k-', 'LineWidth', 2);
    
    % Add dashed lines delimiting the modelled frequency range
    yline(minIdx, 'k--', 'LineWidth', 1);
    yline(maxIdx, 'k--', 'LineWidth', 1);
    
end

% -------------------------------------------------------------------------
% Subject-wise modelling: coefficients for each participant (stored in output)
% -------------------------------------------------------------------------

for subi = 1:nsubs
    
    this_pats = patterns(:,:,subi);         % [nVox x nFrex]
    sub_x = repmat(MNI(:,1), [1, nfrex]);
    sub_y = repmat(MNI(:,2), [1, nfrex]);
    sub_z = repmat(MNI(:,3), [1, nfrex]);
    coordMats_sub = {sub_x, sub_y, sub_z};
    
    for dim = 1:3  % X, Y, Z
        
        coordVec_sub   = [];
        freqIdxVec_sub = [];
        
        % Collect coord/frequency index pairs across selected frequencies
        for frexi = idx_range2model
            this_pat_freq = this_pats(:,frexi);
            mask_sub      = ~isnan(this_pat_freq);  % suprathreshold voxels only
            
            if any(mask_sub)
                coordVec_sub   = [coordVec_sub;   coordMats_sub{dim}(mask_sub,frexi)];
                freqIdxVec_sub = [freqIdxVec_sub; frexi*ones(sum(mask_sub),1)];
            end
        end
        
        nDat_sub = numel(freqIdxVec_sub);
        if nDat_sub < 3
            % Not enough data to fit a quadratic model reliably
            continue
        end
        
        % Total variance for R^2
        sst_sub = sum( (freqIdxVec_sub - mean(freqIdxVec_sub)).^2 );
        
        % ----- Linear model: freqIdx = b1*coord + b0 -----
        p_lin_sub = polyfit(coordVec_sub, freqIdxVec_sub, 1);  % [b1 b0]
        y_lin_sub = polyval(p_lin_sub, coordVec_sub);
        sse_lin_sub = sum((freqIdxVec_sub - y_lin_sub).^2);
        r2_lin_sub  = 1 - sse_lin_sub/sst_sub;
        
        % ----- Quadratic model: freqIdx = b2*coord^2 + b1*coord + b0 -----
        p_quad_sub = polyfit(coordVec_sub, freqIdxVec_sub, 2); % [b2 b1 b0]
        y_quad_sub = polyval(p_quad_sub, coordVec_sub);
        sse_quad_sub = sum((freqIdxVec_sub - y_quad_sub).^2);
        r2_quad_sub  = 1 - sse_quad_sub/sst_sub;
        
        % Store R^2 (subject-wise)
        goodFit.R2_linear(dim,subi)    = r2_lin_sub;
        goodFit.R2_quadratic(dim,subi) = r2_quad_sub;
        
        % ----- Model selection via BIC (subject-wise) -----
        k_lin_sub   = 2;
        k_quad_sub  = 3;
        bic_lin_sub  = nDat_sub*log(sse_lin_sub/nDat_sub)  + k_lin_sub*log(nDat_sub);
        bic_quad_sub = nDat_sub*log(sse_quad_sub/nDat_sub) + k_quad_sub*log(nDat_sub);
        
        if bic_quad_sub < bic_lin_sub
            % Quadratic wins
            gradCoeff(dim,:,subi)       = p_quad_sub;    % [b2 b1 b0]
            goodFit.R2_best(dim,subi)   = r2_quad_sub;
            goodFit.bestOrder(dim,subi) = 2;
        else
            % Linear wins (pad quadratic term with 0)
            gradCoeff(dim,:,subi)       = [0 p_lin_sub]; % [0 b1 b0]
            goodFit.R2_best(dim,subi)   = r2_lin_sub;
            goodFit.bestOrder(dim,subi) = 1;
        end
        
    end
end


%% Optional: individual-subject plots controlled by 'plot_all' flag
% When plot_all = true, produce subject-wise spatial gradient plots and
% overlay subject-specific polynomial fits, using the coefficients stored
% in gradCoeff and the best model order in goodFit.bestOrder.

if plot_all
    
    % Frequency index bounds already computed above:
    minIdx = min(idx_range2model);
    maxIdx = max(idx_range2model);
    
    for subi = 1:nsubs
        
        % Subject-specific patterns: [nVox x nFrex]
        this_pats = patterns(:,:,subi);
        
        % Subject-specific coordinates repeated across frequencies
        sub_x = repmat(MNI(:,1), [1, nfrex]);
        sub_y = repmat(MNI(:,2), [1, nfrex]);
        sub_z = repmat(MNI(:,3), [1, nfrex]);
        
        % Prepare variables for plotting
        [pats2plot_sub, x2plot_sub, y2plot_sub, z2plot_sub, size2plot_sub] = ...
            deal(nan(size(this_pats)));
        
        for frexi = 1:nfrex
            pats2plot_sub(:,frexi) = frexi + (this_pats(:,frexi)./this_pats(:,frexi) - 1);
            x2plot_sub(:,frexi)    = sub_x(:,frexi);
            y2plot_sub(:,frexi)    = sub_y(:,frexi);
            z2plot_sub(:,frexi)    = sub_z(:,frexi);
            size2plot_sub(:,frexi) = 50 * this_pats(:,frexi);
        end
        
        % ----- Subject-specific polynomial curves for plotting, from stored coeffs -----
        coordMats_sub = {sub_x, sub_y, sub_z};
        coordLine_sub = cell(3,1);
        modelIdx_sub  = cell(3,1);
        
        for dim = 1:3  % 1 = X, 2 = Y, 3 = Z
            
            % Skip if no coefficients stored for this subject/dimension
            if all(isnan(gradCoeff(dim,:,subi)))
                coordLine_sub{dim} = [];
                modelIdx_sub{dim}  = [];
                continue
            end
            
            coordVec_sub   = [];
            freqIdxVec_sub = [];
            
            % Collect coord/frequency index pairs across selected frequencies
            for frexi = idx_range2model
                this_pat_freq = this_pats(:,frexi);
                mask_sub      = ~isnan(this_pat_freq);  % suprathreshold voxels only
                
                if any(mask_sub)
                    coordVec_sub   = [coordVec_sub;   coordMats_sub{dim}(mask_sub,frexi)];
                    freqIdxVec_sub = [freqIdxVec_sub; frexi*ones(sum(mask_sub),1)];
                end
            end
            
            nDat_sub = numel(freqIdxVec_sub);
            if nDat_sub < 3
                coordLine_sub{dim} = [];
                modelIdx_sub{dim}  = [];
                continue
            end
            
            % Coordinate range for plotting the subject's model
            coordMin_sub = min(coordVec_sub);
            coordMax_sub = max(coordVec_sub);
            coordLine    = linspace(coordMin_sub, coordMax_sub, 200);
            
            % Evaluate model in frequency-index space using stored coefficients
            coeff_sub     = squeeze(gradCoeff(dim,:,subi));   % [b2 b1 b0]
            bestOrder_sub = goodFit.bestOrder(dim,subi);      % 1 or 2
            
            if bestOrder_sub == 2
                modelIdx = polyval(coeff_sub, coordLine);          % [b2 b1 b0]
            else
                modelIdx = polyval(coeff_sub(2:3), coordLine);     % [b1 b0]
            end
            
            % Bound model within the selected frequency index range
            modelIdx(modelIdx < minIdx | modelIdx > maxIdx) = nan;
            
            coordLine_sub{dim} = coordLine;
            modelIdx_sub{dim}  = modelIdx;
        end
        
        % ----- Produce subject-specific plot with overlaid fits -----
        figure
        % X
        subplot(131), hold on
        scatter(x2plot_sub, pats2plot_sub, size2plot_sub, Cmap, 'filled', ...
                'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2)
        if ~isempty(coordLine_sub{1})
            plot(coordLine_sub{1}, modelIdx_sub{1}, 'k-', 'LineWidth', 2);
        end
        yline(minIdx, 'k--', 'LineWidth', 1);
        yline(maxIdx, 'k--', 'LineWidth', 1);
        xlim([min(x2plot_sub(:)) max(x2plot_sub(:))])
        set(gcf, 'Color', 'w')
        set(gca, 'FontSize', 14)
        set(gca, 'YTick', 1:nfrex, 'YTickLabel', compose('%.4g Hz', frex));
        xlabel(coordlabs{1}, 'FontSize', 18)
        ylabel('Frequency', 'FontSize', 18)
        title('X gradient', 'FontSize', 14)
        
        % Y
        subplot(132), hold on
        scatter(y2plot_sub, pats2plot_sub, size2plot_sub, Cmap, 'filled', ...
                'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2)
        if ~isempty(coordLine_sub{2})
            plot(coordLine_sub{2}, modelIdx_sub{2}, 'k-', 'LineWidth', 2);
        end
        yline(minIdx, 'k--', 'LineWidth', 1);
        yline(maxIdx, 'k--', 'LineWidth', 1);
        xlim([min(y2plot_sub(:)) max(y2plot_sub(:))])
        set(gcf, 'Color', 'w')
        set(gca, 'FontSize', 14)
        set(gca, 'YTick', 1:nfrex, 'YTickLabel', compose('%.4g Hz', frex));
        xlabel(coordlabs{2}, 'FontSize', 18)
        ylabel('Frequency', 'FontSize', 18)
        title('Y gradient', 'FontSize', 14)
        
        % Z
        subplot(133), hold on
        scatter(z2plot_sub, pats2plot_sub, size2plot_sub, Cmap, 'filled', ...
                'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2)
        if ~isempty(coordLine_sub{3})
            plot(coordLine_sub{3}, modelIdx_sub{3}, 'k-', 'LineWidth', 2);
        end
        yline(minIdx, 'k--', 'LineWidth', 1);
        yline(maxIdx, 'k--', 'LineWidth', 1);
        xlim([min(z2plot_sub(:)) max(z2plot_sub(:))])
        set(gcf, 'Color', 'w')
        set(gca, 'FontSize', 14)
        set(gca, 'YTick', 1:nfrex, 'YTickLabel', compose('%.4g Hz', frex));
        xlabel(coordlabs{3}, 'FontSize', 18)
        ylabel('Frequency', 'FontSize', 18)
        title('Z gradient' , 'FontSize', 14)
        sgtitle(['Spatial gradient across frequencies - Component #' num2str(comp2model) ' - Subject #' num2str(subi)],'FontSize', 18)

        
    end
end

%% Reorder polynomial coefficients in the OUTPUT (final step)

% Internal MATLAB convention used for modelling:
%   Quadratic: [b2 b1 b0]
%   Linear:    [0  b1 b0]
%
% OUTPUT convention:
%   Quadratic: [b0 b1 b2]
%   Linear:    [b0 b1 0]

for subi = 1:nsubs
    for dim = 1:3
        
        coeff = squeeze(gradCoeff(dim,:,subi));  % [b2 b1 b0] or [0 b1 b0]
        order = goodFit.bestOrder(dim,subi);     % 1 or 2
        
        if order == 2
            % Quadratic: [b2 b1 b0] -> [b0 b1 b2]
            gradCoeff(dim,:,subi) = [coeff(3) coeff(2) coeff(1)];
            
        elseif order == 1
            % Linear: [0 b1 b0] -> [b0 b1 0]
            gradCoeff(dim,:,subi) = [coeff(3) coeff(2) 0];
        end
        
    end
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
