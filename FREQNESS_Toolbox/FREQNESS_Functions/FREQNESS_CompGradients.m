function [gradCoeff, goodFit] = FREQNESS_CompGradients(FREQ, MNI, varargin)

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
%  activation patterns produced by FREQNESS across COMPONENTS, for a given
%  frequency. This is useful to test whether the activation coefficients in
%  the networks are distributed across the X,Y,Z dimensions as a function
%  of their component index, following spatial gradients. 
%
%  These operations are carried out on the FREQ.pats output produced by 
%  FREQNESS_NetworkEstimation(). The input can consist of either a 2D or
%  3D matrix, depending on whether FREQNESS was run on individual
%  participants or at the group level. 
%
%  If FREQ.pats contains participants as a 3rd/4th dimension, plotting is
%  for all participants aggregated; modelling is for individual
%  participants producing parameters for all of them; the fit is
%  visualized as grand-average in the component range of interest.
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
%  - freq2model  : optional argument to select the frequency (in Hz, if
%                  FREQ.frex is provided) at which you intend to model the
%                  gradient across components. If not specified, the
%                  function will default to the frequency showing the
%                  highest eigenvalue in the first component (i.e., the
%                  most prominent network).
%
%  - comps2model : optional argument to select a range of components you
%                  intend to model. E.g., comps2model = [1 5] will fit the
%                  best polynomial model across components 1 through 5, as
%                  a function of the X,Y,Z coordinates of their dominant
%                  activation coefficients.
%                  When no input is given, the function will fit to ALL
%                  components in the second dimension of FREQ.pats.
%
%  - threshold_sd : number of standard deviations above the mean used to
%                  retain pattern coefficients. Default: 2.
%
%  - plot_all    : logical flag (true/false). When true, the function
%                  produces additional figures for each individual
%                  participant, plotting their spatial activation patterns
%                  across components and the corresponding polynomial fits
%                  (in addition to the grand-average visualization).
%                  By default, plot_all = false and only the aggregated
%                  group-level scatterplots and fits are shown.
%
% ------------------------------------------------------------------------
%  OUTPUT ARGUMENTS:
% ------------------------------------------------------------------------
%
%  - gradCoeff : Polynomial coefficients modelling the component-wise
%                gradient for each individual participant. Size:
%                [3 x 3 x nSubs], where:
%                  • dim 1–3 correspond to X, Y, Z
%                  • the 2nd dimension holds [b0 b1 b2].
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
%  Aarhus (DK), Oxford (UK), 17/11/2025
%
% ========================================================================

%% Map inputs from FREQ 

% Handle optional arguments (name-value pairs)
opts = struct('freq2model', [], 'comps2model', [], 'threshold_sd', 2, 'plot_all', false);
opts = parse_name_value_pairs(opts, varargin{:});

freq2model  = opts.freq2model;
comps2model = opts.comps2model;
plot_all    = opts.plot_all;   % local flag to control individual-subject plotting
threshold_sd = opts.threshold_sd;

% Check main structure
if ~isstruct(FREQ)
    error('Input must be a struct with fields FREQ.pats (and optionally FREQ.frex).');
end

% Assign spatial patterns
if ~isfield(FREQ,'pats') || isempty(FREQ.pats)
    error('FREQ.pats is missing or empty.');
end
if ~isnumeric(FREQ.pats) || ~isreal(FREQ.pats) || any(~isfinite(FREQ.pats(:)))
    error('FREQ.pats must be a real numeric array containing finite values.');
end
[nvoxs, nComp, nfrex, nsubs] = size(FREQ.pats);

% Check MNI coordinates
if ~isnumeric(MNI) || ~isreal(MNI) || any(~isfinite(MNI(:))) || ...
        size(MNI,1) ~= nvoxs || size(MNI,2) ~= 3
    error('The MNI matrix must be [nVox x 3] and match the first dimension of FREQ.pats.');
end
if ~isnumeric(threshold_sd) || ~isscalar(threshold_sd) || ...
        ~isfinite(threshold_sd) || threshold_sd < 0
    error('threshold_sd must be one non-negative finite scalar.');
end

% Assign frequencies (if present)
if isfield(FREQ,'frex') && ~isempty(FREQ.frex)
    frex = FREQ.frex(:); % all frequencies in Hz
    if numel(frex) ~= nfrex
        error('Length of FREQ.frex does not match the 3rd dimension of FREQ.pats.');
    end
else
    frex = (1:nfrex)'; % fallback to indices
    if ~isempty(freq2model)
        warning(['FREQ.frex is missing or empty. The variable "freq2model" will be interpreted as a frequency INDEX, not Hz. ' ...
                 'We recommend explicitly including FREQ.frex in the input structure.']);
    else
        warning(['FREQ.frex is missing or empty. Indices will be used instead of actual frequencies expressed in Hz. ' ...
                 'We recommend explicitly including FREQ.frex in the input structure.']);
    end
end

% Pick frequency (index)
if isempty(freq2model)
    warning('No frequency was defined. Defaulting to the most prominent frequency in the network landscape.');
    if ~isfield(FREQ,'evals') || isempty(FREQ.evals)
        error(['To automatically select the most prominent network when freq2model is not specified, ' ...
               'the input structure must contain FREQ.evals (eigenvalues).']);
    end
    evals = FREQ.evals;
    if size(evals,2) ~= nfrex
        error('The 2nd dimension of FREQ.evals must match the 3rd dimension of FREQ.pats (nfrex).');
    end
    % Assume FREQ.evals is [nComp x nFrex x (nSubs)].
    % Extract eigenvalues for the first component.
    evals_first_comp = squeeze(evals(1,:,:)); % [nFrex x nSubs] or [nFrex x 1]
    % Average across participants if present.
    if isvector(evals_first_comp)
        evals_mean = abs(evals_first_comp(:));
    else
        evals_mean = mean(abs(evals_first_comp),2,'omitnan');
    end
    if numel(evals_mean) ~= nfrex
        error('Inconsistent dimensions between FREQ.evals and FREQ.pats (nfrex).');
    end
    if all(~isfinite(evals_mean))
        error('FREQ.evals does not contain finite values for automatic frequency selection.');
    end
    [~, which_frex] = max(evals_mean,[],'omitnan');
else
    if isscalar(freq2model)
        % If we have proper Hz values (FREQ.frex), treat freq2model as Hz
        if isfield(FREQ,'frex') && ~isempty(FREQ.frex)
            [~, which_frex] = min(abs(frex - freq2model));
            tol = 1e-3; % Hz
            if ~ismembertol(freq2model, frex(:), tol)
                warning('Requested freq2model is not exactly present in FREQ.frex. Using the closest available frequency instead.');
            end
        else
            % Interpret as index
            if ~isfinite(freq2model) || freq2model ~= round(freq2model) || ...
                    freq2model < 1 || freq2model > nfrex
                error(['When FREQ.frex is missing, freq2model is interpreted as an index, ' ...
                       'which must be between 1 and ', num2str(nfrex), '.']);
            end
            which_frex = freq2model;
        end
    else
        error('Input variable "freq2model" must be a scalar.');
    end
end

% Check / pick components range
if isempty(comps2model)
    disp('Component range not specified. Defaulting to all components.');
    idx_comps2model = 1:nComp;
else
    if ~isnumeric(comps2model) || numel(comps2model) ~= 2 || ...
            any(~isfinite(comps2model)) || any(comps2model ~= round(comps2model))
        error(['The variable comps2model must be a 2-element vector containing the boundaries ' ...
               'of the component range to model as a gradient. E.g., [1 5] will model all ' ...
               'components ranging from 1 to 5 in FREQ.pats.']);
    end
    comps2model = sort(comps2model(:)');
    if comps2model(1) < 1 || comps2model(2) > nComp
        error(['Values in comps2model must be between 1 and ', num2str(nComp), '.']);
    end
    idx_comps2model = comps2model(1):comps2model(2);
end

% For convenience, keep the exact Hz of the selected frequency
freq_selected = frex(which_frex);

% Display for the user
fprintf('\nFREQNESS Component Gradients: modelling spatial gradients across components %d to %d at %.1f Hz for %d participants.\n', ...
    idx_comps2model(1), idx_comps2model(end), freq_selected, nsubs);

%% Process spatial activation patterns

% Extract patterns for the selected frequency
% Resulting size: [nVox x nComp x nSubs]
patterns = squeeze( FREQ.pats(:,:,which_frex,:) );

% Ensure patterns has consistent 3D shape (even if nSubs=1)
if ndims(patterns) == 2
    patterns = reshape(patterns, [nvoxs, nComp, 1]);
    nsubs = 1;
end

% Convert to absolute values
patterns = abs(patterns);

% Threshold within subjects across components
for subi = 1:nsubs
    for compi = 1:nComp
        % Assign
        this_pat = patterns(:,compi,subi);
        % Normalize 0-1
        if max(this_pat) ~= 0
            this_pat = this_pat / max(this_pat);
        else
            this_pat(:) = 0;
        end
        % Threshold
        thresh = mean(this_pat) + threshold_sd*std(this_pat);
        this_pat(this_pat<thresh) = nan;
        % Re-assign
        patterns(:,compi,subi) = this_pat;
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

% Concatenate only coordinates along component dimension
cat_x = repmat(cat_x,[1,nComp]);
cat_y = repmat(cat_y,[1,nComp]);
cat_z = repmat(cat_z,[1,nComp]);


%% Visualize spatial distribution across XYZ spatial dimentions (components)

% Produce color map for the components
Cmap = parula(nComp);

% Labels
coordlabs = {'X-coordinates';'Y-coordinates';'Z-coordinates'};

% The following operations are for visualization only; modelling will be
% carried out on the original variables
shuffler = .15;   % jitter for spatial coordinates (horizontal)
jitter_idx = (1:size(cat_x,1))';
jitter_x = shuffler*sin(jitter_idx*sqrt(2));
jitter_y = shuffler*sin(jitter_idx*sqrt(3));
jitter_z = shuffler*sin(jitter_idx*sqrt(5));
[pats2plot,x2plot,y2plot,z2plot,size2plot] = deal(nan(size(cat_pats)));
for compi = 1:nComp

    % Adding component offsets for visual readability
    pats2plot(:,compi) = compi + (cat_pats(:,compi)./cat_pats(:,compi) - 1);  % with component offsets
    x2plot(:,compi)    = cat_x(:,compi) + jitter_x;
    y2plot(:,compi)    = cat_y(:,compi) + jitter_y;
    z2plot(:,compi)    = cat_z(:,compi) + jitter_z;

    % Marker size: voxel amplitudes from cat_pats (after thresholding)
    size2plot(:,compi) = 50*cat_pats(:,compi);
    size2plot(size2plot(:,compi)<=0,compi) = nan;

end

% Frequency label (for title)
if isfield(FREQ,'frex') && ~isempty(FREQ.frex)
    freq_label = sprintf('%.4g Hz', freq_selected);
else
    freq_label = sprintf('Index %d', which_frex);
end

% Produce plot (group-level)
figure
% X
subplot(131), hold on
scatter(x2plot, pats2plot, size2plot, Cmap, 'filled','MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2)
xlim([min(x2plot(:)) max(x2plot(:))])
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)
set(gca, 'YTick', 1:nComp, 'YTickLabel', compose('Comp %d', 1:nComp));
xlabel(coordlabs{1}, 'FontSize', 18)
ylabel('Component', 'FontSize', 18)
title(['X gradient - Group @ ', freq_label], 'FontSize', 18)
% Y
subplot(132), hold on 
scatter(y2plot, pats2plot, size2plot, Cmap, 'filled','MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2)
xlim([min(y2plot(:)) max(y2plot(:))])
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)
set(gca, 'YTick', 1:nComp, 'YTickLabel', compose('Comp %d', 1:nComp));
xlabel(coordlabs{2}, 'FontSize', 18)
ylabel('Component', 'FontSize', 18)
title(['Y gradient - Group @ ', freq_label], 'FontSize', 18)
% Z
subplot(133), hold on
scatter(z2plot, pats2plot, size2plot, Cmap, 'filled','MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2)
xlim([min(z2plot(:)) max(z2plot(:))])
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)
set(gca, 'YTick', 1:nComp, 'YTickLabel', compose('Comp %d', 1:nComp));
xlabel(coordlabs{3}, 'FontSize', 18)
ylabel('Component', 'FontSize', 18)
title(['Z gradient - Group @ ', freq_label], 'FontSize', 18)



%% Model spatial gradients across components and overlay fits on scatterplots

% -------------------------------------------------------------------------
% Here we parametrize spatial gradients by fitting polynomial models
% to the same conceptual variables used in the plots:
%   y = component index (one row per component)
%   x = spatial coordinates (X, Y, Z)
%
% We:
%   1) Select the component range specified by comps2model (if given).
%   2) For each dimension (X, Y, Z), collect all suprathreshold voxels
%      across the selected components.
%   3) Fit 1st- and 2nd-order polynomials:
%           compIdx ~ poly(coord)
%      and compute R^2 for both.
%   4) Use BIC only to select the best model; we then store subject-wise:
%         - gradCoeff(dim,coeff,subi) = [b2 b1 b0] (quadratic coefficient may be 0)
%         - goodFit.R2_* and goodFit.bestOrder per participant.
%   5) Overlay a group-level polynomial (based on concatenated data) on
%      the existing scatterplots, and add dashed y-lines marking the
%      modelled component range.
% -------------------------------------------------------------------------

% Determine component indices to model
idx_range2model = idx_comps2model;
minIdx = min(idx_range2model);
maxIdx = max(idx_range2model);

% -------------------------------------------------------------------------
% Prepare outputs (subject-wise)
% -------------------------------------------------------------------------
% gradCoeff(dim,coeff,subi) = [b2 b1 b0] for X,Y,Z; quadratic term may be 0
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
    compIdxVec = [];
    
    % Collect coord/component index pairs across selected components
    for compi = idx_range2model
        this_pat = cat_pats(:,compi);
        mask     = ~isnan(this_pat);  % suprathreshold voxels only
        
        if any(mask)
            coordVec   = [coordVec;   coordMats_group{dim}(mask,compi)];
            compIdxVec = [compIdxVec; compi*ones(sum(mask),1)];
        end
    end
    
    nDat = numel(compIdxVec);
    if nDat < 3 || numel(unique(coordVec)) < 3 || numel(unique(compIdxVec)) < 2
        % Not enough data to fit a quadratic model reliably
        continue
    end
    
    % Total variance for R^2
    sst = sum( (compIdxVec - mean(compIdxVec)).^2 );
    
    % ----- Linear model: compIdx = b1*coord + b0 -----
    p_lin = polyfit(coordVec, compIdxVec, 1);  % [b1 b0]
    y_lin = polyval(p_lin, coordVec);
    sse_lin = sum((compIdxVec - y_lin).^2);
    r2_lin  = 1 - sse_lin/sst;
    
    % ----- Quadratic model: compIdx = b2*coord^2 + b1*coord + b0 -----
    p_quad = polyfit(coordVec, compIdxVec, 2); % [b2 b1 b0]
    y_quad = polyval(p_quad, coordVec);
    sse_quad = sum((compIdxVec - y_quad).^2);
    r2_quad  = 1 - sse_quad/sst;
    
    % Store R^2 (group-level)
    goodFit_group.R2_linear(dim)    = r2_lin;
    goodFit_group.R2_quadratic(dim) = r2_quad;
    
    % ----- Model selection via BIC (group-level, for visualization only) -----
    k_lin   = 2;  % parameters: slope + intercept
    k_quad  = 3;  % parameters: quad + slope + intercept
    bic_lin  = nDat*log(max(sse_lin,eps)/nDat)  + k_lin*log(nDat);
    bic_quad = nDat*log(max(sse_quad,eps)/nDat) + k_quad*log(nDat);
    
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
    
    % Evaluate model in component-index space
    if goodFit_group.bestOrder(dim) == 2
        modelIdx = polyval(gradCoeff_group(dim,:), coordLine);          % [b2 b1 b0]
    else
        modelIdx = polyval(gradCoeff_group(dim,2:3), coordLine);        % [b1 b0]
    end
    
    % Bound model within the selected component index range
    modelIdx(modelIdx < minIdx | modelIdx > maxIdx) = nan;
    
    % Select subplot (same as your scatter layout)
    subplot(1,3,dim); hold on
    
    % Plot fitted model as a black line
    plot(coordLine, modelIdx, 'k-', 'LineWidth', 2);
    
    % Add dashed lines delimiting the modelled component range
    yline(minIdx, 'k--', 'LineWidth', 1);
    yline(maxIdx, 'k--', 'LineWidth', 1);
    
end

% -------------------------------------------------------------------------
% Subject-wise modelling: coefficients for each participant (stored in output)
% -------------------------------------------------------------------------

for subi = 1:nsubs
    
    this_pats = patterns(:,:,subi);         % [nVox x nComp]
    sub_x = repmat(MNI(:,1), [1, nComp]);
    sub_y = repmat(MNI(:,2), [1, nComp]);
    sub_z = repmat(MNI(:,3), [1, nComp]);
    coordMats_sub = {sub_x, sub_y, sub_z};
    
    for dim = 1:3  % X, Y, Z
        
        coordVec_sub   = [];
        compIdxVec_sub = [];
        
        % Collect coord/component index pairs across selected components
        for compi = idx_range2model
            this_pat_comp = this_pats(:,compi);
            mask_sub      = ~isnan(this_pat_comp);  % suprathreshold voxels only
            
            if any(mask_sub)
                coordVec_sub   = [coordVec_sub;   coordMats_sub{dim}(mask_sub,compi)];
                compIdxVec_sub = [compIdxVec_sub; compi*ones(sum(mask_sub),1)];
            end
        end
        
        nDat_sub = numel(compIdxVec_sub);
        if nDat_sub < 3 || numel(unique(coordVec_sub)) < 3 || ...
                numel(unique(compIdxVec_sub)) < 2
            % Not enough data to fit a quadratic model reliably
            continue
        end
        
        % Total variance for R^2
        sst_sub = sum( (compIdxVec_sub - mean(compIdxVec_sub)).^2 );
        
        % ----- Linear model: compIdx = b1*coord + b0 -----
        p_lin_sub = polyfit(coordVec_sub, compIdxVec_sub, 1);  % [b1 b0]
        y_lin_sub = polyval(p_lin_sub, coordVec_sub);
        sse_lin_sub = sum((compIdxVec_sub - y_lin_sub).^2);
        r2_lin_sub  = 1 - sse_lin_sub/sst_sub;
        
        % ----- Quadratic model: compIdx = b2*coord^2 + b1*coord + b0 -----
        p_quad_sub = polyfit(coordVec_sub, compIdxVec_sub, 2); % [b2 b1 b0]
        y_quad_sub = polyval(p_quad_sub, coordVec_sub);
        sse_quad_sub = sum((compIdxVec_sub - y_quad_sub).^2);
        r2_quad_sub  = 1 - sse_quad_sub/sst_sub;
        
        % Store R^2 (subject-wise)
        goodFit.R2_linear(dim,subi)    = r2_lin_sub;
        goodFit.R2_quadratic(dim,subi) = r2_quad_sub;
        
        % ----- Model selection via BIC (subject-wise) -----
        k_lin_sub   = 2;
        k_quad_sub  = 3;
        bic_lin_sub  = nDat_sub*log(max(sse_lin_sub,eps)/nDat_sub)  + k_lin_sub*log(nDat_sub);
        bic_quad_sub = nDat_sub*log(max(sse_quad_sub,eps)/nDat_sub) + k_quad_sub*log(nDat_sub);
        
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
    
    for subi = 1:nsubs
        
        % Subject-specific patterns: [nVox x nComp]
        this_pats = patterns(:,:,subi);
        
        % Subject-specific coordinates repeated across components
        sub_x = repmat(MNI(:,1), [1, nComp]);
        sub_y = repmat(MNI(:,2), [1, nComp]);
        sub_z = repmat(MNI(:,3), [1, nComp]);
        
        % Prepare variables for plotting
        [pats2plot_sub, x2plot_sub, y2plot_sub, z2plot_sub, size2plot_sub] = ...
            deal(nan(size(this_pats)));
        
        for compi = 1:nComp
            pats2plot_sub(:,compi) = compi + (this_pats(:,compi)./this_pats(:,compi) - 1);
            x2plot_sub(:,compi)    = sub_x(:,compi);
            y2plot_sub(:,compi)    = sub_y(:,compi);
            z2plot_sub(:,compi)    = sub_z(:,compi);
            size2plot_sub(:,compi) = 50 * this_pats(:,compi);
            size2plot_sub(size2plot_sub(:,compi)<=0,compi) = nan;
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
            compIdxVec_sub = [];
            
            % Collect coord/component index pairs across selected components
            for compi = idx_range2model
                this_pat_comp = this_pats(:,compi);
                mask_sub      = ~isnan(this_pat_comp);  % suprathreshold voxels only
                
                if any(mask_sub)
                    coordVec_sub   = [coordVec_sub;   coordMats_sub{dim}(mask_sub,compi)];
                    compIdxVec_sub = [compIdxVec_sub; compi*ones(sum(mask_sub),1)];
                end
            end
            
            nDat_sub = numel(compIdxVec_sub);
            if nDat_sub < 3
                coordLine_sub{dim} = [];
                modelIdx_sub{dim}  = [];
                continue
            end
            
            % Coordinate range for plotting the subject's model
            coordMin_sub = min(coordVec_sub);
            coordMax_sub = max(coordVec_sub);
            coordLine    = linspace(coordMin_sub, coordMax_sub, 200);
            
            % Evaluate model in component-index space using stored coefficients
            coeff_sub     = squeeze(gradCoeff(dim,:,subi));   % [b2 b1 b0]
            bestOrder_sub = goodFit.bestOrder(dim,subi);      % 1 or 2
            
            if bestOrder_sub == 2
                modelIdx = polyval(coeff_sub, coordLine);          % [b2 b1 b0]
            else
                modelIdx = polyval(coeff_sub(2:3), coordLine);     % [b1 b0]
            end
            
            % Bound model within the selected component index range
            modelIdx(modelIdx < minIdx | modelIdx > maxIdx) = nan;
            
            coordLine_sub{dim} = coordLine;
            modelIdx_sub{dim}  = modelIdx;
        end
        
        % Frequency label (for subject figures)
        subj_title_suffix = sprintf('@ %s', freq_label);

        
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
        set(gca, 'YTick', 1:nComp, 'YTickLabel', compose('Comp %d', 1:nComp));
        xlabel(coordlabs{1}, 'FontSize', 18)
        ylabel('Component', 'FontSize', 18)
        title(sprintf('X gradient – subject %d %s', subi, subj_title_suffix), 'FontSize', 18)
        
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
        set(gca, 'YTick', 1:nComp, 'YTickLabel', compose('Comp %d', 1:nComp));
        xlabel(coordlabs{2}, 'FontSize', 18)
        ylabel('Component', 'FontSize', 18)
        title(sprintf('Y gradient – subject %d %s', subi, subj_title_suffix), 'FontSize', 18)
        
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
        set(gca, 'YTick', 1:nComp, 'YTickLabel', compose('Comp %d', 1:nComp));
        xlabel(coordlabs{3}, 'FontSize', 18)
        ylabel('Component', 'FontSize', 18)
        title(sprintf('Z gradient – subject %d %s', subi, subj_title_suffix), 'FontSize', 18)
        
    end
end


%% Reorder polynomial coefficients in the OUTPUT (final step)

% Internal MATLAB convention used for modelling:
%   Quadratic: [b2 b1 b0]
%   Linear:    [0  b1 b0]
%
% OUTPUT convention (for gradCoeff):
%   Quadratic: [b0 b1 b2]
%   Linear:    [b0 b1 0]

for subi = 1:nsubs
    for dim = 1:3
        
        coeff = squeeze(gradCoeff(dim,:,subi));   % [b2 b1 b0] or [0 b1 b0]
        order = goodFit.bestOrder(dim,subi);      % 1 or 2
        
        % Skip if no model was estimated
        if all(isnan(coeff)) || isnan(order)
            continue
        end
        
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
