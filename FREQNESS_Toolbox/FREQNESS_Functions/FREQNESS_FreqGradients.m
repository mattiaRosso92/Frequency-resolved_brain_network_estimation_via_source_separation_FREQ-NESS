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
%  activation patterns produced by FREQNESS. For every frequency and
%  participant, absolute suprathreshold activation coefficients weight one
%  MNI center of mass. The X,Y,Z center coordinates are then modelled as a
%  function of physical frequency.
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
%                 between 8 and 12 Hz, modelling their X,Y,Z
%                 activation-weighted center coordinates.
%                 When no input is given, the function will fit to all 
%                 frequencies in the FREQ.frex field.  
%
%  - comp2model : optional argument to select one component to for all frequencies
%                 you intend to model. E.g., comp2model = 1 will operate on
%                 the 1st component for all frequencies.
%                 When no input is given, the function will default to the
%                 1st component.
%
%  - threshold_sd : number of standard deviations above the mean used to
%                 retain pattern coefficients. Default: 1.
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
%  - gradCoeff : Polynomial coefficients modelling center position from
%                frequency for each participant. Size: [3 x 3 x nSubs], where:
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
opts = struct('frex2model', [], 'comp2model', [], 'threshold_sd', 1, 'plot_all', false);
opts = parse_name_value_pairs(opts, varargin{:});

frex2model = opts.frex2model;
comp2model = opts.comp2model;
plot_all   = opts.plot_all;   % local flag to control individual-subject plotting
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
elseif ~isnumeric(comp2model) || ~isscalar(comp2model) || ~isfinite(comp2model) || ...
        comp2model ~= round(comp2model) || comp2model < 1 || comp2model > nComp
    error('comp2model must be an integer between 1 and %d.',nComp);
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
    if ~isnumeric(frex2model) || numel(frex2model) ~= 2 || any(~isfinite(frex2model))
        error(['The variable frex2model must be a 2-element vector containing the boundaries ' ...
               'of the frequency range to model as a gradient. E.g., [8 12] will model all ' ...
               'frequencies ranging from 8 to 12 Hz in FREQ.frex.']);
    end

    frex2model = sort(frex2model(:),'ascend');
    % Find closest 1st and last frequencies
    [~, idx_first] = min(abs(frex - frex2model(1)));
    [~, idx_last]  = min(abs(frex - frex2model(2)));
    idx_range2model  = min(idx_first,idx_last):max(idx_first,idx_last);

    % Warn if requested boundaries are not exact matches
    tol = 1e-3; % Hz (or "index units" if frex are indices)
    if any( ~ismembertol(frex2model(:), frex(:), tol) )
        warning('Some requested frequencies in frex2model are not present in FREQ.frex. Using closest matches instead.');
    end
end

 
% Display for the user
fprintf('\nFREQNESS Frequency Gradients: modelling spatial gradients from %.1f to %.1f Hz for component %d across %d participants.\n', ...
    frex(idx_range2model(1)),frex(idx_range2model(end)),which_comp,nsubs);

%% Process spatial activation patters

% Convert to absolute values
patterns = abs(patterns);

% Threshold within subjects
for subi = 1:nsubs
    for frexi = 1:nfrex
        % Assign
        this_pat = patterns(:,frexi,subi);
        % Normalize 0-1
        max_pat = max(this_pat);
        if max_pat > 0
            this_pat = this_pat/max_pat;
        else
            this_pat(:) = 0;
        end
        % Threshold
        thresh = mean(this_pat) + threshold_sd*std(this_pat);
        this_pat(this_pat<thresh) = nan;
        % Re-assign
        patterns(:,frexi,subi) = this_pat;
    end
end

% Compute one activation-weighted MNI centroid per frequency and participant
centers = compute_activation_centers(patterns,MNI);

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
jitter_idx = (1:size(cat_x,1))';
jitter_x = shuffler*sin(jitter_idx*sqrt(2));
jitter_y = shuffler*sin(jitter_idx*sqrt(3));
jitter_z = shuffler*sin(jitter_idx*sqrt(5));
[pats2plot,x2plot,y2plot,z2plot,size2plot] = deal(nan(size(cat_pats)));
for frexi = 1:nfrex

    % Adding frequency offsets for visual readability
    pats2plot(:,frexi) = frexi + (cat_pats(:,frexi)./cat_pats(:,frexi) - 1) ;  % with frequency offsets
    x2plot(:,frexi) = cat_x(:,frexi) + jitter_x;
    y2plot(:,frexi) = cat_y(:,frexi) + jitter_y;
    z2plot(:,frexi) = cat_z(:,frexi) + jitter_z;

    % Marker size: voxel amplitudes from cat_pats (after thresholding)
    size2plot(:,frexi) = 50*cat_pats(:,frexi);
    size2plot(size2plot(:,frexi)<=0,frexi) = nan;
    
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
sgtitle(['Spatial gradient across frequencies - Component #' num2str(which_comp) ' - Group level'],'FontSize', 18)



%% Model spatial gradients and overlay fits on scatterplots

% -------------------------------------------------------------------------
% Here we parametrize spatial gradients by fitting polynomial models to one
% activation-weighted MNI centroid per frequency and participant:
%   y = centroid coordinate (X, Y, or Z; mm)
%   x = physical frequency (Hz)
%
% We:
%   1) Select the frequency range specified by frex2model (if given),
%      using FREQ.frex (in Hz) to find the closest indices.
%   2) For each dimension (X, Y, Z), collect the participant centroids
%      across the selected frequencies.
%   3) Fit 1st- and 2nd-order polynomials:
%           centroid ~ poly(frequency)
%      and compute R^2 for both.
%   4) Use BIC only to select the best model; we then store:
%         - gradCoeff(dim,coeff,subi) for each participant
%         - goodFit.R2_* and goodFit.bestOrder per participant
%   5) Overlay a group-level polynomial (with equal participant weighting) on the
%      existing scatterplots, and add dashed y-lines marking the modelled
%      frequency range.
% -------------------------------------------------------------------------

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

for dim = 1:3  % 1 = X, 2 = Y, 3 = Z (group-level modelling)
    nSelected = numel(idx_range2model);
    predictorVec = repmat(frex(idx_range2model),nsubs,1);
    centerVec = nan(nSelected*nsubs,1);
    for subi = 1:nsubs
        this_center = centers(dim,idx_range2model,subi);
        rows = (subi-1)*nSelected+(1:nSelected);
        centerVec(rows) = this_center(:);
    end

    [gradCoeff_group(dim,:), goodFit_group.R2_linear(dim), ...
        goodFit_group.R2_quadratic(dim), goodFit_group.R2_best(dim), ...
        goodFit_group.bestOrder(dim)] = fit_center_gradient(predictorVec,centerVec);

    if isnan(goodFit_group.bestOrder(dim))
        continue
    end
    
    % ----- Overlay selected group-level model on the corresponding subplot -----
    
    % Evaluate centroid position at each selected physical frequency
    modelCoord = polyval(gradCoeff_group(dim,:),frex(idx_range2model));
    minIdx = min(idx_range2model);
    maxIdx = max(idx_range2model);
    
    % Select subplot (same as your scatter layout)
    subplot(1,3,dim); hold on
    
    % Plot fitted model as a black line
    plot(modelCoord, idx_range2model, 'k-', 'LineWidth', 2);
    
    % Add dashed lines delimiting the modelled frequency range
    yline(minIdx, 'k--', 'LineWidth', 1);
    yline(maxIdx, 'k--', 'LineWidth', 1);
    
end

% -------------------------------------------------------------------------
% Subject-wise modelling: coefficients for each participant (stored in output)
% -------------------------------------------------------------------------

for subi = 1:nsubs
    for dim = 1:3  % X, Y, Z
        this_center = centers(dim,idx_range2model,subi);
        centerVec = this_center(:);
        [gradCoeff(dim,:,subi), goodFit.R2_linear(dim,subi), ...
            goodFit.R2_quadratic(dim,subi), goodFit.R2_best(dim,subi), ...
            goodFit.bestOrder(dim,subi)] = ...
            fit_center_gradient(frex(idx_range2model),centerVec);
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
            size2plot_sub(size2plot_sub(:,frexi)<=0,frexi) = nan;
        end
        
        % ----- Subject-specific polynomial curves for plotting, from stored coeffs -----
        coordLine_sub = cell(3,1);
        modelIdx_sub  = cell(3,1);
        
        for dim = 1:3  % 1 = X, 2 = Y, 3 = Z
            
            % Skip if no coefficients stored for this subject/dimension
            if all(isnan(gradCoeff(dim,:,subi)))
                coordLine_sub{dim} = [];
                modelIdx_sub{dim}  = [];
                continue
            end
            
            % Evaluate centroid position from physical frequency
            coeff_sub     = squeeze(gradCoeff(dim,:,subi));   % [b2 b1 b0]
            coordLine_sub{dim} = polyval(coeff_sub,frex(idx_range2model));
            modelIdx_sub{dim}  = idx_range2model;
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
        sgtitle(['Spatial gradient across frequencies - Component #' num2str(which_comp) ' - Subject #' num2str(subi)],'FontSize', 18)

        
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

%% Helper Function: Compute Activation-Weighted Centers
function centers = compute_activation_centers(patterns,MNI)

centers = nan(3,size(patterns,2),size(patterns,3));
for subi = 1:size(patterns,3)
    for leveli = 1:size(patterns,2)
        weights = patterns(:,leveli,subi);
        valid = isfinite(weights) & weights>0;
        total_weight = sum(weights(valid));
        if any(valid) && isfinite(total_weight) && total_weight>0
            centers(:,leveli,subi) = ...
                (MNI(valid,:)'*weights(valid))/total_weight;
        end
    end
end

end

%% Helper Function: Fit Center Trajectory
function [coeff,r2_linear,r2_quadratic,r2_best,best_order] = ...
    fit_center_gradient(predictor,center)

coeff = nan(1,3);
[r2_linear,r2_quadratic,r2_best,best_order] = deal(nan);
predictor = predictor(:);
center = center(:);
valid = isfinite(predictor) & isfinite(center);
predictor = predictor(valid);
center = center(valid);

if numel(predictor)<3 || numel(unique(predictor))<2 || ...
        numel(unique(center))<2
    return
end

sst = sum((center-mean(center)).^2);
if ~isfinite(sst) || sst<=0
    return
end

p_linear = polyfit(predictor,center,1);
predicted_linear = polyval(p_linear,predictor);
sse_linear = sum((center-predicted_linear).^2);
r2_linear = 1-sse_linear/sst;

coeff = [0 p_linear];
r2_best = r2_linear;
best_order = 1;

if numel(predictor)<4 || numel(unique(predictor))<4
    return
end

p_quadratic = polyfit(predictor,center,2);
predicted_quadratic = polyval(p_quadratic,predictor);
sse_quadratic = sum((center-predicted_quadratic).^2);
r2_quadratic = 1-sse_quadratic/sst;

zero_tolerance = eps*max(sst,1)*numel(predictor)*32;
linear_perfect = sse_linear<=zero_tolerance;
quadratic_perfect = sse_quadratic<=zero_tolerance;
if quadratic_perfect && ~linear_perfect
    quadratic_wins = true;
elseif linear_perfect
    quadratic_wins = false;
else
    nDat = numel(predictor);
    bic_linear = nDat*log(sse_linear/nDat)+2*log(nDat);
    bic_quadratic = nDat*log(sse_quadratic/nDat)+3*log(nDat);
    quadratic_wins = bic_quadratic<bic_linear;
end

if quadratic_wins
    coeff = p_quadratic;
    r2_best = r2_quadratic;
    best_order = 2;
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
