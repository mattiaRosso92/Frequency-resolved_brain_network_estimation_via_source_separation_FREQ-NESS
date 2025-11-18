function FREQNESS_Visualizer(FREQ,Landscape,Patterns,varargin)

% ========================================================================
%  FREQNESS VISUALIZER: NETWORK LANDSCAPE & BRAIN ACTIVATION PATTERNS
%
%  Please cite the first FREQNESS paper:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
%  Kringelbach, M. L., & Bonetti, L. (2025).
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain
%  Networks During Auditory Stimulation.
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%
%  This script visualizes:
%
%  - Network prominence across frequencies
%    Quantified as % variance explained based on the associated
%    eigenvalues, producing a frequency-resolved network landscape
%    (see Rosso et al., 2025, Advanced Science, Figure 2).
%
%  - Spatial activation patterns of the networks
%    - All requested networks are plotted together in a 3D visualization of the brain.
%    - NIFTI files are generated for each network, allowing further
%      inspection in FSLeyes or similar software for visualization.
%
%  If multiple participants are included in the input FREQ structure, the
%  function will by default plot only grand-average landscapes and
%  activation patterns. Individual participants can be plotted by setting
%  the optional flag 'plot_all' to true:
%
%      FREQNESS_Visualizer(FREQ,Landscape,Patterns,'plot_all',true)
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - FREQ              : Structure outputted by FREQNESS_NetworkEstimation function.
%                        We recommend not to modify this structure before giving it
%                        as input to the current function (FREQNESS_Visualizer).
%      - FREQ.evals    : Normalized eigenvalues (variance explained in % points).
%      - FREQ.pats     : Spatial activation patterns used to generate NIFTI files
%                        for visualizing network topographies.
%      - FREQ.frex     : Vector of frequencies analyzed; used to select subsets
%                        for visualization.
%
%
%  - Landscape                : Structure containing settings for network landscape visualization
%      - Landscape.frex       : Vector of frequencies to visualize in the network landscape.
%      - Landscape.ncomps     : Number of components (networks) to visualize in the landscape.
%
%
%  - Patterns                 : Structure containing settings for activation patterns visualization
%      - Patterns.MNI_coords  : MNI coordinates provided in the same order as your data
%                               (N x 3, where N is the brain voxel number)
%      - Patterns.frex        : Vector of frequencies to visualize the associated
%                               networks' spatial activation patterns.
%      - Patterns.ncomps      : Number of components (networks) to visualize as
%                               spatial activation patterns.
%      - Patterns.path_output : Output path to save NIFTI images.
%
%  - Optional name-value pair:
%
%      - 'plot_all'        : Logical flag (default = false).
%                            If multiple participants are present in FREQ,
%                            the default behaviour is to plot only grand-average
%                            landscapes and spatial activation patterns.
%                            Setting 'plot_all' to true additionally plots
%                            figures for all individual participants:
%
%                            FREQNESS_Visualizer(FREQ,Landscape,Patterns,'plot_all',true)
%
%                            This flag only affects figure generation and does
%                            not change the numerical results or the NIFTI output.
%
%
%
%  NOTE 1: While the 3D plot produced by this function (Output #2)
%  is a convenient way to quickly inspect the topographies of
%  multiple networks at once, it is recommended to use the
%  NIFTI files (Output #3) for an accurate depiction of the
%  individual networks' topographies.
%  NIFTI files can be visualized using FSLeyes or an equivalent
%  software.
%
%  NOTE 2: Output #2 (3D plot) is supported for any brain MNI space (e.g. 1,2,8mm, etc).
%          Output #3 (nifti image) is currently supported only for 8mm.
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


% NOTE: we acknowledge the NIFTI Toolbox, which is used by FREQNESS to generate nifti images.
% Jimmy Shen (2025). Tools for NIfTI and ANALYZE image
% (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
% MATLAB Central File Exchange. Retrieved February 26, 2025.


%% Optional arguments

plot_all = false; % default: only grand-average plots when multiple subjects

if ~isempty(varargin)
    if mod(numel(varargin),2) ~= 0
        error('Optional arguments must be provided as name-value pairs.');
    end
    for iArg = 1:2:numel(varargin)
        name  = varargin{iArg};
        value = varargin{iArg+1};
        if ischar(name) || isstring(name)
            switch lower(char(name))
                case 'plot_all'
                    plot_all = logical(value);
                otherwise
                    warning('Unknown option "%s" ignored.', char(name));
            end
        end
    end
end


%% Input structure check

% Compute number of dimensions
ndimens = length(size(FREQ.evals));

% Check input size to verify if more than one participant is being analyzed
if ndimens < 3 % single subject in the input
    nsubs = 1;
else
    nsubs = size(FREQ.evals,3); % group of subjects in the input
end
disp(['Input FREQ structure given for ' num2str(nsubs) ' participants.']);



%% Output #1: NETWORK LANDSCAPE

% Define and sort frequencies to plot
nfrex = length(Landscape.frex);
Landscape.frex = sort(Landscape.frex,'ascend');

% Map requested freqs to nearest available in FREQ.frex
idx_frex2plot = zeros(1,nfrex);
for k = 1:nfrex
    [~, idx_temp] = min(abs(FREQ.frex - Landscape.frex(k)));
    idx_frex2plot(k) = idx_temp;
end

% Warn if some requested freqs are not exact matches
tol = 1e-3; % Hz
if ~all(ismembertol(Landscape.frex(:), FREQ.frex(:), tol))
    warning('Some requested frequencies in Landscape.frex are not present in FREQ.frex. Using closest matches instead.');
end

% Colors (red->blue gradient)
numLines = Landscape.ncomps;
colors = flipud([linspace(0,1,numLines)', zeros(numLines,1), linspace(1,0,numLines)']);

% Check requested components
if Landscape.ncomps > size(FREQ.evals,1)
    error('The requested number of components exceeds the number of estimated networks!');
end
comps2plot = 1:Landscape.ncomps;


% -------- Individual participants --------
% By default, plot individual participants only if there is a single subject.
if plot_all || nsubs == 1
    for subi = 1:nsubs
        figure; hold on

        for compi = comps2plot
            y = squeeze(FREQ.evals(compi, idx_frex2plot, subi));   % 1 x nfrex
            if  compi == 1
                plot(FREQ.frex(idx_frex2plot), y, 'Color', colors(compi,:), 'LineWidth', 3);
            else
                plot(FREQ.frex(idx_frex2plot), y, 'Color', colors(compi,:), 'LineWidth', 2);
            end
        end
        xlim([Landscape.frex(1) Landscape.frex(end)]);
        ylim([0 max(FREQ.evals(:))+.1*max(FREQ.evals(:))])
        xticks(FREQ.frex(idx_frex2plot));
        xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Explained variance (%)', 'FontSize', 14, 'FontWeight', 'bold');
        title(['Brain network landscape - Participant #' num2str(subi)], 'FontSize', 16, 'FontWeight', 'bold');


        legend(arrayfun(@(x) sprintf('GED component #%d', x), comps2plot, 'UniformOutput', false), ...
            'Location', 'Northeast', 'FontSize', 12);
        grid minor; set(gca, 'FontSize', 12, 'LineWidth', 1.5); box on; set(gcf, 'Color', 'w');
    end
end


% -------- Grand-average (only if group) --------
if nsubs > 1

    % mean/SEM over subjects (3rd dim)
    avg_evals = squeeze( mean( FREQ.evals(comps2plot, idx_frex2plot, :), 3 ) );   % [ncomps x nfrex]
    sem_evals = squeeze( std ( FREQ.evals(comps2plot, idx_frex2plot, :), 0, 3 ) / sqrt(nsubs) );

    figure; hold on
    for compi = 1:numel(comps2plot)
        yi  = avg_evals(compi,:);    % 1 x nfrex
        sei = sem_evals(compi,:);    % 1 x nfrex
        if compi == 1
            errorbar(FREQ.frex(idx_frex2plot), yi, sei, 'Color', colors(compi,:), 'LineWidth', 3);
        else
            errorbar(FREQ.frex(idx_frex2plot), yi, sei, 'Color', colors(compi,:), 'LineWidth', 2);
        end
    end

    xlim([Landscape.frex(1) Landscape.frex(end)]);
    xticks(FREQ.frex(idx_frex2plot));
    xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Explained variance (%)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Brain network landscape - Grand-average', 'FontSize', 16, 'FontWeight', 'bold');

    legend(arrayfun(@(x) sprintf('GED component #%d', x), comps2plot, 'UniformOutput', false), ...
        'Location', 'Northeast', 'FontSize', 12);
    grid minor; set(gca, 'FontSize', 12, 'LineWidth', 1.5); box on; set(gcf, 'Color', 'w');

end

%% Output #2: NETWORKS' SPATIAL ACTIVATION PATTERNS (3D visualization)

% Update number of frequencies for spatial patterns and sort
nfrex = length(Patterns.frex);
Patterns.frex = sort(Patterns.frex,'ascend');

% Find indices of FREQ.frex that match elements in Patterns.frex
idx_frex2plot = zeros(nfrex,1);
for frexi = 1:nfrex
    [~, idx_temp] = min(abs(FREQ.frex - Patterns.frex(frexi))); % matching the closest elements to user-defined frequency
    idx_frex2plot(frexi) = idx_temp;
end

% Warning if requested freqs do not exactly match FREQ.frex
tol = 1e-3; % tolerance in Hz, adjust as needed
% Warning only if some requested freqs are NOT present in FREQ.frex (within tol)
if ~all(ismembertol(Patterns.frex(:), FREQ.frex(:), tol))
    warning('Some requested frequencies in Patterns.frex are not present in FREQ.frex. Using closest matches instead.');
end


% Plot activation patterns
col_frex = .8*parula(nfrex); % color mapping for frequencies
skipper  = 1;          % parameter for donwsampling the visualization
scale_size = 100;      % scaling factor for activation patterns in the brain
thresh_nsdt = 1; % how many std away from the mean, for thresholding the visualization

% Check the number of requested components
if Patterns.ncomps > size(FREQ.pats,2)
    error('The requested number of components exceeds the number of estimated networks!')
else

    % Pre-allocate matrix to store all patterns and compute grand-average
    nvoxs  = size(FREQ.pats,1);
    ncomps = Patterns.ncomps;
    all_pats = zeros(nvoxs,ncomps,nfrex,nsubs);

    % Plot for individual participants
    for subi = 1:nsubs
        for compi = 1:ncomps

            if plot_all || nsubs == 1
                openfig('BrainTemplate_GT.fig');
                hold on
            end

            for frexi = 1:nfrex

                % Take the absolute value and normalize 0-to-1 (deal with sign ambiguity from source reconstruction)
                this_pat = FREQ.pats(:,compi,idx_frex2plot(frexi),subi);
                this_pat = abs(this_pat);                 % deal with sign ambiguity
                this_pat = this_pat / max(this_pat);     % normalize 0–1

                % Assign to matrix to compute the grand-average later on
                all_pats(:,compi,frexi,subi) = this_pat;

                if plot_all || nsubs == 1
                    % Produce patterns to plot
                    pat2plot = squeeze(this_pat); % assign temporary variable
                    pat2plot( pat2plot < (mean(pat2plot)+thresh_nsdt*std(pat2plot)) ) = nan; % apply threshold
                    pat2plot(isnan(Patterns.MNI_coords(:,1))) = nan;

                    % Assign temporary activation pattern to plot
                    mni2plot = Patterns.MNI_coords;
                    mni2plot(isnan(pat2plot(:,1)),:) = []; % clear from nans
                    pat2plot(isnan(pat2plot)) = [];   % repeat for activation patterns

                    % Plot in 3D (magnitude mapped to transparency)
                    idx_vox     = 1:skipper:length(pat2plot);
                    mni2plot_sub = mni2plot(idx_vox,:);
                    alpha_vals  = pat2plot(idx_vox);   % already normalized 0–1
                    scatter3(mni2plot_sub(:,1), mni2plot_sub(:,2), mni2plot_sub(:,3), ...
                        scale_size, ...                     % constant marker size
                        col_frex(frexi,:), 'filled', ...    % color encodes frequency
                        'MarkerFaceAlpha','flat', ...
                        'MarkerEdgeAlpha','flat', ...
                        'AlphaData',alpha_vals, ...
                        'AlphaDataMapping','none');
                    hold on

                    title(["Network topography - Component #" num2str(compi) ' - Participant #' num2str(subi)],'FontSize',10)
                    %     legend_labels = arrayfun(@(x) sprintf('%.1f Hz', x), Patterns.frex, 'UniformOutput', false);
                    %     legend(legend_labels);
                    % Create legend with fixed-size markers matching the colors in col_frex
                    hold on;
                    legend_handles = gobjects(1, length(Patterns.frex)); % Preallocate legend handles
                    for frexi = 1:length(Patterns.frex)
                        legend_handles(frexi) = plot3(nan, nan, nan, '.', 'Color', col_frex(frexi,:), 'MarkerSize', 12);
                        % Uses NaN to avoid plotting actual points, just storing color info
                    end
                    % Adjust legend properties
                    legend(legend_handles, arrayfun(@(x) sprintf('%.1f Hz', x), FREQ.frex(idx_frex2plot), 'UniformOutput', false), ...
                        'FontSize', 14, ...  % Increase font size
                        'Location', 'northeastoutside'); % Move legend outside the plot area
                    set(legend_handles, 'MarkerSize', 20); % Increase marker size in legend
                    rotate3d on; axis off; axis vis3d; axis equal
                    set(gcf, 'Color', 'w'); % Set figure background to white
                end
            end
        end

    end

    % Group-level visualization
    if nsubs > 1
        for compi = 1:ncomps
            openfig('BrainTemplate_GT.fig');
            hold on
            for frexi = 1:nfrex

                % Compute average activation pattern to plot (NOTE: it must be already in abs values)
                avg_pats2plot = squeeze( mean( all_pats(:,compi,frexi,:), length(size(all_pats)) ) );
                avg_pats2plot = abs( avg_pats2plot / max(avg_pats2plot) );
                avg_pats2plot( avg_pats2plot < (mean(avg_pats2plot)+thresh_nsdt*std(avg_pats2plot)) ) = nan; % apply threshold
                avg_pats2plot(isnan(Patterns.MNI_coords(:,1))) = nan;

                % Assign temporary activation pattern to plot
                mni2plot = Patterns.MNI_coords;
                mni2plot(isnan(avg_pats2plot(:,1)),:) = []; % clear from nans
                avg_pats2plot(isnan(avg_pats2plot)) = [];   % repeat for activation patterns

                % Plot in 3D (magnitude mapped to transparency instead of size)
                idx_vox        = 1:skipper:length(avg_pats2plot);
                mni2plot_sub   = mni2plot(idx_vox,:);
                alpha_vals_avg = avg_pats2plot(idx_vox);   % already normalized 0–1

                scatter3(mni2plot_sub(:,1), mni2plot_sub(:,2), mni2plot_sub(:,3), ...
                    scale_size, ...                     % constant marker size
                    col_frex(frexi,:), 'filled', ...    % color encodes frequency
                    'MarkerFaceAlpha','flat', ...
                    'MarkerEdgeAlpha','flat', ...
                    'AlphaData',alpha_vals_avg, ...
                    'AlphaDataMapping','none');
                hold on

            end
            title(["Network topography - Component #" num2str(compi) ' - GRAND-AVERAGE'],'FontSize',10)
            %     legend_labels = arrayfun(@(x) sprintf('%.1f Hz', x), Patterns.frex, 'UniformOutput', false);
            %     legend(legend_labels);
            % Create legend with fixed-size markers matching the colors in col_frex
            hold on;
            legend_handles = gobjects(1, length(Patterns.frex)); % Preallocate legend handles
            for frexi = 1:length(Patterns.frex)
                legend_handles(frexi) = plot3(nan, nan, nan, '.', 'Color', col_frex(frexi,:), 'MarkerSize', 12);
                % Uses NaN to avoid plotting actual points, just storing color info
            end
            % Adjust legend properties
            legend(legend_handles, arrayfun(@(x) sprintf('%.1f Hz', x), FREQ.frex(idx_frex2plot), 'UniformOutput', false), ...
                'FontSize', 14, ...  % Increase font size
                'Location', 'northeastoutside'); % Move legend outside the plot area
            set(legend_handles, 'MarkerSize', 20); % Increase marker size in legend
            rotate3d on; axis off; axis vis3d; axis equal
            set(gcf, 'Color', 'w'); % Set figure background to white

        end
    end
end


%% Output #3: NIFTI Files (supported only for source reconstruction in

% Create directory for storing FREQNESS NIFTI output files
nifti_path = [Patterns.path_output '/FREQNESS_Output/FREQNESS_nifti'];
mkdir(nifti_path)

% load template
% template_nii = load_nii(['MNI152_T1_' num2str(Patterns.mm) 'mm_Template.nii.gz']);
template_nii     = load_nii('MNI152_8mm_brain_diy.nii.gz');
avg_template_nii = load_nii('MNI152_8mm_brain_diy.nii.gz');

% Get template image data and initialize an empty volume
nii_data = template_nii.img;
nii_data(:) = 0;  % Set all voxels to zero
nii_data = double(nii_data);
% Repeat for average
avg_nii_data = avg_template_nii.img;
avg_nii_data(:) = 0;  % Set all voxels to zero
avg_nii_data = double(avg_nii_data);

% Extract affine transformation matrix from srow_x, srow_y, srow_z
affine = [template_nii.hdr.hist.srow_x;
    template_nii.hdr.hist.srow_y;
    template_nii.hdr.hist.srow_z;
    0 0 0 1]; % Append [0 0 0 1] to make it 4x4

% extract MNI coordinates
MNI_coords = Patterns.MNI_coords;

% Produce NIFTI files
for frexi = 1:nfrex

    for compi = 1:Patterns.ncomps  % Iterate over components given as input

        % Compute average activation pattern to plot (NOTE: it must be already in abs values)
        avg_SO = squeeze( mean( all_pats(:,compi,frexi,:), length(size(all_pats)) ) );

        for subi = 1:nsubs

            % Extract spatial activation pattern (per participant)
            SO = all_pats(:,compi,frexi,subi);

            num_points = size(MNI_coords, 1);
            voxel_coords = zeros(num_points, 3);

            % Convert MNI Coordinates to Voxel Indices
            for ii = 1:num_points
                coord = [MNI_coords(ii, :) 1];  % Add homogeneous coordinate
                voxel = affine\coord';% inv(affine) * coord';  % Convert to voxel space
                voxel_coords(ii, :) = (voxel(1:3)); % Extract rounded voxel indices

                % Adjust for voxel center vs edge (half voxel shift)
                voxel_coords(ii, :) = voxel_coords(ii, :) + 1;  % Subtract 1 voxel (adjust for 8mm shift)

            end

            % Assign activation patterns to the corresponding voxels
            for ii = 1:num_points
                x = voxel_coords(ii, 1);
                y = voxel_coords(ii, 2);
                z = voxel_coords(ii, 3);

                % Ensure indices are within image boundaries
                if all([x, y, z] > 0) && all([x, y, z] <= size(nii_data))
                    nii_data(x, y, z) = SO(ii);
                end
                % Repeat for average
                if all([x, y, z] > 0) && all([x, y, z] <= size(avg_nii_data))
                    avg_nii_data(x, y, z) = avg_SO(ii);
                end

            end

            % Create nii template from individual participants' data
            template_nii.img = nii_data;

            % Create a NIFTI image from the 3D data matrix (8 mm resolution)
            nii = make_nii(nii_data, [8 8 8]);
            nii.img = nii_data;  % Store matrix within image structure
            nii.hdr.hist = template_nii.hdr.hist;  % Copy header information from mask

            % Display saving progress
            disp(['Saving NIFTI images - ' num2str(FREQ.frex(idx_frex2plot(frexi))) ' Hz networks - Comp #' num2str(compi) '_Sub#' num2str(subi)])

            % Save the NIFTI file
            save_nii(nii, [nifti_path '/ActivationPattern_Frex_' num2str(FREQ.frex(idx_frex2plot(frexi))) 'Hz_Comp#' num2str(compi) '_Sub#' num2str(subi) '.nii']);

        end

        % Replicate for grand average
        if nsubs > 1

            % Create nii template
            avg_template_nii.img = avg_nii_data;
            % Create a NIFTI image from the 3D data matrix (8 mm resolution)
            avg_nii = make_nii(avg_nii_data, [8 8 8]);
            avg_nii.img = avg_nii_data;  % Store matrix within image structure
            avg_nii.hdr.hist = avg_template_nii.hdr.hist;  % Copy header information from mask

            % Display saving progress
            disp(['Saving NIFTI images - ' num2str(FREQ.frex(idx_frex2plot(frexi))) ' Hz networks - Comp #' num2str(compi) ' - GRAND-AVERAGE'])

            % Save the NIFTI file
            save_nii(avg_nii, [nifti_path '/ActivationPattern_Frex_' num2str(FREQ.frex(idx_frex2plot(frexi))) 'Hz_Comp#' num2str(compi) '_Sub#' num2str(subi) ' - GRAND-AVERAGE.nii']);

        end

    end

end


end



