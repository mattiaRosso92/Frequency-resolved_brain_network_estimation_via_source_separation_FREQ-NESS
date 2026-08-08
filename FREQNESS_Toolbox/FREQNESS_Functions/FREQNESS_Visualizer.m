function FREQNESS_Visualizer(FREQ,Landscape,Patterns,varargin)

% ========================================================================
%  FREQNESS VISUALIZER: NETWORK LANDSCAPE & BRAIN ACTIVATION PATTERNS
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
%
%  This function visualizes:
%
%  - Network prominence across frequencies
%    Quantified as % variance explained based on the associated
%    eigenvalues, producing a frequency-resolved network landscape
%    (see Rosso et al., 2025, Advanced Science, Figure 2).
%
%  - Spatial activation patterns of the networks
%    - All valid MNI coordinates are shown as visible black reference dots.
%    - Requested frequencies are encoded by a low-to-high viridis gradient.
%    - Normalized activation magnitude is encoded by marker size.
%    - NIFTI files can be generated for each requested network.
%
%  If multiple participants are included in the input FREQ structure, the
%  function will by default plot only grand-average landscapes and
%  activation patterns. Individual participants can be plotted by setting
%  the optional flag 'plot_all' to true.
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - FREQ              : Structure outputted by FREQNESS_NetworkEstimation.
%      - FREQ.evals    : Normalized eigenvalues (variance explained in % points).
%      - FREQ.pats     : Spatial activation patterns.
%      - FREQ.frex     : Vector of frequencies analyzed.
%
%  - Landscape                : Structure containing settings for the
%                               network landscape visualization.
%      - Landscape.frex       : Frequencies to visualize.
%      - Landscape.ncomps     : Number of components to visualize.
%
%  - Patterns                 : Structure containing settings for activation
%                               pattern visualization. Pass [] to generate
%                               only the network landscape.
%      - Patterns.MNI_coords  : MNI coordinates in the same order as the data
%                               [nVoxels x 3].
%      - Patterns.frex        : Frequencies to visualize.
%      - Patterns.ncomps      : Number of components to visualize.
%      - Patterns.path_output : Output path for NIFTI images. Required only
%                               when 'save_nifti' is true.
%
% ------------------------------------------------------------------------
%  OPTIONAL NAME-VALUE PAIRS:
% ------------------------------------------------------------------------
%  - 'plot_all'    : Plot every participant in addition to the group
%                    average. Default: false.
%
%  - 'threshold_sd': Activation threshold expressed as standard deviations
%                    above the mean normalized pattern. Default: 1.
%
%  - 'save_nifti'  : Save individual and grand-average NIFTI maps.
%                    Default: true.
%
%  NOTE 1: the 3D scatter plot is supported for any MNI source space.
%  NOTE 2: NIFTI output is currently supported only for the bundled
%          MNI152 8-mm template.
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

if nargin < 3
    Patterns = [];
end


%% Handle optional arguments

opts = struct('plot_all',false,'threshold_sd',1,'save_nifti',true);
opts = parse_name_value_pairs(opts,varargin{:});

plot_all    = opts.plot_all;
threshold_sd = opts.threshold_sd;
save_nifti  = opts.save_nifti;

if ~islogical(plot_all) || ~isscalar(plot_all)
    error('plot_all must be one logical value (true or false).');
end

if ~isnumeric(threshold_sd) || ~isreal(threshold_sd) || ...
        ~isscalar(threshold_sd) || ~isfinite(threshold_sd) || threshold_sd < 0
    error('threshold_sd must be one non-negative finite scalar.');
end

if ~islogical(save_nifti) || ~isscalar(save_nifti)
    error('save_nifti must be one logical value (true or false).');
end


%% Check FREQ and Landscape inputs

if nargin < 2 || ~isstruct(FREQ) || ~isstruct(Landscape)
    error('FREQ and Landscape must be provided as structures.');
end

required_freq_fields = {'evals','frex'};
for fieldi = 1:numel(required_freq_fields)
    if ~isfield(FREQ,required_freq_fields{fieldi}) || isempty(FREQ.(required_freq_fields{fieldi}))
        error('FREQ.%s is missing or empty.',required_freq_fields{fieldi});
    end
end

if ~isnumeric(FREQ.evals) || ~isreal(FREQ.evals) || ...
        any(~isfinite(FREQ.evals(:))) || ndims(FREQ.evals) > 3
    error(['FREQ.evals must be a real finite numeric array in ' ...
        '[nComp x nFrex x (nSubs)] format.']);
end

frex = FREQ.frex(:)';
if ~isnumeric(frex) || ~isreal(frex) || isempty(frex) || ...
        any(~isfinite(frex)) || any(frex <= 0)
    error('FREQ.frex must contain positive finite frequencies.');
end

if size(FREQ.evals,2) ~= numel(frex)
    error('The frequency dimension of FREQ.evals does not match FREQ.frex.');
end

if ~isfield(Landscape,'frex') || ~isfield(Landscape,'ncomps')
    error('Landscape must contain the fields frex and ncomps.');
end

landscape_frex = validate_frequency_vector(Landscape.frex,'Landscape.frex');
landscape_ncomps = validate_component_number(Landscape.ncomps,'Landscape.ncomps');

if landscape_ncomps > size(FREQ.evals,1)
    error('Landscape.ncomps exceeds the number of estimated networks.');
end

[idx_landscape,landscape_frex_actual] = map_frequencies( ...
    frex,landscape_frex,'Landscape');

nsubs = size(FREQ.evals,3);
fprintf('\nFREQNESS Visualizer: visualizing %d participants.\n',nsubs);


%% Output #1: Network landscape

comps2plot = 1:landscape_ncomps;
colors = .8*parula(landscape_ncomps);
max_eval = max(FREQ.evals(:));
if max_eval <= 0
    max_eval = 1;
end

% Individual participants
if plot_all || nsubs == 1
    for subi = 1:nsubs
        figure('Color','w'); hold on

        for compi = comps2plot
            y = squeeze(FREQ.evals(compi,idx_landscape,subi));
            y = y(:)';
            linewidth = 2;
            if compi == 1
                linewidth = 3;
            end
            plot(landscape_frex_actual,y,'Color',colors(compi,:), ...
                'LineWidth',linewidth);
        end

        set_frequency_limits(landscape_frex_actual)
        ylim([0 max_eval*1.1])
        xticks(landscape_frex_actual)
        xlabel('Frequency (Hz)','FontSize',14,'FontWeight','bold');
        ylabel('Explained variance (%)','FontSize',14,'FontWeight','bold');
        title(['Brain network landscape - Participant #' num2str(subi)], ...
            'FontSize',16,'FontWeight','bold');
        legend(arrayfun(@(x) sprintf('GED component #%d',x),comps2plot, ...
            'UniformOutput',false),'Location','Northeast','FontSize',12);
        grid minor
        set(gca,'FontSize',12,'LineWidth',1.5)
        box on
    end
end

% Grand average
if nsubs > 1
    figure('Color','w'); hold on

    for compi = comps2plot
        this_evals = squeeze(FREQ.evals(compi,idx_landscape,:));
        this_evals = reshape(this_evals,numel(idx_landscape),nsubs);
        avg_evals = mean(this_evals,2,'omitnan');
        sem_evals = std(this_evals,0,2,'omitnan')/sqrt(nsubs);
        linewidth = 2;
        if compi == 1
            linewidth = 3;
        end
        errorbar(landscape_frex_actual,avg_evals,sem_evals, ...
            'Color',colors(compi,:),'LineWidth',linewidth);
    end

    set_frequency_limits(landscape_frex_actual)
    ylim([0 max_eval*1.1])
    xticks(landscape_frex_actual)
    xlabel('Frequency (Hz)','FontSize',14,'FontWeight','bold');
    ylabel('Explained variance (%)','FontSize',14,'FontWeight','bold');
    title('Brain network landscape - Grand-average', ...
        'FontSize',16,'FontWeight','bold');
    legend(arrayfun(@(x) sprintf('GED component #%d',x),comps2plot, ...
        'UniformOutput',false),'Location','Northeast','FontSize',12);
    grid minor
    set(gca,'FontSize',12,'LineWidth',1.5)
    box on
end

% Landscape-only mode
if isempty(Patterns)
    return
end


%% Check Patterns inputs

if ~isstruct(Patterns)
    error('Patterns must be a structure or an empty array.');
end

required_pattern_fields = {'MNI_coords','frex','ncomps'};
for fieldi = 1:numel(required_pattern_fields)
    if ~isfield(Patterns,required_pattern_fields{fieldi}) || isempty(Patterns.(required_pattern_fields{fieldi}))
        error('Patterns.%s is missing or empty.',required_pattern_fields{fieldi});
    end
end

if ~isfield(FREQ,'pats') || isempty(FREQ.pats) || ...
        ~isnumeric(FREQ.pats) || ~isreal(FREQ.pats) || ...
        any(isinf(FREQ.pats(:))) || ndims(FREQ.pats) > 4
    error(['FREQ.pats must be a real numeric array in ' ...
        '[nVoxels x nComp x nFrex x (nSubs)] format.']);
end

if size(FREQ.pats,2) ~= size(FREQ.evals,1) || ...
        size(FREQ.pats,3) ~= numel(frex) || size(FREQ.pats,4) ~= nsubs
    error('FREQ.evals and FREQ.pats have incompatible dimensions.');
end

MNI_coords = Patterns.MNI_coords;
if ~isnumeric(MNI_coords) || ~isreal(MNI_coords) || ...
        ndims(MNI_coords) ~= 2 || size(MNI_coords,2) ~= 3 || ...
        size(MNI_coords,1) ~= size(FREQ.pats,1) || any(isinf(MNI_coords(:)))
    error(['Patterns.MNI_coords must be a real [nVoxels x 3] numeric matrix ' ...
        'with one row per voxel in FREQ.pats.']);
end

if ~any(all(isfinite(MNI_coords),2))
    error('Patterns.MNI_coords contains no valid coordinates.');
end

pattern_frex = validate_frequency_vector(Patterns.frex,'Patterns.frex');
pattern_ncomps = validate_component_number(Patterns.ncomps,'Patterns.ncomps');
if pattern_ncomps > size(FREQ.pats,2)
    error('Patterns.ncomps exceeds the number of estimated networks.');
end

[idx_patterns,pattern_frex_actual] = map_frequencies( ...
    frex,pattern_frex,'Patterns');

if save_nifti && (~isfield(Patterns,'path_output') || isempty(Patterns.path_output) || ...
        ~(ischar(Patterns.path_output) || ...
        (isstring(Patterns.path_output) && isscalar(Patterns.path_output))))
    error('Patterns.path_output must be provided when save_nifti is true.');
end


%% Output #2: Spatial activation patterns

nvoxs = size(FREQ.pats,1);
nfrex_patterns = numel(idx_patterns);
all_pats = nan(nvoxs,pattern_ncomps,nfrex_patterns,nsubs);

% Normalize every participant pattern independently
for subi = 1:nsubs
    for compi = 1:pattern_ncomps
        for frexi = 1:nfrex_patterns
            this_pat = abs(FREQ.pats(:,compi,idx_patterns(frexi),subi));
            valid_pat = isfinite(this_pat);
            this_pat(~valid_pat) = nan;
            if any(valid_pat)
                max_pat = max(this_pat(valid_pat));
                if max_pat > 0
                    this_pat(valid_pat) = this_pat(valid_pat)/max_pat;
                else
                    this_pat(valid_pat) = 0;
                end
            end
            all_pats(:,compi,frexi,subi) = this_pat;
        end
    end
end

% Individual participant patterns
if plot_all || nsubs == 1
    for subi = 1:nsubs
        for compi = 1:pattern_ncomps
            title_text = ['Network topography - Component #' num2str(compi) ...
                ' - Participant #' num2str(subi)];
            plot_activation_patterns(MNI_coords, ...
                all_pats(:,compi,:,subi),pattern_frex_actual, ...
                threshold_sd,title_text)
        end
    end
end

% Grand-average patterns
group_pats = [];
if nsubs > 1
    group_pats = mean(all_pats,4,'omitnan');
    for compi = 1:pattern_ncomps
        for frexi = 1:nfrex_patterns
            this_pat = group_pats(:,compi,frexi);
            valid_pat = isfinite(this_pat);
            if any(valid_pat)
                max_pat = max(this_pat(valid_pat));
                if max_pat > 0
                    this_pat(valid_pat) = this_pat(valid_pat)/max_pat;
                else
                    this_pat(valid_pat) = 0;
                end
            end
            group_pats(:,compi,frexi) = this_pat;
        end

        title_text = ['Network topography - Component #' num2str(compi) ...
            ' - GRAND-AVERAGE'];
        plot_activation_patterns(MNI_coords,group_pats(:,compi,:), ...
            pattern_frex_actual,threshold_sd,title_text)
    end
end


%% Output #3: NIFTI files

if save_nifti
    write_nifti_patterns(MNI_coords,all_pats,group_pats, ...
        pattern_frex_actual,Patterns.path_output)
end

end


%% Helper function: Plot activation patterns
function plot_activation_patterns(MNI_coords,patterns,frex,threshold_sd,title_text)

valid_mni = all(isfinite(MNI_coords),2);
mni_brain = MNI_coords(valid_mni,:);

fig = figure('Color','w','Units','pixels','Position',[100 100 900 700], ...
    'Renderer','opengl');
ax = axes('Parent',fig);
hold(ax,'on')

% Anatomical reference: all valid MNI coordinates as visible black dots
scatter3(ax,mni_brain(:,1),mni_brain(:,2),mni_brain(:,3),12,[0 0 0], ...
    'filled','MarkerFaceAlpha',0.35,'MarkerEdgeColor','none');

% Frequency gradient
cmap = viridis_colormap(256);
if numel(frex) == 1
    color_idx = round(size(cmap,1)/2);
    color_limits = [frex(1)-max(abs(frex(1))*0.05,0.5) ...
        frex(1)+max(abs(frex(1))*0.05,0.5)];
else
    color_idx = 1 + round((frex-min(frex))/(max(frex)-min(frex))*(size(cmap,1)-1));
    color_limits = [min(frex) max(frex)];
end

% Active voxels: frequency is color and normalized magnitude is marker size
active_mni = [];
active_values = [];
active_colors = [];
for frexi = 1:numel(frex)
    this_pat = squeeze(patterns(:,1,frexi));
    valid_pat = isfinite(this_pat);
    if ~any(valid_pat)
        continue
    end
    cutoff = mean(this_pat(valid_pat)) + threshold_sd*std(this_pat(valid_pat));
    idx_active = valid_mni & valid_pat & this_pat >= cutoff & this_pat > 0;

    if any(idx_active)
        this_values = this_pat(idx_active);
        this_mni = MNI_coords(idx_active,:);
        this_colors = repmat(cmap(color_idx(frexi),:),numel(this_values),1);
        active_values = [active_values; this_values]; %#ok<AGROW>
        active_mni = [active_mni; this_mni]; %#ok<AGROW>
        active_colors = [active_colors; this_colors]; %#ok<AGROW>
    end
end

% Plot all frequencies together, ordering points only by activation magnitude
if ~isempty(active_values)
    [active_values,sort_idx] = sort(active_values,'ascend');
    active_mni = active_mni(sort_idx,:);
    active_colors = active_colors(sort_idx,:);
    marker_size = 20 + 100*active_values;

    scatter3(ax,active_mni(:,1),active_mni(:,2),active_mni(:,3), ...
        marker_size,active_colors,'filled', ...
        'MarkerFaceAlpha',0.9,'MarkerEdgeColor','none');
end

colormap(ax,cmap)
clim(ax,color_limits)
cb = colorbar(ax);
cb.Label.String = 'Frequency (Hz)';
cb.Label.FontWeight = 'bold';
cb.Ticks = unique(frex);
cb.FontSize = 11;
cb.Position = [0.87 0.18 0.025 0.64];

% Match Matplotlib's centered equal-scale MNI-space projection
coordinate_min = min(mni_brain,[],1);
coordinate_max = max(mni_brain,[],1);
coordinate_span = coordinate_max-coordinate_min;
coordinate_span(coordinate_span == 0) = 1;
coordinate_padding = 0.05*coordinate_span;
xlim(ax,[coordinate_min(1)-coordinate_padding(1) ...
    coordinate_max(1)+coordinate_padding(1)])
ylim(ax,[coordinate_min(2)-coordinate_padding(2) ...
    coordinate_max(2)+coordinate_padding(2)])
zlim(ax,[coordinate_min(3)-coordinate_padding(3) ...
    coordinate_max(3)+coordinate_padding(3)])
daspect(ax,[1 1 1])
pbaspect(ax,coordinate_span)
axis(ax,'vis3d')
axis(ax,'off')
view(ax,135,25)
camproj(ax,'perspective')
camtarget(ax,(coordinate_min+coordinate_max)/2)
camzoom(ax,0.60)
ax.Position = [0.43 0.10 0.60 0.80];
title(ax,title_text,'FontSize',14,'FontWeight','bold','Visible','on')
rotate3d(fig,'on')

end


%% Helper function: Viridis frequency colormap
function cmap = viridis_colormap(ncolors)

anchor_colors = [ ...
    0.267004 0.004874 0.329415; ...
    0.281412 0.155834 0.469201; ...
    0.244972 0.287675 0.537260; ...
    0.190631 0.407061 0.556089; ...
    0.147607 0.511733 0.557049; ...
    0.119699 0.618490 0.536347; ...
    0.208030 0.718701 0.472873; ...
    0.430983 0.808473 0.346476; ...
    0.709898 0.868751 0.169257; ...
    0.993248 0.906157 0.143936];

anchor_positions = linspace(0,1,size(anchor_colors,1));
color_positions = linspace(0,1,ncolors);
cmap = interp1(anchor_positions,anchor_colors,color_positions,'linear');

end


%% Helper function: Write NIFTI patterns
function write_nifti_patterns(MNI_coords,all_pats,group_pats,frex,path_output)

function_path = fileparts(mfilename('fullpath'));
toolbox_path = fileparts(function_path);
template_path = fullfile(toolbox_path,'FREQNESS_ExternalFunctions', ...
    'nifti_tools','MNI152_8mm_brain_diy.nii.gz');

if ~isfile(template_path)
    error('NIFTI template not found: %s',template_path);
end

nifti_path = fullfile(char(path_output),'FREQNESS_Output','FREQNESS_nifti');
if ~exist(nifti_path,'dir')
    mkdir(nifti_path)
end

template_nii = load_nii(template_path);
image_size = size(template_nii.img);
affine = [template_nii.hdr.hist.srow_x; ...
    template_nii.hdr.hist.srow_y; ...
    template_nii.hdr.hist.srow_z; ...
    0 0 0 1];

% Convert MNI coordinates to one-based MATLAB voxel indices
homogeneous_coords = [MNI_coords ones(size(MNI_coords,1),1)]';
voxel_coords = affine\homogeneous_coords;
voxel_coords = round(voxel_coords(1:3,:)') + 1;
valid_voxels = all(isfinite(MNI_coords),2) & ...
    all(voxel_coords >= 1,2) & ...
    voxel_coords(:,1) <= image_size(1) & ...
    voxel_coords(:,2) <= image_size(2) & ...
    voxel_coords(:,3) <= image_size(3);

nsubs = size(all_pats,4);
ncomps = size(all_pats,2);

for frexi = 1:numel(frex)
    frequency_label = num2str(frex(frexi),'%g');

    for compi = 1:ncomps
        for subi = 1:nsubs
            this_pat = all_pats(:,compi,frexi,subi);
            idx2write = valid_voxels & isfinite(this_pat);
            this_voxels = voxel_coords(idx2write,:);

            nii_data = zeros(image_size);
            linear_idx = sub2ind(image_size,this_voxels(:,1), ...
                this_voxels(:,2),this_voxels(:,3));
            nii_data(linear_idx) = this_pat(idx2write);

            nii = make_nii(nii_data,[8 8 8]);
            nii.hdr.hist = template_nii.hdr.hist;
            output_name = fullfile(nifti_path, ...
                ['ActivationPattern_Frex_' frequency_label 'Hz_Comp#' ...
                num2str(compi) '_Sub#' num2str(subi) '.nii']);
            save_nii(nii,output_name)
            disp(['Saving NIFTI image: ' output_name])
        end

        if nsubs > 1
            this_pat = group_pats(:,compi,frexi);
            idx2write = valid_voxels & isfinite(this_pat);
            this_voxels = voxel_coords(idx2write,:);

            nii_data = zeros(image_size);
            linear_idx = sub2ind(image_size,this_voxels(:,1), ...
                this_voxels(:,2),this_voxels(:,3));
            nii_data(linear_idx) = this_pat(idx2write);

            nii = make_nii(nii_data,[8 8 8]);
            nii.hdr.hist = template_nii.hdr.hist;
            output_name = fullfile(nifti_path, ...
                ['ActivationPattern_Frex_' frequency_label 'Hz_Comp#' ...
                num2str(compi) '_GRAND-AVERAGE.nii']);
            save_nii(nii,output_name)
            disp(['Saving NIFTI image: ' output_name])
        end
    end
end

end


%% Helper function: Validate frequency vector
function frex = validate_frequency_vector(frex,name)

if ~isnumeric(frex) || ~isreal(frex) || ~isvector(frex) || isempty(frex) || ...
        any(~isfinite(frex)) || any(frex <= 0)
    error('%s must be a non-empty vector of positive finite frequencies.',name);
end
frex = sort(frex(:)');

end


%% Helper function: Validate component number
function ncomps = validate_component_number(ncomps,name)

if ~isnumeric(ncomps) || ~isreal(ncomps) || ~isscalar(ncomps) || ...
        ~isfinite(ncomps) || ncomps ~= round(ncomps) || ncomps < 1
    error('%s must be one positive integer.',name);
end

end


%% Helper function: Map frequencies
function [idx_frex,frex_actual] = map_frequencies(frex,frex_requested,label)

idx_frex = zeros(1,numel(frex_requested));
for frexi = 1:numel(frex_requested)
    [~,idx_frex(frexi)] = min(abs(frex-frex_requested(frexi)));
end
frex_actual = frex(idx_frex);

if ~all(ismembertol(frex_requested(:),frex(:),1e-3))
    warning(['Some requested frequencies in ' label ...
        '.frex are not present in FREQ.frex. Using closest matches instead.']);
end

end


%% Helper function: Set frequency limits
function set_frequency_limits(frex)

if numel(frex) == 1
    padding = max(abs(frex(1))*0.05,0.5);
    xlim([frex(1)-padding frex(1)+padding])
else
    xlim([min(frex) max(frex)])
end

end


%% Helper function: Parse name-value pairs
function opts = parse_name_value_pairs(opts,varargin)

if mod(length(varargin),2) ~= 0
    error('Arguments must be given as name-value pairs.');
end

for i = 1:2:length(varargin)
    if ~(ischar(varargin{i}) || (isstring(varargin{i}) && isscalar(varargin{i})))
        error('Optional argument names must be text.');
    end
    name = lower(char(varargin{i}));
    if isfield(opts,name)
        opts.(name) = varargin{i+1};
    else
        error(['Unrecognized argument: ' name]);
    end
end

end
