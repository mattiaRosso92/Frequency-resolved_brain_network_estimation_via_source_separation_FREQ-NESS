function O = GED_SingleSubjects_v3(S)

O = []; 

% Settings (assign from input structure)
srate      = S.srate;       % recorded at 1000 Hz, downsampled by 4
targetfrex = S.targetfrex;  % stimulation frequency: hard-coded 1 / 0.410 s (range is 0.400 and 0.420 but 0.410 is the mode) 
fwhm       = S.fwhm;        % filter width
shr        = S.shr;         % value between 0 and 1 (0 small perturbation, 1 big perturbation) - for Mattia regularisation (until 0.001 it works for restoring the full rank of the matrix)
ncomps     = S.ncomps2keep; % how many components to keep (given the huge amount of components = N voxels)
comps2see  = S.comps2see;   % list components of interest, of which you want to save nifti image
% For permutation testing
permu_flag = S.permu_flag;
permu_name = S.permu_name;
% For data 
path_script = S.path_script;    % path to the script calling this function
path_input  = S.path_source;    % to import source data ;
path_output_pre = S.path_GEDoutput; % to insert folder for GED data, per frequency 
if S.permu_flag > 0
    name_folder = [S.name_folder '_' permu_name{permu_flag} '/'];
else
    name_folder = [S.name_folder '/'];
end

% Create the new folder in the desired output path
mkdir([path_output_pre name_folder])
path_output = [path_output_pre name_folder]; % this is where we save GED output for all SUBJs  

% Prepare subject list
addpath(path_script) % path of script calling this function 
cd(path_input);
list = dir ([path_input 'SUBJ*']);

% Compute covariance matrices
for subi = 1:length(list)
  
    disp(['Loading data and computing covariance of ' list(subi).name])   
    
    % Import broadband source data
    load([path_input list(subi).name])
    broadData = OUT.sources_ERFs(:,:);
    
    % Filter data
    narrowData = filterFGx(broadData,srate,targetfrex,fwhm,0); % turn vis_flag to 0 for suppressing the plot 

    % Permutate broad and narrowband data, if requested
    if permu_flag == 1
        % Randomize labels
        n2permute  = size(broadData,1);     % labels (voxels) | rows
        idx_perm   = randperm(n2permute);
        broadData  = broadData(idx_perm,:); % apply same indexes to broad and narrow
        narrowData = narrowData(idx_perm,:);
        
    elseif permu_flag == 2
        % Randomize labels (independently for each timepoint)
        n2permute  = size(broadData,1);     % labels (voxels) | rows
        for timi = 1:size(broadData,2) % loop over timepoints
            idx_perm   = randperm(n2permute); % new randomization for every timepoint  
            broadData(:,timi)  = broadData(idx_perm,timi); % apply same indexes to broad and narrow
            narrowData(:,timi) = narrowData(idx_perm,timi);
        end
        
    elseif permu_flag == 3
        % Randomize timepoints
        n2permute  = size(broadData,2);     % timepoints | columns
        idx_perm   = randperm(n2permute);
        broadData  = broadData(:,idx_perm); % apply same indexes to broad and narrow
        narrowData = narrowData(:,idx_perm);
        
    elseif permu_flag == 4
        % Randomize all
        n2permute  = length(broadData(:));  % all elements
        idx_perm   = randperm(n2permute);
        broadData  = broadData(:);          % reshape to 1-D
        broadData  = broadData(idx_perm);   % shuffle
        broadData  = reshape(broadData, size(narrowData)); % re-shape back to 2D       
        narrowData = narrowData(:); % repeat for narrowData
        narrowData = narrowData(idx_perm);
        narrowData  = reshape(narrowData, size(broadData));
                
    end
    
    % Compute covariance matrices
    covS = cov(narrowData');
    covR = cov(broadData');

    % Regularisation (to bring covariance matrices to full rank)
    % narrow data
    covS = covS  + 1e-6*eye(size(covS )); %making matrix full rank again by adding a small perturbation/noise (matrix regularisation)
    evalsR = eig(covR );  % for regularization of R covariance matrix
    covR  = (1-shr)*covR  + shr * mean(evalsR) * eye(size(covR )); % regularization 

% Visualize covariance matrices
% figure, clf
% subplot(121)
% imagesc(covS_avg ), axis square
% title('Covariance S (Narrowband)')
% colorbar
% subplot(122)
% imagesc(covR_avg ), axis square
% title('Covariance R (Broadband)')
% colorbar

    % Eigendecomposition 
    [evecs,evals] = eig(covS ,covR);
    [evals,sidx]  = sort( diag(evals),'descend' ); % first output returns sorted evals extracted from diagonal
    evecs = evecs(:,sidx);          % sort eigenvectors
    evals = evals.*100./sum(evals); % normalize eigenvalues to percent variance explained   
    % compute filter forward model and flip sign
    GEDmap = zeros(ncomps,size(broadData,1));
    GEDts  = zeros(ncomps,size(broadData,2));
    disp('Computing time series')
    for compi = 1:ncomps
        GEDmap(compi,:) = evecs(:,compi)'*covS ; % get component
        [~,idxmax] = max(abs(GEDmap(:,1)));     % find max magnitude
        GEDmap(compi,:)  = GEDmap(compi,:)*sign(GEDmap(compi,idxmax)); % possible sign flip
        % compute time series (projections)
        GEDts(compi,:) = evecs(:,compi)'*broadData;
        disp(['Component #' num2str(compi) ' has been computed.'])
    end

    % Save output
    save([path_output, list(subi).name(1:9), '.mat'], 'covS','covR','evecs','evals','GEDmap','GEDts','-v7.3')  


    %plotting weights in the brain (nifti images)
    maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
    for compi = 1:length(comps2see) %over significant PCs
        %using nifti toolbox
        SO = GEDmap(compi,:)';
        %building nifti image
        SS = size(maskk.img);
        dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
        for sourci = 1:size(SO,1) %over brain sources
            dumm = find(maskk.img == sourci); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
            [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
            dumimg(i1,i2,i3,:) = SO(sourci,:); %storing values for all time-points in the image matrix
        end
        nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
        nii.img = dumimg; %storing matrix within image structure
        nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
        disp(['saving nifti images - subject ' list(subi).name(1:9)])
        save_nii(nii,[path_output list(subi).name(1:9) '_Comp_' num2str(compi) '_Var_' num2str(round(evals(compi))) '.nii']); %printing image
    end

end

end

