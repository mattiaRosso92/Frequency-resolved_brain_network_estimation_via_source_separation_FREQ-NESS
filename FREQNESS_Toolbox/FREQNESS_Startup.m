function FREQNESS_Startup(path_home)
% ========================================================================
%  FREQNESS STARTUP FUNCTION
% ========================================================================
%
%  Please cite the first FREQNESS paper:
%  Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P., Kringelbach, M. L., & Bonetti, L. (2025). 
%  FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain Networks During Auditory Stimulation. 
%  Advanced Science, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
% ========================================================================
%
%  This function initializes the FREQNESS toolbox.
%
%  The user must provide `path_home` as a string, keeping the original 
%  subfolder structure of the FREQNESS toolbox unchanged.
%
%  Example usage:
%      FREQNESS_Startup('/Users/mattiarosso/Desktop/FREQNESS_Toolbox/')
%
%  The function will add the relevant subfolders to the user's current MATLAB path.


% ------------------------------------------------------------------------
%  AUTHORS:
%  Mattia Rosso & Leonardo Bonetti
%  mattia.rosso@clin.au.dk
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  Aarhus (DK), Oxford (UK), 24/02/2025
% ========================================================================


% Add main toolbox directories and NIFTI tools subfolder to MATLAB path
addpath(path_home);
addpath([path_home '/FREQNESS_Functions/'])
addpath([path_home '/FREQNESS_ExternalFunctions/'])
addpath([path_home '/FREQNESS_ExternalFunctions/nifti_tools']);


end
