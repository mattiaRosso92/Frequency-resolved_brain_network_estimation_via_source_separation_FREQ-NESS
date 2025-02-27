This repository contains two folders:
-: Code and functions used for the FREQNESS paper published in Advanced Science (Rosso et al., 2025)
-: Matlab FREQNESS toolbox (additional information is provided at the bottom of this page)

FREQ-NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved Brain Networks During Auditory Stimulation

Matlab leading script and functions for the paper entitled: 
FREQ-NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved Brain Networks During Auditory Stimulation.
M. Rosso, G. Fernandez-Rubio, P. E. Keller, E. Brattico, P. Vuust, M. L. Kringelbach, L. Bonetti.
Adv. Sci. 2025, 2413195.
https://doi.org/10.1002/advs.202413195

Additional relevant codes and functions are available here: https://github.com/leonardob92/LBPD-1.0.git

Abstract: The brain is a dynamic system whose network organization is often studied by focusing on specific frequency bands or anatomical regions, leading to fragmented insights, or by employing complex and elaborate methods that hinder straightforward interpretations. To address this issue, a new analytical pipeline named FREQuency-resolved Network Estimation via Source Separation (FREQ-NESS) is introduced. This is designed to estimate the activation and spatial configuration of simultaneous brain networks across frequencies by analyzing the frequency-resolved multivariate covariance between whole-brain voxel time series. FREQ-NESS is applied to source-reconstructed magnetoencephalography (MEG) data during resting state and isochronous auditory stimulation. Results reveal simultaneous, frequency-specific brain networks during resting state, such as the default mode, alpha-band, and motor-beta networks. During auditory stimulation, FREQ-NESS detects: (1) emergence of networks attuned to the stimulation frequency, (2) spatial reorganization of existing networks, such as alpha-band networks shifting from occipital to sensorimotor areas, (3) stability of networks unaffected by auditory stimuli. Furthermore, auditory stimulation significantly enhances cross-frequency coupling, with the phase of attuned auditory networks modulating the gamma band amplitude of medial temporal lobe networks. In conclusion, FREQ-NESS effectively maps the brain’s spatiotemporal dynamics, providing a comprehensive view of brain function by revealing simultaneous, frequency-resolved networks and their interaction.


The conceptualization and code implementation of this work were carried out in collaboration with Prof. Leonardo Bonetti (https://github.com/leonardob92).

Corresponding authors:
Mattia Rosso     - mattia.rosso@clin.au.dk
Leonardo Bonetti - leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk



FREQNESS TOOLBOX

% ========================================================================
%  FREQNESS: EXAMPLE SCRIPT FOR BRAIN NETWORK ESTIMATION
%
%  Please cite the first FREQNESS paper:
%  M. Rosso, G. Fernandez-Rubio, P. E. Keller, E. Brattico, P. Vuust, M. L. Kringelbach, L. Bonetti.
%  FREQ-NESS Reveals the Dynamic Reconfiguration of Frequency-Resolved Brain Networks During Auditory Stimulation.
%  Adv. Sci. 2025, 2413195.
%  https://doi.org/10.1002/advs.202413195
%
%
% ========================================================================
%
%  This script demonstrates the application of FREQNESS 
%  (Frequency-Resolved Brain Network Estimation via Source Separation).
%
%  The example dataset consists of a 3-minute recording with 
%  two conditions:
%  - Resting state ('rest')
%  - Passive listening to an isochronous metronome at 2.4 Hz ('beat')
%  This data is available at the following link (Zenodo repository):
%  https://zenodo.org/records/14922536
%  We recommend to place the data in the folder 'Data_Example' which can be
%  found in the FREQNESS_Toolbox main folder.  
%
%  The script loads the example datasets into the MATLAB workspace and 
%  calls three core functions:
%
% ------------------------------------------------------------------------
%  FUNCTIONS OVERVIEW:
% ------------------------------------------------------------------------
%  - FREQNESS_Startup(path_home) 
%    Initializes the environment. Given `path_home` as input, it sets up 
%    the necessary directories.
%
%  - FREQNESS_NetworkEstimation(...) 
%    Performs Generalized Eigendecomposition (GED) over a user-defined 
%    sample of frequencies, separating frequency-resolved networks.
%
%  - FREQNESS_Visualizer(...) 
%    Takes as input selected outputs from `FREQNESS_NetworkEstimation` 
%    to visualize both the network landscape and brain topographies.
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

