# FREQ-NESS

## Download the toolbox

**[Download the latest FREQ-NESS Toolbox release](https://github.com/mattiaRosso92/Frequency-resolved_brain_network_estimation_via_source_separation_FREQ-NESS/releases/latest/download/FREQNESS.zip)**

The recommended release archive extracts to the short folder name `FREQNESS`,
which avoids Windows path-length problems caused by the repository's long
name. See the [toolbox instructions](FREQNESS_Toolbox/README.md) for MATLAB and
Python setup.

## Citation

If you use this toolbox, please cite:

Rosso, M., Fernández‐Rubio, G., Keller, P. E., Brattico, E., Vuust, P., Kringelbach, M. L., & Bonetti, L. (2025). 
FREQ‐NESS Reveals the Dynamic Reconfiguration of Frequency‐Resolved Brain Networks During Auditory Stimulation. 
Advanced Science, 2413195.
https://doi.org/10.1002/advs.202413195

% ========================================================================

This repository contains two folders:

1) FREQNESS_AdvancedScience_2025: Code and functions used for the FREQNESS paper published in Advanced Science (Rosso et al., 2025)

2) FREQNESS_Toolbox: MATLAB and Python implementations, main pipelines, documentation, examples, spatial coordinates, and supporting functions

% ========================================================================

Abstract: The brain is a dynamic system whose network organization is often studied by focusing on specific frequency bands or anatomical regions, leading to fragmented insights, or by employing complex and elaborate methods that hinder straightforward interpretations. To address this issue, a new analytical pipeline named FREQuency-resolved Network Estimation via Source Separation (FREQ-NESS) is introduced. This is designed to estimate the activation and spatial configuration of simultaneous brain networks across frequencies by analyzing the frequency-resolved multivariate covariance between whole-brain voxel time series. FREQ-NESS is applied to source-reconstructed magnetoencephalography (MEG) data during resting state and isochronous auditory stimulation. Results reveal simultaneous, frequency-specific brain networks during resting state, such as the default mode, alpha-band, and motor-beta networks. During auditory stimulation, FREQ-NESS detects: (1) emergence of networks attuned to the stimulation frequency, (2) spatial reorganization of existing networks, such as alpha-band networks shifting from occipital to sensorimotor areas, (3) stability of networks unaffected by auditory stimuli. Furthermore, auditory stimulation significantly enhances cross-frequency coupling, with the phase of attuned auditory networks modulating the gamma band amplitude of medial temporal lobe networks. In conclusion, FREQ-NESS effectively maps the brain’s spatiotemporal dynamics, providing a comprehensive view of brain function by revealing simultaneous, frequency-resolved networks and their interaction.

% ========================================================================

The conceptualization and code implementation of this work were carried out in collaboration with Prof. Leonardo Bonetti (https://github.com/leonardob92).

Corresponding authors:
Mattia Rosso     - mattia.rosso@clin.au.dk
Leonardo Bonetti - leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk

% ========================================================================



