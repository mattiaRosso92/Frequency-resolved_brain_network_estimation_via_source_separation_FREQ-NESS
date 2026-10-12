"""Complete Python counterpart of ../FREQNESS_MainPipeline.m.

Block 0 defines the settings. Block 1 loads, estimates, saves and releases
one participant at a time. Block 2 combines saved results for optional group
analyses. The CORE block runs by default; group analyses are disabled.

Keep RUN_CORE = True and RUN_GROUP_ANALYSES = False for the initial CORE run
to limit memory use. Once CORE finishes, set RUN_CORE = False and
RUN_GROUP_ANALYSES = True, keep the same estimation settings, and run again.
Use ANALYSES_TO_RUN to select group sections. Command-line --analysis
selections also enable the requested group sections.

If you use this toolbox, please cite:
Rosso, M., Fernández-Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
Kringelbach, M. L., & Bonetti, L. (2025). FREQ-NESS Reveals the Dynamic
Reconfiguration of Frequency-Resolved Brain Networks During Auditory
Stimulation. Advanced Science, 2413195.
https://doi.org/10.1002/advs.202413195

Toolbox authors
---------------
Mattia Rosso - Center for Music in the Brain, Aarhus University -
mattia.rosso@clin.au.dk
Chiara Malvaso - University of Bologna - chiara.malvaso2@unibo.it
Leonardo Bonetti - Center for Music in the Brain, Aarhus University; Centre
for Eudaimonia and Human Flourishing, Linacre College, University of Oxford -
leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
"""

from pathlib import Path
import sys

import numpy as np

from freqness import (
    ALL_ANALYSES,
    FREQNESS_MainPipeline,
    FREQNESSPipelineConfig,
)
from freqness.FREQNESS_MainPipeline import main


# ============================================================================
# BLOCK 0 - STARTUP AND ANALYSIS SETTINGS
# ============================================================================

TOOLBOX_ROOT = Path(__file__).resolve().parent.parent

# Core results are saved outside FREQNESS_Data, in one folder per condition.
RUN_CORE = True
RUN_GROUP_ANALYSES = False
PATH_RESULTS = None  # None uses FREQNESS_Results inside TOOLBOX_ROOT
OVERWRITE_CORE = False

# Optional sections for Block 2. Group arrays remain in memory when enabled.
# ANALYSES_TO_RUN = ("FREQNESS_EntropyLandscape",)
ANALYSES_TO_RUN = ALL_ANALYSES


# ============================================================================
# 1) FREQNESS_NetworkEstimation
# ============================================================================

FREQUENCIES = np.arange(1.2, 24.0 + 0.6, 1.2)
SAMPLING_RATE = 250.0

# Optional arguments can be added here without changing the pipeline engine:
NETWORK_OPTIONS = {
    # "duration": 20,
    # "fwidth": None,  # scalar or one explicit width per frequency
    # "filter": "logarithmic",
    # "regularisation": 0.01,
    "ncomps": 5,  # components retained in evecs, pats and ts; evals retains all
    # "bad_segments": None,
    # "rescale": False,
}


# ============================================================================
# 2) FREQNESS_Visualizer
# ============================================================================

PLOT_ALL = False
LANDSCAPE_FREQUENCIES = FREQUENCIES
LANDSCAPE_COMPONENTS = 10
PATTERN_FREQUENCIES = np.array([2.4, 8.4])
PATTERN_COMPONENTS = 1
SAVE_NIFTI = True


# ============================================================================
# 3) FREQNESS_EntropyLandscape
# ============================================================================

# This section requires no settings beyond the saved FREQ structure.


# ============================================================================
# 4) FREQNESS_ExponentialDK
# ============================================================================

EXPDK_RANGE_TO_FIT = np.array([9.0, 20.0])
EXPDK_COMPONENT = 1


# ============================================================================
# 5) FREQNESS_FreqGradients
# ============================================================================

FREQGRAD_RANGE_TO_MODEL = np.array([8.0, 12.0])
FREQGRAD_COMPONENT = 1


# ============================================================================
# 6) FREQNESS_CompGradients
# ============================================================================

COMPGRAD_FREQUENCY = 10.0
COMPGRAD_COMPONENTS = np.array([1, 5], dtype=np.int64)


# ============================================================================
# 7) FREQNESS_CrossCoupling
# ============================================================================

LFO_FREQUENCY = 2.0


# ============================================================================
# 8) FREQNESS_BackProjection
# ============================================================================

BACKPROJ_FREQUENCY = 8.4
BACKPROJ_COMPONENTS = np.array([1], dtype=np.int64)


# ============================================================================
# 9) FREQNESS_NetworkRemoval
# ============================================================================

NETREM_FREQUENCY = 8.4
NETREM_COMPONENTS = np.array([1], dtype=np.int64)
NETREM_PLOT_LANDSCAPE = True
NETREM_LANDSCAPE_COMPONENTS = 3


# ============================================================================
# OUTPUT AND DISPLAY SETTINGS
# ============================================================================

SHOW_FIGURES = True
VISUALIZER_OUTPUT_DIRECTORY = None


CONFIG = FREQNESSPipelineConfig(
    toolbox_root=TOOLBOX_ROOT,
    run_core=RUN_CORE,
    run_group_analyses=RUN_GROUP_ANALYSES,
    path_results=PATH_RESULTS,
    overwrite_core=OVERWRITE_CORE,
    frex=FREQUENCIES,
    srate=SAMPLING_RATE,
    analyses=ANALYSES_TO_RUN,
    network_options=NETWORK_OPTIONS,
    plot_all=PLOT_ALL,
    show=SHOW_FIGURES,
    save_nifti=SAVE_NIFTI,
    visualizer_output_directory=VISUALIZER_OUTPUT_DIRECTORY,
    landscape_frex=LANDSCAPE_FREQUENCIES,
    landscape_ncomps=LANDSCAPE_COMPONENTS,
    pattern_frex=PATTERN_FREQUENCIES,
    pattern_ncomps=PATTERN_COMPONENTS,
    expdk_range2fit=EXPDK_RANGE_TO_FIT,
    expdk_which_comp=EXPDK_COMPONENT,
    freqgrad_frex2model=FREQGRAD_RANGE_TO_MODEL,
    freqgrad_comp2model=FREQGRAD_COMPONENT,
    compgrad_freq2model=COMPGRAD_FREQUENCY,
    compgrad_comps2model=COMPGRAD_COMPONENTS,
    lfo_freq=LFO_FREQUENCY,
    backproj_freq2project=BACKPROJ_FREQUENCY,
    backproj_comps2project=BACKPROJ_COMPONENTS,
    netrem_freq2remove=NETREM_FREQUENCY,
    netrem_comps2remove=NETREM_COMPONENTS,
    netrem_plot_landscape=NETREM_PLOT_LANDSCAPE,
    netrem_landscape_ncomps=NETREM_LANDSCAPE_COMPONENTS,
)


# ============================================================================
# BLOCK 1 - NETWORK ESTIMATION AND CORE EXPORT
# BLOCK 2 - OPTIONAL GROUP ANALYSES
# ============================================================================
# The package pipeline runs these blocks in order using the settings above.
# Each complete FREQ is saved to CORE/<condition>/<participant>_FREQ.npz.
# Existing files are reused only when input metadata and settings match.


def run_configured_pipeline():
    """Run the configured CORE and optional group-analysis blocks."""
    result = FREQNESS_MainPipeline(CONFIG)
    completed = sum(bool(condition.core_files) for condition in result.conditions)
    print(
        f"FREQ-NESS pipeline completed for {completed} non-empty "
        f"condition(s)."
    )
    return result


if __name__ == "__main__":
    if len(sys.argv) == 1:
        run_configured_pipeline()
    else:
        raise SystemExit(main(default_toolbox_root=TOOLBOX_ROOT, default_config=CONFIG))
