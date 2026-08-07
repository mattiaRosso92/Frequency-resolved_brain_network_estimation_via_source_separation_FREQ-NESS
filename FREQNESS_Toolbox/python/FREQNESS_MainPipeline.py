"""Complete Python counterpart of ../FREQNESS_MainPipeline.m.

The full pipeline is configured below. During MATLAB/Python validation, change
`ANALYSES_TO_RUN` to a one-item tuple so that one complete function section is
executed at a time. Command-line `--analysis` selections provide the same
control without editing this file.
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
# STARTUP AND SECTION SELECTION
# ============================================================================

TOOLBOX_ROOT = Path(__file__).resolve().parent.parent

# Complete pipeline by default. For staged testing, use for example:
# ANALYSES_TO_RUN = ("FREQNESS_NetworkEstimation",)
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
    # "fwidth": 0.1,
    # "filter": "logarithmic",
    # "regularisation": 0.01,
    # "ncomps": 30,
    # "bad_segments": None,
}


# ============================================================================
# 2) FREQNESS_Visualizer
# ============================================================================

PLOT_ALL = False
LANDSCAPE_FREQUENCIES = FREQUENCIES
LANDSCAPE_COMPONENTS = 10
PATTERN_FREQUENCIES = np.array([2.0, 8.0])
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
# OUTPUT AND DISPLAY SETTINGS
# ============================================================================

SHOW_FIGURES = True
SAVE_COMPARISON_OUTPUTS = True
OUTPUT_DIRECTORY = None


CONFIG = FREQNESSPipelineConfig(
    toolbox_root=TOOLBOX_ROOT,
    frex=FREQUENCIES,
    srate=SAMPLING_RATE,
    analyses=ANALYSES_TO_RUN,
    network_options=NETWORK_OPTIONS,
    plot_all=PLOT_ALL,
    show=SHOW_FIGURES,
    save_outputs=SAVE_COMPARISON_OUTPUTS,
    save_nifti=SAVE_NIFTI,
    output_directory=OUTPUT_DIRECTORY,
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
)


def run_configured_pipeline():
    """Run every section listed in `ANALYSES_TO_RUN`."""
    result = FREQNESS_MainPipeline(CONFIG)
    completed = sum(bool(condition.outputs) for condition in result.conditions)
    print(
        f"FREQ-NESS pipeline completed for {completed} non-empty "
        f"condition(s)."
    )
    for condition in result.conditions:
        for analysis, output_path in condition.files.items():
            print(f"{condition.name} | {analysis} | {output_path}")
    return result


if __name__ == "__main__":
    if len(sys.argv) == 1:
        run_configured_pipeline()
    else:
        raise SystemExit(main(default_toolbox_root=TOOLBOX_ROOT))
