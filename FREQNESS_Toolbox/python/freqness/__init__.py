"""Python implementation of the FREQ-NESS toolbox.

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

from .FREQNESS_BackProjection import FREQNESS_BackProjection
from .FREQNESS_CoreIO import (
    FREQNESS_LoadCore,
    FREQNESS_LoadCoreFields,
    FREQNESS_SaveCore,
)
from .FREQNESS_CrossCoupling import (
    FREQNESS_CrossCoupling,
    FREQNESSCrossCouplingResult,
)
from .FREQNESS_CompGradients import FREQNESS_CompGradients
from .FREQNESS_ComputeFilterWidths import FREQNESS_ComputeFilterWidths
from .FREQNESS_EntropyLandscape import (
    FREQNESS_EntropyLandscape,
    FREQNESSEntropyLandscapeResult,
)
from .FREQNESS_ExponentialDK import (
    FREQNESS_ExponentialDK,
    FREQNESSExponentialGoodFit,
)
from .FREQNESS_FreqGradients import (
    FREQNESS_FreqGradients,
    FREQNESSGradientGoodFit,
)
from .FREQNESS_InducedResponses import (
    FREQNESS_InducedResponses,
    FREQNESSInducedResponsesResult,
)
from .FREQNESS_IndividualPeak import (
    FREQNESS_IndividualPeak,
    FREQNESSIndividualPeakResult,
)
from .FREQNESS_LandscapeComparison import (
    FREQNESS_LandscapeComparison,
    FREQNESSLandscapeComparisonResult,
    FREQNESSLandscapeFactorialComparisonResult,
    FREQNESSLandscapeFactorialEffect,
    FREQNESSLandscapeFactorialMeasureResult,
    FREQNESSLandscapeMeasureResult,
)
from .FREQNESS_MainPipeline import (
    ALL_ANALYSES,
    BACK_PROJECTION,
    COMP_GRADIENTS,
    CROSS_COUPLING,
    ENTROPY_LANDSCAPE,
    EXPONENTIAL_DK,
    FREQNESS_MainPipeline,
    FREQNESSPipelineConditionResult,
    FREQNESSPipelineConfig,
    FREQNESSPipelineResult,
    FREQ_GRADIENTS,
    LANDSCAPE_COMPARISON,
    NETWORK_ESTIMATION,
    NETWORK_REMOVAL,
    SPATIAL_COMPARISON,
    STABILITY_INDEX,
    VISUALIZER,
)
from .FREQNESS_NetworkEstimation import (
    FREQNESS_NetworkEstimation,
    FREQNESSResult,
)
from .FREQNESS_NetworkRemoval import FREQNESS_NetworkRemoval
from .FREQNESS_Startup import FREQNESS_Startup
from .FREQNESS_SpatialComparison import (
    FREQNESS_SpatialComparison,
    FREQNESSSpatialCluster,
    FREQNESSSpatialComparisonResult,
    FREQNESSSpatialFactorialCluster,
    FREQNESSSpatialFactorialComparisonResult,
    FREQNESSSpatialFactorialEffect,
)
from .FREQNESS_StabilityIndex import (
    FREQNESS_StabilityIndex,
    FREQNESSStabilityIndexResult,
)
from .FREQNESS_Visualizer import FREQNESS_Visualizer, FREQNESSVisualization
from .filterFGx import filterFGx

estimate_networks = FREQNESS_NetworkEstimation
back_project = FREQNESS_BackProjection
cross_coupling = FREQNESS_CrossCoupling
component_gradients = FREQNESS_CompGradients
entropy_landscape = FREQNESS_EntropyLandscape
exponential_decay = FREQNESS_ExponentialDK
frequency_gradients = FREQNESS_FreqGradients
induced_responses = FREQNESS_InducedResponses
individual_peak = FREQNESS_IndividualPeak
landscape_comparison = FREQNESS_LandscapeComparison
main_pipeline = FREQNESS_MainPipeline
remove_network = FREQNESS_NetworkRemoval
startup = FREQNESS_Startup
spatial_comparison = FREQNESS_SpatialComparison
stability_index = FREQNESS_StabilityIndex
visualize = FREQNESS_Visualizer

__all__ = [
    "FREQNESS_BackProjection",
    "FREQNESS_CrossCoupling",
    "FREQNESSCrossCouplingResult",
    "FREQNESS_CompGradients",
    "FREQNESS_ComputeFilterWidths",
    "FREQNESS_EntropyLandscape",
    "FREQNESSEntropyLandscapeResult",
    "FREQNESS_ExponentialDK",
    "FREQNESSExponentialGoodFit",
    "FREQNESS_FreqGradients",
    "FREQNESSGradientGoodFit",
    "FREQNESS_InducedResponses",
    "FREQNESSInducedResponsesResult",
    "FREQNESS_IndividualPeak",
    "FREQNESSIndividualPeakResult",
    "FREQNESS_LandscapeComparison",
    "FREQNESSLandscapeComparisonResult",
    "FREQNESSLandscapeFactorialComparisonResult",
    "FREQNESSLandscapeFactorialEffect",
    "FREQNESSLandscapeFactorialMeasureResult",
    "FREQNESSLandscapeMeasureResult",
    "FREQNESS_LoadCore",
    "FREQNESS_LoadCoreFields",
    "FREQNESS_SaveCore",
    "FREQNESS_MainPipeline",
    "FREQNESSPipelineConditionResult",
    "FREQNESSPipelineConfig",
    "FREQNESSPipelineResult",
    "FREQNESS_NetworkEstimation",
    "FREQNESS_NetworkRemoval",
    "FREQNESS_Startup",
    "FREQNESS_SpatialComparison",
    "FREQNESSSpatialCluster",
    "FREQNESSSpatialComparisonResult",
    "FREQNESSSpatialFactorialCluster",
    "FREQNESSSpatialFactorialComparisonResult",
    "FREQNESSSpatialFactorialEffect",
    "FREQNESS_StabilityIndex",
    "FREQNESSStabilityIndexResult",
    "FREQNESSResult",
    "FREQNESS_Visualizer",
    "FREQNESSVisualization",
    "back_project",
    "cross_coupling",
    "component_gradients",
    "entropy_landscape",
    "estimate_networks",
    "exponential_decay",
    "frequency_gradients",
    "filterFGx",
    "induced_responses",
    "individual_peak",
    "landscape_comparison",
    "main_pipeline",
    "remove_network",
    "startup",
    "spatial_comparison",
    "stability_index",
    "visualize",
    "ALL_ANALYSES",
    "BACK_PROJECTION",
    "COMP_GRADIENTS",
    "CROSS_COUPLING",
    "ENTROPY_LANDSCAPE",
    "EXPONENTIAL_DK",
    "FREQ_GRADIENTS",
    "LANDSCAPE_COMPARISON",
    "NETWORK_ESTIMATION",
    "NETWORK_REMOVAL",
    "SPATIAL_COMPARISON",
    "STABILITY_INDEX",
    "VISUALIZER",
]

__version__ = "3.2.0"
