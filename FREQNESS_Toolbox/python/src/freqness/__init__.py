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
from .FREQNESS_CrossCoupling import (
    FREQNESS_CrossCoupling,
    FREQNESSCrossCouplingResult,
)
from .FREQNESS_CompGradients import FREQNESS_CompGradients
from .FREQNESS_EntropyLandscape import FREQNESS_EntropyLandscape
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
    NETWORK_ESTIMATION,
    NETWORK_REMOVAL,
    VISUALIZER,
)
from .FREQNESS_NetworkEstimation import (
    FREQNESS_NetworkEstimation,
    FREQNESSResult,
)
from .FREQNESS_NetworkRemoval import FREQNESS_NetworkRemoval
from .FREQNESS_Startup import FREQNESS_Startup
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
main_pipeline = FREQNESS_MainPipeline
remove_network = FREQNESS_NetworkRemoval
startup = FREQNESS_Startup
visualize = FREQNESS_Visualizer

__all__ = [
    "FREQNESS_BackProjection",
    "FREQNESS_CrossCoupling",
    "FREQNESSCrossCouplingResult",
    "FREQNESS_CompGradients",
    "FREQNESS_EntropyLandscape",
    "FREQNESS_ExponentialDK",
    "FREQNESSExponentialGoodFit",
    "FREQNESS_FreqGradients",
    "FREQNESSGradientGoodFit",
    "FREQNESS_InducedResponses",
    "FREQNESSInducedResponsesResult",
    "FREQNESS_MainPipeline",
    "FREQNESSPipelineConditionResult",
    "FREQNESSPipelineConfig",
    "FREQNESSPipelineResult",
    "FREQNESS_NetworkEstimation",
    "FREQNESS_NetworkRemoval",
    "FREQNESS_Startup",
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
    "main_pipeline",
    "remove_network",
    "startup",
    "visualize",
    "ALL_ANALYSES",
    "BACK_PROJECTION",
    "COMP_GRADIENTS",
    "CROSS_COUPLING",
    "ENTROPY_LANDSCAPE",
    "EXPONENTIAL_DK",
    "FREQ_GRADIENTS",
    "NETWORK_ESTIMATION",
    "NETWORK_REMOVAL",
    "VISUALIZER",
]

__version__ = "0.1.0.dev0"
