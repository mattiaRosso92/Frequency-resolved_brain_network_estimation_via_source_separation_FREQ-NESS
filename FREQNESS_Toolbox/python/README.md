# FREQ-NESS for Python

This directory contains the Python implementation of the FREQ-NESS toolbox.
The existing MATLAB implementation remains unchanged in the parent directory.

The first public API preserves the established FREQ-NESS function names:

```python
from freqness import (
    FREQNESS_EntropyLandscape,
    FREQNESS_ExponentialDK,
    FREQNESS_NetworkEstimation,
    FREQNESS_Startup,
    FREQNESS_Visualizer,
)

allData, MNI, path_home = FREQNESS_Startup("/path/to/FREQNESS_Toolbox")
FREQ = FREQNESS_NetworkEstimation(allData[0], frex, srate)

Landscape = {"frex": FREQ.frex, "ncomps": 3}
Patterns = {
    "MNI_coords": MNI,
    "frex": [2.4, 8, 12],
    "ncomps": 3,
    "path_output": "/path/to/results",
}
visualization = FREQNESS_Visualizer(
    FREQ,
    Landscape,
    Patterns,
    plot_all=False,
)

H2, ED = FREQNESS_EntropyLandscape(FREQ)
decayCoeff, goodFit = FREQNESS_ExponentialDK(
    FREQ,
    which_comp=1,
    range2fit=[2, 30],
)
```

The spatial-pattern view plots all valid MNI locations as small black points.
Requested frequencies follow a low-to-high `viridis` colour gradient, while
normalized activation magnitude controls marker size. Pass
`frequency_panels=True` for separate frequency-specific views, or
`save_nifti=False` when only figures are required.

`FREQNESS_EntropyLandscape` computes and returns both quadratic Rényi entropy
(`H2`) and effective dimensionality (`ED`) to preserve MATLAB numerical
compatibility. Its figures intentionally show only `H2`.

`FREQNESS_ExponentialDK` fits `A * exp(-lambda * frequency)` to one component
using MATLAB-compatible one-based component numbering. It returns one decay
coefficient and original-space R² value per participant.

Install visualization dependencies with:

```text
python -m pip install -e ".[visualize]"
```

Python-style aliases are also available:

```python
from freqness import (
    entropy_landscape,
    estimate_networks,
    exponential_decay,
    startup,
    visualize,
)
```

## Development installation

From this directory:

```text
python -m pip install -e ".[test]"
pytest
```

## Citation

Rosso, M., Fernández-Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
Kringelbach, M. L., & Bonetti, L. (2025). FREQ-NESS Reveals the Dynamic
Reconfiguration of Frequency-Resolved Brain Networks During Auditory
Stimulation. *Advanced Science*, 2413195.
https://doi.org/10.1002/advs.202413195
