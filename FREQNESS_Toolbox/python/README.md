# FREQ-NESS for Python

This directory contains the Python implementation of the FREQ-NESS toolbox.
The existing MATLAB implementation remains unchanged in the parent directory.

The first public API preserves the established FREQ-NESS function names:

```python
from freqness import (
    FREQNESS_BackProjection,
    FREQNESS_CrossCoupling,
    FREQNESS_EntropyLandscape,
    FREQNESS_ExponentialDK,
    FREQNESS_InducedResponses,
    FREQNESS_NetworkEstimation,
    FREQNESS_NetworkRemoval,
    FREQNESS_Startup,
    FREQNESS_Visualizer,
)

allData, MNI, path_home = FREQNESS_Startup("/path/to/FREQNESS_Toolbox")
FREQ = FREQNESS_NetworkEstimation(allData[0], frex, srate)
backProj = FREQNESS_BackProjection(
    FREQ,
    freq2project=10,
    comps2project=[1, 2],
)
dataClean, removedActivity = FREQNESS_NetworkRemoval(
    FREQ,
    allData[0],
    freq2remove=10,
    comps2remove=[1, 2],
)

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
CFC = FREQNESS_CrossCoupling(
    FREQ,
    lfo_freq=2.4,
    frex2model=[4, 30],
    which_comp=1,
)
IND = FREQNESS_InducedResponses(
    FREQ,
    events=[501, 1501, 2501],
    epoch_window=[-0.5, 1.0],
    baseline_window=[-0.4, 0.0],
    which_comp=1,
    frex2model=[2, 30],
    plot_avg=True,
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

`FREQNESS_CrossCoupling` computes coupling between the phase of a neural
low-frequency FREQ-NESS component and the power of the same component at
higher carrier frequencies. Phase-binned carrier power is fitted with a fixed
first-harmonic regression. This replaces the external MATLAB `sineFit.m`
dependency with deterministic linear least squares and adds raw and normalized
amplitude, preferred phase, MSE, R², fitted PAC curves, and valid-bin counts.

`FREQNESS_BackProjection` reconstructs the broadband voxel-space contribution
of one or more networks from the complete retained GED filter set and stored
component time series. Component numbers remain one-based for MATLAB
compatibility. The closest analyzed frequency is used when necessary, and
stored rescaling factors are removed so output is expressed in the original
data units.

`FREQNESS_NetworkRemoval` uses that backprojection to subtract selected
broadband network contributions from the original analyzed data. It returns
both the cleaned data and removed activity in matching dimensions and original
input units, with one-based component numbering retained for MATLAB
compatibility.

`FREQNESS_InducedResponses` continuously filters selected broadband component
time series at their corresponding network frequencies, computes Hilbert
power, and applies trial-wise baseline normalization in decibels before trial
averaging. Event samples and component numbers remain one-based. Epoch
endpoints are inclusive, while baseline windows use `[start, end)` so a
prestimulus baseline ending at zero excludes event onset.

Install visualization dependencies with:

```text
python -m pip install -e ".[visualize]"
```

Python-style aliases are also available:

```python
from freqness import (
    back_project,
    cross_coupling,
    entropy_landscape,
    estimate_networks,
    exponential_decay,
    induced_responses,
    remove_network,
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
