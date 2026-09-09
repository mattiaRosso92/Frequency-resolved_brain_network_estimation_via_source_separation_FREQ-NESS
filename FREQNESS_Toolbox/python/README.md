# FREQ-NESS for Python

This directory contains the Python implementation of the FREQ-NESS toolbox.
It follows the MATLAB implementation in the parent directory.

The first public API preserves the established FREQ-NESS function names:

```python
from freqness import (
    FREQNESS_BackProjection,
    FREQNESS_CompGradients,
    FREQNESS_ComputeFilterWidths,
    FREQNESS_CrossCoupling,
    FREQNESS_EntropyLandscape,
    FREQNESS_ExponentialDK,
    FREQNESS_FreqGradients,
    FREQNESS_InducedResponses,
    FREQNESS_MainPipeline,
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
freqGradCoeff, freqGoodFit = FREQNESS_FreqGradients(
    FREQ,
    MNI,
    frex2model=[4, 30],
    comp2model=1,
    show=False,
)
compGradCoeff, compGoodFit = FREQNESS_CompGradients(
    FREQ,
    MNI,
    freq2model=10,
    comps2model=[1, 3],
    show=False,
)

# Run the complete folder-based workflow, or select exact function sections.
pipeline = FREQNESS_MainPipeline()
```

When `FREQNESS_NetworkEstimation` receives no explicit `fwidth`, its default
`logarithmic` schedule varies spectral selectivity smoothly from Q = 7 at the
lowest requested frequency to Q = 3.5 at the highest. The `linear` option uses
constant Q = 5 (`FWHM = frequency / 5`). These automatic widths depend on
frequencies in Hz and their numerical range, not on the number or spacing of
frequency bins. A scalar `fwidth` still applies one constant width to every
filter, while a vector supplies one explicit width per frequency. The same
schedule is independently available through `FREQNESS_ComputeFilterWidths`.

The spatial-pattern view plots all valid MNI locations as small black points.
Requested frequencies follow a low-to-high `viridis` colour gradient, while
normalized activation magnitude controls marker size. Pass
`frequency_panels=True` for separate frequency-specific views, or
`save_nifti=False` when only figures are required.

`FREQNESS_EntropyLandscape` computes and returns both quadratic Rényi entropy
(`H2`) and effective dimensionality (`ED`) to preserve MATLAB numerical
compatibility. `FREQ.evals` always contains the complete normalized
eigenspectrum, independently of the number of eigenvectors, patterns, and time
series retained through `ncomps`. Entropy therefore uses every eigenvalue while
component-dependent outputs remain memory-limited. Its figures intentionally
show only `H2`.

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

`FREQNESS_FreqGradients` and `FREQNESS_CompGradients` model whether dominant
spatial-pattern locations vary linearly or quadratically along the MNI X, Y,
and Z axes. Both return coefficients in `[b0, b1, b2]` order and use BIC to
select model order. The regression targets remain one-based frequency-bin or
component indices for MATLAB compatibility. The legacy thresholds are exposed
as `threshold_sd` (defaults: 1 for frequency gradients and 2 for component
gradients), and retained-voxel counts are included in the fit diagnostics.

`FREQNESS_MainPipeline` mirrors the complete MATLAB pipeline across every
folder in `FREQNESS_Data`. By default it runs the complete workflow and keeps
numerical results in memory. Users can select one or more sections by exact
function name when only a subset of analyses is needed. Network estimation is
automatically computed as the prerequisite for selected secondary analyses.
The BackProjection and NetworkRemoval sections use MATLAB-compatible,
one-based component numbering. NetworkRemoval returns the cleaned data and
removed activity in memory; when requested, the pipeline re-estimates the
cleaned data and plots its first three network components. The re-estimation
is deliberately kept in the pipeline, leaving NetworkRemoval as a focused
subtraction function. The pipeline keeps numerical outputs in memory; only
the Visualizer's established NIfTI export can write files, when enabled.

Long-running Python functions mirror the MATLAB console messages, including
condition, participant, frequency, modelling, backprojection, removal, and
NIfTI-export progress.

For interactive use, open `FREQNESS_NotebookPipeline.ipynb` in Jupyter or VS
Code. The notebook guides one selected condition through every main-pipeline
analysis and keeps intermediate outputs visible. `FREQNESS_MainPipeline.py`
remains the recommended interface for unattended processing of every dataset
folder, while `src/freqness/FREQNESS_MainPipeline.py` provides the reusable
package engine and configuration API.

From the `python` directory, run only the core section with:

```text
python FREQNESS_MainPipeline.py \
    --analysis FREQNESS_NetworkEstimation \
    --no-show
```

Run only one secondary section, for example:

```text
python FREQNESS_MainPipeline.py \
    --analysis FREQNESS_EntropyLandscape \
    --no-show
```

Install visualization dependencies with:

```text
python -m pip install -e ".[visualize]"
```

Python-style aliases are also available:

```python
from freqness import (
    back_project,
    component_gradients,
    cross_coupling,
    entropy_landscape,
    estimate_networks,
    exponential_decay,
    frequency_gradients,
    induced_responses,
    main_pipeline,
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

If you use this toolbox, please cite:

Rosso, M., Fernández-Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
Kringelbach, M. L., & Bonetti, L. (2025). FREQ-NESS Reveals the Dynamic
Reconfiguration of Frequency-Resolved Brain Networks During Auditory
Stimulation. *Advanced Science*, 2413195.
https://doi.org/10.1002/advs.202413195

## Toolbox authors

- Mattia Rosso — Center for Music in the Brain, Aarhus University —
  mattia.rosso@clin.au.dk
- Chiara Malvaso — University of Bologna — chiara.malvaso2@unibo.it
- Leonardo Bonetti — Center for Music in the Brain, Aarhus University; Centre
  for Eudaimonia and Human Flourishing, Linacre College, University of Oxford —
  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
