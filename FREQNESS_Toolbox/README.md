# FREQ-NESS Toolbox

FREQ-NESS extracts frequency-resolved brain networks from EEG, MEG, or iEEG
data using generalized eigendecomposition. MATLAB and Python implementations
are included and follow the same toolbox organization and function names.

## Recommended download

**[Download the latest FREQ-NESS Toolbox release](https://github.com/mattiaRosso92/Frequency-resolved_brain_network_estimation_via_source_separation_FREQ-NESS/releases/latest/download/FREQNESS.zip)**

Use this release asset rather than GitHub's automatically generated **Source
code** archives. The recommended archive extracts to the short top-level folder
`FREQNESS`, avoiding the Windows path-length problem caused by the repository's
long name.

On Windows, extract it near the drive root, for example:

```text
C:\FREQNESS
```

The extracted top-level folder may safely be renamed or moved. Internal paths
are resolved relative to the toolbox rather than from the repository's name.

## Data organization

Put one condition or group in each folder inside `FREQNESS_Data`. Put each
participant in a separate `.mat` file containing exactly one numeric
variables-by-time matrix. Participants within a condition must have matching
matrix dimensions.

```text
FREQNESS/
├── FREQNESS_Data/
│   ├── Condition_1/
│   │   ├── participant_01.mat
│   │   └── participant_02.mat
│   └── Condition_2/
│       └── participant_01.mat
├── FREQNESS_MainPipeline.m
└── python/
    ├── FREQNESS_MainPipeline.py
    ├── FREQNESS_NotebookPipeline.ipynb
    └── freqness/
        └── FREQNESS Python package
```

The bundled `FREQNESS_MNI_Coordinates` file corresponds to the 3559-voxel
MNI152 8-mm grid. Replace it when using a different source space or voxel
ordering. Matching coordinates are required for spatial visualizations and
gradient modelling. `FREQNESS_SpatialComparison` uses the same coordinates to
perform pairwise or multidimensional factorial voxel-wise comparisons with
MNI-adjacency cluster correction. Generalized designs are supplied as a cell
array whose dimensions represent experimental factors; each factor can be
declared as within- or between-subject, while established two-structure calls
remain supported.
`FREQNESS_StabilityIndex` quantifies the temporal frequency stability of
selected network time series without requiring behavioral data.
`FREQNESS_LandscapeComparison` compares participant-level landscapes using
pairwise or multidimensional factorial designs. Generalized MATLAB calls use
the same cell-array factor representation as `FREQNESS_SpatialComparison`
and support within-subject, between-subject, and mixed designs. Internal
Benjamini-Hochberg FDR correction is applied across tested frequencies for
each effect, measure, and component by default.
For both comparison functions, an incompatible non-parametric test selection
is replaced by the design-appropriate test with a warning: `signrank` for
within-subject data and `ranksum` for independent groups.

For group analyses, every participant must use the same standardized source
space: the source count, ordering, and anatomical correspondence must all
match. Equal matrix dimensions alone cannot verify anatomical correspondence.
Source-space transformation or morphing must therefore be completed upstream
with the software used for source reconstruction. FREQ-NESS does not perform
that transformation. Python warns when an input contains 8,000 or more sources,
because its dense covariance matrices and generalized eigendecomposition can
require substantial memory and computation time; ordinary lower-resolution
inputs do not trigger this warning.

## MATLAB quick start

1. Open `FREQNESS_MainPipeline.m` in MATLAB.
2. Edit the settings in `DEFINE ANALYSIS SETTINGS`, including sampling rate,
   frequencies, and the analyses to run.
3. Run the script. It initializes the toolbox automatically.

The MATLAB pipeline is organized into four blocks:

- **Block 0 — Startup and settings:** inspects participant filenames and matrix
  dimensions, loads MNI coordinates, and defines the analysis settings.
- **Block 1 — Network estimation and core export:** loads one participant at a
  time, estimates its networks, saves the complete `FREQ` structure, and releases
  the participant data and result before proceeding to the next file.
- **Block 2 — Optional group analyses:** combines the saved participant results
  and runs the established visualization and secondary-analysis sections,
  including `FREQNESS_StabilityIndex` when `run_stability_index = true`.
- **Block 3 — Selective comparisons:** matches participant identifiers and
  builds compact statistical inputs containing only the requested landscape,
  its frequency axis, and the participant dimension. It never reconstructs
  group-level `ts` or `evecs` arrays.

By default, `run_core = true` and `run_group_analyses = false`. Core results are
saved under `FREQNESS_Results/CORE/<condition>/<participant>_FREQ.mat`. Each file
contains `FREQ`, `PARTICIPANT` source-file metadata, and `NETWORK_SETTINGS`.
Small top-level `FREQ_evals`, `FREQ_pats`, and `FREQ_frex` caches support
selective Block 3 loading; the full component time series remain only inside
`FREQ` for subsequent analyses. Files use MATLAB Version 7.3 to support large
arrays.

The MATLAB MainPipeline retains five components per frequency by default
(`network_ncomps = 5`). This setting limits `FREQ.evecs`, `FREQ.pats`, and
`FREQ.ts`; `FREQ.evals` always contains the complete eigenspectrum. The default
of the standalone `FREQNESS_NetworkEstimation` function is unchanged.

Set `overwrite_core = true` to replace saved results. Otherwise, existing files
are reused only when the saved network settings, source path, dimensions, file
size, and modification time match. These metadata checks do not compare the
contents of the source file. Use another `path_results` to keep separate runs.

Run Block 0 before running either analysis block separately. To run the group
analyses from existing core files, set `run_core = false` and
`run_group_analyses = true`, retaining the network settings used to create the
files. Block 2 restores the participant dimension, keeps group outputs in the
workspace, and requires memory for the combined results. Network removal also
loads the original data for the current condition. NIfTI exports are separated
by condition under `FREQNESS_Results/GROUP`.

Block 3 can run directly from matching CORE files for spatial-pattern and
eigenvalue comparisons. It reads only selected `pats` or `evals` slices and
applies frequency/component selection before concatenating participants.
Entropy (`H2`/`ED`) and cross-coupling amplitude comparisons use the already
compact corresponding Block 2 outputs. Stability Index comparisons similarly
use `STAB.stability_index` or `STAB.mean_frequency` produced in Block 2.
Within-subject comparisons reorder the second condition by participant
identifier; between-group comparisons permit different participant counts.

The original `[allData,MNI,path_home] = FREQNESS_Startup()` call still loads all
datasets. Use `[~,MNI,path_home,folder_data] = FREQNESS_Startup(false)` to inspect
files without loading participant arrays. Matrix contents, including finite
values, are checked when each participant is loaded. The numerical functions
continue to accept both single-participant and multi-participant arrays.

The Python MainPipeline follows the same participant-wise workflow, using
`<participant>_FREQ.npz` files. See the Python quick start below.

Individual functions are available in `FREQNESS_Functions`. See
`README_Toolbox.rtf` and the function headers for detailed MATLAB guidance.

When MATLAB or Python `FREQNESS_NetworkEstimation` receives no explicit
`fwidth`, its default `logarithmic` schedule varies spectral selectivity
smoothly from Q = 7 at the lowest requested frequency to Q = 3.5 at the
highest. The `linear` option uses constant Q = 5 (`FWHM = frequency / 5`).
These automatic widths depend on frequencies in Hz and their numerical range,
not on the number or spacing of frequency bins. A scalar `fwidth` still
applies one constant width to every filter, while a vector supplies one
explicit width per frequency.

## Updating from the previous MainPipeline

Keep your participant files in their existing `FREQNESS_Data/<condition>` folders.
Copy your analysis settings into Block 0 of the updated pipeline, then run it.
The default CORE stage now loads, estimates, saves and releases one participant
at a time. Results are created automatically under
`FREQNESS_Results/CORE/<condition>`; they are no longer collected in the workspace.

The new pipeline default is **five components** for `evecs`, `pats` and `ts`;
`evals` retains all eigenvalues. To reproduce an earlier run, use its original
component count and other settings.

For the previous group analyses, enable `run_group_analyses` in MATLAB or
`RUN_GROUP_ANALYSES` in the Python script. Set `run_core` / `RUN_CORE` to false
to reuse matching CORE files. Group analyses combine results in memory, so
memory use increases with participant count. Existing matching files are skipped;
set `overwrite_core` / `OVERWRITE_CORE` to true to replace them, or choose a
new results directory. Run Block 0 before an individual MATLAB analysis block.

## Python quick start

Python 3.10 or newer is required. Open a terminal in the extracted `python`
folder and create an isolated environment.

Windows:

```text
py -m venv .venv
.venv\Scripts\python -m pip install -e ".[visualize]"
.venv\Scripts\python FREQNESS_MainPipeline.py
```

macOS or Linux:

```text
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[visualize]"
.venv/bin/python FREQNESS_MainPipeline.py
```

For a guided interactive workflow, open
`python/FREQNESS_NotebookPipeline.ipynb` in Jupyter or VS Code and run its
cells in order. It exposes intermediate results for one selected participant
within a chosen condition.

For automated processing of every condition, edit Block 0 near the top of
`python/FREQNESS_MainPipeline.py` and run the script. `RUN_CORE = True` and
`RUN_GROUP_ANALYSES = False` are the defaults. `NETWORK_OPTIONS["ncomps"] = 5`
retains five components while preserving the full eigenspectrum. Python saves
complete FREQ arrays and metadata in `.npz` archives; these are not MATLAB
`.mat` files. Load one with `FREQNESS_LoadCore(path)` for subsequent analyses.
Both interfaces initialize the Python package themselves; users do not run
`FREQNESS_Startup.py` separately. If `FREQNESS_Data` has no condition folders,
the toolbox displays a warning and performs no analysis.

Do not run `python/freqness/FREQNESS_MainPipeline.py` directly. That file is
the internal package engine and uses package-relative imports. Always run the
user-facing `python/FREQNESS_MainPipeline.py` script from the `python` folder as
shown above.

For individual functions, installation options, and Python examples, see the
[Python README](python/README.md).

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
- Leonardo Bonetti — Center for Music in the Brain, Aarhus University; Centre
  for Eudaimonia and Human Flourishing, Linacre College, University of Oxford —
  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk

Python-specific authorship and contact information is provided in
`python/README.md` and the Python source headers.
