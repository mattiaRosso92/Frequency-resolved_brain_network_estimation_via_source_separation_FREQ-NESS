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
    └── FREQNESS_MainPipeline.py
```

The bundled `FREQNESS_MNI_Coordinates` file corresponds to the 3559-voxel
MNI152 8-mm grid. Replace it when using a different source space or voxel
ordering. Matching coordinates are required for spatial visualizations and
gradient modelling.

## MATLAB quick start

1. Open `FREQNESS_MainPipeline.m` in MATLAB.
2. Edit the settings in `DEFINE ANALYSIS SETTINGS`, including sampling rate,
   frequencies, and the analyses to run.
3. Run the script. It initializes the toolbox automatically.

Individual functions are available in `FREQNESS_Functions`. See
`README_Toolbox.rtf` and the function headers for detailed MATLAB guidance.

When MATLAB `FREQNESS_NetworkEstimation` receives no explicit `fwidth`, its
default `logarithmic` schedule varies spectral selectivity smoothly from
Q = 7 at the lowest requested frequency to Q = 3.5 at the highest. The
`linear` option uses constant Q = 5 (`FWHM = frequency / 5`). These automatic
widths depend on frequencies in Hz and their numerical range, not on the
number or spacing of frequency bins. A scalar `fwidth` still applies one
constant width to every filter, while a vector supplies one explicit width
per frequency.

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
cells in order. It exposes intermediate results for one selected condition.

For automated processing of every condition, edit the settings near the top
of `python/FREQNESS_MainPipeline.py` and run the script. Both interfaces
initialize the Python package themselves; users do not run
`FREQNESS_Startup.py` separately. If `FREQNESS_Data` has no condition folders,
the toolbox displays a warning and performs no analysis.

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
