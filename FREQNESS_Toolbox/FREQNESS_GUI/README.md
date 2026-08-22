# FREQ-NESS GUI

This folder is a self-contained graphical layer over the existing
`FREQNESS_Functions`. No GUI file requires edits to the toolbox startup,
pipeline, or numerical functions.

## Launch

Add this folder to the MATLAB path, or make it the current folder, and run:

```matlab
FREQNESS_GUI
```

The launcher adds only the existing function and dependency folders needed by
the current GUI session.

## Two-page workflow

### Page 1 — Network estimation

1. **Data and anatomical import**
   - Drop or browse to one folder named `Dataset*`; dropping its participant
     MAT-files also resolves their common folder.
   - Each top-level `.mat` file is one participant.
   - Each file must contain exactly one real numeric voxels-by-time matrix.
   - The bundled MNI coordinate file loads automatically when present. A custom
     finite N-by-3 coordinate matrix can be dropped or browsed independently,
     and its black-dot 3-D preview can be rotated by dragging.

2. **Core network estimation**
   - Mandatory inputs: frequency range, frequency step, and sampling rate.
   - The two-ended slider snaps to the selected step. Exact endpoint fields and
     an aligned frequency-domain filter-bank plot show the Gaussian filters
     that will be sent to the backend. Curves progress from light to dark as
     centre frequency increases.
   - Optional inputs mirror `FREQNESS_NetworkEstimation`.
   - Participants are processed independently and persisted immediately.

### Page 2 — Secondary analyses

3. **Analysis context and configuration**
   - Drop or browse to the corresponding `FREQ_Networks*` folder.
   - Secondary outputs use one subfolder per backend function inside the
     derived `FREQ_Analyses*` folder.
   - A grouped analysis tree exposes one function at a time.
   - Mandatory fields stay visible; optional settings remain collapsed until
     requested.
   - Every frequency input uses the imported `FREQ.frex` grid. Range sliders
     start at the full available interval and can be narrowed to a contiguous
     subset; single-frequency sliders snap to one exact available frequency.
   - Readiness checks identify missing FREQ results, MNI coordinates, events,
     or original source data before execution.
   - The Run Analysis button combines the selected participant results for
     one group-capable backend call, then persists the group output and each
     participant's extracted numerical output independently.

## Output contract

The text following the `Dataset` prefix is preserved:

```text
Dataset_1/
    sub-001.mat
    sub-002.mat

FREQ_Networks_1/
    sub-001_FREQ.mat
    sub-002_FREQ.mat
    FREQNESS_Manifest.mat

FREQ_Analyses_1/
    Visualizer/
        Group_Visualizer.mat
        sub-001_Visualizer.mat
        sub-002_Visualizer.mat
        FREQNESS_AnalysisManifest.mat
        Figures/
    EntropyLandscape/
        Group_EntropyLandscape.mat
        sub-001_EntropyLandscape.mat
        sub-002_EntropyLandscape.mat
        FREQNESS_AnalysisManifest.mat
        Figures/
    ExponentialDK/
    FreqGradients/
    CompGradients/
    CrossCoupling/
    BackProjection/
    NetworkRemoval/
    InducedResponses/
```

Each participant result contains:

- `FREQ`: the output of `FREQNESS_NetworkEstimation` for that participant;
- `participant`: stable source-file and participant metadata;
- `configuration`: the complete configuration used for the run.

`FREQNESS_Manifest.mat` indexes every participant and records whether its run
completed, failed, or was skipped because an existing output was retained.
This permits safe resumption of long analyses.

Each secondary-analysis function folder contains:

- `Group_<Function>.mat`: the complete backend output for the selected group;
- `<participant>_<Function>.mat`: that participant's extracted numerical
  output, participant metadata, and complete GUI configuration;
- `FREQNESS_AnalysisManifest.mat`: execution status and paths for every output;
- `Figures/`: captured MATLAB figures in editable FIG and PNG formats.

Outputs are written through temporary files and finalized atomically, so an
interrupted or failed run does not replace an existing completed MAT-file with
a partially written result.

## Drag and drop

Native OS file/folder drag-and-drop is enabled by the vendored MIT-licensed
`uiFileDnD` utility. Browse buttons remain available as a platform fallback.

## Current milestone

The data-import, participant-wise network estimation, and secondary-analysis
paths are functional. The Secondary Analyses page imports network-result
context, provides grouped function selection, builds and validates
function-specific configurations, runs the existing numerical functions, and
persists participant, group, figure, and manifest outputs.
