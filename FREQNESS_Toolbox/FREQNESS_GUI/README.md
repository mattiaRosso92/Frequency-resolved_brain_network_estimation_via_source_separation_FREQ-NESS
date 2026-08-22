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

## Three-layer workflow

1. **Data import**
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
     an aligned FWHM plot show the vector that will be sent to the backend.
   - Optional inputs mirror `FREQNESS_NetworkEstimation`.
   - Participants are processed independently and persisted immediately.

3. **Secondary analyses**
   - Drop or browse to the corresponding `FREQ_Networks*` folder.
   - Secondary outputs use one subfolder per backend function inside the
     derived `FREQ_Analyses*` folder.

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
    EntropyLandscape/
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

## Drag and drop

Native OS file/folder drag-and-drop is enabled by the vendored MIT-licensed
`uiFileDnD` utility. Browse buttons remain available as a platform fallback.

## Current milestone

The data-import and participant-wise network-estimation path is functional.
The secondary layer currently defines and validates its input/output contract;
function-specific configuration panels and executors are the next milestone.
