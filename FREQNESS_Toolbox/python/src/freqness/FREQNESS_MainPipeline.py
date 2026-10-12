"""Configurable folder-based FREQ-NESS analysis pipeline.

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

from __future__ import annotations

import argparse
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Mapping, Sequence
import warnings

import numpy as np
from numpy.typing import NDArray

from .FREQNESS_BackProjection import FREQNESS_BackProjection
from .FREQNESS_CompGradients import FREQNESS_CompGradients
from .FREQNESS_CrossCoupling import FREQNESS_CrossCoupling
from .FREQNESS_CoreIO import FREQNESS_LoadCore, FREQNESS_SaveCore, _core_matches
from .FREQNESS_EntropyLandscape import FREQNESS_EntropyLandscape
from .FREQNESS_ExponentialDK import FREQNESS_ExponentialDK
from .FREQNESS_FreqGradients import FREQNESS_FreqGradients
from .FREQNESS_NetworkEstimation import FREQNESS_NetworkEstimation
from .FREQNESS_NetworkRemoval import FREQNESS_NetworkRemoval
from .FREQNESS_Startup import FREQNESS_Startup, _load_participant
from .FREQNESS_Visualizer import FREQNESS_Visualizer


FloatArray = NDArray[np.float64]

NETWORK_ESTIMATION = "FREQNESS_NetworkEstimation"
VISUALIZER = "FREQNESS_Visualizer"
ENTROPY_LANDSCAPE = "FREQNESS_EntropyLandscape"
EXPONENTIAL_DK = "FREQNESS_ExponentialDK"
FREQ_GRADIENTS = "FREQNESS_FreqGradients"
COMP_GRADIENTS = "FREQNESS_CompGradients"
CROSS_COUPLING = "FREQNESS_CrossCoupling"
BACK_PROJECTION = "FREQNESS_BackProjection"
NETWORK_REMOVAL = "FREQNESS_NetworkRemoval"

ALL_ANALYSES = (
    NETWORK_ESTIMATION,
    VISUALIZER,
    ENTROPY_LANDSCAPE,
    EXPONENTIAL_DK,
    FREQ_GRADIENTS,
    COMP_GRADIENTS,
    CROSS_COUPLING,
    BACK_PROJECTION,
    NETWORK_REMOVAL,
)


@dataclass(slots=True)
class FREQNESSPipelineConfig:
    """Settings corresponding to the editable MATLAB pipeline section.

    ``network_options`` is forwarded to ``FREQNESS_NetworkEstimation``. Its
    automatic widths use the physical logarithmic schedule by default; pass
    ``filter='linear'`` for constant Q = 5, or ``fwidth`` for explicit widths.
    The pipeline retains five components by default; evals remains complete.
    CORE export runs by default. Set run_group_analyses=True to combine saved
    participant files and run the sections selected by analyses. Group outputs
    remain in memory. Set run_core=False to use existing matching CORE files.
    """

    toolbox_root: str | Path | None = None
    frex: FloatArray = field(
        default_factory=lambda: np.arange(1.2, 24.0 + 0.6, 1.2)
    )
    srate: float = 250.0
    analyses: tuple[str, ...] = ALL_ANALYSES
    network_options: dict[str, Any] = field(default_factory=lambda: {"ncomps": 5})
    plot_all: bool = False
    show: bool = True
    save_nifti: bool = True
    visualizer_output_directory: str | Path | None = None
    landscape_frex: FloatArray | None = None
    landscape_ncomps: int = 10
    pattern_frex: FloatArray = field(
        default_factory=lambda: np.array([2.4, 8.4])
    )
    pattern_ncomps: int = 1
    expdk_range2fit: FloatArray | None = field(
        default_factory=lambda: np.array([9.0, 20.0])
    )
    expdk_which_comp: int = 1
    freqgrad_frex2model: FloatArray | None = field(
        default_factory=lambda: np.array([8.0, 12.0])
    )
    freqgrad_comp2model: int = 1
    compgrad_freq2model: float | None = 10.0
    compgrad_comps2model: NDArray[np.int64] | None = field(
        default_factory=lambda: np.array([1, 5], dtype=np.int64)
    )
    lfo_freq: float = 2.0
    backproj_freq2project: float = 8.4
    backproj_comps2project: NDArray[np.int64] = field(
        default_factory=lambda: np.array([1], dtype=np.int64)
    )
    netrem_freq2remove: float = 8.4
    netrem_comps2remove: NDArray[np.int64] = field(
        default_factory=lambda: np.array([1], dtype=np.int64)
    )
    netrem_plot_landscape: bool = True
    netrem_landscape_ncomps: int = 3
    run_core: bool = True
    run_group_analyses: bool = False
    path_results: str | Path | None = None
    overwrite_core: bool = False


@dataclass(slots=True)
class FREQNESSPipelineConditionResult:
    """Participant CORE paths and optional in-memory outputs for one folder."""

    name: str
    outputs: dict[str, Any] = field(default_factory=dict)
    core_files: list[Path] = field(default_factory=list)


@dataclass(slots=True)
class FREQNESSPipelineResult:
    """Complete output of :func:`FREQNESS_MainPipeline`."""

    path_home: Path
    MNI: NDArray[np.number] | None
    conditions: list[FREQNESSPipelineConditionResult]


def _validate_config(config: FREQNESSPipelineConfig) -> tuple[str, ...]:
    if not isinstance(config, FREQNESSPipelineConfig):
        raise TypeError("config must be a FREQNESSPipelineConfig")
    analyses = tuple(config.analyses)
    if not analyses:
        raise ValueError("config.analyses cannot be empty")
    unknown = sorted(set(analyses).difference(ALL_ANALYSES))
    if unknown:
        raise ValueError(
            "Unknown analysis name(s): "
            + ", ".join(unknown)
            + ". Use exact FREQ-NESS function names."
        )
    if len(set(analyses)) != len(analyses):
        raise ValueError("config.analyses cannot contain duplicates")
    for value, name in (
        (config.run_core, "run_core"),
        (config.run_group_analyses, "run_group_analyses"),
        (config.overwrite_core, "overwrite_core"),
        (config.plot_all, "plot_all"),
        (config.show, "show"),
        (config.save_nifti, "save_nifti"),
        (config.netrem_plot_landscape, "netrem_plot_landscape"),
    ):
        if not isinstance(value, (bool, np.bool_)):
            raise TypeError(f"config.{name} must be a boolean")
    if not isinstance(config.network_options, Mapping):
        raise TypeError("config.network_options must be a mapping")
    return analyses


def _require_mni(
    MNI: NDArray[np.number] | None,
    analysis: str,
    condition: str,
) -> bool:
    if MNI is not None:
        return True
    warnings.warn(
        f"Skipping {analysis} for {condition}: no valid MNI coordinates loaded.",
        UserWarning,
        stacklevel=3,
    )
    return False


def _visualizer_output_root(
    config: FREQNESSPipelineConfig,
    path_results: Path,
) -> Path:
    if config.visualizer_output_directory is None:
        return path_results / "GROUP"
    return Path(config.visualizer_output_directory).expanduser().resolve()


def FREQNESS_MainPipeline(
    config: FREQNESSPipelineConfig | None = None,
) -> FREQNESSPipelineResult:
    """Estimate and save one participant at a time, then optionally analyse groups.

    Block 0 inspects input files and defines the network settings. Block 1 saves
    complete participant results in FREQNESS_Results/CORE/<condition> and releases
    their arrays. Block 2 is optional: it combines saved FREQ results and retains
    selected secondary outputs in memory. Numerical estimator APIs are unchanged.
    """
    # ========================================================================
    # BLOCK 0 - STARTUP AND ANALYSIS SETTINGS
    # ========================================================================
    if config is None:
        config = FREQNESSPipelineConfig()
    selected = _validate_config(config)
    _, MNI, root, folder_data = FREQNESS_Startup(
        config.toolbox_root, load_data=False, return_folders=True,
    )
    path_results = (
        root / "FREQNESS_Results" if config.path_results is None
        else Path(config.path_results).expanduser().resolve()
    )
    network_options = {
        "duration": None, "fwidth": None, "filter": "logarithmic",
        "regularisation": 0.01, "ncomps": 5, "bad_segments": None, "rescale": False,
        **dict(config.network_options),
    }
    NETWORK_SETTINGS = {"frex": np.asarray(config.frex), "srate": config.srate,
                        **network_options}
    results = [
        FREQNESSPipelineConditionResult(
            name=folder["name"],
            core_files=[
                path_results / "CORE" / folder["name"] / (path.stem + "_FREQ.npz")
                for path in folder["files_mat"]
            ],
        )
        for folder in folder_data
    ]

    # ========================================================================
    # BLOCK 1 - NETWORK ESTIMATION AND CORE EXPORT (one participant at a time)
    # ========================================================================
    if config.run_core:
        for folder, condition_result in zip(folder_data, results, strict=True):
            print(f"\nESTIMATING CONDITION/GROUP: {folder['name']}")
            files_mat = folder["files_mat"]
            if not files_mat:
                continue
            core_path = path_results / "CORE" / folder["name"]
            core_path.mkdir(parents=True, exist_ok=True)

            for subi, source_file in enumerate(files_mat):
                print(f"\nPARTICIPANT #{subi + 1}/{len(files_mat)}: {source_file.name}")
                output_file = condition_result.core_files[subi]
                source_info = source_file.stat()
                PARTICIPANT = {
                    "condition": folder["name"], "name": source_file.stem,
                    "source_file": str(source_file), "source_size": folder["data_size"],
                    "source_bytes": source_info.st_size,
                    "source_mtime_ns": source_info.st_mtime_ns,
                }
                if output_file.exists() and not config.overwrite_core:
                    if not _core_matches(output_file, PARTICIPANT, NETWORK_SETTINGS):
                        raise ValueError(
                            f"Existing core file does not match the current inputs or settings: "
                            f"{output_file}. Set overwrite_core=True or choose another path_results."
                        )
                    print("Core output already exists; skipping participant.")
                    continue

                # Retain the last participant dimension as a singleton.
                this_data = _load_participant(source_file, folder["data_size"])
                FREQ = FREQNESS_NetworkEstimation(
                    this_data[:, :, np.newaxis], np.asarray(config.frex, dtype=float),
                    config.srate, **network_options,
                )
                del this_data
                FREQNESS_SaveCore(
                    output_file, FREQ, PARTICIPANT, NETWORK_SETTINGS,
                    overwrite=config.overwrite_core,
                )
                print(f"Saved core output: {output_file}")
                del FREQ

    # ========================================================================
    # BLOCK 2 - OPTIONAL GROUP ANALYSES (combine saved participant FREQ outputs)
    # ========================================================================
    if config.run_group_analyses:
        for condi, (folder, condition_result) in enumerate(zip(folder_data, results, strict=True)):
            condition = folder["name"]
            print(f"\nANALYSING CONDITION/GROUP #{condi + 1}")
            files_mat = folder["files_mat"]
            nsubs = len(files_mat)
            if not nsubs:
                continue

            for subi, source_file in enumerate(files_mat):
                input_file = condition_result.core_files[subi]
                if not input_file.is_file():
                    raise FileNotFoundError(f"Core output missing: {input_file}. Run Block 1 first.")
                source_info = source_file.stat()
                PARTICIPANT = {
                    "condition": condition, "name": source_file.stem,
                    "source_file": str(source_file), "source_size": folder["data_size"],
                    "source_bytes": source_info.st_size,
                    "source_mtime_ns": source_info.st_mtime_ns,
                }
                if not _core_matches(input_file, PARTICIPANT, NETWORK_SETTINGS):
                    raise ValueError(
                        f"Core output does not match the current inputs or settings: {input_file}. "
                        "Run Block 1 with overwrite_core=True."
                    )
                this_freq = FREQNESS_LoadCore(input_file)
                # Allocate once, then restore participants along the last dimension.
                if subi == 0:
                    FREQ = replace(
                        this_freq,
                        evals=np.empty((*this_freq.evals.shape[:-1], nsubs)),
                        evecs=np.empty((*this_freq.evecs.shape[:-1], nsubs)),
                        pats=np.empty((*this_freq.pats.shape[:-1], nsubs)),
                        ts=np.empty((*this_freq.ts.shape[:-1], nsubs)),
                        scale_factors=np.empty(nsubs),
                    )
                FREQ.evals[..., subi] = this_freq.evals[..., 0]
                FREQ.evecs[..., subi] = this_freq.evecs[..., 0]
                FREQ.pats[..., subi] = this_freq.pats[..., 0]
                FREQ.ts[..., subi] = this_freq.ts[..., 0]
                FREQ.scale_factors[subi] = this_freq.scale_factors[0]
                del this_freq

            if NETWORK_ESTIMATION in selected:
                condition_result.outputs[NETWORK_ESTIMATION] = FREQ

            if VISUALIZER in selected and _require_mni(MNI, VISUALIZER, condition):
                landscape_frex = (
                    np.asarray(config.frex, dtype=float)
                    if config.landscape_frex is None
                    else np.asarray(config.landscape_frex, dtype=float)
                )
                visualization = FREQNESS_Visualizer(
                    FREQ,
                    {
                        "frex": landscape_frex,
                        "ncomps": config.landscape_ncomps,
                    },
                    {
                        "frex": np.asarray(config.pattern_frex, dtype=float),
                        "ncomps": config.pattern_ncomps,
                        "path_output": _visualizer_output_root(config, path_results) / condition,
                        "MNI_coords": MNI,
                    },
                    plot_all=config.plot_all,
                    save_nifti=config.save_nifti,
                    show=config.show,
                )
                condition_result.outputs[VISUALIZER] = visualization

            if ENTROPY_LANDSCAPE in selected:
                H2, ED = FREQNESS_EntropyLandscape(
                    FREQ,
                    plot_all=config.plot_all,
                    show=config.show,
                )
                condition_result.outputs[ENTROPY_LANDSCAPE] = {
                    "H2": H2,
                    "ED": ED,
                }

            if EXPONENTIAL_DK in selected:
                condition_result.outputs[EXPONENTIAL_DK] = FREQNESS_ExponentialDK(
                    FREQ,
                    which_comp=config.expdk_which_comp,
                    range2fit=config.expdk_range2fit,
                    plot_all=config.plot_all,
                    show=config.show,
                )

            if FREQ_GRADIENTS in selected and _require_mni(
                MNI, FREQ_GRADIENTS, condition
            ):
                condition_result.outputs[FREQ_GRADIENTS] = (
                    FREQNESS_FreqGradients(
                        FREQ,
                        MNI,
                        frex2model=config.freqgrad_frex2model,
                        comp2model=config.freqgrad_comp2model,
                        plot_all=config.plot_all,
                        show=config.show,
                    )
                )

            if COMP_GRADIENTS in selected and _require_mni(
                MNI, COMP_GRADIENTS, condition
            ):
                condition_result.outputs[COMP_GRADIENTS] = (
                    FREQNESS_CompGradients(
                        FREQ,
                        MNI,
                        freq2model=config.compgrad_freq2model,
                        comps2model=config.compgrad_comps2model,
                        plot_all=config.plot_all,
                        show=config.show,
                    )
                )

            if CROSS_COUPLING in selected:
                condition_result.outputs[CROSS_COUPLING] = (
                    FREQNESS_CrossCoupling(
                        FREQ,
                        config.lfo_freq,
                        MNI=MNI,
                        plot_all=config.plot_all,
                        show=config.show,
                    )
                )

            if BACK_PROJECTION in selected:
                condition_result.outputs[BACK_PROJECTION] = FREQNESS_BackProjection(
                    FREQ,
                    config.backproj_freq2project,
                    comps2project=config.backproj_comps2project,
                )

            if NETWORK_REMOVAL in selected:
                # Import raw data only for this optional analysis and condition.
                ntime = FREQ.ts.shape[1]
                analyzed_data = np.empty((FREQ.evecs.shape[0], ntime, nsubs))
                for subi, source_file in enumerate(files_mat):
                    this_data = _load_participant(source_file, folder["data_size"])
                    analyzed_data[:, :, subi] = this_data[:, :ntime]
                    del this_data
                data_clean, removed_activity = FREQNESS_NetworkRemoval(
                    FREQ,
                    analyzed_data,
                    config.netrem_freq2remove,
                    comps2remove=config.netrem_comps2remove,
                )
                del analyzed_data
                network_removal_output: dict[str, Any] = {
                    "dataClean": data_clean,
                    "removedActivity": removed_activity,
                }
                if config.netrem_plot_landscape:
                    FREQ_clean = FREQNESS_NetworkEstimation(
                        data_clean,
                        np.asarray(config.frex, dtype=float),
                        config.srate,
                        **network_options,
                    )
                    clean_visualization = FREQNESS_Visualizer(
                        FREQ_clean,
                        {
                            "frex": np.asarray(config.frex, dtype=float),
                            "ncomps": config.netrem_landscape_ncomps,
                        },
                        None,
                        save_nifti=False,
                        show=config.show,
                    )
                    network_removal_output.update(
                        {
                            "FREQ_clean": FREQ_clean,
                            "visualization": clean_visualization,
                        }
                    )
                condition_result.outputs[NETWORK_REMOVAL] = network_removal_output
            del FREQ

    return FREQNESSPipelineResult(
        path_home=root,
        MNI=MNI,
        conditions=results,
    )


def main(
    argv: Sequence[str] | None = None,
    *,
    default_toolbox_root: str | Path | None = None,
    default_config: FREQNESSPipelineConfig | None = None,
) -> int:
    """Command-line entry point used by the mirrored pipeline script."""
    parser = argparse.ArgumentParser(
        description=(
            "Run selected FREQ-NESS analyses over FREQNESS_Data folders. "
            "Repeat --analysis to run multiple sections."
        )
    )
    parser.add_argument(
        "--toolbox-root",
        type=Path,
        default=default_toolbox_root,
        help="Path to FREQNESS_Toolbox.",
    )
    parser.add_argument(
        "--analysis",
        action="append",
        choices=ALL_ANALYSES,
        help="Select group sections by exact function name; repeat as needed.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Create figures without opening interactive windows.",
    )
    parser.add_argument(
        "--no-nifti",
        action="store_true",
        help="Do not export NIfTI files from FREQNESS_Visualizer.",
    )
    parser.add_argument("--group-analyses", action="store_true",
                        help="Combine saved participant files and run group analyses.")
    parser.add_argument("--skip-core", action="store_true",
                        help="Use existing CORE files without estimating networks.")
    parser.add_argument("--overwrite-core", action="store_true",
                        help="Replace existing participant CORE files.")
    parser.add_argument("--results-directory", type=Path,
                        help="Override the FREQNESS_Results output directory.")
    arguments = parser.parse_args(argv)
    config = replace(default_config) if default_config is not None else FREQNESSPipelineConfig()
    if arguments.toolbox_root is not None:
        config.toolbox_root = arguments.toolbox_root
    if arguments.analysis is not None:
        config.analyses = tuple(arguments.analysis)
        config.run_group_analyses = any(name != NETWORK_ESTIMATION for name in config.analyses)
    if arguments.group_analyses:
        config.run_group_analyses = True
    if arguments.skip_core:
        config.run_core = False
    if arguments.overwrite_core:
        config.overwrite_core = True
    if arguments.results_directory is not None:
        config.path_results = arguments.results_directory
    if arguments.no_show:
        config.show = False
    if arguments.no_nifti:
        config.save_nifti = False
    result = FREQNESS_MainPipeline(config)
    completed = sum(bool(condition.core_files) for condition in result.conditions)
    print(
        f"FREQ-NESS pipeline completed for {completed} non-empty "
        f"condition(s)."
    )
    return 0


__all__ = [
    "ALL_ANALYSES",
    "BACK_PROJECTION",
    "COMP_GRADIENTS",
    "CROSS_COUPLING",
    "ENTROPY_LANDSCAPE",
    "EXPONENTIAL_DK",
    "FREQNESS_MainPipeline",
    "FREQNESSPipelineConditionResult",
    "FREQNESSPipelineConfig",
    "FREQNESSPipelineResult",
    "FREQ_GRADIENTS",
    "NETWORK_ESTIMATION",
    "NETWORK_REMOVAL",
    "VISUALIZER",
    "main",
]
