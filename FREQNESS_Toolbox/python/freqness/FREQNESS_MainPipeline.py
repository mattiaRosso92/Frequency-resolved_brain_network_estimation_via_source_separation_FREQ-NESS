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
from types import SimpleNamespace
from typing import Any, Mapping, Sequence
import warnings

if __name__ == "__main__" and not __package__:
    raise SystemExit(
        "Do not run freqness/FREQNESS_MainPipeline.py directly. "
        "Open a terminal in the python folder and run: "
        "python FREQNESS_MainPipeline.py"
    )

import numpy as np
from numpy.typing import NDArray

from .FREQNESS_BackProjection import FREQNESS_BackProjection
from .FREQNESS_CompGradients import FREQNESS_CompGradients
from .FREQNESS_CrossCoupling import FREQNESS_CrossCoupling
from .FREQNESS_CoreIO import (
    FREQNESS_LoadCore,
    FREQNESS_LoadCoreFields,
    FREQNESS_SaveCore,
    _core_matches,
    _core_participant,
)
from .FREQNESS_EntropyLandscape import FREQNESS_EntropyLandscape
from .FREQNESS_ExponentialDK import FREQNESS_ExponentialDK
from .FREQNESS_FreqGradients import FREQNESS_FreqGradients
from .FREQNESS_LandscapeComparison import FREQNESS_LandscapeComparison
from .FREQNESS_NetworkEstimation import FREQNESS_NetworkEstimation
from .FREQNESS_NetworkRemoval import FREQNESS_NetworkRemoval
from .FREQNESS_SpatialComparison import FREQNESS_SpatialComparison
from .FREQNESS_StabilityIndex import FREQNESS_StabilityIndex
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
STABILITY_INDEX = "FREQNESS_StabilityIndex"
BACK_PROJECTION = "FREQNESS_BackProjection"
NETWORK_REMOVAL = "FREQNESS_NetworkRemoval"
SPATIAL_COMPARISON = "FREQNESS_SpatialComparison"
LANDSCAPE_COMPARISON = "FREQNESS_LandscapeComparison"

ALL_ANALYSES = (
    NETWORK_ESTIMATION,
    VISUALIZER,
    ENTROPY_LANDSCAPE,
    EXPONENTIAL_DK,
    FREQ_GRADIENTS,
    COMP_GRADIENTS,
    CROSS_COUPLING,
    STABILITY_INDEX,
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
    stability_which_comp: NDArray[np.int64] = field(
        default_factory=lambda: np.array([1], dtype=np.int64)
    )
    stability_frex2model: FloatArray | None = None
    stability_fwhm: FloatArray | float | None = None
    stability_median_window: float = 0.4
    stability_time_window: FloatArray | None = None
    stability_keep_frequency_sliding: bool = False
    stability_plot_avg: bool = True
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
    run_spatial_comparison: bool = False
    run_landscape_comparison: bool = False
    comparison_conditions: tuple[str, str] | None = None
    comparison_design: str = "within"
    comparison_plot_results: bool = True
    spatialcomparison_test: str | None = None
    spatialcomparison_frex: FloatArray | str = field(
        default_factory=lambda: np.array([2.4, 8.4])
    )
    spatialcomparison_comps: NDArray[np.int64] = field(
        default_factory=lambda: np.array([1], dtype=np.int64)
    )
    spatialcomparison_cluster_forming_alpha: float = 0.001
    spatialcomparison_alpha: float = 0.05
    spatialcomparison_n_permutations: int | str = "auto"
    landscapecomparison_test: str | None = None
    landscapecomparison_measure: str = "evals"
    landscapecomparison_frex: FloatArray | str | None = None
    landscapecomparison_comps: NDArray[np.int64] = field(
        default_factory=lambda: np.array([1], dtype=np.int64)
    )
    landscapecomparison_correction: str = "fdr"
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
    comparisons: dict[str, Any] = field(default_factory=dict)


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
        (config.stability_keep_frequency_sliding, "stability_keep_frequency_sliding"),
        (config.stability_plot_avg, "stability_plot_avg"),
        (config.run_spatial_comparison, "run_spatial_comparison"),
        (config.run_landscape_comparison, "run_landscape_comparison"),
        (config.comparison_plot_results, "comparison_plot_results"),
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


def _comparison_design(design: str) -> str:
    key = str(design).lower()
    if key in {"within", "paired"}:
        return "within"
    if key in {"between", "independent"}:
        return "between"
    raise ValueError("config.comparison_design must be 'within' or 'between'")


def _comparison_condition_indices(
    folder_data: Sequence[Mapping[str, Any]],
    requested: tuple[str, str] | None,
) -> tuple[int, int]:
    if requested is None:
        populated = [
            index for index, folder in enumerate(folder_data) if folder["files_mat"]
        ]
        if len(populated) < 2:
            raise ValueError(
                "At least two non-empty conditions are required for comparisons"
            )
        return populated[0], populated[1]
    if (
        not isinstance(requested, (tuple, list))
        or len(requested) != 2
        or not all(isinstance(name, str) and name for name in requested)
    ):
        raise ValueError(
            "config.comparison_conditions must contain two condition names"
        )
    lookup = {str(folder["name"]).casefold(): index for index, folder in enumerate(folder_data)}
    keys = tuple(name.casefold() for name in requested)
    if keys[0] == keys[1]:
        raise ValueError("comparison_conditions must identify two different conditions")
    missing = [name for name, key in zip(requested, keys, strict=True) if key not in lookup]
    if missing:
        raise ValueError("Unknown comparison condition(s): " + ", ".join(missing))
    return lookup[keys[0]], lookup[keys[1]]


def _comparison_records(
    folder: Mapping[str, Any],
    condition_result: FREQNESSPipelineConditionResult,
    network_settings: Mapping[str, Any],
) -> list[tuple[str, Path]]:
    records: list[tuple[str, Path]] = []
    identifiers: set[str] = set()
    for source_file, core_file in zip(
        folder["files_mat"], condition_result.core_files, strict=True
    ):
        if not core_file.is_file():
            raise FileNotFoundError(
                f"Core output missing: {core_file}. Run Block 1 first."
            )
        source_info = source_file.stat()
        expected = {
            "condition": folder["name"],
            "name": source_file.stem,
            "source_file": str(source_file),
            "source_size": folder["data_size"],
            "source_bytes": source_info.st_size,
            "source_mtime_ns": source_info.st_mtime_ns,
        }
        if not _core_matches(core_file, expected, dict(network_settings)):
            raise ValueError(
                f"Core output does not match the current inputs or settings: {core_file}. "
                "Run Block 1 with overwrite_core=True."
            )
        participant = _core_participant(core_file)
        identifier = str(participant.get("name", ""))
        key = identifier.casefold()
        if not identifier or key in identifiers:
            raise ValueError(
                f"Condition {folder['name']} contains missing or duplicate participant identifiers"
            )
        identifiers.add(key)
        records.append((identifier, core_file))
    return records


def _match_comparison_records(
    records_A: list[tuple[str, Path]],
    records_B: list[tuple[str, Path]],
    design: str,
) -> tuple[list[tuple[str, Path]], list[tuple[str, Path]]]:
    if design == "between":
        return records_A, records_B
    lookup_B = {
        identifier.casefold(): (identifier, core_file)
        for identifier, core_file in records_B
    }
    identifiers_A = [identifier.casefold() for identifier, _ in records_A]
    if set(identifiers_A) != set(lookup_B):
        only_A = sorted(set(identifiers_A).difference(lookup_B))
        only_B = sorted(set(lookup_B).difference(identifiers_A))
        raise ValueError(
            "Within-subject comparisons require identical participant identifiers "
            f"in both conditions. Only in A: {only_A}; only in B: {only_B}"
        )
    ordered_B = [lookup_B[identifier] for identifier in identifiers_A]
    return records_A, ordered_B


def _comparison_frequencies(
    available: NDArray[np.float64],
    requested: Any,
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    frequencies = np.asarray(available, dtype=float).reshape(-1)
    if requested is None or (isinstance(requested, str) and requested.lower() == "all"):
        indices = np.arange(frequencies.size, dtype=np.int64)
        return indices, frequencies.copy()
    if isinstance(requested, str):
        raise ValueError("Comparison frequencies must be numeric, None, or 'all'")
    values = np.asarray(requested, dtype=float).reshape(-1)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("Comparison frequencies must contain finite values")
    if np.any(values < frequencies[0]) or np.any(values > frequencies[-1]):
        raise ValueError("Comparison frequencies fall outside the CORE frequency range")
    tolerance = np.inf if frequencies.size == 1 else np.median(np.diff(frequencies)) / 2.0
    indices = np.empty(values.size, dtype=np.int64)
    for index, value in enumerate(values):
        distance = np.abs(frequencies - value)
        nearest = np.flatnonzero(np.isclose(distance, np.min(distance)))
        if nearest.size != 1 or distance[nearest[0]] > tolerance + np.finfo(float).eps:
            raise ValueError(f"No unique CORE frequency matches {value:g} Hz")
        indices[index] = nearest[0]
    if np.unique(indices).size != indices.size:
        raise ValueError("Two requested comparison frequencies resolve to the same bin")
    resolved = frequencies[indices]
    if not np.allclose(values, resolved, atol=1e-12, rtol=0.0):
        warnings.warn(
            "Some comparison frequencies are unavailable; using closest CORE bins.",
            UserWarning,
            stacklevel=3,
        )
    return indices, resolved


def _comparison_components(requested: Any, available: int) -> NDArray[np.int64]:
    values = np.asarray(requested, dtype=float).reshape(-1)
    if (
        values.size == 0
        or not np.all(np.isfinite(values))
        or not np.all(values == np.round(values))
        or np.any(values < 1)
        or np.any(values > available)
        or np.unique(values).size != values.size
    ):
        raise ValueError(
            f"Comparison components must be unique integers between 1 and {available}"
        )
    return values.astype(np.int64)


def _assemble_core_comparison(
    records: Sequence[tuple[str, Path]],
    field_name: str,
    requested_frex: Any,
    requested_components: Any,
) -> tuple[SimpleNamespace, NDArray[np.int64]]:
    if not records:
        raise ValueError("A comparison condition contains no participant CORE files")
    first = FREQNESS_LoadCoreFields(records[0][1], (field_name, "frex"))
    available_frex = np.asarray(first["frex"], dtype=float).reshape(-1)
    indices, resolved_frex = _comparison_frequencies(available_frex, requested_frex)
    values = np.asarray(first[field_name])
    if field_name == "pats":
        if values.ndim != 4 or values.shape[-1] != 1:
            raise ValueError("Participant CORE pats must end in one participant axis")
        components = _comparison_components(requested_components, values.shape[1])
        selected = values[:, components - 1, :, :][:, :, indices, 0]
    elif field_name == "evals":
        if values.ndim != 3 or values.shape[-1] != 1:
            raise ValueError("Participant CORE evals must end in one participant axis")
        components = _comparison_components(requested_components, values.shape[0])
        selected = values[components - 1, :, :][:, indices, 0]
    else:
        raise ValueError("Selective CORE comparisons support only pats or evals")
    combined = np.empty((*selected.shape, len(records)), dtype=selected.dtype)
    combined[..., 0] = selected
    del first, values, selected

    for subi, (_, core_file) in enumerate(records[1:], start=1):
        participant = FREQNESS_LoadCoreFields(core_file, (field_name, "frex"))
        if not np.allclose(
            np.asarray(participant["frex"], dtype=float).reshape(-1),
            available_frex,
            atol=1e-10,
            rtol=0.0,
        ):
            raise ValueError("Participant CORE frequency grids do not match")
        participant_values = np.asarray(participant[field_name])
        if field_name == "pats":
            combined[..., subi] = participant_values[:, components - 1, :, :][
                :, :, indices, 0
            ]
        else:
            combined[..., subi] = participant_values[components - 1, :, :][
                :, indices, 0
            ]
        del participant, participant_values
    return SimpleNamespace(**{field_name: combined, "frex": resolved_frex}), components


def _assemble_output_comparison(
    output: Any,
    field_name: str,
    axis_name: str,
    requested_frex: Any,
    requested_components: Any,
    participant_order: NDArray[np.int64],
) -> tuple[SimpleNamespace, NDArray[np.int64]]:
    """Build one thin statistical input from a secondary group output."""
    if not hasattr(output, field_name) or not hasattr(output, axis_name):
        raise ValueError(
            f"Secondary output does not contain {field_name!r} and {axis_name!r}"
        )
    available_frex = np.asarray(getattr(output, axis_name), dtype=float).reshape(-1)
    indices, resolved_frex = _comparison_frequencies(available_frex, requested_frex)
    values = np.asarray(getattr(output, field_name))
    if values.ndim not in (2, 3) or values.shape[-1] != participant_order.size:
        raise ValueError(
            f"Secondary output {field_name} does not match the participant count"
        )
    if values.ndim == 2:
        if values.shape[0] != available_frex.size:
            raise ValueError(
                f"Secondary output {field_name} must be [nFrex x nParticipants]"
            )
        components = _comparison_components(requested_components, 1)
        selected = values[indices, :][:, participant_order]
    elif values.ndim == 3:
        if values.shape[1] != available_frex.size:
            raise ValueError(
                f"Secondary output {field_name} must be "
                "[nFeatures x nFrex x nParticipants]"
            )
        components = _comparison_components(requested_components, values.shape[0])
        selected = values[components - 1, :, :][:, indices, :][..., participant_order]
    else:
        raise ValueError(
            f"Secondary output {field_name} must be two- or three-dimensional"
        )
    return (
        SimpleNamespace(**{field_name: selected, axis_name: resolved_frex}),
        components,
    )


def _participant_order(
    original: Sequence[tuple[str, Path]],
    matched: Sequence[tuple[str, Path]],
) -> NDArray[np.int64]:
    """Map matched participant records back to the secondary-output order."""
    lookup = {core_file: index for index, (_, core_file) in enumerate(original)}
    return np.asarray([lookup[core_file] for _, core_file in matched], dtype=np.int64)


def FREQNESS_MainPipeline(
    config: FREQNESSPipelineConfig | None = None,
) -> FREQNESSPipelineResult:
    """Estimate and save one participant at a time, then optionally analyse groups.

    Block 0 inspects input files and defines the network settings. Block 1 saves
    complete participant results in FREQNESS_Results/CORE/<condition> and releases
    their arrays. Block 2 optionally combines saved FREQ results. Block 3 builds
    thin statistical inputs containing only the requested landscape fields,
    frequency axes, and participant dimension. Numerical APIs are unchanged.
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
    comparisons: dict[str, Any] = {}

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
                ENTR = FREQNESS_EntropyLandscape(
                    FREQ,
                    plot_all=config.plot_all,
                    show=config.show,
                )
                condition_result.outputs[ENTROPY_LANDSCAPE] = ENTR

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

            if STABILITY_INDEX in selected:
                condition_result.outputs[STABILITY_INDEX] = (
                    FREQNESS_StabilityIndex(
                        FREQ,
                        which_comp=config.stability_which_comp,
                        frex2model=config.stability_frex2model,
                        fwhm=config.stability_fwhm,
                        median_window=config.stability_median_window,
                        time_window=config.stability_time_window,
                        keep_frequency_sliding=(
                            config.stability_keep_frequency_sliding
                        ),
                        plot_avg=config.stability_plot_avg,
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

    # ========================================================================
    # BLOCK 3 - SELECTIVE CONDITION/GROUP COMPARISONS
    # ========================================================================
    # Block 3 compares selected spatial patterns or frequency landscapes
    # between two conditions or groups. It supports paired and independent
    # designs, matches paired participants by identifier, and returns the
    # statistical results separately from the condition outputs.
    #
    # For low memory use, each sample is represented by one temporary structure
    # containing only the selected landscape, frequency axis, and participant
    # dimension. Complete group ts/evecs arrays are never assembled here.
    if config.run_spatial_comparison or config.run_landscape_comparison:
        condition_A, condition_B = _comparison_condition_indices(
            folder_data, config.comparison_conditions
        )
        folder_A = folder_data[condition_A]
        folder_B = folder_data[condition_B]
        result_A = results[condition_A]
        result_B = results[condition_B]
        design = _comparison_design(config.comparison_design)
        original_records_A = _comparison_records(folder_A, result_A, NETWORK_SETTINGS)
        original_records_B = _comparison_records(folder_B, result_B, NETWORK_SETTINGS)
        records_A, records_B = _match_comparison_records(
            original_records_A, original_records_B, design
        )
        order_A = _participant_order(original_records_A, records_A)
        order_B = _participant_order(original_records_B, records_B)
        sample_labels = (str(folder_A["name"]), str(folder_B["name"]))

        if config.run_spatial_comparison:
            if MNI is None:
                raise ValueError(
                    "FREQNESS_SpatialComparison requires matching MNI coordinates"
                )
            compact_A, original_components = _assemble_core_comparison(
                records_A,
                "pats",
                config.spatialcomparison_frex,
                config.spatialcomparison_comps,
            )
            compact_B, components_B = _assemble_core_comparison(
                records_B,
                "pats",
                config.spatialcomparison_frex,
                config.spatialcomparison_comps,
            )
            if not np.array_equal(original_components, components_B):
                raise RuntimeError("Internal spatial component selection mismatch")
            spatial_stats = FREQNESS_SpatialComparison(
                compact_A,
                compact_B,
                MNI,
                frex2compare=compact_A.frex,
                comps2compare=np.arange(1, original_components.size + 1),
                design=design,
                test=config.spatialcomparison_test,
                condition_labels=sample_labels,
                cluster_forming_alpha=config.spatialcomparison_cluster_forming_alpha,
                alpha=config.spatialcomparison_alpha,
                n_permutations=config.spatialcomparison_n_permutations,
                plot_results=config.comparison_plot_results,
                show=config.show,
            )
            spatial_stats.comps = original_components.copy()
            for cluster in spatial_stats.clusters:
                cluster.comp = int(original_components[cluster.comp - 1])
            comparisons[SPATIAL_COMPARISON] = spatial_stats
            del compact_A, compact_B, spatial_stats

        if config.run_landscape_comparison:
            requested_measure = str(config.landscapecomparison_measure)
            measure_key = requested_measure.lower()
            if measure_key == "evals":
                compact_A, original_components = _assemble_core_comparison(
                    records_A,
                    "evals",
                    config.landscapecomparison_frex,
                    config.landscapecomparison_comps,
                )
                compact_B, components_B = _assemble_core_comparison(
                    records_B,
                    "evals",
                    config.landscapecomparison_frex,
                    config.landscapecomparison_comps,
                )
                measure_name = "evals"
            else:
                secondary_sources = {
                    "h2": (ENTROPY_LANDSCAPE, "H2", "frex"),
                    "ed": (ENTROPY_LANDSCAPE, "ED", "frex"),
                    "amplitude_raw": (
                        CROSS_COUPLING,
                        "amplitude_raw",
                        "carrier_frex",
                    ),
                    "amplitude_normalized": (
                        CROSS_COUPLING,
                        "amplitude_normalized",
                        "carrier_frex",
                    ),
                    "stability_index": (
                        STABILITY_INDEX,
                        "stability_index",
                        "frex",
                    ),
                    "mean_frequency": (
                        STABILITY_INDEX,
                        "mean_frequency",
                        "frex",
                    ),
                }
                if measure_key not in secondary_sources:
                    raise ValueError(
                        "Block 3 can assemble evals, H2, ED, amplitude_raw, "
                        "amplitude_normalized, stability_index, or "
                        "mean_frequency. Pass other compact landscapes directly "
                        "to FREQNESS_LandscapeComparison."
                    )
                analysis_name, measure_name, axis_name = secondary_sources[measure_key]
                if (
                    analysis_name not in result_A.outputs
                    or analysis_name not in result_B.outputs
                ):
                    raise ValueError(
                        f"landscapecomparison_measure={requested_measure!r} requires "
                        f"{analysis_name} to be selected and Block 2 to be run."
                    )
                compact_A, original_components = _assemble_output_comparison(
                    result_A.outputs[analysis_name],
                    measure_name,
                    axis_name,
                    config.landscapecomparison_frex,
                    config.landscapecomparison_comps,
                    order_A,
                )
                compact_B, components_B = _assemble_output_comparison(
                    result_B.outputs[analysis_name],
                    measure_name,
                    axis_name,
                    config.landscapecomparison_frex,
                    config.landscapecomparison_comps,
                    order_B,
                )
            if not np.array_equal(original_components, components_B):
                raise RuntimeError("Internal landscape component selection mismatch")
            frequency_axis = (
                compact_A.carrier_frex
                if hasattr(compact_A, "carrier_frex")
                else compact_A.frex
            )
            landscape_stats = FREQNESS_LandscapeComparison(
                compact_A,
                compact_B,
                measure=measure_name,
                design=design,
                test=config.landscapecomparison_test,
                which_comp=np.arange(1, original_components.size + 1),
                frex2compare=frequency_axis,
                sample_labels=sample_labels,
                correction=config.landscapecomparison_correction,
                plot_results=config.comparison_plot_results,
                show=config.show,
            )
            landscape_stats.results[measure_name].comps = original_components.copy()
            comparisons[LANDSCAPE_COMPARISON] = landscape_stats
            del compact_A, compact_B, landscape_stats

        del original_records_A, original_records_B, records_A, records_B

    return FREQNESSPipelineResult(
        path_home=root,
        MNI=MNI,
        conditions=results,
        comparisons=comparisons,
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
    parser.add_argument("--spatial-comparison", action="store_true",
                        help="Run the selective CORE spatial comparison block.")
    parser.add_argument("--landscape-comparison", action="store_true",
                        help="Run the selective landscape comparison block.")
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
    if arguments.spatial_comparison:
        config.run_spatial_comparison = True
    if arguments.landscape_comparison:
        config.run_landscape_comparison = True
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
    "LANDSCAPE_COMPARISON",
    "NETWORK_ESTIMATION",
    "NETWORK_REMOVAL",
    "SPATIAL_COMPARISON",
    "STABILITY_INDEX",
    "VISUALIZER",
    "main",
]
