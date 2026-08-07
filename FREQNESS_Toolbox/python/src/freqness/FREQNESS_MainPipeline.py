"""Configurable folder-based FREQ-NESS analysis pipeline."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence
import warnings

import numpy as np
from numpy.typing import NDArray

from .FREQNESS_BackProjection import FREQNESS_BackProjection
from .FREQNESS_CompGradients import FREQNESS_CompGradients
from .FREQNESS_CrossCoupling import FREQNESS_CrossCoupling
from .FREQNESS_EntropyLandscape import FREQNESS_EntropyLandscape
from .FREQNESS_ExponentialDK import FREQNESS_ExponentialDK
from .FREQNESS_FreqGradients import FREQNESS_FreqGradients
from .FREQNESS_NetworkEstimation import FREQNESS_NetworkEstimation
from .FREQNESS_NetworkRemoval import FREQNESS_NetworkRemoval
from .FREQNESS_Startup import FREQNESS_Startup
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
    """Settings corresponding to the editable MATLAB pipeline section."""

    toolbox_root: str | Path | None = None
    frex: FloatArray = field(
        default_factory=lambda: np.arange(1.2, 24.0 + 0.6, 1.2)
    )
    srate: float = 250.0
    analyses: tuple[str, ...] = ALL_ANALYSES
    network_options: dict[str, Any] = field(default_factory=dict)
    plot_all: bool = False
    show: bool = True
    save_nifti: bool = True
    visualizer_output_directory: str | Path | None = None
    landscape_frex: FloatArray | None = None
    landscape_ncomps: int = 10
    pattern_frex: FloatArray = field(
        default_factory=lambda: np.array([2.0, 8.0])
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


@dataclass(slots=True)
class FREQNESSPipelineConditionResult:
    """In-memory outputs for one dataset folder."""

    name: str
    outputs: dict[str, Any] = field(default_factory=dict)


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


def _condition_names(root: Path, count: int) -> list[str]:
    data_root = root / "FREQNESS_Data"
    names = sorted(path.name for path in data_root.iterdir() if path.is_dir())
    if len(names) != count:
        raise RuntimeError(
            "Dataset folders changed while FREQNESS_Startup was loading them"
        )
    return names


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
    toolbox_root: Path,
) -> Path:
    if config.visualizer_output_directory is None:
        return toolbox_root
    return Path(config.visualizer_output_directory).expanduser().resolve()


def FREQNESS_MainPipeline(
    config: FREQNESSPipelineConfig | None = None,
) -> FREQNESSPipelineResult:
    """Run selected FREQ-NESS sections across all dataset folders.

    Network estimation is always computed as the prerequisite for every
    non-empty dataset, but is included in the returned outputs only when its
    section is selected. All numerical analysis results remain in memory. The
    only optional filesystem output is the Visualizer's established NIfTI
    export, controlled by save_nifti.
    """
    if config is None:
        config = FREQNESSPipelineConfig()
    selected = _validate_config(config)
    all_data, MNI, root = FREQNESS_Startup(config.toolbox_root)
    names = _condition_names(root, len(all_data))
    results: list[FREQNESSPipelineConditionResult] = []

    for condition, data in zip(names, all_data, strict=True):
        condition_result = FREQNESSPipelineConditionResult(name=condition)
        results.append(condition_result)
        if data is None:
            warnings.warn(
                f"Skipping {condition}: no participant .mat files found.",
                UserWarning,
                stacklevel=2,
            )
            continue

        FREQ = FREQNESS_NetworkEstimation(
            data,
            np.asarray(config.frex, dtype=float),
            config.srate,
            **dict(config.network_options),
        )
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
                    "path_output": _visualizer_output_root(config, root),
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
            analyzed_data = np.asarray(data, dtype=float)[
                :, : FREQ.ts.shape[1], ...
            ]
            data_clean, removed_activity = FREQNESS_NetworkRemoval(
                FREQ,
                analyzed_data,
                config.netrem_freq2remove,
                comps2remove=config.netrem_comps2remove,
            )
            network_removal_output: dict[str, Any] = {
                "dataClean": data_clean,
                "removedActivity": removed_activity,
            }
            if config.netrem_plot_landscape:
                FREQ_clean = FREQNESS_NetworkEstimation(
                    data_clean,
                    np.asarray(config.frex, dtype=float),
                    config.srate,
                    **dict(config.network_options),
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

    return FREQNESSPipelineResult(
        path_home=root,
        MNI=MNI,
        conditions=results,
    )


def main(
    argv: Sequence[str] | None = None,
    *,
    default_toolbox_root: str | Path | None = None,
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
        help="Exact function name to run; repeat as needed. Default: all.",
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
    arguments = parser.parse_args(argv)
    config = FREQNESSPipelineConfig(
        toolbox_root=arguments.toolbox_root,
        analyses=(
            tuple(arguments.analysis)
            if arguments.analysis is not None
            else ALL_ANALYSES
        ),
        show=not arguments.no_show,
        save_nifti=not arguments.no_nifti,
    )
    result = FREQNESS_MainPipeline(config)
    completed = sum(bool(condition.outputs) for condition in result.conditions)
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
