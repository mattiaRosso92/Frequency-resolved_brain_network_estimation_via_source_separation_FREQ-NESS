"""Configurable, comparison-ready FREQ-NESS analysis pipeline."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field, fields as dataclass_fields
from pathlib import Path
from typing import Any, Mapping, Sequence
import warnings

import numpy as np
from numpy.typing import NDArray
from scipy.io import loadmat, savemat

from .FREQNESS_CompGradients import FREQNESS_CompGradients
from .FREQNESS_CrossCoupling import FREQNESS_CrossCoupling
from .FREQNESS_EntropyLandscape import FREQNESS_EntropyLandscape
from .FREQNESS_ExponentialDK import FREQNESS_ExponentialDK
from .FREQNESS_FreqGradients import FREQNESS_FreqGradients
from .FREQNESS_NetworkEstimation import (
    FREQNESS_NetworkEstimation,
    FREQNESSResult,
)
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

ALL_ANALYSES = (
    NETWORK_ESTIMATION,
    VISUALIZER,
    ENTROPY_LANDSCAPE,
    EXPONENTIAL_DK,
    FREQ_GRADIENTS,
    COMP_GRADIENTS,
    CROSS_COUPLING,
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
    save_outputs: bool = True
    save_nifti: bool = True
    output_directory: str | Path | None = None
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


@dataclass(slots=True)
class FREQNESSPipelineConditionResult:
    """Outputs and exported files for one dataset folder."""

    name: str
    outputs: dict[str, Any] = field(default_factory=dict)
    files: dict[str, Path] = field(default_factory=dict)


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
        (config.save_outputs, "save_outputs"),
        (config.save_nifti, "save_nifti"),
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


def _output_root(config: FREQNESSPipelineConfig, root: Path) -> Path:
    if config.output_directory is not None:
        return Path(config.output_directory).expanduser().resolve()
    return root / "python" / "FREQNESS_ComparisonOutputs" / "Python"


def _checkpoint_path(output_root: Path, condition: str) -> Path:
    return output_root / condition / f"{NETWORK_ESTIMATION}_python.mat"


def _freq_payload(FREQ: FREQNESSResult, condition: str) -> dict[str, Any]:
    return {
        "condition": condition,
        "evals": FREQ.evals,
        "evecs": FREQ.evecs,
        "pats": FREQ.pats,
        "ts": FREQ.ts,
        "frex": FREQ.frex,
        "fwhm": FREQ.fwhm,
        "srate": FREQ.srate,
        "duration": FREQ.duration,
        "bad_segments": FREQ.bad_segments,
        "regularisation": FREQ.regularisation,
        "scale_factors": FREQ.scale_factors,
    }


def _scalar(contents: Mapping[str, Any], name: str) -> float:
    values = np.asarray(contents[name], dtype=float)
    if values.size != 1:
        raise ValueError(f"Invalid {name} in saved network-estimation output")
    return float(values.reshape(-1)[0])


def _load_freq_checkpoint(path: Path) -> FREQNESSResult:
    if not path.is_file():
        raise FileNotFoundError(
            f"Required core checkpoint not found: {path}. Run "
            f"{NETWORK_ESTIMATION} for this condition first."
        )
    contents = loadmat(path)
    required = {
        "evals",
        "evecs",
        "pats",
        "ts",
        "frex",
        "fwhm",
        "srate",
        "duration",
        "bad_segments",
        "regularisation",
        "scale_factors",
    }
    missing = sorted(required.difference(contents))
    if missing:
        raise ValueError(
            f"Core checkpoint {path} is missing: {', '.join(missing)}"
        )
    return FREQNESSResult(
        evals=np.asarray(contents["evals"], dtype=float),
        evecs=np.asarray(contents["evecs"], dtype=float),
        pats=np.asarray(contents["pats"], dtype=float),
        ts=np.asarray(contents["ts"], dtype=float),
        frex=np.asarray(contents["frex"], dtype=float).reshape(-1),
        fwhm=np.asarray(contents["fwhm"], dtype=float).reshape(-1),
        srate=_scalar(contents, "srate"),
        duration=_scalar(contents, "duration"),
        bad_segments=np.asarray(
            contents["bad_segments"], dtype=np.int64
        ).reshape(-1),
        regularisation=_scalar(contents, "regularisation"),
        scale_factors=np.asarray(
            contents["scale_factors"], dtype=float
        ).reshape(-1),
    )


def _serializable_dataclass(value: Any) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    for item in dataclass_fields(value):
        if item.name == "figures":
            continue
        field_value = getattr(value, item.name)
        if field_value is None or isinstance(field_value, (str, Path)):
            continue
        if np.isscalar(field_value) or isinstance(field_value, np.ndarray):
            payload[item.name] = field_value
    return payload


def _save_payload(
    output_root: Path,
    condition: str,
    analysis: str,
    payload: Mapping[str, Any],
) -> Path:
    condition_root = output_root / condition
    condition_root.mkdir(parents=True, exist_ok=True)
    path = condition_root / f"{analysis}_python.mat"
    savemat(
        path,
        {"condition": condition, **dict(payload)},
        do_compression=True,
        long_field_names=True,
    )
    return path


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


def FREQNESS_MainPipeline(
    config: FREQNESSPipelineConfig | None = None,
) -> FREQNESSPipelineResult:
    """Run selected MATLAB-equivalent FREQ-NESS stages across dataset folders.

    Analyses are selected using their exact public function names. If a
    secondary analysis is requested without `FREQNESS_NetworkEstimation`, the
    pipeline loads the previously exported Python core checkpoint. This makes
    one-function-at-a-time MATLAB/Python validation deterministic and avoids
    silently recomputing a different prerequisite.
    """
    if config is None:
        config = FREQNESSPipelineConfig()
    selected = _validate_config(config)
    all_data, MNI, root = FREQNESS_Startup(config.toolbox_root)
    names = _condition_names(root, len(all_data))
    output_root = _output_root(config, root)
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

        if NETWORK_ESTIMATION in selected:
            FREQ = FREQNESS_NetworkEstimation(
                data,
                np.asarray(config.frex, dtype=float),
                config.srate,
                **dict(config.network_options),
            )
            condition_result.outputs[NETWORK_ESTIMATION] = FREQ
            if config.save_outputs:
                path = _save_payload(
                    output_root,
                    condition,
                    NETWORK_ESTIMATION,
                    _freq_payload(FREQ, condition),
                )
                condition_result.files[NETWORK_ESTIMATION] = path
        else:
            FREQ = _load_freq_checkpoint(
                _checkpoint_path(output_root, condition)
            )

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
                    "path_output": output_root / condition,
                    "MNI_coords": MNI,
                },
                plot_all=config.plot_all,
                save_nifti=config.save_nifti,
                show=config.show,
            )
            condition_result.outputs[VISUALIZER] = visualization
            if config.save_outputs:
                payload = {
                    "landscape_frequencies": visualization.landscape_frequencies,
                    "pattern_frequencies": visualization.pattern_frequencies,
                    "normalized_patterns": visualization.normalized_patterns,
                }
                if visualization.group_patterns is not None:
                    payload["group_patterns"] = visualization.group_patterns
                condition_result.files[VISUALIZER] = _save_payload(
                    output_root, condition, VISUALIZER, payload
                )

        if ENTROPY_LANDSCAPE in selected:
            H2, ED = FREQNESS_EntropyLandscape(
                FREQ,
                plot_all=config.plot_all,
                show=config.show,
            )
            output = {"H2": H2, "ED": ED}
            condition_result.outputs[ENTROPY_LANDSCAPE] = output
            if config.save_outputs:
                condition_result.files[ENTROPY_LANDSCAPE] = _save_payload(
                    output_root, condition, ENTROPY_LANDSCAPE, output
                )

        if EXPONENTIAL_DK in selected:
            decay_coefficients, good_fit = FREQNESS_ExponentialDK(
                FREQ,
                which_comp=config.expdk_which_comp,
                range2fit=config.expdk_range2fit,
                plot_all=config.plot_all,
                show=config.show,
            )
            output = {
                "decayCoeff": decay_coefficients,
                "goodFit_R2": good_fit.R2,
            }
            condition_result.outputs[EXPONENTIAL_DK] = (
                decay_coefficients,
                good_fit,
            )
            if config.save_outputs:
                condition_result.files[EXPONENTIAL_DK] = _save_payload(
                    output_root, condition, EXPONENTIAL_DK, output
                )

        if FREQ_GRADIENTS in selected and _require_mni(
            MNI, FREQ_GRADIENTS, condition
        ):
            coefficients, good_fit = FREQNESS_FreqGradients(
                FREQ,
                MNI,
                frex2model=config.freqgrad_frex2model,
                comp2model=config.freqgrad_comp2model,
                plot_all=config.plot_all,
                show=config.show,
            )
            condition_result.outputs[FREQ_GRADIENTS] = (
                coefficients,
                good_fit,
            )
            if config.save_outputs:
                payload = {
                    "gradCoeff": coefficients,
                    **{
                        f"goodFit_{name}": value
                        for name, value in _serializable_dataclass(
                            good_fit
                        ).items()
                    },
                }
                condition_result.files[FREQ_GRADIENTS] = _save_payload(
                    output_root, condition, FREQ_GRADIENTS, payload
                )

        if COMP_GRADIENTS in selected and _require_mni(
            MNI, COMP_GRADIENTS, condition
        ):
            coefficients, good_fit = FREQNESS_CompGradients(
                FREQ,
                MNI,
                freq2model=config.compgrad_freq2model,
                comps2model=config.compgrad_comps2model,
                plot_all=config.plot_all,
                show=config.show,
            )
            condition_result.outputs[COMP_GRADIENTS] = (
                coefficients,
                good_fit,
            )
            if config.save_outputs:
                payload = {
                    "gradCoeff": coefficients,
                    **{
                        f"goodFit_{name}": value
                        for name, value in _serializable_dataclass(
                            good_fit
                        ).items()
                    },
                }
                condition_result.files[COMP_GRADIENTS] = _save_payload(
                    output_root, condition, COMP_GRADIENTS, payload
                )

        if CROSS_COUPLING in selected:
            coupling = FREQNESS_CrossCoupling(
                FREQ,
                config.lfo_freq,
                MNI=MNI,
                plot_all=config.plot_all,
                show=config.show,
            )
            condition_result.outputs[CROSS_COUPLING] = coupling
            if config.save_outputs:
                condition_result.files[CROSS_COUPLING] = _save_payload(
                    output_root,
                    condition,
                    CROSS_COUPLING,
                    _serializable_dataclass(coupling),
                )

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
            "Repeat --analysis to run multiple stages."
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
        "--output-directory",
        type=Path,
        help="Override the comparison-output directory.",
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
    parser.add_argument(
        "--no-save",
        action="store_true",
        help="Do not save MATLAB-compatible stage outputs.",
    )
    arguments = parser.parse_args(argv)
    config = FREQNESSPipelineConfig(
        toolbox_root=arguments.toolbox_root,
        analyses=(
            tuple(arguments.analysis)
            if arguments.analysis is not None
            else ALL_ANALYSES
        ),
        output_directory=arguments.output_directory,
        show=not arguments.no_show,
        save_nifti=not arguments.no_nifti,
        save_outputs=not arguments.no_save,
    )
    result = FREQNESS_MainPipeline(config)
    completed = sum(bool(condition.outputs) for condition in result.conditions)
    print(
        f"FREQ-NESS pipeline completed for {completed} non-empty "
        f"condition(s)."
    )
    for condition in result.conditions:
        for analysis, path in condition.files.items():
            print(f"{condition.name} | {analysis} | {path}")
    return 0


__all__ = [
    "ALL_ANALYSES",
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
    "VISUALIZER",
    "main",
]
