"""Shared validation, fitting, and plotting for spatial gradient analyses.

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

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence
import warnings

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(slots=True)
class FREQNESSGradientGoodFit:
    """Diagnostics returned by both FREQ-NESS spatial-gradient functions."""

    R2_linear: FloatArray
    R2_quadratic: FloatArray
    R2_best: FloatArray
    bestOrder: FloatArray
    threshold_sd: float
    retained_voxels: IntArray
    figures: list[Any] = field(default_factory=list)


@dataclass(slots=True)
class _Fit:
    coefficients: FloatArray
    R2_linear: float
    R2_quadratic: float
    R2_best: float
    best_order: float


def field_value(container: Any, name: str, *, required: bool = True) -> Any:
    """Read a field from a mapping or attribute-based FREQ container."""
    if isinstance(container, Mapping):
        if name in container:
            return container[name]
    elif hasattr(container, name):
        return getattr(container, name)
    if required:
        raise ValueError(f"FREQ.{name} is required")
    return None


def pattern_array(FREQ: Any) -> FloatArray:
    """Return patterns as voxels x components x frequencies x participants."""
    raw = np.asarray(field_value(FREQ, "pats"))
    if raw.size == 0:
        raise ValueError("FREQ.pats is empty")
    if not np.issubdtype(raw.dtype, np.number):
        raise TypeError("FREQ.pats must contain numeric values")
    if np.iscomplexobj(raw):
        raise TypeError("FREQ.pats must contain real values")
    patterns = np.asarray(raw, dtype=float)
    if patterns.ndim == 3:
        patterns = patterns[..., np.newaxis]
    elif patterns.ndim != 4:
        raise ValueError(
            "FREQ.pats must have shape (voxels, components, frequencies" 
            "[, participants])"
        )
    if min(patterns.shape) == 0:
        raise ValueError("FREQ.pats is empty")
    if np.any(np.isinf(patterns)):
        raise ValueError("FREQ.pats contains infinite values")
    return patterns


def mni_array(MNI: Any, nvoxels: int) -> FloatArray:
    """Validate MNI coordinates."""
    coordinates = np.asarray(MNI)
    if not np.issubdtype(coordinates.dtype, np.number):
        raise TypeError("MNI must contain numeric values")
    if np.iscomplexobj(coordinates):
        raise TypeError("MNI must contain real values")
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.shape != (nvoxels, 3):
        raise ValueError(
            "MNI must have shape (n_voxels, 3) and match FREQ.pats"
        )
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("MNI must contain only finite coordinates")
    return coordinates


def frequency_axis(
    FREQ: Any,
    nfrequencies: int,
) -> tuple[FloatArray, bool]:
    """Return frequencies or one-based frequency indices."""
    raw = field_value(FREQ, "frex", required=False)
    if raw is None or np.asarray(raw).size == 0:
        return np.arange(1, nfrequencies + 1, dtype=float), False
    frequencies = np.asarray(raw, dtype=float)
    if (
        frequencies.ndim not in (1, 2)
        or (frequencies.ndim == 2 and 1 not in frequencies.shape)
        or frequencies.size != nfrequencies
        or not np.all(np.isfinite(frequencies))
    ):
        raise ValueError(
            "FREQ.frex must be a finite vector matching the frequency "
            "dimension of FREQ.pats"
        )
    return frequencies.reshape(-1).copy(), True


def validate_threshold(threshold_sd: Any) -> float:
    """Validate a non-negative standard-deviation threshold multiplier."""
    if isinstance(threshold_sd, (bool, np.bool_)):
        raise TypeError("threshold_sd must be a non-negative finite number")
    try:
        threshold = float(threshold_sd)
    except (TypeError, ValueError) as error:
        raise TypeError(
            "threshold_sd must be a non-negative finite number"
        ) from error
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("threshold_sd must be a non-negative finite number")
    return threshold


def validate_plot_flags(plot_all: Any, show: Any) -> tuple[bool, bool]:
    """Validate plotting flags without accepting integers as booleans."""
    for value, name in ((plot_all, "plot_all"), (show, "show")):
        if not isinstance(value, (bool, np.bool_)):
            raise TypeError(f"{name} must be a boolean")
    return bool(plot_all), bool(show)


def threshold_patterns(
    patterns: FloatArray,
    threshold_sd: float,
) -> tuple[FloatArray, IntArray]:
    """Normalize and threshold each level independently within participant."""
    thresholded = np.full(patterns.shape, np.nan, dtype=float)
    retained = np.zeros(patterns.shape[1:], dtype=np.int64)
    absolute = np.abs(patterns)

    for subject in range(absolute.shape[2]):
        for level in range(absolute.shape[1]):
            values = absolute[:, level, subject]
            finite = np.isfinite(values)
            if not np.any(finite):
                continue
            maximum = float(np.max(values[finite]))
            if maximum <= 0:
                continue
            normalized = np.full(values.shape, np.nan, dtype=float)
            normalized[finite] = values[finite] / maximum
            finite_values = normalized[finite]
            ddof = 1 if finite_values.size > 1 else 0
            cutoff = float(
                np.mean(finite_values)
                + threshold_sd * np.std(finite_values, ddof=ddof)
            )
            keep = finite & (normalized >= cutoff)
            thresholded[keep, level, subject] = normalized[keep]
            retained[level, subject] = int(np.count_nonzero(keep))
    return thresholded, retained


def _least_squares_fit(
    x: FloatArray,
    y: FloatArray,
    order: int,
) -> tuple[FloatArray, float, float] | None:
    design = np.column_stack([x**power for power in range(order + 1)])
    if np.linalg.matrix_rank(design) < order + 1:
        return None
    coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
    predicted = design @ coefficients
    residual = y - predicted
    sse = float(residual @ residual)
    sst = float(np.sum((y - np.mean(y)) ** 2))
    if not np.isfinite(sse) or sst <= 0:
        return None
    r_squared = 1.0 - sse / sst
    return np.asarray(coefficients, dtype=float), sse, float(r_squared)


def fit_gradient(x: FloatArray, y: FloatArray) -> _Fit:
    """Fit linear/quadratic models and choose order by BIC."""
    empty = _Fit(np.full(3, np.nan), np.nan, np.nan, np.nan, np.nan)
    finite = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x[finite], dtype=float)
    y = np.asarray(y[finite], dtype=float)
    if x.size < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return empty

    linear = _least_squares_fit(x, y, 1)
    if linear is None:
        return empty
    linear_coefficients, linear_sse, linear_r2 = linear

    quadratic = None
    if x.size >= 3 and np.unique(x).size >= 3:
        quadratic = _least_squares_fit(x, y, 2)
    quadratic_r2 = np.nan if quadratic is None else quadratic[2]

    coefficients = np.array(
        [linear_coefficients[0], linear_coefficients[1], 0.0], dtype=float
    )
    best_order = 1.0
    best_r2 = linear_r2
    if quadratic is not None:
        quadratic_coefficients, quadratic_sse, quadratic_r2 = quadratic
        scale = max(float(np.sum((y - np.mean(y)) ** 2)), 1.0)
        zero_tolerance = np.finfo(float).eps * scale * max(x.size, 1) * 32
        linear_perfect = linear_sse <= zero_tolerance
        quadratic_perfect = quadratic_sse <= zero_tolerance
        if quadratic_perfect and not linear_perfect:
            quadratic_wins = True
        elif linear_perfect:
            quadratic_wins = False
        else:
            bic_linear = x.size * np.log(linear_sse / x.size) + 2 * np.log(x.size)
            bic_quadratic = (
                x.size * np.log(quadratic_sse / x.size) + 3 * np.log(x.size)
            )
            quadratic_wins = bool(bic_quadratic < bic_linear)
        if quadratic_wins:
            coefficients = np.asarray(quadratic_coefficients, dtype=float)
            best_order = 2.0
            best_r2 = quadratic_r2

    return _Fit(
        coefficients,
        float(linear_r2),
        float(quadratic_r2),
        float(best_r2),
        best_order,
    )


def collect_level_data(
    patterns: FloatArray,
    coordinates: FloatArray,
    level_indices: IntArray,
    dimension: int,
) -> tuple[FloatArray, FloatArray]:
    """Collect unweighted coordinate/one-based-level pairs."""
    x_parts: list[FloatArray] = []
    y_parts: list[FloatArray] = []
    for level in level_indices:
        valid = np.isfinite(patterns[:, level])
        if np.any(valid):
            x_parts.append(coordinates[valid, dimension])
            y_parts.append(
                np.full(np.count_nonzero(valid), level + 1, dtype=float)
            )
    if not x_parts:
        return np.empty(0, dtype=float), np.empty(0, dtype=float)
    return np.concatenate(x_parts), np.concatenate(y_parts)


def model_subjects(
    patterns: FloatArray,
    MNI: FloatArray,
    level_indices: IntArray,
    threshold_sd: float,
    retained_voxels: IntArray,
) -> tuple[FloatArray, FREQNESSGradientGoodFit]:
    """Fit X/Y/Z gradients independently for every participant."""
    nsubjects = patterns.shape[2]
    coefficients = np.full((3, 3, nsubjects), np.nan, dtype=float)
    linear = np.full((3, nsubjects), np.nan, dtype=float)
    quadratic = np.full((3, nsubjects), np.nan, dtype=float)
    best = np.full((3, nsubjects), np.nan, dtype=float)
    order = np.full((3, nsubjects), np.nan, dtype=float)

    for subject in range(nsubjects):
        for dimension in range(3):
            x, y = collect_level_data(
                patterns[:, :, subject], MNI, level_indices, dimension
            )
            fit = fit_gradient(x, y)
            coefficients[dimension, :, subject] = fit.coefficients
            linear[dimension, subject] = fit.R2_linear
            quadratic[dimension, subject] = fit.R2_quadratic
            best[dimension, subject] = fit.R2_best
            order[dimension, subject] = fit.best_order

    diagnostics = FREQNESSGradientGoodFit(
        R2_linear=linear,
        R2_quadratic=quadratic,
        R2_best=best,
        bestOrder=order,
        threshold_sd=threshold_sd,
        retained_voxels=retained_voxels,
    )
    return coefficients, diagnostics


def group_fits(
    patterns: FloatArray,
    MNI: FloatArray,
    level_indices: IntArray,
) -> list[_Fit]:
    """Fit aggregated participant data for visualization only."""
    group_patterns = np.concatenate(
        [patterns[:, :, subject] for subject in range(patterns.shape[2])],
        axis=0,
    )
    group_coordinates = np.tile(MNI, (patterns.shape[2], 1))
    fits: list[_Fit] = []
    for dimension in range(3):
        x, y = collect_level_data(
            group_patterns, group_coordinates, level_indices, dimension
        )
        fits.append(fit_gradient(x, y))
    return fits


def _evaluate(coefficients: FloatArray, x: FloatArray) -> FloatArray:
    return coefficients[0] + coefficients[1] * x + coefficients[2] * x**2


def plot_gradients(
    patterns: FloatArray,
    MNI: FloatArray,
    level_indices: IntArray,
    level_labels: Sequence[str],
    ylabel: str,
    title: str,
    *,
    plot_all: bool,
) -> list[Any]:
    """Plot the group map and optional participant-specific maps."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "Spatial-gradient plotting requires matplotlib. "
            "Install freqness[plot]."
        ) from error

    figures: list[Any] = []
    group_patterns = np.concatenate(
        [patterns[:, :, subject] for subject in range(patterns.shape[2])],
        axis=0,
    )
    group_coordinates = np.tile(MNI, (patterns.shape[2], 1))
    figures.append(
        _plot_one(
            group_patterns,
            group_coordinates,
            level_indices,
            level_labels,
            ylabel,
            f"{title} - Group level",
            group_fits(patterns, MNI, level_indices),
            jitter=True,
        )
    )
    if plot_all:
        for subject in range(patterns.shape[2]):
            fits = []
            for dimension in range(3):
                x, y = collect_level_data(
                    patterns[:, :, subject], MNI, level_indices, dimension
                )
                fits.append(fit_gradient(x, y))
            figures.append(
                _plot_one(
                    patterns[:, :, subject],
                    MNI,
                    level_indices,
                    level_labels,
                    ylabel,
                    f"{title} - Participant #{subject + 1}",
                    fits,
                    jitter=False,
                )
            )
    return figures


def _plot_one(
    patterns: FloatArray,
    coordinates: FloatArray,
    level_indices: IntArray,
    level_labels: Sequence[str],
    ylabel: str,
    title: str,
    fits: Sequence[_Fit],
    *,
    jitter: bool,
) -> Any:
    import matplotlib.pyplot as plt

    nlevels = patterns.shape[1]
    colors = plt.get_cmap("viridis")(
        np.linspace(0.08, 0.92, max(nlevels, 2))[:nlevels]
    )
    figure, axes = plt.subplots(
        1, 3, figsize=(16, 6), sharey=True, constrained_layout=True
    )
    rng = np.random.default_rng(0)
    coordinate_labels = ("X-coordinate", "Y-coordinate", "Z-coordinate")
    selected_min = int(np.min(level_indices)) + 1
    selected_max = int(np.max(level_indices)) + 1

    for dimension, axis in enumerate(axes):
        for level in range(nlevels):
            amplitudes = patterns[:, level]
            valid = np.isfinite(amplitudes)
            if not np.any(valid):
                continue
            x = coordinates[valid, dimension].copy()
            if jitter:
                x += rng.normal(0.0, 0.15, size=x.size)
            y = np.full(x.size, level + 1, dtype=float)
            sizes = 12.0 + 55.0 * amplitudes[valid]
            axis.scatter(
                x,
                y,
                s=sizes,
                color=colors[level],
                alpha=0.7,
                linewidths=0.2,
                edgecolors="black",
            )

        fit = fits[dimension]
        if np.all(np.isfinite(fit.coefficients)):
            x_fit, _ = collect_level_data(
                patterns, coordinates, level_indices, dimension
            )
            if x_fit.size:
                x_line = np.linspace(float(np.min(x_fit)), float(np.max(x_fit)), 300)
                y_line = _evaluate(fit.coefficients, x_line)
                visible = (y_line >= selected_min) & (y_line <= selected_max)
                y_line = np.where(visible, y_line, np.nan)
                axis.plot(x_line, y_line, color="black", linewidth=2)
        axis.axhline(selected_min, color="black", linestyle="--", linewidth=1)
        if selected_max != selected_min:
            axis.axhline(selected_max, color="black", linestyle="--", linewidth=1)
        axis.set_xlabel(coordinate_labels[dimension])
        axis.set_title(f"{coordinate_labels[dimension][0]} gradient")
        axis.grid(True, alpha=0.18)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        finite_coordinates = coordinates[:, dimension][
            np.isfinite(coordinates[:, dimension])
        ]
        if finite_coordinates.size:
            extent = float(np.ptp(finite_coordinates))
            padding = max(0.04 * extent, 0.5)
            axis.set_xlim(
                float(np.min(finite_coordinates)) - padding,
                float(np.max(finite_coordinates)) + padding,
            )

    axes[0].set_ylabel(ylabel)
    axes[0].set_yticks(np.arange(1, nlevels + 1))
    axes[0].set_yticklabels(level_labels)
    axes[0].set_ylim(0.5, nlevels + 0.5)
    figure.suptitle(title, fontsize=16, fontweight="bold")
    return figure


def warn_missing_frequency_axis() -> None:
    warnings.warn(
        "FREQ.frex is missing or empty. One-based frequency indices will "
        "be used instead of Hz.",
        UserWarning,
        stacklevel=3,
    )


__all__ = ["FREQNESSGradientGoodFit"]
