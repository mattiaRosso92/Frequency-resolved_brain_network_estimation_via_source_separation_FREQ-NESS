"""Model spatial gradients across FREQ-NESS frequencies."""

from __future__ import annotations

from typing import Any
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._gradient_utils import (
    FREQNESSGradientGoodFit,
    frequency_axis,
    mni_array,
    model_subjects,
    pattern_array,
    plot_gradients,
    threshold_patterns,
    validate_plot_flags,
    validate_threshold,
    warn_missing_frequency_axis,
)


FloatArray = NDArray[np.float64]


def _component_index(comp2model: Any, ncomponents: int) -> int:
    if comp2model is None:
        return 0
    valid = (
        isinstance(comp2model, (int, np.integer))
        and not isinstance(comp2model, (bool, np.bool_))
        and 1 <= int(comp2model) <= ncomponents
    )
    if not valid:
        warnings.warn(
            f"comp2model must be an integer between 1 and {ncomponents}. "
            "Defaulting to component 1.",
            UserWarning,
            stacklevel=3,
        )
        return 0
    return int(comp2model) - 1


def _frequency_indices(
    frequencies: FloatArray,
    has_frequencies: bool,
    frex2model: ArrayLike | None,
) -> NDArray[np.int64]:
    if frex2model is None or np.asarray(frex2model).size == 0:
        return np.arange(frequencies.size, dtype=np.int64)
    if not has_frequencies:
        raise ValueError("frex2model was provided in Hz, but FREQ.frex is missing")
    requested = np.asarray(frex2model, dtype=float)
    if requested.ndim != 1 or requested.size != 2:
        raise ValueError("frex2model must be a two-element vector such as [8, 12]")
    if not np.all(np.isfinite(requested)):
        raise ValueError("frex2model must contain finite values")
    requested = np.sort(requested)
    first = int(np.argmin(np.abs(frequencies - requested[0])))
    last = int(np.argmin(np.abs(frequencies - requested[1])))
    first, last = sorted((first, last))
    matched = np.array([frequencies[first], frequencies[last]])
    if not np.all(np.isclose(requested, matched, atol=1e-3, rtol=0.0)):
        warnings.warn(
            "Some requested frequencies in frex2model are not present in "
            "FREQ.frex. Using closest matches instead.",
            UserWarning,
            stacklevel=3,
        )
    return np.arange(first, last + 1, dtype=np.int64)


def FREQNESS_FreqGradients(
    FREQ: Any,
    MNI: ArrayLike,
    *,
    frex2model: ArrayLike | None = None,
    comp2model: int | None = None,
    threshold_sd: float = 1.0,
    plot_all: bool = False,
    show: bool = True,
) -> tuple[FloatArray, FREQNESSGradientGoodFit]:
    """Fit spatial gradients across frequency for one component.

    Spatial patterns are normalized and thresholded independently for every
    frequency and participant. For each MNI dimension, both linear and
    quadratic models predict the one-based frequency-bin index from the
    coordinates of suprathreshold voxels; BIC selects the reported model.
    Coefficients are returned in the explicit order ``[b0, b1, b2]``.

    ``threshold_sd=1`` preserves the MATLAB selection rule. It is exposed as
    a keyword so future sensitivity analyses can change it without modifying
    the implementation.
    """
    plot_all, show = validate_plot_flags(plot_all, show)
    threshold_sd = validate_threshold(threshold_sd)
    all_patterns = pattern_array(FREQ)
    nvoxels, ncomponents, nfrequencies, _ = all_patterns.shape
    coordinates = mni_array(MNI, nvoxels)
    frequencies, has_frequencies = frequency_axis(FREQ, nfrequencies)
    if not has_frequencies:
        warn_missing_frequency_axis()
    component = _component_index(comp2model, ncomponents)
    selected_frequencies = _frequency_indices(
        frequencies, has_frequencies, frex2model
    )

    patterns = all_patterns[:, component, :, :]
    patterns, retained = threshold_patterns(patterns, threshold_sd)
    coefficients, diagnostics = model_subjects(
        patterns,
        coordinates,
        selected_frequencies,
        threshold_sd,
        retained,
    )
    labels = (
        [f"{frequency:.4g} Hz" for frequency in frequencies]
        if has_frequencies
        else [f"Index {index}" for index in range(1, nfrequencies + 1)]
    )
    diagnostics.figures = plot_gradients(
        patterns,
        coordinates,
        selected_frequencies,
        labels,
        "Frequency",
        f"Spatial gradient across frequencies - Component #{component + 1}",
        plot_all=plot_all,
    )
    if show:
        import matplotlib.pyplot as plt

        plt.show()
    return coefficients, diagnostics


frequency_gradients = FREQNESS_FreqGradients


__all__ = [
    "FREQNESS_FreqGradients",
    "FREQNESSGradientGoodFit",
    "frequency_gradients",
]
