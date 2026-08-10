"""Model spatial gradients across FREQ-NESS components.

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

from typing import Any
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._gradient_utils import (
    FREQNESSGradientGoodFit,
    field_value,
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


def _automatic_frequency(FREQ: Any, nfrequencies: int) -> int:
    raw = np.asarray(field_value(FREQ, "evals"))
    if raw.size == 0:
        raise ValueError(
            "FREQ.evals is required when freq2model is not specified"
        )
    if not np.issubdtype(raw.dtype, np.number):
        raise TypeError("FREQ.evals must contain numeric values")
    eigenvalues = np.asarray(raw)
    if eigenvalues.ndim == 2:
        eigenvalues = eigenvalues[..., np.newaxis]
    elif eigenvalues.ndim != 3:
        raise ValueError(
            "FREQ.evals must have shape (eigenvalues, frequencies"
            "[, participants])"
        )
    if eigenvalues.shape[1] != nfrequencies:
        raise ValueError(
            "The frequency dimension of FREQ.evals must match FREQ.pats"
        )
    first_component = np.abs(eigenvalues[0, :, :])
    if not np.any(np.isfinite(first_component)):
        raise ValueError("FREQ.evals contains no finite first-component values")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        average = np.nanmean(first_component, axis=1)
    return int(np.nanargmax(average))


def _frequency_index(
    FREQ: Any,
    frequencies: FloatArray,
    has_frequencies: bool,
    freq2model: Any,
) -> int:
    if freq2model is None:
        warnings.warn(
            "No frequency was defined. Defaulting to the most prominent "
            "first-component frequency in FREQ.evals.",
            UserWarning,
            stacklevel=3,
        )
        return _automatic_frequency(FREQ, frequencies.size)
    if isinstance(freq2model, (bool, np.bool_)) or not np.isscalar(freq2model):
        raise TypeError("freq2model must be a scalar")
    try:
        requested = float(freq2model)
    except (TypeError, ValueError) as error:
        raise TypeError("freq2model must be a scalar") from error
    if not np.isfinite(requested):
        raise ValueError("freq2model must be finite")
    if has_frequencies:
        index = int(np.argmin(np.abs(frequencies - requested)))
        if not np.isclose(requested, frequencies[index], atol=1e-3, rtol=0.0):
            warnings.warn(
                "Requested freq2model is not present in FREQ.frex. Using "
                "the closest available frequency instead.",
                UserWarning,
                stacklevel=3,
            )
        return index
    if not requested.is_integer() or not 1 <= int(requested) <= frequencies.size:
        raise ValueError(
            "Without FREQ.frex, freq2model is a one-based integer index "
            f"between 1 and {frequencies.size}"
        )
    return int(requested) - 1


def _component_indices(
    comps2model: ArrayLike | None,
    ncomponents: int,
) -> NDArray[np.int64]:
    if comps2model is None or np.asarray(comps2model).size == 0:
        return np.arange(ncomponents, dtype=np.int64)
    requested = np.asarray(comps2model)
    if requested.ndim != 1 or requested.size != 2:
        raise ValueError("comps2model must be a two-element vector such as [1, 5]")
    if not np.issubdtype(requested.dtype, np.number):
        raise TypeError("comps2model must contain integer component numbers")
    requested = np.asarray(requested, dtype=float)
    if not np.all(np.isfinite(requested)) or not np.all(
        requested == np.floor(requested)
    ):
        raise ValueError("comps2model must contain integer component numbers")
    first, last = sorted(int(value) for value in requested)
    if first < 1 or last > ncomponents:
        raise ValueError(
            f"comps2model values must be between 1 and {ncomponents}"
        )
    return np.arange(first - 1, last, dtype=np.int64)


def FREQNESS_CompGradients(
    FREQ: Any,
    MNI: ArrayLike,
    *,
    freq2model: float | int | None = None,
    comps2model: ArrayLike | None = None,
    threshold_sd: float = 2.0,
    plot_all: bool = False,
    show: bool = True,
) -> tuple[FloatArray, FREQNESSGradientGoodFit]:
    """Fit spatial gradients across components at one frequency.

    Both linear and quadratic models predict the one-based component index
    from suprathreshold MNI coordinates. BIC selects the reported model and
    coefficients are returned as ``[b0, b1, b2]``. If ``freq2model`` is
    omitted, the peak mean absolute first-component eigenvalue selects the
    frequency, matching the MATLAB workflow.

    ``threshold_sd=2`` preserves the MATLAB selection rule while allowing
    explicit sensitivity analyses through a stable keyword argument.
    """
    plot_all, show = validate_plot_flags(plot_all, show)
    threshold_sd = validate_threshold(threshold_sd)
    all_patterns = pattern_array(FREQ)
    nvoxels, ncomponents, nfrequencies, _ = all_patterns.shape
    coordinates = mni_array(MNI, nvoxels)
    frequencies, has_frequencies = frequency_axis(FREQ, nfrequencies)
    if not has_frequencies:
        warn_missing_frequency_axis()
    frequency = _frequency_index(
        FREQ, frequencies, has_frequencies, freq2model
    )
    selected_components = _component_indices(comps2model, ncomponents)

    patterns = all_patterns[:, :, frequency, :]
    patterns, retained = threshold_patterns(patterns, threshold_sd)
    coefficients, diagnostics = model_subjects(
        patterns,
        coordinates,
        selected_components,
        threshold_sd,
        retained,
    )
    frequency_label = (
        f"{frequencies[frequency]:.4g} Hz"
        if has_frequencies
        else f"frequency index {frequency + 1}"
    )
    diagnostics.figures = plot_gradients(
        patterns,
        coordinates,
        selected_components,
        [f"Comp {component}" for component in range(1, ncomponents + 1)],
        "Component",
        f"Spatial gradient across components at {frequency_label}",
        plot_all=plot_all,
    )
    if show:
        import matplotlib.pyplot as plt

        plt.show()
    return coefficients, diagnostics


component_gradients = FREQNESS_CompGradients


__all__ = [
    "FREQNESS_CompGradients",
    "FREQNESSGradientGoodFit",
    "component_gradients",
]
