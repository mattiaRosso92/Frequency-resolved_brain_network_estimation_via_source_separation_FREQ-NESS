"""Remove selected FREQ-NESS network contributions from broadband data.

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

from .FREQNESS_BackProjection import (
    FREQNESS_BackProjection,
    _canonical_inputs,
    _component_indices,
    _frequency_value,
)


FloatArray = NDArray[np.float64]


def _canonical_data(data: ArrayLike) -> tuple[FloatArray, bool]:
    values = np.asarray(data)
    if values.size == 0:
        raise ValueError("data must be a non-empty numeric array")
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise TypeError("data must contain real numeric values")
    if values.ndim not in (2, 3):
        raise ValueError(
            "data must have shape (voxels, time) or "
            "(voxels, time, participants)"
        )
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("data must contain finite values")
    was_two_dimensional = values.ndim == 2
    if was_two_dimensional:
        values = values[:, :, np.newaxis]
    return values, was_two_dimensional


def FREQNESS_NetworkRemoval(
    FREQ: Any,
    data: ArrayLike,
    freq2remove: float,
    *,
    comps2remove: ArrayLike | None = None,
) -> tuple[FloatArray, FloatArray]:
    """Remove selected broadband network activity from the original data.

    ``data`` must be the original broadband array, with the same voxel,
    analyzed-time, and participant ordering used for ``FREQ``. Component
    indices remain MATLAB-compatible and one-based. The returned tuple is
    ``(dataClean, removedActivity)``, and both arrays preserve the dimensionality
    of ``data``.

    BackProjection restores internally rescaled participants to the original
    input units before subtraction. Consequently, the outputs satisfy
    ``data == dataClean + removedActivity`` within floating-point precision.
    """
    requested_frequency = _frequency_value(freq2remove, "freq2remove")
    frequencies, eigenvectors, time_series, _ = _canonical_inputs(FREQ)
    component_numbers, _ = _component_indices(
        comps2remove,
        eigenvectors.shape[1],
        "comps2remove",
    )
    broadband, was_two_dimensional = _canonical_data(data)

    expected_shape = (
        eigenvectors.shape[0],
        time_series.shape[1],
        eigenvectors.shape[3],
    )
    if broadband.shape[0] != expected_shape[0]:
        raise ValueError(
            "The number of voxels in data does not match FREQ.evecs"
        )
    if broadband.shape[1] != expected_shape[1]:
        raise ValueError(
            "The number of timepoints in data does not match FREQ.ts; "
            "provide the same analyzed data segment used by "
            "FREQNESS_NetworkEstimation"
        )
    if broadband.shape[2] != expected_shape[2]:
        raise ValueError(
            "The number of participants in data does not match the FREQ structure"
        )

    frequency_index = int(np.argmin(np.abs(frequencies - requested_frequency)))
    actual_frequency = float(frequencies[frequency_index])
    if not np.isclose(
        requested_frequency,
        actual_frequency,
        atol=1e-10,
        rtol=0.0,
    ):
        warnings.warn(
            f"Requested freq2remove ({requested_frequency:.3f} Hz) is not "
            f"present in FREQ.frex; using the closest available frequency: "
            f"{actual_frequency:.3f} Hz.",
            UserWarning,
            stacklevel=2,
        )

    removed_activity = FREQNESS_BackProjection(
        FREQ,
        actual_frequency,
        comps2project=component_numbers,
    )
    if removed_activity.ndim == 2:
        removed_activity = removed_activity[:, :, np.newaxis]

    data_clean = broadband - removed_activity
    reconstruction_error = float(
        np.max(np.abs(broadband - data_clean - removed_activity))
    )
    tolerance = (
        10.0
        * np.finfo(float).eps
        * max(1.0, float(np.max(np.abs(broadband))))
    )
    if reconstruction_error > tolerance:
        warnings.warn(
            "The subtraction identity exceeded the expected numerical "
            f"tolerance (maximum error: {reconstruction_error:.3g}).",
            RuntimeWarning,
            stacklevel=2,
        )

    if was_two_dimensional:
        return data_clean[:, :, 0], removed_activity[:, :, 0]
    return data_clean, removed_activity


__all__ = ["FREQNESS_NetworkRemoval"]
