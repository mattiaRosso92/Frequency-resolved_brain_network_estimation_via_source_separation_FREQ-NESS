"""Backproject selected FREQ-NESS components into voxel space.

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

from typing import Any, Mapping
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import pinv


FloatArray = NDArray[np.float64]


def _field(container: Any, name: str, *, required: bool = True) -> Any:
    if isinstance(container, Mapping):
        if name in container:
            return container[name]
    elif hasattr(container, name):
        return getattr(container, name)
    if required:
        raise ValueError(f"FREQ.{name} is missing or empty")
    return None


def _real_finite_array(value: Any, name: str) -> FloatArray:
    array = np.asarray(value)
    if array.size == 0:
        raise ValueError(f"{name} is missing or empty")
    if not np.issubdtype(array.dtype, np.number) or np.iscomplexobj(array):
        raise TypeError(f"{name} must contain real numeric values")
    array = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values")
    return array


def _frequency_value(value_to_check: Any, name: str = "freq2project") -> float:
    value = np.asarray(value_to_check)
    if (
        value.ndim != 0
        or not np.issubdtype(value.dtype, np.number)
        or np.iscomplexobj(value)
        or isinstance(value_to_check, (bool, np.bool_))
    ):
        raise TypeError(f"{name} must be one positive finite scalar frequency")
    frequency = float(value)
    if not np.isfinite(frequency) or frequency <= 0:
        raise ValueError(f"{name} must be one positive finite scalar frequency")
    return frequency


def _component_indices(
    comps2project: ArrayLike | None,
    ncomponents: int,
    name: str = "comps2project",
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    if comps2project is None:
        one_based = np.array([1], dtype=np.int64)
    else:
        raw = np.asarray(comps2project)
        if raw.size == 0:
            one_based = np.array([1], dtype=np.int64)
        else:
            if raw.ndim == 0:
                raw = raw.reshape(1)
            if (
                raw.ndim != 1
                or not np.issubdtype(raw.dtype, np.number)
                or np.iscomplexobj(raw)
                or np.issubdtype(raw.dtype, np.bool_)
                or not np.all(np.isfinite(raw))
                or not np.all(raw == np.round(raw))
            ):
                raise ValueError(
                    f"{name} must contain finite integer component indices"
                )

            # MATLAB uses unique(..., 'stable'). Preserve the public order while
            # retaining MATLAB-compatible one-based component numbering.
            stable_unique = list(dict.fromkeys(raw.astype(np.int64).tolist()))
            one_based = np.asarray(stable_unique, dtype=np.int64)

    if np.any(one_based < 1) or np.any(one_based > ncomponents):
        raise ValueError(
            f"{name} must contain indices between 1 and {ncomponents}"
        )
    return one_based, one_based - 1


def _component_text(component_numbers: NDArray[np.int64]) -> str:
    if component_numbers.size == 1:
        return str(int(component_numbers[0]))
    return "[" + " ".join(str(int(value)) for value in component_numbers) + "]"


def _canonical_inputs(
    FREQ: Any,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    frequencies = _real_finite_array(_field(FREQ, "frex"), "FREQ.frex")
    if frequencies.ndim != 1:
        raise ValueError("FREQ.frex must be a one-dimensional vector")
    if np.any(frequencies <= 0):
        raise ValueError("FREQ.frex must contain positive frequencies")

    eigenvectors = _real_finite_array(_field(FREQ, "evecs"), "FREQ.evecs")
    time_series = _real_finite_array(_field(FREQ, "ts"), "FREQ.ts")

    if eigenvectors.ndim == 3:
        eigenvectors = eigenvectors[:, :, :, np.newaxis]
    if time_series.ndim == 3:
        time_series = time_series[:, :, :, np.newaxis]
    if eigenvectors.ndim != 4:
        raise ValueError(
            "FREQ.evecs must have shape "
            "(voxels, components, frequencies, participants)"
        )
    if time_series.ndim != 4:
        raise ValueError(
            "FREQ.ts must have shape "
            "(components, time, frequencies, participants)"
        )

    nvoxels, ncomponents, nfrequencies, nsubjects = eigenvectors.shape
    if min(nvoxels, ncomponents, nfrequencies, nsubjects) == 0:
        raise ValueError("FREQ.evecs dimensions must be non-empty")
    if frequencies.size != nfrequencies:
        raise ValueError(
            "Length of FREQ.frex does not match the frequency dimension "
            "of FREQ.evecs"
        )
    if (
        time_series.shape[0] != ncomponents
        or time_series.shape[1] == 0
        or time_series.shape[2] != nfrequencies
        or time_series.shape[3] != nsubjects
    ):
        raise ValueError("FREQ.evecs and FREQ.ts have incompatible dimensions")

    raw_scales = _field(FREQ, "scale_factors", required=False)
    if raw_scales is None:
        scale_factors = np.ones(nsubjects, dtype=float)
    else:
        scale_factors = _real_finite_array(
            raw_scales,
            "FREQ.scale_factors",
        ).reshape(-1)
        if scale_factors.size != nsubjects:
            raise ValueError(
                "FREQ.scale_factors must contain one value per participant"
            )
        if np.any(scale_factors <= 0):
            raise ValueError("FREQ.scale_factors must contain positive values")

    return frequencies, eigenvectors, time_series, scale_factors


def FREQNESS_BackProjection(
    FREQ: Any,
    freq2project: float,
    *,
    comps2project: ArrayLike | None = None,
) -> FloatArray:
    """Reconstruct selected FREQ-NESS networks in voxel space.

    The reconstruction uses the dual forward model paired with the stored
    broadband component time series::

        Y = W.T @ X
        A = pinv(W.T)
        X_network = A[:, components] @ Y[components, :]

    ``comps2project`` deliberately uses MATLAB-compatible one-based component
    numbering. When the requested frequency is absent, the closest analyzed
    frequency is used. If network estimation rescaled a participant internally,
    its stored ``scale_factors`` value is removed so the result is returned in
    the units of the original input data.

    A single participant produces ``(voxels, time)`` output; multiple
    participants produce ``(voxels, time, participants)`` output. When fewer
    components than voxels were retained, selecting all retained components
    reconstructs the data's projection onto that retained GED subspace.
    """
    requested_frequency = _frequency_value(freq2project)
    if comps2project is None or np.asarray(comps2project).size == 0:
        print(
            "\nFREQNESS BackProjection: no components specified. Defaulting to component #1."
        )
    frequencies, eigenvectors, time_series, scale_factors = _canonical_inputs(
        FREQ
    )
    component_numbers, component_indices = _component_indices(
        comps2project,
        eigenvectors.shape[1],
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
            f"Requested freq2project ({requested_frequency:.3f} Hz) is not "
            f"present in FREQ.frex; using the closest available frequency: "
            f"{actual_frequency:.3f} Hz.",
            UserWarning,
            stacklevel=2,
        )

    nvoxels = eigenvectors.shape[0]
    ntime = time_series.shape[1]
    nsubjects = eigenvectors.shape[3]
    print(
        "\nFREQNESS BackProjection: backprojecting component(s) "
        f"{_component_text(component_numbers)} at {actual_frequency:.3f} Hz "
        f"for {nsubjects} participants."
    )
    back_projection = np.empty((nvoxels, ntime, nsubjects), dtype=float)

    for subject in range(nsubjects):
        filters = eigenvectors[:, :, frequency_index, subject]
        forward_model = pinv(filters.T, check_finite=False)
        components = time_series[:, :, frequency_index, subject]
        reconstruction = (
            forward_model[:, component_indices]
            @ components[component_indices, :]
        )
        back_projection[:, :, subject] = reconstruction / scale_factors[subject]

    if nsubjects == 1:
        return back_projection[:, :, 0]
    return back_projection


__all__ = ["FREQNESS_BackProjection"]
