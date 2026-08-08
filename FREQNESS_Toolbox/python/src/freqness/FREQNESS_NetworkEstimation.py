"""Frequency-resolved network estimation via generalized eigendecomposition.

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

from dataclasses import dataclass
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import eigh

from .filterFGx import filterFGx


FloatArray = NDArray[np.float64]


@dataclass(slots=True)
class FREQNESSResult:
    """Outputs of :func:`FREQNESS_NetworkEstimation`.

    Array dimensions intentionally mirror the MATLAB implementation. The final
    subject dimension is retained for single-participant inputs.
    """

    evals: FloatArray
    evecs: FloatArray
    pats: FloatArray
    ts: FloatArray
    frex: FloatArray
    fwhm: FloatArray
    srate: float
    duration: float
    bad_segments: NDArray[np.int64]
    regularisation: float
    scale_factors: FloatArray


def _as_subject_array(data: ArrayLike) -> FloatArray:
    values = np.asarray(data)
    if values.ndim == 2:
        values = values[:, :, np.newaxis]
    if values.ndim != 3 or min(values.shape) == 0:
        raise ValueError(
            "data must have shape (variables, time) or "
            "(variables, time, participants)"
        )
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise TypeError("data must contain real numeric values")
    if not np.all(np.isfinite(values)):
        raise ValueError("data contains NaN or infinite values")
    return np.asarray(values, dtype=float).copy()


def _validate_frequencies(frex: ArrayLike, srate: float) -> FloatArray:
    frequencies = np.asarray(frex, dtype=float)
    if frequencies.ndim != 1 or frequencies.size == 0:
        raise ValueError("frex must be a non-empty one-dimensional vector")
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0):
        raise ValueError("frex must contain positive finite values")
    if np.any(frequencies >= srate / 2):
        raise ValueError("all frequencies must be lower than the Nyquist frequency")
    return np.sort(frequencies)


def _filter_widths(
    frequencies: FloatArray,
    fwidth: ArrayLike | None,
    filter_type: str,
) -> FloatArray:
    n_frequencies = frequencies.size
    if fwidth is not None:
        widths = np.asarray(fwidth, dtype=float)
        if widths.ndim == 0:
            widths = np.repeat(widths.item(), n_frequencies)
        elif widths.ndim == 1 and widths.size == n_frequencies:
            widths = widths.copy()
        else:
            raise ValueError(
                "fwidth must be a scalar or contain one value per frequency"
            )
        if not np.all(np.isfinite(widths)) or np.any(widths <= 0):
            raise ValueError("fwidth must contain positive finite values")
        return widths

    reference_frequency = 2.439
    reference_fwhm = 0.35
    n_above = 80
    n_below = 6
    frequencies_above = np.linspace(
        reference_frequency, reference_frequency * (n_above / 2), n_above
    )
    frequencies_below = np.linspace(
        reference_frequency / 2,
        reference_frequency / (2 * n_below),
        n_below,
    )
    reference_frequencies = np.concatenate(
        [frequencies_below[::-1], frequencies_above]
    )
    fwhm_above = np.logspace(
        np.log10(reference_fwhm),
        np.log10(reference_fwhm * n_above),
        n_above,
    )
    fwhm_below = np.logspace(
        np.log10(reference_fwhm),
        np.log10(reference_fwhm / n_below),
        n_below + 1,
    )
    reference_widths = np.concatenate([fwhm_below[:0:-1], fwhm_above])
    lowest_index = int(
        np.argmin(np.abs(reference_frequencies - frequencies[0]))
    )
    lowest_width = reference_widths[lowest_index]

    if filter_type.lower() == "logarithmic":
        return np.logspace(
            np.log10(lowest_width),
            np.log10(lowest_width * n_frequencies),
            n_frequencies,
        )
    if filter_type.lower() == "linear":
        return np.linspace(
            lowest_width,
            lowest_width * n_frequencies,
            n_frequencies,
        )
    raise ValueError("filter must be 'logarithmic' or 'linear'")


def _duration_samples(duration: float | None, n_time: int, srate: float) -> int:
    if duration is None:
        return n_time
    if not np.isscalar(duration) or not np.isfinite(duration) or duration <= 0:
        raise ValueError("duration must be a positive finite scalar in seconds")
    requested = float(duration) * srate
    rounded = int(round(requested))
    if not np.isclose(requested, rounded):
        raise ValueError("duration * srate must define an integer number of samples")
    if rounded > n_time:
        raise ValueError("requested duration exceeds the length of the data")
    if rounded < 2:
        raise ValueError("the requested duration contains fewer than two samples")
    return rounded


def _good_sample_indices(
    bad_segments: ArrayLike | None,
    n_time: int,
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    if bad_segments is None:
        bad_one_based = np.empty(0, dtype=np.int64)
    else:
        raw = np.asarray(bad_segments)
        if raw.ndim != 1 or not np.issubdtype(raw.dtype, np.number):
            raise ValueError("bad_segments must be a vector of integer sample indices")
        if not np.all(np.isfinite(raw)) or not np.all(raw == np.round(raw)):
            raise ValueError("bad_segments must contain finite integer sample indices")
        bad_one_based = np.unique(raw.astype(np.int64))
        if np.any(bad_one_based < 1) or np.any(bad_one_based > n_time):
            raise ValueError(
                "one or more bad_segments indices fall outside the analyzed time range"
            )

    mask = np.ones(n_time, dtype=bool)
    mask[bad_one_based - 1] = False
    good = np.flatnonzero(mask)
    if good.size < 2:
        raise ValueError("fewer than two samples remain for covariance estimation")
    return good, bad_one_based


def _prepare_scale_factors(data: FloatArray, rescale: bool) -> FloatArray:
    factors = np.ones(data.shape[2], dtype=float)
    for subject in range(data.shape[2]):
        scale_reference = float(np.median(np.abs(data[:, :, subject])))
        if scale_reference == 0:
            raise ValueError(
                f"participant {subject + 1} has a median absolute data value of zero"
            )
        order = int(np.floor(np.log10(scale_reference)))
        if order < -2:
            warnings.warn(
                f"Participant {subject + 1} data have a low scale "
                f"(order of magnitude 10e{order}); covariance regularization "
                "may affect network estimation.",
                UserWarning,
                stacklevel=3,
            )
            if rescale:
                factors[subject] = 10.0 ** abs(order - 1)
                data[:, :, subject] *= factors[subject]
    return factors


def FREQNESS_NetworkEstimation(
    data: ArrayLike,
    frex: ArrayLike,
    srate: float,
    *,
    duration: float | None = None,
    fwidth: ArrayLike | None = None,
    filter: str = "logarithmic",
    regularisation: float = 0.01,
    ncomps: int = 30,
    bad_segments: ArrayLike | None = None,
    rescale: bool = False,
) -> FREQNESSResult:
    """Estimate frequency-resolved networks using generalized eigendecomposition.

    Parameters follow the MATLAB function of the same name. ``bad_segments``
    deliberately uses MATLAB-compatible one-based sample indices. Low-amplitude
    data are never rescaled interactively; pass ``rescale=True`` to apply the
    MATLAB function's automatic scale factor after its warning condition.
    """
    if not np.isscalar(srate) or not np.isfinite(srate) or srate <= 0:
        raise ValueError("srate must be a positive finite scalar")
    if (
        not np.isscalar(regularisation)
        or not np.isfinite(regularisation)
        or regularisation < 0
        or regularisation >= 1
    ):
        raise ValueError("regularisation must be in the interval [0, 1)")
    if not isinstance(ncomps, (int, np.integer)) or ncomps < 1:
        raise ValueError("ncomps must be a positive integer")
    if not isinstance(rescale, (bool, np.bool_)):
        raise TypeError("rescale must be a boolean")

    subject_data = _as_subject_array(data)
    n_variables, n_time, n_subjects = subject_data.shape
    if ncomps > n_variables:
        raise ValueError("ncomps cannot exceed the number of input variables")

    sampling_rate = float(srate)
    frequencies = _validate_frequencies(frex, sampling_rate)
    widths = _filter_widths(frequencies, fwidth, filter)
    n_samples = _duration_samples(duration, n_time, sampling_rate)
    subject_data = subject_data[:, :n_samples, :]
    good_samples, bad_one_based = _good_sample_indices(bad_segments, n_samples)
    scale_factors = _prepare_scale_factors(subject_data, bool(rescale))

    n_frequencies = frequencies.size
    eigenvalues = np.zeros((ncomps, n_frequencies, n_subjects), dtype=float)
    eigenvectors = np.zeros(
        (n_variables, ncomps, n_frequencies, n_subjects), dtype=float
    )
    patterns = np.zeros_like(eigenvectors)
    time_series = np.zeros(
        (ncomps, n_samples, n_frequencies, n_subjects), dtype=float
    )

    for subject in range(n_subjects):
        broadband = subject_data[:, :, subject]
        covariance_r = np.atleast_2d(
            np.cov(broadband[:, good_samples], rowvar=True, ddof=1)
        )
        covariance_r_eigenvalues = np.linalg.eigvalsh(covariance_r)
        covariance_r = (
            (1.0 - regularisation) * covariance_r
            + regularisation
            * np.mean(covariance_r_eigenvalues)
            * np.eye(n_variables)
        )

        for frequency_index, (frequency, width) in enumerate(
            zip(frequencies, widths, strict=True)
        ):
            narrowband, _, _ = filterFGx(
                broadband,
                sampling_rate,
                float(frequency),
                float(width),
            )
            covariance_s = np.atleast_2d(
                np.cov(narrowband[:, good_samples], rowvar=True, ddof=1)
            )
            covariance_s = covariance_s + 1e-6 * np.eye(n_variables)

            all_values, all_vectors = eigh(
                covariance_s,
                covariance_r,
                check_finite=False,
            )
            order = np.argsort(all_values)[::-1]
            all_values = np.real(all_values[order])
            all_vectors = np.real(all_vectors[:, order])
            value_sum = np.sum(all_values)
            if np.isclose(value_sum, 0.0):
                raise np.linalg.LinAlgError(
                    "generalized eigenvalues sum to zero and cannot be normalized"
                )
            all_values = all_values * 100.0 / value_sum

            retained_vectors = all_vectors[:, :ncomps]
            eigenvalues[:, frequency_index, subject] = all_values[:ncomps]
            eigenvectors[:, :, frequency_index, subject] = retained_vectors
            time_series[:, :, frequency_index, subject] = (
                retained_vectors.T @ broadband
            )

            retained_patterns = covariance_s @ retained_vectors
            for component in range(ncomps):
                pattern = retained_patterns[:, component]
                largest = int(np.argmax(np.abs(pattern)))
                sign = np.sign(pattern[largest])
                if sign != 0:
                    pattern = pattern * sign
                patterns[:, component, frequency_index, subject] = pattern

    return FREQNESSResult(
        evals=eigenvalues,
        evecs=eigenvectors,
        pats=patterns,
        ts=time_series,
        frex=frequencies,
        fwhm=widths,
        srate=sampling_rate,
        duration=n_samples / sampling_rate,
        bad_segments=bad_one_based,
        regularisation=float(regularisation),
        scale_factors=scale_factors,
    )


__all__ = ["FREQNESS_NetworkEstimation", "FREQNESSResult"]
