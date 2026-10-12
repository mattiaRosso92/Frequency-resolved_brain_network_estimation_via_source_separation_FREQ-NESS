"""Frequency stability of FREQ-NESS network time series.

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
from typing import Any, Mapping, TYPE_CHECKING
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import rank_filter
from scipy.signal import hilbert

from .filterFGx import filterFGx

if TYPE_CHECKING:
    from matplotlib.figure import Figure


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(slots=True)
class FREQNESSStabilityIndexResult:
    """Outputs of :func:`FREQNESS_StabilityIndex`."""

    stability_index: FloatArray
    mean_frequency: FloatArray
    time: FloatArray
    frex: FloatArray
    fwhm: FloatArray
    comps: IntArray
    srate: float
    median_window: float
    median_samples: int
    frequency_sliding: FloatArray | None
    frequency_sliding_raw: FloatArray | None
    figures: list[Figure]


def _field(container: Any, name: str) -> Any:
    if isinstance(container, Mapping):
        if name in container:
            return container[name]
    elif hasattr(container, name):
        return getattr(container, name)
    raise ValueError(f"FREQ.{name} is missing or empty")


def _boolean(value: Any, name: str) -> bool:
    if not isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a boolean")
    return bool(value)


def _matlab_round_positive(value: float) -> int:
    """Round one positive value as MATLAB does."""
    return int(np.floor(float(value) + 0.5))


def _canonical_inputs(
    FREQ: Any,
) -> tuple[FloatArray, FloatArray, FloatArray, float]:
    time_series_raw = np.asarray(_field(FREQ, "ts"))
    if (
        time_series_raw.size == 0
        or not np.issubdtype(time_series_raw.dtype, np.number)
        or np.iscomplexobj(time_series_raw)
    ):
        raise TypeError("FREQ.ts must contain real numeric values")
    time_series = np.asarray(time_series_raw, dtype=float)
    if time_series.ndim == 2:
        time_series = time_series[:, :, np.newaxis, np.newaxis]
    elif time_series.ndim == 3:
        time_series = time_series[:, :, :, np.newaxis]
    if time_series.ndim != 4 or min(time_series.shape) == 0:
        raise ValueError(
            "FREQ.ts must have shape "
            "(components, time, frequencies, participants)"
        )
    if time_series.shape[1] < 3:
        raise ValueError("FREQ.ts must contain at least three time samples")
    if not np.all(np.isfinite(time_series)):
        raise ValueError("FREQ.ts contains NaN or infinite values")

    frequencies_raw = np.asarray(_field(FREQ, "frex"))
    widths_raw = np.asarray(_field(FREQ, "fwhm"))
    for values, name in (
        (frequencies_raw, "FREQ.frex"),
        (widths_raw, "FREQ.fwhm"),
    ):
        if (
            values.size == 0
            or not np.issubdtype(values.dtype, np.number)
            or np.iscomplexobj(values)
        ):
            raise TypeError(f"{name} must contain real numeric values")
    frequencies = np.asarray(frequencies_raw, dtype=float).reshape(-1)
    widths = np.asarray(widths_raw, dtype=float).reshape(-1)
    if frequencies.size != time_series.shape[2]:
        raise ValueError(
            "Length of FREQ.frex does not match the frequency dimension of FREQ.ts"
        )
    if widths.size != time_series.shape[2]:
        raise ValueError("FREQ.fwhm must contain one value per frequency")

    sampling_rate_raw = np.asarray(_field(FREQ, "srate"))
    if (
        sampling_rate_raw.ndim != 0
        or not np.issubdtype(sampling_rate_raw.dtype, np.number)
        or np.iscomplexobj(sampling_rate_raw)
    ):
        raise TypeError("FREQ.srate must be one positive finite scalar")
    sampling_rate = float(sampling_rate_raw)
    if not np.isfinite(sampling_rate) or sampling_rate <= 0:
        raise ValueError("FREQ.srate must be one positive finite scalar")
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0):
        raise ValueError("FREQ.frex must contain positive finite frequencies")
    if np.any(frequencies >= sampling_rate / 2.0):
        raise ValueError("FREQ.frex must be lower than the Nyquist frequency")
    if not np.all(np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError("FREQ.fwhm must contain positive finite filter widths")

    return time_series, frequencies, widths, sampling_rate


def _component_indices(
    which_comp: ArrayLike | None,
    ncomponents: int,
) -> tuple[IntArray, IntArray]:
    if which_comp is None or np.asarray(which_comp).size == 0:
        one_based = np.array([1], dtype=np.int64)
    else:
        raw = np.asarray(which_comp)
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
            raise ValueError("which_comp must contain finite integer indices")
        one_based = np.asarray(
            list(dict.fromkeys(raw.astype(np.int64).tolist())),
            dtype=np.int64,
        )
    if np.any(one_based < 1) or np.any(one_based > ncomponents):
        raise ValueError(
            f"which_comp must contain indices between 1 and {ncomponents}"
        )
    return one_based, one_based - 1


def _frequency_indices(
    frequencies: FloatArray,
    frex2model: ArrayLike | None,
) -> IntArray:
    if frex2model is None or np.asarray(frex2model).size == 0:
        return np.arange(frequencies.size, dtype=np.int64)

    requested = np.asarray(frex2model)
    if (
        requested.ndim > 1
        or not np.issubdtype(requested.dtype, np.number)
        or np.iscomplexobj(requested)
        or not np.all(np.isfinite(requested))
    ):
        raise ValueError(
            "frex2model must be one frequency or a two-element range"
        )
    requested = np.asarray(requested, dtype=float).reshape(-1)
    if requested.size == 1:
        index = int(np.argmin(np.abs(frequencies - requested[0])))
        if not np.isclose(requested[0], frequencies[index], atol=1e-10, rtol=0.0):
            warnings.warn(
                f"Requested frex2model ({requested[0]:.3f} Hz) is not present "
                f"in FREQ.frex; using the closest available frequency: "
                f"{frequencies[index]:.3f} Hz.",
                UserWarning,
                stacklevel=3,
            )
        return np.array([index], dtype=np.int64)
    if requested.size != 2:
        raise ValueError(
            "frex2model must contain either one frequency or a two-element range"
        )

    requested.sort()
    first = int(np.argmin(np.abs(frequencies - requested[0])))
    last = int(np.argmin(np.abs(frequencies - requested[1])))
    lower, upper = sorted((first, last))
    if not np.all(
        [
            np.any(np.isclose(value, frequencies, atol=1e-10, rtol=0.0))
            for value in requested
        ]
    ):
        warnings.warn(
            "One or both frex2model limits are unavailable; using the closest "
            f"range: {frequencies[lower]:.3f} to {frequencies[upper]:.3f} Hz.",
            UserWarning,
            stacklevel=3,
        )
    return np.arange(lower, upper + 1, dtype=np.int64)


def _moving_median(values: FloatArray, window: int) -> FloatArray:
    """Match MATLAB ``movmedian(values, window, 2)`` for finite 2D data."""
    data = np.asarray(values, dtype=float)
    if data.ndim != 2:
        raise ValueError("values must have shape (signals, time)")
    if not isinstance(window, (int, np.integer)) or isinstance(
        window, (bool, np.bool_)
    ):
        raise TypeError("window must be an integer")
    window = int(window)
    if window < 1 or window > data.shape[1]:
        raise ValueError("window must lie between 1 and the number of time samples")
    if window == 1:
        return data.copy()

    filter_size = (1, window)
    if window % 2:
        filtered = rank_filter(
            data,
            rank=window // 2,
            size=filter_size,
            mode="nearest",
        )
    else:
        lower = rank_filter(
            data,
            rank=window // 2 - 1,
            size=filter_size,
            mode="nearest",
        )
        upper = rank_filter(
            data,
            rank=window // 2,
            size=filter_size,
            mode="nearest",
        )
        filtered = (lower + upper) / 2.0

    left = window // 2
    right = window - left - 1
    for index in range(left):
        filtered[:, index] = np.median(
            data[:, : index + right + 1],
            axis=1,
        )
    for index in range(data.shape[1] - right, data.shape[1]):
        filtered[:, index] = np.median(
            data[:, index - left :],
            axis=1,
        )
    return np.asarray(filtered, dtype=float)


def _plot_stability(
    frequencies: FloatArray,
    values: FloatArray,
    errors: FloatArray | None,
    title: str,
) -> Figure:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "Stability-index plotting requires the optional 'plot' dependencies."
        ) from exc

    figure, axis = plt.subplots()
    if errors is None:
        axis.plot(
            frequencies,
            values,
            "o-",
            color=(0.15, 0.35, 0.65),
            markerfacecolor=(0.15, 0.35, 0.65),
            linewidth=1.5,
        )
    else:
        axis.errorbar(
            frequencies,
            values,
            yerr=errors,
            fmt="o-",
            color=(0.15, 0.35, 0.65),
            markerfacecolor=(0.15, 0.35, 0.65),
            linewidth=1.5,
            capsize=4,
        )
    axis.set_xlabel("Network frequency (Hz)")
    axis.set_ylabel("Stability Index (Hz)")
    axis.set_title(title)
    axis.grid(True, alpha=0.3)
    figure.tight_layout()
    return figure


def FREQNESS_StabilityIndex(
    FREQ: Any,
    *,
    which_comp: ArrayLike | None = None,
    frex2model: ArrayLike | None = None,
    median_window: float = 0.4,
    keep_frequency_sliding: bool = False,
    plot_avg: bool = False,
    plot_all: bool = False,
    show: bool = True,
) -> FREQNESSStabilityIndexResult:
    """Compute the Stability Index of frequency-resolved network activity.

    Component numbers deliberately remain one-based for MATLAB compatibility.
    The Stability Index is the sample standard deviation, in hertz, of the
    moving-median-smoothed instantaneous frequency. Lower values indicate more
    stable oscillatory frequency over time.
    """
    keep_frequency_sliding = _boolean(
        keep_frequency_sliding,
        "keep_frequency_sliding",
    )
    plot_avg = _boolean(plot_avg, "plot_avg")
    plot_all = _boolean(plot_all, "plot_all")
    show = _boolean(show, "show")

    if (
        isinstance(median_window, (bool, np.bool_))
        or not np.isscalar(median_window)
        or not np.isreal(median_window)
    ):
        raise TypeError("median_window must be one positive finite scalar")
    median_window = float(median_window)
    if not np.isfinite(median_window) or median_window <= 0:
        raise ValueError("median_window must be one positive finite duration")

    time_series, frequencies, widths, sampling_rate = _canonical_inputs(FREQ)
    if which_comp is None or np.asarray(which_comp).size == 0:
        print(
            "\nFREQNESS StabilityIndex: no components specified. "
            "Defaulting to component #1."
        )
    component_numbers, component_indices = _component_indices(
        which_comp,
        time_series.shape[0],
    )
    frequency_indices = _frequency_indices(frequencies, frex2model)
    selected_frequencies = frequencies[frequency_indices]
    selected_widths = widths[frequency_indices]

    median_samples = max(
        1,
        _matlab_round_positive(median_window * sampling_rate),
    )
    ntime_frequency = time_series.shape[1] - 1
    if median_samples > ntime_frequency:
        raise ValueError(
            f"median_window corresponds to {median_samples} samples, but "
            f"instantaneous frequency contains only {ntime_frequency} samples"
        )
    effective_window = median_samples / sampling_rate

    ncomponents = component_numbers.size
    nfrequencies = frequency_indices.size
    nsubjects = time_series.shape[3]
    stability_index = np.full(
        (ncomponents, nfrequencies, nsubjects),
        np.nan,
        dtype=float,
    )
    mean_frequency = np.full_like(stability_index, np.nan)
    time = np.arange(1, time_series.shape[1], dtype=float) / sampling_rate

    if keep_frequency_sliding:
        frequency_sliding = np.full(
            (ncomponents, ntime_frequency, nfrequencies, nsubjects),
            np.nan,
            dtype=float,
        )
        frequency_sliding_raw = np.full_like(frequency_sliding, np.nan)
    else:
        frequency_sliding = None
        frequency_sliding_raw = None

    print(
        "\nFREQNESS StabilityIndex: computing frequency stability for "
        f"{ncomponents} component(s), {nfrequencies} frequencies, and "
        f"{nsubjects} participants."
    )
    print(
        f"Moving-median window: {median_samples} samples "
        f"({effective_window:.4g} seconds)."
    )

    for subject in range(nsubjects):
        print(
            f"FREQNESS StabilityIndex: participant "
            f"#{subject + 1}/{nsubjects}."
        )
        for output_frequency, source_frequency in enumerate(frequency_indices):
            component_time_series = time_series[
                component_indices,
                :,
                source_frequency,
                subject,
            ]
            filtered, _, _ = filterFGx(
                component_time_series,
                sampling_rate,
                float(frequencies[source_frequency]),
                float(widths[source_frequency]),
            )
            phase = np.unwrap(np.angle(hilbert(filtered, axis=1)), axis=1)
            raw = sampling_rate * np.diff(phase, axis=1) / (2.0 * np.pi)
            smoothed = _moving_median(raw, median_samples)

            stability_index[:, output_frequency, subject] = np.std(
                smoothed,
                axis=1,
                ddof=1,
            )
            mean_frequency[:, output_frequency, subject] = np.mean(
                smoothed,
                axis=1,
            )
            if frequency_sliding is not None:
                frequency_sliding[:, :, output_frequency, subject] = smoothed
            if frequency_sliding_raw is not None:
                frequency_sliding_raw[:, :, output_frequency, subject] = raw

    figures: list[Figure] = []
    if plot_avg:
        group_mean = np.mean(stability_index, axis=2)
        if nsubjects > 1:
            group_sem = np.std(stability_index, axis=2, ddof=1) / np.sqrt(
                nsubjects
            )
        else:
            group_sem = np.zeros_like(group_mean)
        for component in range(ncomponents):
            figures.append(
                _plot_stability(
                    selected_frequencies,
                    group_mean[component],
                    group_sem[component],
                    f"Stability Index - component #{component_numbers[component]}",
                )
            )
    if plot_all:
        for subject in range(nsubjects):
            for component in range(ncomponents):
                figures.append(
                    _plot_stability(
                        selected_frequencies,
                        stability_index[component, :, subject],
                        None,
                        f"Participant #{subject + 1} - component "
                        f"#{component_numbers[component]}",
                    )
                )
    if figures and show:
        import matplotlib.pyplot as plt

        plt.show()

    return FREQNESSStabilityIndexResult(
        stability_index=stability_index,
        mean_frequency=mean_frequency,
        time=time,
        frex=selected_frequencies.copy(),
        fwhm=selected_widths.copy(),
        comps=component_numbers.copy(),
        srate=sampling_rate,
        median_window=effective_window,
        median_samples=median_samples,
        frequency_sliding=frequency_sliding,
        frequency_sliding_raw=frequency_sliding_raw,
        figures=figures,
    )


__all__ = ["FREQNESS_StabilityIndex", "FREQNESSStabilityIndexResult"]

