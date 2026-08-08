"""Event-related induced responses of frequency-resolved networks.

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
from scipy.signal import hilbert

from .filterFGx import filterFGx

if TYPE_CHECKING:
    from matplotlib.figure import Figure


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(slots=True)
class FREQNESSInducedResponsesResult:
    """Outputs of :func:`FREQNESS_InducedResponses`."""

    power: FloatArray
    time: FloatArray
    frex: FloatArray
    fwhm: FloatArray
    comps: IntArray
    events: list[IntArray]
    ntrials: IntArray
    epoch_window: FloatArray
    baseline_window: FloatArray
    srate: float
    power_trials: list[FloatArray] | None
    power_scaled: FloatArray | None
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


def _time_window(value: ArrayLike, name: str) -> FloatArray:
    window = np.asarray(value)
    if (
        window.size != 2
        or not np.issubdtype(window.dtype, np.number)
        or np.iscomplexobj(window)
    ):
        raise ValueError(f"{name} must contain two finite increasing values")
    window = np.asarray(window, dtype=float).reshape(-1)
    if not np.all(np.isfinite(window)) or window[0] >= window[1]:
        raise ValueError(f"{name} must contain two finite increasing values")
    return window


def _matlab_round(value: float) -> int:
    """Round halves away from zero, matching MATLAB rather than NumPy."""
    magnitude = np.floor(abs(float(value)) + 0.5)
    return int(np.copysign(magnitude, value))


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
    if time_series.ndim == 3:
        time_series = time_series[:, :, :, np.newaxis]
    if time_series.ndim != 4 or min(time_series.shape) == 0:
        raise ValueError(
            "FREQ.ts must have shape "
            "(components, time, frequencies, participants)"
        )

    frequencies = np.asarray(_field(FREQ, "frex"))
    widths = np.asarray(_field(FREQ, "fwhm"))
    for values, name in ((frequencies, "FREQ.frex"), (widths, "FREQ.fwhm")):
        if (
            values.size == 0
            or not np.issubdtype(values.dtype, np.number)
            or np.iscomplexobj(values)
        ):
            raise TypeError(f"{name} must contain real numeric values")
    frequencies = np.asarray(frequencies, dtype=float).reshape(-1)
    widths = np.asarray(widths, dtype=float).reshape(-1)
    nfrequencies = time_series.shape[2]
    if frequencies.size != nfrequencies:
        raise ValueError(
            "Length of FREQ.frex does not match the frequency dimension of FREQ.ts"
        )
    if widths.size != nfrequencies:
        raise ValueError("FREQ.fwhm must contain one value per frequency")
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0):
        raise ValueError("FREQ.frex must contain positive finite frequencies")
    if not np.all(np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError("FREQ.fwhm must contain positive finite filter widths")

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
    if np.any(frequencies >= sampling_rate / 2.0):
        raise ValueError("FREQ.frex must be lower than the Nyquist frequency")

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


def _event_vector(value: Any, subject: int, ntime: int) -> IntArray:
    raw = np.asarray(value)
    if raw.ndim == 0:
        raw = raw.reshape(1)
    if (
        raw.size == 0
        or raw.ndim != 1
        or not np.issubdtype(raw.dtype, np.number)
        or np.iscomplexobj(raw)
        or np.issubdtype(raw.dtype, np.bool_)
        or not np.all(np.isfinite(raw))
        or not np.all(raw == np.round(raw))
    ):
        raise ValueError(
            f"events for participant {subject + 1} must contain finite integer "
            "sample indices"
        )
    events = raw.astype(np.int64)
    if np.any(events < 1) or np.any(events > ntime):
        raise ValueError(
            f"events for participant {subject + 1} contain samples outside "
            "the FREQ.ts time range"
        )
    stable_unique = np.asarray(
        list(dict.fromkeys(events.tolist())),
        dtype=np.int64,
    )
    if stable_unique.size != events.size:
        warnings.warn(
            f"Participant {subject + 1} events contain duplicate samples; "
            "duplicates were removed.",
            UserWarning,
            stacklevel=3,
        )
    return stable_unique


def _normalize_events(events: Any, nsubjects: int, ntime: int) -> list[IntArray]:
    if nsubjects == 1:
        candidate = np.asarray(events)
        if candidate.ndim == 1 and np.issubdtype(candidate.dtype, np.number):
            values = [events]
        elif isinstance(events, (list, tuple)) and len(events) == 1:
            values = [events[0]]
        else:
            raise ValueError(
                "events must be a numeric vector for one participant"
            )
    else:
        if not isinstance(events, (list, tuple)) or len(events) != nsubjects:
            raise ValueError(
                f"events must contain one vector per participant ({nsubjects} vectors)"
            )
        values = list(events)
    return [
        _event_vector(value, subject, ntime)
        for subject, value in enumerate(values)
    ]


def _nanmean(values: FloatArray, axis: int) -> FloatArray:
    finite = np.isfinite(values)
    counts = np.sum(finite, axis=axis)
    totals = np.sum(np.where(finite, values, 0.0), axis=axis)
    result = np.full(np.shape(totals), np.nan, dtype=float)
    np.divide(totals, counts, out=result, where=counts > 0)
    return result


def _scaled_power(power: FloatArray, components: IntArray) -> FloatArray:
    scaled = np.full_like(power, np.nan)
    for subject in range(power.shape[3]):
        for component in range(power.shape[0]):
            values = power[component, :, :, subject]
            finite_values = values[np.isfinite(values)]
            maximum = float(np.max(finite_values)) if finite_values.size else np.nan
            if np.isfinite(maximum) and maximum != 0.0:
                scaled[component, :, :, subject] = values / maximum
            else:
                warnings.warn(
                    f"Participant {subject + 1}, component "
                    f"{components[component]}: maximum scaling was not applied "
                    "because the maximum power was invalid or zero.",
                    UserWarning,
                    stacklevel=3,
                )
    return scaled


def _centers_to_edges(values: FloatArray, fallback_width: float) -> FloatArray:
    if values.size == 1:
        half_width = fallback_width / 2.0
        return np.array([values[0] - half_width, values[0] + half_width])
    midpoints = (values[:-1] + values[1:]) / 2.0
    return np.concatenate(
        [
            [values[0] - (midpoints[0] - values[0])],
            midpoints,
            [values[-1] + (values[-1] - midpoints[-1])],
        ]
    )


def _plot_response(
    values: FloatArray,
    time: FloatArray,
    frequencies: FloatArray,
    sampling_rate: float,
    title: str,
) -> Figure:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "Induced-response plotting requires the optional 'plot' dependencies."
        ) from exc

    figure, axis = plt.subplots()
    finite_values = values[np.isfinite(values)]
    limit = float(np.max(np.abs(finite_values))) if finite_values.size else 1.0
    if not np.isfinite(limit) or limit == 0.0:
        limit = 1.0
    time_edges = _centers_to_edges(time, 1.0 / sampling_rate)
    frequency_fallback = max(float(frequencies[0]) * 0.1, 0.5)
    frequency_edges = _centers_to_edges(frequencies, frequency_fallback)
    image = axis.pcolormesh(
        time_edges,
        frequency_edges,
        values.T,
        shading="flat",
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
    )
    axis.axvline(0.0, color="black", linestyle="--", linewidth=1.0)
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Frequency (Hz)")
    axis.set_title(title)
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label("Power change (dB)")
    figure.tight_layout()
    return figure


def FREQNESS_InducedResponses(
    FREQ: Any,
    events: Any,
    epoch_window: ArrayLike,
    baseline_window: ArrayLike,
    *,
    which_comp: ArrayLike | None = None,
    frex2model: ArrayLike | None = None,
    keep_trials: bool = False,
    scale_max: bool = False,
    plot_avg: bool = False,
    plot_all: bool = False,
    show: bool = True,
) -> FREQNESSInducedResponsesResult:
    """Compute trial-wise baseline-normalized induced network power.

    Event samples and component numbers deliberately remain one-based for
    MATLAB compatibility. Epoch endpoints are inclusive, while the baseline
    uses the half-open interval ``[start, end)``. Continuous matched-frequency
    filtering precedes epoching, and trial power is converted to decibels before
    averaging so non-phase-locked induced activity is retained.
    """
    keep_trials = _boolean(keep_trials, "keep_trials")
    scale_max = _boolean(scale_max, "scale_max")
    plot_avg = _boolean(plot_avg, "plot_avg")
    plot_all = _boolean(plot_all, "plot_all")
    show = _boolean(show, "show")

    epoch = _time_window(epoch_window, "epoch_window")
    baseline = _time_window(baseline_window, "baseline_window")
    if baseline[0] < epoch[0] or baseline[1] > epoch[1]:
        raise ValueError("baseline_window must fall entirely within epoch_window")

    time_series, frequencies, widths, sampling_rate = _canonical_inputs(FREQ)
    component_numbers, component_indices = _component_indices(
        which_comp,
        time_series.shape[0],
    )
    frequency_indices = _frequency_indices(frequencies, frex2model)
    selected_frequencies = frequencies[frequency_indices]
    selected_widths = widths[frequency_indices]
    nsubjects = time_series.shape[3]
    ntime_continuous = time_series.shape[1]
    all_events = _normalize_events(events, nsubjects, ntime_continuous)

    first_offset = _matlab_round(epoch[0] * sampling_rate)
    last_offset = _matlab_round(epoch[1] * sampling_rate)
    epoch_offsets = np.arange(first_offset, last_offset + 1, dtype=np.int64)
    time = epoch_offsets.astype(float) / sampling_rate
    baseline_indices = (time >= baseline[0]) & (time < baseline[1])
    if not np.any(baseline_indices):
        raise ValueError(
            "baseline_window does not contain any samples at the current "
            "sampling rate"
        )

    ncomponents = component_indices.size
    nfrequencies = frequency_indices.size
    nepoch = epoch_offsets.size
    power_average = np.full(
        (ncomponents, nepoch, nfrequencies, nsubjects),
        np.nan,
        dtype=float,
    )
    valid_events: list[IntArray] = []
    ntrials = np.zeros(nsubjects, dtype=np.int64)
    power_trials: list[FloatArray] | None = [] if keep_trials else None

    for subject in range(nsubjects):
        subject_events = all_events[subject]
        valid = (
            (subject_events + first_offset >= 1)
            & (subject_events + last_offset <= ntime_continuous)
        )
        kept_events = subject_events[valid]
        excluded = int(np.sum(~valid))
        if excluded:
            warnings.warn(
                f"Participant {subject + 1}: excluding {excluded} event(s) "
                "because the requested epoch extends outside FREQ.ts.",
                UserWarning,
                stacklevel=2,
            )
        if kept_events.size == 0:
            raise ValueError(
                f"Participant {subject + 1} has no valid events after "
                "epoch-boundary checks"
            )
        valid_events.append(kept_events)
        ntrials[subject] = kept_events.size
        epoch_indices = (
            kept_events[:, np.newaxis] - 1 + epoch_offsets[np.newaxis, :]
        )
        subject_trials = (
            np.full(
                (ncomponents, nepoch, nfrequencies, kept_events.size),
                np.nan,
                dtype=float,
            )
            if keep_trials
            else None
        )

        for output_frequency, source_frequency in enumerate(frequency_indices):
            continuous = time_series[
                component_indices,
                :,
                source_frequency,
                subject,
            ]
            if not np.all(np.isfinite(continuous)):
                raise ValueError(
                    f"FREQ.ts contains non-finite values for participant "
                    f"{subject + 1} at {frequencies[source_frequency]:.3f} Hz"
                )
            filtered, _, _ = filterFGx(
                continuous,
                sampling_rate,
                float(frequencies[source_frequency]),
                float(widths[source_frequency]),
            )
            continuous_power = np.abs(hilbert(filtered, axis=1)) ** 2
            epoched_power = continuous_power[:, epoch_indices]
            baseline_power = _nanmean(
                epoched_power[:, :, baseline_indices],
                axis=2,
            )
            invalid_baselines = ~np.isfinite(baseline_power) | (baseline_power <= 0)
            with np.errstate(divide="ignore", invalid="ignore"):
                normalized = 10.0 * np.log10(
                    epoched_power / baseline_power[:, :, np.newaxis]
                )
            normalized[invalid_baselines] = np.nan

            for component in range(ncomponents):
                invalid_count = int(np.sum(invalid_baselines[component]))
                if invalid_count:
                    warnings.warn(
                        f"Participant {subject + 1}, component "
                        f"{component_numbers[component]}, frequency "
                        f"{frequencies[source_frequency]:.3f} Hz: "
                        f"{invalid_count} trial(s) had invalid baseline power "
                        "and were retained as NaN.",
                        UserWarning,
                        stacklevel=2,
                    )

            power_average[:, :, output_frequency, subject] = _nanmean(
                normalized,
                axis=1,
            )
            if subject_trials is not None:
                subject_trials[:, :, output_frequency, :] = np.transpose(
                    normalized,
                    (0, 2, 1),
                )

        if power_trials is not None and subject_trials is not None:
            power_trials.append(subject_trials)

    power_scaled = (
        _scaled_power(power_average, component_numbers) if scale_max else None
    )
    figures: list[Figure] = []
    if plot_avg:
        group_average = _nanmean(power_average, axis=3)
        for component in range(ncomponents):
            figures.append(
                _plot_response(
                    group_average[component],
                    time,
                    selected_frequencies,
                    sampling_rate,
                    f"Induced responses - component #{component_numbers[component]}",
                )
            )
    if plot_all:
        for subject in range(nsubjects):
            for component in range(ncomponents):
                figures.append(
                    _plot_response(
                        power_average[component, :, :, subject],
                        time,
                        selected_frequencies,
                        sampling_rate,
                        f"Participant #{subject + 1} - component "
                        f"#{component_numbers[component]}",
                    )
                )
    if figures and show:
        import matplotlib.pyplot as plt

        plt.show()

    return FREQNESSInducedResponsesResult(
        power=power_average,
        time=time,
        frex=selected_frequencies.copy(),
        fwhm=selected_widths.copy(),
        comps=component_numbers.copy(),
        events=valid_events,
        ntrials=ntrials,
        epoch_window=epoch,
        baseline_window=baseline,
        srate=sampling_rate,
        power_trials=power_trials,
        power_scaled=power_scaled,
        figures=figures,
    )


__all__ = ["FREQNESS_InducedResponses", "FREQNESSInducedResponsesResult"]
