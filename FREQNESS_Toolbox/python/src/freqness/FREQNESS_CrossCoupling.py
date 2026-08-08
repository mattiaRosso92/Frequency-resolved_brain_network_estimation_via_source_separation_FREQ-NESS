"""Neural phase-amplitude cross-frequency coupling for FREQ-NESS networks.

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
from typing import Any, Mapping
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal import hilbert

from .filterFGx import filterFGx


FloatArray = NDArray[np.float64]


@dataclass(slots=True)
class FREQNESSCrossCouplingResult:
    """Outputs of :func:`FREQNESS_CrossCoupling`."""

    PAC_all: FloatArray
    PAC_avg: FloatArray
    sAmpl: FloatArray
    pShift: FloatArray
    dcOff: FloatArray
    mFrex: FloatArray
    goodFit: FloatArray
    carrier_frex: FloatArray
    lfo_freq: float
    comp: int
    phase_edges: FloatArray
    phase_centers: FloatArray
    fitted_PAC: FloatArray
    coefficients: FloatArray
    amplitude_raw: FloatArray
    amplitude_normalized: FloatArray
    preferred_phase: FloatArray
    dc_offset: FloatArray
    mse: FloatArray
    r2: FloatArray
    valid_bins: NDArray[np.int64]
    lfo_phase: FloatArray
    figures: list[Any]


@dataclass(slots=True)
class _HarmonicFit:
    coefficients: FloatArray
    fitted: FloatArray
    amplitude: float
    amplitude_normalized: float
    preferred_phase: float
    offset: float
    mse: float
    r2: float
    nvalid: int


def _field(container: Any, name: str, *, required: bool = True) -> Any:
    if isinstance(container, Mapping):
        if name in container:
            return container[name]
    elif hasattr(container, name):
        return getattr(container, name)
    if required:
        raise ValueError(f"FREQ.{name} is required")
    return None


def _vector(values: Any, name: str) -> FloatArray:
    vector = np.asarray(values, dtype=float)
    if vector.ndim not in (1, 2) or (vector.ndim == 2 and 1 not in vector.shape):
        raise ValueError(f"{name} must be a vector")
    vector = vector.reshape(-1)
    if vector.size == 0 or not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain finite values")
    return vector


def _canonical_inputs(
    FREQ: Any,
) -> tuple[FloatArray, FloatArray, float, FloatArray]:
    frequencies = _vector(_field(FREQ, "frex"), "FREQ.frex")
    if np.any(frequencies <= 0) or np.any(np.diff(frequencies) <= 0):
        raise ValueError("FREQ.frex must be positive and strictly increasing")

    time_series = np.asarray(_field(FREQ, "ts"))
    if not np.issubdtype(time_series.dtype, np.number) or np.iscomplexobj(time_series):
        raise TypeError("FREQ.ts must contain real numeric values")
    time_series = np.asarray(time_series, dtype=float)
    if time_series.ndim == 3:
        time_series = time_series[:, :, :, np.newaxis]
    if time_series.ndim != 4:
        raise ValueError(
            "FREQ.ts must have shape "
            "(components, time, frequencies[, participants])"
        )
    if min(time_series.shape) == 0 or time_series.shape[1] < 2:
        raise ValueError("FREQ.ts is empty or contains fewer than two time samples")
    if time_series.shape[2] != frequencies.size:
        raise ValueError("FREQ.ts and FREQ.frex have incompatible dimensions")
    if not np.all(np.isfinite(time_series)):
        raise ValueError("FREQ.ts contains NaN or infinite values")

    sampling_rate = float(_field(FREQ, "srate"))
    if not np.isfinite(sampling_rate) or sampling_rate <= 0:
        raise ValueError("FREQ.srate must be a positive finite scalar")
    if np.any(frequencies >= sampling_rate / 2):
        raise ValueError("FREQ.frex must be below the Nyquist frequency")

    widths = _vector(_field(FREQ, "fwhm"), "FREQ.fwhm")
    if widths.size != frequencies.size:
        raise ValueError("FREQ.fwhm must contain one value per FREQ.frex entry")
    if np.any(widths <= 0):
        raise ValueError("FREQ.fwhm must contain positive values")
    return frequencies, time_series, sampling_rate, widths


def _component_index(which_comp: Any, ncomponents: int) -> int:
    valid = (
        isinstance(which_comp, (int, np.integer))
        and not isinstance(which_comp, (bool, np.bool_))
        and 1 <= int(which_comp) <= ncomponents
    )
    if not valid:
        raise ValueError(
            f"which_comp must be an integer between 1 and {ncomponents}"
        )
    return int(which_comp) - 1


def _carrier_indices(
    frequencies: FloatArray,
    lfo_index: int,
    frex2model: ArrayLike | None,
) -> NDArray[np.int64]:
    if frex2model is None or np.asarray(frex2model).size == 0:
        bounds = np.array([frequencies[0], frequencies[-1]], dtype=float)
    else:
        bounds = np.asarray(frex2model, dtype=float)
        if bounds.ndim != 1 or bounds.size != 2:
            raise ValueError("frex2model must be a two-element frequency range")
        if not np.all(np.isfinite(bounds)) or bounds[0] > bounds[1]:
            raise ValueError("frex2model must contain finite ascending boundaries")
    mask = (frequencies >= bounds[0]) & (frequencies <= bounds[1])
    mask[: lfo_index + 1] = False
    indices = np.flatnonzero(mask).astype(np.int64)
    if indices.size == 0:
        raise ValueError("No carrier frequencies remain above the selected LFO")
    return indices


def _first_harmonic_fit(
    phase: ArrayLike,
    power: ArrayLike,
    min_valid_bins: int,
) -> _HarmonicFit:
    phase_values = np.asarray(phase, dtype=float).reshape(-1)
    power_values = np.asarray(power, dtype=float).reshape(-1)
    if phase_values.shape != power_values.shape:
        raise ValueError("phase and power must have matching shapes")
    valid = np.isfinite(phase_values) & np.isfinite(power_values)
    nvalid = int(np.count_nonzero(valid))
    empty = _HarmonicFit(
        coefficients=np.full(3, np.nan),
        fitted=np.full(power_values.shape, np.nan),
        amplitude=np.nan,
        amplitude_normalized=np.nan,
        preferred_phase=np.nan,
        offset=np.nan,
        mse=np.nan,
        r2=np.nan,
        nvalid=nvalid,
    )
    if nvalid < min_valid_bins:
        return empty

    valid_phase = phase_values[valid]
    observed = power_values[valid]
    design = np.column_stack(
        [np.ones(nvalid), np.cos(valid_phase), np.sin(valid_phase)]
    )
    coefficients = np.linalg.lstsq(design, observed, rcond=None)[0]
    fitted = (
        coefficients[0]
        + coefficients[1] * np.cos(phase_values)
        + coefficients[2] * np.sin(phase_values)
    )
    residuals = observed - design @ coefficients
    sse = float(np.sum(residuals**2))
    sst = float(np.sum((observed - np.mean(observed)) ** 2))
    amplitude = float(np.hypot(coefficients[1], coefficients[2]))
    preferred_phase = float(np.arctan2(coefficients[2], coefficients[1]))
    offset = float(coefficients[0])
    normalized = (
        float(100.0 * amplitude / abs(offset))
        if abs(offset) > np.finfo(float).eps
        else np.nan
    )
    return _HarmonicFit(
        coefficients=np.asarray(coefficients, dtype=float),
        fitted=np.asarray(fitted, dtype=float),
        amplitude=amplitude,
        amplitude_normalized=normalized,
        preferred_phase=preferred_phase,
        offset=offset,
        mse=float(np.mean(residuals**2)),
        r2=1.0 - sse / sst if sst > np.finfo(float).eps else np.nan,
        nvalid=nvalid,
    )


def _nanmean(values: FloatArray, axis: int) -> FloatArray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(values, axis=axis)


def _equal_3d_axes(axis: Any, coordinates: FloatArray) -> None:
    spans = np.ptp(coordinates, axis=0)
    spans[spans == 0] = 1.0
    axis.set_box_aspect(spans)
    axis.view_init(elev=30, azim=135)


def _plot_pac_summary(
    phase_centers: FloatArray,
    PAC_avg: FloatArray,
    carrier_frequencies: FloatArray,
    amplitude: FloatArray,
    lfo_frequency: float,
    component: int,
) -> Any:
    import matplotlib.pyplot as plt

    colors = plt.get_cmap("viridis")(
        np.linspace(0.05, 0.95, carrier_frequencies.size)
    )
    figure, axes = plt.subplots(2, 1, figsize=(9, 8), constrained_layout=True)
    for carrier, frequency in enumerate(carrier_frequencies):
        axes[0].plot(
            phase_centers,
            PAC_avg[carrier],
            color=colors[carrier],
            linewidth=1.6,
            label=f"{frequency:g} Hz",
        )
    axes[0].set_xlabel("LFO phase (rad)", fontweight="bold")
    axes[0].set_ylabel("Carrier power (a.u.)", fontweight="bold")
    axes[0].set_title(
        f"PAC across frequencies - Group average: LFO {lfo_frequency:g} Hz "
        f"(component {component})",
        fontweight="bold",
    )
    axes[0].legend(loc="center left", bbox_to_anchor=(1.01, 0.5))

    axes[1].plot(
        carrier_frequencies,
        _nanmean(amplitude, axis=1),
        color=colors[0],
        linewidth=1.8,
    )
    axes[1].set_xlabel("Carrier frequency (Hz)", fontweight="bold")
    axes[1].set_ylabel("First-harmonic amplitude (a.u.)", fontweight="bold")
    axes[1].set_title(
        f"Power modulation: LFO {lfo_frequency:g} Hz (component {component})",
        fontweight="bold",
    )
    for axis in axes:
        axis.grid(True, which="both", alpha=0.25)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    return figure


def _plot_participant_pac(
    phase_centers: FloatArray,
    PAC_all: FloatArray,
    carrier_frequencies: FloatArray,
    amplitude: FloatArray,
    lfo_frequency: float,
    component: int,
) -> list[Any]:
    import matplotlib.pyplot as plt

    colors = plt.get_cmap("viridis")(
        np.linspace(0.05, 0.95, carrier_frequencies.size)
    )
    figures: list[Any] = []
    for subject in range(PAC_all.shape[1]):
        figure, axes = plt.subplots(2, 1, figsize=(9, 8), constrained_layout=True)
        for carrier, frequency in enumerate(carrier_frequencies):
            axes[0].plot(
                phase_centers,
                PAC_all[carrier, subject],
                color=colors[carrier],
                linewidth=1.6,
                label=f"{frequency:g} Hz",
            )
        axes[0].set_xlabel("LFO phase (rad)", fontweight="bold")
        axes[0].set_ylabel("Carrier power (a.u.)", fontweight="bold")
        axes[0].set_title(
            f"PAC - Participant #{subject + 1}: LFO {lfo_frequency:g} Hz "
            f"(component {component})",
            fontweight="bold",
        )
        axes[0].legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
        axes[1].plot(
            carrier_frequencies,
            amplitude[:, subject],
            color=colors[0],
            linewidth=1.8,
        )
        axes[1].set_xlabel("Carrier frequency (Hz)", fontweight="bold")
        axes[1].set_ylabel("First-harmonic amplitude (a.u.)", fontweight="bold")
        axes[1].set_title(
            f"Power modulation - Participant #{subject + 1}",
            fontweight="bold",
        )
        for axis in axes:
            axis.grid(True, which="both", alpha=0.25)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
        figures.append(figure)
    return figures


def _normalise_pattern(pattern: FloatArray) -> FloatArray:
    values = np.abs(np.asarray(pattern, dtype=float))
    finite = np.isfinite(values)
    maximum = float(np.max(values[finite])) if np.any(finite) else 0.0
    if maximum > 0:
        values[finite] /= maximum
    else:
        values[finite] = 0.0
    return values


def _plot_spatial_patterns(
    FREQ: Any,
    MNI: ArrayLike,
    component_index: int,
    lfo_index: int,
    carrier_indices: NDArray[np.int64],
    PAC_avg: FloatArray,
    frequencies: FloatArray,
) -> Any | None:
    import matplotlib.pyplot as plt

    raw_patterns = _field(FREQ, "pats", required=False)
    if raw_patterns is None:
        warnings.warn(
            "MNI coordinates were provided but FREQ.pats is unavailable; "
            "skipping spatial plots.",
            UserWarning,
            stacklevel=3,
        )
        return None
    patterns = np.asarray(raw_patterns, dtype=float)
    if patterns.ndim == 3:
        patterns = patterns[:, :, :, np.newaxis]
    if patterns.ndim != 4 or patterns.shape[1] <= component_index:
        warnings.warn("FREQ.pats has incompatible dimensions; skipping spatial plots.")
        return None
    if patterns.shape[2] != frequencies.size:
        warnings.warn("FREQ.pats has incompatible frequencies; skipping spatial plots.")
        return None

    coordinates = np.asarray(MNI, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape != (patterns.shape[0], 3):
        raise ValueError("MNI must have shape (voxels, 3) matching FREQ.pats")
    valid = np.all(np.isfinite(coordinates), axis=1)
    if not np.any(valid):
        raise ValueError("MNI contains no valid coordinates")

    lfo_pattern = _nanmean(patterns[:, component_index, lfo_index, :], axis=1)
    carrier_mean_power = _nanmean(PAC_avg, axis=1)
    if not np.any(np.isfinite(carrier_mean_power)):
        warnings.warn("PAC estimates are unavailable; skipping spatial plots.")
        return None
    peak_local = int(np.nanargmax(carrier_mean_power))
    peak_index = int(carrier_indices[peak_local])
    peak_pattern = _nanmean(patterns[:, component_index, peak_index, :], axis=1)
    visual_patterns = [_normalise_pattern(lfo_pattern), _normalise_pattern(peak_pattern)]
    titles = [
        f"LFO network ({frequencies[lfo_index]:g} Hz, component {component_index + 1})",
        f"Peak carrier ({frequencies[peak_index]:g} Hz, component {component_index + 1})",
    ]

    figure = plt.figure(figsize=(12, 5.5), constrained_layout=True)
    for index, (pattern, title) in enumerate(zip(visual_patterns, titles, strict=True)):
        axis = figure.add_subplot(1, 2, index + 1, projection="3d")
        scatter = axis.scatter(
            coordinates[valid, 0],
            coordinates[valid, 1],
            coordinates[valid, 2],
            s=20,
            c=pattern[valid],
            cmap="viridis",
            vmin=0,
            vmax=1,
            linewidths=0,
            depthshade=False,
        )
        figure.colorbar(scatter, ax=axis, shrink=0.7, label="Normalized magnitude")
        axis.set_title(title, fontweight="bold")
        axis.set_xlabel("MNI X")
        axis.set_ylabel("MNI Y")
        axis.set_zlabel("MNI Z")
        _equal_3d_axes(axis, coordinates[valid])
    return figure


def FREQNESS_CrossCoupling(
    FREQ: Any,
    lfo_freq: float,
    *,
    MNI: ArrayLike | None = None,
    frex2model: ArrayLike | None = None,
    which_comp: int = 1,
    plot_all: bool = False,
    nbins: int = 37,
    min_valid_bins: int | None = None,
    show: bool = True,
) -> FREQNESSCrossCouplingResult:
    """Compute neural LFO phase-to-neural carrier power coupling.

    The FREQ-NESS signal path mirrors MATLAB: the selected component is
    re-filtered at the LFO and carrier frequencies, LFO phase and carrier power
    are obtained with the analytic signal, and mean carrier power is binned by
    LFO phase. The former external ``sineFit.m`` stage is intentionally replaced
    by a deterministic fixed first-harmonic regression over phase-bin centres.

    ``which_comp`` uses MATLAB-compatible one-based numbering. The legacy
    fields ``sAmpl``, ``pShift``, ``dcOff``, and ``goodFit`` alias raw harmonic
    amplitude, preferred phase, DC offset, and MSE. ``mFrex`` is fixed to
    ``1 / nbins`` cycles per bin, equivalent to one cycle per LFO phase cycle.
    """
    for value, name in ((plot_all, "plot_all"), (show, "show")):
        if not isinstance(value, (bool, np.bool_)):
            raise TypeError(f"{name} must be a boolean")
    if not isinstance(nbins, (int, np.integer)) or isinstance(nbins, (bool, np.bool_)):
        raise TypeError("nbins must be an integer")
    if nbins < 6:
        raise ValueError("nbins must be at least 6")
    if min_valid_bins is None:
        minimum_bins = max(6, int(np.ceil(nbins / 2)))
    else:
        if (
            not isinstance(min_valid_bins, (int, np.integer))
            or isinstance(min_valid_bins, (bool, np.bool_))
            or not 3 <= int(min_valid_bins) <= nbins
        ):
            raise ValueError("min_valid_bins must be an integer from 3 to nbins")
        minimum_bins = int(min_valid_bins)

    frequencies, time_series, sampling_rate, widths = _canonical_inputs(FREQ)
    if (
        not np.isscalar(lfo_freq)
        or not np.isfinite(lfo_freq)
        or float(lfo_freq) <= 0
    ):
        raise ValueError("lfo_freq must be a positive finite scalar")
    component_index = _component_index(which_comp, time_series.shape[0])
    lfo_index = int(np.argmin(np.abs(frequencies - float(lfo_freq))))
    actual_lfo = float(frequencies[lfo_index])
    if not np.isclose(actual_lfo, float(lfo_freq), atol=1e-3, rtol=0.0):
        warnings.warn(
            f"lfo_freq {float(lfo_freq):g} Hz is unavailable; using "
            f"{actual_lfo:g} Hz.",
            UserWarning,
            stacklevel=2,
        )
    carrier_indices = _carrier_indices(frequencies, lfo_index, frex2model)
    carrier_frequencies = frequencies[carrier_indices]
    ncarriers = carrier_indices.size
    nsubjects = time_series.shape[3]
    ntime = time_series.shape[1]

    lfo_phase = np.empty((ntime, nsubjects), dtype=float)
    for subject in range(nsubjects):
        lfo_signal, _, _ = filterFGx(
            time_series[component_index, :, lfo_index, subject],
            sampling_rate,
            actual_lfo,
            float(widths[lfo_index]),
        )
        lfo_phase[:, subject] = np.angle(hilbert(lfo_signal))

    phase_edges = np.linspace(-np.pi, np.pi, nbins + 1)
    phase_centers = phase_edges[:-1] + np.diff(phase_edges) / 2.0
    PAC_all = np.full((ncarriers, nsubjects, nbins), np.nan, dtype=float)
    for carrier, frequency_index in enumerate(carrier_indices):
        for subject in range(nsubjects):
            carrier_signal, _, _ = filterFGx(
                time_series[component_index, :, frequency_index, subject],
                sampling_rate,
                float(frequencies[frequency_index]),
                float(widths[frequency_index]),
            )
            carrier_power = np.abs(hilbert(carrier_signal)) ** 2
            phase = lfo_phase[:, subject]
            for phase_bin in range(nbins):
                in_bin = (phase > phase_edges[phase_bin]) & (
                    phase <= phase_edges[phase_bin + 1]
                )
                if np.any(in_bin):
                    PAC_all[carrier, subject, phase_bin] = float(
                        np.mean(carrier_power[in_bin])
                    )

    PAC_avg = _nanmean(PAC_all, axis=1)
    coefficients = np.full((ncarriers, nsubjects, 3), np.nan, dtype=float)
    fitted_PAC = np.full_like(PAC_all, np.nan)
    amplitude_raw = np.full((ncarriers, nsubjects), np.nan, dtype=float)
    amplitude_normalized = np.full_like(amplitude_raw, np.nan)
    preferred_phase = np.full_like(amplitude_raw, np.nan)
    dc_offset = np.full_like(amplitude_raw, np.nan)
    mse = np.full_like(amplitude_raw, np.nan)
    r2 = np.full_like(amplitude_raw, np.nan)
    valid_bins = np.zeros((ncarriers, nsubjects), dtype=np.int64)
    for carrier in range(ncarriers):
        for subject in range(nsubjects):
            fit = _first_harmonic_fit(
                phase_centers,
                PAC_all[carrier, subject],
                minimum_bins,
            )
            coefficients[carrier, subject] = fit.coefficients
            fitted_PAC[carrier, subject] = fit.fitted
            amplitude_raw[carrier, subject] = fit.amplitude
            amplitude_normalized[carrier, subject] = fit.amplitude_normalized
            preferred_phase[carrier, subject] = fit.preferred_phase
            dc_offset[carrier, subject] = fit.offset
            mse[carrier, subject] = fit.mse
            r2[carrier, subject] = fit.r2
            valid_bins[carrier, subject] = fit.nvalid

    mFrex = np.full_like(amplitude_raw, np.nan)
    mFrex[np.isfinite(amplitude_raw)] = 1.0 / nbins
    figures: list[Any] = [
        _plot_pac_summary(
            phase_centers,
            PAC_avg,
            carrier_frequencies,
            amplitude_raw,
            actual_lfo,
            int(which_comp),
        )
    ]
    if MNI is not None and np.asarray(MNI).size > 0:
        spatial_figure = _plot_spatial_patterns(
            FREQ,
            MNI,
            component_index,
            lfo_index,
            carrier_indices,
            PAC_avg,
            frequencies,
        )
        if spatial_figure is not None:
            figures.append(spatial_figure)
    if plot_all:
        figures.extend(
            _plot_participant_pac(
                phase_centers,
                PAC_all,
                carrier_frequencies,
                amplitude_raw,
                actual_lfo,
                int(which_comp),
            )
        )
    if show:
        import matplotlib.pyplot as plt

        plt.show()

    return FREQNESSCrossCouplingResult(
        PAC_all=PAC_all,
        PAC_avg=PAC_avg,
        sAmpl=amplitude_raw,
        pShift=preferred_phase,
        dcOff=dc_offset,
        mFrex=mFrex,
        goodFit=mse,
        carrier_frex=carrier_frequencies,
        lfo_freq=actual_lfo,
        comp=int(which_comp),
        phase_edges=phase_edges,
        phase_centers=phase_centers,
        fitted_PAC=fitted_PAC,
        coefficients=coefficients,
        amplitude_raw=amplitude_raw,
        amplitude_normalized=amplitude_normalized,
        preferred_phase=preferred_phase,
        dc_offset=dc_offset,
        mse=mse,
        r2=r2,
        valid_bins=valid_bins,
        lfo_phase=lfo_phase,
        figures=figures,
    )


__all__ = ["FREQNESS_CrossCoupling", "FREQNESSCrossCouplingResult"]
