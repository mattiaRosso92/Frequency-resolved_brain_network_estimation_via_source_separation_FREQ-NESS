"""Model exponential decay of FREQ-NESS eigenvalues across frequency."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


@dataclass(slots=True)
class FREQNESSExponentialGoodFit:
    """Goodness-of-fit information returned by FREQNESS_ExponentialDK."""

    R2: FloatArray


def _field(container: Any, name: str, *, required: bool = True) -> Any:
    if isinstance(container, Mapping):
        if name in container:
            return container[name]
    elif hasattr(container, name):
        return getattr(container, name)
    if required:
        raise ValueError(f"FREQ.{name} is required")
    return None


def _eigenspectrum_array(FREQ: Any) -> FloatArray:
    values = np.asarray(_field(FREQ, "evals"))
    if values.size == 0:
        raise ValueError("FREQ.evals is empty")
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise TypeError("FREQ.evals must contain real numeric values")
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, np.newaxis, np.newaxis]
    elif values.ndim == 2:
        values = values[:, :, np.newaxis]
    elif values.ndim != 3:
        raise ValueError(
            "FREQ.evals must be 1D, 2D, or 3D "
            "(components, frequencies[, participants])"
        )
    if min(values.shape) == 0:
        raise ValueError("FREQ.evals is empty")
    return values


def _frequency_axis(FREQ: Any, nfrequencies: int) -> tuple[FloatArray, str, bool]:
    raw = _field(FREQ, "frex", required=False)
    if raw is not None:
        frequencies = np.asarray(raw, dtype=float)
        if (
            frequencies.ndim in (1, 2)
            and (frequencies.ndim == 1 or 1 in frequencies.shape)
            and frequencies.size == nfrequencies
            and np.all(np.isfinite(frequencies))
        ):
            return frequencies.reshape(-1).copy(), "Frequency (Hz)", True
    return (
        np.arange(1, nfrequencies + 1, dtype=float),
        "Frequency index",
        False,
    )


def _component_index(which_comp: Any, ncomponents: int) -> int:
    if which_comp is None:
        return 0
    valid = (
        isinstance(which_comp, (int, np.integer))
        and not isinstance(which_comp, (bool, np.bool_))
        and 1 <= int(which_comp) <= ncomponents
    )
    if not valid:
        warnings.warn(
            f"which_comp must be an integer between 1 and {ncomponents}. "
            "Defaulting to component 1.",
            UserWarning,
            stacklevel=3,
        )
        return 0
    return int(which_comp) - 1


def _fit_indices(
    x_axis: FloatArray,
    has_frequencies: bool,
    range2fit: ArrayLike | None,
) -> NDArray[np.int64]:
    if range2fit is None:
        return np.arange(x_axis.size, dtype=np.int64)
    requested = np.asarray(range2fit, dtype=float)
    if requested.size == 0:
        return np.arange(x_axis.size, dtype=np.int64)
    if not has_frequencies:
        raise ValueError("range2fit was provided in Hz, but FREQ.frex is missing")
    if requested.ndim != 1 or requested.size != 2:
        raise ValueError("range2fit must be a two-element vector such as [8, 12]")
    if not np.all(np.isfinite(requested)):
        raise ValueError("range2fit must contain finite values")
    if requested[0] > requested[1]:
        raise ValueError("range2fit boundaries must be in ascending order")

    first = int(np.argmin(np.abs(x_axis - requested[0])))
    last = int(np.argmin(np.abs(x_axis - requested[1])))
    if first > last:
        raise ValueError("range2fit does not map to an ascending frequency range")
    matched = np.array([x_axis[first], x_axis[last]])
    if not np.all(np.isclose(requested, matched, atol=1e-3, rtol=0.0)):
        warnings.warn(
            "Some requested frequencies in range2fit are not present in "
            "FREQ.frex. Using closest matches instead.",
            UserWarning,
            stacklevel=3,
        )
    return np.arange(first, last + 1, dtype=np.int64)


def _fit_exponential(
    x_values: FloatArray,
    y_values: FloatArray,
) -> tuple[float, float, float]:
    mask = (y_values > 0) & np.isfinite(y_values) & np.isfinite(x_values)
    if np.count_nonzero(mask) < 3:
        return np.nan, np.nan, np.nan

    x_fit = x_values[mask]
    y_fit = y_values[mask]
    design = np.column_stack([x_fit, np.ones(x_fit.size)])
    slope, intercept = np.linalg.lstsq(design, np.log(y_fit), rcond=None)[0]
    decay = float(-slope)
    amplitude = float(np.exp(intercept))
    predicted = amplitude * np.exp(-decay * x_fit)
    sse = float(np.sum((y_fit - predicted) ** 2))
    sst = float(np.sum((y_fit - np.mean(y_fit)) ** 2))
    r_squared = 1.0 - sse / sst if sst > 0 else np.nan
    return decay, amplitude, r_squared


def _nanmean(values: FloatArray, axis: int) -> FloatArray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(values, axis=axis)


def _set_x_limits(axis: Any, x_axis: FloatArray) -> None:
    if x_axis.size == 1:
        padding = max(abs(float(x_axis[0])) * 0.05, 0.5)
        axis.set_xlim(x_axis[0] - padding, x_axis[0] + padding)
    else:
        axis.set_xlim(float(np.min(x_axis)), float(np.max(x_axis)))


def _style_axis(axis: Any, x_axis: FloatArray, xlabel: str, title: str) -> None:
    axis.set_xlabel(xlabel, fontweight="bold")
    axis.set_ylabel("Eigenvalue", fontweight="bold")
    axis.set_title(title, fontweight="bold")
    axis.grid(True, which="both", alpha=0.25)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    _set_x_limits(axis, x_axis)


def _plot_decay(
    eigenspectrum: FloatArray,
    x_axis: FloatArray,
    xlabel: str,
    fit_indices: NDArray[np.int64],
    component_index: int,
    decay_coefficients: FloatArray,
    amplitudes: FloatArray,
    group_decay: float,
    group_amplitude: float,
    *,
    plot_all: bool,
) -> list[Any]:
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "FREQNESS_ExponentialDK plotting requires matplotlib. "
            "Install freqness[plot]."
        ) from error

    nsubjects = eigenspectrum.shape[1]
    colormap = plt.get_cmap("viridis")
    spectrum_color = colormap(0.12)
    fit_color = colormap(0.78)
    figures: list[Any] = []

    average = _nanmean(eigenspectrum, axis=1)
    if nsubjects > 1:
        sem = np.nanstd(eigenspectrum, axis=1, ddof=1) / np.sqrt(nsubjects)
    else:
        sem = None

    figure, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
    if sem is not None and not np.all(np.isnan(sem)):
        axis.fill_between(
            x_axis,
            average - sem,
            average + sem,
            color="black",
            alpha=0.1,
            linewidth=0,
        )
    spectrum_line = axis.plot(
        x_axis,
        average,
        color=spectrum_color,
        linewidth=2,
        label="Mean eigenspectrum",
    )[0]
    handles = [spectrum_line]
    if np.isfinite(group_decay) and np.isfinite(group_amplitude):
        fit_line = axis.plot(
            x_axis[fit_indices],
            group_amplitude * np.exp(-group_decay * x_axis[fit_indices]),
            color=fit_color,
            linewidth=2.4,
            label="Exponential fit",
        )[0]
        handles.append(fit_line)
    range_line = axis.axvline(
        x_axis[fit_indices[0]],
        color="black",
        linestyle="--",
        linewidth=1,
        label="Fit range",
    )
    axis.axvline(
        x_axis[fit_indices[-1]],
        color="black",
        linestyle="--",
        linewidth=1,
    )
    handles.append(range_line)
    _style_axis(
        axis,
        x_axis,
        xlabel,
        f"Exponential decay fit - Component #{component_index + 1}",
    )
    axis.legend(handles=handles, loc="upper right")
    figures.append(figure)

    if plot_all:
        for subject in range(nsubjects):
            figure, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
            spectrum_line = axis.plot(
                x_axis,
                eigenspectrum[:, subject],
                color=spectrum_color,
                linewidth=1.5,
                label="Eigenspectrum",
            )[0]
            handles = [spectrum_line]
            if np.isfinite(decay_coefficients[subject, 0]) and np.isfinite(
                amplitudes[subject, 0]
            ):
                fit_line = axis.plot(
                    x_axis[fit_indices],
                    amplitudes[subject, 0]
                    * np.exp(
                        -decay_coefficients[subject, 0] * x_axis[fit_indices]
                    ),
                    color=fit_color,
                    linewidth=2.4,
                    label="Exponential fit",
                )[0]
                handles.append(fit_line)
            range_line = axis.axvline(
                x_axis[fit_indices[0]],
                color="black",
                linestyle="--",
                linewidth=1,
                label="Fit range",
            )
            axis.axvline(
                x_axis[fit_indices[-1]],
                color="black",
                linestyle="--",
                linewidth=1,
            )
            handles.append(range_line)
            _style_axis(
                axis,
                x_axis,
                xlabel,
                f"Participant #{subject + 1} - Exponential decay fit "
                f"(component {component_index + 1})",
            )
            axis.legend(handles=handles, loc="upper right")
            figures.append(figure)
    return figures


def FREQNESS_ExponentialDK(
    FREQ: Any,
    *,
    which_comp: int | None = None,
    range2fit: ArrayLike | None = None,
    plot_all: bool = False,
    show: bool = True,
) -> tuple[FloatArray, FREQNESSExponentialGoodFit]:
    """Fit exponential eigenvalue decay across frequencies.

    The model is ``y = A * exp(-lambda * x)``. ``which_comp`` deliberately
    uses MATLAB-compatible one-based component numbering. Fits require at
    least three positive finite eigenvalues. R² is evaluated in the original,
    non-logarithmic eigenvalue space, matching the MATLAB implementation.

    Returns
    -------
    decayCoeff, goodFit
        ``decayCoeff`` and ``goodFit.R2`` both have shape
        ``(participants, 1)``.
    """
    for value, name in ((plot_all, "plot_all"), (show, "show")):
        if not isinstance(value, (bool, np.bool_)):
            raise TypeError(f"{name} must be a boolean")

    values = _eigenspectrum_array(FREQ)
    ncomponents, nfrequencies, nsubjects = values.shape
    component_index = _component_index(which_comp, ncomponents)
    x_axis, xlabel, has_frequencies = _frequency_axis(FREQ, nfrequencies)
    fit_indices = _fit_indices(x_axis, has_frequencies, range2fit)
    component_spectrum = values[component_index, :, :]

    decay_coefficients = np.full((nsubjects, 1), np.nan, dtype=float)
    amplitudes = np.full((nsubjects, 1), np.nan, dtype=float)
    r_squared = np.full((nsubjects, 1), np.nan, dtype=float)
    for subject in range(nsubjects):
        decay, amplitude, subject_r2 = _fit_exponential(
            x_axis[fit_indices],
            component_spectrum[fit_indices, subject],
        )
        decay_coefficients[subject, 0] = decay
        amplitudes[subject, 0] = amplitude
        r_squared[subject, 0] = subject_r2

    group_spectrum = _nanmean(component_spectrum, axis=1)
    group_decay, group_amplitude, _ = _fit_exponential(
        x_axis[fit_indices], group_spectrum[fit_indices]
    )
    _plot_decay(
        component_spectrum,
        x_axis,
        xlabel,
        fit_indices,
        component_index,
        decay_coefficients,
        amplitudes,
        group_decay,
        group_amplitude,
        plot_all=bool(plot_all),
    )
    if show:
        import matplotlib.pyplot as plt

        plt.show()

    return decay_coefficients, FREQNESSExponentialGoodFit(R2=r_squared)


__all__ = ["FREQNESS_ExponentialDK", "FREQNESSExponentialGoodFit"]
