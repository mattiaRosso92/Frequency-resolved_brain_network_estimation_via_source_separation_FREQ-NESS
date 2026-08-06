"""Frequency-domain Gaussian filtering used by FREQ-NESS."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal import hilbert


def _nearest_index(values: NDArray[np.floating], target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


def filterFGx(
    data: ArrayLike,
    srate: float,
    f: float,
    fwhm: float,
    showplot: bool = False,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Apply a narrow-band frequency-domain Gaussian filter.

    This is an independent Python translation of the ``filterFGx`` algorithm
    used by the FREQ-NESS MATLAB toolbox. The original filter design and
    implementations are by Mike X Cohen; see
    https://github.com/mikexcohen/GED_tutorial.

    Parameters
    ----------
    data
        Array with shape ``(channels, time)`` or a one-dimensional time series.
    srate
        Sampling rate in Hz.
    f
        Peak frequency of the filter in Hz.
    fwhm
        Full width at half maximum in Hz.
    showplot
        Plot the frequency-domain filter when ``True``.

    Returns
    -------
    filtdat
        Filtered data with the same dimensionality as the input.
    empVals
        Empirical peak frequency, frequency-domain FWHM in Hz, and temporal
        FWHM in milliseconds.
    fx
        Gain-normalized frequency-domain Gaussian.
    """
    values = np.asarray(data)
    input_was_1d = values.ndim == 1
    if input_was_1d:
        values = values[np.newaxis, :]

    if values.ndim != 2 or values.shape[1] < 2:
        raise ValueError("data must have shape (channels, time) with at least two samples")
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise TypeError("data must contain real numeric values")
    if not np.all(np.isfinite(values)):
        raise ValueError("data contains NaN or infinite values")
    if not np.isfinite(srate) or srate <= 0:
        raise ValueError("srate must be a positive finite scalar")
    if not np.isfinite(f) or f <= 0 or f >= srate / 2:
        raise ValueError("f must be positive and lower than the Nyquist frequency")
    if not np.isfinite(fwhm) or fwhm <= 0:
        raise ValueError("fwhm must be a positive finite scalar")

    n_time = values.shape[1]
    hz = np.linspace(0.0, float(srate), n_time)

    normalized_width = fwhm * (2.0 * np.pi - 1.0) / (4.0 * np.pi)
    fx = np.exp(-0.5 * ((hz - f) / normalized_width) ** 2)
    fx /= np.max(fx)

    spectrum = np.fft.fft(values, axis=1)
    filtered = 2.0 * np.real(np.fft.ifft(spectrum * fx[np.newaxis, :], axis=1))

    peak_index = _nearest_index(hz, f)
    left_index = _nearest_index(fx[: peak_index + 1], 0.5)
    right_index = peak_index + _nearest_index(fx[peak_index:], 0.5)

    impulse_envelope = np.abs(
        hilbert(np.real(np.fft.fftshift(np.fft.ifft(fx))))
    )
    impulse_envelope /= np.max(impulse_envelope)
    time = np.arange(n_time, dtype=float) / float(srate)
    temporal_peak = int(np.argmax(impulse_envelope))
    temporal_left = _nearest_index(impulse_envelope[: temporal_peak + 1], 0.5)
    temporal_right = temporal_peak + _nearest_index(
        impulse_envelope[temporal_peak:], 0.5
    )

    empirical_values = np.asarray(
        [
            hz[peak_index],
            hz[right_index] - hz[left_index],
            (time[temporal_right] - time[temporal_left]) * 1000.0,
        ],
        dtype=float,
    )

    if showplot:
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError(
                "Plotting filterFGx requires the optional 'plot' dependencies."
            ) from exc

        _, axis = plt.subplots()
        axis.plot(hz, fx, color="black", linewidth=1.3, label="Filter kernel")
        axis.set_xlim(max(f - 10.0, 0.0), f + 10.0)
        axis.set_title(
            f"Empirical filter: {empirical_values[0]:g}, "
            f"{empirical_values[1]:g} Hz"
        )
        axis.set_xlabel("Frequency (Hz)")
        axis.set_ylabel("Filter gain")
        axis.legend()

    if input_was_1d:
        filtered = filtered[0]
    return filtered, empirical_values, fx


__all__ = ["filterFGx"]

