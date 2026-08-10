"""Quadratic Rényi entropy and effective dimensionality across frequencies.

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

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]


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
            "(eigenvalues, frequencies[, participants])"
        )
    if min(values.shape) == 0:
        raise ValueError("FREQ.evals is empty")
    if np.any(np.isinf(values)):
        raise ValueError("FREQ.evals contains infinite values")
    if np.any(values < 0):
        raise ValueError("FREQ.evals cannot contain negative eigenvalues")
    return values


def _frequency_axis(FREQ: Any, nfrequencies: int) -> tuple[FloatArray, str]:
    raw = _field(FREQ, "frex", required=False)
    if raw is not None:
        frequencies = np.asarray(raw, dtype=float)
        if (
            frequencies.ndim in (1, 2)
            and (frequencies.ndim == 1 or 1 in frequencies.shape)
            and frequencies.size == nfrequencies
            and np.all(np.isfinite(frequencies))
        ):
            return frequencies.reshape(-1).copy(), "Frequency (Hz)"
    return np.arange(1, nfrequencies + 1, dtype=float), "Frequency index"


def _matlab_output(values: FloatArray) -> FloatArray | float:
    if values.size == 1:
        return float(values[0, 0])
    return values


def _set_frequency_limits(axis: Any, frequencies: FloatArray) -> None:
    if frequencies.size == 1:
        padding = max(abs(float(frequencies[0])) * 0.05, 0.5)
        axis.set_xlim(frequencies[0] - padding, frequencies[0] + padding)
    else:
        axis.set_xlim(float(np.min(frequencies)), float(np.max(frequencies)))


def _style_axis(axis: Any, frequencies: FloatArray, xlabel: str, title: str) -> None:
    axis.set_xlabel(xlabel, fontweight="bold")
    axis.set_ylabel(r"$H_2$", fontweight="bold")
    axis.set_title(title, fontweight="bold")
    axis.grid(True, which="both", alpha=0.25)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    _set_frequency_limits(axis, frequencies)


def _plot_h2(
    H2: FloatArray,
    frequencies: FloatArray,
    xlabel: str,
    *,
    plot_all: bool,
) -> list[Any]:
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "FREQNESS_EntropyLandscape plotting requires matplotlib. "
            "Install freqness[plot]."
        ) from error

    nsubjects = H2.shape[1]
    color = plt.get_cmap("viridis")(0.12)
    figures: list[Any] = []

    with np.errstate(invalid="ignore"):
        average = np.nanmean(H2, axis=1)
    if nsubjects > 1:
        sem = np.nanstd(H2, axis=1, ddof=1) / np.sqrt(nsubjects)
    else:
        sem = None

    figure, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
    if sem is not None and not np.all(np.isnan(sem)):
        axis.fill_between(
            frequencies,
            average - sem,
            average + sem,
            color="black",
            alpha=0.1,
            linewidth=0,
        )
    axis.plot(frequencies, average, color=color, linewidth=2)
    _style_axis(
        axis,
        frequencies,
        xlabel,
        "Quadratic Rényi Entropy - Grand Average",
    )
    figures.append(figure)

    if plot_all:
        for subject in range(nsubjects):
            figure, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
            axis.plot(
                frequencies,
                H2[:, subject],
                color=color,
                linewidth=1.5,
            )
            _style_axis(
                axis,
                frequencies,
                xlabel,
                f"Quadratic Rényi Entropy - Participant #{subject + 1}",
            )
            figures.append(figure)
    return figures


def FREQNESS_EntropyLandscape(
    FREQ: Any,
    *,
    plot_all: bool = False,
    show: bool = True,
) -> tuple[FloatArray | float, FloatArray | float]:
    """Compute quadratic Rényi entropy and effective dimensionality.

    The numerical outputs mirror the MATLAB function. ``FREQ.evals`` may be
    shaped as ``(eigenvalues,)``, ``(eigenvalues, frequencies)``, or
    ``(eigenvalues, frequencies, participants)``. NetworkEstimation supplies
    the complete eigenspectrum independently of ``ncomps``. The visualization
    deliberately includes only quadratic Rényi entropy (H2); effective
    dimensionality (ED) is returned but is not plotted.

    Parameters
    ----------
    FREQ
        A :class:`~freqness.FREQNESSResult`, mapping, or object exposing an
        ``evals`` field and, optionally, a ``frex`` field.
    plot_all
        Produce an additional H2 figure for every participant.
    show
        Call :func:`matplotlib.pyplot.show` after creating the figures.

    Returns
    -------
    H2, ED
        Arrays shaped ``(frequencies, participants)``. When there is exactly
        one frequency and one participant, both outputs are scalars, matching
        MATLAB behaviour.
    """
    for value, name in ((plot_all, "plot_all"), (show, "show")):
        if not isinstance(value, (bool, np.bool_)):
            raise TypeError(f"{name} must be a boolean")

    eigenspectrum = _eigenspectrum_array(FREQ)
    nfrequencies = eigenspectrum.shape[1]
    frequencies, xlabel = _frequency_axis(FREQ, nfrequencies)

    tiny = np.finfo(float).tiny
    sums = np.sum(eigenspectrum, axis=0, keepdims=True)
    probabilities = eigenspectrum / np.maximum(sums, tiny)
    squared_probability_sum = np.sum(probabilities**2, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        H2_matrix = -np.log(squared_probability_sum)

    numerator = np.sum(eigenspectrum, axis=0) ** 2
    denominator = np.maximum(np.sum(eigenspectrum**2, axis=0), tiny)
    ED_matrix = numerator / denominator

    _plot_h2(
        H2_matrix,
        frequencies,
        xlabel,
        plot_all=bool(plot_all),
    )
    if show:
        import matplotlib.pyplot as plt

        plt.show()

    return _matlab_output(H2_matrix), _matlab_output(ED_matrix)


__all__ = ["FREQNESS_EntropyLandscape"]
