"""Individual maxima of FREQ-NESS network landscapes within a frequency range.

If you use this toolbox, please cite:
Rosso, M., Fernández-Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
Kringelbach, M. L., & Bonetti, L. (2025). FREQ-NESS Reveals the Dynamic
Reconfiguration of Frequency-Resolved Brain Networks During Auditory
Stimulation. Advanced Science, 2413195.
https://doi.org/10.1002/advs.202413195

Method: individual alpha-network peak analysis described in
https://doi.org/10.64898/2026.09.07.749858

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

from dataclasses import dataclass, field
from typing import Any, Mapping, TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING:
    from matplotlib.figure import Figure


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
BoolArray = NDArray[np.bool_]


@dataclass(slots=True)
class FREQNESSIndividualPeakResult:
    """Individual maxima; numerical arrays retain (components, participants).

    ``which_comp`` is one-based. ``idx_peak`` indexes the original frequency
    axis and is zero-based for NumPy; missing indices are NaN.
    """

    peak_frex: FloatArray
    peak_value: FloatArray
    idx_peak: FloatArray
    which_comp: IntArray
    frex_range: FloatArray
    frex: FloatArray
    at_boundary: BoolArray
    is_tied: BoolArray
    figures: list[Figure] = field(default_factory=list)


def _field(container: Any, name: str) -> Any:
    if isinstance(container, Mapping) and name in container:
        return container[name]
    if hasattr(container, name):
        return getattr(container, name)
    raise ValueError(f"FREQ.{name} is required")


def FREQNESS_IndividualPeak(
    FREQ: Any,
    frex_range: ArrayLike,
    *,
    which_comp: ArrayLike = 1,
    plot_avg: bool = False,
    plot_all: bool = False,
    show: bool = True,
) -> FREQNESSIndividualPeakResult:
    """Find each participant's largest landscape value within a Hz interval.

    Parameters
    ----------
    FREQ
        Mapping or object with ``evals`` shaped (components, frequencies,
        participants), or (components, frequencies) for one participant,
        and positive, strictly increasing ``frex`` in Hz.
    frex_range
        Required [lower, upper] limits within the stored frequency axis.
        Only existing bins inside this inclusive range are searched; at
        least two stored bins must be included.
    which_comp
        One-based component index or vector of indices. Default: 1.
    plot_avg
        Plot the group distribution of individual peak frequencies, with
        participant points, a boxplot, and group mean. Default: False.
    plot_all
        Plot each participant's full selected-component landscape, search
        interval, and selected maximum. Default: False.
    show
        Display generated figures. Default: True.

    Returns
    -------
    FREQNESSIndividualPeakResult
        Peak frequencies, stored values, and indices in the original axis,
        all shaped (selected components, participants). Boundary and tied
        maxima are flagged. The requested range and searched bins are kept.

    Notes
    -----
    Values are used as supplied, without normalization, smoothing, or
    interpolation. Ties select the lowest-frequency maximum. Non-finite
    bins are ignored; all-missing windows return NaN and false flags.
    A range maximum does not establish an interior oscillatory peak, and a
    network-landscape peak is distinct from a power-spectrum peak.
    """
    raw_evals = np.asarray(_field(FREQ, "evals"))
    if (
        raw_evals.size == 0
        or raw_evals.ndim not in (2, 3)
        or not np.issubdtype(raw_evals.dtype, np.number)
        or np.iscomplexobj(raw_evals)
    ):
        raise ValueError(
            "FREQ.evals must be real numeric "
            "(components, frequencies[, participants])"
        )
    evals = np.asarray(raw_evals, dtype=float)
    if evals.ndim == 2:
        evals = evals[:, :, np.newaxis]

    raw_frex = np.asarray(_field(FREQ, "frex"))
    if (
        not np.issubdtype(raw_frex.dtype, np.number)
        or np.iscomplexobj(raw_frex)
        or raw_frex.ndim not in (1, 2)
        or (raw_frex.ndim == 2 and 1 not in raw_frex.shape)
    ):
        raise ValueError("FREQ.frex must be a real numeric frequency vector")
    frex = np.asarray(raw_frex, dtype=float).reshape(-1)
    if (
        frex.size != evals.shape[1]
        or not np.all(np.isfinite(frex))
        or np.any(frex <= 0)
        or np.any(np.diff(frex) <= 0)
    ):
        raise ValueError(
            "FREQ.frex must contain positive, strictly increasing frequencies "
            "matching the frequency dimension of FREQ.evals"
        )

    raw_range = np.asarray(frex_range)
    if (
        raw_range.size != 2
        or raw_range.ndim not in (1, 2)
        or (raw_range.ndim == 2 and 1 not in raw_range.shape)
        or not np.issubdtype(raw_range.dtype, np.number)
        or np.iscomplexobj(raw_range)
    ):
        raise ValueError("frex_range must contain two finite ascending limits")
    search_range = np.asarray(raw_range, dtype=float).reshape(-1)
    if not np.all(np.isfinite(search_range)) or search_range[0] >= search_range[1]:
        raise ValueError("frex_range must contain two finite ascending limits")
    if search_range[0] < frex[0] or search_range[1] > frex[-1]:
        raise ValueError("frex_range must lie within the analyzed FREQ.frex range")
    idx_frex = np.flatnonzero((frex >= search_range[0]) & (frex <= search_range[1]))
    if idx_frex.size < 2:
        raise ValueError("frex_range must include at least two stored frequency bins")

    raw_comps = np.asarray(which_comp)
    if (
        raw_comps.size == 0
        or raw_comps.ndim > 2
        or (raw_comps.ndim == 2 and 1 not in raw_comps.shape)
        or not np.issubdtype(raw_comps.dtype, np.number)
        or np.iscomplexobj(raw_comps)
        or not np.all(np.isfinite(raw_comps))
        or not np.all(raw_comps == np.round(raw_comps))
        or np.any(raw_comps < 1)
        or np.any(raw_comps > evals.shape[0])
    ):
        raise ValueError("which_comp must contain valid one-based component indices")
    comps = np.asarray(list(dict.fromkeys(raw_comps.reshape(-1))), dtype=np.int64)
    for value, name in ((plot_avg, "plot_avg"), (plot_all, "plot_all"), (show, "show")):
        if not isinstance(value, (bool, np.bool_)):
            raise TypeError(f"{name} must be a boolean")

    ncomps = comps.size
    nsubs = evals.shape[2]
    PEAK = FREQNESSIndividualPeakResult(
        peak_frex=np.full((ncomps, nsubs), np.nan),
        peak_value=np.full((ncomps, nsubs), np.nan),
        idx_peak=np.full((ncomps, nsubs), np.nan),
        which_comp=comps,
        frex_range=search_range.copy(),
        frex=frex[idx_frex].copy(),
        at_boundary=np.zeros((ncomps, nsubs), dtype=bool),
        is_tied=np.zeros((ncomps, nsubs), dtype=bool),
    )
    print(
        f"\nFREQNESS Individual Peak: searching {search_range[0]:.4g}-"
        f"{search_range[1]:.4g} Hz ({idx_frex.size} bins) for "
        f"{ncomps} component(s) and {nsubs} participant(s)."
    )
    for compi, component in enumerate(comps):
        for subi in range(nsubs):
            this_landscape = evals[component - 1, idx_frex, subi]
            mask_valid = np.isfinite(this_landscape)
            if not np.any(mask_valid):
                continue
            this_landscape = np.where(mask_valid, this_landscape, -np.inf)
            idx_local = int(np.argmax(this_landscape))
            peak_value = this_landscape[idx_local]
            idx_peak = idx_frex[idx_local]
            PEAK.peak_frex[compi, subi] = frex[idx_peak]
            PEAK.peak_value[compi, subi] = peak_value
            PEAK.idx_peak[compi, subi] = idx_peak
            PEAK.at_boundary[compi, subi] = idx_local in (0, idx_frex.size - 1)
            PEAK.is_tied[compi, subi] = np.count_nonzero(this_landscape == peak_value) > 1

    if plot_avg or plot_all:
        PEAK.figures.extend(_plot_peak_results(evals, frex, PEAK, plot_avg, plot_all))
        if show:
            import matplotlib.pyplot as plt

            plt.show()
    return PEAK


def _plot_peak_results(
    evals: FloatArray,
    frex: FloatArray,
    PEAK: FREQNESSIndividualPeakResult,
    plot_avg: bool,
    plot_all: bool,
) -> list[Figure]:
    """Plot individual estimates without modifying them or the global RNG."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "FREQNESS_IndividualPeak plotting requires matplotlib. "
            "Install freqness[plot]."
        ) from error

    col_line = "#0072B2"
    col_peak = "#D55E00"
    figures: list[Figure] = []
    for compi, component in enumerate(PEAK.which_comp):
        if plot_avg:
            figure, axis = plt.subplots(figsize=(4.8, 5.2), constrained_layout=True)
            this_peaks = PEAK.peak_frex[compi]
            this_peaks = this_peaks[np.isfinite(this_peaks)]
            if this_peaks.size:
                # Match MATLAB boxchart's midpoint-CDF quartile convention.
                lower, median, upper = np.percentile(
                    this_peaks, [25, 50, 75], method="hazen"
                )
                spread = 1.5 * (upper - lower)
                axis.bxp(
                    [{
                        "q1": lower, "med": median, "q3": upper,
                        "whislo": np.min(this_peaks[this_peaks >= lower - spread]),
                        "whishi": np.max(this_peaks[this_peaks <= upper + spread]),
                    }],
                    positions=[1],
                    widths=0.35,
                    showfliers=False,
                    patch_artist=True,
                    boxprops={"facecolor": col_line, "alpha": 0.15},
                    medianprops={"color": col_line, "linewidth": 1.8},
                )
                offsets = 0.12 * np.sin(np.arange(1, this_peaks.size + 1) * 2.399963)
                axis.scatter(
                    1 + offsets, this_peaks, s=40, color=col_line,
                    alpha=0.75, edgecolor="white", zorder=3,
                )
                axis.scatter(
                    1, np.mean(this_peaks), s=90, color=col_peak,
                    marker="D", edgecolor="white", zorder=4,
                )
            else:
                axis.text(1, np.mean(PEAK.frex_range), "No valid peaks", ha="center")
            padding = 0.08 * np.diff(PEAK.frex_range)[0]
            axis.set_ylim(PEAK.frex_range + [-padding, padding])
            axis.set_xlim(0.6, 1.4)
            axis.set_xticks([1], ["Participants"])
            axis.set_ylabel("Individual peak frequency (Hz)", fontweight="bold")
            axis.set_title(
                f"Individual Peak - Component {component}\n"
                f"n = {this_peaks.size} | Diamond: group mean"
            )
            figures.append(figure)
        if plot_all:
            for subi in range(evals.shape[2]):
                figure, axis = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
                this_landscape = evals[component - 1, :, subi]
                axis.plot(
                    frex, np.where(np.isfinite(this_landscape), this_landscape, np.nan),
                    color=col_line, linewidth=2,
                )
                axis.axvspan(*PEAK.frex_range, color=col_line, alpha=0.07)
                for boundary in PEAK.frex_range:
                    axis.axvline(boundary, linestyle="--", color="0.5", linewidth=1)
                if np.isfinite(PEAK.peak_frex[compi, subi]):
                    axis.scatter(
                        PEAK.peak_frex[compi, subi], PEAK.peak_value[compi, subi],
                        s=70, color=col_peak, edgecolor="white", zorder=4,
                    )
                    subtitle = f"Selected maximum: {PEAK.peak_frex[compi, subi]:.4g} Hz"
                else:
                    subtitle = "No valid values in the search interval"
                axis.set_xlim(frex[0], frex[-1])
                axis.set_xlabel("Frequency (Hz)", fontweight="bold")
                axis.set_ylabel("Stored landscape value", fontweight="bold")
                axis.set_title(
                    f"Individual Peak - Participant {subi + 1} - "
                    f"Component {component}\n{subtitle}"
                )
                figures.append(figure)
    for figure in figures:
        axis = figure.axes[0]
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.15)
        axis.set_axisbelow(True)
    return figures


individual_peak = FREQNESS_IndividualPeak

__all__ = ["FREQNESS_IndividualPeak", "FREQNESSIndividualPeakResult", "individual_peak"]
