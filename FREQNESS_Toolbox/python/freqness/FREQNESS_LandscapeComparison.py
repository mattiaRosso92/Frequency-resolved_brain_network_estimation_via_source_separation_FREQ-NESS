"""Pairwise or factorial comparison of participant-level landscapes.

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

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence, TYPE_CHECKING
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import stats

from ._factorial_utils import (
    anova_statistics,
    effect_definitions,
    is_factorial_grid,
    long_data,
    match_subjects,
    normalize_grid,
    permutation_anova_p_values,
    permutation_count,
    rank_statistics,
    validate_design,
    validate_factor_names,
    validate_level_labels,
    validate_test,
)

if TYPE_CHECKING:
    from matplotlib.figure import Figure


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
BoolArray = NDArray[np.bool_]

_AUTO_MEASURES = (
    "evals",
    "H2",
    "ED",
    "amplitude_raw",
    "amplitude_normalized",
    "stability_index",
    "mean_frequency",
)


@dataclass(slots=True)
class FREQNESSLandscapeMeasureResult:
    """Statistics for one landscape field."""

    measure: str
    requested_frex: FloatArray
    frex: FloatArray
    idx_frex: IntArray
    comps: IntArray
    values_A: FloatArray
    values_B: FloatArray
    mean_A: FloatArray
    mean_B: FloatArray
    median_A: FloatArray
    median_B: FloatArray
    std_A: FloatArray
    std_B: FloatArray
    sem_A: FloatArray
    sem_B: FloatArray
    mean_difference: FloatArray
    median_difference: FloatArray
    statistic: FloatArray
    statistic_name: str
    effect_size: FloatArray
    effect_size_name: str
    p_uncorrected: FloatArray
    q_fdr: FloatArray
    significant: BoolArray
    fdr_threshold: FloatArray
    n_A: IntArray
    n_B: IntArray
    figures: list[Figure]


@dataclass
class FREQNESSLandscapeComparisonResult:
    """Outputs of :func:`FREQNESS_LandscapeComparison`.

    Individual measure results are available both through ``results`` and as
    attributes preserving the input field name, for example ``STATS.evals``.
    """

    analysis: str
    design: str
    test: str
    sample_labels: tuple[str, str]
    alpha: float
    correction: str
    correction_scope: str
    measures: tuple[str, ...]
    results: dict[str, FREQNESSLandscapeMeasureResult]
    figures: list[Figure] = field(default_factory=list)

    def __getattr__(self, name: str) -> FREQNESSLandscapeMeasureResult:
        results = self.__dict__.get("results", {})
        if name in results:
            return results[name]
        raise AttributeError(name)


@dataclass(slots=True)
class FREQNESSLandscapeFactorialMeasureResult:
    """Frequency-wise statistics for one measure and factorial effect."""

    measure: str
    name: str
    key: str
    factors: tuple[int, ...]
    requested_frex: FloatArray
    frex: FloatArray
    idx_frex: IntArray
    comps: IntArray
    values: NDArray[np.object_]
    statistic: FloatArray
    statistic_name: str
    effect_size: FloatArray
    effect_size_name: str
    p_uncorrected: FloatArray
    q_fdr: FloatArray
    significant: BoolArray
    fdr_threshold: FloatArray
    df1: FloatArray
    df2: FloatArray
    figures: list[Figure] = field(default_factory=list)


@dataclass
class FREQNESSLandscapeFactorialEffect:
    """All analyzed landscape measures for one factorial effect."""

    name: str
    key: str
    factors: tuple[int, ...]
    results: dict[str, FREQNESSLandscapeFactorialMeasureResult]

    def __getattr__(self, name: str) -> FREQNESSLandscapeFactorialMeasureResult:
        results = self.__dict__.get("results", {})
        if name in results:
            return results[name]
        raise AttributeError(name)


@dataclass(slots=True)
class FREQNESSLandscapeFactorialComparisonResult:
    """Outputs for a generalized factorial landscape comparison."""

    analysis: str
    factor_names: tuple[str, ...]
    level_labels: tuple[tuple[str, ...], ...]
    factor_design: tuple[str, ...]
    design: tuple[str, ...]
    input_shape: tuple[int, ...]
    test: str
    test_requested: str
    tail: str
    alpha: float
    correction: str
    correction_scope: str
    measures: tuple[str, ...]
    effects: dict[str, FREQNESSLandscapeFactorialEffect]
    effect_names: tuple[str, ...]
    cell_means: dict[str, NDArray[np.object_]]
    cell_n_participants: IntArray
    subject_ids: NDArray[np.object_]
    n_permutations: int
    random_seed: int
    figures: list[Figure] = field(default_factory=list)


def _field_name(container: Any, requested: str) -> str | None:
    if isinstance(container, Mapping):
        names = tuple(str(name) for name in container.keys())
    elif isinstance(container, np.ndarray):
        return None
    else:
        names = tuple(name for name in dir(container) if not name.startswith("_"))
    matches = [name for name in names if name.lower() == requested.lower()]
    return matches[0] if matches else None


def _field(container: Any, requested: str) -> Any:
    name = _field_name(container, requested)
    if name is None:
        raise ValueError(f"Input does not contain the requested field '{requested}'")
    if isinstance(container, Mapping):
        return container[name]
    return getattr(container, name)


def _is_numeric_input(container: Any) -> bool:
    """Return true for raw numeric array-like inputs, but not structures."""
    if isinstance(container, Mapping) or hasattr(container, "__dict__"):
        return False
    try:
        values = np.asarray(container)
    except (TypeError, ValueError):
        return False
    return values.ndim > 0 and np.issubdtype(values.dtype, np.number)


def _measure_selection(
    landscape_A: Any,
    landscape_B: Any,
    measure: str | Sequence[str] | None,
) -> tuple[str, ...]:
    numeric_A = _is_numeric_input(landscape_A)
    numeric_B = _is_numeric_input(landscape_B)
    if numeric_A != numeric_B:
        raise TypeError("LANDSCAPE_A and LANDSCAPE_B must use the same input format")
    if numeric_A:
        if measure is None:
            return ("values",)
        requested = (measure,) if isinstance(measure, str) else tuple(measure)
        if len(requested) != 1 or requested[0].lower() != "values":
            raise ValueError("Raw numeric inputs support only measure='values'")
        return ("values",)

    if measure is None:
        selected = tuple(
            name
            for name in _AUTO_MEASURES
            if _field_name(landscape_A, name) is not None
            and _field_name(landscape_B, name) is not None
        )
        if not selected:
            raise ValueError(
                "No shared landscape fields were recognized; supply measure explicitly"
            )
        return selected

    requested = (measure,) if isinstance(measure, str) else tuple(measure)
    if not requested or not all(isinstance(name, str) and name for name in requested):
        raise ValueError("measure must be text or a non-empty sequence of field names")
    selected: list[str] = []
    for name in requested:
        actual_A = _field_name(landscape_A, name)
        actual_B = _field_name(landscape_B, name)
        if actual_A is None or actual_B is None:
            raise ValueError(f"Both inputs must contain the measure '{name}'")
        if actual_A.lower() != actual_B.lower():
            raise ValueError(f"Measure names do not match for '{name}'")
        if actual_A not in selected:
            selected.append(actual_A)
    return tuple(selected)


def _vector(values: Any, name: str) -> FloatArray:
    vector = np.asarray(values, dtype=float)
    if vector.ndim not in (1, 2) or (vector.ndim == 2 and 1 not in vector.shape):
        raise ValueError(f"{name} must be a vector")
    vector = vector.reshape(-1)
    if vector.size == 0 or not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain finite values")
    if np.any(np.diff(vector) <= 0):
        raise ValueError(f"{name} must be strictly increasing")
    return vector


def _frequency_axis(
    container: Any,
    measure: str,
    supplied_frex: ArrayLike | None,
) -> FloatArray:
    if supplied_frex is not None:
        return _vector(supplied_frex, "frex")
    if _is_numeric_input(container):
        raise ValueError("frex is required for raw numeric landscape inputs")
    axis_name = (
        "carrier_frex"
        if measure.lower() in {"amplitude_raw", "amplitude_normalized"}
        else "frex"
    )
    return _vector(_field(container, axis_name), f"LANDSCAPE.{axis_name}")


def _canonical_values(
    container: Any,
    measure: str,
    frequencies: FloatArray,
    label: str,
) -> FloatArray:
    raw = container if measure == "values" else _field(container, measure)
    values = np.asarray(raw)
    if not np.issubdtype(values.dtype, np.number) or np.iscomplexobj(values):
        raise TypeError(f"{label}.{measure} must contain real numeric values")
    values = np.asarray(values, dtype=float)
    if np.any(np.isinf(values)):
        raise ValueError(f"{label}.{measure} contains infinite values")
    if values.ndim == 1:
        if values.size != frequencies.size:
            raise ValueError(f"{label}.{measure} does not match its frequency axis")
        values = values[np.newaxis, :, np.newaxis]
    elif values.ndim == 2:
        if values.shape[0] != frequencies.size:
            raise ValueError(
                f"{label}.{measure} must be [nFrex x nParticipants]"
            )
        values = values[np.newaxis, :, :]
    elif values.ndim == 3:
        if values.shape[1] != frequencies.size:
            raise ValueError(
                f"{label}.{measure} must be "
                "[nFeatures x nFrex x nParticipants]"
            )
    else:
        raise ValueError(
            f"{label}.{measure} must be one-, two-, or three-dimensional"
        )
    if min(values.shape) == 0:
        raise ValueError(f"{label}.{measure} is empty")
    return values


def _component_selection(which_comp: ArrayLike | None, nfeatures: int) -> IntArray:
    requested = np.array([1]) if which_comp is None else np.asarray(which_comp)
    if requested.ndim == 0:
        requested = requested.reshape(1)
    if requested.ndim not in (1, 2) or (
        requested.ndim == 2 and 1 not in requested.shape
    ):
        raise ValueError("which_comp must be a scalar or vector")
    requested = requested.reshape(-1)
    if (
        requested.size == 0
        or not np.issubdtype(requested.dtype, np.number)
        or np.iscomplexobj(requested)
        or not np.all(np.isfinite(requested))
        or not np.all(requested == np.round(requested))
    ):
        raise ValueError("which_comp must contain finite integer component numbers")
    components = requested.astype(np.int64)
    if np.any(components < 1) or np.any(components > nfeatures):
        raise ValueError(f"which_comp must lie between 1 and {nfeatures}")
    if np.unique(components).size != components.size:
        raise ValueError("which_comp cannot contain duplicates")
    return components


def _frequency_selection(
    frequencies: FloatArray,
    frex2compare: ArrayLike | str | None,
    tolerance: float | None,
) -> tuple[FloatArray, FloatArray, IntArray]:
    if frex2compare is None or (
        isinstance(frex2compare, str) and frex2compare.lower() == "all"
    ):
        indices = np.arange(frequencies.size, dtype=np.int64)
        return frequencies.copy(), frequencies.copy(), indices
    if isinstance(frex2compare, str):
        raise ValueError("frex2compare must be numeric, 'all', or None")
    requested = np.asarray(frex2compare, dtype=float).reshape(-1)
    if requested.size == 0 or not np.all(np.isfinite(requested)):
        raise ValueError("frex2compare must contain finite frequencies")
    if tolerance is None:
        if frequencies.size == 1:
            tolerance = max(abs(float(frequencies[0])) * 1e-10, 1e-10)
        else:
            tolerance = float(np.median(np.diff(frequencies)) / 2.0)
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("frex_tolerance must be a non-negative finite scalar")
    indices = np.array(
        [int(np.argmin(np.abs(frequencies - value))) for value in requested],
        dtype=np.int64,
    )
    distances = np.abs(frequencies[indices] - requested)
    if np.any(distances > tolerance):
        raise ValueError("One or more requested frequencies exceed frex_tolerance")
    if np.unique(indices).size != indices.size:
        raise ValueError("frex2compare maps more than one request to the same bin")
    return requested, frequencies[indices], indices


def _validate_labels(sample_labels: Sequence[str]) -> tuple[str, str]:
    if isinstance(sample_labels, str):
        raise ValueError("sample_labels must contain two non-empty strings")
    try:
        labels = tuple(sample_labels)
    except TypeError as error:
        raise ValueError(
            "sample_labels must contain two non-empty strings"
        ) from error
    if len(labels) != 2 or not all(
        isinstance(label, str) and label for label in labels
    ):
        raise ValueError("sample_labels must contain two non-empty strings")
    return labels[0], labels[1]


def _analysis_options(
    design: str,
    test: str | None,
    tail: str,
) -> tuple[str, str, str]:
    design_key = str(design).lower()
    if design_key in {"within", "paired"}:
        design_key = "within"
    elif design_key in {"between", "independent"}:
        design_key = "between"
    else:
        raise ValueError("design must be 'within' or 'between'")

    if test is None:
        test_key = "signrank" if design_key == "within" else "ranksum"
    else:
        test_key = str(test).lower()
    if test_key not in {"ttest", "signrank", "ranksum"}:
        raise ValueError("test must be 'ttest', 'signrank', or 'ranksum'")
    if design_key == "within" and test_key == "ranksum":
        warnings.warn(
            "test='ranksum' is an independent-samples test and is incompatible "
            "with design='within'. Using test='signrank' instead.",
            UserWarning,
            stacklevel=3,
        )
        test_key = "signrank"
    elif design_key == "between" and test_key == "signrank":
        warnings.warn(
            "test='signrank' is a paired-samples test and is incompatible with "
            "design='between'. Using test='ranksum' instead.",
            UserWarning,
            stacklevel=3,
        )
        test_key = "ranksum"

    tail_key = str(tail).lower()
    if tail_key not in {"both", "left", "right"}:
        raise ValueError("tail must be 'both', 'left', or 'right'")
    return design_key, test_key, tail_key


def _tail_name(tail: str) -> str:
    return {"both": "two-sided", "left": "less", "right": "greater"}[tail]


def _paired_test(
    values_A: FloatArray,
    values_B: FloatArray,
    test: str,
    tail: str,
) -> tuple[float, str, float, str, float]:
    differences = values_A - values_B
    if np.all(differences == 0):
        if test == "ttest":
            return 0.0, "t", 0.0, "Cohen's dz", 1.0
        return 0.0, "W+", 0.0, "matched rank-biserial r", 1.0
    if test == "ttest":
        result = stats.ttest_rel(values_A, values_B, alternative=_tail_name(tail))
        deviation = float(np.std(differences, ddof=1))
        effect = (
            float(np.mean(differences) / deviation)
            if deviation > 0
            else np.sign(np.mean(differences)) * np.inf
        )
        return float(result.statistic), "t", effect, "Cohen's dz", float(result.pvalue)

    nonzero = differences[differences != 0]
    ranks = stats.rankdata(np.abs(nonzero), method="average")
    w_positive = float(np.sum(ranks[nonzero > 0]))
    w_negative = float(np.sum(ranks[nonzero < 0]))
    denominator = w_positive + w_negative
    effect = (w_positive - w_negative) / denominator if denominator else 0.0
    result = stats.wilcoxon(
        values_A,
        values_B,
        zero_method="wilcox",
        correction=False,
        alternative=_tail_name(tail),
        method="approx",
    )
    return (
        w_positive,
        "W+",
        float(effect),
        "matched rank-biserial r",
        float(result.pvalue),
    )


def _independent_test(
    values_A: FloatArray,
    values_B: FloatArray,
    test: str,
    tail: str,
) -> tuple[float, str, float, str, float]:
    if np.all(values_A == values_A[0]) and np.all(values_B == values_A[0]):
        if test == "ttest":
            return 0.0, "t", 0.0, "Hedges' g", 1.0
        return (
            values_A.size * values_B.size / 2.0,
            "U",
            0.0,
            "rank-biserial r",
            1.0,
        )
    if test == "ttest":
        result = stats.ttest_ind(
            values_A,
            values_B,
            equal_var=False,
            alternative=_tail_name(tail),
        )
        n_A, n_B = values_A.size, values_B.size
        variance_A = float(np.var(values_A, ddof=1))
        variance_B = float(np.var(values_B, ddof=1))
        pooled = np.sqrt(
            ((n_A - 1) * variance_A + (n_B - 1) * variance_B)
            / (n_A + n_B - 2)
        )
        if pooled > 0:
            cohens_d = float((np.mean(values_A) - np.mean(values_B)) / pooled)
            correction = 1.0 - 3.0 / (4.0 * (n_A + n_B) - 9.0)
            effect = correction * cohens_d
        else:
            mean_difference = float(np.mean(values_A) - np.mean(values_B))
            effect = 0.0 if mean_difference == 0 else np.sign(mean_difference) * np.inf
        return float(result.statistic), "t", effect, "Hedges' g", float(result.pvalue)

    result = stats.mannwhitneyu(
        values_A,
        values_B,
        alternative=_tail_name(tail),
        method="asymptotic",
    )
    effect = 2.0 * float(result.statistic) / (values_A.size * values_B.size) - 1.0
    return float(result.statistic), "U", effect, "rank-biserial r", float(result.pvalue)


def _bh_fdr(p_values: FloatArray, alpha: float) -> tuple[FloatArray, float]:
    q_values = np.full(p_values.shape, np.nan)
    valid = np.isfinite(p_values)
    values = p_values[valid]
    if values.size == 0:
        return q_values, np.nan
    order = np.argsort(values)
    sorted_p = values[order]
    ranks = np.arange(1, sorted_p.size + 1, dtype=float)
    sorted_q = sorted_p * sorted_p.size / ranks
    sorted_q = np.minimum.accumulate(sorted_q[::-1])[::-1]
    sorted_q = np.minimum(sorted_q, 1.0)
    restored = np.empty_like(sorted_q)
    restored[order] = sorted_q
    q_values[valid] = restored
    accepted = sorted_p <= alpha * ranks / sorted_p.size
    threshold = (
        float(sorted_p[np.flatnonzero(accepted)[-1]])
        if np.any(accepted)
        else np.nan
    )
    return q_values, threshold


def _empty_matrix(ncomponents: int, nfrequencies: int) -> FloatArray:
    return np.full((ncomponents, nfrequencies), np.nan)


def _compare_measure(
    measure: str,
    values_A: FloatArray,
    values_B: FloatArray,
    requested_frex: FloatArray,
    resolved_frex: FloatArray,
    frequency_indices: IntArray,
    components: IntArray,
    design: str,
    test: str,
    tail: str,
    alpha: float,
    correction: str,
) -> FREQNESSLandscapeMeasureResult:
    selected_A = values_A[components - 1][:, frequency_indices, :]
    selected_B = values_B[components - 1][:, frequency_indices, :]
    ncomponents, nfrequencies = selected_A.shape[:2]
    matrices = [_empty_matrix(ncomponents, nfrequencies) for _ in range(16)]
    (
        mean_A,
        mean_B,
        median_A,
        median_B,
        std_A,
        std_B,
        sem_A,
        sem_B,
        mean_difference,
        median_difference,
        statistic,
        effect_size,
        p_uncorrected,
        q_fdr,
        n_A_float,
        n_B_float,
    ) = matrices
    statistic_name = ""
    effect_size_name = ""

    for component in range(ncomponents):
        for frequency in range(nfrequencies):
            sample_A = selected_A[component, frequency]
            sample_B = selected_B[component, frequency]
            if design == "within":
                valid = np.isfinite(sample_A) & np.isfinite(sample_B)
                sample_A = sample_A[valid]
                sample_B = sample_B[valid]
            else:
                sample_A = sample_A[np.isfinite(sample_A)]
                sample_B = sample_B[np.isfinite(sample_B)]
            n_A_float[component, frequency] = sample_A.size
            n_B_float[component, frequency] = sample_B.size
            if sample_A.size:
                mean_A[component, frequency] = np.mean(sample_A)
                median_A[component, frequency] = np.median(sample_A)
                if sample_A.size > 1:
                    std_A[component, frequency] = np.std(sample_A, ddof=1)
                    sem_A[component, frequency] = (
                        std_A[component, frequency] / np.sqrt(sample_A.size)
                    )
            if sample_B.size:
                mean_B[component, frequency] = np.mean(sample_B)
                median_B[component, frequency] = np.median(sample_B)
                if sample_B.size > 1:
                    std_B[component, frequency] = np.std(sample_B, ddof=1)
                    sem_B[component, frequency] = (
                        std_B[component, frequency] / np.sqrt(sample_B.size)
                    )
            if sample_A.size and sample_B.size:
                mean_difference[component, frequency] = np.mean(sample_A) - np.mean(
                    sample_B
                )
                median_difference[component, frequency] = np.median(
                    sample_A
                ) - np.median(sample_B)
            if sample_A.size < 2 or sample_B.size < 2:
                continue
            if design == "within":
                output = _paired_test(sample_A, sample_B, test, tail)
            else:
                output = _independent_test(sample_A, sample_B, test, tail)
            (
                statistic[component, frequency],
                statistic_name,
                effect_size[component, frequency],
                effect_size_name,
                p_uncorrected[component, frequency],
            ) = output

    fdr_threshold = np.full(ncomponents, np.nan)
    if correction == "fdr":
        for component in range(ncomponents):
            q_fdr[component], fdr_threshold[component] = _bh_fdr(
                p_uncorrected[component], alpha
            )
    else:
        q_fdr[:] = p_uncorrected
        fdr_threshold[:] = alpha

    return FREQNESSLandscapeMeasureResult(
        measure=measure,
        requested_frex=requested_frex,
        frex=resolved_frex,
        idx_frex=frequency_indices,
        comps=components,
        values_A=selected_A,
        values_B=selected_B,
        mean_A=mean_A,
        mean_B=mean_B,
        median_A=median_A,
        median_B=median_B,
        std_A=std_A,
        std_B=std_B,
        sem_A=sem_A,
        sem_B=sem_B,
        mean_difference=mean_difference,
        median_difference=median_difference,
        statistic=statistic,
        statistic_name=statistic_name,
        effect_size=effect_size,
        effect_size_name=effect_size_name,
        p_uncorrected=p_uncorrected,
        q_fdr=q_fdr,
        significant=np.isfinite(q_fdr) & (q_fdr <= alpha),
        fdr_threshold=fdr_threshold,
        n_A=n_A_float.astype(np.int64),
        n_B=n_B_float.astype(np.int64),
        figures=[],
    )


def _plot_measure(
    result: FREQNESSLandscapeMeasureResult,
    labels: tuple[str, str],
    correction: str,
    alpha: float,
) -> list[Figure]:
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "FREQNESS_LandscapeComparison plotting requires matplotlib. "
            "Install freqness[plot]."
        ) from error

    figures: list[Figure] = []
    colors = plt.get_cmap("viridis")([0.15, 0.75])
    for component_index, component in enumerate(result.comps):
        figure, axes = plt.subplots(
            2,
            1,
            figsize=(9, 7),
            sharex=True,
            constrained_layout=True,
            gridspec_kw={"height_ratios": [2.0, 1.0]},
        )
        for mean, sem, label, color in (
            (
                result.mean_A[component_index],
                result.sem_A[component_index],
                labels[0],
                colors[0],
            ),
            (
                result.mean_B[component_index],
                result.sem_B[component_index],
                labels[1],
                colors[1],
            ),
        ):
            axes[0].plot(result.frex, mean, linewidth=2, color=color, label=label)
            if np.any(np.isfinite(sem)):
                axes[0].fill_between(
                    result.frex,
                    mean - sem,
                    mean + sem,
                    color=color,
                    alpha=0.18,
                    linewidth=0,
                )
        axes[0].set_ylabel(result.measure, fontweight="bold")
        axes[0].set_title(
            f"{result.measure} landscape comparison - component #{component}",
            fontweight="bold",
        )
        axes[0].legend(frameon=False)
        axes[0].grid(True, alpha=0.25)
        axes[0].spines[["top", "right"]].set_visible(False)

        difference = result.mean_difference[component_index]
        axes[1].axhline(0.0, color="black", linewidth=1, linestyle="--")
        axes[1].plot(result.frex, difference, color="black", linewidth=1.8)
        significant = result.significant[component_index]
        if np.any(significant):
            axes[1].scatter(
                result.frex[significant],
                difference[significant],
                s=42,
                color=plt.get_cmap("viridis")(0.95),
                zorder=3,
                label=(
                    f"FDR q ≤ {alpha:g}"
                    if correction == "fdr"
                    else f"p ≤ {alpha:g}"
                ),
            )
            axes[1].legend(frameon=False)
        axes[1].set_xlabel("Frequency (Hz)", fontweight="bold")
        axes[1].set_ylabel(f"{labels[0]} - {labels[1]}", fontweight="bold")
        axes[1].grid(True, alpha=0.25)
        axes[1].spines[["top", "right"]].set_visible(False)
        figures.append(figure)
    return figures


def _pairwise_landscape_comparison(
    LANDSCAPE_A: Any,
    LANDSCAPE_B: Any,
    *,
    measure: str | Sequence[str] | None = None,
    design: str = "within",
    test: str | None = None,
    which_comp: ArrayLike | None = None,
    frex: ArrayLike | None = None,
    frex2compare: ArrayLike | str | None = None,
    sample_labels: tuple[str, str] = ("A", "B"),
    tail: str = "both",
    alpha: float = 0.05,
    correction: str = "fdr",
    frex_tolerance: float | None = None,
    plot_results: bool = True,
    show: bool = True,
) -> FREQNESSLandscapeComparisonResult:
    """Compare two within-subject or between-group frequency landscapes.

    Component numbers use the MATLAB-compatible one-based convention. When
    ``measure`` is omitted, every shared recognized landscape field is tested
    separately and exposed with its original field name under ``STATS``.
    """
    design_key, test_key, tail_key = _analysis_options(design, test, tail)
    labels = _validate_labels(sample_labels)
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between 0 and 1")
    correction_key = str(correction).lower()
    if correction_key not in {"fdr", "none"}:
        raise ValueError("correction must be 'fdr' or 'none'")
    if not isinstance(plot_results, (bool, np.bool_)) or not isinstance(
        show, (bool, np.bool_)
    ):
        raise TypeError("plot_results and show must be Boolean")

    measures = _measure_selection(LANDSCAPE_A, LANDSCAPE_B, measure)
    results: dict[str, FREQNESSLandscapeMeasureResult] = {}
    all_figures: list[Figure] = []
    print(
        "\nFREQNESS Landscape Comparison: comparing "
        f"{labels[0]} and {labels[1]} using a {design_key}-subject "
        f"{test_key} analysis."
    )
    for measure_name in measures:
        frequency_axis_A = _frequency_axis(LANDSCAPE_A, measure_name, frex)
        frequency_axis_B = _frequency_axis(LANDSCAPE_B, measure_name, frex)
        if not np.allclose(frequency_axis_A, frequency_axis_B, atol=1e-10, rtol=0.0):
            raise ValueError(
                f"Frequency grids do not match for measure '{measure_name}'"
            )
        values_A = _canonical_values(
            LANDSCAPE_A, measure_name, frequency_axis_A, "LANDSCAPE_A"
        )
        values_B = _canonical_values(
            LANDSCAPE_B, measure_name, frequency_axis_B, "LANDSCAPE_B"
        )
        if values_A.shape[0] != values_B.shape[0]:
            raise ValueError(
                f"Feature counts do not match for measure '{measure_name}'"
            )
        if design_key == "within" and values_A.shape[2] != values_B.shape[2]:
            raise ValueError(
                "Within-subject comparisons require matching participant counts"
            )
        if values_A.shape[2] < 2 or values_B.shape[2] < 2:
            raise ValueError("Each sample must contain at least two participants")
        components = _component_selection(which_comp, values_A.shape[0])
        requested, resolved, indices = _frequency_selection(
            frequency_axis_A, frex2compare, frex_tolerance
        )
        print(
            f"  {measure_name}: {resolved.size} frequencies, "
            f"component(s) {components.tolist()}."
        )
        result = _compare_measure(
            measure_name,
            values_A,
            values_B,
            requested,
            resolved,
            indices,
            components,
            design_key,
            test_key,
            tail_key,
            float(alpha),
            correction_key,
        )
        if plot_results:
            result.figures.extend(
                _plot_measure(result, labels, correction_key, float(alpha))
            )
            all_figures.extend(result.figures)
        results[measure_name] = result

    if all_figures and show:
        import matplotlib.pyplot as plt

        plt.show()

    return FREQNESSLandscapeComparisonResult(
        analysis="FREQ-NESS landscape comparison",
        design=design_key,
        test=test_key,
        sample_labels=labels,
        alpha=float(alpha),
        correction=correction_key,
        correction_scope="frequencies separately for each measure and component",
        measures=measures,
        results=results,
        figures=all_figures,
    )


def _factorial_measure_selection(
    grid: NDArray[np.object_],
    measure: str | Sequence[str] | None,
) -> tuple[str, ...]:
    """Select landscape fields shared by every factorial cell."""
    numeric = tuple(_is_numeric_input(item) for item in grid.flat)
    if any(numeric) and not all(numeric):
        raise TypeError("All factorial cells must use the same input format")
    if all(numeric):
        requested = () if measure is None else (
            (measure,) if isinstance(measure, str) else tuple(measure)
        )
        if requested and (
            len(requested) != 1 or str(requested[0]).lower() != "values"
        ):
            raise ValueError("Raw numeric factorial cells support measure='values'")
        return ("values",)
    if measure is None:
        selected = tuple(
            name
            for name in _AUTO_MEASURES
            if all(_field_name(item, name) is not None for item in grid.flat)
        )
        if not selected:
            raise ValueError(
                "No shared landscape fields were recognized; supply measure explicitly"
            )
        return selected
    requested = (measure,) if isinstance(measure, str) else tuple(measure)
    if not requested or not all(isinstance(name, str) and name for name in requested):
        raise ValueError("measure must be text or a non-empty sequence")
    selected: list[str] = []
    for requested_name in requested:
        actual = _field_name(grid.flat[0], requested_name)
        if actual is None or not all(
            _field_name(item, requested_name) is not None for item in grid.flat
        ):
            raise ValueError(
                f"Every factorial cell must contain measure '{requested_name}'"
            )
        if actual not in selected:
            selected.append(actual)
    return tuple(selected)


def _factorial_landscape_values(
    grid: NDArray[np.object_],
    measures: tuple[str, ...],
    supplied_frex: ArrayLike | None,
) -> tuple[dict[str, NDArray[np.object_]], dict[str, FloatArray], IntArray]:
    """Validate every selected measure and participant dimension."""
    values: dict[str, NDArray[np.object_]] = {}
    frequency_axes: dict[str, FloatArray] = {}
    participant_counts = np.zeros(grid.shape, dtype=np.int64)
    for measure_index, measure in enumerate(measures):
        current_grid = np.empty(grid.shape, dtype=object)
        reference_frex: FloatArray | None = None
        reference_features: int | None = None
        for index in np.ndindex(grid.shape):
            frequencies = _frequency_axis(grid[index], measure, supplied_frex)
            current = _canonical_values(
                grid[index], measure, frequencies, f"LANDSCAPE_GRID{index}"
            )
            if reference_frex is None:
                reference_frex = frequencies
                reference_features = current.shape[0]
            elif current.shape[0] != reference_features or not np.allclose(
                frequencies, reference_frex, atol=1e-10, rtol=0.0
            ):
                raise ValueError(
                    f"All factorial cells must share the feature dimension and "
                    f"frequency grid for measure '{measure}'"
                )
            if current.shape[2] < 2:
                raise ValueError(
                    "Every factorial cell must contain at least two participants"
                )
            if measure_index == 0:
                participant_counts[index] = current.shape[2]
            elif current.shape[2] != participant_counts[index]:
                raise ValueError(
                    "Selected measures within each factorial cell must have "
                    "matching participant counts"
                )
            current_grid[index] = current
        assert reference_frex is not None
        values[measure] = current_grid
        frequency_axes[measure] = reference_frex
    return values, frequency_axes, participant_counts


def _factorial_fdr(
    p_values: FloatArray,
    correction: str,
    alpha: float,
) -> tuple[FloatArray, FloatArray, BoolArray]:
    """Correct frequencies independently within every component."""
    q_values = np.full(p_values.shape, np.nan)
    thresholds = np.full(p_values.shape[0], np.nan)
    if correction == "fdr":
        for component in range(p_values.shape[0]):
            q_values[component], thresholds[component] = _bh_fdr(
                p_values[component], alpha
            )
    else:
        q_values[:] = p_values
        thresholds[:] = alpha
    return q_values, thresholds, np.isfinite(q_values) & (q_values <= alpha)


def _plot_factorial_measure(
    result: FREQNESSLandscapeFactorialMeasureResult,
    correction: str,
    alpha: float,
) -> list[Figure]:
    """Plot one omnibus statistic landscape per selected component."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "FREQNESS_LandscapeComparison plotting requires matplotlib. "
            "Install freqness[plot]."
        ) from error

    figures: list[Figure] = []
    for component_index, component in enumerate(result.comps):
        figure, axis = plt.subplots(figsize=(9, 5), constrained_layout=True)
        axis.plot(
            result.frex,
            result.statistic[component_index],
            color=plt.get_cmap("viridis")(0.45),
            linewidth=2,
            label=result.statistic_name,
        )
        significant = result.significant[component_index]
        if np.any(significant):
            axis.scatter(
                result.frex[significant],
                result.statistic[component_index, significant],
                s=48,
                color=plt.get_cmap("viridis")(0.95),
                label=(
                    f"FDR q ≤ {alpha:g}"
                    if correction == "fdr"
                    else f"p ≤ {alpha:g}"
                ),
                zorder=3,
            )
        axis.set_xlabel("Frequency (Hz)", fontweight="bold")
        axis.set_ylabel(result.statistic_name, fontweight="bold")
        axis.set_title(
            f"{result.measure} - {result.name} - component #{component}",
            fontweight="bold",
        )
        axis.grid(True, alpha=0.25)
        axis.spines[["top", "right"]].set_visible(False)
        axis.legend(frameon=False)
        figures.append(figure)
    return figures


def _factorial_landscape_comparison(
    LANDSCAPE_GRID: Any,
    *,
    measure: str | Sequence[str] | None,
    design: str | Sequence[str],
    test: str | None,
    which_comp: ArrayLike | None,
    frex: ArrayLike | None,
    frex2compare: ArrayLike | str | None,
    sample_labels: tuple[str, str] | None,
    factor_names: Sequence[str] | None,
    level_labels: Sequence[Sequence[str]] | None,
    subject_ids: Any,
    effects: str | Sequence[str],
    tail: str,
    alpha: float,
    correction: str,
    frex_tolerance: float | None,
    n_permutations: int | str,
    random_seed: int,
    plot_results: bool,
    show: bool,
) -> FREQNESSLandscapeFactorialComparisonResult:
    grid, factor_sizes, cell_levels = normalize_grid(LANDSCAPE_GRID)
    names = validate_factor_names(factor_names, len(factor_sizes))
    labels = validate_level_labels(level_labels, names, factor_sizes, sample_labels)
    factor_design = validate_design(design, len(factor_sizes))
    test_key, test_requested = validate_test(
        test, factor_design, factor_sizes, pairwise_default="rank"
    )
    definitions = effect_definitions(names, effects)
    tail_key = str(tail).lower()
    if tail_key not in {"both", "left", "right"}:
        raise ValueError("tail must be 'both', 'left', or 'right'")
    if tail_key != "both" and test_key not in {"ttest", "signrank", "ranksum"}:
        warnings.warn(
            "Factorial omnibus tests are non-directional. Using tail='both' instead.",
            UserWarning,
            stacklevel=3,
        )
        tail_key = "both"
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between 0 and 1")
    correction_key = str(correction).lower()
    if correction_key not in {"fdr", "none"}:
        raise ValueError("correction must be 'fdr' or 'none'")
    if not isinstance(random_seed, (int, np.integer)) or isinstance(
        random_seed, (bool, np.bool_)
    ):
        raise ValueError("random_seed must be an integer")
    if not isinstance(plot_results, (bool, np.bool_)) or not isinstance(
        show, (bool, np.bool_)
    ):
        raise TypeError("plot_results and show must be Boolean")

    measures = _factorial_measure_selection(grid, measure)
    value_grids, frequency_axes, participant_counts = _factorial_landscape_values(
        grid, measures, frex
    )
    orders, identifiers = match_subjects(
        participant_counts, cell_levels, factor_design, subject_ids
    )
    for measure_name in measures:
        for index in np.ndindex(grid.shape):
            value_grids[measure_name][index] = value_grids[measure_name][index][
                :, :, orders[index]
            ]

    print(
        "\nFREQNESS Landscape Comparison: testing a "
        f"{'/'.join(factor_design)} design with {len(factor_sizes)} factors, "
        f"{grid.size} factorial cells, and {len(measures)} measure(s) "
        f"using {test_key}."
    )
    effect_outputs = {
        definition.key: FREQNESSLandscapeFactorialEffect(
            name=definition.name,
            key=definition.key,
            factors=definition.factors,
            results={},
        )
        for definition in definitions
    }
    cell_means: dict[str, NDArray[np.object_]] = {}
    figures: list[Figure] = []
    permutations_used = 0

    for measure_name in measures:
        frequencies = frequency_axes[measure_name]
        components = _component_selection(
            which_comp, value_grids[measure_name].flat[0].shape[0]
        )
        requested, resolved, frequency_indices = _frequency_selection(
            frequencies, frex2compare, frex_tolerance
        )
        selected = np.empty(grid.shape, dtype=object)
        for index in np.ndindex(grid.shape):
            selected[index] = value_grids[measure_name][index][
                components - 1
            ][:, frequency_indices, :]
        Y, observation_levels, participant_index, means, _ = long_data(
            selected, cell_levels, factor_design
        )
        cell_means[measure_name] = means
        print(
            f"  {measure_name}: {resolved.size} frequencies, "
            f"component(s) {components.tolist()}."
        )

        pairwise_tests = {"ttest", "signrank", "ranksum"}
        if (
            len(factor_sizes) == 1
            and factor_sizes[0] == 2
            and test_key in pairwise_tests
        ):
            pairwise = _compare_measure(
                measure_name,
                value_grids[measure_name].flat[0],
                value_grids[measure_name].flat[1],
                requested,
                resolved,
                frequency_indices,
                components,
                factor_design[0],
                test_key,
                tail_key,
                float(alpha),
                correction_key,
            )
            definition = definitions[0]
            output = FREQNESSLandscapeFactorialMeasureResult(
                measure=measure_name,
                name=definition.name,
                key=definition.key,
                factors=definition.factors,
                requested_frex=requested,
                frex=resolved,
                idx_frex=frequency_indices,
                comps=components,
                values=selected,
                statistic=pairwise.statistic,
                statistic_name=pairwise.statistic_name,
                effect_size=pairwise.effect_size,
                effect_size_name=pairwise.effect_size_name,
                p_uncorrected=pairwise.p_uncorrected,
                q_fdr=pairwise.q_fdr,
                significant=pairwise.significant,
                fdr_threshold=pairwise.fdr_threshold,
                df1=np.full(pairwise.statistic.shape, np.nan),
                df2=np.full(pairwise.statistic.shape, np.nan),
            )
            if plot_results:
                pairwise.figures.extend(
                    _plot_measure(
                        pairwise,
                        (labels[0][0], labels[0][1]),
                        correction_key,
                        float(alpha),
                    )
                )
                output.figures.extend(pairwise.figures)
                figures.extend(pairwise.figures)
            effect_outputs[definition.key].results[measure_name] = output
            continue

        if test_key in {"friedman", "kruskalwallis"}:
            raw_results, models = rank_statistics(
                Y, observation_levels, participant_index, test_key, definitions
            )
        else:
            raw_results, models = anova_statistics(
                Y,
                observation_levels,
                participant_index,
                factor_sizes,
                factor_design,
                definitions,
            )
        if test_key == "permutation_anova":
            permutations_used = permutation_count(n_permutations)
            print(
                f"    Computing {permutations_used} design-aware "
                "Monte Carlo permutations."
            )
            rng = np.random.default_rng(random_seed)
            for raw_result, model in zip(raw_results, models, strict=True):
                raw_result["p_uncorrected"] = permutation_anova_p_values(
                    model,
                    raw_result["statistic"],
                    permutations_used,
                    rng,
                )

        for raw_result in raw_results:
            definition = raw_result["definition"]
            shape = (components.size, resolved.size)
            p_uncorrected = np.reshape(raw_result["p_uncorrected"], shape)
            q_values, thresholds, significant = _factorial_fdr(
                p_uncorrected, correction_key, float(alpha)
            )
            output = FREQNESSLandscapeFactorialMeasureResult(
                measure=measure_name,
                name=definition.name,
                key=definition.key,
                factors=definition.factors,
                requested_frex=requested,
                frex=resolved,
                idx_frex=frequency_indices,
                comps=components,
                values=selected,
                statistic=np.reshape(raw_result["statistic"], shape),
                statistic_name=raw_result["statistic_name"],
                effect_size=np.reshape(raw_result["effect_size"], shape),
                effect_size_name=raw_result["effect_size_name"],
                p_uncorrected=p_uncorrected,
                q_fdr=q_values,
                significant=significant,
                fdr_threshold=thresholds,
                df1=np.reshape(raw_result["df1"], shape),
                df2=np.reshape(raw_result["df2"], shape),
            )
            if plot_results:
                output.figures.extend(
                    _plot_factorial_measure(output, correction_key, float(alpha))
                )
                figures.extend(output.figures)
            effect_outputs[definition.key].results[measure_name] = output

    if figures and show:
        import matplotlib.pyplot as plt

        plt.show()
    return FREQNESSLandscapeFactorialComparisonResult(
        analysis=(
            "pairwise"
            if len(factor_sizes) == 1
            and factor_sizes[0] == 2
            and test_key in {"ttest", "signrank", "ranksum"}
            else "factorial"
        ),
        factor_names=names,
        level_labels=labels,
        factor_design=factor_design,
        design=factor_design,
        input_shape=factor_sizes,
        test=test_key,
        test_requested=test_requested,
        tail=tail_key,
        alpha=float(alpha),
        correction=correction_key,
        correction_scope=(
            "frequencies separately for each effect, measure, and component"
        ),
        measures=measures,
        effects=effect_outputs,
        effect_names=tuple(effect_outputs),
        cell_means=cell_means,
        cell_n_participants=participant_counts,
        subject_ids=identifiers,
        n_permutations=permutations_used,
        random_seed=int(random_seed),
        figures=figures,
    )


def FREQNESS_LandscapeComparison(
    LANDSCAPE_A: Any,
    LANDSCAPE_B: Any = None,
    *,
    measure: str | Sequence[str] | None = None,
    design: str | Sequence[str] = "within",
    test: str | None = None,
    which_comp: ArrayLike | None = None,
    frex: ArrayLike | None = None,
    frex2compare: ArrayLike | str | None = None,
    sample_labels: tuple[str, str] | None = None,
    factor_names: Sequence[str] | None = None,
    level_labels: Sequence[Sequence[str]] | None = None,
    subject_ids: Any = None,
    effects: str | Sequence[str] = "all",
    tail: str = "both",
    alpha: float = 0.05,
    correction: str = "fdr",
    frex_tolerance: float | None = None,
    n_permutations: int | str = 10_000,
    random_seed: int = 0,
    plot_results: bool = True,
    show: bool = True,
) -> FREQNESSLandscapeComparisonResult | FREQNESSLandscapeFactorialComparisonResult:
    """Compare frequency landscapes using pairwise or factorial designs.

    Existing calls pass two landscapes. Generalized designs pass one nested
    grid whose axes are experimental factors and whose cells contain the
    corresponding participant-level landscapes.
    """
    if is_factorial_grid(LANDSCAPE_A):
        if LANDSCAPE_B is not None:
            raise TypeError(
                "A factorial landscape grid is followed directly by keyword options"
            )
        return _factorial_landscape_comparison(
            LANDSCAPE_A,
            measure=measure,
            design=design,
            test=test,
            which_comp=which_comp,
            frex=frex,
            frex2compare=frex2compare,
            sample_labels=sample_labels,
            factor_names=factor_names,
            level_labels=level_labels,
            subject_ids=subject_ids,
            effects=effects,
            tail=tail,
            alpha=alpha,
            correction=correction,
            frex_tolerance=frex_tolerance,
            n_permutations=n_permutations,
            random_seed=random_seed,
            plot_results=plot_results,
            show=show,
        )
    if LANDSCAPE_B is None:
        raise TypeError("LANDSCAPE_A and LANDSCAPE_B are required")
    if not isinstance(design, str):
        raise ValueError("Pairwise design must be 'within' or 'between'")
    labels = ("A", "B") if sample_labels is None else sample_labels
    return _pairwise_landscape_comparison(
        LANDSCAPE_A,
        LANDSCAPE_B,
        measure=measure,
        design=design,
        test=test,
        which_comp=which_comp,
        frex=frex,
        frex2compare=frex2compare,
        sample_labels=labels,
        tail=tail,
        alpha=alpha,
        correction=correction,
        frex_tolerance=frex_tolerance,
        plot_results=plot_results,
        show=show,
    )


landscape_comparison = FREQNESS_LandscapeComparison


__all__ = [
    "FREQNESS_LandscapeComparison",
    "FREQNESSLandscapeComparisonResult",
    "FREQNESSLandscapeFactorialComparisonResult",
    "FREQNESSLandscapeFactorialEffect",
    "FREQNESSLandscapeFactorialMeasureResult",
    "FREQNESSLandscapeMeasureResult",
    "landscape_comparison",
]
