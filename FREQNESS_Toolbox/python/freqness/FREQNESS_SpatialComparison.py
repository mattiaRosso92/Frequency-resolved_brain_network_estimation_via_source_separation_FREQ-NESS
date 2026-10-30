"""Pairwise or factorial comparison of spatial activation patterns.

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
from scipy.spatial import cKDTree

from ._factorial_utils import (
    anova_statistics,
    effect_definitions,
    is_factorial_grid,
    long_data,
    match_subjects,
    normalize_grid,
    permutation_count,
    permute_model,
    rank_statistics,
    validate_design,
    validate_factor_names,
    validate_level_labels,
    validate_test,
)
from ._repeated_spatial import repeated_spatial_models, repeated_spatial_statistic

if TYPE_CHECKING:
    from matplotlib.figure import Figure


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
BoolArray = NDArray[np.bool_]


@dataclass(slots=True)
class FREQNESSSpatialCluster:
    """One observed supra-threshold cluster; ``peak_voxel`` is zero-based."""

    frex: float
    comp: int
    polarity: int
    size: int
    mass: float
    p_corrected: float
    peak_t: float
    peak_statistic: float
    peak_voxel: int
    peak_mni: FloatArray


@dataclass(slots=True)
class FREQNESSSpatialComparisonResult:
    """Outputs of :func:`FREQNESS_SpatialComparison`."""

    requested_frex: FloatArray
    frex: FloatArray
    idx_frex: IntArray
    comps: IntArray
    condition_labels: tuple[str, str]
    design: str
    test: str
    n_participants: int | tuple[int, int]
    n_A: int
    n_B: int
    mean_difference: FloatArray
    statistic: FloatArray
    statistic_name: str
    t_values: FloatArray
    effect_size: FloatArray
    effect_size_name: str
    effect_size_dz: FloatArray
    p_uncorrected: FloatArray
    p_cluster_corrected: FloatArray
    significant_mask: BoolArray
    cluster_labels: IntArray
    clusters: list[FREQNESSSpatialCluster]
    map_mode: str
    normalization: str
    cluster_forming_alpha: float
    alpha: float
    connectivity: int
    grid_spacing: FloatArray
    correction_scope: str
    n_permutations: int
    exact_permutations: bool
    random_seed: int
    sign_flips_a: FloatArray
    sign_flips_b: FloatArray


@dataclass(slots=True)
class FREQNESSSpatialFactorialCluster:
    """One cluster associated with a factorial main effect or interaction."""

    effect: str
    frex: float
    comp: int
    size: int
    mass: float
    p_corrected: float
    peak_f: float
    peak_statistic: float
    peak_voxel: int
    peak_mni: FloatArray


@dataclass(slots=True)
class FREQNESSSpatialFactorialEffect:
    """Spatial statistics for one factorial effect."""

    name: str
    key: str
    factors: tuple[int, ...]
    statistic: FloatArray
    statistic_name: str
    effect_size: FloatArray
    effect_size_name: str
    p_uncorrected: FloatArray
    p_cluster_corrected: FloatArray
    significant_mask: BoolArray
    cluster_labels: IntArray
    clusters: list[FREQNESSSpatialFactorialCluster]
    df1: FloatArray
    df2: FloatArray
    epsilon_gg: FloatArray | None = None
    df1_gg: FloatArray | None = None
    df2_gg: FloatArray | None = None
    contrast_matrix: FloatArray | None = None
    ss_effect: FloatArray | None = None
    ss_error: FloatArray | None = None
    p_uncorrected_method: str = "parametric reference"


@dataclass(slots=True)
class FREQNESSSpatialFactorialComparisonResult:
    """Outputs for a generalized factorial spatial comparison."""

    analysis: str
    requested_frex: FloatArray
    frex: FloatArray
    idx_frex: IntArray
    comps: IntArray
    factor_names: tuple[str, ...]
    level_labels: tuple[tuple[str, ...], ...]
    factor_design: tuple[str, ...]
    design: tuple[str, ...]
    input_shape: tuple[int, ...]
    test: str
    test_requested: str
    effects: dict[str, FREQNESSSpatialFactorialEffect]
    effect_names: tuple[str, ...]
    cell_means: NDArray[np.object_]
    cell_n_participants: IntArray
    subject_ids: NDArray[np.object_]
    map_mode: str
    normalization: str
    cluster_forming_alpha: float
    alpha: float
    connectivity: int
    grid_spacing: FloatArray
    correction_scope: str
    n_permutations: int
    exact_permutations: bool
    random_seed: int
    sign_flips: NDArray[np.object_]
    figures: list[Figure] = field(default_factory=list)
    inference_method: str | None = None
    permutation_unit: str | None = None
    permutation_assumption: str | None = None
    n_participants: int | None = None
    contrast_cell_levels: IntArray | None = None


def _field(container: Any, name: str) -> Any:
    if isinstance(container, Mapping) and name in container:
        return container[name]
    if hasattr(container, name):
        return getattr(container, name)
    raise ValueError(f"FREQ.{name} is required")


def _canonical_freq(freq: Any, label: str) -> tuple[FloatArray, FloatArray]:
    patterns = np.asarray(_field(freq, "pats"))
    if not np.issubdtype(patterns.dtype, np.number) or np.iscomplexobj(patterns):
        raise TypeError(f"{label}.pats must contain real numeric values")
    patterns = np.asarray(patterns, dtype=float)
    if patterns.ndim == 2:
        patterns = patterns[:, :, np.newaxis, np.newaxis]
    elif patterns.ndim == 3:
        patterns = patterns[..., np.newaxis]
    if patterns.ndim != 4 or min(patterns.shape) == 0:
        raise ValueError(
            f"{label}.pats must have shape "
            "(voxels, components[, frequencies, participants])"
        )
    if not np.all(np.isfinite(patterns)):
        raise ValueError(f"{label}.pats contains NaN or infinite values")

    frequency_axis = np.asarray(_field(freq, "frex"), dtype=float).reshape(-1)
    if frequency_axis.size != patterns.shape[2]:
        raise ValueError(f"{label}.frex and {label}.pats have incompatible dimensions")
    if not np.all(np.isfinite(frequency_axis)) or np.any(np.diff(frequency_axis) <= 0):
        raise ValueError(f"{label}.frex must be finite and strictly increasing")
    return patterns, frequency_axis


def _frequency_selection(
    frequency_axis: FloatArray,
    requested: ArrayLike | str,
    tolerance: float | None,
) -> tuple[FloatArray, FloatArray, IntArray]:
    if isinstance(requested, str):
        if requested.lower() != "all":
            raise ValueError("frex2compare must be numeric or the string 'all'")
        indices = np.arange(frequency_axis.size, dtype=np.int64)
        return frequency_axis.copy(), frequency_axis.copy(), indices

    requested_values = np.asarray(requested, dtype=float).reshape(-1)
    if requested_values.size == 0 or not np.all(np.isfinite(requested_values)):
        raise ValueError("frex2compare must contain at least one finite value")
    if np.any(requested_values < frequency_axis[0]) or np.any(
        requested_values > frequency_axis[-1]
    ):
        raise ValueError("Requested frex2compare values fall outside the available range")

    if tolerance is None:
        tolerance_value = (
            np.inf
            if frequency_axis.size == 1
            else float(np.median(np.diff(frequency_axis)) / 2.0)
        )
    else:
        tolerance_value = float(tolerance)
        if not np.isfinite(tolerance_value) or tolerance_value < 0:
            raise ValueError("frex_tolerance must be a non-negative finite scalar")

    indices = np.empty(requested_values.size, dtype=np.int64)
    for index, frequency in enumerate(requested_values):
        distances = np.abs(frequency_axis - frequency)
        nearest = np.flatnonzero(np.isclose(distances, np.min(distances)))
        if nearest.size != 1:
            raise ValueError(f"Frequency {frequency:g} Hz is equally close to two bins")
        indices[index] = nearest[0]
        if distances[nearest[0]] > tolerance_value + np.finfo(float).eps:
            raise ValueError(
                f"No analyzed frequency lies within {tolerance_value:g} Hz of "
                f"{frequency:g} Hz"
            )
    if np.unique(indices).size != indices.size:
        raise ValueError("Two frex2compare values resolve to the same analyzed bin")
    resolved = frequency_axis[indices]
    if not np.allclose(requested_values, resolved, atol=1e-12, rtol=0.0):
        warnings.warn(
            "Some requested frex2compare values are not present in FREQ.frex. "
            "Using closest matches instead.",
            UserWarning,
            stacklevel=3,
        )
    return requested_values, resolved, indices


def _component_selection(components: ArrayLike | None, ncomponents: int) -> IntArray:
    if components is None:
        print("comps2compare not specified. Defaulting to the 1st component.")
        return np.array([1], dtype=np.int64)
    values = np.asarray(components)
    if values.ndim > 1 or values.size == 0:
        raise ValueError("comps2compare must be a non-empty vector")
    numeric = values.astype(float).reshape(-1)
    if (
        not np.all(np.isfinite(numeric))
        or not np.all(numeric == np.round(numeric))
        or np.any(numeric < 1)
        or np.any(numeric > ncomponents)
    ):
        raise ValueError(
            f"comps2compare must contain integers between 1 and {ncomponents}"
        )
    selected = numeric.astype(np.int64)
    if np.unique(selected).size != selected.size:
        raise ValueError("comps2compare must not contain duplicates")
    return selected


def _analysis_options(design: str, test: str | None) -> tuple[str, str]:
    design_key = str(design).lower()
    if design_key in {"within", "paired"}:
        design_key = "within"
    elif design_key in {"between", "independent"}:
        design_key = "between"
    else:
        raise ValueError("design must be 'within' or 'between'")

    test_key = "ttest" if test is None else str(test).lower()
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
    return design_key, test_key


def _prepare_maps(
    maps_a: FloatArray,
    maps_b: FloatArray,
    map_mode: str,
    normalization: str,
    design: str,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    prepared_a = maps_a.copy()
    prepared_b = maps_b.copy()
    flips_a = np.ones(maps_a.shape[1:], dtype=float)
    flips_b = np.ones(maps_b.shape[1:], dtype=float)

    if map_mode == "magnitude":
        prepared_a = np.abs(prepared_a)
        prepared_b = np.abs(prepared_b)
    elif map_mode == "signed_aligned":
        for freqi in range(prepared_a.shape[1]):
            for compi in range(prepared_a.shape[2]):
                if design == "within":
                    for subi in range(prepared_a.shape[3]):
                        if np.dot(
                            prepared_a[:, freqi, compi, subi],
                            prepared_b[:, freqi, compi, subi],
                        ) < 0:
                            prepared_b[:, freqi, compi, subi] *= -1
                            flips_b[freqi, compi, subi] *= -1
                    reference_maps = 0.5 * (
                        prepared_a[:, freqi, compi, :]
                        + prepared_b[:, freqi, compi, :]
                    )
                else:
                    reference_maps = np.concatenate(
                        (
                            prepared_a[:, freqi, compi, :],
                            prepared_b[:, freqi, compi, :],
                        ),
                        axis=1,
                    )
                reference = np.linalg.svd(
                    reference_maps, full_matrices=False
                )[0][:, 0]
                peak = int(np.argmax(np.abs(reference)))
                if reference[peak] < 0:
                    reference *= -1
                for subi in range(prepared_a.shape[3]):
                    if np.dot(
                        prepared_a[:, freqi, compi, subi], reference
                    ) < 0:
                        prepared_a[:, freqi, compi, subi] *= -1
                        flips_a[freqi, compi, subi] *= -1
                for subi in range(prepared_b.shape[3]):
                    if np.dot(
                        prepared_b[:, freqi, compi, subi], reference
                    ) < 0:
                        prepared_b[:, freqi, compi, subi] *= -1
                        flips_b[freqi, compi, subi] *= -1
    else:
        raise ValueError("map_mode must be 'magnitude' or 'signed_aligned'")

    if normalization == "l2":
        for array, label in ((prepared_a, "FREQ_A"), (prepared_b, "FREQ_B")):
            norms = np.linalg.norm(array, axis=0, keepdims=True)
            if np.any(norms <= np.finfo(float).eps):
                raise ValueError(f"{label} contains a zero-norm selected pattern")
            array /= norms
    elif normalization != "none":
        raise ValueError("normalization must be 'l2' or 'none'")
    return prepared_a, prepared_b, flips_a, flips_b


def _paired_moments(
    differences: FloatArray,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    nsubjects = differences.shape[-1]
    mean_difference = np.mean(differences, axis=-1)
    standard_deviation = np.std(differences, axis=-1, ddof=1)
    denominator = standard_deviation / np.sqrt(nsubjects)
    t_values = np.divide(
        mean_difference,
        denominator,
        out=np.zeros_like(mean_difference),
        where=denominator > np.finfo(float).eps,
    )
    zero_variance_effect = (denominator <= np.finfo(float).eps) & (
        np.abs(mean_difference) > np.finfo(float).eps
    )
    t_values[zero_variance_effect] = np.sign(mean_difference[zero_variance_effect]) * np.inf
    return mean_difference, standard_deviation, t_values


def _paired_statistics(differences: FloatArray) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    nsubjects = differences.shape[-1]
    mean_difference, standard_deviation, t_values = _paired_moments(differences)
    zero_variance_effect = (standard_deviation <= np.finfo(float).eps) & (
        np.abs(mean_difference) > np.finfo(float).eps
    )
    effect_size = np.divide(
        mean_difference,
        standard_deviation,
        out=np.zeros_like(mean_difference),
        where=standard_deviation > np.finfo(float).eps,
    )
    effect_size[zero_variance_effect] = np.sign(mean_difference[zero_variance_effect]) * np.inf
    p_values = 2.0 * stats.t.sf(np.abs(t_values), nsubjects - 1)
    return mean_difference, t_values, effect_size, p_values


def _within_statistics(
    differences: FloatArray,
    test: str,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, str, str]:
    if test == "ttest":
        mean_difference, statistic, effect_size, p_values = _paired_statistics(
            differences
        )
        return (
            mean_difference,
            statistic,
            effect_size,
            p_values,
            "t",
            "Cohen's dz",
        )

    mean_difference = np.mean(differences, axis=-1)
    nonzero = differences != 0
    absolute = np.where(nonzero, np.abs(differences), np.nan)
    ranks = stats.rankdata(
        absolute,
        axis=-1,
        method="average",
        nan_policy="omit",
    )
    ranks = np.nan_to_num(ranks, nan=0.0)
    signed_rank = np.sum(np.sign(differences) * ranks, axis=-1)
    n_nonzero = np.sum(nonzero, axis=-1)
    variance = n_nonzero * (n_nonzero + 1) * (2 * n_nonzero + 1) / 6.0
    statistic = np.divide(
        signed_rank,
        np.sqrt(variance),
        out=np.zeros_like(signed_rank),
        where=variance > 0,
    )
    total_rank = n_nonzero * (n_nonzero + 1) / 2.0
    effect_size = np.divide(
        signed_rank,
        total_rank,
        out=np.zeros_like(signed_rank),
        where=total_rank > 0,
    )
    p_values = 2.0 * stats.norm.sf(np.abs(statistic))
    p_values[n_nonzero == 0] = 1.0
    return (
        mean_difference,
        statistic,
        effect_size,
        p_values,
        "z",
        "matched rank-biserial r",
    )


def _between_statistics(
    maps_a: FloatArray,
    maps_b: FloatArray,
    test: str,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, str, str]:
    n_a = maps_a.shape[-1]
    n_b = maps_b.shape[-1]
    mean_a = np.mean(maps_a, axis=-1)
    mean_b = np.mean(maps_b, axis=-1)
    mean_difference = mean_b - mean_a

    if test == "ttest":
        variance_a = np.var(maps_a, axis=-1, ddof=1)
        variance_b = np.var(maps_b, axis=-1, ddof=1)
        term_a = variance_a / n_a
        term_b = variance_b / n_b
        denominator = np.sqrt(term_a + term_b)
        statistic = np.divide(
            mean_difference,
            denominator,
            out=np.zeros_like(mean_difference),
            where=denominator > np.finfo(float).eps,
        )
        zero_variance_effect = (
            denominator <= np.finfo(float).eps
        ) & (np.abs(mean_difference) > np.finfo(float).eps)
        statistic[zero_variance_effect] = (
            np.sign(mean_difference[zero_variance_effect]) * np.inf
        )
        degrees_freedom = np.divide(
            (term_a + term_b) ** 2,
            term_a**2 / (n_a - 1) + term_b**2 / (n_b - 1),
            out=np.full_like(mean_difference, np.inf),
            where=(term_a + term_b) > np.finfo(float).eps,
        )
        p_values = 2.0 * stats.t.sf(np.abs(statistic), degrees_freedom)
        p_values[(term_a + term_b) == 0] = np.where(
            mean_difference[(term_a + term_b) == 0] == 0,
            1.0,
            0.0,
        )

        pooled = np.sqrt(
            ((n_a - 1) * variance_a + (n_b - 1) * variance_b)
            / (n_a + n_b - 2)
        )
        cohens_d = np.divide(
            mean_difference,
            pooled,
            out=np.zeros_like(mean_difference),
            where=pooled > np.finfo(float).eps,
        )
        cohens_d[zero_variance_effect] = (
            np.sign(mean_difference[zero_variance_effect]) * np.inf
        )
        correction = 1.0 - 3.0 / (4.0 * (n_a + n_b) - 9.0)
        return (
            mean_difference,
            statistic,
            correction * cohens_d,
            p_values,
            "t",
            "Hedges' g",
        )

    combined = np.concatenate((maps_b, maps_a), axis=-1)
    ranks = stats.rankdata(combined, axis=-1, method="average")
    u_b = np.sum(ranks[..., :n_b], axis=-1) - n_b * (n_b + 1) / 2.0
    centered_u = u_b - n_a * n_b / 2.0
    variance = n_a * n_b * (n_a + n_b + 1) / 12.0
    statistic = centered_u / np.sqrt(variance)
    effect_size = 2.0 * u_b / (n_a * n_b) - 1.0
    p_values = 2.0 * stats.norm.sf(np.abs(statistic))
    return (
        mean_difference,
        statistic,
        effect_size,
        p_values,
        "z",
        "rank-biserial r",
    )


def _spatial_adjacency(
    coordinates: FloatArray, connectivity: int
) -> tuple[list[IntArray], FloatArray]:
    if connectivity not in (6, 18, 26):
        raise ValueError("connectivity must be 6, 18, or 26")
    spacing = np.full(3, np.nan)
    scaled = np.zeros_like(coordinates)
    for dimension in range(3):
        values = np.unique(coordinates[:, dimension])
        gaps = np.diff(values)
        gaps = gaps[gaps > np.finfo(float).eps]
        if gaps.size:
            spacing[dimension] = float(np.min(gaps))
            scaled[:, dimension] = (
                coordinates[:, dimension] - values[0]
            ) / spacing[dimension]
    radius = {6: 1.01, 18: np.sqrt(2.0) + 0.01, 26: np.sqrt(3.0) + 0.01}[connectivity]
    maximum_axes = {6: 1, 18: 2, 26: 3}[connectivity]
    tree = cKDTree(scaled)
    adjacency: list[IntArray] = []
    for index, candidates in enumerate(tree.query_ball_point(scaled, radius)):
        candidates_array = np.asarray(candidates, dtype=np.int64)
        candidates_array = candidates_array[candidates_array != index]
        delta = np.abs(scaled[candidates_array] - scaled[index])
        on_grid = np.all(np.abs(delta - np.round(delta)) < 1e-6, axis=1)
        immediate = np.all(delta <= 1.0 + 1e-6, axis=1)
        changed_axes = np.sum(delta > 1e-6, axis=1)
        candidates_array = candidates_array[
            on_grid
            & immediate
            & (changed_axes >= 1)
            & (changed_axes <= maximum_axes)
        ]
        adjacency.append(candidates_array)
    return adjacency, spacing


def _clusters(
    mask: BoolArray,
    adjacency: list[IntArray],
) -> list[IntArray]:
    visited = np.zeros(mask.size, dtype=bool)
    components: list[IntArray] = []
    for seed in np.flatnonzero(mask):
        if visited[seed]:
            continue
        stack = [int(seed)]
        visited[seed] = True
        members: list[int] = []
        while stack:
            current = stack.pop()
            members.append(current)
            for neighbor in adjacency[current]:
                if mask[neighbor] and not visited[neighbor]:
                    visited[neighbor] = True
                    stack.append(int(neighbor))
        components.append(np.asarray(members, dtype=np.int64))
    return components


def _maximum_cluster_mass(
    statistic: FloatArray,
    p_values: FloatArray,
    cluster_forming_alpha: float,
    adjacency: list[IntArray],
) -> float:
    maximum = 0.0
    for freqi in range(statistic.shape[1]):
        for compi in range(statistic.shape[2]):
            spatial_map = statistic[:, freqi, compi]
            significant = p_values[:, freqi, compi] <= cluster_forming_alpha
            for mask in (
                significant & (spatial_map > 0),
                significant & (spatial_map < 0),
            ):
                for members in _clusters(mask, adjacency):
                    maximum = max(maximum, float(np.sum(np.abs(spatial_map[members]))))
    return maximum


def _permutation_plan(
    nsubjects: int,
    n_permutations: int | str,
    random_seed: int,
) -> tuple[FloatArray, bool]:
    possible = 2**nsubjects
    if isinstance(n_permutations, str):
        mode = n_permutations.lower()
        if mode not in ("auto", "exact"):
            raise ValueError("n_permutations must be a positive integer, 'auto', or 'exact'")
        if mode == "exact" and possible > 1_048_576:
            raise ValueError("Exact enumeration is too large; use Monte Carlo permutations")
        exact = mode == "exact" or possible <= 65_536
        count = possible if exact else 10_000
    else:
        count = int(n_permutations)
        if count <= 0 or count != n_permutations:
            raise ValueError("n_permutations must be a positive integer")
        exact = count >= possible and possible <= 65_536
        count = possible if exact else count
    if exact:
        integers = np.arange(count, dtype=np.uint64)[:, None]
        bits = (integers >> np.arange(nsubjects, dtype=np.uint64)) & 1
        signs = 2.0 * bits.astype(float) - 1.0
    else:
        rng = np.random.default_rng(random_seed)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(count, nsubjects))
    return signs, exact


def _between_permutation_plan(
    n_a: int,
    n_b: int,
    n_permutations: int | str,
    random_seed: int,
) -> tuple[BoolArray, bool]:
    from itertools import combinations
    from math import comb

    total = n_a + n_b
    possible = comb(total, n_a)
    if isinstance(n_permutations, str):
        mode = n_permutations.lower()
        if mode not in ("auto", "exact"):
            raise ValueError(
                "n_permutations must be a positive integer, 'auto', or 'exact'"
            )
        if mode == "exact" and possible > 1_048_576:
            raise ValueError(
                "Exact enumeration is too large; use Monte Carlo permutations"
            )
        exact = mode == "exact" or possible <= 65_536
        count = possible if exact else 10_000
    else:
        count = int(n_permutations)
        if count <= 0 or count != n_permutations:
            raise ValueError("n_permutations must be a positive integer")
        exact = count >= possible and possible <= 65_536
        count = possible if exact else count

    assignments = np.zeros((count, total), dtype=bool)
    if exact:
        for permutation, indices in enumerate(combinations(range(total), n_a)):
            assignments[permutation, np.asarray(indices, dtype=np.int64)] = True
    else:
        rng = np.random.default_rng(random_seed)
        for permutation in range(count):
            indices = rng.choice(total, size=n_a, replace=False)
            assignments[permutation, indices] = True
    return assignments, exact


def _publication_spatial_colormap():
    """Return the sequential palette shared by comparison figures."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "freqness_spatial",
        (
            (0.2670, 0.0049, 0.3294),
            (0.2297, 0.3224, 0.5457),
            (0.1276, 0.5669, 0.5506),
            (0.3692, 0.7889, 0.3829),
            (0.9932, 0.9062, 0.1439),
        ),
    )


def _format_publication_spatial_axis(axis, coordinates: FloatArray) -> None:
    """Apply a consistent orthographic MNI view without decorative axes."""
    extent = np.ptp(coordinates, axis=0)
    axis.set_box_aspect(np.maximum(extent, 1.0))
    axis.set_proj_type("ortho")
    axis.view_init(elev=18, azim=-62)
    axis.set_axis_off()


def _plot_results(
    coordinates: FloatArray,
    result: FREQNESSSpatialComparisonResult,
    show: bool,
) -> None:
    import matplotlib.pyplot as plt

    for freqi, frequency in enumerate(result.frex):
        for compi, component in enumerate(result.comps):
            difference = result.mean_difference[:, freqi, compi]
            significant = result.significant_mask[:, freqi, compi]
            limit = max(float(np.max(np.abs(difference))), np.finfo(float).eps)
            figure = plt.figure(figsize=(11.5, 5.2), dpi=140, facecolor="white")
            axes = []
            points = None
            for panel, (values, mask, title) in enumerate(
                (
                    (difference, np.ones(difference.size, dtype=bool), "Mean difference"),
                    (difference, significant, "Cluster-corrected difference"),
                ),
                start=1,
            ):
                axis = figure.add_subplot(1, 2, panel, projection="3d")
                axes.append(axis)
                axis.scatter(
                    *coordinates.T,
                    s=9,
                    c="#CBD2D9",
                    alpha=0.65,
                    linewidths=0,
                    depthshade=False,
                )
                if np.any(mask):
                    points = axis.scatter(
                        *coordinates[mask].T,
                        c=values[mask],
                        cmap="coolwarm",
                        vmin=-limit,
                        vmax=limit,
                        s=34 if panel == 1 else 50,
                        alpha=0.95,
                        edgecolors="#17212B" if panel == 2 else "none",
                        linewidths=0.7 if panel == 2 else 0,
                        depthshade=False,
                    )
                else:
                    axis.text2D(
                        0.5,
                        0.04,
                        "No corrected clusters",
                        transform=axis.transAxes,
                        ha="center",
                        color="#66717D",
                        fontsize=9,
                    )
                axis.set_title(title, fontsize=12, fontweight="bold", pad=8)
                _format_publication_spatial_axis(axis, coordinates)
            if points is not None:
                colorbar = figure.colorbar(
                    points, ax=axes, shrink=0.68, aspect=24, pad=0.015
                )
                colorbar.set_label("Mean difference", fontsize=10)
                colorbar.outline.set_linewidth(0.6)
            figure.suptitle(
                f"{result.condition_labels[1]} − {result.condition_labels[0]}  ·  "
                f"{frequency:g} Hz  ·  Component {component}",
                fontsize=15,
                fontweight="bold",
                y=0.96,
            )
            figure.text(
                0.5,
                0.025,
                f"Outlined markers: cluster-corrected p ≤ {result.alpha:g}",
                ha="center",
                fontsize=9,
                color="#52606D",
            )
            figure.subplots_adjust(left=0.01, right=0.91, top=0.86, bottom=0.08)
    if show:
        plt.show()


def _pairwise_spatial_comparison(
    FREQ_A: Any,
    FREQ_B: Any,
    MNI: ArrayLike,
    *,
    frex2compare: ArrayLike | str,
    comps2compare: ArrayLike | None = None,
    design: str = "within",
    test: str | None = None,
    condition_labels: tuple[str, str] = ("A", "B"),
    map_mode: str = "magnitude",
    normalization: str = "l2",
    cluster_forming_alpha: float = 0.001,
    alpha: float = 0.05,
    connectivity: int = 6,
    frex_tolerance: float | None = None,
    correction_scope: str = "all_maps",
    n_permutations: int | str = "auto",
    random_seed: int = 0,
    plot_results: bool = True,
    show: bool = True,
) -> FREQNESSSpatialComparisonResult:
    """Compare FREQ-NESS spatial patterns with cluster correction.

    For a within-subject design, participant order must match across
    ``FREQ_A`` and ``FREQ_B``. Between-group designs may have unequal sample
    sizes. Incompatible non-parametric test names are corrected with a warning.
    Component numbers retain the MATLAB convention and are one-based, while
    ``idx_frex`` is zero-based for direct NumPy indexing.
    """
    design_key, test_key = _analysis_options(design, test)
    patterns_a, frequency_axis_a = _canonical_freq(FREQ_A, "FREQ_A")
    patterns_b, frequency_axis_b = _canonical_freq(FREQ_B, "FREQ_B")
    if patterns_a.shape[:3] != patterns_b.shape[:3]:
        raise ValueError(
            "FREQ_A.pats and FREQ_B.pats must have matching voxel, component, "
            "and frequency dimensions"
        )
    if not np.allclose(frequency_axis_a, frequency_axis_b, atol=1e-10, rtol=0.0):
        raise ValueError("FREQ_A.frex and FREQ_B.frex must contain the same grid")
    n_a = patterns_a.shape[3]
    n_b = patterns_b.shape[3]
    if design_key == "within" and n_a != n_b:
        raise ValueError(
            "Within-subject comparisons require matching participant counts"
        )
    if n_a < 2 or n_b < 2:
        raise ValueError("Each sample must contain at least two participants")

    coordinates = np.asarray(MNI, dtype=float)
    if coordinates.shape != (patterns_a.shape[0], 3) or not np.all(np.isfinite(coordinates)):
        raise ValueError("MNI must have shape (n_voxels, 3) and contain finite values")
    if len(condition_labels) != 2 or not all(isinstance(label, str) for label in condition_labels):
        raise ValueError("condition_labels must contain two strings")
    if correction_scope != "all_maps":
        raise ValueError("correction_scope currently supports only 'all_maps'")
    for value, name in ((cluster_forming_alpha, "cluster_forming_alpha"), (alpha, "alpha")):
        if not np.isfinite(value) or not 0 < value < 1:
            raise ValueError(f"{name} must lie strictly between 0 and 1")
    if not isinstance(random_seed, (int, np.integer)) or isinstance(random_seed, (bool, np.bool_)):
        raise ValueError("random_seed must be an integer")
    if not isinstance(plot_results, (bool, np.bool_)) or not isinstance(show, (bool, np.bool_)):
        raise TypeError("plot_results and show must be Boolean")

    requested, resolved, frequency_indices = _frequency_selection(
        frequency_axis_a, frex2compare, frex_tolerance
    )
    selected_components = _component_selection(comps2compare, patterns_a.shape[1])
    component_indices = selected_components - 1
    selected_a = patterns_a[:, component_indices, :, :][:, :, frequency_indices, :]
    selected_b = patterns_b[:, component_indices, :, :][:, :, frequency_indices, :]
    selected_a = np.transpose(selected_a, (0, 2, 1, 3))
    selected_b = np.transpose(selected_b, (0, 2, 1, 3))
    selected_a, selected_b, flips_a, flips_b = _prepare_maps(
        selected_a, selected_b, map_mode, normalization, design_key
    )
    adjacency, grid_spacing = _spatial_adjacency(coordinates, connectivity)

    if design_key == "within":
        differences = selected_b - selected_a
        (
            mean_difference,
            statistic,
            effect_size,
            p_uncorrected,
            statistic_name,
            effect_size_name,
        ) = _within_statistics(differences, test_key)
        permutation_plan, exact = _permutation_plan(
            n_a, n_permutations, int(random_seed)
        )
        permutation_description = "paired sign-flip"
    else:
        (
            mean_difference,
            statistic,
            effect_size,
            p_uncorrected,
            statistic_name,
            effect_size_name,
        ) = _between_statistics(selected_a, selected_b, test_key)
        permutation_plan, exact = _between_permutation_plan(
            n_a, n_b, n_permutations, int(random_seed)
        )
        permutation_description = "independent label"

    print(
        "\nFREQNESS Spatial Comparison: comparing "
        f"{condition_labels[0]} and {condition_labels[1]} at "
        f"{resolved.size} frequencies and {selected_components.size} components "
        f"using a {design_key}-subject {test_key} analysis "
        f"(n={n_a}, n={n_b})."
    )
    print(
        f"Computing {permutation_plan.shape[0]} "
        f"{'exact' if exact else 'Monte Carlo'} "
        f"{permutation_description} permutations..."
    )
    null_maxima = np.zeros(permutation_plan.shape[0])
    progress_step = max(1, permutation_plan.shape[0] // 10)
    if design_key == "between":
        pooled_maps = np.concatenate((selected_a, selected_b), axis=-1)
    for permutation, assignment in enumerate(permutation_plan):
        if design_key == "within":
            _, permuted_statistic, _, permuted_p, _, _ = _within_statistics(
                differences * assignment.reshape(1, 1, 1, -1),
                test_key,
            )
        else:
            _, permuted_statistic, _, permuted_p, _, _ = _between_statistics(
                pooled_maps[..., assignment],
                pooled_maps[..., ~assignment],
                test_key,
            )
        null_maxima[permutation] = _maximum_cluster_mass(
            permuted_statistic,
            permuted_p,
            cluster_forming_alpha,
            adjacency,
        )
        if (
            (permutation + 1) % progress_step == 0
            or permutation + 1 == permutation_plan.shape[0]
        ):
            print(f"  Permutation {permutation + 1}/{permutation_plan.shape[0]}")

    corrected = np.ones_like(statistic)
    significant = np.zeros_like(statistic, dtype=bool)
    labels = np.zeros_like(statistic, dtype=np.int64)
    cluster_records: list[FREQNESSSpatialCluster] = []
    for freqi, frequency in enumerate(resolved):
        for compi, component in enumerate(selected_components):
            spatial_map = statistic[:, freqi, compi]
            spatial_p = p_uncorrected[:, freqi, compi]
            supra_threshold = spatial_p <= cluster_forming_alpha
            positive_label = 0
            negative_label = 0
            for polarity, mask in (
                (1, supra_threshold & (spatial_map > 0)),
                (-1, supra_threshold & (spatial_map < 0)),
            ):
                for members in _clusters(mask, adjacency):
                    mass = float(np.sum(np.abs(spatial_map[members])))
                    if exact:
                        p_corrected = float(np.mean(null_maxima >= mass))
                    else:
                        p_corrected = float(
                            (1 + np.count_nonzero(null_maxima >= mass))
                            / (permutation_plan.shape[0] + 1)
                        )
                    corrected[members, freqi, compi] = p_corrected
                    significant[members, freqi, compi] = p_corrected <= alpha
                    if polarity > 0:
                        positive_label += 1
                        cluster_label = positive_label
                    else:
                        negative_label += 1
                        cluster_label = -negative_label
                    labels[members, freqi, compi] = cluster_label
                    peak_local = int(np.argmax(np.abs(spatial_map[members])))
                    peak_voxel = int(members[peak_local])
                    cluster_records.append(
                        FREQNESSSpatialCluster(
                            frex=float(frequency),
                            comp=int(component),
                            polarity=polarity,
                            size=int(members.size),
                            mass=mass,
                            p_corrected=p_corrected,
                            peak_t=float(spatial_map[peak_voxel]),
                            peak_statistic=float(spatial_map[peak_voxel]),
                            peak_voxel=peak_voxel,
                            peak_mni=coordinates[peak_voxel].copy(),
                        )
                    )

    result = FREQNESSSpatialComparisonResult(
        requested_frex=requested,
        frex=resolved,
        idx_frex=frequency_indices,
        comps=selected_components,
        condition_labels=(condition_labels[0], condition_labels[1]),
        design=design_key,
        test=test_key,
        n_participants=n_a if design_key == "within" else (n_a, n_b),
        n_A=n_a,
        n_B=n_b,
        mean_difference=mean_difference,
        statistic=statistic,
        statistic_name=statistic_name,
        t_values=(
            statistic.copy()
            if test_key == "ttest"
            else np.full_like(statistic, np.nan)
        ),
        effect_size=effect_size,
        effect_size_name=effect_size_name,
        effect_size_dz=(
            effect_size.copy()
            if design_key == "within" and test_key == "ttest"
            else np.full_like(effect_size, np.nan)
        ),
        p_uncorrected=p_uncorrected,
        p_cluster_corrected=corrected,
        significant_mask=significant,
        cluster_labels=labels,
        clusters=cluster_records,
        map_mode=map_mode,
        normalization=normalization,
        cluster_forming_alpha=float(cluster_forming_alpha),
        alpha=float(alpha),
        connectivity=connectivity,
        grid_spacing=grid_spacing,
        correction_scope=correction_scope,
        n_permutations=permutation_plan.shape[0],
        exact_permutations=exact,
        random_seed=int(random_seed),
        sign_flips_a=flips_a,
        sign_flips_b=flips_b,
    )
    if plot_results:
        _plot_results(coordinates, result, bool(show))
    return result


def _factorial_maps(
    patterns: NDArray[np.object_],
    components: IntArray,
    frequency_indices: IntArray,
    map_mode: str,
    normalization: str,
) -> tuple[NDArray[np.object_], NDArray[np.object_]]:
    """Select, align, and normalize maps from every factorial cell."""
    maps = np.empty(patterns.shape, dtype=object)
    flips = np.empty(patterns.shape, dtype=object)
    for index in np.ndindex(patterns.shape):
        current = patterns[index][:, components - 1][:, :, frequency_indices, :]
        current = np.transpose(current, (0, 2, 1, 3)).copy()
        maps[index] = current
        flips[index] = np.ones(current.shape[1:], dtype=float)

    if map_mode == "magnitude":
        for index in np.ndindex(maps.shape):
            maps[index] = np.abs(maps[index])
    elif map_mode == "signed_aligned":
        reference = maps.flat[0][..., 0].copy()
        for freqi in range(reference.shape[1]):
            for compi in range(reference.shape[2]):
                if np.linalg.norm(reference[:, freqi, compi]) <= np.finfo(float).eps:
                    reference[:, freqi, compi] = 1.0
                for index in np.ndindex(maps.shape):
                    for subject in range(maps[index].shape[3]):
                        if np.dot(
                            maps[index][:, freqi, compi, subject],
                            reference[:, freqi, compi],
                        ) < 0:
                            maps[index][:, freqi, compi, subject] *= -1.0
                            flips[index][freqi, compi, subject] = -1.0
    else:
        raise ValueError("map_mode must be 'magnitude' or 'signed_aligned'")

    if normalization == "l2":
        for index in np.ndindex(maps.shape):
            norms = np.linalg.norm(maps[index], axis=0, keepdims=True)
            if np.any(norms <= np.finfo(float).eps):
                raise ValueError("A selected factorial pattern has zero L2 norm")
            maps[index] = maps[index] / norms
    elif normalization != "none":
        raise ValueError("normalization must be 'l2' or 'none'")
    return maps, flips


def _correct_factorial_clusters(
    result: dict[str, Any],
    coordinates: FloatArray,
    adjacency: list[IntArray],
    null_maxima: FloatArray,
    cluster_forming_alpha: float,
    alpha: float,
    n_permutations: int,
    exact: bool,
    frequencies: FloatArray,
    components: IntArray,
) -> FREQNESSSpatialFactorialEffect:
    """Apply maximum-cluster correction to one factorial effect."""
    statistic = result["statistic"]
    p_uncorrected = result["p_uncorrected"]
    corrected = np.ones(statistic.shape)
    significant = np.zeros(statistic.shape, dtype=bool)
    labels = np.zeros(statistic.shape, dtype=np.int64)
    records: list[FREQNESSSpatialFactorialCluster] = []
    definition = result["definition"]
    for freqi, frequency in enumerate(frequencies):
        for compi, component in enumerate(components):
            spatial_map = statistic[:, freqi, compi]
            mask = (
                (p_uncorrected[:, freqi, compi] <= cluster_forming_alpha)
                & (spatial_map > 0)
            )
            for cluster_index, members in enumerate(_clusters(mask, adjacency), start=1):
                mass = float(np.sum(spatial_map[members]))
                p_corrected = (
                    float(np.mean(null_maxima >= mass))
                    if exact
                    else float(
                        (1 + np.count_nonzero(null_maxima >= mass))
                        / (n_permutations + 1)
                    )
                )
                corrected[members, freqi, compi] = p_corrected
                significant[members, freqi, compi] = p_corrected <= alpha
                labels[members, freqi, compi] = cluster_index
                peak = int(members[np.argmax(spatial_map[members])])
                records.append(
                    FREQNESSSpatialFactorialCluster(
                        effect=definition.key,
                        frex=float(frequency),
                        comp=int(component),
                        size=int(members.size),
                        mass=mass,
                        p_corrected=p_corrected,
                        peak_f=float(spatial_map[peak]),
                        peak_statistic=float(spatial_map[peak]),
                        peak_voxel=peak,
                        peak_mni=coordinates[peak].copy(),
                    )
                )
    return FREQNESSSpatialFactorialEffect(
        name=definition.name,
        key=definition.key,
        factors=definition.factors,
        statistic=statistic,
        statistic_name=result["statistic_name"],
        effect_size=result["effect_size"],
        effect_size_name=result["effect_size_name"],
        p_uncorrected=p_uncorrected,
        p_cluster_corrected=corrected,
        significant_mask=significant,
        cluster_labels=labels,
        clusters=records,
        df1=result["df1"],
        df2=result["df2"],
        **{key: result[key] for key in (
            "epsilon_gg", "df1_gg", "df2_gg", "contrast_matrix",
            "ss_effect", "ss_error", "p_uncorrected_method"
        ) if key in result},
    )


def _plot_factorial_results(
    coordinates: FloatArray,
    result: FREQNESSSpatialFactorialComparisonResult,
    show: bool,
) -> list[Figure]:
    """Plot omnibus spatial statistics and corrected clusters."""
    import matplotlib.pyplot as plt

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    figures: list[Figure] = []
    cmap = _publication_spatial_colormap()
    for effect_name in result.effect_names:
        effect = result.effects[effect_name]
        for compi, component in enumerate(result.comps):
            figure = plt.figure(
                figsize=(4.5 * result.frex.size + 0.8, 5.0),
                dpi=140,
                facecolor="white",
            )
            limit = max(float(np.nanmax(effect.statistic[:, :, compi])), 1.0)
            axes = []
            for freqi, frequency in enumerate(result.frex):
                axis = figure.add_subplot(
                    1, result.frex.size, freqi + 1, projection="3d"
                )
                axes.append(axis)
                values = effect.statistic[:, freqi, compi]
                axis.scatter(
                    *coordinates.T,
                    s=10,
                    c="#CBD2D9",
                    alpha=0.68,
                    linewidths=0,
                    depthshade=False,
                )
                significant = effect.significant_mask[:, freqi, compi]
                if np.any(significant):
                    axis.scatter(
                        *coordinates[significant].T,
                        s=54,
                        c=values[significant],
                        cmap=cmap,
                        vmin=0,
                        vmax=limit,
                        edgecolors="#17212B",
                        linewidths=0.75,
                        depthshade=False,
                    )
                else:
                    axis.text2D(
                        0.5,
                        0.04,
                        "No corrected clusters",
                        transform=axis.transAxes,
                        ha="center",
                        color="#66717D",
                        fontsize=9,
                    )
                axis.set_title(
                    f"{frequency:g} Hz", fontsize=12, fontweight="bold", pad=6
                )
                _format_publication_spatial_axis(axis, coordinates)
            mapper = ScalarMappable(norm=Normalize(vmin=0, vmax=limit), cmap=cmap)
            colorbar = figure.colorbar(
                mapper, ax=axes, shrink=0.68, aspect=24, pad=0.018
            )
            colorbar.set_label(f"{effect.statistic_name} statistic", fontsize=10)
            colorbar.outline.set_linewidth(0.6)
            figure.suptitle(
                f"Spatial comparison — {effect.name}",
                fontsize=15,
                fontweight="bold",
                y=0.97,
            )
            figure.text(
                0.5,
                0.91,
                f"Component {component}  ·  Outlined markers: "
                f"cluster-corrected p ≤ {result.alpha:g}",
                ha="center",
                fontsize=9.5,
                color="#52606D",
            )
            figure.subplots_adjust(left=0.01, right=0.90, top=0.84, bottom=0.04)
            figures.append(figure)
    if figures and show:
        plt.show()
    return figures


def _factorial_spatial_comparison(
    FREQ_GRID: Any,
    MNI: ArrayLike,
    *,
    frex2compare: ArrayLike | str,
    comps2compare: ArrayLike | None,
    design: str | Sequence[str],
    test: str | None,
    condition_labels: tuple[str, str] | None,
    factor_names: Sequence[str] | None,
    level_labels: Sequence[Sequence[str]] | None,
    subject_ids: Any,
    effects: str | Sequence[str],
    map_mode: str,
    normalization: str,
    cluster_forming_alpha: float,
    alpha: float,
    connectivity: int,
    frex_tolerance: float | None,
    correction_scope: str,
    n_permutations: int | str,
    random_seed: int,
    plot_results: bool,
    show: bool,
) -> FREQNESSSpatialFactorialComparisonResult:
    grid, factor_sizes, cell_levels = normalize_grid(FREQ_GRID)
    names = validate_factor_names(factor_names, len(factor_sizes))
    labels = validate_level_labels(
        level_labels, names, factor_sizes, condition_labels
    )
    factor_design = validate_design(design, len(factor_sizes))
    test_key, test_requested = validate_test(
        test, factor_design, factor_sizes, pairwise_default="ttest"
    )
    definitions = effect_definitions(names, effects)
    if correction_scope not in {"all_maps", "all_effects_and_maps"}:
        raise ValueError(
            "correction_scope must be 'all_maps' or 'all_effects_and_maps'"
        )
    for value, name in (
        (cluster_forming_alpha, "cluster_forming_alpha"),
        (alpha, "alpha"),
    ):
        if not np.isfinite(value) or not 0 < value < 1:
            raise ValueError(f"{name} must lie strictly between 0 and 1")
    if not isinstance(random_seed, (int, np.integer)) or isinstance(
        random_seed, (bool, np.bool_)
    ):
        raise ValueError("random_seed must be an integer")

    patterns = np.empty(grid.shape, dtype=object)
    frequency_axis: FloatArray | None = None
    reference_shape: tuple[int, int, int] | None = None
    participant_counts = np.zeros(grid.shape, dtype=np.int64)
    for index in np.ndindex(grid.shape):
        current, current_frex = _canonical_freq(grid[index], f"FREQ_GRID{index}")
        if frequency_axis is None:
            frequency_axis = current_frex
            reference_shape = current.shape[:3]
        elif current.shape[:3] != reference_shape or not np.allclose(
            current_frex, frequency_axis, atol=1e-10, rtol=0.0
        ):
            raise ValueError(
                "All factorial FREQ inputs must share voxel, component, and "
                "frequency dimensions and the same FREQ.frex grid"
            )
        if current.shape[3] < 2:
            raise ValueError("Every factorial cell must contain at least two participants")
        patterns[index] = current
        participant_counts[index] = current.shape[3]
    assert frequency_axis is not None and reference_shape is not None
    coordinates = np.asarray(MNI, dtype=float)
    if coordinates.shape != (reference_shape[0], 3) or not np.all(
        np.isfinite(coordinates)
    ):
        raise ValueError("MNI must have shape (n_voxels, 3) and contain finite values")

    orders, identifiers = match_subjects(
        participant_counts, cell_levels, factor_design, subject_ids
    )
    for index in np.ndindex(grid.shape):
        patterns[index] = patterns[index][..., orders[index]]
    requested, resolved, frequency_indices = _frequency_selection(
        frequency_axis, frex2compare, frex_tolerance
    )
    components = _component_selection(comps2compare, reference_shape[1])

    pairwise_tests = {"ttest", "signrank", "ranksum"}
    if len(factor_sizes) == 1 and factor_sizes[0] == 2 and test_key in pairwise_tests:
        reordered = [
            {"pats": patterns.flat[index], "frex": frequency_axis}
            for index in range(2)
        ]
        pairwise = _pairwise_spatial_comparison(
            reordered[0],
            reordered[1],
            coordinates,
            frex2compare=frex2compare,
            comps2compare=comps2compare,
            design=factor_design[0],
            test=test_key,
            condition_labels=(labels[0][0], labels[0][1]),
            map_mode=map_mode,
            normalization=normalization,
            cluster_forming_alpha=cluster_forming_alpha,
            alpha=alpha,
            connectivity=connectivity,
            frex_tolerance=frex_tolerance,
            correction_scope="all_maps",
            n_permutations=n_permutations,
            random_seed=random_seed,
            plot_results=plot_results,
            show=show,
        )
        definition = definitions[0]
        effect = FREQNESSSpatialFactorialEffect(
            name=definition.name,
            key=definition.key,
            factors=definition.factors,
            statistic=pairwise.statistic,
            statistic_name=pairwise.statistic_name,
            effect_size=pairwise.effect_size,
            effect_size_name=pairwise.effect_size_name,
            p_uncorrected=pairwise.p_uncorrected,
            p_cluster_corrected=pairwise.p_cluster_corrected,
            significant_mask=pairwise.significant_mask,
            cluster_labels=pairwise.cluster_labels,
            clusters=[
                FREQNESSSpatialFactorialCluster(
                    effect=definition.key,
                    frex=cluster.frex,
                    comp=cluster.comp,
                    size=cluster.size,
                    mass=cluster.mass,
                    p_corrected=cluster.p_corrected,
                    peak_f=cluster.peak_statistic,
                    peak_statistic=cluster.peak_statistic,
                    peak_voxel=cluster.peak_voxel,
                    peak_mni=cluster.peak_mni,
                )
                for cluster in pairwise.clusters
            ],
            df1=np.full(pairwise.statistic.shape, np.nan),
            df2=np.full(pairwise.statistic.shape, np.nan),
        )
        prepared, flips = _factorial_maps(
            patterns, components, frequency_indices, map_mode, normalization
        )
        _, _, _, cell_means, cell_counts = long_data(
            prepared, cell_levels, factor_design
        )
        return FREQNESSSpatialFactorialComparisonResult(
            analysis="pairwise",
            requested_frex=pairwise.requested_frex,
            frex=pairwise.frex,
            idx_frex=pairwise.idx_frex,
            comps=pairwise.comps,
            factor_names=names,
            level_labels=labels,
            factor_design=factor_design,
            design=factor_design,
            input_shape=factor_sizes,
            test=test_key,
            test_requested=test_requested,
            effects={definition.key: effect},
            effect_names=(definition.key,),
            cell_means=cell_means,
            cell_n_participants=cell_counts,
            subject_ids=identifiers,
            map_mode=map_mode,
            normalization=normalization,
            cluster_forming_alpha=float(cluster_forming_alpha),
            alpha=float(alpha),
            connectivity=connectivity,
            grid_spacing=pairwise.grid_spacing,
            correction_scope="all_maps",
            n_permutations=pairwise.n_permutations,
            exact_permutations=pairwise.exact_permutations,
            random_seed=int(random_seed),
            sign_flips=flips,
        )

    maps, flips = _factorial_maps(
        patterns, components, frequency_indices, map_mode, normalization
    )
    Y, observation_levels, subject_index, cell_means, cell_counts = long_data(
        maps, cell_levels, factor_design
    )
    repeated_anova = all(value == "within" for value in factor_design) and test_key in {
        "anova", "permutation_anova"
    }
    if repeated_anova:
        models = repeated_spatial_models(maps, cell_levels, factor_sizes, definitions)
        raw_results = []
        for definition, model in zip(definitions, models):
            result = repeated_spatial_statistic(model)
            result.update(
                definition=definition, statistic_name="F",
                effect_size_name="partial eta squared",
                contrast_matrix=model.contrast_matrix,
                p_uncorrected_method="Greenhouse-Geisser F approximation",
            )
            raw_results.append(result)
    elif test_key in {"friedman", "kruskalwallis"}:
        raw_results, models = rank_statistics(
            Y, observation_levels, subject_index, test_key, definitions
        )
    else:
        if "within" in factor_design:
            warnings.warn(
                "Mixed-design spatial ANOVA retains the legacy residual-permutation "
                "method and requires within-block exchangeability; the corrected "
                "participant-contrast method applies only to all-within designs.",
                RuntimeWarning, stacklevel=2,
            )
        raw_results, models = anova_statistics(
            Y,
            observation_levels,
            subject_index,
            factor_sizes,
            factor_design,
            definitions,
        )
    nvoxels = reference_shape[0]
    nfrequencies = resolved.size
    ncomponents = components.size
    for result in raw_results:
        for key in ("statistic", "effect_size", "p_uncorrected", "df1", "df2",
                    "epsilon_gg", "df1_gg", "df2_gg", "ss_effect", "ss_error"):
            if key not in result:
                continue
            result[key] = np.reshape(
                result[key], (nvoxels, nfrequencies, ncomponents)
            )

    adjacency, spacing = _spatial_adjacency(coordinates, connectivity)
    exact = False
    if repeated_anova:
        sign_plan, exact = _permutation_plan(
            int(cell_counts.flat[0]), n_permutations, random_seed
        )
        count = sign_plan.shape[0]
    else:
        count = permutation_count(n_permutations)
    print(
        "\nFREQNESS Spatial Comparison: testing a "
        f"{'/'.join(factor_design)} design with {len(factor_sizes)} factors, "
        f"{grid.size} factorial cells, {resolved.size} frequencies, and "
        f"{components.size} components using {test_key}."
    )
    permutation_description = "exact participant sign-flip" if exact else "Monte Carlo"
    print(f"Computing {count} {permutation_description} cluster permutations...")
    rng = np.random.default_rng(random_seed)
    null_maxima = np.zeros(count)
    progress_step = max(1, count // 10)
    for permutation in range(count):
        maximum = 0.0
        for model in models:
            if repeated_anova:
                permuted = repeated_spatial_statistic(model, sign_plan[permutation])
                statistic, p_values = permuted["statistic"], permuted["p_uncorrected"]
            else:
                statistic, p_values = permute_model(model, rng)
            statistic = statistic.reshape(nvoxels, nfrequencies, ncomponents)
            p_values = p_values.reshape(nvoxels, nfrequencies, ncomponents)
            maximum = max(
                maximum,
                _maximum_cluster_mass(
                    statistic, p_values, cluster_forming_alpha, adjacency
                ),
            )
        null_maxima[permutation] = maximum
        if (permutation + 1) % progress_step == 0 or permutation + 1 == count:
            print(f"  Permutation {permutation + 1}/{count}")

    effect_results: dict[str, FREQNESSSpatialFactorialEffect] = {}
    for raw_result in raw_results:
        effect = _correct_factorial_clusters(
            raw_result,
            coordinates,
            adjacency,
            null_maxima,
            cluster_forming_alpha,
            alpha,
            count,
            exact,
            resolved,
            components,
        )
        effect_results[effect.key] = effect
    output = FREQNESSSpatialFactorialComparisonResult(
        analysis="factorial",
        requested_frex=requested,
        frex=resolved,
        idx_frex=frequency_indices,
        comps=components,
        factor_names=names,
        level_labels=labels,
        factor_design=factor_design,
        design=factor_design,
        input_shape=factor_sizes,
        test=test_key,
        test_requested=test_requested,
        effects=effect_results,
        effect_names=tuple(effect_results),
        cell_means=cell_means,
        cell_n_participants=cell_counts,
        subject_ids=identifiers,
        map_mode=map_mode,
        normalization=normalization,
        cluster_forming_alpha=float(cluster_forming_alpha),
        alpha=float(alpha),
        connectivity=connectivity,
        grid_spacing=spacing,
        correction_scope="all_effects_and_maps",
        n_permutations=count,
        exact_permutations=exact,
        random_seed=int(random_seed),
        sign_flips=flips,
        inference_method=("within_subject_contrast_signflip" if repeated_anova
                          else "legacy_factorial"),
        permutation_unit="participant" if repeated_anova else "design-dependent",
        permutation_assumption=(
            "independent participants; joint central symmetry of null participant contrast fields"
            if repeated_anova else "legacy design-specific exchangeability"
        ),
        n_participants=int(cell_counts.flat[0]) if repeated_anova else None,
        contrast_cell_levels=cell_levels.copy() if repeated_anova else None,
    )
    if plot_results:
        output.figures.extend(_plot_factorial_results(coordinates, output, bool(show)))
    return output


def FREQNESS_SpatialComparison(
    FREQ_A: Any,
    FREQ_B: Any = None,
    MNI: ArrayLike | None = None,
    *,
    frex2compare: ArrayLike | str,
    comps2compare: ArrayLike | None = None,
    design: str | Sequence[str] = "within",
    test: str | None = None,
    condition_labels: tuple[str, str] | None = None,
    factor_names: Sequence[str] | None = None,
    level_labels: Sequence[Sequence[str]] | None = None,
    subject_ids: Any = None,
    effects: str | Sequence[str] = "all",
    map_mode: str = "magnitude",
    normalization: str = "l2",
    cluster_forming_alpha: float = 0.001,
    alpha: float = 0.05,
    connectivity: int = 6,
    frex_tolerance: float | None = None,
    correction_scope: str = "all_maps",
    n_permutations: int | str = "auto",
    random_seed: int = 0,
    plot_results: bool = True,
    show: bool = True,
) -> FREQNESSSpatialComparisonResult | FREQNESSSpatialFactorialComparisonResult:
    """Compare spatial patterns with pairwise or factorial cluster inference.

    Existing calls pass ``FREQ_A, FREQ_B, MNI``. Generalized designs pass a
    nested grid as the first argument and ``MNI`` as the second; grid axes are
    experimental factors and grid cells contain their level combinations.

    Complete all-within ANOVA grids use orthonormal participant contrasts
    with effect-specific error terms (nominal df q, q*(n-1)). A common sign
    per participant is shared across all contrasts, effects and maps.
    Inference requires independent participants and joint central symmetry
    of null contrast fields after the original cell maps are prepared.
    GG-adjusted F p values form clusters and are recomputed in each draw;
    corrected p values use the maximum cluster mass over the whole family.
    F clusters are non-directional. Mixed designs retain legacy residual
    permutations, require within-block exchangeability, and issue a warning.

    ``n_permutations="auto"`` enumerates up to 65536 participant sign patterns
    (16384 for n=14), otherwise using 10000 draws; ``"exact"`` forces
    enumeration up to the existing safety cap. No new settings are required.
    Effect outputs add ``ss_effect``, ``ss_error``, ``epsilon_gg``, ``df1_gg``,
    ``df2_gg`` and ``contrast_matrix`` (cells x q). Map fields have shape
    (voxels, selected frequencies, selected components). ``df1/df2`` remain
    nominal; ``p_uncorrected`` is the GG approximation, not a permutation p.
    ``contrast_cell_levels`` gives one-based levels for contrast-matrix rows;
    ``inference_method``, ``permutation_unit`` and ``permutation_assumption``
    document the inference used.
    """
    if is_factorial_grid(FREQ_A):
        coordinates = MNI if MNI is not None else FREQ_B
        if coordinates is None:
            raise TypeError("A factorial FREQ grid must be followed by MNI")
        return _factorial_spatial_comparison(
            FREQ_A,
            coordinates,
            frex2compare=frex2compare,
            comps2compare=comps2compare,
            design=design,
            test=test,
            condition_labels=condition_labels,
            factor_names=factor_names,
            level_labels=level_labels,
            subject_ids=subject_ids,
            effects=effects,
            map_mode=map_mode,
            normalization=normalization,
            cluster_forming_alpha=cluster_forming_alpha,
            alpha=alpha,
            connectivity=connectivity,
            frex_tolerance=frex_tolerance,
            correction_scope=correction_scope,
            n_permutations=n_permutations,
            random_seed=random_seed,
            plot_results=plot_results,
            show=show,
        )
    if FREQ_B is None or MNI is None:
        raise TypeError("FREQ_A, FREQ_B, and MNI are required")
    if not isinstance(design, str):
        raise ValueError("Pairwise design must be 'within' or 'between'")
    labels = ("A", "B") if condition_labels is None else condition_labels
    return _pairwise_spatial_comparison(
        FREQ_A,
        FREQ_B,
        MNI,
        frex2compare=frex2compare,
        comps2compare=comps2compare,
        design=design,
        test=test,
        condition_labels=labels,
        map_mode=map_mode,
        normalization=normalization,
        cluster_forming_alpha=cluster_forming_alpha,
        alpha=alpha,
        connectivity=connectivity,
        frex_tolerance=frex_tolerance,
        correction_scope=correction_scope,
        n_permutations=n_permutations,
        random_seed=random_seed,
        plot_results=plot_results,
        show=show,
    )


spatial_comparison = FREQNESS_SpatialComparison


__all__ = [
    "FREQNESS_SpatialComparison",
    "FREQNESSSpatialCluster",
    "FREQNESSSpatialComparisonResult",
    "FREQNESSSpatialFactorialCluster",
    "FREQNESSSpatialFactorialComparisonResult",
    "FREQNESSSpatialFactorialEffect",
    "spatial_comparison",
]
