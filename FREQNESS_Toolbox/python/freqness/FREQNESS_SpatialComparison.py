"""Paired spatial comparison of FREQ-NESS activation patterns.

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
from scipy import stats
from scipy.spatial import cKDTree


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
    n_participants: int
    mean_difference: FloatArray
    t_values: FloatArray
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


def _prepare_maps(
    maps_a: FloatArray,
    maps_b: FloatArray,
    map_mode: str,
    normalization: str,
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
                for subi in range(prepared_a.shape[3]):
                    if np.dot(
                        prepared_a[:, freqi, compi, subi],
                        prepared_b[:, freqi, compi, subi],
                    ) < 0:
                        prepared_b[:, freqi, compi, subi] *= -1
                        flips_b[freqi, compi, subi] *= -1
                paired_mean = 0.5 * (
                    prepared_a[:, freqi, compi, :]
                    + prepared_b[:, freqi, compi, :]
                )
                reference = np.linalg.svd(paired_mean, full_matrices=False)[0][:, 0]
                peak = int(np.argmax(np.abs(reference)))
                if reference[peak] < 0:
                    reference *= -1
                for subi in range(prepared_a.shape[3]):
                    if np.dot(paired_mean[:, subi], reference) < 0:
                        prepared_a[:, freqi, compi, subi] *= -1
                        prepared_b[:, freqi, compi, subi] *= -1
                        flips_a[freqi, compi, subi] *= -1
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
    t_values: FloatArray,
    threshold: float,
    adjacency: list[IntArray],
) -> float:
    maximum = 0.0
    for freqi in range(t_values.shape[1]):
        for compi in range(t_values.shape[2]):
            spatial_map = t_values[:, freqi, compi]
            for mask in (spatial_map >= threshold, spatial_map <= -threshold):
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
            figure = plt.figure(figsize=(12, 5))
            for panel, (values, mask, title) in enumerate(
                (
                    (difference, np.ones(difference.size, dtype=bool), "Mean difference"),
                    (difference, significant, "Cluster-corrected difference"),
                ),
                start=1,
            ):
                axis = figure.add_subplot(1, 2, panel, projection="3d")
                axis.scatter(*coordinates.T, s=6, c="black", alpha=0.18, linewidths=0)
                if np.any(mask):
                    points = axis.scatter(
                        *coordinates[mask].T,
                        c=values[mask],
                        cmap="coolwarm",
                        vmin=-limit,
                        vmax=limit,
                        s=24,
                        alpha=0.9,
                        linewidths=0,
                    )
                    figure.colorbar(points, ax=axis, shrink=0.65, pad=0.08)
                axis.set_title(title)
                axis.set_xlabel("X (mm)")
                axis.set_ylabel("Y (mm)")
                axis.set_zlabel("Z (mm)")
                axis.set_box_aspect(np.ptp(coordinates, axis=0))
            figure.suptitle(
                f"{result.condition_labels[1]} - {result.condition_labels[0]} | "
                f"{frequency:g} Hz | Component #{component}"
            )
            figure.tight_layout()
    if show:
        plt.show()


def FREQNESS_SpatialComparison(
    FREQ_A: Any,
    FREQ_B: Any,
    MNI: ArrayLike,
    *,
    frex2compare: ArrayLike | str,
    comps2compare: ArrayLike | None = None,
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
    """Compare paired FREQ-NESS spatial patterns with cluster correction.

    The participant order must match across ``FREQ_A`` and ``FREQ_B``.
    Component numbers retain the MATLAB convention and are one-based, while
    ``idx_frex`` is zero-based for direct NumPy indexing.
    """
    patterns_a, frequency_axis_a = _canonical_freq(FREQ_A, "FREQ_A")
    patterns_b, frequency_axis_b = _canonical_freq(FREQ_B, "FREQ_B")
    if patterns_a.shape != patterns_b.shape:
        raise ValueError("FREQ_A.pats and FREQ_B.pats must have matching dimensions")
    if not np.allclose(frequency_axis_a, frequency_axis_b, atol=1e-10, rtol=0.0):
        raise ValueError("FREQ_A.frex and FREQ_B.frex must contain the same grid")
    if patterns_a.shape[3] < 2:
        raise ValueError("At least two paired participants are required")

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
        selected_a, selected_b, map_mode, normalization
    )
    differences = selected_b - selected_a
    mean_difference, t_values, effect_size, p_uncorrected = _paired_statistics(differences)
    adjacency, grid_spacing = _spatial_adjacency(coordinates, connectivity)
    threshold = float(stats.t.ppf(1.0 - cluster_forming_alpha / 2.0, patterns_a.shape[3] - 1))
    signs, exact = _permutation_plan(patterns_a.shape[3], n_permutations, int(random_seed))
    print(
        "\nFREQNESS Spatial Comparison: comparing "
        f"{condition_labels[0]} and {condition_labels[1]} at "
        f"{resolved.size} frequencies and {selected_components.size} components "
        f"across {patterns_a.shape[3]} paired participants."
    )
    print(
        f"Computing {signs.shape[0]} "
        f"{'exact' if exact else 'Monte Carlo'} paired sign-flip permutations..."
    )
    null_maxima = np.zeros(signs.shape[0])
    progress_step = max(1, signs.shape[0] // 10)
    for permutation, permutation_signs in enumerate(signs):
        _, _, permuted_t = _paired_moments(
            differences * permutation_signs.reshape(1, 1, 1, -1)
        )
        null_maxima[permutation] = _maximum_cluster_mass(
            permuted_t, threshold, adjacency
        )
        if (permutation + 1) % progress_step == 0 or permutation + 1 == signs.shape[0]:
            print(f"  Permutation {permutation + 1}/{signs.shape[0]}")

    corrected = np.ones_like(t_values)
    significant = np.zeros_like(t_values, dtype=bool)
    labels = np.zeros_like(t_values, dtype=np.int64)
    cluster_records: list[FREQNESSSpatialCluster] = []
    for freqi, frequency in enumerate(resolved):
        for compi, component in enumerate(selected_components):
            spatial_map = t_values[:, freqi, compi]
            positive_label = 0
            negative_label = 0
            for polarity, mask in (
                (1, spatial_map >= threshold),
                (-1, spatial_map <= -threshold),
            ):
                for members in _clusters(mask, adjacency):
                    mass = float(np.sum(np.abs(spatial_map[members])))
                    if exact:
                        p_corrected = float(np.mean(null_maxima >= mass))
                    else:
                        p_corrected = float((1 + np.count_nonzero(null_maxima >= mass)) / (signs.shape[0] + 1))
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
        design="paired",
        n_participants=patterns_a.shape[3],
        mean_difference=mean_difference,
        t_values=t_values,
        effect_size_dz=effect_size,
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
        n_permutations=signs.shape[0],
        exact_permutations=exact,
        random_seed=int(random_seed),
        sign_flips_a=flips_a,
        sign_flips_b=flips_b,
    )
    if plot_results:
        _plot_results(coordinates, result, bool(show))
    return result


spatial_comparison = FREQNESS_SpatialComparison


__all__ = [
    "FREQNESS_SpatialComparison",
    "FREQNESSSpatialCluster",
    "FREQNESSSpatialComparisonResult",
    "spatial_comparison",
]
