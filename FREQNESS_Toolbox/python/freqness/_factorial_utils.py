"""Shared validation and modelling for factorial comparison analyses.

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
from itertools import combinations
from typing import Any, Mapping, Sequence
import warnings

import numpy as np
from numpy.typing import NDArray
from scipy import stats


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True, slots=True)
class FactorialEffect:
    """Definition of one main effect or interaction."""

    name: str
    key: str
    factors: tuple[int, ...]


@dataclass(slots=True)
class FactorialModel:
    """Internal model retained for design-aware permutations."""

    kind: str
    Y: FloatArray
    X_full: FloatArray | None = None
    X_reduced: FloatArray | None = None
    permutation_mode: str = ""
    block_index: IntArray | None = None
    group_levels: IntArray | None = None
    nlevels: int = 0


def _is_container(value: Any) -> bool:
    """Return true for mapping or attribute-based result structures."""
    return isinstance(value, Mapping) or any(
        hasattr(value, name)
        for name in (
            "pats",
            "evals",
            "H2",
            "ED",
            "amplitude_raw",
            "stability_index",
        )
    )


def is_factorial_grid(value: Any) -> bool:
    """Identify nested structure grids without confusing numeric landscapes."""
    if isinstance(value, np.ndarray) and value.dtype == object:
        return value.size >= 2 and all(
            _is_container(item) or isinstance(item, np.ndarray)
            for item in value.flat
        )
    if not isinstance(value, (list, tuple)) or len(value) < 2:
        return False

    def leaves(item: Any) -> list[Any]:
        if _is_container(item) or isinstance(item, np.ndarray):
            return [item]
        if isinstance(item, (list, tuple)):
            output: list[Any] = []
            for child in item:
                output.extend(leaves(child))
            return output
        return [item]

    items = leaves(value)
    return bool(items) and all(_is_container(item) for item in items)


def normalize_grid(value: Any) -> tuple[NDArray[np.object_], tuple[int, ...], IntArray]:
    """Return a rectangular object grid and its one-based factor levels."""
    if isinstance(value, np.ndarray) and value.dtype == object:
        grid = value.copy()
    else:
        try:
            grid = np.asarray(value, dtype=object)
        except (TypeError, ValueError) as error:
            raise ValueError("Factorial input must be a rectangular nested grid") from error
    if grid.size < 2 or grid.ndim == 0 or any(item is None for item in grid.flat):
        raise ValueError("Factorial input must contain at least two non-empty cells")
    active_axes = tuple(axis for axis, size in enumerate(grid.shape) if size > 1)
    if not active_axes:
        raise ValueError("Factorial input must contain at least two cells")
    factor_sizes = tuple(int(grid.shape[axis]) for axis in active_axes)
    grid = np.squeeze(grid, axis=tuple(
        axis for axis, size in enumerate(grid.shape) if size == 1
    ))
    if grid.ndim == 0:
        raise ValueError("Factorial input must contain at least two cells")
    if grid.shape != factor_sizes:
        grid = np.reshape(grid, factor_sizes)
    cell_levels = np.asarray(
        [tuple(index + 1 for index in indices) for indices in np.ndindex(grid.shape)],
        dtype=np.int64,
    )
    return grid, factor_sizes, cell_levels


def validate_factor_names(
    factor_names: Sequence[str] | None,
    nfactors: int,
) -> tuple[str, ...]:
    """Validate or generate factor names."""
    names = (
        tuple(f"Factor{index}" for index in range(1, nfactors + 1))
        if factor_names is None
        else tuple(factor_names)
    )
    if len(names) != nfactors or not all(
        isinstance(name, str) and name.strip() for name in names
    ):
        raise ValueError("factor_names must contain one non-empty name per factor")
    keys = tuple(valid_key(name) for name in names)
    if len({key.lower() for key in keys}) != nfactors:
        raise ValueError("factor_names must remain unique as result keys")
    return names


def validate_level_labels(
    level_labels: Sequence[Sequence[str]] | None,
    factor_names: tuple[str, ...],
    factor_sizes: tuple[int, ...],
    pairwise_labels: Sequence[str] | None = None,
) -> tuple[tuple[str, ...], ...]:
    """Validate or generate labels for every factor level."""
    if level_labels is None:
        labels = tuple(
            tuple(f"Level{index}" for index in range(1, size + 1))
            for size in factor_sizes
        )
        if pairwise_labels is not None:
            if len(factor_sizes) != 1 or factor_sizes[0] != 2:
                raise ValueError(
                    "Pairwise labels apply only to a one-factor, two-level grid"
                )
            labels = (tuple(pairwise_labels),)
    else:
        labels = tuple(tuple(values) for values in level_labels)
    if len(labels) != len(factor_sizes):
        raise ValueError("level_labels must contain one sequence per factor")
    for index, (values, size, name) in enumerate(
        zip(labels, factor_sizes, factor_names, strict=True),
        start=1,
    ):
        if len(values) != size or not all(
            isinstance(value, str) and value.strip() for value in values
        ):
            raise ValueError(
                f"level_labels[{index - 1}] must contain {size} non-empty "
                f"labels for {name}"
            )
    return labels


def validate_design(
    design: str | Sequence[str],
    nfactors: int,
) -> tuple[str, ...]:
    """Return one canonical within/between declaration per factor."""
    values = (design,) * nfactors if isinstance(design, str) else tuple(design)
    if len(values) != nfactors:
        raise ValueError("design must contain one value per factor")
    normalized: list[str] = []
    for value in values:
        key = str(value).lower()
        if key in {"within", "paired"}:
            normalized.append("within")
        elif key in {"between", "independent"}:
            normalized.append("between")
        else:
            raise ValueError("design values must be 'within' or 'between'")
    return tuple(normalized)


def validate_test(
    test: str | None,
    design: tuple[str, ...],
    factor_sizes: tuple[int, ...],
    *,
    pairwise_default: str,
) -> tuple[str, str]:
    """Resolve test/design combinations with explicit user warnings."""
    requested = "auto" if test is None else str(test).lower()
    supported = {
        "auto",
        "ttest",
        "signrank",
        "ranksum",
        "anova",
        "friedman",
        "kruskalwallis",
        "permutation_anova",
    }
    if requested not in supported:
        raise ValueError(f"test must be one of {sorted(supported)}")
    resolved = requested
    one_pair = len(factor_sizes) == 1 and factor_sizes[0] == 2
    if resolved == "auto":
        if one_pair:
            if pairwise_default == "rank":
                resolved = "signrank" if design[0] == "within" else "ranksum"
            else:
                resolved = pairwise_default
        else:
            resolved = "anova"
        return resolved, requested

    if one_pair:
        if design[0] == "within" and resolved == "ranksum":
            warnings.warn(
                "test='ranksum' is incompatible with a within-subject factor. "
                "Using test='signrank' instead.",
                UserWarning,
                stacklevel=3,
            )
            resolved = "signrank"
        elif design[0] == "between" and resolved == "signrank":
            warnings.warn(
                "test='signrank' is incompatible with a between-subject factor. "
                "Using test='ranksum' instead.",
                UserWarning,
                stacklevel=3,
            )
            resolved = "ranksum"
        elif resolved == "friedman":
            warnings.warn(
                "Two-level Friedman analysis is equivalent to signrank here.",
                UserWarning,
                stacklevel=3,
            )
            resolved = "signrank"
        elif resolved == "kruskalwallis":
            warnings.warn(
                "Two-level Kruskal-Wallis analysis is equivalent to ranksum here.",
                UserWarning,
                stacklevel=3,
            )
            resolved = "ranksum"
        return resolved, requested

    if len(factor_sizes) == 1:
        replacements = {
            "ttest": "anova",
            "signrank": "friedman",
            "ranksum": "kruskalwallis",
        }
        if resolved in replacements:
            replacement = replacements[resolved]
            warnings.warn(
                f"{resolved} supports two levels only. Using {replacement} instead.",
                UserWarning,
                stacklevel=3,
            )
            resolved = replacement
        if design[0] == "within" and resolved == "kruskalwallis":
            warnings.warn(
                "kruskalwallis is an independent-samples test. Using friedman "
                "for this within-subject factor instead.",
                UserWarning,
                stacklevel=3,
            )
            resolved = "friedman"
        elif design[0] == "between" and resolved == "friedman":
            warnings.warn(
                "friedman is a repeated-measures test. Using kruskalwallis "
                "for this between-subject factor instead.",
                UserWarning,
                stacklevel=3,
            )
            resolved = "kruskalwallis"
    elif resolved in {"signrank", "ranksum", "friedman", "kruskalwallis"}:
        warnings.warn(
            "The selected rank test does not model factorial interactions. "
            "Using permutation_anova instead.",
            UserWarning,
            stacklevel=3,
        )
        resolved = "permutation_anova"
    elif resolved == "ttest":
        warnings.warn(
            "ttest does not model factorial effects. Using anova instead.",
            UserWarning,
            stacklevel=3,
        )
        resolved = "anova"
    return resolved, requested


def valid_key(value: str) -> str:
    """Convert a user label to a stable Python result key."""
    output = "".join(character if character.isalnum() else "_" for character in value)
    output = "_".join(part for part in output.split("_") if part)
    if not output:
        output = "Effect"
    if output[0].isdigit():
        output = f"Effect_{output}"
    return output


def effect_definitions(
    factor_names: tuple[str, ...],
    requested: str | Sequence[str] = "all",
) -> tuple[FactorialEffect, ...]:
    """Define requested main effects and interactions."""
    effects: list[FactorialEffect] = []
    for length in range(1, len(factor_names) + 1):
        for factors in combinations(range(len(factor_names)), length):
            names = tuple(factor_names[index] for index in factors)
            effects.append(
                FactorialEffect(
                    name=" x ".join(names),
                    key="_x_".join(valid_key(name) for name in names),
                    factors=factors,
                )
            )
    if isinstance(requested, str) and requested.lower() == "all":
        return tuple(effects)
    values = (requested,) if isinstance(requested, str) else tuple(requested)
    if not values or not all(isinstance(value, str) and value for value in values):
        raise ValueError("effects must be 'all' or a non-empty sequence")
    selected: list[FactorialEffect] = []
    for value in values:
        normalized = valid_key(value.replace("*", "_x_").replace(":", "_x_"))
        matches = [
            effect
            for effect in effects
            if value.lower() in {effect.name.lower(), effect.key.lower()}
            or normalized.lower() == effect.key.lower()
        ]
        if not matches:
            raise ValueError(f"Unknown factorial effect: {value}")
        if matches[0] not in selected:
            selected.append(matches[0])
    return tuple(selected)


def _subject_grid(
    subject_ids: Any,
    shape: tuple[int, ...],
) -> NDArray[np.object_]:
    """Read nested participant-ID sequences without expanding their contents."""
    output = np.empty(shape, dtype=object)

    def assign(value: Any, depth: int, prefix: tuple[int, ...]) -> None:
        if depth == len(shape):
            output[prefix] = value
            return
        try:
            items = list(value)
        except TypeError as error:
            raise ValueError("subject_ids must match the factorial grid") from error
        if len(items) != shape[depth]:
            raise ValueError("subject_ids must match the factorial grid")
        for index, item in enumerate(items):
            assign(item, depth + 1, prefix + (index,))

    assign(subject_ids, 0, ())
    return output


def match_subjects(
    participant_counts: IntArray,
    cell_levels: IntArray,
    design: tuple[str, ...],
    subject_ids: Any,
) -> tuple[NDArray[np.object_], NDArray[np.object_]]:
    """Match repeated observations by identifier and return reorder indices."""
    counts = np.asarray(participant_counts, dtype=np.int64)
    orders = np.empty(counts.shape, dtype=object)
    identifiers = np.empty(counts.shape, dtype=object)
    if subject_ids is None:
        for index in np.ndindex(counts.shape):
            count = int(counts[index])
            orders[index] = np.arange(count, dtype=np.int64)
            identifiers[index] = tuple(f"S{value}" for value in range(1, count + 1))
        if "within" in design:
            warnings.warn(
                "subject_ids were not supplied. Repeated observations are assumed "
                "to have identical participant order within each between-subject "
                "stratum.",
                UserWarning,
                stacklevel=3,
            )
    else:
        supplied = _subject_grid(subject_ids, counts.shape)
        for index in np.ndindex(counts.shape):
            values = tuple(str(value) for value in supplied[index])
            if len(values) != int(counts[index]) or any(not value for value in values):
                raise ValueError(
                    f"subject_ids{index} must contain one identifier per participant"
                )
            if len(set(values)) != len(values):
                raise ValueError(f"subject_ids{index} must contain unique identifiers")
            identifiers[index] = values
            orders[index] = np.arange(len(values), dtype=np.int64)

    within = tuple(index for index, value in enumerate(design) if value == "within")
    between = tuple(index for index, value in enumerate(design) if value == "between")
    if not within:
        return orders, identifiers
    if between:
        _, strata = np.unique(cell_levels[:, between], axis=0, return_inverse=True)
    else:
        strata = np.zeros(cell_levels.shape[0], dtype=np.int64)
    flat_ids = list(identifiers.flat)
    flat_orders = list(orders.flat)
    for stratum in np.unique(strata):
        cells = np.flatnonzero(strata == stratum)
        reference = flat_ids[int(cells[0])]
        for cell in cells:
            current = flat_ids[int(cell)]
            if len(current) != len(reference) or set(current) != set(reference):
                raise ValueError(
                    "Within-subject cells in each between-subject stratum must "
                    "contain the same participant identifiers"
                )
            lookup = {identifier: index for index, identifier in enumerate(current)}
            flat_orders[int(cell)] = np.asarray(
                [lookup[identifier] for identifier in reference], dtype=np.int64
            )
            flat_ids[int(cell)] = reference
    for flat_index, index in enumerate(np.ndindex(counts.shape)):
        orders[index] = flat_orders[flat_index]
        identifiers[index] = flat_ids[flat_index]
    return orders, identifiers


def long_data(
    values: NDArray[np.object_],
    cell_levels: IntArray,
    design: tuple[str, ...],
) -> tuple[FloatArray, IntArray, IntArray, NDArray[np.object_], IntArray]:
    """Assemble participants as rows and maps/landscapes as columns."""
    flat_values = list(values.flat)
    counts = np.asarray([value.shape[-1] for value in flat_values], dtype=np.int64)
    ncolumns = int(np.prod(flat_values[0].shape[:-1]))
    Y = np.full((int(np.sum(counts)), ncolumns), np.nan)
    levels = np.zeros((Y.shape[0], cell_levels.shape[1]), dtype=np.int64)
    subject_index = np.zeros(Y.shape[0], dtype=np.int64)
    cell_means = np.empty(values.shape, dtype=object)

    within = tuple(index for index, value in enumerate(design) if value == "within")
    between = tuple(index for index, value in enumerate(design) if value == "between")
    if not within:
        strata = np.arange(values.size, dtype=np.int64)
    elif between:
        _, strata = np.unique(cell_levels[:, between], axis=0, return_inverse=True)
    else:
        strata = np.zeros(values.size, dtype=np.int64)
    offsets: dict[int, int] = {}
    next_subject = 0
    if within:
        for stratum in np.unique(strata):
            first = int(np.flatnonzero(strata == stratum)[0])
            offsets[int(stratum)] = next_subject
            next_subject += int(counts[first])

    row = 0
    for flat_index, (index, current) in enumerate(
        zip(np.ndindex(values.shape), flat_values, strict=True)
    ):
        count = current.shape[-1]
        rows = slice(row, row + count)
        Y[rows] = np.moveaxis(current, -1, 0).reshape(count, -1)
        levels[rows] = cell_levels[flat_index]
        cell_means[index] = np.nanmean(current, axis=-1)
        if within:
            subject_index[rows] = offsets[int(strata[flat_index])] + np.arange(count)
        else:
            subject_index[rows] = next_subject + np.arange(count)
            next_subject += count
        row += count
    return Y, levels, subject_index, cell_means, counts.reshape(values.shape)


def _effect_coding(levels: IntArray, nlevels: int) -> FloatArray:
    coding = np.zeros((levels.size, nlevels - 1))
    for level in range(1, nlevels):
        coding[levels == level, level - 1] = 1.0
    coding[levels == nlevels] = -1.0
    return coding


def _numeric_effects(nfactors: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        factors
        for length in range(1, nfactors + 1)
        for factors in combinations(range(nfactors), length)
    )


def _interaction(matrices: Sequence[FloatArray]) -> FloatArray:
    output = np.ones((matrices[0].shape[0], 1))
    for current in matrices:
        output = np.column_stack(
            [output[:, left] * current[:, right]
             for left in range(output.shape[1])
             for right in range(current.shape[1])]
        )
    return output


def _subject_average(
    Y: FloatArray,
    levels: IntArray,
    subject_index: IntArray,
) -> tuple[FloatArray, IntArray, IntArray]:
    subjects = np.unique(subject_index)
    averaged = np.full((subjects.size, Y.shape[1]), np.nan)
    subject_levels = np.zeros((subjects.size, levels.shape[1]), dtype=np.int64)
    for row, subject in enumerate(subjects):
        selected = subject_index == subject
        averaged[row] = np.nanmean(Y[selected], axis=0)
        subject_levels[row] = levels[np.flatnonzero(selected)[0]]
    return averaged, subject_levels, subjects


def _partial_f_single(
    Y: FloatArray,
    X_full: FloatArray,
    X_reduced: FloatArray,
) -> tuple[float, float, float, float, float]:
    if Y.size < 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    rank_full = np.linalg.matrix_rank(X_full)
    rank_reduced = np.linalg.matrix_rank(X_reduced)
    df1 = rank_full - rank_reduced
    df2 = X_full.shape[0] - rank_full
    if df1 < 1 or df2 < 1:
        return np.nan, np.nan, np.nan, float(df1), float(df2)
    residual_full = Y - X_full @ (np.linalg.pinv(X_full) @ Y)
    residual_reduced = Y - X_reduced @ (np.linalg.pinv(X_reduced) @ Y)
    sse_full = float(residual_full @ residual_full)
    sse_reduced = float(residual_reduced @ residual_reduced)
    numerator = max(sse_reduced - sse_full, 0.0) / df1
    denominator = sse_full / df2
    statistic = numerator / denominator if denominator > np.finfo(float).eps else (
        np.inf if numerator > np.finfo(float).eps else 0.0
    )
    p_value = float(stats.f.sf(statistic, df1, df2))
    eta_denominator = statistic * df1 + df2
    effect_size = statistic * df1 / eta_denominator if eta_denominator > 0 else 0.0
    return float(statistic), p_value, float(effect_size), float(df1), float(df2)


def _partial_f_columns(
    Y: FloatArray,
    X_full: FloatArray,
    X_reduced: FloatArray,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    outputs = [np.full(Y.shape[1], np.nan) for _ in range(5)]
    for column in range(Y.shape[1]):
        valid = np.isfinite(Y[:, column])
        values = _partial_f_single(
            Y[valid, column], X_full[valid], X_reduced[valid]
        )
        for output, value in zip(outputs, values, strict=True):
            output[column] = value
    return tuple(outputs)  # type: ignore[return-value]


def anova_statistics(
    Y: FloatArray,
    levels: IntArray,
    subject_index: IntArray,
    factor_sizes: tuple[int, ...],
    design: tuple[str, ...],
    effects: tuple[FactorialEffect, ...],
) -> tuple[list[dict[str, Any]], list[FactorialModel]]:
    """Compute partial-F tests for requested factorial effects."""
    numeric_effects = _numeric_effects(len(factor_sizes))
    codings = tuple(
        _effect_coding(levels[:, index], factor_sizes[index])
        for index in range(len(factor_sizes))
    )
    matrices = {
        factors: _interaction(tuple(codings[index] for index in factors))
        for factors in numeric_effects
    }
    results: list[dict[str, Any]] = []
    models: list[FactorialModel] = []
    for effect in effects:
        between_only = all(design[index] == "between" for index in effect.factors)
        if between_only:
            current_Y, current_levels, subjects = _subject_average(
                Y, levels, subject_index
            )
            current_codings = tuple(
                _effect_coding(current_levels[:, index], factor_sizes[index])
                for index in range(len(factor_sizes))
            )
            included = tuple(
                factors
                for factors in numeric_effects
                if all(design[index] == "between" for index in factors)
            )
            X_full = np.ones((current_Y.shape[0], 1))
            target_columns: list[int] = []
            for factors in included:
                columns = _interaction(
                    tuple(current_codings[index] for index in factors)
                )
                start = X_full.shape[1]
                X_full = np.column_stack((X_full, columns))
                if factors == effect.factors:
                    target_columns = list(range(start, X_full.shape[1]))
            mode = "between"
            blocks = subjects
        else:
            current_Y = Y
            subjects = np.unique(subject_index)
            X_full = np.column_stack(
                [(subject_index == subject).astype(float) for subject in subjects]
            )
            target_columns = []
            for factors in numeric_effects:
                if not any(design[index] == "within" for index in factors):
                    continue
                columns = matrices[factors]
                start = X_full.shape[1]
                X_full = np.column_stack((X_full, columns))
                if factors == effect.factors:
                    target_columns = list(range(start, X_full.shape[1]))
            mode = "within"
            blocks = subject_index
        if not target_columns:
            raise ValueError(f"The requested effect {effect.name} is not estimable")
        X_reduced = np.delete(X_full, target_columns, axis=1)
        statistic, p_values, effect_size, df1, df2 = _partial_f_columns(
            current_Y, X_full, X_reduced
        )
        results.append(
            {
                "definition": effect,
                "statistic": statistic,
                "statistic_name": "F",
                "effect_size": effect_size,
                "effect_size_name": "partial eta squared",
                "p_uncorrected": p_values,
                "df1": df1,
                "df2": df2,
            }
        )
        models.append(
            FactorialModel(
                kind="anova",
                Y=current_Y,
                X_full=X_full,
                X_reduced=X_reduced,
                permutation_mode=mode,
                block_index=np.asarray(blocks, dtype=np.int64),
            )
        )
    return results, models


def _rankdata_with_ties(values: FloatArray) -> tuple[FloatArray, float]:
    ranks = stats.rankdata(values, method="average")
    _, counts = np.unique(values, return_counts=True)
    ties = float(np.sum(counts[counts > 1] ** 3 - counts[counts > 1]))
    return np.asarray(ranks, dtype=float), ties


def _kruskal_columns(
    Y: FloatArray,
    groups: IntArray,
    nlevels: int,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    statistic = np.full(Y.shape[1], np.nan)
    p_values = np.full(Y.shape[1], np.nan)
    effect_size = np.full(Y.shape[1], np.nan)
    df1 = np.full(Y.shape[1], nlevels - 1.0)
    df2 = np.full(Y.shape[1], np.inf)
    for column in range(Y.shape[1]):
        valid = np.isfinite(Y[:, column])
        values = Y[valid, column]
        current_groups = groups[valid]
        counts = np.bincount(current_groups, minlength=nlevels + 1)[1:]
        if values.size <= nlevels or np.any(counts == 0):
            continue
        ranks, ties = _rankdata_with_ties(values)
        rank_sums = np.asarray(
            [np.sum(ranks[current_groups == level]) for level in range(1, nlevels + 1)]
        )
        H = (
            12.0 / (values.size * (values.size + 1.0))
            * np.sum(rank_sums**2 / counts)
            - 3.0 * (values.size + 1.0)
        )
        correction = 1.0 - ties / (values.size**3 - values.size)
        H = H / correction if correction > np.finfo(float).eps else 0.0
        statistic[column] = max(H, 0.0)
        p_values[column] = stats.chi2.sf(statistic[column], nlevels - 1)
        effect_size[column] = max(
            (statistic[column] - nlevels + 1.0) / (values.size - nlevels), 0.0
        )
    return statistic, p_values, effect_size, df1, df2


def _friedman_columns(
    repeated: FloatArray,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    nlevels = repeated.shape[1]
    statistic = np.full(repeated.shape[2], np.nan)
    p_values = np.full(repeated.shape[2], np.nan)
    effect_size = np.full(repeated.shape[2], np.nan)
    df1 = np.full(repeated.shape[2], nlevels - 1.0)
    df2 = np.full(repeated.shape[2], np.inf)
    for column in range(repeated.shape[2]):
        current = repeated[:, :, column]
        current = current[np.all(np.isfinite(current), axis=1)]
        nsubjects = current.shape[0]
        if nsubjects < 2:
            continue
        rank_sums = np.zeros(nlevels)
        ties = 0.0
        for row in current:
            ranks, current_ties = _rankdata_with_ties(row)
            rank_sums += ranks
            ties += current_ties
        value = (
            12.0 / (nsubjects * nlevels * (nlevels + 1.0))
            * np.sum(rank_sums**2)
            - 3.0 * nsubjects * (nlevels + 1.0)
        )
        correction = 1.0 - ties / (
            nsubjects * nlevels * (nlevels**2 - 1.0)
        )
        value = value / correction if correction > np.finfo(float).eps else 0.0
        statistic[column] = max(value, 0.0)
        p_values[column] = stats.chi2.sf(statistic[column], nlevels - 1)
        effect_size[column] = statistic[column] / (nsubjects * (nlevels - 1.0))
    return statistic, p_values, effect_size, df1, df2


def rank_statistics(
    Y: FloatArray,
    levels: IntArray,
    subject_index: IntArray,
    test: str,
    effects: tuple[FactorialEffect, ...],
) -> tuple[list[dict[str, Any]], list[FactorialModel]]:
    """Compute one-factor Friedman or Kruskal-Wallis statistics."""
    if levels.shape[1] != 1 or len(effects) != 1:
        raise ValueError("friedman and kruskalwallis support one effect only")
    nlevels = int(np.max(levels[:, 0]))
    if test == "kruskalwallis":
        output = _kruskal_columns(Y, levels[:, 0], nlevels)
        statistic_name = "H"
        effect_name = "epsilon squared"
        model = FactorialModel(
            kind="kruskalwallis",
            Y=Y,
            group_levels=levels[:, 0],
            nlevels=nlevels,
        )
    else:
        subjects = np.unique(subject_index)
        repeated = np.full((subjects.size, nlevels, Y.shape[1]), np.nan)
        for subject_row, subject in enumerate(subjects):
            for level in range(1, nlevels + 1):
                rows = np.flatnonzero(
                    (subject_index == subject) & (levels[:, 0] == level)
                )
                if rows.size != 1:
                    raise ValueError(
                        "Friedman analysis requires one observation per "
                        "participant at every factor level"
                    )
                repeated[subject_row, level - 1] = Y[rows[0]]
        output = _friedman_columns(repeated)
        statistic_name = "chi-square"
        effect_name = "Kendall's W"
        model = FactorialModel(kind="friedman", Y=repeated, nlevels=nlevels)
    statistic, p_values, effect_size, df1, df2 = output
    result = {
        "definition": effects[0],
        "statistic": statistic,
        "statistic_name": statistic_name,
        "effect_size": effect_size,
        "effect_size_name": effect_name,
        "p_uncorrected": p_values,
        "df1": df1,
        "df2": df2,
    }
    return [result], [model]


def permutation_count(value: int | str) -> int:
    """Validate the Monte Carlo count for generalized designs."""
    if isinstance(value, str):
        if value.lower() != "auto":
            raise ValueError(
                "Generalized factorial designs support an integer or 'auto'"
            )
        return 10_000
    if isinstance(value, (bool, np.bool_)):
        raise ValueError("n_permutations must be a positive integer")
    count = int(value)
    if count < 1 or count != value:
        raise ValueError("n_permutations must be a positive integer")
    return count


def permute_model(
    model: FactorialModel,
    rng: np.random.Generator,
) -> tuple[FloatArray, FloatArray]:
    """Generate one design-aware permutation of an omnibus statistic."""
    if model.kind == "anova":
        if np.any(~np.isfinite(model.Y)):
            raise ValueError(
                "Spatial factorial permutations require finite selected patterns"
            )
        assert model.X_full is not None and model.X_reduced is not None
        residual = model.Y - model.X_reduced @ (
            np.linalg.pinv(model.X_reduced) @ model.Y
        )
        fitted = model.Y - residual
        if model.permutation_mode == "between":
            order = rng.permutation(model.Y.shape[0])
        else:
            assert model.block_index is not None
            order = np.arange(model.Y.shape[0])
            for subject in np.unique(model.block_index):
                rows = np.flatnonzero(model.block_index == subject)
                order[rows] = rng.permutation(rows)
        permuted = fitted + residual[order]
        return _partial_f_columns(permuted, model.X_full, model.X_reduced)[:2]
    if model.kind == "kruskalwallis":
        assert model.group_levels is not None
        groups = rng.permutation(model.group_levels)
        return _kruskal_columns(model.Y, groups, model.nlevels)[:2]
    if model.kind == "friedman":
        permuted = model.Y.copy()
        for subject in range(permuted.shape[0]):
            permuted[subject] = permuted[subject, rng.permutation(model.nlevels)]
        return _friedman_columns(permuted)[:2]
    raise RuntimeError(f"Unknown factorial model: {model.kind}")


def permutation_anova_p_values(
    model: FactorialModel,
    observed: FloatArray,
    n_permutations: int,
    rng: np.random.Generator,
) -> FloatArray:
    """Compute frequency-wise empirical p values using Freedman-Lane residuals."""
    assert model.X_full is not None and model.X_reduced is not None
    p_values = np.full(observed.shape, np.nan)
    for column, observed_value in enumerate(observed):
        if np.isnan(observed_value):
            continue
        valid = np.isfinite(model.Y[:, column])
        Y = model.Y[valid, column]
        X_full = model.X_full[valid]
        X_reduced = model.X_reduced[valid]
        blocks = None if model.block_index is None else model.block_index[valid]
        _, _, _, df1, df2 = _partial_f_single(Y, X_full, X_reduced)
        if not np.isfinite(df1) or not np.isfinite(df2):
            continue
        residual = Y - X_reduced @ (np.linalg.pinv(X_reduced) @ Y)
        fitted = Y - residual
        count = 0
        for _ in range(n_permutations):
            if model.permutation_mode == "between":
                order = rng.permutation(Y.size)
            else:
                assert blocks is not None
                order = np.arange(Y.size)
                for subject in np.unique(blocks):
                    rows = np.flatnonzero(blocks == subject)
                    order[rows] = rng.permutation(rows)
            permuted = fitted + residual[order]
            statistic = _partial_f_single(permuted, X_full, X_reduced)[0]
            count += int(statistic >= observed_value)
        p_values[column] = (count + 1.0) / (n_permutations + 1.0)
    return p_values


__all__ = [
    "FactorialEffect",
    "FactorialModel",
    "anova_statistics",
    "effect_definitions",
    "is_factorial_grid",
    "long_data",
    "match_subjects",
    "normalize_grid",
    "permutation_anova_p_values",
    "permutation_count",
    "permute_model",
    "rank_statistics",
    "validate_design",
    "validate_factor_names",
    "validate_level_labels",
    "validate_test",
]
