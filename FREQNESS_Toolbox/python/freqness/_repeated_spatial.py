"""Participant-contrast inference for complete within-subject spatial designs.

This module is deliberately separate from the shared landscape/mixed-design
GLM code. Sign flipping requires independent participants and joint symmetry
of the analysed null contrast vectors, including their spatial dependence.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy import stats

from ._factorial_utils import FactorialEffect


@dataclass
class RepeatedSpatialModel:
    scores: NDArray[np.float64]  # participants x contrasts x maps
    second_moment: NDArray[np.float64]  # contrasts x contrasts x maps
    contrast_matrix: NDArray[np.float64]  # factorial cells x contrasts


def helmert_matrix(nlevels: int) -> NDArray[np.float64]:
    """Orthonormal contrasts: each later level minus preceding levels."""
    matrix = np.zeros((nlevels, nlevels - 1))
    for index in range(1, nlevels):
        scale = np.sqrt(index * (index + 1))
        matrix[:index, index - 1] = -1.0 / scale
        matrix[index, index - 1] = index / scale
    return matrix


def repeated_spatial_models(
    maps: NDArray[np.object_],
    cell_levels: NDArray[np.int64],
    factor_sizes: tuple[int, ...],
    definitions: tuple[FactorialEffect, ...],
) -> list[RepeatedSpatialModel]:
    """Project the prepared maps; never normalize a mean or contrast again."""
    # C-order cell traversal agrees with normalize_grid/long_data in Python.
    nsubs = maps.flat[0].shape[-1]
    data = np.stack(
        [value.reshape(-1, nsubs).T for value in maps.flat], axis=1
    )  # participants x cells x maps
    # Remove each participant's grand mean before contrasts for numerical
    # stability. All effect bases are orthogonal to this nuisance term.
    data = data - np.mean(data, axis=1, keepdims=True)
    # Roundoff from projecting an exactly absent effect must not become
    # evidence through division by its equally tiny estimated variance.
    tolerance = (16 * np.finfo(float).eps * maps.size
                 * np.max(np.abs(data), axis=(0, 1)))
    models = []
    for definition in definitions:
        weights = np.ones((maps.size, 1))
        for factori, nlevels in enumerate(factor_sizes):
            coding = (
                helmert_matrix(nlevels)
                if factori in definition.factors
                else np.ones((nlevels, 1)) / np.sqrt(nlevels)
            )
            values = coding[cell_levels[:, factori] - 1]  # public levels are one-based
            weights = (weights[:, :, None] * values[:, None, :]).reshape(maps.size, -1)
        scores = np.einsum("scm,cq->sqm", data, weights)
        scores[np.abs(scores) <= tolerance[None, None, :]] = 0.0
        second = np.einsum("sqm,srm->qrm", scores, scores)
        models.append(RepeatedSpatialModel(scores, second, weights))
    return models


def repeated_spatial_statistic(
    model: RepeatedSpatialModel, signs: NDArray[np.float64] | None = None
) -> dict[str, NDArray[np.float64]]:
    """Effect-specific RM F and GG screening probabilities for one sign draw.

    Only participant means change under sign flips. Retain second moments,
    not a maps-by-permutations tensor. Recompute near-degenerate scatter
    directly to avoid cancellation in Q - n * mean * mean'.
    """
    scores = model.scores
    nsubs, ncontrasts, nmaps = scores.shape
    if signs is None:
        signs = np.ones(nsubs)
    means = np.einsum("s,sqm->qm", signs, scores) / nsubs
    scatter = model.second_moment - nsubs * np.einsum("qm,rm->qrm", means, means)
    ss_error = np.einsum("qqm->m", scatter)
    total = np.einsum("qqm->m", model.second_moment)
    unstable = ss_error <= 64 * np.finfo(float).eps * total
    if np.any(unstable):
        values = scores[:, :, unstable] * signs[:, None, None]
        centered = values - means[:, unstable][None, :, :]
        constant = np.all(np.max(values, axis=0) == np.min(values, axis=0), axis=0)
        centered[:, :, constant] = 0.0
        scatter[:, :, unstable] = np.einsum("sqm,srm->qrm", centered, centered)
        ss_error = np.einsum("qqm->m", scatter)
    ss_effect = nsubs * np.sum(means**2, axis=0)
    statistic = np.divide(
        (nsubs - 1) * ss_effect, ss_error,
        out=np.zeros(nmaps), where=ss_error > 0,
    )
    statistic[(ss_error == 0) & (ss_effect > 0)] = np.inf
    eta = np.divide(
        ss_effect, ss_effect + ss_error,
        out=np.zeros(nmaps), where=(ss_effect + ss_error) > 0,
    )
    epsilon = np.ones(nmaps)
    if ncontrasts > 1:
        frobenius = np.sum(scatter**2, axis=(0, 1))
        np.divide(
            ss_error**2, ncontrasts * frobenius,
            out=epsilon, where=frobenius > 0,
        )
        epsilon = np.clip(epsilon, 1.0 / ncontrasts, 1.0)
    df1 = np.full(nmaps, ncontrasts, dtype=float)
    df2 = np.full(nmaps, ncontrasts * (nsubs - 1), dtype=float)
    return dict(
        statistic=statistic,
        p_uncorrected=stats.f.sf(statistic, df1 * epsilon, df2 * epsilon),
        effect_size=eta, df1=df1, df2=df2, epsilon_gg=epsilon,
        df1_gg=df1 * epsilon, df2_gg=df2 * epsilon,
        ss_effect=ss_effect, ss_error=ss_error,
    )
