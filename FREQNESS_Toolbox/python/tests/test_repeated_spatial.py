"""Independent repeated-measures references, exact inference, and contracts."""
from itertools import product
import numpy as np
import pytest
from scipy import stats
from scipy.sparse.csgraph import connected_components

from freqness import FREQNESS_SpatialComparison
from freqness._factorial_utils import effect_definitions, normalize_grid
from freqness._repeated_spatial import repeated_spatial_models, repeated_spatial_statistic


def dataset(shape=(2, 4), nsubs=8, nfrex=2, ncomps=2):
    rng = np.random.default_rng(17)
    coordinates = np.array(list(np.ndindex(2, 2, 2))) * 8.0
    grid, ids = np.empty(shape, object), np.empty(shape, object)
    baseline = rng.normal(20, 1, (8, ncomps, nfrex, nsubs))
    for index in np.ndindex(shape):
        values = baseline + rng.normal(0, .4, baseline.shape)
        values[:4] += index[0] * .8
        grid[index] = dict(pats=values, frex=np.arange(nfrex) + 8.0)
        ids[index] = list(range(nsubs))
    return grid, ids, coordinates


def run(grid, ids, coordinates, **kwargs):
    options = dict(frex2compare="all", design="within", test="anova",
                   subject_ids=ids, plot_results=False, n_permutations=7,
                   comps2compare=list(range(1, grid.flat[0]["pats"].shape[1] + 1)))
    options.update(kwargs)
    return FREQNESS_SpatialComparison(grid, coordinates, **options)


def prepared(grid, normalize=False):
    maps = np.empty(grid.shape, object)
    for index in np.ndindex(grid.shape):
        values = np.abs(grid[index]["pats"]).transpose(0, 2, 1, 3)
        if normalize:
            values = values / np.linalg.norm(values, axis=0, keepdims=True)
        maps[index] = values
    return maps


def marginal_projection(data, shape, factors):
    """Independent reference using cell means/centering, not Helmert coding."""
    value = data.reshape((data.shape[0],) + shape + (data.shape[-1],))
    for factori in range(len(shape)):
        mean = value.mean(axis=factori + 1, keepdims=True)
        value = value - mean if factori in factors else np.broadcast_to(mean, value.shape)
    return value.reshape(data.shape)


def reference(effect, q):
    nsubs = effect.shape[0]
    mean = effect.mean(axis=0)
    centered = effect - mean
    ss_effect = nsubs * np.sum(mean**2, axis=0)
    ss_error = np.sum(centered**2, axis=(0, 1))
    statistic = (nsubs - 1) * ss_effect / ss_error
    covariance = np.einsum("scm,sdm->cdm", centered, centered)
    epsilon = np.clip(ss_error**2 / (q * np.sum(covariance**2, axis=(0, 1))), 1/q, 1)
    return statistic, epsilon, ss_effect, ss_error


@pytest.mark.parametrize("shape", [(2, 4), (3,), (2, 3, 2)])
def test_effect_statistics_against_independent_marginal_anova(shape):
    grid, ids, coordinates = dataset(shape)
    result = run(grid, ids, coordinates, normalization="none")
    maps = prepared(grid)
    data = np.stack([v.reshape(-1, 8).T for v in maps.flat], axis=1)
    for effect in result.effects.values():
        q = np.prod([shape[f] - 1 for f in effect.factors])
        expected = reference(marginal_projection(data, shape, effect.factors), q)
        for actual, target in zip((effect.statistic, effect.epsilon_gg,
                                   effect.ss_effect, effect.ss_error), expected):
            np.testing.assert_allclose(actual.ravel(), target, rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(effect.effect_size.ravel(), expected[2]/(expected[2]+expected[3]))
        np.testing.assert_allclose(effect.p_uncorrected.ravel(),
                                   stats.f.sf(expected[0], q*expected[1], q*7*expected[1]))
        assert np.all(effect.df1 == q) and np.all(effect.df2 == q*7)
        assert effect.statistic.shape == (8, 2, 2)
        np.testing.assert_allclose(effect.contrast_matrix.T @ effect.contrast_matrix,
                                   np.eye(q), atol=1e-14)


def test_drug_equals_t_squared_after_original_cell_normalization():
    grid, ids, coordinates = dataset(nsubs=14)
    result = run(grid, ids, coordinates, factor_names=["Drug", "Condition"])
    maps = prepared(grid, normalize=True)
    difference = np.mean(np.stack(list(maps[1])), axis=0) - np.mean(np.stack(list(maps[0])), axis=0)
    expected = stats.ttest_1samp(difference, 0, axis=-1).statistic**2
    np.testing.assert_allclose(result.effects["Drug"].statistic, expected, atol=1e-10)
    assert np.all(result.effects["Drug"].df2 == 13)
    assert np.all(result.effects["Condition"].df2 == 39)
    assert np.all(result.effects["Drug_x_Condition"].df2 == 39)
    assert result.inference_method == "within_subject_contrast_signflip"
    assert result.permutation_unit == "participant"


def test_repeated_condition_counterexample_and_zero_effects():
    grid, ids, coordinates = dataset(nsubs=14, nfrex=1, ncomps=1)
    difference = np.linspace(-1, 2, 14)
    for drug, condition in np.ndindex(grid.shape):
        grid[drug, condition]["pats"][:] = 20 + drug*difference
    result = run(grid, ids, coordinates, normalization="none")
    expected = stats.ttest_1samp(difference, 0).statistic**2
    np.testing.assert_allclose(result.effects["Factor1"].statistic, expected)
    for key in ("Factor2", "Factor1_x_Factor2"):
        np.testing.assert_array_equal(result.effects[key].statistic, 0)
        np.testing.assert_array_equal(result.effects[key].p_uncorrected, 1)


def cluster_masses(statistic, p_values, coordinates, alpha):
    """Independent graph-component reference over each frequency/component."""
    adjacency = np.sum(np.abs(coordinates[:, None] - coordinates), axis=2) == 8
    masses = []
    for mapi in range(statistic.shape[1]):
        selected = np.flatnonzero(p_values[:, mapi] <= alpha)
        if selected.size:
            count, labels = connected_components(adjacency[np.ix_(selected, selected)])
            masses.extend(np.sum(statistic[selected[labels == i], mapi]) for i in range(count))
    return masses


def test_exact_joint_family_against_brute_force_cell_projection():
    grid, ids, coordinates = dataset(nsubs=6)
    result = run(grid, ids, coordinates, normalization="none",
                 n_permutations="exact", cluster_forming_alpha=.2)
    data = np.stack([v.reshape(-1, 6).T for v in prepared(grid).flat], axis=1)
    fields = [(marginal_projection(data, grid.shape, effect.factors),
               int(np.prod([grid.shape[f]-1 for f in effect.factors])))
              for effect in result.effects.values()]
    maxima = []
    for signs in product((-1, 1), repeat=6):
        maximum = 0
        for field, q in fields:
            statistic, epsilon, _, _ = reference(field*np.array(signs)[:, None, None], q)
            p_values = stats.f.sf(statistic, q*epsilon, q*5*epsilon)
            masses = cluster_masses(statistic.reshape(8, -1), p_values.reshape(8, -1), coordinates, .2)
            maximum = max([maximum] + masses)
        maxima.append(maximum)
    assert result.exact_permutations and result.n_permutations == 64
    assert sum(len(effect.clusters) for effect in result.effects.values()) > 0
    for effect in result.effects.values():
        for cluster in effect.clusters:
            expected = np.mean(np.asarray(maxima) >= cluster.mass*(1-1e-12))
            assert cluster.p_corrected == expected


def test_reordering_invariance_and_monte_carlo_reproducibility():
    grid, ids, coordinates = dataset()
    baseline = run(grid, ids, coordinates, random_seed=29)
    rng = np.random.default_rng(44)
    for index in np.ndindex(grid.shape):
        order = rng.permutation(8)
        grid[index]["pats"] = -3.5*grid[index]["pats"][..., order]
        ids[index] = order.tolist()
    changed = run(grid, ids, coordinates, random_seed=29)
    again = run(grid, ids, coordinates, random_seed=29)
    for key in baseline.effects:
        np.testing.assert_allclose(baseline.effects[key].statistic, changed.effects[key].statistic, atol=1e-10)
        # Reordering the reference cell changes which seeded sign is assigned
        # to each person; compare the repeated run for Monte Carlo identity.
        np.testing.assert_array_equal(changed.effects[key].p_cluster_corrected,
                                      again.effects[key].p_cluster_corrected)


def test_constant_contrasts_and_basis_rotation():
    grid, _, _ = dataset(nfrex=1, ncomps=1)
    maps = prepared(grid)
    _, sizes, levels = normalize_grid(grid)
    models = repeated_spatial_models(maps, levels, sizes,
                                     effect_definitions(("A", "B"), "all"))
    model = models[1]
    original = repeated_spatial_statistic(model)
    rotation, _ = np.linalg.qr(np.random.default_rng(7).normal(size=(3, 3)))
    model.scores = np.einsum("sqm,qr->srm", model.scores, rotation)
    model.second_moment = np.einsum("sqm,srm->qrm", model.scores, model.scores)
    rotated = repeated_spatial_statistic(model)
    for field in ("statistic", "epsilon_gg", "p_uncorrected"):
        np.testing.assert_allclose(original[field], rotated[field], atol=1e-12)
    model.scores[:] = 2
    model.second_moment = np.einsum("sqm,srm->qrm", model.scores, model.scores)
    constant = repeated_spatial_statistic(model)
    assert np.all(np.isinf(constant["statistic"]))
    assert np.all(constant["p_uncorrected"] == 0)
    assert np.all(constant["effect_size"] == 1)


def test_incomplete_ids_and_zero_norm_are_rejected():
    grid, ids, coordinates = dataset()
    ids[1, 3][-1] = 99
    with pytest.raises(ValueError, match="same participant"):
        run(grid, ids, coordinates)
    ids[1, 3][-1] = 7
    grid[1, 3]["pats"][:] = 0
    with pytest.raises(ValueError, match="zero L2 norm"):
        run(grid, ids, coordinates)


def test_selected_effects_and_permutation_modes_preserve_statistics():
    grid, ids, coordinates = dataset(nsubs=6)
    all_effects = run(grid, ids, coordinates, n_permutations="auto",
                      cluster_forming_alpha=.2)
    selected = run(grid, ids, coordinates, effects=["Factor1"],
                   n_permutations="exact", cluster_forming_alpha=.2)
    assert tuple(selected.effects) == ("Factor1",)
    assert all_effects.exact_permutations and all_effects.n_permutations == 64
    np.testing.assert_array_equal(selected.effects["Factor1"].statistic,
                                  all_effects.effects["Factor1"].statistic)
    assert np.all(selected.effects["Factor1"].p_cluster_corrected <=
                  all_effects.effects["Factor1"].p_cluster_corrected)
    monte_carlo = run(grid, ids, coordinates, n_permutations=11,
                      cluster_forming_alpha=.2)
    assert not monte_carlo.exact_permutations and monte_carlo.n_permutations == 11
    for effect in monte_carlo.effects.values():
        for cluster in effect.clusters:
            assert cluster.p_corrected >= 1/12
            assert np.isclose(cluster.p_corrected*12, round(cluster.p_corrected*12))
