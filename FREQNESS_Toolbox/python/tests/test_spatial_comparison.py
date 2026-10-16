from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from freqness import (
    FREQNESS_SpatialComparison,
    FREQNESSSpatialComparisonResult,
    FREQNESSSpatialFactorialComparisonResult,
    spatial_comparison,
)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def paired_patterns():
    coordinates = np.asarray(list(np.ndindex(3, 3, 3)), dtype=float) * 8.0
    frequencies = np.array([4.0, 8.0, 12.0])
    rng = np.random.default_rng(12)
    condition_a = rng.normal(0.2, 0.03, (27, 2, 3, 6))
    condition_b = condition_a + rng.normal(0.0, 0.015, condition_a.shape)
    central_cluster = np.array([4, 5, 13, 14])
    condition_b[central_cluster, 0, 1, :] += 0.25
    return (
        SimpleNamespace(pats=condition_a, frex=frequencies),
        SimpleNamespace(pats=condition_b, frex=frequencies),
        coordinates,
        central_cluster,
    )


def test_spatial_comparison_returns_structured_cluster_corrected_results(
    paired_patterns,
):
    condition_a, condition_b, coordinates, central_cluster = paired_patterns

    result = FREQNESS_SpatialComparison(
        condition_a,
        condition_b,
        coordinates,
        frex2compare=[8],
        normalization="none",
        cluster_forming_alpha=0.05,
        n_permutations="exact",
        plot_results=False,
    )

    assert isinstance(result, FREQNESSSpatialComparisonResult)
    assert result.design == "within"
    assert result.test == "ttest"
    assert result.n_participants == 6
    assert result.n_A == result.n_B == 6
    assert result.exact_permutations
    assert result.n_permutations == 64
    assert result.mean_difference.shape == (27, 1, 1)
    assert result.comps.tolist() == [1]
    assert result.idx_frex.tolist() == [1]
    assert np.any(result.significant_mask[central_cluster, 0, 0])
    assert any(cluster.p_corrected <= 0.05 for cluster in result.clusters)


def test_nearest_frequency_resolution_is_reported(paired_patterns):
    condition_a, condition_b, coordinates, _ = paired_patterns

    with pytest.warns(UserWarning, match="closest matches"):
        result = FREQNESS_SpatialComparison(
            condition_a,
            condition_b,
            coordinates,
            frex2compare=[7.8, 12],
            comps2compare=[2],
            n_permutations=16,
            plot_results=False,
        )

    np.testing.assert_allclose(result.requested_frex, [7.8, 12])
    np.testing.assert_allclose(result.frex, [8, 12])
    np.testing.assert_array_equal(result.idx_frex, [1, 2])
    np.testing.assert_array_equal(result.comps, [2])
    assert result.correction_scope == "all_maps"


def test_all_frequency_selection_and_figures(paired_patterns):
    condition_a, condition_b, coordinates, _ = paired_patterns

    result = FREQNESS_SpatialComparison(
        condition_a,
        condition_b,
        coordinates,
        frex2compare="all",
        comps2compare=[1, 2],
        n_permutations=8,
        show=False,
    )

    np.testing.assert_allclose(result.frex, condition_a.frex)
    assert result.mean_difference.shape == (27, 3, 2)
    assert len(plt.get_fignums()) == 6


def test_magnitude_l2_comparison_is_sign_and_scale_invariant(paired_patterns):
    condition_a, condition_b, coordinates, _ = paired_patterns
    baseline = FREQNESS_SpatialComparison(
        condition_a,
        condition_b,
        coordinates,
        frex2compare=[8],
        n_permutations=8,
        plot_results=False,
    )
    transformed_a = SimpleNamespace(
        pats=-3.0 * condition_a.pats,
        frex=condition_a.frex,
    )
    transformed_b = SimpleNamespace(
        pats=7.0 * condition_b.pats,
        frex=condition_b.frex,
    )
    transformed = FREQNESS_SpatialComparison(
        transformed_a,
        transformed_b,
        coordinates,
        frex2compare=[8],
        n_permutations=8,
        plot_results=False,
    )

    np.testing.assert_allclose(transformed.mean_difference, baseline.mean_difference)
    np.testing.assert_allclose(transformed.t_values, baseline.t_values)
    np.testing.assert_allclose(
        transformed.p_cluster_corrected, baseline.p_cluster_corrected
    )


def test_signed_alignment_removes_arbitrary_condition_sign(paired_patterns):
    condition_a, _, coordinates, _ = paired_patterns
    condition_b = SimpleNamespace(pats=-condition_a.pats, frex=condition_a.frex)

    result = FREQNESS_SpatialComparison(
        condition_a,
        condition_b,
        coordinates,
        frex2compare=[4],
        map_mode="signed_aligned",
        n_permutations=8,
        plot_results=False,
    )

    np.testing.assert_allclose(result.mean_difference, 0.0, atol=1e-14)
    assert np.all(result.sign_flips_b == -1)
    assert not np.any(result.significant_mask)


def test_frequency_and_input_safeguards(paired_patterns):
    condition_a, condition_b, coordinates, _ = paired_patterns
    with pytest.raises(TypeError):
        FREQNESS_SpatialComparison(condition_a, condition_b, coordinates)
    with pytest.raises(ValueError, match="outside"):
        FREQNESS_SpatialComparison(
            condition_a,
            condition_b,
            coordinates,
            frex2compare=[30],
            plot_results=False,
        )
    with pytest.raises(ValueError, match="same analyzed bin"):
        FREQNESS_SpatialComparison(
            condition_a,
            condition_b,
            coordinates,
            frex2compare=[7.9, 8.1],
            plot_results=False,
        )
    with pytest.raises(ValueError, match="matching participant counts"):
        FREQNESS_SpatialComparison(
            condition_a,
            SimpleNamespace(pats=condition_b.pats[..., :-1], frex=condition_b.frex),
            coordinates,
            frex2compare=[8],
            plot_results=False,
        )


def test_between_subject_ttest_supports_unequal_group_sizes():
    coordinates = np.asarray(list(np.ndindex(3, 3, 3)), dtype=float) * 8.0
    rng = np.random.default_rng(44)
    group_a = rng.normal(0.2, 0.03, (27, 1, 1, 4))
    group_b = rng.normal(0.2, 0.03, (27, 1, 1, 5))
    group_b[[4, 5, 13, 14], 0, 0, :] += 0.3

    result = FREQNESS_SpatialComparison(
        SimpleNamespace(pats=group_a, frex=[8.0]),
        SimpleNamespace(pats=group_b, frex=[8.0]),
        coordinates,
        frex2compare=[8],
        design="between",
        test="ttest",
        normalization="none",
        cluster_forming_alpha=0.05,
        n_permutations="exact",
        plot_results=False,
    )

    assert result.design == "between"
    assert result.test == "ttest"
    assert result.n_participants == (4, 5)
    assert result.n_A == 4
    assert result.n_B == 5
    assert result.statistic_name == "t"
    assert result.effect_size_name == "Hedges' g"
    assert result.exact_permutations
    assert result.n_permutations == 126
    assert result.sign_flips_a.shape[-1] == 4
    assert result.sign_flips_b.shape[-1] == 5


def test_incompatible_nonparametric_tests_warn_and_are_corrected(paired_patterns):
    condition_a, condition_b, coordinates, _ = paired_patterns
    with pytest.warns(UserWarning, match="Using test='signrank'"):
        within = FREQNESS_SpatialComparison(
            condition_a,
            condition_b,
            coordinates,
            frex2compare=[8],
            design="within",
            test="ranksum",
            cluster_forming_alpha=0.05,
            n_permutations=8,
            plot_results=False,
        )
    assert within.test == "signrank"
    assert within.statistic_name == "z"
    assert within.effect_size_name == "matched rank-biserial r"
    assert np.all(np.isnan(within.t_values))

    between_b = SimpleNamespace(
        pats=condition_b.pats[..., :5],
        frex=condition_b.frex,
    )
    with pytest.warns(UserWarning, match="Using test='ranksum'"):
        between = FREQNESS_SpatialComparison(
            condition_a,
            between_b,
            coordinates,
            frex2compare=[8],
            design="between",
            test="signrank",
            cluster_forming_alpha=0.05,
            n_permutations=8,
            plot_results=False,
        )
    assert between.test == "ranksum"
    assert between.statistic_name == "z"
    assert between.effect_size_name == "rank-biserial r"


def test_public_alias_references_matlab_named_function():
    assert spatial_comparison is FREQNESS_SpatialComparison


def test_factorial_within_subject_anova_uses_cluster_permutations():
    coordinates = np.asarray(list(np.ndindex(3, 3, 3)), dtype=float) * 8.0
    rng = np.random.default_rng(91)
    frex = np.array([6.0, 10.0])
    grid = []
    identifiers = []
    baseline = rng.normal(0.0, 0.08, (27, 1, 2, 7))
    for level in range(3):
        patterns = baseline + rng.normal(0.0, 0.03, baseline.shape)
        patterns[[4, 5, 13, 14], 0, 1] += level * 0.35
        grid.append(SimpleNamespace(pats=patterns, frex=frex))
        identifiers.append([f"P{index}" for index in range(7)])

    result = FREQNESS_SpatialComparison(
        grid,
        coordinates,
        frex2compare="all",
        design="within",
        factor_names=["Condition"],
        level_labels=[["Low", "Medium", "High"]],
        subject_ids=identifiers,
        normalization="none",
        cluster_forming_alpha=0.05,
        n_permutations=19,
        show=False,
    )

    assert isinstance(result, FREQNESSSpatialFactorialComparisonResult)
    assert result.analysis == "factorial"
    assert result.test == "anova"
    assert result.factor_names == ("Condition",)
    assert result.effect_names == ("Condition",)
    assert result.correction_scope == "all_effects_and_maps"
    assert result.n_permutations == 19
    assert result.effects["Condition"].statistic.shape == (27, 2, 1)
    assert result.effects["Condition"].p_cluster_corrected.shape == (27, 2, 1)
    assert len(result.figures) == 1


def test_factorial_mixed_design_exposes_main_effects_and_interaction():
    coordinates = np.asarray(list(np.ndindex(2, 2, 2)), dtype=float) * 8.0
    rng = np.random.default_rng(92)
    grid = [[None, None], [None, None]]
    identifiers = [[None, None], [None, None]]
    for group in range(2):
        nsubjects = 5 + group
        group_ids = [f"G{group}P{index}" for index in range(nsubjects)]
        baseline = rng.normal(0.0, 0.1, (8, 1, 1, nsubjects))
        for condition in range(2):
            patterns = baseline + rng.normal(0.0, 0.03, baseline.shape)
            patterns[:4] += 0.2 * group + 0.3 * condition
            grid[group][condition] = SimpleNamespace(pats=patterns, frex=[8.0])
            identifiers[group][condition] = group_ids

    result = FREQNESS_SpatialComparison(
        grid,
        coordinates,
        frex2compare=[8],
        design=["between", "within"],
        factor_names=["Group", "Condition"],
        subject_ids=identifiers,
        normalization="none",
        cluster_forming_alpha=0.1,
        n_permutations=9,
        plot_results=False,
    )

    assert result.factor_design == ("between", "within")
    assert result.effect_names == ("Group", "Condition", "Group_x_Condition")
    np.testing.assert_array_equal(result.cell_n_participants, [[5, 5], [6, 6]])
    for effect in result.effects.values():
        assert effect.statistic.shape == (8, 1, 1)


def test_factorial_two_level_grid_preserves_pairwise_statistics(paired_patterns):
    condition_a, condition_b, coordinates, _ = paired_patterns
    direct = FREQNESS_SpatialComparison(
        condition_a,
        condition_b,
        coordinates,
        frex2compare=[8],
        normalization="none",
        n_permutations=16,
        plot_results=False,
    )
    order = np.array([2, 0, 5, 1, 4, 3])
    reordered_b = SimpleNamespace(
        pats=condition_b.pats[..., order], frex=condition_b.frex
    )
    generalized = FREQNESS_SpatialComparison(
        [condition_a, reordered_b],
        coordinates,
        frex2compare=[8],
        normalization="none",
        subject_ids=[
            [f"P{index}" for index in range(6)],
            [f"P{index}" for index in order],
        ],
        n_permutations=16,
        plot_results=False,
    )

    effect = generalized.effects["Factor1"]
    np.testing.assert_allclose(effect.statistic, direct.statistic)
    np.testing.assert_allclose(
        effect.p_cluster_corrected, direct.p_cluster_corrected
    )
