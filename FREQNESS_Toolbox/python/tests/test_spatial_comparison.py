from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from freqness import (
    FREQNESS_SpatialComparison,
    FREQNESSSpatialComparisonResult,
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
    assert result.design == "paired"
    assert result.n_participants == 6
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
    with pytest.raises(ValueError, match="matching dimensions"):
        FREQNESS_SpatialComparison(
            condition_a,
            SimpleNamespace(pats=condition_b.pats[..., :-1], frex=condition_b.frex),
            coordinates,
            frex2compare=[8],
            plot_results=False,
        )


def test_public_alias_references_matlab_named_function():
    assert spatial_comparison is FREQNESS_SpatialComparison
