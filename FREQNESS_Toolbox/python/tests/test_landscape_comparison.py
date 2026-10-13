from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from freqness import (
    FREQNESS_LandscapeComparison,
    FREQNESSLandscapeComparisonResult,
    landscape_comparison,
)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _paired_freq(seed: int = 31):
    rng = np.random.default_rng(seed)
    frex = np.array([4.0, 8.0, 12.0, 16.0])
    values_a = rng.normal(0.3, 0.04, (2, frex.size, 14))
    values_b = values_a + rng.normal(0.0, 0.01, values_a.shape)
    values_b[0, 1:3] += 0.15
    return (
        SimpleNamespace(evals=values_a, frex=frex),
        SimpleNamespace(evals=values_b, frex=frex),
    )


def test_within_subject_ttest_returns_componentwise_fdr_results():
    condition_a, condition_b = _paired_freq()
    result = FREQNESS_LandscapeComparison(
        condition_a,
        condition_b,
        measure="evals",
        test="ttest",
        which_comp=[1, 2],
        plot_results=False,
    )

    assert isinstance(result, FREQNESSLandscapeComparisonResult)
    assert result.design == "within"
    assert result.test == "ttest"
    assert result.correction == "fdr"
    assert result.measures == ("evals",)
    assert result.evals is result.results["evals"]
    assert result.evals.statistic.shape == (2, 4)
    np.testing.assert_array_equal(result.evals.comps, [1, 2])
    assert np.all(result.evals.significant[0, 1:3])
    assert not np.any(result.evals.significant[1])


def test_omitted_measure_preserves_each_entropy_field():
    rng = np.random.default_rng(8)
    frex = np.array([2.0, 4.0, 6.0])
    entropy_a = SimpleNamespace(
        H2=rng.normal(1.0, 0.03, (3, 10)),
        ED=rng.normal(3.0, 0.08, (3, 10)),
        frex=frex,
    )
    entropy_b = SimpleNamespace(
        H2=entropy_a.H2 + 0.2,
        ED=entropy_a.ED - 0.3,
        frex=frex,
    )

    result = FREQNESS_LandscapeComparison(
        entropy_a,
        entropy_b,
        plot_results=False,
    )

    assert result.measures == ("H2", "ED")
    assert result.H2.measure == "H2"
    assert result.ED.measure == "ED"
    np.testing.assert_array_equal(result.H2.comps, [1])
    np.testing.assert_allclose(result.H2.frex, frex)


def test_between_subject_design_allows_different_group_sizes():
    rng = np.random.default_rng(19)
    frex = np.array([5.0, 10.0, 15.0])
    group_a = SimpleNamespace(
        stability_index=rng.normal(0.4, 0.04, (2, 3, 9)),
        frex=frex,
    )
    group_b = SimpleNamespace(
        stability_index=rng.normal(0.65, 0.04, (2, 3, 13)),
        frex=frex,
    )

    nonparametric = FREQNESS_LandscapeComparison(
        group_a,
        group_b,
        measure="stability_index",
        design="between",
        which_comp=2,
        plot_results=False,
    )
    parametric = FREQNESS_LandscapeComparison(
        group_a,
        group_b,
        measure="stability_index",
        design="between",
        test="ttest",
        which_comp=2,
        plot_results=False,
    )

    assert nonparametric.test == "ranksum"
    assert nonparametric.stability_index.statistic_name == "U"
    assert parametric.stability_index.effect_size_name == "Hedges' g"
    np.testing.assert_array_equal(nonparametric.stability_index.n_A, [[9, 9, 9]])
    np.testing.assert_array_equal(nonparametric.stability_index.n_B, [[13, 13, 13]])


def test_cross_coupling_measures_use_carrier_frequency_axis():
    rng = np.random.default_rng(22)
    carrier_frex = np.array([12.0, 16.0, 20.0])
    coupling_a = SimpleNamespace(
        amplitude_raw=rng.normal(0.1, 0.01, (3, 8)),
        amplitude_normalized=rng.normal(10.0, 1.0, (3, 8)),
        carrier_frex=carrier_frex,
    )
    coupling_b = SimpleNamespace(
        amplitude_raw=coupling_a.amplitude_raw + 0.03,
        amplitude_normalized=coupling_a.amplitude_normalized + 2.0,
        carrier_frex=carrier_frex,
    )

    result = FREQNESS_LandscapeComparison(
        coupling_a,
        coupling_b,
        plot_results=False,
    )

    assert result.measures == ("amplitude_raw", "amplitude_normalized")
    np.testing.assert_allclose(result.amplitude_raw.frex, carrier_frex)
    np.testing.assert_allclose(result.amplitude_normalized.frex, carrier_frex)


def test_raw_numeric_landscapes_require_frequencies_and_accept_lists():
    values_a = [[1.0, 1.1, 0.9], [2.0, 2.1, 1.9]]
    values_b = [[1.3, 1.4, 1.2], [2.3, 2.4, 2.2]]

    with pytest.raises(ValueError, match="frex is required"):
        FREQNESS_LandscapeComparison(
            values_a,
            values_b,
            plot_results=False,
        )

    result = FREQNESS_LandscapeComparison(
        values_a,
        values_b,
        frex=[6, 12],
        plot_results=False,
    )
    assert result.measures == ("values",)
    np.testing.assert_allclose(result.values.frex, [6, 12])


def test_design_test_recovery_and_plot_count():
    condition_a, condition_b = _paired_freq()
    with pytest.warns(UserWarning, match="Using test='ranksum'"):
        between = FREQNESS_LandscapeComparison(
            condition_a,
            condition_b,
            design="between",
            test="signrank",
            plot_results=False,
        )
    assert between.test == "ranksum"

    with pytest.warns(UserWarning, match="Using test='signrank'"):
        within = FREQNESS_LandscapeComparison(
            condition_a,
            condition_b,
            design="within",
            test="ranksum",
            plot_results=False,
        )
    assert within.test == "signrank"

    result = FREQNESS_LandscapeComparison(
        condition_a,
        condition_b,
        measure="evals",
        which_comp=[1, 2],
        show=False,
    )
    assert len(result.figures) == 2
    assert len(plt.get_fignums()) == 2
    assert landscape_comparison is FREQNESS_LandscapeComparison
