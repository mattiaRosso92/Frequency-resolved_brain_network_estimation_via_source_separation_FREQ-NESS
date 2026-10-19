from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from freqness import (
    FREQNESS_IndividualPeak,
    FREQNESSIndividualPeakResult,
    individual_peak,
)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _fixture():
    frex = np.array([4.0, 7.0, 8.0, 9.5, 11.0, 13.0, 16.0])
    values = np.array([
        [999, 1, 2, 9, 4, 2, 999],
        [999, 9, 1, 9, 1, 0, 999],
        [999, np.nan, np.inf, -np.inf, np.nan, np.nan, 999],
        [999, 0, 1, 2, 3, 10, 999],
        [999, -5, -2, -1, -2, -4, 999],
    ]).T
    evals = np.stack((values, values[::-1]), axis=0)
    return {"evals": evals, "frex": frex}


def test_individual_maxima_are_within_range_and_use_the_full_grid_indices():
    FREQ = _fixture()
    original = FREQ["evals"].copy()
    result = FREQNESS_IndividualPeak(FREQ, [7, 13])

    assert isinstance(result, FREQNESSIndividualPeakResult)
    np.testing.assert_allclose(result.peak_frex, [[9.5, 7, np.nan, 13, 9.5]])
    np.testing.assert_allclose(result.peak_value, [[9, 9, np.nan, 10, -1]])
    np.testing.assert_allclose(result.idx_peak, [[3, 1, np.nan, 5, 3]])
    np.testing.assert_array_equal(result.at_boundary, [[False, True, False, True, False]])
    np.testing.assert_array_equal(result.is_tied, [[False, True, False, False, False]])
    np.testing.assert_array_equal(result.frex, [7, 8, 9.5, 11, 13])
    np.testing.assert_array_equal(result.frex_range, [7, 13])
    np.testing.assert_array_equal(result.which_comp, [1])
    np.testing.assert_array_equal(FREQ["evals"], original)
    assert not result.figures
    assert not plt.get_fignums()
    assert individual_peak is FREQNESS_IndividualPeak


def test_non_grid_limits_select_only_interior_bins_without_snapping():
    result = FREQNESS_IndividualPeak(_fixture(), [7.1, 12.9])
    np.testing.assert_array_equal(result.frex, [8, 9.5, 11])
    np.testing.assert_allclose(result.peak_frex, [[9.5, 9.5, np.nan, 11, 9.5]])


def test_partial_missing_values_and_flat_landscape_have_explicit_results():
    FREQ = {"evals": np.array([[[np.inf, 2], [5, 2], [np.nan, 2]]]), "frex": [7, 9, 13]}
    result = FREQNESS_IndividualPeak(FREQ, [7, 13])
    np.testing.assert_allclose(result.peak_frex, [[9, 7]])
    np.testing.assert_array_equal(result.is_tied, [[False, True]])
    np.testing.assert_array_equal(result.at_boundary, [[False, True]])


def test_object_input_singletons_and_component_order_are_preserved():
    FREQ = _fixture()
    single = SimpleNamespace(evals=FREQ["evals"][:, :, 0], frex=FREQ["frex"][:, None])
    result = FREQNESS_IndividualPeak(single, [7, 13], which_comp=[2, 1, 2])
    assert result.peak_frex.shape == (2, 1)
    np.testing.assert_array_equal(result.which_comp, [2, 1])
    np.testing.assert_allclose(result.peak_frex[:, 0], [9.5, 9.5])
    one = FREQNESS_IndividualPeak(
        {"evals": FREQ["evals"][0:1, :, 0], "frex": FREQ["frex"]}, [7, 13]
    )
    assert one.peak_frex.shape == (1, 1)


@pytest.mark.parametrize("plot_avg,plot_all,nfigures", [(True, False, 2), (False, True, 4), (True, True, 6)])
def test_plot_flags_are_independent_and_do_not_change_results(plot_avg, plot_all, nfigures):
    FREQ = _fixture()
    FREQ["evals"] = FREQ["evals"][:, :, :2]
    result = FREQNESS_IndividualPeak(
        FREQ, [7, 13], which_comp=[1, 2],
        plot_avg=plot_avg, plot_all=plot_all, show=False,
    )
    plain = FREQNESS_IndividualPeak(FREQ, [7, 13], which_comp=[1, 2])
    assert len(result.figures) == nfigures
    np.testing.assert_array_equal(result.peak_frex, plain.peak_frex)
    assert any(axis.get_ylabel() == "Individual peak frequency (Hz)"
               for figure in result.figures for axis in figure.axes) == plot_avg


def test_missing_results_can_be_plotted_and_plotting_preserves_rng():
    FREQ = {"evals": np.full((1, 3, 1), np.nan), "frex": [7, 9, 13]}
    state_before = np.random.get_state()
    result = FREQNESS_IndividualPeak(FREQ, [7, 13], plot_avg=True, plot_all=True, show=False)
    state_after = np.random.get_state()
    assert len(result.figures) == 2
    np.testing.assert_array_equal(state_before[1], state_after[1])
    assert state_before[2:] == state_after[2:]


def test_group_box_uses_matlab_quartiles():
    frex = np.array([7, 8, 8.75, 9, 9.5, 9.75, 10, 10.25, 10.5, 11, 11.5, 12, 13])
    peaks = [8, 8.75, 9, 9.5, 9.75, 10, 10, 10.25, 10.5, 11, 11.5, 12]
    evals = np.zeros((1, frex.size, len(peaks)))
    for subi, peak in enumerate(peaks):
        evals[0, np.flatnonzero(frex == peak)[0], subi] = 1
    result = FREQNESS_IndividualPeak(
        {"evals": evals, "frex": frex}, [7, 13], plot_avg=True, show=False
    )
    box_y = result.figures[0].axes[0].patches[0].get_path().vertices[:, 1]
    np.testing.assert_allclose([np.min(box_y), np.max(box_y)], [9.25, 10.75])


@pytest.mark.parametrize("limits,message", [
    ([13, 7], "ascending"), ([7, 7], "ascending"), ([np.nan, 13], "ascending"),
    ([7], "two finite"), ([3, 13], "within"), ([7, 17], "within"),
    ([8.1, 9.6], "at least two"), (["7", "13"], "two finite"),
])
def test_invalid_ranges_are_rejected(limits, message):
    with pytest.raises(ValueError, match=message):
        FREQNESS_IndividualPeak(_fixture(), limits)


@pytest.mark.parametrize("which_comp", [0, 3, 1.5, [], True, np.nan, [1, np.inf]])
def test_invalid_components_are_rejected(which_comp):
    with pytest.raises(ValueError, match="one-based"):
        FREQNESS_IndividualPeak(_fixture(), [7, 13], which_comp=which_comp)


@pytest.mark.parametrize("field,value,message", [
    ("evals", np.zeros((2, 7, 1, 1)), "FREQ.evals"),
    ("evals", np.zeros((2, 7), dtype=complex), "FREQ.evals"),
    ("frex", [4, 7, 8], "frequency dimension"),
    ("frex", [4, 7, 8, 8, 11, 13, 16], "increasing"),
    ("frex", [0, 7, 8, 9.5, 11, 13, 16], "positive"),
])
def test_invalid_frequency_axes_and_shapes_are_rejected(field, value, message):
    FREQ = _fixture()
    FREQ[field] = value
    with pytest.raises(ValueError, match=message):
        FREQNESS_IndividualPeak(FREQ, [7, 13])


def test_invalid_flags_missing_fields_and_unknown_options_are_rejected():
    with pytest.raises(ValueError, match="FREQ.frex is required"):
        FREQNESS_IndividualPeak({"evals": [[1, 2]]}, [7, 13])
    with pytest.raises(TypeError, match="plot_avg"):
        FREQNESS_IndividualPeak(_fixture(), [7, 13], plot_avg=1)
    with pytest.raises(TypeError):
        FREQNESS_IndividualPeak(_fixture(), [7, 13], unknown=True)
