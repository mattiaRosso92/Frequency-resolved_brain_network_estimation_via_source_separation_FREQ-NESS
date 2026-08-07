from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from freqness import (
    FREQNESS_CompGradients,
    FREQNESS_FreqGradients,
    FREQNESSGradientGoodFit,
    component_gradients,
    frequency_gradients,
)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def gradient_data():
    x = np.linspace(-45.0, 45.0, 181)
    MNI = np.column_stack([x, 0.5 * x + 2.0, -0.25 * x + 5.0])
    ncomponents = 4
    frequencies = np.array([2.0, 4.0, 8.0, 12.0, 20.0])
    patterns = np.empty((x.size, ncomponents, frequencies.size, 2))
    for subject in range(2):
        for component in range(ncomponents):
            for frequency in range(frequencies.size):
                center = -25.0 + 11.0 * frequency + 8.0 * component + 0.5 * subject
                patterns[:, component, frequency, subject] = (
                    np.exp(-0.5 * ((x - center) / 2.5) ** 2) + 0.01
                )
    eigenvalues = np.ones((ncomponents, frequencies.size, 2))
    eigenvalues[0, :, :] = np.array([1.0, 8.0, 3.0, 2.0, 1.0])[:, None]
    FREQ = SimpleNamespace(
        pats=patterns,
        frex=frequencies,
        evals=eigenvalues,
    )
    return FREQ, MNI


def test_frequency_gradients_recover_linear_x_gradient(gradient_data):
    FREQ, MNI = gradient_data

    coefficients, good_fit = FREQNESS_FreqGradients(
        FREQ,
        MNI,
        comp2model=1,
        show=False,
    )

    assert coefficients.shape == (3, 3, 2)
    assert isinstance(good_fit, FREQNESSGradientGoodFit)
    assert good_fit.R2_linear.shape == (3, 2)
    assert good_fit.threshold_sd == 1.0
    assert good_fit.retained_voxels.shape == (5, 2)
    np.testing.assert_array_equal(good_fit.bestOrder[0], [1.0, 1.0])
    np.testing.assert_allclose(coefficients[0, 1], 1.0 / 11.0, atol=0.005)
    np.testing.assert_allclose(coefficients[0, 2], 0.0, atol=1e-12)
    assert np.all(good_fit.R2_best[0] > 0.98)
    assert len(good_fit.figures) == 1


def test_component_gradients_recover_linear_x_gradient(gradient_data):
    FREQ, MNI = gradient_data

    coefficients, good_fit = FREQNESS_CompGradients(
        FREQ,
        MNI,
        freq2model=8.0,
        show=False,
    )

    assert coefficients.shape == (3, 3, 2)
    assert good_fit.threshold_sd == 2.0
    assert good_fit.retained_voxels.shape == (4, 2)
    np.testing.assert_array_equal(good_fit.bestOrder[0], [1.0, 1.0])
    np.testing.assert_allclose(coefficients[0, 1], 1.0 / 8.0, atol=0.01)
    np.testing.assert_allclose(coefficients[0, 2], 0.0, atol=1e-12)
    assert np.all(good_fit.R2_best[0] > 0.96)


def test_coefficients_are_b0_b1_b2_and_quadratic_can_win():
    x = np.array([0.0, 1.0, np.sqrt(2.0), np.sqrt(3.0), 2.0])
    MNI = np.column_stack([x, np.zeros_like(x), np.zeros_like(x)])
    patterns = np.eye(5)[:, None, :, None]
    FREQ = {"pats": patterns, "frex": np.arange(1.0, 6.0)}

    coefficients, good_fit = FREQNESS_FreqGradients(
        FREQ,
        MNI,
        threshold_sd=0,
        show=False,
    )

    assert good_fit.bestOrder[0, 0] == 2
    np.testing.assert_allclose(coefficients[0, :, 0], [1.0, 0.0, 1.0], atol=1e-12)
    assert good_fit.R2_best[0, 0] == pytest.approx(1.0)
    assert np.isnan(good_fit.bestOrder[1:, 0]).all()


def test_frequency_range_is_sorted_and_nearest_values_warn(gradient_data):
    FREQ, MNI = gradient_data

    with pytest.warns(UserWarning, match="closest matches"):
        _, good_fit = FREQNESS_FreqGradients(
            FREQ,
            MNI,
            frex2model=[12.2, 3.9],
            show=False,
        )

    axis = good_fit.figures[0].axes[0]
    dashed_levels = sorted(
        float(line.get_ydata()[0])
        for line in axis.lines
        if line.get_linestyle() == "--"
    )
    assert dashed_levels == [2.0, 4.0]


def test_component_range_is_sorted_and_inclusive(gradient_data):
    FREQ, MNI = gradient_data

    _, good_fit = FREQNESS_CompGradients(
        FREQ,
        MNI,
        freq2model=8,
        comps2model=[4, 2],
        show=False,
    )

    axis = good_fit.figures[0].axes[0]
    dashed_levels = sorted(
        float(line.get_ydata()[0])
        for line in axis.lines
        if line.get_linestyle() == "--"
    )
    assert dashed_levels == [2.0, 4.0]


def test_component_gradient_automatically_selects_peak_frequency(gradient_data):
    FREQ, MNI = gradient_data

    with pytest.warns(UserWarning, match="most prominent"):
        _, good_fit = FREQNESS_CompGradients(FREQ, MNI, show=False)

    assert "4 Hz" in good_fit.figures[0]._suptitle.get_text()


def test_missing_frequency_axis_uses_indices_for_component_gradient(gradient_data):
    FREQ, MNI = gradient_data
    no_frequencies = {"pats": FREQ.pats, "evals": FREQ.evals}

    with pytest.warns(UserWarning, match="One-based frequency indices"):
        _, good_fit = FREQNESS_CompGradients(
            no_frequencies,
            MNI,
            freq2model=3,
            show=False,
        )

    assert "frequency index 3" in good_fit.figures[0]._suptitle.get_text()


def test_missing_frequency_axis_rejects_hz_range(gradient_data):
    FREQ, MNI = gradient_data

    with pytest.warns(UserWarning, match="One-based frequency indices"):
        with pytest.raises(ValueError, match="FREQ.frex is missing"):
            FREQNESS_FreqGradients(
                {"pats": FREQ.pats},
                MNI,
                frex2model=[4, 12],
                show=False,
            )


def test_all_zero_maps_return_nan_fits_and_zero_retained_counts():
    patterns = np.zeros((20, 3, 4, 1))
    MNI = np.column_stack(
        [np.arange(20.0), np.arange(20.0), np.arange(20.0)]
    )

    coefficients, good_fit = FREQNESS_FreqGradients(
        {"pats": patterns, "frex": [2, 4, 6, 8]},
        MNI,
        show=False,
    )

    assert np.isnan(coefficients).all()
    assert np.isnan(good_fit.R2_best).all()
    assert not np.any(good_fit.retained_voxels)


def test_plot_all_adds_one_figure_per_participant(gradient_data):
    FREQ, MNI = gradient_data

    _, good_fit = FREQNESS_FreqGradients(
        FREQ,
        MNI,
        plot_all=True,
        show=False,
    )

    assert len(good_fit.figures) == 3
    assert all(len(figure.axes) == 3 for figure in good_fit.figures)
    assert good_fit.figures[1]._suptitle.get_text().endswith("Participant #1")


def test_public_aliases_reference_matlab_named_functions():
    assert frequency_gradients is FREQNESS_FreqGradients
    assert component_gradients is FREQNESS_CompGradients


@pytest.mark.parametrize("function", [FREQNESS_FreqGradients, FREQNESS_CompGradients])
def test_gradient_functions_validate_mni_and_threshold(function, gradient_data):
    FREQ, MNI = gradient_data
    with pytest.raises(ValueError, match="shape"):
        function(FREQ, MNI[:-1], show=False)
    with pytest.raises(ValueError, match="threshold_sd"):
        function(FREQ, MNI, threshold_sd=-1, show=False)


def test_component_gradient_requires_eigenvalues_for_automatic_selection(gradient_data):
    FREQ, MNI = gradient_data
    with pytest.warns(UserWarning, match="most prominent"):
        with pytest.raises(ValueError, match="FREQ.evals"):
            FREQNESS_CompGradients(
                {"pats": FREQ.pats, "frex": FREQ.frex},
                MNI,
                show=False,
            )
