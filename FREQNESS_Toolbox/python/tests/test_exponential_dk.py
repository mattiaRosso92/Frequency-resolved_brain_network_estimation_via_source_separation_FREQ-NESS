import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from freqness import FREQNESS_ExponentialDK, FREQNESSExponentialGoodFit


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _exact_exponential_input():
    frequencies = np.arange(1.0, 7.0)
    decay = np.array([0.2, 0.45])
    amplitude = np.array([12.0, 8.0])
    component_one = np.column_stack(
        [
            amplitude[subject] * np.exp(-decay[subject] * frequencies)
            for subject in range(2)
        ]
    )
    component_two = np.column_stack(
        [
            5.0 * np.exp(-0.1 * frequencies),
            6.0 * np.exp(-0.15 * frequencies),
        ]
    )
    eigenvalues = np.stack([component_one, component_two], axis=0)
    return {"evals": eigenvalues, "frex": frequencies}, decay


def test_exact_exponential_decay_and_original_space_r_squared():
    FREQ, expected_decay = _exact_exponential_input()

    decay_coefficients, good_fit = FREQNESS_ExponentialDK(FREQ, show=False)

    assert isinstance(good_fit, FREQNESSExponentialGoodFit)
    assert decay_coefficients.shape == (2, 1)
    assert good_fit.R2.shape == (2, 1)
    np.testing.assert_allclose(decay_coefficients[:, 0], expected_decay, atol=1e-12)
    np.testing.assert_allclose(good_fit.R2, 1.0, atol=1e-12)


def test_one_based_component_selection():
    FREQ, _ = _exact_exponential_input()

    decay_coefficients, _ = FREQNESS_ExponentialDK(
        FREQ,
        which_comp=2,
        show=False,
    )

    np.testing.assert_allclose(decay_coefficients[:, 0], [0.1, 0.15], atol=1e-12)


def test_invalid_component_warns_and_defaults_to_first():
    FREQ, expected_decay = _exact_exponential_input()

    with pytest.warns(UserWarning, match="Defaulting to component 1"):
        decay_coefficients, _ = FREQNESS_ExponentialDK(
            FREQ,
            which_comp=99,
            show=False,
        )

    np.testing.assert_allclose(decay_coefficients[:, 0], expected_decay, atol=1e-12)


def test_range_uses_closest_frequencies_and_marks_actual_boundaries():
    FREQ, expected_decay = _exact_exponential_input()

    with pytest.warns(UserWarning, match="closest matches"):
        decay_coefficients, _ = FREQNESS_ExponentialDK(
            FREQ,
            range2fit=[1.2, 5.8],
            show=False,
        )

    np.testing.assert_allclose(decay_coefficients[:, 0], expected_decay, atol=1e-12)
    axis = plt.gcf().axes[0]
    vertical_lines = [line for line in axis.lines if line.get_linestyle() == "--"]
    assert len(vertical_lines) == 2
    np.testing.assert_array_equal(vertical_lines[0].get_xdata(), [1, 1])
    np.testing.assert_array_equal(vertical_lines[1].get_xdata(), [6, 6])


def test_missing_frequency_axis_uses_indices_but_rejects_hz_range():
    x = np.arange(1.0, 6.0)
    values = 9.0 * np.exp(-0.3 * x)

    decay_coefficients, _ = FREQNESS_ExponentialDK(
        {"evals": np.vstack([values, values * 0.5])},
        show=False,
    )
    np.testing.assert_allclose(decay_coefficients, [[0.3]], atol=1e-12)
    assert plt.gcf().axes[0].get_xlabel() == "Frequency index"

    with pytest.raises(ValueError, match="FREQ.frex is missing"):
        FREQNESS_ExponentialDK(
            {"evals": np.vstack([values, values * 0.5])},
            range2fit=[1, 4],
            show=False,
        )


def test_empty_range_sequence_fits_all_frequencies():
    FREQ, expected_decay = _exact_exponential_input()

    decay_coefficients, _ = FREQNESS_ExponentialDK(
        FREQ,
        range2fit=[],
        show=False,
    )

    np.testing.assert_allclose(decay_coefficients[:, 0], expected_decay, atol=1e-12)


def test_insufficient_positive_points_returns_nan_without_fit_curve():
    FREQ = {
        "evals": np.array([[4.0, 0.0, np.nan, 2.0], [1.0, 1.0, 1.0, 1.0]]),
        "frex": [1, 2, 3, 4],
    }

    decay_coefficients, good_fit = FREQNESS_ExponentialDK(FREQ, show=False)

    assert np.isnan(decay_coefficients[0, 0])
    assert np.isnan(good_fit.R2[0, 0])
    labels = [line.get_label() for line in plt.gcf().axes[0].lines]
    assert "Exponential fit" not in labels


def test_plot_all_creates_group_and_participant_figures():
    FREQ, _ = _exact_exponential_input()

    FREQNESS_ExponentialDK(FREQ, plot_all=True, show=False)

    figures = [plt.figure(number) for number in plt.get_fignums()]
    assert len(figures) == 3
    assert figures[0].axes[0].get_title() == "Exponential decay fit - Component #1"
    assert "Participant #1" in figures[1].axes[0].get_title()
    assert "Participant #2" in figures[2].axes[0].get_title()


@pytest.mark.parametrize(
    ("range2fit", "message"),
    [
        ([1], "two-element"),
        ([4, 2], "ascending"),
        ([1, np.inf], "finite"),
    ],
)
def test_invalid_fit_ranges_are_rejected(range2fit, message):
    FREQ, _ = _exact_exponential_input()
    with pytest.raises(ValueError, match=message):
        FREQNESS_ExponentialDK(FREQ, range2fit=range2fit, show=False)
