from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from freqness import FREQNESS_EntropyLandscape


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_entropy_and_effective_dimensionality_match_definitions():
    eigenvalues = np.array(
        [
            [[4.0, 3.0], [2.0, 5.0]],
            [[2.0, 3.0], [2.0, 3.0]],
            [[1.0, 2.0], [2.0, 2.0]],
        ]
    )
    FREQ = SimpleNamespace(evals=eigenvalues, frex=np.array([[4.0], [8.0]]))

    ENTR = FREQNESS_EntropyLandscape(FREQ, show=False)

    probabilities = eigenvalues / np.sum(eigenvalues, axis=0, keepdims=True)
    expected_h2 = -np.log(np.sum(probabilities**2, axis=0))
    expected_ed = np.sum(eigenvalues, axis=0) ** 2 / np.sum(
        eigenvalues**2, axis=0
    )
    np.testing.assert_allclose(ENTR.H2, expected_h2)
    np.testing.assert_allclose(ENTR.ED, expected_ed)
    np.testing.assert_allclose(np.exp(ENTR.H2), ENTR.ED)
    np.testing.assert_array_equal(ENTR.frex, [4.0, 8.0])
    assert ENTR.H2.shape == (2, 2)
    assert len(ENTR.figures) == 1


def test_vector_input_retains_frequency_and_participant_dimensions():
    ENTR = FREQNESS_EntropyLandscape(
        {"evals": np.array([4.0, 2.0, 2.0])},
        show=False,
    )

    expected_ed = 64.0 / 24.0
    assert ENTR.H2.shape == (1, 1)
    assert ENTR.ED.shape == (1, 1)
    assert ENTR.H2[0, 0] == pytest.approx(np.log(expected_ed))
    assert ENTR.ED[0, 0] == pytest.approx(expected_ed)
    np.testing.assert_array_equal(ENTR.frex, [1.0])


def test_group_figure_contains_h2_but_not_effective_dimensionality():
    eigenvalues = np.array(
        [
            [[4.0, 3.0], [2.0, 5.0]],
            [[2.0, 3.0], [2.0, 3.0]],
            [[1.0, 2.0], [2.0, 2.0]],
        ]
    )

    FREQNESS_EntropyLandscape(
        {"evals": eigenvalues, "frex": [4, 8]},
        show=False,
    )

    figures = [plt.figure(number) for number in plt.get_fignums()]
    assert len(figures) == 1
    assert len(figures[0].axes) == 1
    axis = figures[0].axes[0]
    assert "Rényi Entropy" in axis.get_title()
    assert axis.get_ylabel() == r"$H_2$"
    all_text = " ".join(text.get_text() for text in figures[0].findobj(plt.Text))
    assert "Effective Dimensionality" not in all_text
    np.testing.assert_array_equal(axis.lines[0].get_xdata(), [4, 8])


def test_plot_all_adds_one_h2_figure_per_participant():
    eigenvalues = np.ones((3, 2, 3), dtype=float)

    FREQNESS_EntropyLandscape(
        {"evals": eigenvalues, "frex": [2, 6]},
        plot_all=True,
        show=False,
    )

    figures = [plt.figure(number) for number in plt.get_fignums()]
    assert len(figures) == 4
    titles = [figure.axes[0].get_title() for figure in figures]
    assert titles[0].endswith("Grand Average")
    assert titles[1:] == [
        "Quadratic Rényi Entropy - Participant #1",
        "Quadratic Rényi Entropy - Participant #2",
        "Quadratic Rényi Entropy - Participant #3",
    ]


def test_frequency_index_fallback_and_output_shape():
    ENTR = FREQNESS_EntropyLandscape(
        {"evals": np.ones((4, 3))},
        show=False,
    )

    assert ENTR.H2.shape == (3, 1)
    assert ENTR.ED.shape == (3, 1)
    np.testing.assert_array_equal(ENTR.frex, [1, 2, 3])
    axis = plt.gcf().axes[0]
    np.testing.assert_array_equal(axis.lines[0].get_xdata(), [1, 2, 3])
    assert axis.get_xlabel() == "Frequency index"


@pytest.mark.parametrize(
    ("eigenvalues", "message"),
    [
        (np.ones((2, 2, 2, 2)), "1D, 2D, or 3D"),
        (np.array([[1.0, -1.0], [2.0, 3.0]]), "negative"),
    ],
)
def test_entropy_landscape_rejects_invalid_eigenspectra(eigenvalues, message):
    with pytest.raises(ValueError, match=message):
        FREQNESS_EntropyLandscape({"evals": eigenvalues}, show=False)
