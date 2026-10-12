import importlib

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from freqness import (
    FREQNESS_StabilityIndex,
    FREQNESSStabilityIndexResult,
    stability_index,
)
from freqness.FREQNESS_StabilityIndex import _moving_median


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _synthetic_freq(*, nsubjects: int = 1, modulated: bool = False):
    sampling_rate = 100.0
    ntime = 4000
    time = np.arange(ntime) / sampling_rate
    frequencies = np.array([5.0, 10.0])
    widths = np.array([4.0, 5.0])
    time_series = np.zeros((2, ntime, frequencies.size, nsubjects), dtype=float)

    for subject in range(nsubjects):
        for frequency_index, frequency in enumerate(frequencies):
            for component in range(2):
                phase = 2.0 * np.pi * frequency * time + 0.2 * subject
                if modulated:
                    modulation_depth = 0.15 if component == 0 else 0.75
                    phase += modulation_depth * np.sin(2.0 * np.pi * 0.4 * time)
                time_series[component, :, frequency_index, subject] = np.sin(
                    phase + 0.1 * component
                )

    return {
        "ts": time_series,
        "frex": frequencies,
        "fwhm": widths,
        "srate": sampling_rate,
    }


def test_constant_sine_has_near_zero_stability_and_no_padded_zero():
    FREQ = _synthetic_freq()

    result = FREQNESS_StabilityIndex(
        FREQ,
        frex2model=5,
        keep_frequency_sliding=True,
        show=False,
    )

    assert isinstance(result, FREQNESSStabilityIndexResult)
    assert result.stability_index.shape == (1, 1, 1)
    assert result.mean_frequency.shape == (1, 1, 1)
    assert result.frequency_sliding is not None
    assert result.frequency_sliding_raw is not None
    assert result.frequency_sliding.shape == (1, 3999, 1, 1)
    assert result.frequency_sliding_raw.shape == (1, 3999, 1, 1)
    assert result.time.shape == (3999,)
    assert result.time[0] == pytest.approx(0.01)
    assert result.mean_frequency[0, 0, 0] == pytest.approx(5.0, abs=1e-10)
    assert result.stability_index[0, 0, 0] < 1e-9
    assert result.frequency_sliding_raw[0, 0, 0, 0] == pytest.approx(
        5.0,
        abs=1e-9,
    )


def test_frequency_modulation_increases_stability_index():
    FREQ = _synthetic_freq(modulated=True)

    result = FREQNESS_StabilityIndex(
        FREQ,
        which_comp=[1, 2],
        frex2model=[5, 10],
        median_window=0.4,
        show=False,
    )

    assert result.stability_index.shape == (2, 2, 1)
    np.testing.assert_array_equal(result.comps, [1, 2])
    np.testing.assert_array_equal(result.frex, [5, 10])
    np.testing.assert_array_equal(result.fwhm, [4, 5])
    assert result.median_samples == 40
    assert result.median_window == pytest.approx(0.4)
    assert np.all(result.stability_index[1] > result.stability_index[0])
    assert result.frequency_sliding is None
    assert result.frequency_sliding_raw is None


def test_even_moving_median_matches_matlab_window_and_endpoints():
    values = np.array([[1.0, 100.0, 2.0, 3.0, 4.0, 5.0]])

    actual = _moving_median(values, 4)

    np.testing.assert_allclose(actual, [[50.5, 2.0, 2.5, 3.5, 3.5, 4.0]])


def test_sampling_rate_controls_window_in_seconds():
    FREQ = _synthetic_freq()
    result = FREQNESS_StabilityIndex(
        FREQ,
        median_window=0.405,
        show=False,
    )

    assert result.median_samples == 41
    assert result.median_window == pytest.approx(0.41)


def test_closest_frequency_warns_and_component_order_is_preserved():
    FREQ = _synthetic_freq(nsubjects=2)

    with pytest.warns(UserWarning, match="closest available frequency: 10.000"):
        result = FREQNESS_StabilityIndex(
            FREQ,
            which_comp=[2, 1, 2],
            frex2model=9.6,
            show=False,
        )

    np.testing.assert_array_equal(result.comps, [2, 1])
    np.testing.assert_array_equal(result.frex, [10])
    assert result.stability_index.shape == (2, 1, 2)


def test_plot_flags_return_group_and_participant_figures():
    FREQ = _synthetic_freq(nsubjects=2, modulated=True)

    result = FREQNESS_StabilityIndex(
        FREQ,
        which_comp=[1, 2],
        plot_avg=True,
        plot_all=True,
        show=False,
    )

    assert len(result.figures) == 6
    assert all(figure.axes[0].get_ylabel() == "Stability Index (Hz)" for figure in result.figures)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"which_comp": 3}, "between 1 and 2"),
        ({"frex2model": [5, 7, 10]}, "one frequency or a two-element range"),
        ({"median_window": 0}, "positive finite duration"),
        ({"median_window": 100}, "instantaneous frequency contains"),
        ({"keep_frequency_sliding": 1}, "must be a boolean"),
    ],
)
def test_invalid_options_raise_clear_errors(kwargs, message):
    FREQ = _synthetic_freq()

    with pytest.raises((TypeError, ValueError), match=message):
        FREQNESS_StabilityIndex(FREQ, show=False, **kwargs)


def test_python_alias_and_module_are_public():
    assert stability_index is FREQNESS_StabilityIndex
    module = importlib.import_module("freqness.FREQNESS_StabilityIndex")
    assert module.FREQNESS_StabilityIndex is FREQNESS_StabilityIndex

