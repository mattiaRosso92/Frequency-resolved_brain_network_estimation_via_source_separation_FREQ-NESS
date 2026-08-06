import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from freqness import FREQNESS_CrossCoupling, FREQNESSCrossCouplingResult
from freqness.FREQNESS_CrossCoupling import _first_harmonic_fit


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _synthetic_freq(nsubjects: int = 2, *, include_patterns: bool = False):
    sampling_rate = 100.0
    duration = 24.0
    time = np.arange(int(sampling_rate * duration)) / sampling_rate
    frequencies = np.array([2.0, 10.0, 16.0])
    time_series = np.empty((1, time.size, frequencies.size, nsubjects), dtype=float)
    for subject in range(nsubjects):
        lfo_phase = 2 * np.pi * frequencies[0] * time - np.pi / 2
        lfo = np.sin(2 * np.pi * frequencies[0] * time)
        time_series[0, :, 0, subject] = lfo
        for index, carrier in enumerate(frequencies[1:], start=1):
            preferred_phase = 0.4 + 0.2 * subject + 0.15 * index
            envelope = 1.2 + 0.55 * np.cos(lfo_phase - preferred_phase)
            time_series[0, :, index, subject] = envelope * np.cos(
                2 * np.pi * carrier * time
            )
    result = {
        "frex": frequencies,
        "fwhm": np.array([0.6, 8.0, 8.0]),
        "srate": sampling_rate,
        "ts": time_series,
    }
    if include_patterns:
        voxels = 8
        patterns = np.empty((voxels, 1, frequencies.size, nsubjects), dtype=float)
        base = np.linspace(0.1, 1.0, voxels)
        for frequency in range(frequencies.size):
            for subject in range(nsubjects):
                patterns[:, 0, frequency, subject] = np.roll(
                    base, frequency + subject
                )
        result["pats"] = patterns
    return result


def test_first_harmonic_fit_recovers_known_parameters():
    phase = np.linspace(-np.pi, np.pi, 37, endpoint=False) + np.pi / 37
    offset = 4.0
    amplitude = 1.5
    preferred_phase = 0.6
    power = offset + amplitude * np.cos(phase - preferred_phase)

    fit = _first_harmonic_fit(phase, power, min_valid_bins=6)

    np.testing.assert_allclose(
        fit.coefficients,
        [
            offset,
            amplitude * np.cos(preferred_phase),
            amplitude * np.sin(preferred_phase),
        ],
        atol=1e-12,
    )
    assert fit.amplitude == pytest.approx(amplitude)
    assert fit.amplitude_normalized == pytest.approx(37.5)
    assert fit.preferred_phase == pytest.approx(preferred_phase)
    assert fit.offset == pytest.approx(offset)
    assert fit.mse == pytest.approx(0, abs=1e-28)
    assert fit.r2 == pytest.approx(1)
    np.testing.assert_allclose(fit.fitted, power, atol=1e-12)


def test_first_harmonic_fit_requires_minimum_valid_bins():
    phase = np.linspace(-np.pi, np.pi, 12, endpoint=False)
    power = 3 + np.cos(phase)
    power[:7] = np.nan

    fit = _first_harmonic_fit(phase, power, min_valid_bins=6)

    assert fit.nvalid == 5
    assert np.isnan(fit.amplitude)
    assert np.all(np.isnan(fit.fitted))


def test_cross_coupling_uses_neural_lfo_and_returns_harmonic_outputs():
    result = FREQNESS_CrossCoupling(
        _synthetic_freq(),
        2.0,
        show=False,
    )

    assert isinstance(result, FREQNESSCrossCouplingResult)
    np.testing.assert_array_equal(result.carrier_frex, [10, 16])
    assert result.lfo_freq == 2.0
    assert result.comp == 1
    assert result.PAC_all.shape == (2, 2, 37)
    assert result.PAC_avg.shape == (2, 37)
    assert result.coefficients.shape == (2, 2, 3)
    assert result.fitted_PAC.shape == result.PAC_all.shape
    assert result.lfo_phase.shape == (2400, 2)
    assert np.all(result.valid_bins == 37)
    assert np.all(result.amplitude_raw > 0)
    assert np.all(result.amplitude_normalized > 0)
    assert np.all(result.r2 > 0.75)
    np.testing.assert_allclose(result.mFrex, 1 / 37)
    assert result.sAmpl is result.amplitude_raw
    assert result.pShift is result.preferred_phase
    assert result.dcOff is result.dc_offset
    assert result.goodFit is result.mse
    assert len(result.figures) == 1


def test_nearest_lfo_warning_and_carrier_range_selection():
    with pytest.warns(UserWarning, match="using 2 Hz"):
        result = FREQNESS_CrossCoupling(
            _synthetic_freq(nsubjects=1),
            2.2,
            frex2model=[9, 12],
            show=False,
        )

    np.testing.assert_array_equal(result.carrier_frex, [10])


def test_empty_optional_arrays_use_defaults():
    result = FREQNESS_CrossCoupling(
        _synthetic_freq(nsubjects=1),
        2,
        MNI=[],
        frex2model=[],
        show=False,
    )

    np.testing.assert_array_equal(result.carrier_frex, [10, 16])
    assert len(result.figures) == 1


def test_plot_all_and_mni_patterns_create_expected_figures():
    FREQ = _synthetic_freq(nsubjects=1, include_patterns=True)
    coordinates = np.column_stack(
        [
            np.linspace(-20, 20, 8),
            np.linspace(-10, 10, 8),
            np.linspace(0, 30, 8),
        ]
    )

    result = FREQNESS_CrossCoupling(
        FREQ,
        2,
        MNI=coordinates,
        plot_all=True,
        show=False,
    )

    assert len(result.figures) == 3
    assert len(result.figures[0].axes) == 2
    assert len(result.figures[1].axes) == 4
    assert "Participant #1" in result.figures[2].axes[0].get_title()


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"which_comp": 2}, "which_comp"),
        ({"frex2model": [1, 2]}, "No carrier"),
        ({"frex2model": [12, 8]}, "ascending"),
        ({"nbins": 5}, "at least 6"),
        ({"nbins": 10, "min_valid_bins": 11}, "min_valid_bins"),
    ],
)
def test_cross_coupling_rejects_invalid_options(kwargs, message):
    with pytest.raises((ValueError, TypeError), match=message):
        FREQNESS_CrossCoupling(
            _synthetic_freq(nsubjects=1),
            2,
            show=False,
            **kwargs,
        )


def test_cross_coupling_rejects_nonpositive_lfo():
    with pytest.raises(ValueError, match="positive finite"):
        FREQNESS_CrossCoupling(
            _synthetic_freq(nsubjects=1),
            0,
            show=False,
        )
