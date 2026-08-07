import importlib

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from freqness import (
    FREQNESS_InducedResponses,
    FREQNESSInducedResponsesResult,
    induced_responses,
)
from freqness.FREQNESS_InducedResponses import _matlab_round


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _synthetic_freq(*, nsubjects: int = 1, zero: bool = False):
    sampling_rate = 100.0
    ntime = 2000
    time = np.arange(ntime) / sampling_rate
    frequencies = np.array([5.0, 10.0, 15.0])
    widths = np.array([2.0, 2.5, 3.0])
    events = [501, 1001, 1501]
    time_series = np.zeros((2, ntime, frequencies.size, nsubjects), dtype=float)

    if not zero:
        for subject in range(nsubjects):
            for frequency_index, frequency in enumerate(frequencies):
                for component in range(2):
                    envelope = np.ones(ntime)
                    for event in events:
                        onset = event - 1
                        envelope[onset : onset + 35] = 2.5 + 0.25 * component
                    phase = 0.2 * subject + 0.1 * component
                    time_series[component, :, frequency_index, subject] = (
                        envelope
                        * np.sin(2 * np.pi * frequency * time + phase)
                    )

    return {
        "ts": time_series,
        "frex": frequencies,
        "fwhm": widths,
        "srate": sampling_rate,
    }, events


def test_induced_response_detects_post_event_power_increase():
    FREQ, events = _synthetic_freq()

    result = FREQNESS_InducedResponses(
        FREQ,
        events,
        [-0.5, 0.6],
        [-0.4, -0.1],
        which_comp=1,
        frex2model=[5, 15],
        show=False,
    )

    assert isinstance(result, FREQNESSInducedResponsesResult)
    assert result.power.shape == (1, 111, 3, 1)
    np.testing.assert_array_equal(result.comps, [1])
    np.testing.assert_array_equal(result.events[0], events)
    np.testing.assert_array_equal(result.ntrials, [3])
    post = (result.time >= 0.1) & (result.time <= 0.25)
    baseline = (result.time >= -0.4) & (result.time < -0.1)
    assert float(np.mean(result.power[0, post, :, 0])) > 4.0
    assert abs(float(np.mean(result.power[0, baseline, :, 0]))) < 1.0
    assert result.power_trials is None
    assert result.power_scaled is None
    assert result.figures == []


def test_keep_trials_preserves_participant_specific_trial_counts():
    FREQ, _ = _synthetic_freq(nsubjects=2)
    events = [[501, 1001, 1501], [501, 1501]]

    result = FREQNESS_InducedResponses(
        FREQ,
        events,
        [-0.2, 0.3],
        [-0.2, 0],
        which_comp=[2, 1],
        frex2model=10,
        keep_trials=True,
        show=False,
    )

    assert result.power.shape == (2, 51, 1, 2)
    assert result.power_trials is not None
    assert result.power_trials[0].shape == (2, 51, 1, 3)
    assert result.power_trials[1].shape == (2, 51, 1, 2)
    np.testing.assert_array_equal(result.comps, [2, 1])
    np.testing.assert_array_equal(result.ntrials, [3, 2])
    np.testing.assert_allclose(
        result.power[:, :, :, 0],
        np.mean(result.power_trials[0], axis=3),
        atol=1e-12,
    )


def test_matlab_half_away_from_zero_rounding_sets_epoch_endpoints():
    assert _matlab_round(2.5) == 3
    assert _matlab_round(-2.5) == -3
    FREQ, events = _synthetic_freq()

    result = FREQNESS_InducedResponses(
        FREQ,
        events,
        [-0.025, 0.025],
        [-0.025, 0],
        frex2model=5,
        show=False,
    )

    np.testing.assert_array_equal(result.time, np.arange(-3, 4) / 100)


def test_baseline_is_half_open_and_excludes_event_onset(monkeypatch):
    module = importlib.import_module("freqness.FREQNESS_InducedResponses")
    sampling_rate = 10.0
    time_series = np.ones((1, 30, 1, 1))
    time_series[0, 10, 0, 0] = 10.0
    FREQ = {"ts": time_series, "frex": [2], "fwhm": [1], "srate": sampling_rate}

    def identity_filter(data, *_args, **_kwargs):
        values = np.asarray(data, dtype=float)
        return values, np.ones(3), np.ones(values.shape[-1])

    monkeypatch.setattr(module, "filterFGx", identity_filter)
    monkeypatch.setattr(module, "hilbert", lambda data, axis: np.asarray(data))

    result = FREQNESS_InducedResponses(
        FREQ,
        [11],
        [-0.2, 0.2],
        [-0.2, 0],
        show=False,
    )

    onset = int(np.flatnonzero(result.time == 0)[0])
    assert result.power[0, onset, 0, 0] == pytest.approx(20.0)


def test_boundary_events_are_excluded_and_duplicates_removed():
    FREQ, _ = _synthetic_freq()

    with pytest.warns(UserWarning) as warnings_seen:
        result = FREQNESS_InducedResponses(
            FREQ,
            [10, 501, 501, 1990],
            [-0.5, 0.5],
            [-0.4, 0],
            frex2model=5,
            show=False,
        )

    messages = [str(item.message) for item in warnings_seen]
    assert any("duplicate" in message for message in messages)
    assert any("excluding 2 event" in message for message in messages)
    np.testing.assert_array_equal(result.events[0], [501])
    np.testing.assert_array_equal(result.ntrials, [1])


def test_reversed_frequency_range_is_sorted_and_nearest_limits_warn():
    FREQ, events = _synthetic_freq()

    with pytest.warns(UserWarning, match="closest range"):
        result = FREQNESS_InducedResponses(
            FREQ,
            events,
            [-0.2, 0.2],
            [-0.2, 0],
            frex2model=[14.2, 5.4],
            show=False,
        )

    np.testing.assert_array_equal(result.frex, [5, 10, 15])
    np.testing.assert_array_equal(result.fwhm, [2, 2.5, 3])


def test_scalar_frequency_uses_closest_available_frequency():
    FREQ, events = _synthetic_freq()

    with pytest.warns(UserWarning, match="closest available frequency: 10.000"):
        result = FREQNESS_InducedResponses(
            FREQ,
            events,
            [-0.2, 0.2],
            [-0.2, 0],
            frex2model=9.6,
            show=False,
        )

    np.testing.assert_array_equal(result.frex, [10])


def test_scale_max_matches_matlab_replication_scaling():
    FREQ, events = _synthetic_freq(nsubjects=2)

    result = FREQNESS_InducedResponses(
        FREQ,
        [events, events],
        [-0.2, 0.4],
        [-0.2, 0],
        which_comp=[1, 2],
        scale_max=True,
        show=False,
    )

    assert result.power_scaled is not None
    for subject in range(2):
        for component in range(2):
            maximum = np.nanmax(result.power[component, :, :, subject])
            np.testing.assert_allclose(
                result.power_scaled[component, :, :, subject],
                result.power[component, :, :, subject] / maximum,
                atol=1e-12,
            )


def test_zero_signal_retains_trials_as_nan_after_invalid_baseline():
    FREQ, events = _synthetic_freq(zero=True)

    with pytest.warns(UserWarning, match="invalid baseline power"):
        result = FREQNESS_InducedResponses(
            FREQ,
            events,
            [-0.2, 0.2],
            [-0.2, 0],
            frex2model=5,
            keep_trials=True,
            show=False,
        )

    assert np.all(np.isnan(result.power))
    assert result.power_trials is not None
    assert np.all(np.isnan(result.power_trials[0]))


def test_group_and_participant_plots_use_frequency_coordinates():
    FREQ, events = _synthetic_freq(nsubjects=2)

    result = FREQNESS_InducedResponses(
        FREQ,
        [events, events],
        [-0.2, 0.3],
        [-0.2, 0],
        which_comp=[1, 2],
        plot_avg=True,
        plot_all=True,
        show=False,
    )

    assert len(result.figures) == 6
    assert result.figures[0].axes[0].get_ylabel() == "Frequency (Hz)"
    assert result.figures[0].axes[1].get_ylabel() == "Power change (dB)"
    assert "component #1" in result.figures[0].axes[0].get_title()
    assert "Participant #2" in result.figures[-1].axes[0].get_title()


def test_python_alias_matches_primary_function():
    FREQ, events = _synthetic_freq()
    primary = FREQNESS_InducedResponses(
        FREQ,
        events,
        [-0.1, 0.2],
        [-0.1, 0],
        frex2model=5,
        show=False,
    )
    alias = induced_responses(
        FREQ,
        events,
        [-0.1, 0.2],
        [-0.1, 0],
        frex2model=5,
        show=False,
    )
    np.testing.assert_allclose(primary.power, alias.power, atol=1e-12)


@pytest.mark.parametrize(
    ("epoch", "baseline", "message"),
    [
        ([0, 0], [-0.1, 0], "epoch_window"),
        ([-0.2, 0.2], [0, -0.1], "baseline_window"),
        ([-0.2, 0.2], [-0.3, 0], "entirely within"),
        ([-0.01, 0.01], [-0.001, 0], "does not contain any samples"),
    ],
)
def test_invalid_time_windows_are_rejected(epoch, baseline, message):
    FREQ, events = _synthetic_freq()
    with pytest.raises(ValueError, match=message):
        FREQNESS_InducedResponses(
            FREQ,
            events,
            epoch,
            baseline,
            frex2model=5,
            show=False,
        )


@pytest.mark.parametrize(
    ("kwargs", "error", "message"),
    [
        ({"which_comp": [0]}, ValueError, "which_comp"),
        ({"which_comp": [1.5]}, ValueError, "which_comp"),
        ({"frex2model": [1, 2, 3]}, ValueError, "frex2model"),
        ({"keep_trials": 1}, TypeError, "keep_trials"),
        ({"plot_avg": "yes"}, TypeError, "plot_avg"),
        ({"show": 0}, TypeError, "show"),
    ],
)
def test_invalid_options_are_rejected(kwargs, error, message):
    FREQ, events = _synthetic_freq()
    with pytest.raises(error, match=message):
        FREQNESS_InducedResponses(
            FREQ,
            events,
            [-0.2, 0.2],
            [-0.2, 0],
            **kwargs,
        )


@pytest.mark.parametrize(
    ("events", "message"),
    [
        ([], "integer"),
        ([0, 100], "outside"),
        ([100.5], "integer"),
        ([np.nan], "integer"),
        ([True], "integer"),
    ],
)
def test_invalid_single_participant_events_are_rejected(events, message):
    FREQ, _ = _synthetic_freq()
    with pytest.raises(ValueError, match=message):
        FREQNESS_InducedResponses(
            FREQ,
            events,
            [-0.2, 0.2],
            [-0.2, 0],
            frex2model=5,
            show=False,
        )


def test_multiple_participants_require_one_event_vector_each():
    FREQ, events = _synthetic_freq(nsubjects=2)
    with pytest.raises(ValueError, match="one vector per participant"):
        FREQNESS_InducedResponses(
            FREQ,
            events,
            [-0.2, 0.2],
            [-0.2, 0],
            frex2model=5,
            show=False,
        )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda f: f.pop("fwhm"), "FREQ.fwhm"),
        (lambda f: f.update(frex=[5, 10]), "frequency dimension"),
        (lambda f: f.update(fwhm=[1, 2]), "one value per frequency"),
        (lambda f: f.update(srate=0), "positive finite scalar"),
        (lambda f: f.update(frex=[5, 10, 50]), "Nyquist"),
    ],
)
def test_invalid_freq_fields_are_rejected(mutation, message):
    FREQ, events = _synthetic_freq()
    mutation(FREQ)
    with pytest.raises((ValueError, TypeError), match=message):
        FREQNESS_InducedResponses(
            FREQ,
            events,
            [-0.2, 0.2],
            [-0.2, 0],
            show=False,
        )


def test_nonfinite_unselected_component_is_allowed_but_selected_is_rejected():
    FREQ, events = _synthetic_freq()
    FREQ["ts"][1, 100, 0, 0] = np.nan

    valid = FREQNESS_InducedResponses(
        FREQ,
        events,
        [-0.2, 0.2],
        [-0.2, 0],
        which_comp=1,
        frex2model=5,
        show=False,
    )
    assert np.all(np.isfinite(valid.power))

    with pytest.raises(ValueError, match="non-finite"):
        FREQNESS_InducedResponses(
            FREQ,
            events,
            [-0.2, 0.2],
            [-0.2, 0],
            which_comp=2,
            frex2model=5,
            show=False,
        )
