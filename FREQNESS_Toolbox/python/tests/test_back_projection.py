from types import SimpleNamespace

import numpy as np
import pytest
from scipy.linalg import pinv

from freqness import FREQNESS_BackProjection, back_project


def _square_freq(*, nsubjects: int = 1, scale_factors=None):
    frequencies = np.array([4.0, 8.0])
    nvoxels = ncomponents = 3
    ntime = 7
    eigenvectors = np.empty(
        (nvoxels, ncomponents, frequencies.size, nsubjects),
        dtype=float,
    )
    time_series = np.empty(
        (ncomponents, ntime, frequencies.size, nsubjects),
        dtype=float,
    )
    original_data = np.empty((nvoxels, ntime, nsubjects), dtype=float)
    scales = (
        np.ones(nsubjects, dtype=float)
        if scale_factors is None
        else np.asarray(scale_factors, dtype=float)
    )

    for subject in range(nsubjects):
        original = (
            np.arange(nvoxels * ntime, dtype=float).reshape(nvoxels, ntime)
            + 0.25
            + subject
        )
        original_data[:, :, subject] = original
        for frequency in range(frequencies.size):
            filters = np.array(
                [
                    [1.0 + 0.1 * frequency, 0.2, 0.1],
                    [0.1, 1.2 + 0.1 * subject, 0.3],
                    [0.2, 0.1, 0.9 + 0.05 * frequency],
                ]
            )
            eigenvectors[:, :, frequency, subject] = filters
            analyzed_data = scales[subject] * original
            time_series[:, :, frequency, subject] = filters.T @ analyzed_data

    FREQ = {
        "frex": frequencies,
        "evecs": eigenvectors,
        "ts": time_series,
        "scale_factors": scales,
    }
    return FREQ, original_data


def test_all_components_exactly_reconstruct_square_system():
    FREQ, original = _square_freq(nsubjects=2)

    result = FREQNESS_BackProjection(
        FREQ,
        4,
        comps2project=[1, 2, 3],
    )

    assert result.shape == original.shape
    np.testing.assert_allclose(result, original, atol=1e-12)


def test_single_participant_output_is_two_dimensional():
    FREQ, original = _square_freq()

    result = FREQNESS_BackProjection(
        FREQ,
        8,
        comps2project=[1, 2, 3],
    )

    assert result.shape == original[:, :, 0].shape
    np.testing.assert_allclose(result, original[:, :, 0], atol=1e-12)


def test_original_units_are_restored_after_internal_rescaling():
    FREQ, original = _square_freq(nsubjects=2, scale_factors=[100.0, 0.25])

    result = FREQNESS_BackProjection(
        FREQ,
        4,
        comps2project=[1, 2, 3],
    )

    np.testing.assert_allclose(result, original, atol=1e-11)


def test_selected_components_use_complete_forward_model_and_one_based_indices():
    rng = np.random.default_rng(12)
    filters = rng.normal(size=(5, 3))
    component_ts = rng.normal(size=(3, 20))
    FREQ = {
        "frex": [6.0],
        "evecs": filters[:, :, np.newaxis],
        "ts": component_ts[:, :, np.newaxis],
    }

    result = FREQNESS_BackProjection(
        FREQ,
        6,
        comps2project=[3, 1, 3],
    )

    expected = pinv(filters.T)[:, [2, 0]] @ component_ts[[2, 0], :]
    np.testing.assert_allclose(result, expected, atol=1e-12)


def test_all_retained_components_reconstruct_retained_subspace_projection():
    rng = np.random.default_rng(21)
    filters = rng.normal(size=(6, 3))
    data = rng.normal(size=(6, 30))
    FREQ = {
        "frex": [5.0],
        "evecs": filters[:, :, np.newaxis],
        "ts": (filters.T @ data)[:, :, np.newaxis],
    }

    result = FREQNESS_BackProjection(FREQ, 5, comps2project=[1, 2, 3])

    expected = pinv(filters.T) @ filters.T @ data
    np.testing.assert_allclose(result, expected, atol=1e-12)
    assert not np.allclose(result, data)


def test_reconstruction_is_invariant_to_paired_filter_and_time_series_sign_flip():
    FREQ, _ = _square_freq()
    baseline = FREQNESS_BackProjection(FREQ, 4, comps2project=[2])

    flipped = {name: np.array(value, copy=True) for name, value in FREQ.items()}
    flipped["evecs"][:, 1, :, :] *= -1
    flipped["ts"][1, :, :, :] *= -1

    result = FREQNESS_BackProjection(flipped, 4, comps2project=[2])

    np.testing.assert_allclose(result, baseline, atol=1e-12)


def test_closest_frequency_warns_and_uses_expected_solution():
    FREQ, _ = _square_freq()

    with pytest.warns(UserWarning, match="closest available frequency: 8.000 Hz"):
        result = FREQNESS_BackProjection(
            FREQ,
            7.7,
            comps2project=[1, 2, 3],
        )

    exact = FREQNESS_BackProjection(FREQ, 8, comps2project=[1, 2, 3])
    np.testing.assert_allclose(result, exact, atol=1e-12)


def test_dataclass_style_input_defaults_to_first_component_and_alias_matches():
    FREQ, _ = _square_freq()
    container = SimpleNamespace(**FREQ)

    result = FREQNESS_BackProjection(container, 4, comps2project=[])
    aliased = back_project(container, 4)

    np.testing.assert_allclose(result, aliased, atol=1e-12)


@pytest.mark.parametrize(
    ("frequency", "error", "message"),
    [
        (0, ValueError, "positive finite"),
        (np.inf, ValueError, "positive finite"),
        ([4], TypeError, "scalar"),
        (True, TypeError, "scalar"),
    ],
)
def test_invalid_requested_frequencies_are_rejected(frequency, error, message):
    FREQ, _ = _square_freq()
    with pytest.raises(error, match=message):
        FREQNESS_BackProjection(FREQ, frequency)


@pytest.mark.parametrize(
    "components",
    [[0], [4], [1.5], [np.nan], [[1, 2]], [True]],
)
def test_invalid_component_selections_are_rejected(components):
    FREQ, _ = _square_freq()
    with pytest.raises(ValueError, match="comps2project"):
        FREQNESS_BackProjection(FREQ, 4, comps2project=components)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda f: f.pop("evecs"), "FREQ.evecs"),
        (lambda f: f.update(frex=[4]), "frequency dimension"),
        (lambda f: f.update(ts=np.repeat(f["ts"], 2, axis=3)), "incompatible"),
        (lambda f: f.update(scale_factors=[1, 2]), "one value per participant"),
        (lambda f: f.update(scale_factors=[0]), "positive"),
    ],
)
def test_incompatible_freq_fields_are_rejected(mutation, message):
    FREQ, _ = _square_freq()
    mutation(FREQ)
    with pytest.raises((ValueError, TypeError), match=message):
        FREQNESS_BackProjection(FREQ, 4)
