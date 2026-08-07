from types import SimpleNamespace

import numpy as np
import pytest

from freqness import (
    FREQNESS_BackProjection,
    FREQNESS_NetworkRemoval,
    remove_network,
)


def _freq_and_data(*, nsubjects: int = 1, scale_factors=None):
    frequencies = np.array([4.0, 8.0])
    nvoxels = ncomponents = 3
    ntime = 9
    scales = (
        np.ones(nsubjects, dtype=float)
        if scale_factors is None
        else np.asarray(scale_factors, dtype=float)
    )
    eigenvectors = np.empty((3, 3, 2, nsubjects), dtype=float)
    time_series = np.empty((3, ntime, 2, nsubjects), dtype=float)
    original = np.empty((3, ntime, nsubjects), dtype=float)

    for subject in range(nsubjects):
        rng = np.random.default_rng(100 + subject)
        original[:, :, subject] = rng.normal(size=(nvoxels, ntime))
        for frequency in range(frequencies.size):
            filters = np.array(
                [
                    [1.0 + 0.1 * frequency, 0.2, 0.1],
                    [0.1, 1.1 + 0.1 * subject, 0.3],
                    [0.2, 0.1, 0.9 + 0.05 * frequency],
                ]
            )
            eigenvectors[:, :, frequency, subject] = filters
            analyzed = scales[subject] * original[:, :, subject]
            time_series[:, :, frequency, subject] = filters.T @ analyzed

    FREQ = {
        "frex": frequencies,
        "evecs": eigenvectors,
        "ts": time_series,
        "scale_factors": scales,
    }
    return FREQ, original


def test_removing_all_components_from_square_system_leaves_zero():
    FREQ, original = _freq_and_data()
    data = original[:, :, 0]

    clean, removed = FREQNESS_NetworkRemoval(
        FREQ,
        data,
        4,
        comps2remove=[1, 2, 3],
    )

    np.testing.assert_allclose(removed, data, atol=1e-12)
    np.testing.assert_allclose(clean, 0.0, atol=1e-12)
    np.testing.assert_allclose(clean + removed, data, atol=1e-12)


def test_selected_component_matches_back_projection():
    FREQ, original = _freq_and_data()
    data = original[:, :, 0]

    clean, removed = FREQNESS_NetworkRemoval(
        FREQ,
        data,
        8,
        comps2remove=[2],
    )

    expected = FREQNESS_BackProjection(FREQ, 8, comps2project=[2])
    np.testing.assert_allclose(removed, expected, atol=1e-12)
    np.testing.assert_allclose(clean, data - expected, atol=1e-12)
    np.testing.assert_allclose(clean + removed, data, atol=1e-12)


def test_multiple_participants_and_rescaling_return_original_units():
    FREQ, data = _freq_and_data(
        nsubjects=2,
        scale_factors=[100.0, 0.01],
    )

    clean, removed = FREQNESS_NetworkRemoval(
        FREQ,
        data,
        4,
        comps2remove=[1, 2, 3],
    )

    assert clean.shape == data.shape
    assert removed.shape == data.shape
    np.testing.assert_allclose(removed, data, atol=1e-11)
    np.testing.assert_allclose(clean, 0.0, atol=1e-11)


def test_three_dimensional_single_participant_shape_is_preserved():
    FREQ, data = _freq_and_data()

    clean, removed = FREQNESS_NetworkRemoval(FREQ, data, 4)

    assert clean.shape == data.shape
    assert removed.shape == data.shape


def test_empty_component_selection_defaults_to_first_and_alias_matches():
    FREQ, data = _freq_and_data()
    container = SimpleNamespace(**FREQ)

    result = FREQNESS_NetworkRemoval(
        container,
        data[:, :, 0],
        4,
        comps2remove=[],
    )
    aliased = remove_network(container, data[:, :, 0], 4)

    np.testing.assert_allclose(result[0], aliased[0], atol=1e-12)
    np.testing.assert_allclose(result[1], aliased[1], atol=1e-12)


def test_closest_frequency_warns_with_network_removal_argument_name():
    FREQ, data = _freq_and_data()

    with pytest.warns(UserWarning, match="freq2remove.*8.000 Hz"):
        result = FREQNESS_NetworkRemoval(
            FREQ,
            data[:, :, 0],
            7.6,
            comps2remove=[2],
        )

    exact = FREQNESS_NetworkRemoval(
        FREQ,
        data[:, :, 0],
        8,
        comps2remove=[2],
    )
    np.testing.assert_allclose(result[0], exact[0], atol=1e-12)
    np.testing.assert_allclose(result[1], exact[1], atol=1e-12)


@pytest.mark.parametrize(
    ("frequency", "error"),
    [(0, ValueError), (np.nan, ValueError), ([4], TypeError), (True, TypeError)],
)
def test_invalid_frequencies_are_rejected_with_public_name(frequency, error):
    FREQ, data = _freq_and_data()
    with pytest.raises(error, match="freq2remove"):
        FREQNESS_NetworkRemoval(FREQ, data[:, :, 0], frequency)


@pytest.mark.parametrize(
    "components",
    [[0], [4], [1.5], [np.inf], [[1, 2]], [False]],
)
def test_invalid_components_are_rejected_with_public_name(components):
    FREQ, data = _freq_and_data()
    with pytest.raises(ValueError, match="comps2remove"):
        FREQNESS_NetworkRemoval(
            FREQ,
            data[:, :, 0],
            4,
            comps2remove=components,
        )


@pytest.mark.parametrize(
    ("data_factory", "message"),
    [
        (lambda d: d[:2, :, 0], "voxels"),
        (lambda d: d[:, :-1, 0], "timepoints"),
        (lambda d: np.repeat(d, 2, axis=2), "participants"),
        (lambda d: d[:, :, :, np.newaxis], "shape"),
        (lambda d: d[:, :, 0].astype(complex), "real numeric"),
        (lambda d: np.full_like(d[:, :, 0], np.nan), "finite"),
        (lambda d: np.empty((0, 0)), "non-empty"),
    ],
)
def test_invalid_or_mismatched_data_are_rejected(data_factory, message):
    FREQ, data = _freq_and_data()
    invalid = data_factory(data)
    with pytest.raises((ValueError, TypeError), match=message):
        FREQNESS_NetworkRemoval(FREQ, invalid, 4)
