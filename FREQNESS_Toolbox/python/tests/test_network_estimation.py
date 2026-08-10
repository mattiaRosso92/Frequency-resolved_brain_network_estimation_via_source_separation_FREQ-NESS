import numpy as np
import pytest

from freqness import FREQNESS_EntropyLandscape, FREQNESS_NetworkEstimation


def test_network_estimation_mirrors_matlab_output_dimensions() -> None:
    rng = np.random.default_rng(7)
    data = rng.standard_normal((4, 256, 2))

    result = FREQNESS_NetworkEstimation(
        data,
        [10, 5],
        64,
        fwidth=[2, 1],
        ncomps=3,
        bad_segments=[1, 25],
    )

    assert result.evals.shape == (4, 2, 2)
    assert result.evecs.shape == (4, 3, 2, 2)
    assert result.pats.shape == (4, 3, 2, 2)
    assert result.ts.shape == (3, 256, 2, 2)
    np.testing.assert_array_equal(result.frex, [5, 10])
    np.testing.assert_array_equal(result.bad_segments, [1, 25])
    np.testing.assert_allclose(np.sum(result.evals, axis=0), 100.0)
    assert np.all(np.diff(result.evals, axis=0) <= 0)

    for subject in range(2):
        for frequency in range(2):
            expected_ts = (
                result.evecs[:, :, frequency, subject].T @ data[:, :, subject]
            )
            np.testing.assert_allclose(
                result.ts[:, :, frequency, subject], expected_ts, atol=1e-12
            )
            largest = np.argmax(
                np.abs(result.pats[:, :, frequency, subject]), axis=0
            )
            oriented = result.pats[
                largest, np.arange(3), frequency, subject
            ]
            assert np.all(oriented >= 0)


def test_single_participant_keeps_subject_dimension() -> None:
    rng = np.random.default_rng(11)
    result = FREQNESS_NetworkEstimation(
        rng.standard_normal((3, 128)),
        [6],
        64,
        fwidth=1,
        ncomps=2,
    )

    assert result.evals.shape == (3, 1, 1)
    assert result.evecs.shape == (3, 2, 1, 1)
    assert result.ts.shape == (2, 128, 1, 1)


def test_complete_eigenspectrum_is_independent_of_ncomps() -> None:
    rng = np.random.default_rng(19)
    data = rng.standard_normal((4, 256, 2))
    options = {"fwidth": [1, 2], "bad_segments": [1, 25]}

    retained_two = FREQNESS_NetworkEstimation(
        data,
        [5, 10],
        64,
        ncomps=2,
        **options,
    )
    retained_four = FREQNESS_NetworkEstimation(
        data,
        [5, 10],
        64,
        ncomps=4,
        **options,
    )

    np.testing.assert_allclose(retained_two.evals, retained_four.evals)
    assert retained_two.evals.shape == retained_four.evals.shape == (4, 2, 2)
    assert retained_two.evecs.shape[1] == retained_two.pats.shape[1] == 2
    assert retained_two.ts.shape[0] == 2
    assert retained_four.evecs.shape[1] == retained_four.pats.shape[1] == 4
    assert retained_four.ts.shape[0] == 4

    H2_two, ED_two = FREQNESS_EntropyLandscape(retained_two, show=False)
    H2_four, ED_four = FREQNESS_EntropyLandscape(retained_four, show=False)
    np.testing.assert_allclose(H2_two, H2_four)
    np.testing.assert_allclose(ED_two, ED_four)


def test_bad_segments_are_matlab_one_based() -> None:
    data = np.arange(60, dtype=float).reshape(3, 20)

    with pytest.raises(ValueError, match="outside"):
        FREQNESS_NetworkEstimation(
            data,
            [2],
            20,
            fwidth=1,
            ncomps=2,
            bad_segments=[0],
        )
