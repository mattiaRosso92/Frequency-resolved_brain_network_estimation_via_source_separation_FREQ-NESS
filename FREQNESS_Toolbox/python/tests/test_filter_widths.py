import numpy as np

from freqness import (
    FREQNESS_ComputeFilterWidths,
    FREQNESS_CrossCoupling,
    FREQNESS_InducedResponses,
    FREQNESS_NetworkEstimation,
)


def _deterministic_data() -> np.ndarray:
    sampling_rate = 64
    time = np.arange(512) / sampling_rate
    return np.vstack(
        [
            np.sin(2 * np.pi * 2 * time) + 0.20 * np.cos(2 * np.pi * 11 * time),
            np.cos(2 * np.pi * 3 * time) + 0.10 * np.sin(2 * np.pi * 13 * time),
            0.70 * np.sin(2 * np.pi * 4 * time + 0.3)
            + 0.30 * np.cos(2 * np.pi * 7 * time),
            0.40 * np.cos(2 * np.pi * 6 * time - 0.2)
            + 0.25 * np.sin(2 * np.pi * 9 * time),
        ]
    )


def test_logarithmic_reference_calibration() -> None:
    actual = FREQNESS_ComputeFilterWidths([1, 10, 100])
    expected = [
        0.142857142857143,
        2.020305089104422,
        28.571428571428573,
    ]

    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=0)


def test_automatic_widths_are_frequency_density_invariant() -> None:
    sparse = np.arange(1.0, 100.0 + 1.0, 1.0)
    dense = np.arange(1.0, 100.0 + 0.25, 0.25)

    sparse_widths = FREQNESS_ComputeFilterWidths(sparse)
    dense_widths = FREQNESS_ComputeFilterWidths(dense)

    np.testing.assert_array_equal(sparse_widths, dense_widths[::4])


def test_linear_mode_uses_constant_q() -> None:
    frequencies = np.array([1.0, 10.0, 100.0])

    actual = FREQNESS_ComputeFilterWidths(frequencies, "linear")

    np.testing.assert_array_equal(actual, frequencies / 5)


def test_single_frequency_uses_logarithmic_midpoint_q() -> None:
    actual = FREQNESS_ComputeFilterWidths(10)

    np.testing.assert_allclose(actual, [10 / np.sqrt(7 * 3.5)], atol=1e-12)


def test_network_estimation_is_frequency_density_invariant() -> None:
    data = _deterministic_data()
    sparse = FREQNESS_NetworkEstimation(data, [2, 4, 6], 64, ncomps=3)
    dense = FREQNESS_NetworkEstimation(
        data,
        [2, 3, 4, 5, 6],
        64,
        ncomps=3,
    )
    shared = slice(None, None, 2)

    np.testing.assert_array_equal(sparse.fwhm, dense.fwhm[shared])
    np.testing.assert_allclose(sparse.evals, dense.evals[:, shared, :], atol=1e-12)
    np.testing.assert_allclose(
        sparse.evecs,
        dense.evecs[:, :, shared, :],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        sparse.pats,
        dense.pats[:, :, shared, :],
        atol=1e-12,
    )
    np.testing.assert_allclose(sparse.ts, dense.ts[:, :, shared, :], atol=1e-12)


def test_explicit_widths_override_automatic_mode() -> None:
    data = _deterministic_data()
    frequencies = [2, 4, 6]
    widths = [0.4, 0.8, 1.2]

    logarithmic = FREQNESS_NetworkEstimation(
        data,
        frequencies,
        64,
        fwidth=widths,
        filter="logarithmic",
        ncomps=3,
    )
    linear = FREQNESS_NetworkEstimation(
        data,
        frequencies,
        64,
        fwidth=widths,
        filter="linear",
        ncomps=3,
    )
    scalar = FREQNESS_NetworkEstimation(
        data,
        frequencies,
        64,
        fwidth=0.8,
        ncomps=1,
    )

    np.testing.assert_array_equal(logarithmic.fwhm, widths)
    np.testing.assert_array_equal(linear.fwhm, widths)
    np.testing.assert_array_equal(scalar.fwhm, [0.8, 0.8, 0.8])
    np.testing.assert_array_equal(logarithmic.evals, linear.evals)
    np.testing.assert_array_equal(logarithmic.evecs, linear.evecs)
    np.testing.assert_array_equal(logarithmic.pats, linear.pats)
    np.testing.assert_array_equal(logarithmic.ts, linear.ts)


def test_downstream_functions_reuse_physical_widths() -> None:
    result = FREQNESS_NetworkEstimation(
        _deterministic_data(),
        [2, 4, 6],
        64,
        ncomps=2,
    )

    coupling = FREQNESS_CrossCoupling(
        result,
        2,
        frex2model=[4, 6],
        which_comp=1,
        plot_all=False,
        nbins=12,
        min_valid_bins=6,
        show=False,
    )
    np.testing.assert_array_equal(coupling.carrier_frex, [4, 6])

    induced = FREQNESS_InducedResponses(
        result,
        events=[160, 320],
        epoch_window=[-0.5, 0.5],
        baseline_window=[-0.5, 0],
        which_comp=1,
        frex2model=[2, 6],
        plot_avg=False,
        plot_all=False,
        show=False,
    )
    np.testing.assert_array_equal(induced.frex, [2, 4, 6])
    np.testing.assert_array_equal(induced.fwhm, result.fwhm)
