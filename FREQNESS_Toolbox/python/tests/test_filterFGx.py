import numpy as np

from freqness import filterFGx


def test_filter_matches_direct_frequency_domain_formula() -> None:
    rng = np.random.default_rng(42)
    data = rng.standard_normal((3, 256))

    filtered, empirical, kernel = filterFGx(data, 128, 10, 2)
    expected = 2 * np.real(
        np.fft.ifft(np.fft.fft(data, axis=1) * kernel[np.newaxis, :], axis=1)
    )

    np.testing.assert_allclose(filtered, expected, rtol=0, atol=1e-12)
    assert filtered.shape == data.shape
    assert empirical.shape == (3,)
    assert kernel.shape == (data.shape[1],)
    assert np.max(kernel) == 1


def test_filter_preserves_one_dimensional_input() -> None:
    time = np.arange(512) / 128
    signal = np.sin(2 * np.pi * 8 * time)

    filtered, _, _ = filterFGx(signal, 128, 8, 1)

    assert filtered.ndim == 1
    assert filtered.shape == signal.shape

