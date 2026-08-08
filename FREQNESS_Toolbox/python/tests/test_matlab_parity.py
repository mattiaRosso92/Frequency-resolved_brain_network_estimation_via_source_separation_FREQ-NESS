from pathlib import Path

import numpy as np
from scipy.io import loadmat

from freqness import FREQNESS_NetworkEstimation, filterFGx


FIXTURE = Path(__file__).parent / "fixtures" / "matlab_core_reference.mat"


def _align_component_signs(
    actual_vectors: np.ndarray,
    reference_vectors: np.ndarray,
) -> np.ndarray:
    aligned = actual_vectors.copy()
    for subject in range(actual_vectors.shape[3]):
        for frequency in range(actual_vectors.shape[2]):
            for component in range(actual_vectors.shape[1]):
                actual = actual_vectors[:, component, frequency, subject]
                reference = reference_vectors[:, component, frequency, subject]
                if np.dot(actual, reference) < 0:
                    aligned[:, component, frequency, subject] *= -1
    return aligned


def test_filter_matches_matlab_reference() -> None:
    fixture = loadmat(FIXTURE, simplify_cells=True)

    filtered, empirical, kernel = filterFGx(
        fixture["data"][:, :, 0], fixture["srate"], 8, 1.5
    )

    np.testing.assert_allclose(filtered, fixture["filter_data"], atol=1e-11)
    np.testing.assert_allclose(empirical, fixture["filter_empirical"], atol=1e-12)
    np.testing.assert_allclose(kernel, fixture["filter_kernel"], atol=1e-13)


def test_core_outputs_match_matlab_reference() -> None:
    fixture = loadmat(FIXTURE, simplify_cells=True)
    matlab = fixture["FREQ"]

    result = FREQNESS_NetworkEstimation(
        fixture["data"],
        fixture["frex"],
        fixture["srate"],
        fwidth=fixture["fwidth"],
        regularisation=0.01,
        ncomps=3,
        bad_segments=fixture["bad_segments"],
    )

    np.testing.assert_allclose(result.frex, matlab["frex"], atol=0)
    np.testing.assert_allclose(result.fwhm, matlab["fwhm"], atol=0)
    np.testing.assert_allclose(result.evals, matlab["evals"], rtol=1e-10, atol=1e-10)

    aligned_vectors = _align_component_signs(result.evecs, matlab["evecs"])
    np.testing.assert_allclose(
        aligned_vectors,
        matlab["evecs"],
        rtol=1e-9,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        result.pats,
        matlab["pats"],
        rtol=1e-9,
        atol=1e-10,
    )

    aligned_time_series = result.ts.copy()
    for subject in range(result.evecs.shape[3]):
        for frequency in range(result.evecs.shape[2]):
            for component in range(result.evecs.shape[1]):
                if np.dot(
                    result.evecs[:, component, frequency, subject],
                    matlab["evecs"][:, component, frequency, subject],
                ) < 0:
                    aligned_time_series[component, :, frequency, subject] *= -1

    np.testing.assert_allclose(
        aligned_time_series,
        matlab["ts"],
        rtol=1e-9,
        atol=1e-10,
    )

