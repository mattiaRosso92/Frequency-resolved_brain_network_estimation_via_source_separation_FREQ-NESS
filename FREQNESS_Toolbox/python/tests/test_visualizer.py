from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pytest

from freqness import FREQNESS_Visualizer, FREQNESSVisualization


def _example_inputs(nsubs: int = 2):
    frequencies = np.array([2.0, 5.0, 10.0])
    eigenvalues = np.empty((2, 3, nsubs), dtype=float)
    patterns = np.empty((6, 2, 3, nsubs), dtype=float)
    for subject in range(nsubs):
        eigenvalues[:, :, subject] = np.array(
            [[30, 25, 20], [15, 14, 12]], dtype=float
        ) + subject
        for component in range(2):
            for frequency in range(3):
                patterns[:, component, frequency, subject] = (
                    np.arange(1, 7, dtype=float)
                    * (component + 1)
                    * (frequency + 1)
                    * (subject + 1)
                )
    coordinates = np.array(
        [
            [-16, -8, 0],
            [-8, 0, 8],
            [0, 8, 16],
            [8, 0, 8],
            [16, -8, 0],
            [np.nan, 0, 0],
        ],
        dtype=float,
    )
    result = SimpleNamespace(
        frex=frequencies,
        evals=eigenvalues,
        pats=patterns,
    )
    landscape = {"frex": [2, 10], "ncomps": 2}
    pattern_settings = {
        "MNI_coords": coordinates,
        "frex": [2, 10],
        "ncomps": 2,
    }
    return result, landscape, pattern_settings


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_visualizer_creates_group_landscape_and_frequency_colored_patterns():
    FREQ, Landscape, Patterns = _example_inputs(nsubs=2)

    result = FREQNESS_Visualizer(
        FREQ,
        Landscape,
        Patterns,
        save_nifti=False,
        show=False,
        threshold_sd=None,
    )

    assert isinstance(result, FREQNESSVisualization)
    assert len(result.landscape_figures) == 1
    assert len(result.pattern_figures) == 2
    assert result.frequency_panel_figures == []
    np.testing.assert_array_equal(result.landscape_frequencies, [2, 10])
    np.testing.assert_array_equal(result.pattern_frequencies, [2, 10])
    assert result.normalized_patterns.shape == (6, 2, 2, 2)
    assert result.group_patterns is not None
    np.testing.assert_allclose(
        np.nanmax(result.group_patterns, axis=0), np.ones((2, 2))
    )

    pattern_axis = result.pattern_figures[0].axes[0]
    assert len(pattern_axis.collections) == 2
    background, overlay = pattern_axis.collections
    np.testing.assert_allclose(background.get_facecolors()[0, :3], [0, 0, 0])
    assert np.unique(overlay.get_facecolors()[:, :3], axis=0).shape[0] == 2
    assert np.ptp(overlay.get_sizes()) > 0
    colorbar_axis = result.pattern_figures[0].axes[1]
    assert colorbar_axis.get_ylabel() == "Frequency (Hz)"


def test_visualizer_supports_landscape_only_without_mni_or_nifti():
    FREQ, Landscape, _ = _example_inputs(nsubs=2)
    landscape_only = SimpleNamespace(frex=FREQ.frex, evals=FREQ.evals)

    result = FREQNESS_Visualizer(
        landscape_only,
        Landscape,
        None,
        save_nifti=True,
        show=False,
    )

    assert len(result.landscape_figures) == 1
    assert result.pattern_figures == []
    assert result.frequency_panel_figures == []
    assert result.nifti_paths == []
    assert result.pattern_frequencies is None
    assert result.normalized_patterns is None
    assert result.group_patterns is None


def test_visualizer_plot_all_and_frequency_panels_include_each_subject():
    FREQ, Landscape, Patterns = _example_inputs(nsubs=2)

    result = FREQNESS_Visualizer(
        FREQ,
        Landscape,
        Patterns,
        plot_all=True,
        frequency_panels=True,
        save_nifti=False,
        show=False,
    )

    assert len(result.landscape_figures) == 3
    assert len(result.pattern_figures) == 6
    assert len(result.frequency_panel_figures) == 6


def test_visualizer_warns_and_uses_nearest_frequencies():
    FREQ, Landscape, Patterns = _example_inputs(nsubs=1)
    Landscape["frex"] = [2.2]
    Patterns["frex"] = [9.8]

    with pytest.warns(UserWarning) as recorded:
        result = FREQNESS_Visualizer(
            FREQ,
            Landscape,
            Patterns,
            save_nifti=False,
            show=False,
        )

    assert len(recorded) == 2
    np.testing.assert_array_equal(result.landscape_frequencies, [2])
    np.testing.assert_array_equal(result.pattern_frequencies, [10])


def test_visualizer_validates_mni_voxel_count():
    FREQ, Landscape, Patterns = _example_inputs(nsubs=1)
    Patterns["MNI_coords"] = Patterns["MNI_coords"][:-1]

    with pytest.raises(ValueError, match="one row per voxel"):
        FREQNESS_Visualizer(
            FREQ,
            Landscape,
            Patterns,
            save_nifti=False,
            show=False,
        )


def test_visualizer_exports_individual_and_group_nifti_files(tmp_path):
    FREQ, Landscape, Patterns = _example_inputs(nsubs=2)
    Landscape["frex"] = [2]
    Landscape["ncomps"] = 1
    Patterns["frex"] = [2]
    Patterns["ncomps"] = 1
    Patterns["path_output"] = tmp_path
    Patterns["MNI_coords"] = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
            [0, 1, 0],
            [0, 2, 0],
            [0, 0, 2],
        ],
        dtype=float,
    )
    template_path = tmp_path / "template.nii.gz"
    nib.save(
        nib.Nifti1Image(np.zeros((4, 4, 4), dtype=np.float32), np.eye(4)),
        template_path,
    )

    result = FREQNESS_Visualizer(
        FREQ,
        Landscape,
        Patterns,
        template_path=template_path,
        show=False,
    )

    assert len(result.nifti_paths) == 3
    assert all(path.is_file() for path in result.nifti_paths)
    individual = nib.load(result.nifti_paths[0])
    assert individual.shape == (4, 4, 4)
    assert individual.get_fdata()[2, 0, 0] > 0
