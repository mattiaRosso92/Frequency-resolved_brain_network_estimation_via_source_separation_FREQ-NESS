from importlib import import_module
from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat

from freqness import (
    BACK_PROJECTION,
    ENTROPY_LANDSCAPE,
    FREQNESS_MainPipeline,
    FREQNESSPipelineConfig,
    FREQNESSResult,
    NETWORK_ESTIMATION,
    NETWORK_REMOVAL,
)


pipeline_module = import_module("freqness.FREQNESS_MainPipeline")


def _toolbox(tmp_path: Path, *, empty_second: bool = False) -> Path:
    first = tmp_path / "FREQNESS_Data" / "Dataset_1"
    second = tmp_path / "FREQNESS_Data" / "Dataset_2"
    mni_root = tmp_path / "FREQNESS_MNI_Coordinates"
    first.mkdir(parents=True)
    second.mkdir()
    mni_root.mkdir()
    savemat(first / "participant.mat", {"data": np.ones((3, 20))})
    if not empty_second:
        savemat(second / "participant.mat", {"data": np.full((3, 20), 2.0)})
    savemat(
        mni_root / "MNI.mat",
        {"MNI": np.column_stack([np.arange(3.0)] * 3)},
    )
    return tmp_path


def _freq_result(scale: float = 1.0) -> FREQNESSResult:
    return FREQNESSResult(
        evals=np.full((2, 2, 1), scale),
        evecs=np.full((3, 2, 2, 1), scale),
        pats=np.full((3, 2, 2, 1), scale),
        ts=np.full((2, 20, 2, 1), scale),
        frex=np.array([2.0, 4.0]),
        fwhm=np.array([0.2, 0.4]),
        srate=10.0,
        duration=2.0,
        bad_segments=np.array([], dtype=np.int64),
        regularisation=0.01,
        scale_factors=np.array([1.0]),
    )


def test_pipeline_defaults_mirror_matlab_settings():
    config = FREQNESSPipelineConfig()

    np.testing.assert_allclose(config.frex, np.arange(1.2, 24.0 + 0.6, 1.2))
    assert config.srate == 250
    assert config.landscape_ncomps == 10
    np.testing.assert_array_equal(config.pattern_frex, [2, 8])
    np.testing.assert_array_equal(config.expdk_range2fit, [9, 20])
    np.testing.assert_array_equal(config.freqgrad_frex2model, [8, 12])
    np.testing.assert_array_equal(config.compgrad_comps2model, [1, 5])
    assert config.lfo_freq == 2
    assert config.backproj_freq2project == 8.4
    np.testing.assert_array_equal(config.backproj_comps2project, [1])
    assert config.netrem_freq2remove == 8.4
    np.testing.assert_array_equal(config.netrem_comps2remove, [1])
    assert config.netrem_plot_landscape is True
    assert config.netrem_landscape_ncomps == 3


def test_core_section_processes_every_dataset_and_keeps_outputs_in_memory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    root = _toolbox(tmp_path)
    original_mat_files = set(root.rglob("*.mat"))
    calls: list[float] = []

    def fake_estimation(data, frex, srate, **options):
        calls.append(float(data[0, 0, 0]))
        assert srate == 250
        assert not options
        return _freq_result(calls[-1])

    monkeypatch.setattr(
        pipeline_module, "FREQNESS_NetworkEstimation", fake_estimation
    )
    result = FREQNESS_MainPipeline(
        FREQNESSPipelineConfig(
            toolbox_root=root,
            analyses=(NETWORK_ESTIMATION,),
            show=False,
        )
    )

    assert calls == [1.0, 2.0]
    assert [item.name for item in result.conditions] == [
        "Dataset_1",
        "Dataset_2",
    ]
    assert all(
        isinstance(item.outputs[NETWORK_ESTIMATION], FREQNESSResult)
        for item in result.conditions
    )
    assert set(root.rglob("*.mat")) == original_mat_files


def test_secondary_section_computes_core_prerequisite_without_exporting_it(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    root = _toolbox(tmp_path)
    estimation_calls = 0
    seen_shapes: list[tuple[int, ...]] = []

    def fake_estimation(*args, **kwargs):
        nonlocal estimation_calls
        estimation_calls += 1
        return _freq_result()

    def fake_entropy(FREQ, **kwargs):
        seen_shapes.append(FREQ.evals.shape)
        return np.array([[0.5]]), np.array([[2.0]])

    monkeypatch.setattr(
        pipeline_module, "FREQNESS_NetworkEstimation", fake_estimation
    )
    monkeypatch.setattr(
        pipeline_module, "FREQNESS_EntropyLandscape", fake_entropy
    )
    result = FREQNESS_MainPipeline(
        FREQNESSPipelineConfig(
            toolbox_root=root,
            analyses=(ENTROPY_LANDSCAPE,),
            show=False,
        )
    )

    assert estimation_calls == 2
    assert seen_shapes == [(2, 2, 1), (2, 2, 1)]
    for condition in result.conditions:
        assert set(condition.outputs) == {ENTROPY_LANDSCAPE}
        np.testing.assert_array_equal(
            condition.outputs[ENTROPY_LANDSCAPE]["H2"], [[0.5]]
        )
        np.testing.assert_array_equal(
            condition.outputs[ENTROPY_LANDSCAPE]["ED"], [[2.0]]
        )


def test_empty_dataset_is_retained_and_skipped(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    root = _toolbox(tmp_path, empty_second=True)
    monkeypatch.setattr(
        pipeline_module,
        "FREQNESS_NetworkEstimation",
        lambda *args, **kwargs: _freq_result(),
    )

    with pytest.warns(UserWarning) as warning_record:
        result = FREQNESS_MainPipeline(
            FREQNESSPipelineConfig(
                toolbox_root=root,
                analyses=(NETWORK_ESTIMATION,),
                show=False,
            )
        )

    assert len(warning_record) == 2
    assert result.conditions[0].outputs
    assert not result.conditions[1].outputs


def test_pipeline_warns_when_data_root_has_no_dataset_folders(
    tmp_path: Path,
) -> None:
    (tmp_path / "FREQNESS_Data").mkdir()
    (tmp_path / "FREQNESS_MNI_Coordinates").mkdir()

    with pytest.warns(UserWarning, match="No dataset folders found"):
        result = FREQNESS_MainPipeline(
            FREQNESSPipelineConfig(
                toolbox_root=tmp_path,
                analyses=(NETWORK_ESTIMATION,),
                show=False,
            )
        )

    assert result.conditions == []


def test_backprojection_and_network_removal_sections_are_independent(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    root = _toolbox(tmp_path)
    estimation_inputs: list[np.ndarray] = []
    backprojection_settings: list[tuple[float, list[int]]] = []
    removal_settings: list[tuple[float, list[int]]] = []
    visualization_settings: list[tuple[list[float], int, object]] = []

    def fake_estimation(data, *args, **kwargs):
        estimation_inputs.append(np.asarray(data).copy())
        return _freq_result(float(np.asarray(data)[0, 0, 0]))

    def fake_backprojection(FREQ, frequency, *, comps2project):
        backprojection_settings.append(
            (float(frequency), np.asarray(comps2project).tolist())
        )
        return np.full((3, 20), 0.25)

    def fake_removal(FREQ, data, frequency, *, comps2remove):
        removal_settings.append(
            (float(frequency), np.asarray(comps2remove).tolist())
        )
        removed = np.full_like(data, 0.25, dtype=float)
        return np.asarray(data, dtype=float) - removed, removed

    def fake_visualizer(FREQ, Landscape, Patterns, **kwargs):
        visualization_settings.append(
            (list(Landscape["frex"]), Landscape["ncomps"], Patterns)
        )
        assert kwargs["save_nifti"] is False
        return "clean landscape"

    monkeypatch.setattr(
        pipeline_module, "FREQNESS_NetworkEstimation", fake_estimation
    )
    monkeypatch.setattr(
        pipeline_module, "FREQNESS_BackProjection", fake_backprojection
    )
    monkeypatch.setattr(
        pipeline_module, "FREQNESS_NetworkRemoval", fake_removal
    )
    monkeypatch.setattr(
        pipeline_module, "FREQNESS_Visualizer", fake_visualizer
    )

    result = FREQNESS_MainPipeline(
        FREQNESSPipelineConfig(
            toolbox_root=root,
            frex=np.array([2.0, 4.0]),
            analyses=(BACK_PROJECTION, NETWORK_REMOVAL),
            show=False,
        )
    )

    assert len(estimation_inputs) == 4
    assert backprojection_settings == [(8.4, [1]), (8.4, [1])]
    assert removal_settings == [(8.4, [1]), (8.4, [1])]
    assert visualization_settings == [([2.0, 4.0], 3, None)] * 2
    for condition in result.conditions:
        assert set(condition.outputs) == {BACK_PROJECTION, NETWORK_REMOVAL}
        removal = condition.outputs[NETWORK_REMOVAL]
        assert removal["visualization"] == "clean landscape"
        assert isinstance(removal["FREQ_clean"], FREQNESSResult)
        np.testing.assert_allclose(
            removal["dataClean"] + removal["removedActivity"],
            estimation_inputs[0]
            if condition.name == "Dataset_1"
            else estimation_inputs[2],
        )


@pytest.mark.parametrize(
    "analyses",
    [(), ("UnknownFunction",), (NETWORK_ESTIMATION, NETWORK_ESTIMATION)],
)
def test_pipeline_rejects_invalid_section_selection(analyses):
    with pytest.raises(ValueError):
        FREQNESS_MainPipeline(FREQNESSPipelineConfig(analyses=analyses))


def test_executable_pipeline_is_complete_and_user_focused():
    script = Path(__file__).resolve().parents[1] / "FREQNESS_MainPipeline.py"
    source = script.read_text()
    assert "ANALYSES_TO_RUN = ALL_ANALYSES" in source
    assert "# 1) FREQNESS_NetworkEstimation" in source
    assert "# 7) FREQNESS_CrossCoupling" in source
    assert "# 8) FREQNESS_BackProjection" in source
    assert "# 9) FREQNESS_NetworkRemoval" in source
    assert "comparison" not in source.lower()
    assert "checkpoint" not in source.lower()
