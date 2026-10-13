from importlib import import_module
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.io import savemat

from freqness import (
    BACK_PROJECTION,
    ENTROPY_LANDSCAPE,
    FREQNESSEntropyLandscapeResult,
    FREQNESS_MainPipeline,
    FREQNESS_LoadCore,
    FREQNESSPipelineConfig,
    FREQNESSResult,
    LANDSCAPE_COMPARISON,
    NETWORK_ESTIMATION,
    NETWORK_REMOVAL,
    SPATIAL_COMPARISON,
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
        evals=np.full((3, 2, 1), scale),
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
    assert config.network_options["ncomps"] == 5
    assert config.run_core is True
    assert config.run_group_analyses is False
    assert config.run_spatial_comparison is False
    assert config.run_landscape_comparison is False
    assert config.overwrite_core is False
    assert config.landscape_ncomps == 10
    np.testing.assert_array_equal(config.pattern_frex, [2.4, 8.4])
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


def test_core_section_saves_every_participant_and_returns_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    root = _toolbox(tmp_path)
    original_mat_files = set(root.rglob("*.mat"))
    calls: list[float] = []

    def fake_estimation(data, frex, srate, **options):
        calls.append(float(data[0, 0, 0]))
        assert srate == 250
        assert options["ncomps"] == 5
        assert data.shape[-1] == 1
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
    assert all(not item.outputs for item in result.conditions)
    for item in result.conditions:
        assert len(item.core_files) == 1
        assert item.core_files[0].is_file()
        assert isinstance(FREQNESS_LoadCore(item.core_files[0]), FREQNESSResult)
    assert set(root.rglob("*.mat")) == original_mat_files


def test_secondary_section_uses_saved_core_prerequisite(
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
        return FREQNESSEntropyLandscapeResult(
            H2=np.array([[0.5]]),
            ED=np.array([[2.0]]),
            frex=np.array([8.0]),
            figures=[],
        )

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
            run_group_analyses=True,
            show=False,
        )
    )

    assert estimation_calls == 2
    assert seen_shapes == [(3, 2, 1), (3, 2, 1)]
    for condition in result.conditions:
        assert set(condition.outputs) == {ENTROPY_LANDSCAPE}
        np.testing.assert_array_equal(
            condition.outputs[ENTROPY_LANDSCAPE].H2, [[0.5]]
        )
        np.testing.assert_array_equal(
            condition.outputs[ENTROPY_LANDSCAPE].ED, [[2.0]]
        )
        np.testing.assert_array_equal(
            condition.outputs[ENTROPY_LANDSCAPE].frex, [8.0]
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

    assert len(warning_record) == 1
    assert result.conditions[0].core_files
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
            run_group_analyses=True,
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
            else estimation_inputs[1],
        )


@pytest.mark.parametrize(
    "analyses",
    [(), ("UnknownFunction",), (NETWORK_ESTIMATION, NETWORK_ESTIMATION)],
)
def test_pipeline_rejects_invalid_section_selection(analyses):
    with pytest.raises(ValueError):
        FREQNESS_MainPipeline(FREQNESSPipelineConfig(analyses=analyses))


def test_selective_comparison_block_builds_thin_matched_structures(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    root = tmp_path
    mni_root = root / "FREQNESS_MNI_Coordinates"
    mni_root.mkdir(parents=True)
    savemat(mni_root / "MNI.mat", {"MNI": np.column_stack([np.arange(3.0)] * 3)})
    for condition, scales in (
        ("Condition_A", {"participant_1": 1.0, "participant_2": 2.0}),
        ("Condition_B", {"participant_1": 10.0, "participant_2": 20.0}),
    ):
        folder = root / "FREQNESS_Data" / condition
        folder.mkdir(parents=True)
        for participant, scale in scales.items():
            savemat(folder / f"{participant}.mat", {"data": np.full((3, 20), scale)})

    monkeypatch.setattr(
        pipeline_module,
        "FREQNESS_NetworkEstimation",
        lambda data, *args, **kwargs: _freq_result(float(data[0, 0, 0])),
    )
    common = dict(
        toolbox_root=root,
        frex=np.array([2.0, 4.0]),
        analyses=(NETWORK_ESTIMATION,),
        show=False,
    )
    FREQNESS_MainPipeline(FREQNESSPipelineConfig(**common))

    original_records = pipeline_module._comparison_records

    def reversed_second_condition(folder, *args):
        records = original_records(folder, *args)
        return list(reversed(records)) if folder["name"] == "Condition_B" else records

    monkeypatch.setattr(
        pipeline_module,
        "_comparison_records",
        reversed_second_condition,
    )
    monkeypatch.setattr(
        pipeline_module,
        "FREQNESS_LoadCore",
        lambda *args, **kwargs: pytest.fail("Block 3 must not load a complete FREQ"),
    )
    captured: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    def fake_spatial(freq_A, freq_B, _mni, **kwargs):
        captured["pats"] = (freq_A.pats.copy(), freq_B.pats.copy())
        assert not hasattr(freq_A, "ts")
        assert kwargs["condition_labels"] == ("Condition_A", "Condition_B")
        return SimpleNamespace(comps=np.array([1]), clusters=[])

    def fake_landscape(freq_A, freq_B, **kwargs):
        captured["evals"] = (freq_A.evals.copy(), freq_B.evals.copy())
        assert not hasattr(freq_A, "ts")
        assert kwargs["sample_labels"] == ("Condition_A", "Condition_B")
        measure = SimpleNamespace(comps=np.array([1]))
        return SimpleNamespace(evals=measure, results={"evals": measure})

    monkeypatch.setattr(pipeline_module, "FREQNESS_SpatialComparison", fake_spatial)
    monkeypatch.setattr(
        pipeline_module,
        "FREQNESS_LandscapeComparison",
        fake_landscape,
    )
    result = FREQNESS_MainPipeline(
        FREQNESSPipelineConfig(
            **common,
            run_core=False,
            run_spatial_comparison=True,
            run_landscape_comparison=True,
            comparison_conditions=("Condition_A", "Condition_B"),
            spatialcomparison_frex=np.array([4.0]),
            spatialcomparison_comps=np.array([2]),
            landscapecomparison_frex=np.array([4.0]),
            landscapecomparison_comps=np.array([2]),
            comparison_plot_results=False,
        )
    )

    assert set(result.comparisons) == {SPATIAL_COMPARISON, LANDSCAPE_COMPARISON}
    pats_A, pats_B = captured["pats"]
    evals_A, evals_B = captured["evals"]
    assert pats_A.shape == (3, 1, 1, 2)
    assert evals_A.shape == (1, 1, 2)
    np.testing.assert_array_equal(pats_A[0, 0, 0], [1.0, 2.0])
    np.testing.assert_array_equal(pats_B[0, 0, 0], [10.0, 20.0])
    np.testing.assert_array_equal(evals_A[0, 0], [1.0, 2.0])
    np.testing.assert_array_equal(evals_B[0, 0], [10.0, 20.0])


def test_selective_comparison_block_uses_thin_secondary_landscape(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    root = tmp_path
    (root / "FREQNESS_MNI_Coordinates").mkdir(parents=True)
    for condition, offset in (("Condition_A", 0.0), ("Condition_B", 10.0)):
        folder = root / "FREQNESS_Data" / condition
        folder.mkdir(parents=True)
        for participant, scale in (("participant_1", 1.0), ("participant_2", 2.0)):
            savemat(
                folder / f"{participant}.mat",
                {"data": np.full((3, 20), scale + offset)},
            )

    monkeypatch.setattr(
        pipeline_module,
        "FREQNESS_NetworkEstimation",
        lambda data, *args, **kwargs: _freq_result(float(data[0, 0, 0])),
    )

    def fake_entropy(FREQ, **kwargs):
        scales = np.asarray(FREQ.evals[0, 0, :], dtype=float)
        return FREQNESSEntropyLandscapeResult(
            H2=np.vstack((scales, scales + 0.5)),
            ED=np.vstack((scales + 1.0, scales + 1.5)),
            frex=np.array([2.0, 4.0]),
            figures=[],
        )

    captured: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    def fake_landscape(entr_A, entr_B, **kwargs):
        assert not hasattr(entr_A, "ts")
        assert set(vars(entr_A)) == {"H2", "frex"}
        captured["H2"] = (entr_A.H2.copy(), entr_B.H2.copy())
        measure = SimpleNamespace(comps=np.array([1]))
        return SimpleNamespace(H2=measure, results={"H2": measure})

    monkeypatch.setattr(pipeline_module, "FREQNESS_EntropyLandscape", fake_entropy)
    monkeypatch.setattr(
        pipeline_module,
        "FREQNESS_LandscapeComparison",
        fake_landscape,
    )
    result = FREQNESS_MainPipeline(
        FREQNESSPipelineConfig(
            toolbox_root=root,
            frex=np.array([2.0, 4.0]),
            analyses=(ENTROPY_LANDSCAPE,),
            run_group_analyses=True,
            run_landscape_comparison=True,
            comparison_conditions=("Condition_A", "Condition_B"),
            landscapecomparison_measure="H2",
            landscapecomparison_frex=np.array([4.0]),
            comparison_plot_results=False,
            show=False,
        )
    )

    assert LANDSCAPE_COMPARISON in result.comparisons
    H2_A, H2_B = captured["H2"]
    assert H2_A.shape == (1, 2)
    np.testing.assert_array_equal(H2_A[0], [1.5, 2.5])
    np.testing.assert_array_equal(H2_B[0], [11.5, 12.5])


def test_executable_pipeline_is_complete_and_user_focused():
    script = Path(__file__).resolve().parents[1] / "FREQNESS_MainPipeline.py"
    source = script.read_text()
    assert "ANALYSES_TO_RUN = ALL_ANALYSES" in source
    assert "# 1) FREQNESS_NetworkEstimation" in source
    assert "# 7) FREQNESS_CrossCoupling" in source
    assert "# 8) FREQNESS_BackProjection" in source
    assert "# 9) FREQNESS_NetworkRemoval" in source
    assert "# 10) SELECTIVE CONDITION/GROUP COMPARISONS" in source
    assert "checkpoint" not in source.lower()
