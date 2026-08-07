from importlib import import_module
from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat

from freqness import (
    ENTROPY_LANDSCAPE,
    FREQNESS_MainPipeline,
    FREQNESSPipelineConfig,
    FREQNESSResult,
    NETWORK_ESTIMATION,
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
    assert "comparison" not in source.lower()
    assert "checkpoint" not in source.lower()
