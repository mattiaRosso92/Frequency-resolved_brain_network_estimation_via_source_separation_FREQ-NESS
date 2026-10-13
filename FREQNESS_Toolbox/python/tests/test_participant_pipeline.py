from dataclasses import fields, replace
from importlib import import_module
import json
from pathlib import Path
import weakref

import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy.io import loadmat, savemat

from freqness import (
    ALL_ANALYSES, BACK_PROJECTION, CROSS_COUPLING, ENTROPY_LANDSCAPE,
    NETWORK_ESTIMATION, NETWORK_REMOVAL,
    FREQNESS_BackProjection, FREQNESS_CrossCoupling, FREQNESS_EntropyLandscape,
    FREQNESS_LoadCore, FREQNESS_MainPipeline, FREQNESS_NetworkEstimation,
    FREQNESS_SaveCore, FREQNESS_Startup, FREQNESSPipelineConfig,
)

pipeline_module = import_module("freqness.FREQNESS_MainPipeline")
startup_module = import_module("freqness.FREQNESS_Startup")
core_module = import_module("freqness.FREQNESS_CoreIO")


@pytest.fixture
def study(tmp_path):
    rng = np.random.default_rng(15)
    time = np.arange(512) / 64
    data = np.sin(2 * np.pi * np.arange(1, 7)[:, None] * time)
    data += 0.15 * rng.standard_normal(data.shape)
    for condition in ("Dataset_1", "Dataset_2", "Empty"):
        (tmp_path / "FREQNESS_Data" / condition).mkdir(parents=True)
    first = tmp_path / "FREQNESS_Data" / "Dataset_1"
    savemat(first / "participant_a.mat", {"data": (1e-5 * data).astype(np.float32)})
    savemat(first / "Participant_B.mat", {"data": data + 0.2 * rng.standard_normal(data.shape)})
    savemat(tmp_path / "FREQNESS_Data" / "Dataset_2" / "participant_a.mat", {"data": 2 * data})
    (tmp_path / "FREQNESS_MNI_Coordinates").mkdir()
    savemat(tmp_path / "FREQNESS_MNI_Coordinates" / "MNI.mat", {"MNI": rng.normal(size=(6, 3))})
    yield FREQNESSPipelineConfig(
        toolbox_root=tmp_path, frex=np.array([4., 8., 12.]), srate=64.,
        network_options={"duration": 4, "bad_segments": [1, 25], "rescale": True},
        show=False, save_nifti=False, landscape_ncomps=5, pattern_frex=np.array([4., 8.]),
        expdk_range2fit=np.array([4., 12.]), freqgrad_frex2model=np.array([4., 12.]),
        compgrad_freq2model=8., lfo_freq=4., backproj_freq2project=8., netrem_freq2remove=8.,
    )
    plt.close("all")


def _batch(config, folder):
    files = sorted(folder.glob("*.mat"), key=lambda path: path.name.casefold())
    data = np.stack([loadmat(path)["data"].astype(float) for path in files], axis=-1)
    freq = FREQNESS_NetworkEstimation(data, config.frex, config.srate,
                                    **{"ncomps": 5, **config.network_options})
    return data, freq


def _assert_same(actual, expected, subi=None):
    for field in fields(actual):
        value = getattr(expected, field.name)
        if subi is not None and field.name in ("evals", "evecs", "pats", "ts", "scale_factors"):
            value = value[..., subi:subi + 1]
        if field.name in ("evecs", "ts"):
            continue
        np.testing.assert_allclose(getattr(actual, field.name), value, rtol=1e-9, atol=1e-10)
    expected_vectors = expected.evecs if subi is None else expected.evecs[..., subi:subi + 1]
    expected_ts = expected.ts if subi is None else expected.ts[..., subi:subi + 1]
    signs = np.sign(np.sum(actual.evecs * expected_vectors, axis=0))
    np.testing.assert_allclose(actual.evecs * signs[None, ...], expected_vectors, rtol=1e-9, atol=1e-10)
    np.testing.assert_allclose(actual.ts * signs[:, None, ...], expected_ts, rtol=1e-9, atol=1e-10)


def test_discovery_loads_only_coordinates_and_keeps_legacy_startup(study, monkeypatch):
    original = startup_module._mat_variables
    reads = []

    def tracked(path):
        reads.append(path)
        assert path.parent.name == "FREQNESS_MNI_Coordinates"
        return original(path)

    monkeypatch.setattr(startup_module, "_mat_variables", tracked)
    all_data, mni, root, folders = FREQNESS_Startup(study.toolbox_root, load_data=False, return_folders=True)
    assert all_data == [None, None, None]
    assert mni.shape == (6, 3)
    assert len(reads) == 1
    assert [path.name for path in folders[0]["files_mat"]] == ["participant_a.mat", "Participant_B.mat"]
    assert folders[0]["data_size"] == (6, 512)
    assert not (root / "FREQNESS_Results").exists()
    monkeypatch.setattr(startup_module, "_mat_variables", original)
    legacy, _, _, legacy_folders = FREQNESS_Startup(root, return_folders=True)
    assert legacy[0].shape == (6, 512, 2)
    assert legacy_folders[0]["files_mat"][0].name == "Participant_B.mat"
    assert legacy[2] is None


@pytest.mark.parametrize("options", [{}, {"filter": "linear"}, {"fwidth": np.array([0.8, 1.6, 2.4])}])
def test_core_matches_batch_and_releases_arrays_between_loads(study, monkeypatch, options):
    study.network_options.update(options)
    refs = []
    original_load = pipeline_module._load_participant
    original_save = pipeline_module.FREQNESS_SaveCore

    def tracked_load(*args):
        assert all(reference() is None for reference in refs)
        value = original_load(*args)
        refs.append(weakref.ref(value))
        return value

    def tracked_save(path, FREQ, *args, **kwargs):
        assert refs[-1]() is None  # raw data released before serialization
        refs.extend(weakref.ref(getattr(FREQ, name)) for name in ("evals", "evecs", "pats", "ts"))
        return original_save(path, FREQ, *args, **kwargs)

    monkeypatch.setattr(pipeline_module, "_load_participant", tracked_load)
    monkeypatch.setattr(pipeline_module, "FREQNESS_SaveCore", tracked_save)
    result = FREQNESS_MainPipeline(study)
    assert all(reference() is None for reference in refs)
    for condition in result.conditions[:2]:
        assert condition.outputs == {}
        _, expected = _batch(study, Path(study.toolbox_root) / "FREQNESS_Data" / condition.name)
        for subi, path in enumerate(condition.core_files):
            actual = FREQNESS_LoadCore(path)
            assert actual.evals.shape == (6, 3, 1)
            assert actual.evecs.shape == actual.pats.shape == (6, 5, 3, 1)
            assert actual.ts.shape == (5, 256, 3, 1)
            _assert_same(actual, expected, subi)
            with np.load(path, allow_pickle=False) as saved:
                assert json.loads(saved["NETWORK_SETTINGS"].item())["ncomps"] == 5
                assert json.loads(saved["PARTICIPANT"].item())["condition"] == condition.name
    assert result.conditions[2].core_files == []


def test_reuse_does_not_load_raw_or_freq_and_settings_require_overwrite(study, monkeypatch):
    result = FREQNESS_MainPipeline(study)
    output = result.conditions[0].core_files[0]
    timestamp = output.stat().st_mtime_ns

    def unexpected(*args, **kwargs):
        pytest.fail("Matching CORE files must be reused without importing participant arrays")

    with monkeypatch.context() as patch:
        patch.setattr(pipeline_module, "_load_participant", unexpected)
        patch.setattr(pipeline_module, "FREQNESS_LoadCore", unexpected)
        FREQNESS_MainPipeline(study)
    assert output.stat().st_mtime_ns == timestamp
    study.network_options["ncomps"] = 2
    with pytest.raises(ValueError, match="overwrite_core=True"):
        FREQNESS_MainPipeline(study)
    FREQNESS_MainPipeline(replace(study, overwrite_core=True))
    actual = FREQNESS_LoadCore(output)
    assert actual.ts.shape[0] == 2
    assert actual.evals.shape[0] == 6


def test_changed_source_and_missing_core_are_rejected(study):
    result = FREQNESS_MainPipeline(study)
    source = Path(study.toolbox_root) / "FREQNESS_Data/Dataset_1/participant_a.mat"
    savemat(source, {"data": np.ones((6, 512))})
    with pytest.raises(ValueError, match="does not match"):
        FREQNESS_MainPipeline(study)
    with pytest.raises(ValueError, match="does not match"):
        FREQNESS_MainPipeline(replace(study, run_core=False, run_group_analyses=True))
    result.conditions[0].core_files[0].unlink()
    with pytest.raises(FileNotFoundError, match="Run Block 1 first"):
        FREQNESS_MainPipeline(replace(study, run_core=False, run_group_analyses=True))


@pytest.mark.parametrize("bad", [np.full((6, 512), np.nan), np.full((6, 512), 1j)])
def test_real_and_finite_checks_are_deferred_until_loading(study, bad):
    source = Path(study.toolbox_root) / "FREQNESS_Data/Dataset_1/participant_a.mat"
    savemat(source, {"data": bad})
    FREQNESS_Startup(study.toolbox_root, load_data=False)
    with pytest.raises(ValueError, match="real finite"):
        FREQNESS_MainPipeline(study)
    assert not list(Path(study.toolbox_root).rglob("*_FREQ.npz"))


def test_group_results_match_batch_and_only_removal_loads_raw(study, monkeypatch):
    FREQNESS_MainPipeline(study)
    original_load = pipeline_module._load_participant
    loaded = []

    def tracked(*args):
        loaded.append(args[0])
        return original_load(*args)

    monkeypatch.setattr(pipeline_module, "_load_participant", tracked)
    config = replace(study, run_core=False, run_group_analyses=True, analyses=ALL_ANALYSES)
    result = FREQNESS_MainPipeline(config)
    assert len(loaded) == 3
    for condition in result.conditions[:2]:
        data, expected = _batch(study, Path(study.toolbox_root) / "FREQNESS_Data" / condition.name)
        outputs = condition.outputs
        _assert_same(outputs[NETWORK_ESTIMATION], expected)
        entropy = FREQNESS_EntropyLandscape(expected, show=False)
        np.testing.assert_allclose(
            outputs[ENTROPY_LANDSCAPE].H2, entropy.H2, atol=1e-12
        )
        np.testing.assert_allclose(
            outputs[ENTROPY_LANDSCAPE].ED, entropy.ED, atol=1e-12
        )
        np.testing.assert_allclose(
            outputs[ENTROPY_LANDSCAPE].frex, entropy.frex, atol=1e-12
        )
        back = FREQNESS_BackProjection(expected, 8., comps2project=[1])
        np.testing.assert_allclose(outputs[BACK_PROJECTION], back, atol=1e-10)
        removal = outputs[NETWORK_REMOVAL]
        np.testing.assert_allclose(removal["dataClean"] + removal["removedActivity"], data[:, :256, :], atol=1e-12)
        cleaned = FREQNESS_NetworkEstimation(removal["dataClean"], config.frex, config.srate,
                                           **{"ncomps": 5, **config.network_options})
        _assert_same(removal["FREQ_clean"], cleaned)
        coupling = FREQNESS_CrossCoupling(expected, 4., show=False)
        np.testing.assert_allclose(outputs[CROSS_COUPLING].PAC_all, coupling.PAC_all, atol=1e-10)
    loaded.clear()
    FREQNESS_MainPipeline(replace(config, analyses=(ENTROPY_LANDSCAPE,)))
    assert loaded == []


def test_failed_save_preserves_previous_file_and_removes_temporary(study, monkeypatch):
    output = FREQNESS_MainPipeline(study).conditions[0].core_files[0]
    freq = FREQNESS_LoadCore(output)
    previous = output.read_bytes()
    with pytest.raises(FileExistsError):
        FREQNESS_SaveCore(output, freq, {}, {})

    def failed_save(path, **contents):
        Path(path).write_bytes(b"interrupted save")
        raise OSError("simulated disk failure")

    monkeypatch.setattr(core_module.np, "savez_compressed", failed_save)
    with pytest.raises(OSError, match="simulated disk failure"):
        FREQNESS_SaveCore(output, freq, {}, {}, overwrite=True)
    assert output.read_bytes() == previous
    assert sorted(path.name for path in output.parent.iterdir()) == [
        "Participant_B_FREQ.npz", "participant_a_FREQ.npz",
    ]


def test_command_line_preserves_config_and_enables_requested_group_block(study, monkeypatch):
    calls = []
    original = pipeline_module.FREQNESS_MainPipeline

    def tracked(config):
        calls.append(config)
        return original(config)

    monkeypatch.setattr(pipeline_module, "FREQNESS_MainPipeline", tracked)
    output_root = Path(study.toolbox_root) / "custom_results"
    assert pipeline_module.main(["--results-directory", str(output_root)], default_config=study) == 0
    assert (output_root / "CORE/Dataset_1/participant_a_FREQ.npz").is_file()
    assert calls[-1].network_options == study.network_options
    assert pipeline_module.main([
        "--analysis", ENTROPY_LANDSCAPE, "--skip-core", "--results-directory", str(output_root),
    ], default_config=study) == 0
    assert calls[-1].run_group_analyses is True
    assert calls[-1].run_core is False
    assert study.run_core is True
    assert study.run_group_analyses is False


def test_notebook_discovers_then_loads_only_selected_participant(study, monkeypatch):
    notebook = json.loads((Path(__file__).parents[1] / "FREQNESS_NotebookPipeline.ipynb").read_text())
    code = ["".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code"]
    namespace = {}
    exec(code[0], namespace)
    namespace["config"] = replace(study, network_options={"ncomps": 5, **study.network_options})
    raw_reads = []

    def tracked(path):
        raw_reads.append(path)
        return loadmat(path)

    namespace["loadmat"] = tracked
    original = startup_module._mat_variables

    def coordinates_only(path):
        assert path.parent.name == "FREQNESS_MNI_Coordinates"
        return original(path)

    monkeypatch.setattr(startup_module, "_mat_variables", coordinates_only)
    exec(next(source for source in code if "folder_data = FREQNESS_Startup" in source), namespace)
    assert raw_reads == []
    exec(next(source for source in code if source.startswith("CONDITION_INDEX")), namespace)
    exec(next(source for source in code if source.startswith("FREQ = FREQNESS_NetworkEstimation")), namespace)
    assert len(raw_reads) == 1
    assert raw_reads[0].name == "participant_a.mat"
    assert namespace["data"].shape == (6, 512, 1)
    assert namespace["FREQ"].ts.shape == (5, 256, 3, 1)
