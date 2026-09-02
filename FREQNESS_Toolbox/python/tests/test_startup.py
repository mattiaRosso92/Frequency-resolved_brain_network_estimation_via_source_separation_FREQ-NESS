from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat

from freqness import FREQNESS_Startup


def test_startup_loads_groups_and_transposes_mni(tmp_path: Path) -> None:
    data_root = tmp_path / "FREQNESS_Data"
    group_a = data_root / "Group_A"
    group_b = data_root / "Group_B"
    mni_root = tmp_path / "FREQNESS_MNI_Coordinates"
    group_a.mkdir(parents=True)
    group_b.mkdir()
    mni_root.mkdir()

    first = np.arange(12, dtype=float).reshape(3, 4)
    second = first + 100
    savemat(group_a / "participant_01.mat", {"data": first})
    savemat(group_a / "participant_02.mat", {"data": second})
    savemat(mni_root / "coordinates.mat", {"MNI": np.arange(15).reshape(3, 5)})

    with pytest.warns(UserWarning, match="No .mat files"):
        all_data, mni, path_home = FREQNESS_Startup(tmp_path)

    assert path_home == tmp_path
    assert len(all_data) == 2
    np.testing.assert_array_equal(all_data[0][:, :, 0], first)
    np.testing.assert_array_equal(all_data[0][:, :, 1], second)
    assert all_data[1] is None
    assert mni is not None
    assert mni.shape == (5, 3)


def test_startup_rejects_multiple_variables(tmp_path: Path) -> None:
    group = tmp_path / "FREQNESS_Data" / "Group_A"
    group.mkdir(parents=True)
    (tmp_path / "FREQNESS_MNI_Coordinates").mkdir()
    savemat(group / "participant.mat", {"first": np.ones((2, 3)), "second": 1})

    with pytest.raises(ValueError, match="exactly one data matrix"):
        FREQNESS_Startup(tmp_path)


def test_startup_warns_when_data_root_has_no_dataset_folders(
    tmp_path: Path,
) -> None:
    (tmp_path / "FREQNESS_Data").mkdir()
    (tmp_path / "FREQNESS_MNI_Coordinates").mkdir()

    with pytest.warns(UserWarning, match="No dataset folders found"):
        all_data, mni, path_home = FREQNESS_Startup(tmp_path)

    assert all_data == []
    assert mni is None
    assert path_home == tmp_path

