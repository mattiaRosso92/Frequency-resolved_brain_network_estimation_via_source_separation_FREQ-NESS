"""Dataset and MNI-coordinate discovery for FREQ-NESS.

If you use this toolbox, please cite:
Rosso, M., Fernández-Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
Kringelbach, M. L., & Bonetti, L. (2025). FREQ-NESS Reveals the Dynamic
Reconfiguration of Frequency-Resolved Brain Networks During Auditory
Stimulation. Advanced Science, 2413195.
https://doi.org/10.1002/advs.202413195

Toolbox authors
---------------
Mattia Rosso - Center for Music in the Brain, Aarhus University -
mattia.rosso@clin.au.dk
Chiara Malvaso - University of Bologna - chiara.malvaso2@unibo.it
Leonardo Bonetti - Center for Music in the Brain, Aarhus University; Centre
for Eudaimonia and Human Flourishing, Linacre College, University of Oxford -
leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
"""

from __future__ import annotations

from pathlib import Path
from typing import TypeAlias
import warnings

import numpy as np
from numpy.typing import NDArray
from scipy.io import loadmat, whosmat


NumericArray: TypeAlias = NDArray[np.number]


def _mat_variables(path: Path) -> dict[str, NumericArray]:
    """Load user variables from a non-v7.3 MATLAB file."""
    try:
        contents = loadmat(path)
    except NotImplementedError as exc:
        raise ValueError(
            f"{path.name} uses MATLAB v7.3/HDF5 format, which is not yet "
            "supported by FREQNESS_Startup. Save it as MATLAB v7 or earlier."
        ) from exc

    return {
        name: value
        for name, value in contents.items()
        if not name.startswith("__")
    }


def _resolve_toolbox_root(path_home: str | Path | None) -> Path:
    if path_home is not None:
        root = Path(path_home).expanduser().resolve()
        if not root.is_dir():
            raise FileNotFoundError(f"FREQ-NESS toolbox directory not found: {root}")
        return root

    candidates = [Path.cwd(), *Path(__file__).resolve().parents]
    for candidate in candidates:
        if (candidate / "FREQNESS_Data").is_dir() and (
            candidate / "FREQNESS_MNI_Coordinates"
        ).is_dir():
            return candidate

    raise FileNotFoundError(
        "Could not locate FREQNESS_Data and FREQNESS_MNI_Coordinates. "
        "Pass the FREQNESS_Toolbox directory to FREQNESS_Startup(path_home)."
    )


def _load_groups(data_directory: Path) -> list[NumericArray | None]:
    if not data_directory.is_dir():
        raise FileNotFoundError(f"FREQ-NESS data directory not found: {data_directory}")

    groups: list[NumericArray | None] = []
    group_directories = sorted(path for path in data_directory.iterdir() if path.is_dir())

    if not group_directories:
        warnings.warn(
            f"No dataset folders found in {data_directory}. Add one folder per "
            "condition containing participant .mat files.",
            UserWarning,
            stacklevel=2,
        )

    for group_directory in group_directories:
        participant_files = sorted(group_directory.glob("*.mat"))
        if not participant_files:
            warnings.warn(
                f"No .mat files found in folder: {group_directory}",
                UserWarning,
                stacklevel=2,
            )
            groups.append(None)
            continue

        participants: list[NumericArray] = []
        reference_shape: tuple[int, int] | None = None

        for participant_file in participant_files:
            variables = _mat_variables(participant_file)
            if len(variables) != 1:
                raise ValueError(
                    f"File {participant_file.name} must contain exactly one data matrix."
                )

            participant = np.asarray(next(iter(variables.values())))
            if participant.ndim != 2 or not np.issubdtype(participant.dtype, np.number):
                raise ValueError(
                    f"File {participant_file.name} must contain one numeric "
                    "variables-by-time matrix."
                )

            if reference_shape is None:
                reference_shape = participant.shape
            elif participant.shape != reference_shape:
                raise ValueError(
                    f"Size mismatch in file {participant_file.name} relative to "
                    f"the first participant in {group_directory}."
                )

            participants.append(participant)

        groups.append(np.stack(participants, axis=2))

    return groups


def _load_mni_coordinates(mni_directory: Path) -> NumericArray | None:
    if not mni_directory.is_dir():
        raise FileNotFoundError(
            f"FREQ-NESS MNI-coordinate directory not found: {mni_directory}"
        )

    mni_files = sorted(mni_directory.glob("*.mat"))
    if not mni_files:
        return None
    if len(mni_files) > 1:
        warnings.warn(
            "Multiple MNI .mat files found. None loaded.",
            UserWarning,
            stacklevel=2,
        )
        return None

    variables = _mat_variables(mni_files[0])
    if len(variables) != 1:
        warnings.warn(
            f"MNI file {mni_files[0].name} ignored: the file must contain "
            "exactly one coordinate matrix.",
            UserWarning,
            stacklevel=2,
        )
        return None

    coordinates = np.asarray(next(iter(variables.values())))
    if coordinates.ndim != 2 or not np.issubdtype(coordinates.dtype, np.number):
        warnings.warn(
            f"MNI file {mni_files[0].name} ignored: coordinates must be numeric.",
            UserWarning,
            stacklevel=2,
        )
        return None

    if coordinates.shape[1] == 3:
        return coordinates
    if coordinates.shape[0] == 3 and coordinates.shape[1] != 3:
        return coordinates.T

    warnings.warn(
        f"MNI file {mni_files[0].name} ignored: coordinate matrix does not "
        "have three columns.",
        UserWarning,
        stacklevel=2,
    )
    return None


def _inspect_groups(data_directory: Path) -> list[dict]:
    """Inspect participant names and dimensions without loading their arrays."""
    if not data_directory.is_dir():
        raise FileNotFoundError(f"FREQ-NESS data directory not found: {data_directory}")
    directories = sorted(
        (path for path in data_directory.iterdir() if path.is_dir()),
        key=lambda path: path.name.casefold(),
    )
    if not directories:
        warnings.warn(
            f"No dataset folders found in {data_directory}. Add one folder per "
            "condition containing participant .mat files.",
            UserWarning,
            stacklevel=2,
        )

    folder_data = []
    numeric_classes = {
        "double", "single", "int8", "uint8", "int16", "uint16",
        "int32", "uint32", "int64", "uint64",
    }
    for this_folder in directories:
        files_mat = sorted(this_folder.glob("*.mat"), key=lambda path: path.name.casefold())
        data_size = None
        if not files_mat:
            warnings.warn(f"No .mat files found in folder: {this_folder}", stacklevel=2)
        for source_file in files_mat:
            try:
                file_info = whosmat(source_file)
            except NotImplementedError as exc:
                raise ValueError(
                    f"{source_file.name} uses MATLAB v7.3/HDF5 format. "
                    "Save input data as MATLAB v7 or earlier."
                ) from exc
            if len(file_info) != 1:
                raise ValueError(f"File {source_file.name} must contain exactly one data matrix.")
            _, shape, numeric_class = file_info[0]
            if numeric_class not in numeric_classes or len(shape) != 2 or min(shape) == 0:
                raise ValueError(
                    f"File {source_file.name} must contain one non-empty numeric "
                    "variables-by-time matrix."
                )
            if data_size is not None and shape != data_size:
                raise ValueError(
                    f"Size mismatch in file {source_file.name} relative to "
                    f"the first participant in {this_folder}."
                )
            data_size = shape
        folder_data.append({
            "name": this_folder.name, "folder": this_folder,
            "files_mat": files_mat, "data_size": data_size,
        })
    return folder_data


def _load_participant(source_file: Path, data_size: tuple[int, int]) -> NDArray[np.float64]:
    """Load and validate one participant, preserving startup's matrix contract."""
    variables = _mat_variables(source_file)
    if len(variables) != 1:
        raise ValueError(f"File {source_file.name} must contain exactly one data matrix.")
    this_data = np.asarray(next(iter(variables.values())))
    if (
        not np.issubdtype(this_data.dtype, np.number)
        or np.iscomplexobj(this_data)
        or this_data.ndim != 2
        or this_data.shape != data_size
        or this_data.size == 0
        or not np.all(np.isfinite(this_data))
    ):
        raise ValueError(
            f"File {source_file.name} must contain one real finite variables-by-time "
            "matrix with the dimensions inspected during startup."
        )
    return np.asarray(this_data, dtype=float)


def FREQNESS_Startup(
    path_home: str | Path | None = None,
    *,
    load_data: bool = True,
    return_folders: bool = False,
) -> tuple:
    """Load or inspect FREQ-NESS datasets and optional MNI coordinates.

    Parameters
    ----------
    path_home
        Path to ``FREQNESS_Toolbox``. When omitted, the function searches the
        current directory and the installed module's parent directories.
    load_data
        False inspects participant files without importing their data arrays.
        The default True preserves the original eager-loading workflow.
    return_folders
        True adds a fourth output containing folder names, ``files_mat`` paths,
        and common ``data_size`` values. The default three-output API is retained.

    Returns
    -------
    allData
        One entry per group or condition. Non-empty entries have shape
        ``(variables, time, participants)``. Empty group folders are ``None``.
        When load_data is False, every entry is None.
    MNI
        MNI coordinates with shape ``(voxels, 3)``, or ``None`` when no valid
        coordinate file is available.
    path_home
        Resolved toolbox directory as a :class:`pathlib.Path`.
    folder_data
        Folder metadata, returned only when return_folders is True.

    Notes
    -----
    Unlike MATLAB startup, this function does not modify import paths. Python
    package installation provides that functionality without runtime side
    effects.
    Output directories are created by MainPipeline when its CORE block runs.
    Discovery checks dimensions and numeric classes; real and finite values
    are checked when each participant is loaded. Inputs remain MATLAB v7 or
    earlier. The discovery workflow sorts names without regard to case.
    """
    if not isinstance(load_data, (bool, np.bool_)) or not isinstance(
        return_folders, (bool, np.bool_)
    ):
        raise TypeError("load_data and return_folders must be booleans")
    root = _resolve_toolbox_root(path_home)
    if load_data:
        print("\nFREQNESS Startup: loading data and MNI coordinates.")
    else:
        print("\nFREQNESS Startup: inspecting data files and loading MNI coordinates.")
    print("\nFREQNESS successfully initialized.")
    print(f"Base directory: {root}\n")
    folder_data = None
    if load_data:
        all_data = _load_groups(root / "FREQNESS_Data")
        if return_folders:
            folder_data = _inspect_groups(root / "FREQNESS_Data")
            folder_data.sort(key=lambda folder: folder["name"])
            for folder in folder_data:
                folder["files_mat"].sort(key=lambda path: path.name)
    else:
        folder_data = _inspect_groups(root / "FREQNESS_Data")
        all_data = [None] * len(folder_data)
    mni = _load_mni_coordinates(root / "FREQNESS_MNI_Coordinates")
    if not load_data and mni is not None:
        if np.iscomplexobj(mni) or np.any(np.isinf(mni)) or any(
            folder["data_size"] is not None and folder["data_size"][0] != mni.shape[0]
            for folder in folder_data
        ):
            warnings.warn("MNI coordinates ignored: invalid coordinates or voxel count mismatch.",
                          stacklevel=2)
            mni = None
    if return_folders:
        return all_data, mni, root, folder_data
    return all_data, mni, root


__all__ = ["FREQNESS_Startup"]
