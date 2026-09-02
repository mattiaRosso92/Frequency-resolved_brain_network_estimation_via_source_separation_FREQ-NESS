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
from scipy.io import loadmat


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


def FREQNESS_Startup(
    path_home: str | Path | None = None,
) -> tuple[list[NumericArray | None], NumericArray | None, Path]:
    """Load convention-based FREQ-NESS datasets and optional MNI coordinates.

    Parameters
    ----------
    path_home
        Path to ``FREQNESS_Toolbox``. When omitted, the function searches the
        current directory and the installed module's parent directories.

    Returns
    -------
    allData
        One entry per group or condition. Non-empty entries have shape
        ``(variables, time, participants)``. Empty group folders are ``None``.
    MNI
        MNI coordinates with shape ``(voxels, 3)``, or ``None`` when no valid
        coordinate file is available.
    path_home
        Resolved toolbox directory as a :class:`pathlib.Path`.

    Notes
    -----
    Unlike MATLAB startup, this function does not modify import paths. Python
    package installation provides that functionality without runtime side
    effects.
    """
    root = _resolve_toolbox_root(path_home)
    all_data = _load_groups(root / "FREQNESS_Data")
    mni = _load_mni_coordinates(root / "FREQNESS_MNI_Coordinates")
    return all_data, mni, root


__all__ = ["FREQNESS_Startup"]
