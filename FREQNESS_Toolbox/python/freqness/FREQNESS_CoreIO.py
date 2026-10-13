"""Participant CORE-file storage for FREQ-NESS.

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

from dataclasses import fields
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping, Sequence

import numpy as np

from .FREQNESS_NetworkEstimation import FREQNESSResult


def _metadata_json(metadata: Mapping[str, Any]) -> str:
    """Normalize small configuration values without copying scientific arrays."""
    values = {
        name: value.tolist() if isinstance(value, (np.ndarray, np.generic)) else value
        for name, value in metadata.items()
    }
    return json.dumps(values, sort_keys=True, allow_nan=False)


def _core_matches(path: Path, participant: dict, network_settings: dict) -> bool:
    """Check saved metadata without importing FREQ arrays."""
    with np.load(path, allow_pickle=False) as saved:
        required = {"PARTICIPANT", "NETWORK_SETTINGS"} | {
            "FREQ_" + item.name for item in fields(FREQNESSResult)
        }
        return (
            required.issubset(saved.files)
            and saved["PARTICIPANT"].item() == _metadata_json(participant)
            and saved["NETWORK_SETTINGS"].item() == _metadata_json(network_settings)
        )


def FREQNESS_SaveCore(
    path: str | Path,
    FREQ: FREQNESSResult,
    PARTICIPANT: dict,
    NETWORK_SETTINGS: dict,
    *,
    overwrite: bool = False,
) -> None:
    """Save a complete participant FREQ and its input/settings metadata.

    NumPy NPZ files preserve every field and singleton dimension, support large
    arrays, and require no pickle or additional dependency. The final file is
    replaced only after a successful temporary-file save in the same directory.
    Existing results require overwrite=True. This is a Python format, not MAT.
    """
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Core output already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    # getattr retains array references; dataclasses.asdict would copy them.
    contents = {"FREQ_" + item.name: getattr(FREQ, item.name) for item in fields(FREQ)}
    contents["PARTICIPANT"] = _metadata_json(PARTICIPANT)
    contents["NETWORK_SETTINGS"] = _metadata_json(NETWORK_SETTINGS)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".npz", delete=False) as stream:
        temporary_file = Path(stream.name)
    try:
        np.savez_compressed(temporary_file, **contents)
        os.replace(temporary_file, path)
    finally:
        temporary_file.unlink(missing_ok=True)


def FREQNESS_LoadCore(path: str | Path) -> FREQNESSResult:
    """Load one complete Python CORE result for any FREQ-NESS analysis.

    Returns the original FREQNESSResult with all eigenvalues and the retained
    eigenvectors, patterns and time series. Participant axes remain explicit.
    Input and settings metadata are stored separately as JSON strings named
    PARTICIPANT and NETWORK_SETTINGS inside the NPZ archive.
    """
    with np.load(path, allow_pickle=False) as saved:
        values = {item.name: saved["FREQ_" + item.name] for item in fields(FREQNESSResult)}
    for name in ("srate", "duration", "regularisation"):
        values[name] = float(values[name])
    return FREQNESSResult(**values)


def FREQNESS_LoadCoreFields(
    path: str | Path,
    field_names: Sequence[str],
) -> dict[str, np.ndarray]:
    """Load selected FREQ fields without reconstructing the complete result.

    This loader is intended for memory-conscious comparisons across participant
    CORE files. Only the requested NPZ members are decompressed; callers can
    slice frequencies and components before concatenating participant arrays.
    """
    requested = tuple(field_names)
    valid = {item.name for item in fields(FREQNESSResult)}
    if not requested or any(name not in valid for name in requested):
        raise ValueError(
            "field_names must contain one or more valid FREQNESSResult fields"
        )
    if len(set(requested)) != len(requested):
        raise ValueError("field_names must not contain duplicates")
    with np.load(Path(path), allow_pickle=False) as saved:
        missing = [name for name in requested if "FREQ_" + name not in saved.files]
        if missing:
            raise ValueError(
                "CORE file is missing FREQ field(s): " + ", ".join(missing)
            )
        return {name: saved["FREQ_" + name] for name in requested}


def _core_participant(path: str | Path) -> dict[str, Any]:
    """Read only the participant metadata stored in one CORE file."""
    with np.load(Path(path), allow_pickle=False) as saved:
        if "PARTICIPANT" not in saved.files:
            raise ValueError(f"CORE file is missing PARTICIPANT metadata: {path}")
        return json.loads(saved["PARTICIPANT"].item())


__all__ = [
    "FREQNESS_LoadCore",
    "FREQNESS_LoadCoreFields",
    "FREQNESS_SaveCore",
]
