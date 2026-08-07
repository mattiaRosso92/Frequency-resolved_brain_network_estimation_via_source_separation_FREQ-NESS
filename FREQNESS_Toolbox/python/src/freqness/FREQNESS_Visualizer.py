"""Visualize FREQ-NESS network landscapes and spatial activation patterns."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, TYPE_CHECKING
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING:
    from matplotlib.figure import Figure


FloatArray = NDArray[np.float64]


@dataclass(slots=True)
class FREQNESSVisualization:
    """Figures and files produced by :func:`FREQNESS_Visualizer`."""

    landscape_figures: list[Figure]
    pattern_figures: list[Figure]
    frequency_panel_figures: list[Figure]
    nifti_paths: list[Path]
    landscape_frequencies: FloatArray
    pattern_frequencies: FloatArray | None
    normalized_patterns: FloatArray | None
    group_patterns: FloatArray | None


def _field(container: Any, name: str, *, required: bool = True) -> Any:
    if isinstance(container, Mapping):
        if name in container:
            return container[name]
    elif hasattr(container, name):
        return getattr(container, name)
    if required:
        raise ValueError(f"{name} is required")
    return None


def _positive_integer(value: Any, name: str) -> int:
    if not isinstance(value, (int, np.integer)) or isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a positive integer")
    if value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _frequency_vector(
    values: ArrayLike,
    name: str,
    *,
    sort: bool = True,
) -> FloatArray:
    frequencies = np.asarray(values, dtype=float)
    if frequencies.ndim != 1 or frequencies.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional vector")
    if not np.all(np.isfinite(frequencies)):
        raise ValueError(f"{name} must contain finite values")
    return np.sort(frequencies) if sort else frequencies.copy()


def _nearest_frequency_indices(
    available: FloatArray,
    requested: ArrayLike,
    label: str,
) -> tuple[NDArray[np.int64], FloatArray]:
    requested_frequencies = _frequency_vector(requested, f"{label}.frex")
    distances = np.abs(available[:, np.newaxis] - requested_frequencies[np.newaxis, :])
    indices = np.argmin(distances, axis=0).astype(np.int64)
    matched = available[indices]
    if not np.all(np.isclose(requested_frequencies, matched, atol=1e-3, rtol=0.0)):
        warnings.warn(
            f"Some requested frequencies in {label}.frex are not present in "
            "FREQ.frex. Using closest matches instead.",
            UserWarning,
            stacklevel=3,
        )
    return indices, matched


def _canonical_inputs(
    FREQ: Any,
    Landscape: Any,
    Patterns: Any,
) -> tuple[
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    int,
    int,
    NDArray[np.int64],
    FloatArray,
    NDArray[np.int64],
    FloatArray,
]:
    frequencies = _frequency_vector(
        _field(FREQ, "frex"), "FREQ.frex", sort=False
    )
    eigenvalues = np.asarray(_field(FREQ, "evals"), dtype=float)
    if eigenvalues.ndim == 2:
        eigenvalues = eigenvalues[:, :, np.newaxis]
    if eigenvalues.ndim != 3:
        raise ValueError(
            "FREQ.evals must have shape (components, frequencies[, participants])"
        )
    patterns = np.asarray(_field(FREQ, "pats"), dtype=float)
    if patterns.ndim == 3:
        patterns = patterns[:, :, :, np.newaxis]
    if patterns.ndim != 4:
        raise ValueError(
            "FREQ.pats must have shape "
            "(voxels, components, frequencies[, participants])"
        )

    if eigenvalues.shape[1] != frequencies.size or patterns.shape[2] != frequencies.size:
        raise ValueError("FREQ frequency dimensions do not match FREQ.frex")
    if eigenvalues.shape[0] != patterns.shape[1]:
        raise ValueError("FREQ.evals and FREQ.pats have different component counts")
    if eigenvalues.shape[2] != patterns.shape[3]:
        raise ValueError("FREQ.evals and FREQ.pats have different participant counts")

    mni_coordinates = np.asarray(_field(Patterns, "MNI_coords"), dtype=float)
    if mni_coordinates.ndim != 2 or mni_coordinates.shape[1] != 3:
        raise ValueError("Patterns.MNI_coords must have shape (voxels, 3)")
    if mni_coordinates.shape[0] != patterns.shape[0]:
        raise ValueError(
            "Patterns.MNI_coords must contain one row per voxel in FREQ.pats"
        )
    if not np.any(np.all(np.isfinite(mni_coordinates), axis=1)):
        raise ValueError("Patterns.MNI_coords contains no valid coordinates")

    landscape_ncomps = _positive_integer(
        _field(Landscape, "ncomps"), "Landscape.ncomps"
    )
    pattern_ncomps = _positive_integer(
        _field(Patterns, "ncomps"), "Patterns.ncomps"
    )
    if landscape_ncomps > eigenvalues.shape[0]:
        raise ValueError(
            "Landscape.ncomps exceeds the number of estimated networks"
        )
    if pattern_ncomps > patterns.shape[1]:
        raise ValueError("Patterns.ncomps exceeds the number of estimated networks")

    landscape_indices, landscape_frequencies = _nearest_frequency_indices(
        frequencies, _field(Landscape, "frex"), "Landscape"
    )
    pattern_indices, pattern_frequencies = _nearest_frequency_indices(
        frequencies, _field(Patterns, "frex"), "Patterns"
    )
    return (
        frequencies,
        eigenvalues,
        patterns,
        mni_coordinates,
        landscape_ncomps,
        pattern_ncomps,
        landscape_indices,
        landscape_frequencies,
        pattern_indices,
        pattern_frequencies,
    )


def _canonical_landscape_inputs(
    FREQ: Any,
    Landscape: Any,
) -> tuple[FloatArray, FloatArray, int, NDArray[np.int64], FloatArray]:
    """Validate the inputs required for a landscape-only visualization."""
    frequencies = _frequency_vector(
        _field(FREQ, "frex"), "FREQ.frex", sort=False
    )
    eigenvalues = np.asarray(_field(FREQ, "evals"), dtype=float)
    if eigenvalues.ndim == 2:
        eigenvalues = eigenvalues[:, :, np.newaxis]
    if eigenvalues.ndim != 3:
        raise ValueError(
            "FREQ.evals must have shape (components, frequencies[, participants])"
        )
    if eigenvalues.shape[1] != frequencies.size:
        raise ValueError("FREQ.evals frequency dimension does not match FREQ.frex")
    landscape_ncomps = _positive_integer(
        _field(Landscape, "ncomps"), "Landscape.ncomps"
    )
    if landscape_ncomps > eigenvalues.shape[0]:
        raise ValueError(
            "Landscape.ncomps exceeds the number of estimated networks"
        )
    landscape_indices, landscape_frequencies = _nearest_frequency_indices(
        frequencies, _field(Landscape, "frex"), "Landscape"
    )
    return (
        frequencies,
        eigenvalues,
        landscape_ncomps,
        landscape_indices,
        landscape_frequencies,
    )


def _normalise_patterns(patterns: FloatArray) -> FloatArray:
    absolute = np.abs(patterns)
    finite = np.isfinite(absolute)
    safe = np.where(finite, absolute, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        maxima = np.nanmax(safe, axis=0, keepdims=True)
    normalised = np.zeros_like(absolute, dtype=float)
    np.divide(
        absolute,
        maxima,
        out=normalised,
        where=finite & np.isfinite(maxima) & (maxima > 0),
    )
    normalised[~finite] = np.nan
    return normalised


def _normalise_vector(values: FloatArray) -> FloatArray:
    result = np.abs(np.asarray(values, dtype=float)).copy()
    finite = np.isfinite(result)
    maximum = float(np.max(result[finite])) if np.any(finite) else 0.0
    if maximum > 0:
        result[finite] /= maximum
    else:
        result[finite] = 0.0
    return result


def _threshold_mask(values: FloatArray, threshold_sd: float | None) -> NDArray[np.bool_]:
    finite = np.isfinite(values)
    if threshold_sd is None:
        return finite & (values > 0)
    if not np.any(finite):
        return finite
    cutoff = float(np.mean(values[finite]) + threshold_sd * np.std(values[finite]))
    return finite & (values >= cutoff) & (values > 0)


def _equal_3d_axes(ax: Any, coordinates: FloatArray) -> None:
    spans = np.ptp(coordinates, axis=0)
    spans[spans == 0] = 1.0
    ax.set_box_aspect(spans)
    ax.view_init(elev=25, azim=135)
    ax.set_axis_off()


def _landscape_plots(
    eigenvalues: FloatArray,
    frequencies: FloatArray,
    ncomps: int,
    *,
    plot_all: bool,
) -> list[Figure]:
    import matplotlib.pyplot as plt

    nsubs = eigenvalues.shape[2]
    component_colors = plt.get_cmap("viridis")(
        np.linspace(0.12, 0.88, ncomps)
    )
    maximum = float(np.nanmax(eigenvalues))
    figures: list[Figure] = []

    if plot_all or nsubs == 1:
        for subject in range(nsubs):
            figure, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
            for component in range(ncomps):
                axis.plot(
                    frequencies,
                    eigenvalues[component, :, subject],
                    color=component_colors[component],
                    linewidth=3 if component == 0 else 2,
                    label=f"GED component #{component + 1}",
                )
            _style_landscape_axis(
                axis,
                frequencies,
                maximum,
                f"Brain network landscape - Participant #{subject + 1}",
            )
            figures.append(figure)

    if nsubs > 1:
        average = np.mean(eigenvalues[:ncomps], axis=2)
        sem = np.std(eigenvalues[:ncomps], axis=2, ddof=0) / np.sqrt(nsubs)
        figure, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
        for component in range(ncomps):
            axis.errorbar(
                frequencies,
                average[component],
                yerr=sem[component],
                color=component_colors[component],
                linewidth=3 if component == 0 else 2,
                capsize=2,
                label=f"GED component #{component + 1}",
            )
        _style_landscape_axis(
            axis,
            frequencies,
            maximum,
            "Brain network landscape - Grand-average",
        )
        figures.append(figure)
    return figures


def _style_landscape_axis(
    axis: Any,
    frequencies: FloatArray,
    maximum: float,
    title: str,
) -> None:
    if frequencies.size == 1:
        padding = max(abs(float(frequencies[0])) * 0.05, 0.5)
        axis.set_xlim(
            float(frequencies[0]) - padding,
            float(frequencies[0]) + padding,
        )
    else:
        axis.set_xlim(float(frequencies[0]), float(frequencies[-1]))
    axis.set_ylim(0, maximum * 1.1 if maximum > 0 else 1.0)
    axis.set_xticks(frequencies)
    axis.set_xlabel("Frequency (Hz)", fontweight="bold")
    axis.set_ylabel("Explained variance (%)", fontweight="bold")
    axis.set_title(title, fontweight="bold")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(loc="upper right")


def _frequency_colors(frequencies: FloatArray) -> tuple[Any, Any, FloatArray]:
    from matplotlib import colors
    import matplotlib.pyplot as plt

    minimum = float(np.min(frequencies))
    maximum = float(np.max(frequencies))
    if np.isclose(minimum, maximum):
        padding = max(abs(minimum) * 0.05, 0.5)
        normalizer = colors.Normalize(minimum - padding, maximum + padding)
    else:
        normalizer = colors.Normalize(minimum, maximum)
    colormap = plt.get_cmap("viridis")
    return colormap, normalizer, colormap(normalizer(frequencies))


def _pattern_figure(
    coordinates: FloatArray,
    patterns: FloatArray,
    frequencies: FloatArray,
    component: int,
    title: str,
    threshold_sd: float | None,
) -> Figure:
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable

    valid_coordinates = np.all(np.isfinite(coordinates), axis=1)
    brain_coordinates = coordinates[valid_coordinates]
    colormap, normalizer, frequency_colors = _frequency_colors(frequencies)

    figure = plt.figure(figsize=(9, 7), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    axis.scatter(
        brain_coordinates[:, 0],
        brain_coordinates[:, 1],
        brain_coordinates[:, 2],
        s=4,
        c="black",
        alpha=0.15,
        linewidths=0,
        depthshade=False,
    )

    point_coordinates: list[FloatArray] = []
    point_values: list[FloatArray] = []
    point_colors: list[FloatArray] = []
    for frequency_index in range(frequencies.size):
        values = patterns[:, component, frequency_index]
        mask = valid_coordinates & _threshold_mask(values, threshold_sd)
        if not np.any(mask):
            continue
        point_coordinates.append(coordinates[mask])
        point_values.append(values[mask])
        point_colors.append(
            np.repeat(
                frequency_colors[frequency_index][np.newaxis, :],
                np.count_nonzero(mask),
                axis=0,
            )
        )

    if point_values:
        overlay_coordinates = np.concatenate(point_coordinates, axis=0)
        overlay_values = np.concatenate(point_values)
        overlay_colors = np.concatenate(point_colors, axis=0)
        order = np.argsort(overlay_values, kind="stable")
        marker_areas = 15.0 + 85.0 * overlay_values[order]
        axis.scatter(
            overlay_coordinates[order, 0],
            overlay_coordinates[order, 1],
            overlay_coordinates[order, 2],
            s=marker_areas,
            c=overlay_colors[order],
            alpha=0.85,
            linewidths=0,
            depthshade=False,
        )

    scalar_mappable = ScalarMappable(norm=normalizer, cmap=colormap)
    scalar_mappable.set_array(frequencies)
    colorbar = figure.colorbar(scalar_mappable, ax=axis, shrink=0.7, pad=0.02)
    colorbar.set_label("Frequency (Hz)", fontweight="bold")
    colorbar.set_ticks(np.unique(frequencies))
    axis.set_title(title, fontweight="bold")
    _equal_3d_axes(axis, brain_coordinates)
    return figure


def _frequency_panel_figure(
    coordinates: FloatArray,
    patterns: FloatArray,
    frequencies: FloatArray,
    component: int,
    title: str,
    threshold_sd: float | None,
) -> Figure:
    import matplotlib.pyplot as plt

    valid_coordinates = np.all(np.isfinite(coordinates), axis=1)
    brain_coordinates = coordinates[valid_coordinates]
    _, _, frequency_colors = _frequency_colors(frequencies)
    columns = min(3, frequencies.size)
    rows = int(np.ceil(frequencies.size / columns))
    figure = plt.figure(
        figsize=(5.2 * columns, 4.6 * rows),
        constrained_layout=True,
    )
    figure.suptitle(title, fontweight="bold")

    for frequency_index, frequency in enumerate(frequencies):
        axis = figure.add_subplot(rows, columns, frequency_index + 1, projection="3d")
        axis.scatter(
            brain_coordinates[:, 0],
            brain_coordinates[:, 1],
            brain_coordinates[:, 2],
            s=4,
            c="black",
            alpha=0.15,
            linewidths=0,
            depthshade=False,
        )
        values = patterns[:, component, frequency_index]
        mask = valid_coordinates & _threshold_mask(values, threshold_sd)
        if np.any(mask):
            order = np.argsort(values[mask], kind="stable")
            active_coordinates = coordinates[mask][order]
            active_values = values[mask][order]
            axis.scatter(
                active_coordinates[:, 0],
                active_coordinates[:, 1],
                active_coordinates[:, 2],
                s=15.0 + 85.0 * active_values,
                c=np.repeat(
                    frequency_colors[frequency_index][np.newaxis, :],
                    active_values.size,
                    axis=0,
                ),
                alpha=0.85,
                linewidths=0,
                depthshade=False,
            )
        axis.set_title(f"{frequency:g} Hz")
        _equal_3d_axes(axis, brain_coordinates)
    return figure


def _pattern_plots(
    coordinates: FloatArray,
    normalized_patterns: FloatArray,
    group_patterns: FloatArray | None,
    frequencies: FloatArray,
    ncomps: int,
    *,
    plot_all: bool,
    frequency_panels: bool,
    threshold_sd: float | None,
) -> tuple[list[Figure], list[Figure]]:
    nsubs = normalized_patterns.shape[3]
    combined: list[Figure] = []
    panels: list[Figure] = []

    if plot_all or nsubs == 1:
        for subject in range(nsubs):
            subject_patterns = normalized_patterns[:, :, :, subject]
            for component in range(ncomps):
                title = (
                    f"Network topography - Component #{component + 1} - "
                    f"Participant #{subject + 1}"
                )
                combined.append(
                    _pattern_figure(
                        coordinates,
                        subject_patterns,
                        frequencies,
                        component,
                        title,
                        threshold_sd,
                    )
                )
                if frequency_panels:
                    panels.append(
                        _frequency_panel_figure(
                            coordinates,
                            subject_patterns,
                            frequencies,
                            component,
                            f"{title} - Frequency-specific views",
                            threshold_sd,
                        )
                    )

    if group_patterns is not None:
        for component in range(ncomps):
            title = f"Network topography - Component #{component + 1} - Grand-average"
            combined.append(
                _pattern_figure(
                    coordinates,
                    group_patterns,
                    frequencies,
                    component,
                    title,
                    threshold_sd,
                )
            )
            if frequency_panels:
                panels.append(
                    _frequency_panel_figure(
                        coordinates,
                        group_patterns,
                        frequencies,
                        component,
                        f"{title} - Frequency-specific views",
                        threshold_sd,
                    )
                )
    return combined, panels


def _default_template_path() -> Path:
    toolbox_root = Path(__file__).resolve().parents[3]
    return (
        toolbox_root
        / "FREQNESS_ExternalFunctions"
        / "nifti_tools"
        / "MNI152_8mm_brain_diy.nii.gz"
    )


def _nifti_output(
    coordinates: FloatArray,
    normalized_patterns: FloatArray,
    frequencies: FloatArray,
    ncomps: int,
    output_path: Path,
    template_path: Path,
) -> list[Path]:
    try:
        import nibabel as nib
    except ImportError as error:
        raise ImportError(
            "NIfTI export requires nibabel. Install freqness[visualize]."
        ) from error

    if not template_path.is_file():
        raise FileNotFoundError(f"NIfTI template not found: {template_path}")
    template = nib.load(str(template_path))
    shape = tuple(int(value) for value in template.shape[:3])
    inverse_affine = np.linalg.inv(template.affine)
    valid_coordinates = np.all(np.isfinite(coordinates), axis=1)
    voxel_indices = np.rint(
        nib.affines.apply_affine(inverse_affine, coordinates[valid_coordinates])
    ).astype(np.int64)
    in_bounds = np.all(voxel_indices >= 0, axis=1)
    in_bounds &= np.all(voxel_indices < np.asarray(shape), axis=1)
    coordinate_rows = np.flatnonzero(valid_coordinates)[in_bounds]
    voxel_indices = voxel_indices[in_bounds]

    nifti_directory = output_path / "FREQNESS_Output" / "FREQNESS_nifti"
    nifti_directory.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    nsubs = normalized_patterns.shape[3]

    for frequency_index, frequency in enumerate(frequencies):
        frequency_label = f"{frequency:g}"
        for component in range(ncomps):
            for subject in range(nsubs):
                values = normalized_patterns[:, component, frequency_index, subject]
                volume = np.zeros(shape, dtype=np.float32)
                finite = np.isfinite(values[coordinate_rows])
                voxels = voxel_indices[finite]
                volume[voxels[:, 0], voxels[:, 1], voxels[:, 2]] = values[
                    coordinate_rows[finite]
                ]
                path = nifti_directory / (
                    f"ActivationPattern_Frex_{frequency_label}Hz_"
                    f"Comp#{component + 1}_Sub#{subject + 1}.nii.gz"
                )
                image = nib.Nifti1Image(volume, template.affine, template.header.copy())
                image.set_data_dtype(np.float32)
                nib.save(image, str(path))
                paths.append(path)

            if nsubs > 1:
                group_values = np.nanmean(
                    normalized_patterns[:, component, frequency_index, :], axis=1
                )
                volume = np.zeros(shape, dtype=np.float32)
                finite = np.isfinite(group_values[coordinate_rows])
                voxels = voxel_indices[finite]
                volume[voxels[:, 0], voxels[:, 1], voxels[:, 2]] = group_values[
                    coordinate_rows[finite]
                ]
                path = nifti_directory / (
                    f"ActivationPattern_Frex_{frequency_label}Hz_"
                    f"Comp#{component + 1}_GRAND-AVERAGE.nii.gz"
                )
                image = nib.Nifti1Image(volume, template.affine, template.header.copy())
                image.set_data_dtype(np.float32)
                nib.save(image, str(path))
                paths.append(path)
    return paths


def FREQNESS_Visualizer(
    FREQ: Any,
    Landscape: Any,
    Patterns: Any | None,
    *,
    plot_all: bool = False,
    threshold_sd: float | None = 1.0,
    frequency_panels: bool = False,
    save_nifti: bool = True,
    template_path: str | Path | None = None,
    show: bool = True,
) -> FREQNESSVisualization:
    """Visualize network landscapes and MNI-space activation patterns.

    ``FREQ`` may be a :class:`~freqness.FREQNESSResult`, a mapping, or an
    object exposing ``evals``, ``pats``, and ``frex`` attributes. ``Landscape``
    and ``Patterns`` accept MATLAB-like mappings or objects with the same field
    names. Pass ``Patterns=None`` to produce only the network landscape, with
    no MNI requirement or NIfTI output. Pattern colour represents the requested
    frequency; marker size represents normalized activation magnitude.

    Parameters
    ----------
    plot_all
        Plot every participant in addition to group averages. A single
        participant is always plotted.
    threshold_sd
        Only overlay values at least this many standard deviations above each
        pattern's mean. Pass ``None`` to show all positive values.
    frequency_panels
        Additionally produce one frequency-specific subplot per requested
        frequency.
    save_nifti
        Save participant maps and, when applicable, grand-average maps. The
        output root is ``Patterns.path_output``.
    template_path
        NIfTI template override. By default, the MATLAB toolbox's 8-mm template
        is reused as data, without calling any MATLAB NIfTI functions.
    show
        Call :func:`matplotlib.pyplot.show` after creating all outputs.
    """
    for value, name in (
        (plot_all, "plot_all"),
        (frequency_panels, "frequency_panels"),
        (save_nifti, "save_nifti"),
        (show, "show"),
    ):
        if not isinstance(value, (bool, np.bool_)):
            raise TypeError(f"{name} must be a boolean")
    if threshold_sd is not None:
        if not np.isscalar(threshold_sd) or not np.isfinite(threshold_sd):
            raise ValueError("threshold_sd must be a finite scalar or None")
        if threshold_sd < 0:
            raise ValueError("threshold_sd cannot be negative")
        threshold_sd = float(threshold_sd)

    if Patterns is None:
        (
            _frequencies,
            eigenvalues,
            landscape_ncomps,
            landscape_indices,
            landscape_frequencies,
        ) = _canonical_landscape_inputs(FREQ, Landscape)
        selected_eigenvalues = eigenvalues[
            :landscape_ncomps, landscape_indices, :
        ]
        landscape_figures = _landscape_plots(
            selected_eigenvalues,
            landscape_frequencies,
            landscape_ncomps,
            plot_all=bool(plot_all),
        )
        if show:
            import matplotlib.pyplot as plt

            plt.show()
        return FREQNESSVisualization(
            landscape_figures=landscape_figures,
            pattern_figures=[],
            frequency_panel_figures=[],
            nifti_paths=[],
            landscape_frequencies=landscape_frequencies,
            pattern_frequencies=None,
            normalized_patterns=None,
            group_patterns=None,
        )

    (
        _frequencies,
        eigenvalues,
        patterns,
        coordinates,
        landscape_ncomps,
        pattern_ncomps,
        landscape_indices,
        landscape_frequencies,
        pattern_indices,
        pattern_frequencies,
    ) = _canonical_inputs(FREQ, Landscape, Patterns)

    selected_eigenvalues = eigenvalues[
        :landscape_ncomps, landscape_indices, :
    ]
    selected_patterns = patterns[
        :, :pattern_ncomps, pattern_indices, :
    ]
    normalized_patterns = _normalise_patterns(selected_patterns)
    if normalized_patterns.shape[3] > 1:
        group_patterns = np.nanmean(normalized_patterns, axis=3)
        for component in range(pattern_ncomps):
            for frequency_index in range(pattern_frequencies.size):
                group_patterns[:, component, frequency_index] = _normalise_vector(
                    group_patterns[:, component, frequency_index]
                )
    else:
        group_patterns = None

    landscape_figures = _landscape_plots(
        selected_eigenvalues,
        landscape_frequencies,
        landscape_ncomps,
        plot_all=bool(plot_all),
    )
    pattern_figures, frequency_panel_figures = _pattern_plots(
        coordinates,
        normalized_patterns,
        group_patterns,
        pattern_frequencies,
        pattern_ncomps,
        plot_all=bool(plot_all),
        frequency_panels=bool(frequency_panels),
        threshold_sd=threshold_sd,
    )

    nifti_paths: list[Path] = []
    if save_nifti:
        path_output = _field(Patterns, "path_output")
        if not isinstance(path_output, (str, Path)):
            raise TypeError("Patterns.path_output must be a string or pathlib.Path")
        resolved_template = (
            Path(template_path).expanduser()
            if template_path is not None
            else _default_template_path()
        )
        nifti_paths = _nifti_output(
            coordinates,
            normalized_patterns,
            pattern_frequencies,
            pattern_ncomps,
            Path(path_output).expanduser(),
            resolved_template,
        )

    if show:
        import matplotlib.pyplot as plt

        plt.show()

    return FREQNESSVisualization(
        landscape_figures=landscape_figures,
        pattern_figures=pattern_figures,
        frequency_panel_figures=frequency_panel_figures,
        nifti_paths=nifti_paths,
        landscape_frequencies=landscape_frequencies,
        pattern_frequencies=pattern_frequencies,
        normalized_patterns=normalized_patterns,
        group_patterns=group_patterns,
    )


__all__ = ["FREQNESS_Visualizer", "FREQNESSVisualization"]
