"""Physically calibrated Gaussian-filter widths for FREQ-NESS.

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

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


def FREQNESS_ComputeFilterWidths(
    frex: ArrayLike,
    filter: str = "logarithmic",
) -> FloatArray:
    """Compute automatic FWHM values from physical frequencies in Hz.

    The default logarithmic schedule varies selectivity from Q = 7 at the
    lowest requested frequency to Q = 3.5 at the highest. The linear schedule
    uses constant Q = 5. A single unique frequency in logarithmic mode uses
    the geometric midpoint Q = sqrt(7 * 3.5). Output order matches ``frex``.
    """
    raw_frequencies = np.asarray(frex)
    if raw_frequencies.ndim == 0:
        raw_frequencies = raw_frequencies.reshape(1)
    if (
        raw_frequencies.ndim != 1
        or raw_frequencies.size == 0
        or not np.issubdtype(raw_frequencies.dtype, np.number)
        or np.iscomplexobj(raw_frequencies)
    ):
        raise ValueError(
            "frex must be a non-empty vector of positive finite frequencies"
        )

    frequencies = np.asarray(raw_frequencies, dtype=float)
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0):
        raise ValueError(
            "frex must be a non-empty vector of positive finite frequencies"
        )
    if not isinstance(filter, str) or filter.lower() not in {
        "logarithmic",
        "linear",
    }:
        raise ValueError("filter must be either 'logarithmic' or 'linear'")

    if filter.lower() == "linear":
        q_values = np.full(frequencies.shape, 5.0)
    else:
        q_low = 7.0
        q_high = 3.5
        minimum = float(np.min(frequencies))
        maximum = float(np.max(frequencies))

        if minimum == maximum:
            q_values = np.full(
                frequencies.shape,
                np.sqrt(q_low * q_high),
            )
        else:
            log_position = np.log(frequencies / minimum) / np.log(
                maximum / minimum
            )
            q_values = q_low * (q_high / q_low) ** log_position

    return frequencies / q_values


__all__ = ["FREQNESS_ComputeFilterWidths"]
