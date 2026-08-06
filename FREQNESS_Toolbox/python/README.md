# FREQ-NESS for Python

This directory contains the Python implementation of the FREQ-NESS toolbox.
The existing MATLAB implementation remains unchanged in the parent directory.

The first public API preserves the established FREQ-NESS function names:

```python
from freqness import FREQNESS_NetworkEstimation, FREQNESS_Startup

allData, MNI, path_home = FREQNESS_Startup("/path/to/FREQNESS_Toolbox")
FREQ = FREQNESS_NetworkEstimation(allData[0], frex, srate)
```

Python-style aliases are also available:

```python
from freqness import estimate_networks, startup
```

## Development installation

From this directory:

```text
python -m pip install -e ".[test]"
pytest
```

## Citation

Rosso, M., Fernández-Rubio, G., Keller, P. E., Brattico, E., Vuust, P.,
Kringelbach, M. L., & Bonetti, L. (2025). FREQ-NESS Reveals the Dynamic
Reconfiguration of Frequency-Resolved Brain Networks During Auditory
Stimulation. *Advanced Science*, 2413195.
https://doi.org/10.1002/advs.202413195

