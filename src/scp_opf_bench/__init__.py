"""SCP-OPF Benchmark — reproducibility package for the SCP-OPF paper."""

from .networks import RADIAL, MESHED, ALL_CASES, build_net, NETWORK_POLICY
from .nr import nr_validate, voltage_deviation_pct

__all__ = [
    "RADIAL",
    "MESHED",
    "ALL_CASES",
    "build_net",
    "NETWORK_POLICY",
    "nr_validate",
    "voltage_deviation_pct",
]

__version__ = "0.1.0"
