"""Per-method ``measure(name)`` entry points.

Each module exposes a single ``measure(name)`` (and for SCP, ``measure(name, eta)``)
that returns a ``MethodResult``.
"""

from .dopf import measure as measure_dopf
from .socp import measure as measure_socp
from .sdp import measure as measure_sdp
from .lm_opf import measure as measure_lm_opf
from .nlp import measure as measure_nlp
from .scp import measure as measure_scp

__all__ = [
    "measure_dopf",
    "measure_socp",
    "measure_sdp",
    "measure_lm_opf",
    "measure_nlp",
    "measure_scp",
]
