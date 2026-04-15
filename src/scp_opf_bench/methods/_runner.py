"""Internal: load the legacy ``measure_runtime`` module, wrapping each
exposed function so it runs with CWD set to the ``_legacy`` directory
(required by ``utils.InitializeConstants`` which reads/writes a pickle
cache at ``networkdata/initialize/<name>.pkl``).
"""
from __future__ import annotations

from functools import wraps

from ..import _legacy  # side-effect: adds _legacy dir to sys.path, creates cache dir
import measure_runtime as _mr  # noqa: E402  (post-path-munging import)


def _in_legacy_cwd(fn):
    @wraps(fn)
    def wrapper(*args, **kwargs):
        with _legacy.in_legacy_cwd():
            return fn(*args, **kwargs)
    return wrapper


measure_dcopf = _in_legacy_cwd(_mr.measure_dcopf)
measure_socp = _in_legacy_cwd(_mr.measure_socp)
measure_sdp = _in_legacy_cwd(_mr.measure_sdp)
measure_lm = _in_legacy_cwd(_mr.measure_lm)
measure_nlp = _in_legacy_cwd(_mr.measure_nlp)
measure_scp = _in_legacy_cwd(_mr.measure_scp)
RADIAL = tuple(_mr.RADIAL)
MESHED = tuple(_mr.MESHED)
