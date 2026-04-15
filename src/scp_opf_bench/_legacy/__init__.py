"""Legacy source files from the original research codebase.

Kept unmodified to preserve numerical reproducibility. The clean wrappers in
``scp_opf_bench.methods`` import the public functions defined here.

On first import this module:
  1. Inserts the ``_legacy`` directory into ``sys.path`` so the legacy
     modules can import each other by flat module name
     (``from utils import ...`` etc.).
  2. Changes the working directory to ``_legacy`` for the duration of any
     `build_net` / `measure_*` call, because `utils.InitializeConstants`
     reads/writes a pickle cache under the relative path
     ``networkdata/initialize/<name>.pkl``.
"""

from __future__ import annotations

import os as _os
import sys as _sys
from contextlib import contextmanager as _contextmanager

_LEGACY_DIR = _os.path.dirname(_os.path.abspath(__file__))
if _LEGACY_DIR not in _sys.path:
    _sys.path.insert(0, _LEGACY_DIR)

# Make sure the pickle-cache directory exists.
_os.makedirs(_os.path.join(_LEGACY_DIR, "networkdata", "initialize"), exist_ok=True)


@_contextmanager
def in_legacy_cwd():
    """Context manager: run a block with CWD set to the _legacy directory."""
    prev = _os.getcwd()
    try:
        _os.chdir(_LEGACY_DIR)
        yield _LEGACY_DIR
    finally:
        _os.chdir(prev)


__all__ = ["in_legacy_cwd"]
