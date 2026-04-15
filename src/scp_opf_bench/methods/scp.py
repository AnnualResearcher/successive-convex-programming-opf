"""SCP — Proposed contraction-certified successive LP/QP in the
voltage-current representation. ``eta=0`` is the Jacobian-free Picard
iteration; ``eta=1`` uses the full sensitivity update.
"""

from __future__ import annotations

from ..result import MethodResult
from ._runner import measure_scp as _measure


def measure(name: str, eta: float = 1.0) -> MethodResult:
    r = _measure(name, float(eta), 1.0)
    return MethodResult(
        method=f"SCP(η={eta})",
        network=name,
        est_mw=float(r["est_mw"]),
        eval_mw=float(r["loss_mw"]),
        delta_v_pct=float(r["delta_v"]),
        runtime_s=float(r["runtime"]),
        status=str(r["status"]),
    )
