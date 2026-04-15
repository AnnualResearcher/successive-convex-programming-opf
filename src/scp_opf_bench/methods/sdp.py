"""SDP — BIM semidefinite relaxation via CVXPY/MOSEK (small networks only)."""

from __future__ import annotations

from ..result import MethodResult
from ._runner import measure_sdp as _measure


def measure(name: str) -> MethodResult:
    r = _measure(name, 1.0)
    return MethodResult(
        method="SDP",
        network=name,
        est_mw=float(r["est_mw"]),
        eval_mw=float(r["loss_mw"]),
        delta_v_pct=float(r["delta_v"]),
        runtime_s=float(r["runtime"]),
        status=str(r["status"]),
    )
