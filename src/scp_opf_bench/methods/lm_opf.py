"""LM-OPF — Linearized meshed OPF (Wang, Shrestha, Dubey, 2024/2025).

Linearized branch-flow active+reactive KCL with cycle-angle recovery:
``sum_{(i,j) in C} (pi/4)(X_ij P_ij - R_ij Q_ij) = 0`` for each cycle C.
"""

from __future__ import annotations

from ..result import MethodResult
from ._runner import measure_lm as _measure


def measure(name: str) -> MethodResult:
    r = _measure(name, 1.0)
    return MethodResult(
        method="LM-OPF",
        network=name,
        est_mw=float(r["est_mw"]),
        eval_mw=float(r["loss_mw"]),
        delta_v_pct=float(r["delta_v"]),
        runtime_s=float(r["runtime"]),
        status=str(r["status"]),
    )
