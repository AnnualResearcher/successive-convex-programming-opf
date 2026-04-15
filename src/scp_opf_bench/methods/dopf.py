"""DOPF — Decoupled OPF.

Stage 1: active-power LP (Simplified DistFlow on radial / B-theta on meshed).
Stage 2: reactive-power QP at fixed P*, minimizing R(P*^2 + Q^2).

Returns a :class:`~scp_opf_bench.result.MethodResult`.
"""

from __future__ import annotations

from ..result import MethodResult
from ._runner import measure_dcopf as _measure


def measure(name: str) -> MethodResult:
    r = _measure(name, 1.0)
    return MethodResult(
        method="DOPF",
        network=name,
        est_mw=float(r["est_mw"]),
        eval_mw=float(r["loss_mw"]),
        delta_v_pct=float(r["delta_v"]),
        runtime_s=float(r["runtime"]),
        status=str(r["status"]),
    )
