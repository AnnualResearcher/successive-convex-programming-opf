"""NR validation utilities — the single source of truth for what 'Eval' means.

Every method's result is validated by running pandapower's Newton–Raphson
power flow with the method's recovered Q dispatch (and any P fixings). The
NR solver tries multiple init strategies because stressed transmission
systems (e.g. IEEE-118, 1888) can be sensitive to the starting point.
"""

from __future__ import annotations

import cmath
import math
from typing import Iterable

import pandapower as pp


def nr_validate(net, buses: Iterable[int]):
    """Run NR with multiple init strategies, return (loss_mw, V_dict, ok).

    The loss accounts for both line and transformer losses. ``V_dict`` maps
    bus index to the complex per-unit voltage (rect form) recovered by NR.
    ``ok`` is False if NR diverges with all init strategies.
    """
    n_buses = sum(1 for _ in buses)
    base_tol = 1e-5 if n_buses < 500 else 1e-4
    inits = [
        ("auto", base_tol),
        ("dc", base_tol),
        ("flat", base_tol),
        ("auto", base_tol * 10),
        ("dc", base_tol * 10),
    ]
    for init, tol in inits:
        try:
            pp.runpp(
                net,
                init=init,
                max_iteration=100,
                tolerance_mva=tol,
                numba=True,
                enforce_q_lims=False,
            )
        except Exception:
            continue
        loss = float(net.res_line.pl_mw.sum() + net.res_trafo.pl_mw.sum())
        V = {
            bus: cmath.rect(
                float(net.res_bus.vm_pu[bus]),
                math.radians(float(net.res_bus.va_degree[bus])),
            )
            for bus in buses
        }
        return loss, V, True
    return float("nan"), {}, False


def voltage_deviation_pct(V_method, V_nr, buses: Iterable[int]) -> float:
    """Average per-bus magnitude deviation ||V_method| - |V_nr|| / |V_nr|, in percent."""
    if not V_method or not V_nr:
        return float("nan")
    err = 0.0
    n = 0
    for bus in buses:
        if bus in V_method and bus in V_nr:
            num = abs(abs(V_method[bus]) - abs(V_nr[bus]))
            err += num / max(abs(V_nr[bus]), 1e-9) * 100.0
            n += 1
    return err / n if n else float("nan")
