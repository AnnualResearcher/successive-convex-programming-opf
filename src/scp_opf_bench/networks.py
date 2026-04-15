"""Single source of truth for per-network configuration.

All benchmark methods import `build_net(name)` from here. This avoids the
classic bug where each method silently uses different load scaling, V_slack,
or sgen Q caps and produces non-comparable results.

Per-network policy
------------------
Each network has a `NetPolicy` describing:
  * `v_slack`         : substation voltage [pu] (matches native pandapower case)
  * `load_scale`      : scalar applied to net.load.p_mw and net.load.q_mvar
  * `q_cap_per_sgen`  : if not None, each sgen's apparent-power rating is
                        reset to sqrt(p^2 + q_cap^2) so its reactive
                        capability is bounded by ±q_cap [MVAR]

Why these knobs are needed
--------------------------
The factory networks for IEEE-14 / 118 / 1888 are defined with PV-bus
generators (V-controlled). Our benchmark uses qdispatch mode, which converts
those PV buses to PQ buses (sgens with dispatchable Q). PV-bus voltage
regulation disappears, so transmission cases at native loading become
NR-infeasible. We document the resulting load scaling and V_slack per
network rather than hide it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class NetPolicy:
    v_slack: float = 1.0
    load_scale: float = 1.0
    q_cap_per_sgen: Optional[float] = None  # MVAR; None = no cap


RADIAL = ("ieee33", "ieee123", "ieee8500")
MESHED = ("ieee14", "ieee118", "ieee1888")
ALL_CASES = RADIAL + MESHED


# Centralized per-network configuration. Edit here, not in each method file.
NETWORK_POLICY: dict[str, NetPolicy] = {
    # Radial: native pu, no Q cap
    "ieee33":   NetPolicy(v_slack=1.0,   load_scale=1.0,  q_cap_per_sgen=None),
    "ieee123":  NetPolicy(v_slack=1.0,   load_scale=1.0,  q_cap_per_sgen=None),
    "ieee8500": NetPolicy(v_slack=1.0,   load_scale=1.0,  q_cap_per_sgen=None),
    # Meshed: native V_slack from pandapower; load scaled to ensure
    # NR-feasibility under qdispatch (PV->PQ removes voltage regulation);
    # per-sgen Q capped at ±30 MVAR to mirror realistic inverter limits.
    "ieee14":   NetPolicy(v_slack=1.06,  load_scale=1.0,  q_cap_per_sgen=30.0),
    "ieee118":  NetPolicy(v_slack=1.035, load_scale=0.5,  q_cap_per_sgen=30.0),
    "ieee1888": NetPolicy(v_slack=1.0,   load_scale=0.3,  q_cap_per_sgen=30.0),
}


def build_net(name: str, *, v_slack_override: Optional[float] = None):
    """Build a pandapower net with the benchmark policy applied.

    Parameters
    ----------
    name : str
        One of the entries in ``ALL_CASES``.
    v_slack_override : float, optional
        If provided, overrides the policy's V_slack (used by the OLTC sweep).

    Returns
    -------
    pandapowerNet
        Configured network ready for any baseline method.
    """
    if name not in NETWORK_POLICY:
        raise ValueError(f"Unknown network '{name}'. Choose from {ALL_CASES}.")
    policy = NETWORK_POLICY[name]
    v_slack = v_slack_override if v_slack_override is not None else policy.v_slack

    from . import _legacy  # ensures sys.path + cache dir are set up
    with _legacy.in_legacy_cwd():
        if name == "ieee33":
            from ._legacy.networks import ieee33
            net = ieee33(ppf=True, use_sgen=True)
        elif name == "ieee123":
            from ._legacy.networks import ieee123
            net = ieee123(ppf=True, use_sgen=True)
        elif name == "ieee8500":
            from ._legacy.networks import ieee8500
            net = ieee8500(ppf=True, use_sgen=True)
        elif name == "ieee14":
            from ._legacy.ieee14 import ieee14
            net = ieee14(mode="qdispatch", sgen=True)
        elif name == "ieee118":
            from ._legacy.ieee118 import ieee118
            net = ieee118(mode="qdispatch", sgen=True)
        elif name == "ieee1888":
            from ._legacy.ieee1888 import ieee1888
            net = ieee1888(mode="qdispatch", sgen=True)
        else:  # pragma: no cover
            raise AssertionError("unreachable")

    # Load scaling
    if policy.load_scale != 1.0:
        net.load.p_mw *= policy.load_scale
        net.load.q_mvar *= policy.load_scale
        net.sgen.p_mw *= policy.load_scale

    # Per-sgen Q cap (preserve current sgen.p_mw, set sn = sqrt(p^2 + q_cap^2))
    if policy.q_cap_per_sgen is not None:
        cap = float(policy.q_cap_per_sgen)
        for idx in net.sgen.index:
            p = float(net.sgen.p_mw[idx])
            net.sgen.loc[idx, "sn_mva"] = (p**2 + cap**2) ** 0.5

    # Slack voltage and tap/shunt cleanup
    net.ext_grid.vm_pu = v_slack
    net.ext_grid.va_degree = 0
    if len(net.trafo) > 0:
        net.trafo["tap_pos"] = 0
    if len(net.shunt) > 0:
        net.shunt["in_service"] = False
    return net


def is_radial(name: str) -> bool:
    return name in RADIAL
