# MILP-SCP OPF with discrete OLTC — multi-network driver
# Runs the tap-integer formulation on several test systems for reviewer response.
#%%
import networks
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import pandapower as pp
import numpy as np
import math
import cmath
import copy
import time
import warnings
warnings.filterwarnings('ignore')
from utils import InitializeConstants, RefreshPQ


def compute_bus_shunts(net, Sbase=1e6):
    """Return per-bus real shunt admittance (line charging B/2, trafo magnetizing /
    Π-equivalent B/2, explicit shunts) read directly from the pandapower ppc branch
    matrix column BR_B (index 4).

    For each branch with non-zero B, MATPOWER attaches (j·B/2) at each end (scaled
    by 1/|tap|^2 at the from-end of a tapped trafo).  This function sums those
    contributions at each bus.  The "tap asymmetry" of Yff vs. Ytt is *not* a
    physical shunt and is already handled by the Ohm's law formulation, so it is
    not included here.

    Returns a dict: bus -> complex shunt admittance in system p.u. at `Sbase`.
    """
    import pandapower as pp
    import numpy as np
    bus_shunt = {}
    try:
        pp.runpp(net, init='flat')
    except Exception:
        try:
            pp.runpp(net)
        except Exception:
            return bus_shunt
    if '_ppc' not in net or net._ppc is None:
        return bus_shunt
    branch = net._ppc['branch']
    bus_lookup = net._pd2ppc_lookups.get('bus', None)
    if bus_lookup is None:
        return bus_shunt

    # ppc_sbase vs target Sbase (our convention Y_new = Y_old * (S_old/S_new))
    ppc_sbase_mva = float(net.sn_mva)
    y_scale = (ppc_sbase_mva * 1e6) / Sbase  # scale=1 when system_sbase = net.sn_mva*1e6

    # Build reverse map: ppc_idx -> pp bus index
    ppc_to_pp = {}
    for pp_bus in net.bus.index:
        try:
            ppc_to_pp[int(bus_lookup[pp_bus])] = int(pp_bus)
        except Exception:
            pass

    # Collect per-branch shunt contributions.
    # ppc branch columns: F_BUS=0, T_BUS=1, BR_R=2, BR_X=3, BR_B=4, ..., TAP=8, SHIFT=9
    n_line = len(net.line)
    n_trafo = len(net.trafo)
    for i in range(n_line + n_trafo):
        B = complex(0, branch[i, 4])  # j*B (B is susceptance)
        if abs(B) < 1e-14:
            continue
        f = int(branch[i, 0])
        t = int(branch[i, 1])
        tap = float(branch[i, 8])
        if abs(tap) < 1e-12:
            tap = 1.0
        # At from-bus: (j*B/2) / |tap|^2
        # At to-bus: (j*B/2)
        y_from = B / 2.0 / (tap * tap)
        y_to = B / 2.0
        pp_f = ppc_to_pp.get(f)
        pp_t = ppc_to_pp.get(t)
        if pp_f is not None:
            bus_shunt[pp_f] = bus_shunt.get(pp_f, complex(0, 0)) + y_from * y_scale
        if pp_t is not None:
            bus_shunt[pp_t] = bus_shunt.get(pp_t, complex(0, 0)) + y_to * y_scale

    # Also add any explicit net.shunt elements that are in_service
    try:
        for sidx in net.shunt.index:
            if not net.shunt.loc[sidx, 'in_service']:
                continue
            bus = int(net.shunt.loc[sidx, 'bus'])
            # Shunt Y = (P - jQ) / |V_n|^2 in per-unit; pandapower stores P, Q in MW/MVAr at 1 pu
            # Y_pu = (P + jQ)_mva / S_base_mva  (conductance positive, susceptance negative for capacitive Q<0)
            p_mw = float(net.shunt.loc[sidx, 'p_mw'])
            q_mvar = float(net.shunt.loc[sidx, 'q_mvar'])
            sbase_mva = Sbase / 1e6
            g = p_mw / sbase_mva
            b = -q_mvar / sbase_mva  # negative: capacitive (Q<0 → b>0 susceptance)
            bus_shunt[bus] = bus_shunt.get(bus, complex(0, 0)) + complex(g, b)
    except Exception:
        pass

    return bus_shunt


def compute_multi_vbase_rx(net, Sbase=1e6):
    """Return (R, X) dicts in system-base per-unit, read directly from pandapower's
    internal ppc branch matrix after a Newton-Raphson solve.  This delegates ALL trafo
    parameter handling (vk_percent, vkr_percent, shift_degree, negative leakage,
    multi-voltage ratios, tap_side logic) to pandapower and avoids re-implementing them.

    pandapower stores per-unit R, X on a 1 MVA system base when net.sn_mva = 1
    (our convention), or on net.sn_mva MVA otherwise; we rescale to `Sbase`.
    """
    import pandapower as pp
    import numpy as np

    R, X = {}, {}

    # pandapower needs a runpp call to populate _ppc
    try:
        pp.runpp(net, init='flat')
    except Exception:
        # Try run without flat init in case the solve diverges at tap settings
        try:
            pp.runpp(net)
        except Exception:
            pass

    if '_ppc' not in net or net._ppc is None or 'branch' not in net._ppc:
        # fallback: return empty, caller will keep old values
        return R, X

    branch = net._ppc['branch']
    bus_lookup = net._pd2ppc_lookups.get('bus', None)
    # scaling factor: ppc branch R,X are in per-unit on net.sn_mva system base (MVA)
    # our target system base is Sbase [VA]
    ppc_sbase_mva = float(net.sn_mva) if hasattr(net, 'sn_mva') else 1.0
    scale = (ppc_sbase_mva * 1e6) / Sbase  # Sbase ohm / ppc ohm

    def _get_rx(row):
        r = float(branch[row, 2]) * scale
        x = float(branch[row, 3]) * scale
        return r, x

    # Lines: first N_line rows of branch correspond to net.line in index order
    line_rows = list(range(len(net.line)))
    for pos, lidx in enumerate(net.line.index):
        if not net.line.loc[lidx, 'in_service']:
            continue
        fr = int(net.line.loc[lidx, 'from_bus'])
        to = int(net.line.loc[lidx, 'to_bus'])
        r, x = _get_rx(pos)
        R[(fr, to)] = r; R[(to, fr)] = r
        X[(fr, to)] = x; X[(to, fr)] = x

    # Bus-bus switches
    try:
        for sw in net.switch.loc[(net.switch.et == 'b') & (net.switch.closed)].index:
            fr = int(net.switch.loc[sw, 'bus'])
            to = int(net.switch.loc[sw, 'element'])
            R[(fr, to)] = 0.0; R[(to, fr)] = 0.0
            X[(fr, to)] = 0.0; X[(to, fr)] = 0.0
    except Exception:
        pass

    # Trafos: rows immediately after lines
    n_line = len(net.line)
    for pos, tfidx in enumerate(net.trafo.index):
        if not net.trafo.loc[tfidx, 'in_service']:
            continue
        hv = int(net.trafo.loc[tfidx, 'hv_bus'])
        lv = int(net.trafo.loc[tfidx, 'lv_bus'])
        row = n_line + pos
        r, x = _get_rx(row)
        R[(hv, lv)] = r; R[(lv, hv)] = r
        X[(hv, lv)] = x; X[(lv, hv)] = x

    return R, X


def solve_scp_oltc(net, max_iter=10, eta=1.0, tol=1e-3, V_min=0.90, V_max=1.10,
                   with_oltc=True, verbose=False, damping=0.0, rho=1e-1):
    """Run the SCP MILP loop with discrete OLTC tap variables.

    Returns dict with keys: converged, iterations, obj, runtime_total, runtime_last,
    taps, Vtilde, err_vs_nr, nr_loss_mw, nr_v_range.
    """
    # Pick the "right" system MVA base. For distribution feeders (single low voltage)
    # 1 MVA is reasonable; for transmission (net.sn_mva = 100 usually) use net.sn_mva.
    system_sbase = float(net.sn_mva) * 1e6  # Pa

    result = InitializeConstants(net, Sbase=system_sbase, update=True)
    (Vbase, Ibase, Sbase, Zbase, Buses, SubstationBuses, NonsubstationBuses,
     BidirectionalLines, UnidirectionalLines, UnidirectionalTrafo,
     NonTrafoUnidirectionalLines, TrafoTappos, TrafoTapStepPercent,
     LineTupleToIndex, TrafoTupleToIndex, R, X, AdjacencyList,
     DirectedLinesForEachTree, UnidirectionalLineTupleToIndex) = result

    # Override R, X by reading pandapower's internal ppc branch matrix directly.
    # This delegates all parameter handling (vk, shift_degree, negative reactance,
    # multi-voltage trafo ratio, tap_side) to pandapower.
    R_mv, X_mv = compute_multi_vbase_rx(net, Sbase=system_sbase)
    for key in R:
        if key in R_mv:
            R[key] = R_mv[key]
            X[key] = X_mv[key]

    # Per-bus shunt admittance from Ybus diagonal (captures line charging, trafo
    # magnetizing, MATPOWER tap-side shunt corrections, explicit shunts).
    bus_shunt = compute_bus_shunts(net, Sbase=system_sbase)

    P, Q = RefreshPQ(net, system_sbase)

    Parents, Childs = dict(), dict()
    for bus in Buses:
        Parents[bus], Childs[bus] = [], []
    for line in UnidirectionalLines:
        Parents[line[1]].append(line[0])
        Childs[line[0]].append(line[1])

    _slack_vm0 = float(net.ext_grid.vm_pu.iloc[0])
    Vtilde = {bus: complex(_slack_vm0, 0) for bus in Buses}
    Isgentilde = {bus: 0 for bus in Buses}
    trafo_set = set(UnidirectionalTrafo)

    total_runtime = 0.0
    current_taps = {trafo: 0 for trafo in UnidirectionalTrafo}
    prev_taps = {trafo: None for trafo in UnidirectionalTrafo}
    last_runtime = 0.0
    converged = False
    iteration = 0
    loss_true = float('nan')
    last_Isgenr = {bus: 0.0 for bus in Buses}
    last_Isgeni = {bus: 0.0 for bus in Buses}

    for iteration in range(max_iter):
        Vtilde_prev = copy.deepcopy(Vtilde)
        model = gp.Model(f'SCP_OLTC_iter{iteration}')

        # --- Voltage vars ---
        Vr, Vi = dict(), dict()
        for bus in Buses:
            Vr[bus] = model.addVar(lb=-float('inf'))
            Vi[bus] = model.addVar(lb=-float('inf'))

        # --- Branch currents ---
        Jr, Ji = dict(), dict()
        for line in UnidirectionalLines:
            Jr[line] = model.addVar(lb=-float('inf'))
            Ji[line] = model.addVar(lb=-float('inf'))

        # --- Sgen currents ---
        Isgenr, Isgeni = dict(), dict()
        for bus in Buses:
            Isgenr[bus], Isgeni[bus] = 0, 0
        for sgen in net.sgen.index:
            bus = net.sgen.bus[sgen]
            Isgenr[bus] = model.addVar(lb=-float('inf'))
            Isgeni[bus] = model.addVar(lb=-float('inf'))

        # --- OLTC binaries + J_hv_out auxiliaries ---
        delta = {}
        tap_info = {}
        J_hv_out_r, J_hv_out_i = dict(), dict()
        for trafo in UnidirectionalTrafo:
            tidx = TrafoTupleToIndex[trafo]
            tap_min_val = int(net.trafo.loc[tidx, 'tap_min'])
            tap_max_val = int(net.trafo.loc[tidx, 'tap_max'])
            step = TrafoTapStepPercent[trafo]
            if step is None or np.isnan(step) or step == 0:
                # No tap changer: fix at neutral
                tap_min_val = 0
                tap_max_val = 0
                step = 0
            tap_side = net.trafo.loc[tidx, 'tap_side']
            if tap_side is None or (isinstance(tap_side, float) and np.isnan(tap_side)):
                tap_side = 'hv'  # default
            tap_info[trafo] = (tap_min_val, tap_max_val, step, tap_side)

            if not with_oltc:
                # Force tap=0 (neutral)
                tap_min_val = 0
                tap_max_val = 0
                tap_info[trafo] = (0, 0, step, tap_side)

            is_fixed = (tap_min_val == tap_max_val)
            if not is_fixed:
                for k in range(tap_min_val, tap_max_val + 1):
                    delta[trafo, k] = model.addVar(vtype=GRB.BINARY)
                    if iteration > 0 and current_taps[trafo] == k:
                        delta[trafo, k].Start = 1.0
                    elif iteration > 0:
                        delta[trafo, k].Start = 0.0
                model.addConstr(
                    gp.quicksum(delta[trafo, k] for k in range(tap_min_val, tap_max_val + 1)) == 1
                )
            J_hv_out_r[trafo] = model.addVar(lb=-float('inf'))
            J_hv_out_i[trafo] = model.addVar(lb=-float('inf'))

        # --- Substation voltage fix (read from net.ext_grid.vm_pu) ---
        slack_vm = float(net.ext_grid.vm_pu.iloc[0])
        for bus in SubstationBuses:
            model.addConstr(Vr[bus] == slack_vm)
            model.addConstr(Vi[bus] == 0)

        # --- Linearized current injection ---
        Ir, Ii = dict(), dict()
        for bus in Buses:
            Ir[bus] = model.addVar(lb=-float('inf'))
            Ii[bus] = model.addVar(lb=-float('inf'))

            sq = abs(Vtilde[bus])**2
            orig_Ir = P[bus]*Vtilde[bus].real/sq + Q[bus]*Vtilde[bus].imag/sq
            orig_Ii = P[bus]*Vtilde[bus].imag/sq - Q[bus]*Vtilde[bus].real/sq
            rcoef_vr = P[bus]*(1/sq - 2*Vtilde[bus].real**2/sq**2) - Q[bus]*(2*Vtilde[bus].real*Vtilde[bus].imag/sq**2)
            rcoef_vi = P[bus]*(-2*Vtilde[bus].real*Vtilde[bus].imag/sq**2) + Q[bus]*(1/sq - 2*Vtilde[bus].imag**2/sq**2)
            rresi = Vtilde[bus].real*rcoef_vr + Vtilde[bus].imag*rcoef_vi
            icoef_vr = P[bus]*(2*Vtilde[bus].real*Vtilde[bus].imag/sq**2) + Q[bus]*(1/sq - 2*Vtilde[bus].real**2/sq**2)
            icoef_vi = P[bus]*(2*Vtilde[bus].imag**2/sq**2 - 1/sq) - Q[bus]*(2*Vtilde[bus].real*Vtilde[bus].imag/sq**2)
            iresi = Vtilde[bus].real*icoef_vr + Vtilde[bus].imag*icoef_vi

            # Bus shunt: y_sh * V drawn from bus (acts like extra load → subtracted from Ir,Ii)
            y_sh = bus_shunt.get(int(bus), complex(0, 0))
            g_sh, b_sh = y_sh.real, y_sh.imag

            model.addConstr(Ir[bus] == orig_Ir + eta*(rcoef_vr*Vr[bus] + rcoef_vi*Vi[bus] - rresi) + Isgenr[bus]
                            - g_sh*Vr[bus] + b_sh*Vi[bus])
            model.addConstr(Ii[bus] == orig_Ii + eta*(icoef_vr*Vr[bus] + icoef_vi*Vi[bus] - iresi) + Isgeni[bus]
                            - g_sh*Vi[bus] - b_sh*Vr[bus])

        # --- KCL ---
        for bus in NonsubstationBuses:
            lhs_r, lhs_i = 0, 0
            for parent in Parents[bus]:
                edge = (parent, bus)
                lhs_r = lhs_r + Jr[edge]
                lhs_i = lhs_i + Ji[edge]
            rhs_r, rhs_i = 0, 0
            for child in Childs[bus]:
                edge = (bus, child)
                if edge in trafo_set:
                    rhs_r = rhs_r + J_hv_out_r[edge]
                    rhs_i = rhs_i + J_hv_out_i[edge]
                else:
                    rhs_r = rhs_r + Jr[edge]
                    rhs_i = rhs_i + Ji[edge]
            model.addConstr(lhs_r == -Ir[bus] + rhs_r)
            model.addConstr(lhs_i == -Ii[bus] + rhs_i)

        # --- Ohm for non-trafo lines ---
        for line in NonTrafoUnidirectionalLines:
            model.addConstr(Vr[line[0]] - Vr[line[1]] == R[line]*Jr[line] - X[line]*Ji[line])
            model.addConstr(Vi[line[0]] - Vi[line[1]] == R[line]*Ji[line] + X[line]*Jr[line])

        # --- Ohm + current ratio for trafos ---
        for trafo in UnidirectionalTrafo:
            hv, lv = trafo
            tap_min_val, tap_max_val, step, tap_side = tap_info[trafo]
            is_fixed = (tap_min_val == tap_max_val)  # single tap position → no binaries needed

            for k in range(tap_min_val, tap_max_val + 1):
                t_k = 1 + k * step / 100
                if t_k <= 0:
                    continue
                if tap_side == 'hv':
                    TAP = t_k
                else:
                    TAP = 1 / t_k
                inv_TAP = 1 / TAP

                if is_fixed:
                    # Fixed tap: plain linear constraints (no binary, no indicator)
                    model.addConstr(Vr[hv] - TAP*Vr[lv] - TAP*(R[trafo]*Jr[trafo] - X[trafo]*Ji[trafo]) == 0)
                    model.addConstr(Vi[hv] - TAP*Vi[lv] - TAP*(R[trafo]*Ji[trafo] + X[trafo]*Jr[trafo]) == 0)
                    model.addConstr(J_hv_out_r[trafo] - inv_TAP*Jr[trafo] == 0)
                    model.addConstr(J_hv_out_i[trafo] - inv_TAP*Ji[trafo] == 0)
                else:
                    model.addGenConstrIndicator(
                        delta[trafo, k], 1,
                        Vr[hv] - TAP*Vr[lv] - TAP*(R[trafo]*Jr[trafo] - X[trafo]*Ji[trafo]), GRB.EQUAL, 0.0)
                    model.addGenConstrIndicator(
                        delta[trafo, k], 1,
                        Vi[hv] - TAP*Vi[lv] - TAP*(R[trafo]*Ji[trafo] + X[trafo]*Jr[trafo]), GRB.EQUAL, 0.0)
                    model.addGenConstrIndicator(
                        delta[trafo, k], 1,
                        J_hv_out_r[trafo] - inv_TAP*Jr[trafo], GRB.EQUAL, 0.0)
                    model.addGenConstrIndicator(
                        delta[trafo, k], 1,
                        J_hv_out_i[trafo] - inv_TAP*Ji[trafo], GRB.EQUAL, 0.0)

        # --- Voltage magnitude limits (linearized via Vr only) ---
        for bus in NonsubstationBuses:
            model.addConstr(Vr[bus] >= V_min)
            model.addConstr(Vr[bus] <= V_max)

        # --- Sgen constraints ---
        for sgen in net.sgen.index:
            bus = net.sgen.bus[sgen]
            sn_mva = net.sgen.loc[sgen, 'sn_mva']
            p_mw = net.sgen.loc[sgen, 'p_mw']
            if sn_mva**2 - p_mw**2 <= 0:
                qlim = 0
            else:
                qlim = (sn_mva**2 - p_mw**2)**0.5 * 1e6 / Sbase
            model.addConstr(Vtilde[bus].real*Isgenr[bus] + Vtilde[bus].imag*Isgeni[bus] == 0)
            model.addConstr(Vtilde[bus].imag*Isgenr[bus] - Vtilde[bus].real*Isgeni[bus] <= qlim)
            model.addConstr(Vtilde[bus].imag*Isgenr[bus] - Vtilde[bus].real*Isgeni[bus] >= -qlim)

        # Objective: total series loss + proximal regularization on V
        # The proximal term ρ·|V-V_prev|² prevents oscillation of tap decisions
        # caused by linearization error. ρ is small so original loss is not distorted.
        obj = gp.quicksum(R[line]*(Jr[line]*Jr[line] + Ji[line]*Ji[line]) for line in UnidirectionalLines)
        if iteration > 0 and rho > 0:
            # Sparse proximal: only penalize buses adjacent to OLTC trafos
            # For large networks (1000+ buses) the full-bus quadratic blows up memory.
            prox_buses = set()
            for trafo in UnidirectionalTrafo:
                tap_min_val = tap_info[trafo][0]
                tap_max_val = tap_info[trafo][1]
                if tap_min_val < tap_max_val:  # OLTC-active
                    prox_buses.add(trafo[0])
                    prox_buses.add(trafo[1])
            if not prox_buses:
                prox_buses = set(NonsubstationBuses)  # fallback if no OLTC
            prox = gp.quicksum(
                (Vr[bus] - Vtilde[bus].real)*(Vr[bus] - Vtilde[bus].real) +
                (Vi[bus] - Vtilde[bus].imag)*(Vi[bus] - Vtilde[bus].imag)
                for bus in prox_buses
            )
            obj = obj + rho * prox
        model.setObjective(obj)

        model.setParam('OutputFlag', 1 if verbose else 0)
        model.setParam('MIPGap', 5e-3)  # 0.5% gap — good enough for SCP outer loop
        model.setParam('TimeLimit', 600)  # 10 min per iteration max
        model.optimize()

        if model.status not in [GRB.OPTIMAL, GRB.SUBOPTIMAL, GRB.TIME_LIMIT]:
            if model.SolCount == 0:
                print(f'  iter {iteration}: MIP failed (status={model.status})')
                break
        if model.SolCount == 0:
            print(f'  iter {iteration}: no incumbent (status={model.status})')
            break

        last_runtime = model.Runtime
        total_runtime += model.Runtime

        for bus in Buses:
            V_new = complex(Vr[bus].X, Vi[bus].X)
            if damping > 0:
                Vtilde[bus] = damping * Vtilde[bus] + (1 - damping) * V_new
            else:
                Vtilde[bus] = V_new

        prev_taps_copy = dict(current_taps)
        for trafo in UnidirectionalTrafo:
            tap_min_val, tap_max_val, _, _ = tap_info[trafo]
            if tap_min_val == tap_max_val:
                current_taps[trafo] = tap_min_val  # frozen
                continue
            for k in range(tap_min_val, tap_max_val + 1):
                if delta[trafo, k].X > 0.5:
                    current_taps[trafo] = k
                    break
        taps_unchanged = all(current_taps[t] == prev_taps_copy[t] for t in UnidirectionalTrafo)

        # Grab sgen currents for final validation
        last_Isgenr = {}
        last_Isgeni = {}
        for bus in Buses:
            try:
                last_Isgenr[bus] = Isgenr[bus].X
                last_Isgeni[bus] = Isgeni[bus].X
            except Exception:
                last_Isgenr[bus] = 0
                last_Isgeni[bus] = 0

        # Evaluate true loss (without proximal term)
        loss_true = sum(
            R[line]*(Jr[line].X**2 + Ji[line].X**2) for line in UnidirectionalLines
        )
        v_change = sum(abs(Vtilde[bus] - Vtilde_prev[bus]) for bus in Buses) / len(Buses)
        print(f'  iter {iteration}: loss={loss_true:.6f}, ΔV={v_change:.2e}, taps_same={taps_unchanged}, runtime={model.Runtime:.1f}s')

        if iteration > 0 and v_change < tol and taps_unchanged:
            converged = True
            break

    final_obj = loss_true  # report true loss (without proximal penalty)

    # --- Apply sgen Q and tap positions to pandapower for validation ---
    # Ssg is in per-unit on system_sbase; convert to MVAr for pandapower.
    sbase_mva = system_sbase / 1e6
    for bus in Buses:
        Isgen = complex(last_Isgenr[bus], last_Isgeni[bus])
        Ssg = Vtilde[bus] * Isgen.conjugate()
        Isgentilde[bus] = Isgen
    for idx in net.sgen.index:
        bus = net.sgen.loc[idx, 'bus']
        Isgen = complex(last_Isgenr[bus], last_Isgeni[bus])
        Ssg = Vtilde[bus] * Isgen.conjugate()
        net.sgen.loc[idx, 'q_mvar'] = Ssg.imag * sbase_mva

    for trafo in UnidirectionalTrafo:
        tidx = TrafoTupleToIndex[trafo]
        net.trafo.loc[tidx, 'tap_pos'] = current_taps[trafo]

    # Warm-start continuation for NR validation
    net.trafo.tap_pos = 0
    try:
        pp.runpp(net)
    except Exception:
        pass
    target_taps = {TrafoTupleToIndex[trafo]: current_taps[trafo] for trafo in UnidirectionalTrafo}
    max_abs = max((abs(k) for k in target_taps.values()), default=0)
    nr_converged = True
    for si in range(1, max_abs + 1):
        for tidx, tgt in target_taps.items():
            sign = 1 if tgt > 0 else (-1 if tgt < 0 else 0)
            net.trafo.loc[tidx, 'tap_pos'] = sign * min(si, abs(tgt))
        try:
            pp.runpp(net, init='results')
        except Exception:
            nr_converged = False
            break

    nr_loss = np.nan  # actual branch+trafo loss (MW)
    nr_supply = np.nan  # total ext-grid supply
    err = np.nan
    v_range = (np.nan, np.nan)
    if nr_converged:
        try:
            pp.runpp(net, init='results')
            Veval = {bus: cmath.rect(net.res_bus.vm_pu[bus], math.radians(net.res_bus.va_degree[bus])) for bus in Buses}
            err = sum(abs(Veval[bus] - Vtilde[bus])/abs(Veval[bus]) for bus in Buses) / len(Buses) * 100
            nr_loss = float(net.res_line.pl_mw.sum() + net.res_trafo.pl_mw.sum())
            nr_supply = float(net.res_ext_grid.p_mw.sum())
            v_range = (float(net.res_bus.vm_pu.min()), float(net.res_bus.vm_pu.max()))
        except Exception:
            nr_converged = False

    return {
        'converged': converged,
        'iterations': iteration + 1,
        'obj': final_obj,
        'runtime_total': total_runtime,
        'runtime_last': last_runtime,
        'taps': dict(current_taps),
        'Vtilde': Vtilde,
        'err_vs_nr_pct': err,
        'nr_converged': nr_converged,
        'nr_loss_mw': nr_loss,
        'nr_supply_mw': nr_supply,
        'nr_v_range': v_range,
        'n_buses': len(Buses),
        'n_lines': len(NonTrafoUnidirectionalLines),
        'n_trafos': len(UnidirectionalTrafo),
        'n_sgens': len(net.sgen),
    }


def build_net(name):
    """Return a freshly-built pandapower network for SCP experiments."""
    if name == 'ieee33':
        net = networks.ieee33(vvc=True, use_sgen=True)
        net.ext_grid.vm_pu = 1
        if 'shunt' in net and len(net.shunt) > 0:
            net.shunt['in_service'] = False
        # Enable OLTC: step=0 in factory → 1.5%
        net.trafo['tap_step_percent'] = 1.5
    elif name == 'ieee123':
        net = networks.ieee123(vvc=True, use_sgen=True)
        net.ext_grid.vm_pu = 1
        net.shunt['in_service'] = False
        # Unrealistic vk=0.08% → 6%
        net.trafo['vk_percent'] = 6.0
        net.trafo['vkr_percent'] = 0.6
    elif name == 'ieee118':
        from ieee118 import ieee118
        net = ieee118(mode='vvc', sgen=True)
        net.ext_grid.vm_pu = 1
        if 'shunt' in net and len(net.shunt) > 0:
            net.shunt['in_service'] = False
        net.line.c_nf_per_km = 0
        net.line.g_us_per_km = 0
        # 3 OLTCs (the rest frozen at tap=0) — MIQP with all 13 becomes dominated by
        # linearization noise when multiple cross-voltage taps are moved simultaneously.
        keep = list(net.trafo.index)[:3]
        for idx in net.trafo.index:
            if idx in keep:
                continue
            net.trafo.loc[idx, 'tap_pos'] = 0
            net.trafo.loc[idx, 'tap_min'] = 0
            net.trafo.loc[idx, 'tap_max'] = 0
            net.trafo.loc[idx, 'tap_step_percent'] = 0
    elif name == 'ieee14':
        from ieee14 import ieee14
        net = ieee14(mode='vvc', sgen=True)
        net.ext_grid.vm_pu = 1
        if 'shunt' in net and len(net.shunt) > 0:
            net.shunt['in_service'] = False
        net.line.c_nf_per_km = 0
        net.line.g_us_per_km = 0
        # Keep only 1 OLTC (first trafo); freeze the rest at tap=0 with fixed range
        keep = list(net.trafo.index)[:1]
        for idx in net.trafo.index:
            if idx in keep:
                continue
            net.trafo.loc[idx, 'tap_pos'] = 0
            net.trafo.loc[idx, 'tap_min'] = 0
            net.trafo.loc[idx, 'tap_max'] = 0
            net.trafo.loc[idx, 'tap_step_percent'] = 0
    elif name == 'ieee8500':
        net = networks.ieee8500(vvc=True, use_sgen=True)
        net.ext_grid.vm_pu = 1
        if 'shunt' in net and len(net.shunt) > 0:
            net.shunt['in_service'] = False
        net.line.c_nf_per_km = 0
        net.line.g_us_per_km = 0
        # All 20 trafos as OLTC
    elif name == 'ieee1888':
        from ieee1888 import ieee1888
        net = ieee1888(mode='vvc', sgen=True)
        net.ext_grid.vm_pu = 1
        if 'shunt' in net and len(net.shunt) > 0:
            net.shunt['in_service'] = False
        net.line.c_nf_per_km = 0
        net.line.g_us_per_km = 0
        # 555 trafos — keep 10 as OLTC, freeze the rest
        keep = list(net.trafo.index)[:10]
        for idx in net.trafo.index:
            if idx in keep:
                continue
            net.trafo.loc[idx, 'tap_pos'] = 0
            net.trafo.loc[idx, 'tap_min'] = 0
            net.trafo.loc[idx, 'tap_max'] = 0
            net.trafo.loc[idx, 'tap_step_percent'] = 0
    else:
        raise ValueError(f'Unknown network name: {name}')

    # Normalize tap parameters for all trafos
    # case118/case14 have NaN tap_min/tap_max/tap_step → assign realistic defaults
    net.trafo['tap_pos'] = net.trafo['tap_pos'].fillna(0).astype(int)
    net.trafo['tap_step_percent'] = net.trafo['tap_step_percent'].fillna(1.5)
    net.trafo['tap_side'] = net.trafo['tap_side'].fillna('hv')
    # Force a ±4 tap range where NaN — narrower than ±9 to keep binary count manageable
    net.trafo['tap_min'] = net.trafo['tap_min'].fillna(-4).astype(int)
    net.trafo['tap_max'] = net.trafo['tap_max'].fillna(4).astype(int)
    net.trafo['tap_neutral'] = net.trafo['tap_neutral'].fillna(0).astype(int) if 'tap_neutral' in net.trafo.columns else 0
    # Clamp existing ranges to ±4 to limit MIQP binary count (9 positions instead of 19)
    net.trafo['tap_min'] = net.trafo['tap_min'].clip(lower=-4)
    net.trafo['tap_max'] = net.trafo['tap_max'].clip(upper=4)

    net.ext_grid.va_degree = 0
    return net


def run_case(name):
    print(f'\n======== {name} ========')
    # With OLTC
    net = build_net(name)
    try:
        pp.runpp(net)
        base_loss = float(net.res_ext_grid.p_mw.sum())
    except Exception:
        base_loss = np.nan
    print(f'Base: buses={len(net.bus)}, lines={len(net.line)}, trafos={len(net.trafo)}, sgens={len(net.sgen)}')
    print(f'Base NR loss (tap=0): {base_loss:.4f} MW')

    print('--- Solving WITH OLTC ---')
    r_on = solve_scp_oltc(net, max_iter=10, with_oltc=True)
    r_on['base_loss_mw'] = base_loss

    # Without OLTC (tap forced to 0) — fresh net
    net2 = build_net(name)
    print('--- Solving WITHOUT OLTC (tap fixed at 0) ---')
    r_off = solve_scp_oltc(net2, max_iter=10, with_oltc=False)

    print(f'\nResults for {name}:')
    print(f'  converged: ON={r_on["converged"]}, OFF={r_off["converged"]}')
    print(f'  obj (p.u.):      ON={r_on["obj"]:.6f}, OFF={r_off["obj"]:.6f}')
    print(f'  iterations:      ON={r_on["iterations"]}, OFF={r_off["iterations"]}')
    print(f'  total runtime:   ON={r_on["runtime_total"]:.1f}s, OFF={r_off["runtime_total"]:.1f}s')
    print(f'  NR loss (MW):    ON={r_on["nr_loss_mw"]:.4f}, OFF={r_off["nr_loss_mw"]:.4f}')
    print(f'  ΔV vs NR (%):    ON={r_on["err_vs_nr_pct"]:.2e}, OFF={r_off["err_vs_nr_pct"]:.2e}')
    print(f'  V range (NR):    ON={r_on["nr_v_range"]}, OFF={r_off["nr_v_range"]}')
    print(f'  Taps found:      {list(r_on["taps"].values())}')
    if r_on["nr_converged"] and r_off["nr_converged"]:
        improvement = (r_off['nr_loss_mw'] - r_on['nr_loss_mw']) / r_off['nr_loss_mw'] * 100
        print(f'  Improvement:     {improvement:.3f} %')

    return {'name': name, 'on': r_on, 'off': r_off}


#%%
if __name__ == '__main__':
    import sys
    cases = sys.argv[1:] if len(sys.argv) > 1 else ['ieee33', 'ieee14', 'ieee123', 'ieee118']
    summary = {}
    for case in cases:
        try:
            summary[case] = run_case(case)
        except Exception as e:
            print(f'\n{case} FAILED: {type(e).__name__}: {e}')
            import traceback; traceback.print_exc()
            summary[case] = {'error': str(e)}

    print('\n\n======== SUMMARY TABLE ========')
    print(f'{"Network":<12}{"buses":<8}{"trafos":<8}{"iter":<6}{"obj_ON":<12}{"obj_OFF":<12}{"NR_ON(MW)":<12}{"NR_OFF(MW)":<12}{"impr%":<8}{"ΔV%":<10}{"rt(s)":<8}')
    for case in cases:
        s = summary.get(case, {})
        if 'on' in s:
            on, off = s['on'], s['off']
            improvement = '-'
            if on['nr_converged'] and off['nr_converged']:
                improvement = f"{(off['nr_loss_mw']-on['nr_loss_mw'])/off['nr_loss_mw']*100:.3f}"
            print(f"{case:<12}{on['n_buses']:<8}{on['n_trafos']:<8}{on['iterations']:<6}"
                  f"{on['obj']:<12.6f}{off['obj']:<12.6f}"
                  f"{on['nr_loss_mw']:<12.4f}{off['nr_loss_mw']:<12.4f}"
                  f"{improvement:<8}{on['err_vs_nr_pct']:<10.2e}{on['runtime_total']:<8.1f}")
        else:
            print(f"{case:<12}FAILED")
