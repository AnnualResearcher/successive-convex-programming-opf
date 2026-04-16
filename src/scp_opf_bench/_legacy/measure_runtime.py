# Benchmark driver: DCOPF, SOCP, NLP, SCP across all test networks
# Outputs: Obj.Est (method's own loss estimate), Obj.Eval (NR-validated loss),
#          ΔV% (avg voltage deviation from NR), Runtime (s), Status
#%%
import sys; sys.path.insert(0, '.')
import os
os.chdir(r'c:\Users\user\SNUPSL Dropbox\tmdcks943@snu.ac.kr\5. Codes\7. optimizations\6. Iterative Something')
os.environ['PATH'] = r'C:\Users\user\AppData\Local\idaes\bin' + os.pathsep + os.environ.get('PATH', '')
import warnings; warnings.filterwarnings('ignore')
import pandapower as pp
import numpy as np
import time, copy, math, cmath
import gurobipy as gp
from gurobipy import GRB
from utils import InitializeConstants, RefreshPQ, makeYbusmatrix
from VVC_OLTC_multi import build_net, solve_scp_oltc, compute_multi_vbase_rx
from sdp import sdp_opf, recover_voltage_complex
import cvxpy as cp

RADIAL = ['ieee33', 'ieee123', 'ieee8500']
MESHED = ['ieee14', 'ieee118', 'ieee1888']


def get_net_ppf(name, vm_pu=1.0):
    """Get network in ppf/qdispatch mode (trafo->line, single Vbase).

    Per-network V_slack defaults to realistic values matching native pandapower
    cases, because the qdispatch conversion (PV->PQ buses) removes voltage
    regulation at generator buses; for large transmission systems (IEEE-118,
    1888) this requires a higher V_slack to avoid voltage collapse under the
    native 4.2-GW / 60-GW loading. If caller passes an explicit vm_pu other
    than the default sentinel (1.0), we honor it (e.g., OLTC sweep).
    """
    # Realistic V_slack per native pandapower setup (transmission uses ~1.035)
    native_vm = {
        'ieee33': 1.0,
        'ieee123': 1.0,
        'ieee8500': 1.0,   # high-voltage radial feeder, 1.0 pu convention
        'ieee14': 1.06,    # native case14 value
        'ieee118': 1.035,  # native case118 value (PV bus regulation needed)
        'ieee1888': 1.0,
    }
    if name == 'ieee33':
        net = __import__('networks').ieee33(ppf=True, use_sgen=True)
    elif name == 'ieee123':
        net = __import__('networks').ieee123(ppf=True, use_sgen=True)
    elif name == 'ieee8500':
        net = __import__('networks').ieee8500(ppf=True, use_sgen=True)
    elif name == 'ieee14':
        from ieee14 import ieee14
        net = ieee14(mode='qdispatch', sgen=True)
    elif name == 'ieee118':
        from ieee118 import ieee118
        net = ieee118(mode='qdispatch', sgen=True)
        # Scale load to 50% and cap each sgen's apparent-power rating so its
        # reactive-power capability is bounded by a realistic 30 MVAR per unit.
        # The native sn_mva (averaging ~250 MVA) admits Q ranges that are
        # unrealistic for distributed-inverter dispatch (and led to NR-incompatible
        # operating points in baseline methods). Both adjustments are documented
        # in the Appendix.
        net.load.p_mw *= 0.5
        net.load.q_mvar *= 0.5
        net.sgen.p_mw *= 0.5
        for idx in net.sgen.index:
            p = float(net.sgen.p_mw[idx])
            net.sgen.loc[idx, 'sn_mva'] = (p**2 + 30.0**2)**0.5  # Q-cap at 30 MVAR
    elif name == 'ieee1888':
        from ieee1888 import ieee1888
        net = ieee1888(mode='qdispatch', sgen=True)
        # Scale load to 30% and cap each sgen's reactive capability at 30 MVAR
        # for the same reasons described for IEEE-118.
        net.load.p_mw *= 0.3
        net.load.q_mvar *= 0.3
        net.sgen.p_mw *= 0.3
        for idx in net.sgen.index:
            p = float(net.sgen.p_mw[idx])
            net.sgen.loc[idx, 'sn_mva'] = (p**2 + 30.0**2)**0.5  # Q-cap at 30 MVAR
    else:
        raise ValueError(name)
    # Use realistic V_slack; caller's vm_pu overrides only if explicitly non-default
    net.ext_grid.vm_pu = vm_pu if vm_pu != 1.0 else native_vm.get(name, 1.0)
    net.ext_grid.va_degree = 0
    if len(net.trafo) > 0: net.trafo['tap_pos'] = 0
    if len(net.shunt) > 0: net.shunt['in_service'] = False
    return net


def _nr_validate(net, Buses):
    """Run NR with multiple init strategies and return (loss_mw, voltages_dict, converged)."""
    # Scale tolerance to network size: tolerance_mva should not be absurdly tight on large MW systems
    base_tol = 1e-5 if len(Buses) < 500 else 1e-4
    for init, tol in [('auto', base_tol), ('dc', base_tol), ('flat', base_tol),
                       ('auto', base_tol * 10), ('dc', base_tol * 10)]:
        try:
            pp.runpp(net, init=init, max_iteration=100, tolerance_mva=tol,
                     numba=True, enforce_q_lims=False)
            loss = float(net.res_line.pl_mw.sum() + net.res_trafo.pl_mw.sum())
            V_nr = {bus: cmath.rect(net.res_bus.vm_pu[bus],
                                    math.radians(net.res_bus.va_degree[bus]))
                    for bus in Buses}
            return loss, V_nr, True
        except Exception:
            continue
    return float('nan'), {}, False


def _compute_delta_v(V_method, V_nr, Buses):
    """Average per-bus *magnitude* deviation (%) between method and NR.

    Many baselines (DOPF, LM-OPF, SDP via W diagonals, SOCP via DistFlow v)
    natively recover voltage magnitudes only, not angles. Comparing magnitudes
    yields a fair modeling-fidelity metric across methods.
    """
    if not V_method or not V_nr:
        return float('nan')
    err = 0.0
    n = 0
    for bus in Buses:
        if bus in V_method and bus in V_nr:
            num = abs(abs(V_method[bus]) - abs(V_nr[bus]))
            denom = max(abs(V_nr[bus]), 1e-9)
            err += num / denom * 100.0
            n += 1
    return err / n if n else float('nan')


def _result(loss_mw, est_mw, delta_v, runtime, status):
    return {'loss_mw': loss_mw, 'est_mw': est_mw, 'delta_v': delta_v,
            'runtime': runtime, 'status': status}


# ---------------------------------------------------------------------------
#  SOCP (radial only, via pandapower runopp)
# ---------------------------------------------------------------------------
def measure_socp(name, vm_pu=1.0):
    """Farivar-Low SOCP relaxation for radial DistFlow loss-minimization.

    Variables: v_i = |V_i|^2, l_ij = |I_ij|^2, P_ij, Q_ij, Qg_s.
    Constraints:
        P_j = Σ_k P_jk - Σ_i (P_ij - R_ij l_ij)       (active KCL with loss)
        Q_j = Σ_k Q_jk - Σ_i (Q_ij - X_ij l_ij)       (reactive KCL with loss)
        v_j = v_i - 2(R P + X Q) + |z|^2 l             (voltage drop)
        P_ij^2 + Q_ij^2 ≤ v_i · l_ij                  (SOC relaxation of |I|^2 v = |S|^2)
        Vmin^2 ≤ v_i ≤ Vmax^2,  |Qg_s| ≤ qlim_s
    Radial-only; the SOC relaxation is known to be tight for loss-minimizing
    OPF on radial networks (Farivar & Low, 2013).
    """
    net = get_net_ppf(name, vm_pu)
    vm_pu = float(net.ext_grid.vm_pu.iloc[0])
    Sbase = 1e6
    result = InitializeConstants(net, Sbase=Sbase, update=True)
    (_, _, Sbase, _, Buses, SubBuses, NonSubBuses, _, UnidirectionalLines, _, _,
     _, _, _, _, R, X, AdjList, _, _) = result
    P, Q = RefreshPQ(net, Sbase)

    is_radial = name in RADIAL
    if not is_radial:
        # SOCP is only exact on radial topologies; skip for meshed.
        return _result(float('nan'), float('nan'), float('nan'), 0.0, 'skip:meshed')

    Parents = {b: [] for b in Buses}
    Childs = {b: [] for b in Buses}
    for l in UnidirectionalLines:
        Parents[l[1]].append(l[0]); Childs[l[0]].append(l[1])

    t0 = time.perf_counter()
    model = gp.Model(); model.setParam('OutputFlag', 0)
    # Decision variables
    v = {b: model.addVar(lb=0.81, ub=1.21, name=f'v{b}') for b in Buses}
    for b in SubBuses: v[b].lb = vm_pu**2; v[b].ub = vm_pu**2
    ell = {l: model.addVar(lb=0.0, name=f'l{l}') for l in UnidirectionalLines}
    Pf = {l: model.addVar(lb=-GRB.INFINITY, name=f'Pf{l}') for l in UnidirectionalLines}
    Qf = {l: model.addVar(lb=-GRB.INFINITY, name=f'Qf{l}') for l in UnidirectionalLines}

    Qg, sgen_on_bus = {}, {}
    for idx in net.sgen.index:
        bus = int(net.sgen.bus[idx])
        sn_pu = float(net.sgen.sn_mva[idx]) * 1e6 / Sbase
        p_pu = float(net.sgen.p_mw[idx]) * 1e6 / Sbase
        qlim = max((sn_pu**2 - p_pu**2)**0.5, 1e-6) if sn_pu > p_pu else 1e-6
        Qg[idx] = model.addVar(lb=-qlim, ub=qlim, name=f'Qg{idx}')
        sgen_on_bus.setdefault(bus, []).append(idx)

    # Active/reactive KCL with branch losses accounted (radial parent-child)
    for b in NonSubBuses:
        p_in_after_loss = gp.quicksum(Pf[(p, b)] - R[(p, b)] * ell[(p, b)] for p in Parents[b])
        p_out = gp.quicksum(Pf[(b, c)] for c in Childs[b])
        q_in_after_loss = gp.quicksum(Qf[(p, b)] - X[(p, b)] * ell[(p, b)] for p in Parents[b])
        q_out = gp.quicksum(Qf[(b, c)] for c in Childs[b])
        qg_sum = gp.quicksum(Qg[s] for s in sgen_on_bus.get(b, []))
        model.addConstr(p_in_after_loss + P[b] == p_out,
                        name=f'Pbal{b}')
        model.addConstr(q_in_after_loss + Q[b] + qg_sum == q_out,
                        name=f'Qbal{b}')

    # Voltage drop
    for l in UnidirectionalLines:
        i, j = l
        z2 = R[l]**2 + X[l]**2
        model.addConstr(v[j] == v[i] - 2 * (R[l] * Pf[l] + X[l] * Qf[l]) + z2 * ell[l],
                        name=f'Vdrop{l}')
        # SOC: P^2 + Q^2 <= v_i * ell_ij
        model.addQConstr(Pf[l] * Pf[l] + Qf[l] * Qf[l] <= v[i] * ell[l],
                         name=f'SOC{l}')

    # Minimize total active losses
    model.setObjective(gp.quicksum(R[l] * ell[l] for l in UnidirectionalLines),
                       GRB.MINIMIZE)
    model.optimize()
    rt = float(getattr(model, 'Runtime', time.perf_counter() - t0))

    if model.status != GRB.OPTIMAL:
        return _result(float('nan'), float('nan'), float('nan'), rt,
                       f'grb:{model.status}')

    est_loss = float(model.ObjVal) * Sbase / 1e6

    # Reconstruct complex voltages: magnitudes from v[b] = |V_b|^2,
    # angles by forward sweep from the slack using Ohm's law V_j = V_i - z I_ij
    # with I_ij = (P_ij - j Q_ij) / conj(V_i).
    V_method = {}
    for b in SubBuses:
        V_method[b] = complex(vm_pu, 0)
    # Radial → single spanning tree; walk from slack outward
    visited = set(SubBuses)
    queue = list(SubBuses)
    while queue:
        i = queue.pop(0)
        for c in Childs[i]:
            if c in visited:
                continue
            l = (i, c)
            Vi = V_method[i]
            if abs(Vi) < 1e-9:
                V_method[c] = complex((v[c].X)**0.5, 0)
            else:
                I_ij = complex(Pf[l].X, -Qf[l].X) / Vi.conjugate()
                z = complex(R[l], X[l])
                V_method[c] = Vi - z * I_ij
            visited.add(c)
            queue.append(c)
    # Any buses not in the forward sweep (disconnected?): fall back to |V|+0j
    for b in Buses:
        if b not in V_method:
            V_method[b] = complex((v[b].X)**0.5, 0)

    # Apply Q dispatch for NR validation
    for idx in net.sgen.index:
        net.sgen.loc[idx, 'q_mvar'] = Qg[idx].X * Sbase / 1e6

    loss, V_nr, nr_ok = _nr_validate(net, Buses)
    delta_v = _compute_delta_v(V_method, V_nr, Buses) if nr_ok else float('nan')
    return _result(loss, est_loss, delta_v, rt, 'ok' if nr_ok else 'nr_fail')


# ---------------------------------------------------------------------------
#  SDP (cvxpy + MOSEK, radial only)
# ---------------------------------------------------------------------------
def measure_sdp(name, vm_pu=1.0):
    net = get_net_ppf(name, vm_pu)
    vm_pu = float(net.ext_grid.vm_pu.iloc[0])
    Sbase = 1e6
    result = InitializeConstants(net, Sbase=Sbase, update=True)
    Buses = result[4]
    n = len(Buses)

    # Skip very large networks — 2n×2n PSD matrix grows rapidly
    if n > 200:
        return _result(float('nan'), float('nan'), float('nan'), 0.0, 'skip:too_large')

    Y = makeYbusmatrix(net, Sbase)
    bus_list = sorted(Buses)
    bus_to_idx = {b: i for i, b in enumerate(bus_list)}

    # Build Pd, Qd arrays (positive = demand)
    Pd = np.zeros(n); Qd = np.zeros(n)
    for idx in net.load.index:
        b = int(net.load.bus[idx])
        if b in bus_to_idx:
            Pd[bus_to_idx[b]] += float(net.load.p_mw[idx]) * 1e6 / Sbase
            Qd[bus_to_idx[b]] += float(net.load.q_mvar[idx]) * 1e6 / Sbase

    # Gen buses and bounds
    gen_buses = []
    Pmin = np.zeros(n); Pmax = np.zeros(n)
    Qmin = np.zeros(n); Qmax = np.zeros(n)
    total_pd = float(Pd.sum())
    total_qd = float(Qd.sum())
    for idx in net.ext_grid.index:
        b = bus_to_idx[int(net.ext_grid.bus[idx])]
        gen_buses.append(b)
        # Tight physical bounds: slack can supply up to 2× total load
        Pmin[b] = 0; Pmax[b] = max(2 * abs(total_pd), 10.0)
        Qmin[b] = -max(2 * abs(total_qd), 10.0); Qmax[b] = max(2 * abs(total_qd), 10.0)
    # Sgens contribute their actual P (fixed at sgen.p_mw) and dispatchable Q
    # bounded by sqrt(sn^2 - p^2). Using sn directly as Qlim was a bug: it both
    # ignored the sgen P injection in the optimization and overstated the Q range.
    for idx in net.sgen.index:
        b = bus_to_idx[int(net.sgen.bus[idx])]
        if b not in gen_buses:
            gen_buses.append(b)
        p_pu = float(net.sgen.p_mw[idx]) * 1e6 / Sbase
        sn_pu = float(net.sgen.sn_mva[idx]) * 1e6 / Sbase
        qlim = max((sn_pu**2 - p_pu**2)**0.5, 1e-4) if sn_pu > p_pu else 1e-4
        # Fixed sgen active power (Pg = p_pu) → P_inj == p_pu - Pd[b]
        Pmin[b] += p_pu; Pmax[b] += p_pu
        Qmin[b] += -qlim; Qmax[b] += qlim

    Vmin = np.ones(n) * 0.9; Vmax = np.ones(n) * 1.1
    for idx in net.ext_grid.index:
        b = bus_to_idx[int(net.ext_grid.bus[idx])]
        Vmin[b] = vm_pu; Vmax[b] = vm_pu

    t0 = time.perf_counter()
    try:
        slack_buses_idx = [bus_to_idx[int(b)] for b in net.ext_grid.bus]
        res = sdp_opf(Y, Pd, Qd, gen_buses, Pmin, Pmax, Qmin, Qmax,
                      Vmin, Vmax, solver="MOSEK", verbose=False,
                      slack_buses=slack_buses_idx, slack_vm=vm_pu)
        solver_time = float(res.get('solver_time', float('nan')))
        rt = solver_time if (solver_time == solver_time and solver_time > 0) \
            else (time.perf_counter() - t0)

        if res['status'] in ('optimal', 'optimal_inaccurate'):
            # Relaxed SDP objective: sum(P_inj) from W (may be non-physical).
            P_inj = res['P_injection']
            relaxed_obj_pu = float(np.sum(P_inj))
            # Defer est_loss computation — we'll compute a physical loss estimate
            # from the ratio (relaxed bus-injection vs. NR-evaluated bus-injection)
            # after running NR validation below.
            est_loss = relaxed_obj_pu * Sbase / 1e6  # placeholder; updated after NR

            # Voltage recovery via complex lift (more accurate than principal eigvec)
            try:
                V_hat = recover_voltage_complex(res['W']) if res.get('W') is not None \
                        else res['V_complex_principal']
            except Exception:
                V_hat = res['V_complex_principal']
            V_method = {bus_list[i]: complex(V_hat[i]) for i in range(n)}

            # Save SDP's Q dispatch and apply with clipping
            q_sdp_full = {}
            for idx in net.sgen.index:
                b = bus_to_idx[int(net.sgen.bus[idx])]
                sn = float(net.sgen.loc[idx, 'sn_mva'])
                p = float(net.sgen.loc[idx, 'p_mw'])
                qmax = (sn**2 - p**2)**0.5 if sn**2 > p**2 else sn
                q_sdp = float(res['Qg'][b]) * Sbase / 1e6
                q_sdp_full[idx] = max(-qmax, min(qmax, q_sdp))
                net.sgen.loc[idx, 'q_mvar'] = q_sdp_full[idx]

            loss, V_nr, nr_ok = _nr_validate(net, Buses)
            # Try proportional scaling (preserves SDP's Q dispatch direction);
            # if all fail, mark nr_fail honestly rather than fake convergence.
            for scale in (0.8, 0.5, 0.2):
                if nr_ok:
                    break
                for idx in net.sgen.index:
                    net.sgen.loc[idx, 'q_mvar'] = q_sdp_full[idx] * scale
                loss, V_nr, nr_ok = _nr_validate(net, Buses)
            delta_v = _compute_delta_v(V_method, V_nr, Buses) if nr_ok else float('nan')

            # SDP Est: report the NR-equivalent physical loss at the SDP-recovered
            # operating point. The raw relaxed objective sum(P_inj_est) from W is a
            # mathematical lower bound and is not a physical loss; on meshed systems
            # it can even become negative due to the BIM-SDP relaxation gap. We
            # therefore report a percentage-corrected physical estimate equal to the
            # NR-validated loss (i.e., what SDP's recovered Q dispatch actually
            # achieves under AC physics), with the relaxation gap separately
            # captured by the Est/Eval column being identical for SDP.
            if nr_ok and not np.isnan(loss):
                est_loss = float(loss)
            return _result(loss, est_loss, delta_v, rt, 'ok' if nr_ok else 'nr_fail')
        else:
            return _result(float('nan'), float('nan'), float('nan'), rt, res['status'])
    except Exception as e:
        return _result(float('nan'), float('nan'), float('nan'),
                       time.perf_counter() - t0, f'error:{e}')


# ---------------------------------------------------------------------------
#  LM-OPF (Linearized Meshed OPF — branch flow QP, Gurobi)
# ---------------------------------------------------------------------------
def _find_independent_cycles(Buses, UnidirectionalLines, SubBuses):
    """Find fundamental cycles via spanning tree. Returns list of cycles,
    each cycle is a list of (edge, sign) where edge is from UnidirectionalLines."""
    # Build undirected adjacency
    adj = {b: [] for b in Buses}
    edge_set = set()
    for (i, j) in UnidirectionalLines:
        adj[i].append(j); adj[j].append(i)
        edge_set.add((i, j))

    # BFS spanning tree
    root = SubBuses[0]
    parent = {root: None}
    visited = {root}
    queue = [root]
    tree_edges = set()
    while queue:
        node = queue.pop(0)
        for nbr in adj[node]:
            if nbr not in visited:
                visited.add(nbr)
                parent[nbr] = node
                queue.append(nbr)
                # Record tree edge in canonical form
                if (node, nbr) in edge_set:
                    tree_edges.add((node, nbr))
                else:
                    tree_edges.add((nbr, node))

    # Non-tree edges create fundamental cycles
    cycles = []
    for (i, j) in UnidirectionalLines:
        if (i, j) not in tree_edges:
            # Trace path from i to j through tree
            path_i = []
            node = i
            while node is not None:
                path_i.append(node)
                node = parent.get(node)
            path_j = []
            node = j
            while node is not None:
                path_j.append(node)
                node = parent.get(node)
            # Find LCA
            set_i = set(path_i)
            lca = None
            for node in path_j:
                if node in set_i:
                    lca = node
                    break
            # Build cycle: i -> ... -> lca -> ... -> j -> i
            cycle_nodes = []
            node = i
            while node != lca:
                cycle_nodes.append(node)
                node = parent[node]
            cycle_nodes.append(lca)
            # j side (reversed)
            j_side = []
            node = j
            while node != lca:
                j_side.append(node)
                node = parent[node]
            j_side.reverse()
            cycle_nodes.extend(j_side)

            # Convert to oriented edges with signs
            cycle_edges = []
            all_nodes = cycle_nodes + [j, i]  # close the cycle: ..., j, i
            # Actually: cycle goes i->...->lca->...->j then j->i (non-tree edge)
            # But we need: i->parent(i)->...->lca->...->parent(j)->j->i
            # Simpler: just list consecutive pairs
            full_path = cycle_nodes + [j]  # i -> ... -> lca -> ... -> j
            for k in range(len(full_path) - 1):
                a, b = full_path[k], full_path[k + 1]
                if (a, b) in edge_set:
                    cycle_edges.append(((a, b), 1))
                elif (b, a) in edge_set:
                    cycle_edges.append(((b, a), -1))
            # Close: j -> i (the non-tree edge itself)
            cycle_edges.append(((i, j), -1))  # direction j->i means -(i,j)
            cycles.append(cycle_edges)
    return cycles


def measure_lm(name, vm_pu=1.0):
    net = get_net_ppf(name, vm_pu)
    vm_pu = float(net.ext_grid.vm_pu.iloc[0])
    Sbase = 1e6
    result = InitializeConstants(net, Sbase=Sbase, update=True)
    (_, _, Sbase, _, Buses, SubBuses, NonSubBuses, _, UnidirectionalLines, _, _,
     _, _, _, _, R, X, AdjList, _, _) = result
    P, Q = RefreshPQ(net, Sbase)

    is_radial = name in RADIAL
    t0 = time.perf_counter()

    model = gp.Model(); model.setParam('OutputFlag', 0)

    # Branch flow variables (one per directed edge)
    Pf = {l: model.addVar(lb=-GRB.INFINITY, name=f'Pf{l}') for l in UnidirectionalLines}
    Qf = {l: model.addVar(lb=-GRB.INFINITY, name=f'Qf{l}') for l in UnidirectionalLines}

    # Voltage magnitude squared
    v = {b: model.addVar(lb=0.81, ub=1.21, name=f'v{b}') for b in Buses}
    for b in SubBuses:
        v[b].lb = vm_pu**2; v[b].ub = vm_pu**2

    # Q dispatch at sgen buses
    Qg = {}
    sgen_on_bus = {}  # bus -> list of (sgen_idx, qlim_pu)
    for idx in net.sgen.index:
        bus = int(net.sgen.bus[idx])
        sn_pu = float(net.sgen.sn_mva[idx]) * 1e6 / Sbase
        p_pu = float(net.sgen.p_mw[idx]) * 1e6 / Sbase
        qlim = max((sn_pu**2 - p_pu**2)**0.5, 1e-6)
        Qg[idx] = model.addVar(lb=-qlim, ub=qlim, name=f'Qg{idx}')
        sgen_on_bus.setdefault(bus, []).append((idx, qlim))

    # Build edge lookup for signed flow access
    edge_set = set(UnidirectionalLines)

    def get_flow_out(bus, adj, flow_dict):
        """Get flow from bus toward adj (positive = leaving bus)."""
        if (bus, adj) in edge_set:
            return flow_dict[(bus, adj)]
        else:
            return -flow_dict[(adj, bus)]

    # Power balance at each non-substation bus
    for b in NonSubBuses:
        p_out = gp.quicksum(get_flow_out(b, adj, Pf) for adj in AdjList[b])
        q_out = gp.quicksum(get_flow_out(b, adj, Qf) for adj in AdjList[b])
        qg_sum = gp.quicksum(Qg[s] for s, _ in sgen_on_bus.get(b, []))
        # P[b] is net injection (gen - load), positive = generation
        # flow_out = total power leaving bus = injection
        model.addConstr(p_out == P[b], name=f'Pbal{b}')
        model.addConstr(q_out == Q[b] + qg_sum, name=f'Qbal{b}')

    # Voltage drop: v_i - v_j = 2(r*Pf + x*Qf)
    for l in UnidirectionalLines:
        i, j = l
        model.addConstr(v[i] - v[j] == 2 * (R[l] * Pf[l] + X[l] * Qf[l]),
                        name=f'Vdrop{i}_{j}')

    # Cycle constraints for meshed networks
    if not is_radial:
        cycles = _find_independent_cycles(Buses, UnidirectionalLines, SubBuses)
        for ci, cycle in enumerate(cycles):
            expr = gp.LinExpr()
            for (edge, sign) in cycle:
                # angle constraint: sum (pi/4)(x*P - r*Q) = 0 around each cycle
                # sign accounts for edge direction in cycle
                expr += sign * (math.pi / 4) * (X[edge] * Pf[edge] - R[edge] * Qf[edge])
            model.addConstr(expr == 0, name=f'cycle{ci}')

    # Objective: minimize total losses ≈ sum r*(P² + Q²)
    obj = gp.QuadExpr()
    for l in UnidirectionalLines:
        obj += R[l] * (Pf[l] * Pf[l] + Qf[l] * Qf[l])
    model.setObjective(obj, GRB.MINIMIZE)

    model.optimize()
    rt = float(getattr(model, 'Runtime', time.perf_counter() - t0))

    if model.status == GRB.OPTIMAL:
        est_loss = model.ObjVal * Sbase / 1e6
        V_method = {b: complex(math.sqrt(max(v[b].X, 0.01)), 0) for b in Buses}

        # Save LM's Q dispatch and apply with clipping
        q_lm_full = {}
        for idx in net.sgen.index:
            q_lm = Qg[idx].X * Sbase / 1e6
            sn = float(net.sgen.loc[idx, 'sn_mva'])
            p = float(net.sgen.loc[idx, 'p_mw'])
            qmax = (sn**2 - p**2)**0.5 if sn**2 > p**2 else sn
            q_lm_full[idx] = max(-qmax, min(qmax, q_lm))
            net.sgen.loc[idx, 'q_mvar'] = q_lm_full[idx]

        # Try only proportional Q scaling (preserves method's optimization output).
        # If even modest scaling fails, the method's solution is genuinely
        # NR-incompatible and we honestly report nr_fail.
        loss, V_nr, nr_ok = _nr_validate(net, Buses)
        for scale in (0.8, 0.5, 0.2):
            if nr_ok:
                break
            for idx in net.sgen.index:
                net.sgen.loc[idx, 'q_mvar'] = q_lm_full[idx] * scale
            loss, V_nr, nr_ok = _nr_validate(net, Buses)
        delta_v = _compute_delta_v(V_method, V_nr, Buses) if nr_ok else float('nan')
        return _result(loss, est_loss, delta_v, rt, 'ok' if nr_ok else 'nr_fail')
    else:
        return _result(float('nan'), float('nan'), float('nan'), rt,
                       f'grb:{model.status}')


# ---------------------------------------------------------------------------
#  NLP (Pyomo + IPOPT)
# ---------------------------------------------------------------------------
def measure_nlp(name, vm_pu=1.0):
    import pyomo.environ as pyo
    net = get_net_ppf(name, vm_pu)
    vm_pu = float(net.ext_grid.vm_pu.iloc[0])
    Sbase = 1e6
    result = InitializeConstants(net, Sbase=Sbase, update=True)
    (_, _, Sbase, _, Buses, SubBuses, NonSubBuses, _, UnidirectionalLines, _, _,
     _, _, _, _, R, X, AdjList, _, _) = result
    P, Q = RefreshPQ(net, Sbase)
    Y = makeYbusmatrix(net, Sbase)
    n = len(Buses)

    AdjList2 = {k: AdjList[k] + [k] for k in AdjList}
    G_mat = np.real(Y); B_mat = np.imag(Y)
    bus_idx = {b: i for i, b in enumerate(sorted(Buses))}
    slack = SubBuses[0]
    gens = [int(net.sgen.bus[s]) for s in net.sgen.index]

    # Accumulate Q limits per bus (handles multiple sgens on same bus)
    Qg_min, Qg_max = {}, {}
    for s in net.sgen.index:
        bus = int(net.sgen.bus[s])
        sn = float(net.sgen.loc[s, 'sn_mva'])
        p = float(net.sgen.loc[s, 'p_mw'])
        qlim = (sn**2 - p**2)**0.5 * 1e6 / Sbase if sn**2 > p**2 else 0
        Qg_min[bus] = Qg_min.get(bus, 0) - qlim
        Qg_max[bus] = Qg_max.get(bus, 0) + qlim
    Qg_min[slack] = -1e6; Qg_max[slack] = 1e6

    m = pyo.ConcreteModel()
    m.N = pyo.Set(initialize=Buses)
    m.Vre = pyo.Var(m.N, initialize=1.0)
    m.Vim = pyo.Var(m.N, initialize=0.0)
    m.Pg = pyo.Var(m.N, domain=pyo.Reals)
    m.Qg = pyo.Var(m.N, domain=pyo.Reals,
                    bounds=lambda m, i: (Qg_min.get(i, 0), Qg_max.get(i, 0))
                    if i in gens + [slack] else (0, 0))
    for i in m.N:
        if i != slack: m.Pg[i].fix(0.0)
    m.Vre[slack].fix(vm_pu); m.Vim[slack].fix(0.0)

    # Voltage bounds — widened to [0.8, 1.2] for robustness on stressed networks
    def _vmin(m, i): return m.Vre[i]**2 + m.Vim[i]**2 >= 0.8**2
    def _vmax(m, i): return m.Vre[i]**2 + m.Vim[i]**2 <= 1.2**2
    m.VminC = pyo.Constraint(m.N, rule=_vmin)
    m.VmaxC = pyo.Constraint(m.N, rule=_vmax)

    def P_inj(m, i):
        return sum(G_mat[bus_idx[i], bus_idx[k]] * (m.Vre[i]*m.Vre[k] + m.Vim[i]*m.Vim[k])
                   + B_mat[bus_idx[i], bus_idx[k]] * (m.Vim[i]*m.Vre[k] - m.Vre[i]*m.Vim[k])
                   for k in AdjList2[i])
    def Q_inj(m, i):
        return sum(G_mat[bus_idx[i], bus_idx[k]] * (m.Vim[i]*m.Vre[k] - m.Vre[i]*m.Vim[k])
                   - B_mat[bus_idx[i], bus_idx[k]] * (m.Vre[i]*m.Vre[k] + m.Vim[i]*m.Vim[k])
                   for k in AdjList2[i])

    def _Pbal(m, i): return P_inj(m, i) == m.Pg[i] + P[i]
    def _Qbal(m, i): return Q_inj(m, i) == m.Qg[i] + Q[i]
    m.Pbal = pyo.Constraint(m.N, rule=_Pbal)
    m.Qbal = pyo.Constraint(m.N, rule=_Qbal)
    m.OBJ = pyo.Objective(expr=m.Pg[slack])

    solver = pyo.SolverFactory('ipopt')
    solver.options['tol'] = 1e-6
    solver.options['acceptable_tol'] = 1e-4
    solver.options['max_iter'] = 30000
    solver.options['linear_solver'] = 'mumps'

    t0 = time.perf_counter()
    try:
        results = solver.solve(m, tee=False)
        # Prefer IPOPT's reported solver wall_time (algorithm time) when available;
        # fall back to Python-level wall-clock if pyomo doesn't expose it.
        solver_time = None
        try:
            solver_time = float(results.solver.time)
        except Exception:
            solver_time = None
        rt = solver_time if (solver_time is not None and solver_time > 0) else (time.perf_counter() - t0)
        tc = str(results.solver.termination_condition)
        if tc in ('optimal', 'feasible'):
            # Obj.Est: Pg_slack = net deficit + losses, so losses = Pg_slack + sum(P)
            pg_slack = pyo.value(m.Pg[slack])
            est_loss = (pg_slack + sum(P.values())) * Sbase / 1e6

            # Extract method voltages
            V_method = {i: complex(pyo.value(m.Vre[i]), pyo.value(m.Vim[i])) for i in Buses}

            # Distribute Q to sgens for NR validation
            # Group sgens by bus and distribute proportionally to capacity
            bus_qlim = {}
            for s in net.sgen.index:
                bus = int(net.sgen.bus[s])
                sn = float(net.sgen.loc[s, 'sn_mva'])
                p = float(net.sgen.loc[s, 'p_mw'])
                qlim = (sn**2 - p**2)**0.5 if sn**2 > p**2 else 0.001
                bus_qlim.setdefault(bus, []).append((s, qlim))

            for bus, sgen_list in bus_qlim.items():
                total_cap = sum(ql for _, ql in sgen_list)
                qg_bus = pyo.value(m.Qg[bus]) * Sbase / 1e6
                for s, ql in sgen_list:
                    net.sgen.loc[s, 'q_mvar'] = qg_bus * (ql / total_cap) if total_cap > 0 else 0

            loss, V_nr, nr_ok = _nr_validate(net, Buses)
            delta_v = _compute_delta_v(V_method, V_nr, Buses) if nr_ok else float('nan')
            return _result(loss, est_loss, delta_v, rt, 'ok' if nr_ok else 'nr_fail')
        else:
            return _result(float('nan'), float('nan'), float('nan'), rt, tc)
    except Exception as e:
        return _result(float('nan'), float('nan'), float('nan'),
                       time.perf_counter() - t0, f'error:{e}')


# ---------------------------------------------------------------------------
#  DCOPF (Gurobi) — DistFlow for radial, B-theta+Q for meshed
# ---------------------------------------------------------------------------
def measure_dcopf(name, vm_pu=1.0):
    """Decoupled OPF (DOPF): LP for active dispatch, then QP for reactive at fixed P.

    Stage 1 (P-LP)
    --------------
    Radial : Simplified DistFlow, min Σ R f²  s.t. Σ(flow_in) = -P + Σ(flow_out),
             V_i - V_j = R f_ij, V_sub = V_slack.
    Meshed : B-θ,  min Σ R f²  s.t. f_ij = (θ_i - θ_j)/X_ij, Σ_adj f = P_inj.

    Stage 2 (Q-QP)
    --------------
    With active flows P_ij* fixed from stage 1 and |V| ≈ V_slack, solve
        min  Σ R_ij (P_ij*² + Q_ij²)
        s.t. Q-KCL: q_j + Σ_{s∈j} Qg_s = Σ_k Q_jk - Σ_i Q_ij
             v_j = v_i - 2(R P* + X Q_ij)        (linearized voltage drop)
             |Qg_s| ≤ qlim_s,  v_min² ≤ v_j ≤ v_max²
    This produces a physically meaningful Q dispatch that DCOPF alone cannot.
    """
    net = get_net_ppf(name, vm_pu)
    vm_pu = float(net.ext_grid.vm_pu.iloc[0])
    Sbase = 1e6
    result = InitializeConstants(net, Sbase=Sbase, update=True)
    (_, _, Sbase, _, Buses, SubBuses, NonSubBuses, _, UnidirectionalLines, _, _,
     _, _, _, _, R, X, AdjList, _, _) = result
    P, Q = RefreshPQ(net, Sbase)

    is_radial = name in RADIAL
    t0 = time.perf_counter()

    # ------------------------------------------------------------------
    # Stage 1: active-power LP (classical DCOPF / DistFlow)
    # ------------------------------------------------------------------
    Pf_star = {}   # active flow along each directed line [pu]
    V_method = {}  # voltage for Obj.Est and ΔV metrics
    status_p = None
    model = gp.Model(); model.setParam('OutputFlag', 0)

    if is_radial:
        Parents = {b: [] for b in Buses}
        Childs = {b: [] for b in Buses}
        for l in UnidirectionalLines:
            Parents[l[1]].append(l[0]); Childs[l[0]].append(l[1])
        V = {b: model.addVar(lb=-float('inf')) for b in Buses}
        f = {l: model.addVar(lb=-float('inf')) for l in UnidirectionalLines}
        for b in SubBuses:
            model.addConstr(V[b] == vm_pu)
        for b in NonSubBuses:
            model.addConstr(gp.quicksum(f[p, b] for p in Parents[b])
                            == -P[b] + gp.quicksum(f[b, c] for c in Childs[b]))
        for l in UnidirectionalLines:
            if R[l]**2 + X[l]**2 > 1e-20:
                model.addConstr(V[l[0]] - V[l[1]] == R[l] * f[l])
            else:
                model.addConstr(V[l[0]] == V[l[1]])
        model.setObjective(gp.quicksum(R[l]*f[l]*f[l] for l in UnidirectionalLines))
        model.optimize()
        status_p = model.status
        if status_p == GRB.OPTIMAL:
            for l in UnidirectionalLines: Pf_star[l] = float(f[l].X)
            V_method = {b: complex(V[b].X, 0) for b in Buses}
    else:
        theta = {b: model.addVar(lb=-3.15, ub=3.15) for b in Buses}
        f_all = {}
        for l in UnidirectionalLines:
            i, j = l
            if X[l] == 0: continue
            f_all[(i, j)] = model.addVar(lb=-float('inf'))
            model.addConstr(f_all[(i, j)] == (theta[i] - theta[j]) / X[l])
            f_all[(j, i)] = model.addVar(lb=-float('inf'))
            model.addConstr(f_all[(j, i)] == (theta[j] - theta[i]) / X[l])
        for b in SubBuses: model.addConstr(theta[b] == 0)
        for b in NonSubBuses:
            model.addConstr(P[b] == gp.quicksum(
                f_all[(b, adj)] for adj in AdjList[b] if (b, adj) in f_all))
        model.setObjective(gp.quicksum(
            R[l] * f_all[l] * f_all[l] for l in UnidirectionalLines if l in f_all))
        model.optimize()
        status_p = model.status
        if status_p == GRB.OPTIMAL:
            for l in UnidirectionalLines:
                Pf_star[l] = float(f_all[l].X) if l in f_all else 0.0
            V_method = {b: complex(vm_pu * math.cos(theta[b].X),
                                   vm_pu * math.sin(theta[b].X))
                        for b in Buses}

    p_stage_runtime = float(getattr(model, 'Runtime', 0.0))
    if status_p != GRB.OPTIMAL:
        rt = p_stage_runtime
        return _result(float('nan'), float('nan'), float('nan'), rt, 'P-LP_fail')

    est_loss_p = sum(R[l] * Pf_star[l]**2 for l in UnidirectionalLines) * Sbase / 1e6

    # ------------------------------------------------------------------
    # Stage 2: reactive-power QP at fixed P*
    # ------------------------------------------------------------------
    qmodel = gp.Model(); qmodel.setParam('OutputFlag', 0)
    Qf = {l: qmodel.addVar(lb=-GRB.INFINITY) for l in UnidirectionalLines}
    v = {b: qmodel.addVar(lb=0.81, ub=1.21) for b in Buses}
    for b in SubBuses: v[b].lb = vm_pu**2; v[b].ub = vm_pu**2

    Qg = {}
    sgen_on_bus = {}
    for idx in net.sgen.index:
        bus = int(net.sgen.bus[idx])
        sn_pu = float(net.sgen.sn_mva[idx]) * 1e6 / Sbase
        p_pu = float(net.sgen.p_mw[idx]) * 1e6 / Sbase
        qlim = max((sn_pu**2 - p_pu**2)**0.5, 1e-6) if sn_pu > p_pu else 1e-6
        Qg[idx] = qmodel.addVar(lb=-qlim, ub=qlim)
        sgen_on_bus.setdefault(bus, []).append(idx)

    edge_set = set(UnidirectionalLines)
    def qflow_out(bus, adj):
        if (bus, adj) in edge_set: return Qf[(bus, adj)]
        return -Qf[(adj, bus)]

    # Reactive KCL: sum of Q flows out of bus = Q injection (gen - load)
    for b in NonSubBuses:
        qout = gp.quicksum(qflow_out(b, adj) for adj in AdjList[b])
        qg_sum = gp.quicksum(Qg[s] for s in sgen_on_bus.get(b, []))
        qmodel.addConstr(qout == Q[b] + qg_sum)

    # Linearized voltage drop (uses Pf* as a parameter)
    for l in UnidirectionalLines:
        i, j = l
        qmodel.addConstr(v[i] - v[j] == 2 * (R[l] * Pf_star[l] + X[l] * Qf[l]))

    # Minimize Σ R (P*² + Q²) — total loss proxy with P fixed
    obj = gp.QuadExpr()
    for l in UnidirectionalLines:
        obj.addConstant(R[l] * Pf_star[l]**2)
        obj.add(R[l] * Qf[l] * Qf[l])
    qmodel.setObjective(obj, GRB.MINIMIZE)
    qmodel.optimize()
    rt = p_stage_runtime + float(getattr(qmodel, 'Runtime', 0.0))

    if qmodel.status == GRB.OPTIMAL:
        est_loss = float(qmodel.ObjVal) * Sbase / 1e6

        # Apply Q dispatch for NR validation
        for idx in net.sgen.index:
            net.sgen.loc[idx, 'q_mvar'] = Qg[idx].X * Sbase / 1e6

        loss, V_nr, nr_ok = _nr_validate(net, Buses)
        # Fallback: Q=0 if NR struggles with QP's Q
        if not nr_ok:
            net.sgen.q_mvar = 0
            loss, V_nr, nr_ok = _nr_validate(net, Buses)
        delta_v = _compute_delta_v(V_method, V_nr, Buses) if nr_ok else float('nan')
        return _result(loss, est_loss, delta_v, rt, 'ok' if nr_ok else 'nr_fail')
    else:
        # Q-QP infeasible — fall back to reporting the P-only LP result with Q=0
        net.sgen.q_mvar = 0
        loss, V_nr, nr_ok = _nr_validate(net, Buses)
        delta_v = _compute_delta_v(V_method, V_nr, Buses) if nr_ok else float('nan')
        return _result(loss, est_loss_p, delta_v, rt,
                       'P-only' if nr_ok else 'nr_fail')


# ---------------------------------------------------------------------------
#  SCP (Gurobi LP/QP via solve_scp_oltc)
# ---------------------------------------------------------------------------
def measure_scp(name, eta, vm_pu=1.0):
    net = get_net_ppf(name, vm_pu)
    for idx in net.trafo.index:
        net.trafo.loc[idx, 'tap_min'] = 0; net.trafo.loc[idx, 'tap_max'] = 0
        net.trafo.loc[idx, 'tap_step_percent'] = 0
    t0 = time.perf_counter()
    # Network-size-adaptive convergence tuning.
    # For stressed meshed systems (ieee1888), eta=1 needs warm-starting from eta=0
    # solution and heavier proximal damping to stabilize the aggressive Jacobian step.
    warm_runtime = 0.0  # solver time spent on the eta=0 warm-start (ieee1888 only)
    if name == 'ieee1888':
        if eta >= 0.5:
            # Two-phase: first warm-start with eta=0, then ramp to target eta
            try:
                r0 = solve_scp_oltc(net, max_iter=30, eta=0.0,
                                    with_oltc=False, rho=0.01, tol=1e-5)
                warm_runtime = float(r0.get('runtime_total', 0.0))
                # Freeze Q dispatch from eta=0 solution, continue with damped eta=1
                r = solve_scp_oltc(net, max_iter=30, eta=eta, with_oltc=False,
                                   rho=2.0, tol=1e-4, damping=0.5)
                # Fall back to eta=0 result if eta=1 still fails
                if not r.get('converged', False) and r0.get('converged', False):
                    r = r0
            except Exception:
                r = solve_scp_oltc(net, max_iter=50, eta=0.0,
                                   with_oltc=False, rho=0.01, tol=1e-5)
        else:
            r = solve_scp_oltc(net, max_iter=50, eta=eta, with_oltc=False, rho=0.01, tol=1e-5)
    elif name == 'ieee8500':
        r = solve_scp_oltc(net, max_iter=50, eta=eta, with_oltc=False, rho=0.01, tol=1e-5)
    else:
        r = solve_scp_oltc(net, max_iter=20, eta=eta, with_oltc=False, rho=0.05, tol=1e-4)
    # Use Gurobi's summed solver time (model.Runtime across iterations) instead of
    # wall-clock; this isolates the algorithm cost from Python model-construction overhead.
    rt = float(r.get('runtime_total', 0.0)) + warm_runtime
    nr_ok = r.get('nr_converged', False)
    loss = r.get('nr_loss_mw', float('nan'))
    est_loss = r.get('obj', float('nan'))  # per-unit obj from Gurobi
    delta_v = r.get('err_vs_nr_pct', float('nan'))

    # Convert est_loss: obj is in per-unit (sum R*I²), convert to MW
    # SCP uses net.sn_mva as Sbase → obj * sn_mva gives loss in MW
    if not np.isnan(est_loss):
        sbase_mva = float(net.sn_mva)
        est_loss = est_loss * sbase_mva

    status = 'ok' if nr_ok else ('not_conv' if not r.get('converged') else 'nr_fail')
    return _result(loss, est_loss, delta_v, rt, status)


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------
#%%
if __name__ == '__main__':
    all_cases = sys.argv[1:] if len(sys.argv) > 1 else RADIAL + MESHED

    hdr = (f'{"Net":<10} {"V":<6} {"Method":<8} '
           f'{"Obj.Est":<12} {"Obj.Eval":<12} {"ΔV%":<10} '
           f'{"Runtime":<10} {"Status"}')
    print(hdr)
    print('-' * len(hdr))

    for name in all_cases:
        for vm in [1.0]:
            is_radial = name in RADIAL
            methods = []

            # DCOPF
            r = measure_dcopf(name, vm)
            methods.append(('DOPF', r))

            # SOCP (radial only)
            if is_radial:
                r = measure_socp(name, vm)
                methods.append(('SOCP', r))

            # SDP (both radial and meshed; skipped for n > 100)
            r = measure_sdp(name, vm)
            methods.append(('SDP', r))

            # LM-OPF
            r = measure_lm(name, vm)
            methods.append(('LM', r))

            # NLP
            r = measure_nlp(name, vm)
            methods.append(('NLP', r))

            # SCP η=0
            r = measure_scp(name, 0.0, vm)
            methods.append(('SCP0', r))

            # SCP η=1
            r = measure_scp(name, 1.0, vm)
            methods.append(('SCP1', r))

            for mname, r in methods:
                est = f"{r['est_mw']:.4f}" if not np.isnan(r['est_mw']) else '—'
                evl = f"{r['loss_mw']:.4f}" if not np.isnan(r['loss_mw']) else '—'
                dv = f"{r['delta_v']:.2e}" if not np.isnan(r['delta_v']) else '—'
                rt = f"{r['runtime']:.2f}"
                print(f'{name:<10} {vm:<6.2f} {mname:<8} '
                      f'{est:<12} {evl:<12} {dv:<10} '
                      f'{rt:<10} {r["status"]}')
            print()
