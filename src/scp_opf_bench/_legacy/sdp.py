#%%
# sdp_opf.py
# SDP relaxation of AC-OPF with W = [Re(V); Im(V)] [Re(V); Im(V)]^T  (2n x 2n)
# Requirements: pip install cvxpy numpy

import numpy as np
import cvxpy as cp


def build_power_matrices(Ybus: np.ndarray):
    """
    Build linear maps for active/reactive power injections:
        P_i = trace(Hp[i] @ W),  Q_i = trace(Hq[i] @ W)
    where W = [vr; vi] [vr; vi]^T ∈ S_+^{2n}
    Ybus = G + jB (n x n)
    """
    assert Ybus.shape[0] == Ybus.shape[1], "Ybus must be square."
    n = Ybus.shape[0]
    G = np.asarray(np.real(Ybus))
    B = np.asarray(np.imag(Ybus))

    Hp = []  # list of 2n x 2n matrices
    Hq = []

    Mk = []
    Mlm = []
    yMlm = []

    I_n = np.eye(n)

    for i in range(n):
        E = np.zeros((n, n))
        E[i, i] = 1.0

        # ----- Active power P_i -----
        # Blocks for z = [vr; vi]
        H11_P = E @ G                                 # vr^T E G vr
        H22_P = E @ G                                 # vi^T E G vi
        # Cross terms: - vr^T E B vi + vi^T E B vr = vr^T(-E B + B^T E) vi
        H12_P = 0.5 * (-E @ B + B.T @ E)

        Hp_i = np.block([
            [H11_P,  H12_P     ],
            [H12_P.T, H22_P    ]
        ])
        Hp.append(Hp_i)

        # ----- Reactive power Q_i -----
        # Quadratic terms: - vr^T E B vr  and  - vi^T E B vi
        sym_EB = 0.5 * (E @ B + B.T @ E)             # symmetrized
        H11_Q = -sym_EB
        H22_Q = -sym_EB
        # Cross terms:  vi^T E G vr - vr^T E G vi = vr^T( G^T E - E G ) vi
        H12_Q = 0.5 * (G.T @ E - E @ G)

        Hq_i = np.block([
            [H11_Q,  H12_Q     ],
            [H12_Q.T, H22_Q    ]
        ])
        Hq.append(Hq_i)

        # Mk
        e = np.zeros((n,1))     
        e[i][0] =1
        
        M11_k = e@e.T
        Mk_i = np.block([[M11_k, np.zeros_like(M11_k)], [np.zeros_like(M11_k), M11_k]])
        Mk.append(Mk_i)



    for l in range(n):
        for m in range(n):
            if abs(Ybus[l][m]) > 1e-8 and m > l:
                el, em = np.zeros((n,1)), np.zeros((n,1))
                M12_lm = el@el.T
                el[l][0], em[m][0] = 1, 1
                elm = el - em
                M11_lm  = elm@elm.T

                Mlm_i = np.block([[M11_lm, M12_lm], [M12_lm, M11_lm]])
                Mlm.append(Mlm_i)

                yMlm.append(-Ybus[l][m].real*Mlm_i)



    # Voltage magnitude selector J_i: |V_i|^2 = vr_i^2 + vi_i^2 = trace(J_i W)
    J = []
    for i in range(n):
        Ji = np.zeros((2*n, 2*n))
        Ji[i, i] = 1.0            # vr_i^2
        Ji[n+i, n+i] = 1.0        # vi_i^2
        J.append(Ji)

    return Hp, Hq, Mk, Mlm, yMlm, J # 논문의 Y_k의 real part : HP, 논문의 Y_k의 imaginary part: Hq


def sdp_opf(
    Ybus: np.ndarray,
    Pd: np.ndarray, Qd: np.ndarray,        # load demand (+) at each bus
    gen_buses: list,                       # list of generator bus indices (0-based)
    Pmin: np.ndarray, Pmax: np.ndarray,    # only meaningful at gen_buses; ignored elsewhere
    Qmin: np.ndarray, Qmax: np.ndarray,
    Vmin: np.ndarray, Vmax: np.ndarray,    # voltage mag bounds [pu]
    cost_quad: np.ndarray = None,          # c2 for Pg^2, length n; use 0 for non-gens or if not used
    cost_lin:  np.ndarray = None,          # c1 for Pg
    cost_const: np.ndarray = None,         # c0
    solver: str = "MOSEK",                 # "MOSEK" recommended; fallback "SCS"
    verbose: bool = False,
    slack_buses: list = None,              # list of slack bus indices (ext_grid)
    slack_vm: float = 1.0,                 # slack bus voltage magnitude
):
    """
    Solve SDP-OPF (relaxation). Returns dict with results.
    Conventions:
      - Injections:  P_inj[i] = Pg[i] - Pd[i],  Q_inj[i] = Qg[i] - Qd[i]
      - For non-gen buses, Pg=Qg=0 so P_inj = -Pd, Q_inj = -Qd (equality).
      - For gen buses i∈G:  Pmin[i] ≤ Pg[i] ≤ Pmax[i], Qmin[i] ≤ Qg[i] ≤ Qmax[i]
                            => Pmin-Pd ≤ P_inj ≤ Pmax-Pd (same for Q)
    """
    n = Ybus.shape[0]
    Pd = Pd.astype(float).reshape(n)
    Qd = Qd.astype(float).reshape(n)
    Vmin = Vmin.astype(float).reshape(n)
    Vmax = Vmax.astype(float).reshape(n)

    # Defaults for cost: minimize total generation (≈ losses) if not given
    if cost_quad is None: cost_quad = np.zeros(n)
    if cost_lin  is None: cost_lin  = np.ones(n) * 1.0
    if cost_const is None: cost_const = np.zeros(n)
    cost_quad = cost_quad.reshape(n)
    cost_lin  = cost_lin.reshape(n)
    cost_const= cost_const.reshape(n)

    Hp, Hq, Mk, Mlm, yMlm, J = build_power_matrices(Ybus)

    # Decision: W ≽ 0  (2n x 2n)
    W = cp.Variable((2*n, 2*n), PSD=True)

    # Linear injection functionals
    P_inj = [cp.trace(Hp[i] @ W) for i in range(n)]
    Q_inj = [cp.trace(Hq[i] @ W) for i in range(n)]

    constr = []

    # Voltage magnitude bounds
    for i in range(n):
        constr += [ Vmin[i]**2 <= cp.trace(J[i] @ W),
                    cp.trace(J[i] @ W) <= Vmax[i]**2 ]

    # Power balance & gen limits
    genset = set(gen_buses or [])
    slackset = set(slack_buses or [])
    for i in range(n):
        if (i in genset) and (i not in slackset):
            # Sgen bus: fixed P (via Pmin==Pmax), variable Q
            constr += [ Pmin[i] - Pd[i] <= P_inj[i],
                        P_inj[i] <= Pmax[i] - Pd[i] ]
            constr += [ Qmin[i] - Qd[i] <= Q_inj[i],
                        Q_inj[i] <= Qmax[i] - Qd[i] ]
        elif i in slackset:
            # Slack bus: variable P and Q, fixed voltage
            constr += [ Pmin[i] - Pd[i] <= P_inj[i],
                        P_inj[i] <= Pmax[i] - Pd[i] ]
            constr += [ Qmin[i] - Qd[i] <= Q_inj[i],
                        Q_inj[i] <= Qmax[i] - Qd[i] ]
            constr += [ cp.trace(Mk[i] @ W) == slack_vm**2 ]
        else:
            # No generator: injections fixed to -load
            constr += [ P_inj[i] == -Pd[i],
                        Q_inj[i] == -Qd[i] ]

    

    # Objective: minimize total network losses = sum of nodal active-power injections
    # (Kirchhoff: sum P_inj across all buses = total real power dissipation in lines)
    obj = cp.Minimize(cp.sum(P_inj))

    prob = cp.Problem(obj, constr)

    # Solve
    try:
        prob.solve(solver=getattr(cp, solver), verbose=verbose)
    except Exception as e:
        # Fallback to SCS if MOSEK unavailable or fails
        print(f'{solver} failed: {e}, falling back to SCS')
        try:
            prob.solve(solver=cp.SCS, verbose=verbose, max_iters=5000, eps=1e-6)
        except Exception as e2:
            print(f'SCS also failed: {e2}')

    status = prob.status
    optval = prob.value

    # Extract numbers
    P_inj_val = np.array([p.value if hasattr(p, "value") else None for p in P_inj], dtype=float)
    Q_inj_val = np.array([q.value if hasattr(q, "value") else None for q in Q_inj], dtype=float)
    V2_val    = np.array([cp.trace(J[i] @ W).value for i in range(n)], dtype=float)

    # Simple (heuristic) rank-1 recovery: leading eigenvector of W
    W_val = W.value
    if W_val is None:
        return {
            "status": status,
            "objective": optval,
            "W": None,
            "P_injection": P_inj_val,
            "Q_injection": Q_inj_val,
            "Pg": P_inj_val + Pd if P_inj_val[0] is not None else np.full(n, np.nan),
            "Qg": Q_inj_val + Qd if Q_inj_val[0] is not None else np.full(n, np.nan),
            "V2": V2_val,
            "Vmag_principal": np.full(n, np.nan),
            "V_complex_principal": np.full(n, np.nan, dtype=complex),
        }
    evals, evecs = np.linalg.eigh((W_val + W_val.T) * 0.5)
    idx = np.argmax(evals)
    z_hat = evecs[:, idx] * np.sqrt(max(evals[idx], 0.0))   # principal component
    vr_hat = z_hat[:n]
    vi_hat = z_hat[n:]
    V_hat = vr_hat + 1j*vi_hat
    Vmag_hat = np.abs(V_hat)

    return {
        "status": status,
        "objective": optval,
        "W": W_val,
        "P_injection": P_inj_val,                 # Pg - Pd  (pu)
        "Q_injection": Q_inj_val,                 # Qg - Qd  (pu)
        "Pg": P_inj_val + Pd,                     # inferred generation
        "Qg": Q_inj_val + Qd,
        "V2": V2_val,                             # |V|^2 from W-diagonals
        "Vmag_principal": Vmag_hat,               # heuristic voltage recovery
        "V_complex_principal": V_hat,             # heuristic complex voltage
    }

#%%
if __name__ == "__main__":
  import networks
  import gurobipy as gp
  from gurobipy import GRB
  import pandas as pd
  import pandapower as pp
  import pandapower.networks as pn
  import numpy as np
  import math
  import cmath
  from utils import InitializeConstants, RefreshPQ, makeYbusmatrix
  net = networks.ieee123(ppf= True, use_sgen = True)
#net.line.r_ohm_per_km = 0.0
#net.line.x_ohm_per_km = 0.00001
#net = networks.three_bus_example()
#net.line.x_ohm_per_km = 0.054
#net.line.r_ohm_per_km = 0.001

  net.line.r_ohm_per_km = 2.139135
  net.line.x_ohm_per_km = 2.168586
  net.line.length_km = 0.05334
  net.load.p_mw = 0.04
  net.load.q_mvar =0.02
  Sbase = 1e8
  Vbase, Ibase, Sbase, Zbase, Buses, SubstationBuses, NonsubstationBuses, BidirectionalLines, UnidirectionalLines, UnidirectionalTrafo, NonTrafoUnidirectionalLines, TrafoTappos, TrafoTapStepPercent, LineTupleToIndex, TrafoTupleToIndex, R, X, AdjacencyList, DirectedLinesForEachTree, UnidirectionalLineTupleToIndex = InitializeConstants(net, Sbase = Sbase)
  for key, value in zip(R.keys(), R.values()):
      a= (0,1)
      if value < 1e-10:
          R[key] = R[a]
          X[key] = X[a]
      else:
          a = key
  P, Q = RefreshPQ(net, Sbase)

#%%
# --------------------------
# Minimal usage example
# --------------------------
if __name__ == "__main__":
    # 3-bus toy Ybus (pu). Replace with your own system data.
    # Here we make a simple symmetric example (no shunts/line charging).
    Y = makeYbusmatrix(net, Sbase)

    n = Y.shape[0]
    Pd, Qd = [], [] 
    for bus in net.bus.index:
        Pd.append(-P[bus]), Qd.append(-Q[bus])
    
    Pd = np.array(Pd)
    Qd = np.array(Qd)


    gen_buses = list(net.ext_grid.bus) + list(net.sgen.bus)                 # bus 0 (slack), bus 1 (PV)
    gen_buses.sort()
    Pmin = np.zeros(n); Pmax = np.ones(n)*Pd.sum()*2
    Qmin = -np.ones(n)*Qd.sum()*2.; Qmax = np.ones(n)*Qd.sum()*2

    Vmin = np.ones(n)*0.9
    Vmax = np.ones(n)*1.1

    # Quadratic gen cost c2*Pg^2 + c1*Pg + c0 (only used on gen buses)
    c2 = np.zeros(n); c1 = np.zeros(n); c0 = np.zeros(n)
    c2[0] = 0.01; c1[0] = 1.0
    c2[1] = 0.02; c1[1] = 0.8

    res = sdp_opf(
        Ybus=Y, Pd=Pd, Qd=Qd,
        gen_buses=gen_buses,
        Pmin=Pmin, Pmax=Pmax, Qmin=Qmin, Qmax=Qmax,
        Vmin=Vmin, Vmax=Vmax,
        cost_quad=c2, cost_lin=c1, cost_const=c0,
        solver="MOSEK", verbose=True
    )

    print("Status:", res["status"])
    print("Objective:", res["objective"])
    print("Pg (pu):", np.round(res["Pg"], 4))
    print("Qg (pu):", np.round(res["Qg"], 4))
    print("|V| (principal):", np.round(res["Vmag_principal"], 4))

def recover_voltage_complex(W):
    n = W.shape[0]//2
    Wrr, Wri = W[:n,:n], W[:n,n:]
    Wir, Wii = W[n:,:n], W[n:,n:]
    Wc = Wrr + Wii + 1j*(Wir - Wri)         # complex lift

    m = np.sqrt(np.maximum(np.real(np.diag(Wc)), 0.0))
    H = np.zeros_like(Wc, dtype=complex)
    nz = m > 1e-12
    for i in range(n):
        for j in range(n):
            if nz[i] and nz[j]:
                H[i,j] = Wc[i,j]/(m[i]*m[j])

    # Hermitian symmetrization & leading eigenvector
    evals, evecs = np.linalg.eigh((H + H.conj().T)/2)
    u = evecs[:, np.argmax(evals)]
    phi = np.angle(u)
    s = 0  # slack index
    V = m * np.exp(1j*(phi - phi[s]))
    return V
