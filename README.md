# Successive Convex Programming for Loss-Minimizing OPF

Reference implementation and benchmark suite for the paper

> *A Successive Linear/Convex Programming Framework under a Voltage–Current
> Representation for Loss-Minimizing Optimal Power Flow*
> — Seungchan Jo, Jae-Young Oh, and Gyu-Sub Lee
> (*International Journal of Electrical Power & Energy Systems*, 2026,
> manuscript IJEPES-D-26-00259).

Compares the proposed contraction-certified Successive Convex Programming
(SCP) against five established OPF baselines on six standard IEEE /
pandapower test systems.

## Methods

| Method | Solver | Formulation | Notes |
|--------|--------|-------------|-------|
| **Lin-DistFlow** | Gurobi (LP + QP) | Simplified DistFlow LP(P) → QP(Q\|P*) | Radial only; decoupled active LP then reactive QP |
| **DOPF** | Gurobi (LP + QP) | B-θ LP(P) → QP(Q\|P*) | Meshed only; decoupled active LP then reactive QP |
| **SOCP** | Gurobi (QCQP) | Farivar–Low DistFlow | Radial only; exact relaxation under loss-min |
| **SDP** | CVXPY + MOSEK | BIM with `W ≽ 0` | Limited to networks with ≤ 200 buses |
| **LM-OPF** | Gurobi (QP) | Linear BFM + cycle recovery | Wang, Shrestha, Dubey (2026, IEEE TSG early access) |
| **NLP** | Pyomo + IPOPT | Full rectangular AC | |
| **SCP** | Gurobi (LP / QP) | Proposed iterative scheme | η ∈ {0, 1} (Picard / full sensitivity) |

### Mathematical formulations

Let

- `V_i = v_r_i + j v_im_i` : complex bus voltage (pu)
- `v_i = |V_i|²`             : squared magnitude (pu)
- `P_ij, Q_ij`                : branch active/reactive flow (pu)
- `ℓ_ij = |I_ij|²`            : squared branch current magnitude (pu)
- `R_ij, X_ij`                : branch resistance/reactance (pu)
- `p_i, q_i`                  : net nodal active/reactive injection (gen − load)
- `Q_g^s`                     : controllable reactive dispatch at sgen s (pu)
- `q_lim^s = √(sn_s² − p_s²)` : reactive capability

#### 1. Lin-DistFlow / DOPF — Decoupled OPF

**Stage 1 (active LP).**

- Radial (Simplified DistFlow — referred to as **Lin-DistFlow**):

  ```
  min   Σ R_ij f_ij²
  s.t.  Σ_{p: p→i} f_pi + p_i = Σ_{c: i→c} f_ic       ∀ i ∉ slack
        V_i − V_j = R_ij f_ij                          ∀ (i, j)
        V_slack = V*_slack
  ```

- Meshed (B–θ — referred to as **DOPF**):

  ```
  min   Σ R_ij f_ij²
  s.t.  f_ij = (θ_i − θ_j) / X_ij                      ∀ (i, j)
        Σ_{j∈adj(i)} f_ij = p_i                        ∀ i ∉ slack,  θ_slack = 0
  ```

**Stage 2 (reactive QP with `P_ij* ` fixed from Stage 1).**

```
min   Σ R_ij (P_ij*² + Q_ij²)
s.t.  q_i + Σ_{s ∈ i} Q_g^s = Σ_{k: i→k} Q_ik − Σ_{j: j→i} Q_ji   (reactive KCL)
      v_i − v_j = 2 (R_ij P_ij* + X_ij Q_ij)                       (lin. voltage drop)
      v_min² ≤ v_i ≤ v_max²,   |Q_g^s| ≤ q_lim^s
```

#### 2. SOCP — Farivar-Low DistFlow (radial)

```
min   Σ R_ij ℓ_ij
s.t.  p_j = Σ_{k: j→k} P_jk − Σ_{i: i→j} (P_ij − R_ij ℓ_ij)
      q_j = Σ_{k: j→k} Q_jk − Σ_{i: i→j} (Q_ij − X_ij ℓ_ij) + Σ_{s∈j} Q_g^s
      v_j = v_i − 2 (R_ij P_ij + X_ij Q_ij) + |z_ij|² ℓ_ij
      P_ij² + Q_ij² ≤ v_i · ℓ_ij                                  (SOC relaxation)
      v_min² ≤ v_i ≤ v_max²,   |Q_g^s| ≤ q_lim^s,   v_slack = V*_slack²
```

Tight for loss-minimizing OPF on radial networks (Farivar & Low, 2013).

#### 3. SDP — BIM relaxation

Lift `W = [v_r; v_im][v_r; v_im]ᵀ ∈ S_+^{2n}`:

```
min   Σ_i (P_g^i − p_d^i)                                         (total loss)
s.t.  tr(H_p^i · W) = P_g^i − p_d^i                               ∀ i    (active KCL)
      tr(H_q^i · W) = Q_g^i − q_d^i                               ∀ i    (reactive KCL)
      v_min² ≤ tr(J_i · W) ≤ v_max²
      P_min^i ≤ P_g^i ≤ P_max^i,  Q_min^i ≤ Q_g^i ≤ Q_max^i
      tr(M_slack · W) = V*_slack²
      W ≽ 0
```

If `rank(W) = 1` the relaxation is exact. On heavily meshed systems
rank-1 is not guaranteed, so the reported "Est." for SDP is the NR-
validated loss at the SDP-recovered operating point (the raw relaxed
objective is a lower bound, not a physical loss).

#### 4. LM-OPF (Wang et al. 2026) — meshed linearized BFM

```
min   Σ R_ij (P_ij² + Q_ij²)
s.t.  p_j = Σ_{k: j→k} P_jk − Σ_{i: i→j} P_ij                        (active KCL, losses dropped)
      q_j = Σ_{k: j→k} Q_jk − Σ_{i: i→j} Q_ij + Σ_{s∈j} Q_g^s
      v_j = v_i − 2 (R_ij P_ij + X_ij Q_ij)                          (linear voltage drop)
      Σ_{(i,j)∈C}  (π/4) (X_ij P_ij − R_ij Q_ij) = 0                 ∀ cycle C
      v_min² ≤ v_i ≤ v_max²,  |Q_g^s| ≤ q_lim^s
```

The cycle constraint is the linearized angle-recovery condition obtained
from `arctan(t) ≈ (π/4)·t` for `|t| ≤ 1`.

#### 5. NLP — full AC-OPF (rectangular)

```
min   P_g^slack
s.t.  P_i(V) = P_g^i − p_d^i                                      (nonlin. AC KCL)
      Q_i(V) = Q_g^i − q_d^i
      v_min² ≤ v_r_i² + v_im_i² ≤ v_max²
      bounds on P_g, Q_g
```

Solved with IPOPT via Pyomo. May converge to local optima.

#### 6. SCP (proposed) — successive LP/QP in current-injection domain

At iterate k with reference voltage Ṽ^(k):

```
min   Σ R_ij |I_ij|²
s.t.  Ŷ V − η · J_load(Ṽ^(k)) (V − Ṽ^(k)) − fI(Ṽ^(k)) = 0           (linearized KCL)
      V_slack = V*_slack + j·0
      |I_ij|² ≤ I_max²,   voltage + Q-gen limits
Update Ṽ^(k+1) ← V*;  iterate until ‖ΔṼ‖ < tol.
```

- `η = 0` : Jacobian-free Picard iteration (tightest contraction certificate)
- `η = 1` : full sensitivity update (faster local convergence)

## Test systems

| Network | Buses | Branches | Loads | Sgens | Total P_d | Type | V_slack | Load scaling | Q cap/sgen |
|---------|-------|----------|-------|-------|-----------|------|---------|--------------|------------|
| IEEE 33 | 33 | 32 | 32 | 4 | 3.16 MW | Radial | 1.00 pu | 100% | — |
| IEEE 123 | 126 | 120 | 85 | 16 | 3.49 MW | Radial | 1.00 pu | 100% | — |
| IEEE 8500 | 4,876 | 4,875 | 1,177 | 1,213 | 16.16 MW | Radial | 1.00 pu | 100% | — |
| IEEE 14 | 14 | 20 | 11 | 4 | 259 MW | Meshed | 1.06 pu | 100% | ±30 MVAR |
| IEEE 118 | 118 | 179 | 99 | 53 | 2,121 MW | Meshed | 1.035 pu | 50% | ±30 MVAR |
| French 1888 | 1,888 | 2,308 | 943 | 271 | 17,882 MW | Meshed | 1.00 pu | 30% | ±30 MVAR |

### Data sources

All networks are constructed from standard pandapower 3.0 test cases:

- **IEEE 33**: `pandapower.networks.case33bw()` — 12.66 kV radial distribution feeder (Baran & Wu, 1989)
- **IEEE 123**: reconstructed from IEEE PES Test Feeder data (`IEEE123Node/`) — 4.16 kV three-phase feeder reduced to balanced single-phase equivalent
- **IEEE 8500**: pandapower `case_ieee8500()` — 7.2 kV large radial feeder (balanced reduction), cached as `ieee8500.p` since the native constructor requires OpenDSS/comtypes
- **IEEE 14**: `pandapower.networks.case14()` — 5-generator meshed transmission system
- **IEEE 118**: `pandapower.networks.case118()` — 118-bus transmission system; loads scaled to 50% of nominal to maintain NR feasibility after PV→PQ bus conversion
- **French 1888**: `pandapower.networks.case1888rte()` — French transmission grid (RTE); loads scaled to 30% of nominal for the same PV→PQ feasibility reason

### Preprocessing

- Originally unbalanced feeders are reduced to balanced single-phase equivalents
- Distributed generators (sgens) are placed at buses selected by each network's factory function
- Line shunt capacitances and transformer magnetizing branches are set to zero
- Bus voltage limits: [0.9, 1.1] pu throughout
- For meshed cases (IEEE 14/118, French 1888): each sgen's apparent-power rating is set to `sn = √(p² + 30²)` MVA, bounding reactive capability to ±30 MVAR per unit. This matches realistic distributed-inverter ratings and is applied uniformly to all baselines
- IEEE-118 and French-1888 loads are scaled because the qdispatch conversion (PV→PQ) removes generator voltage regulation; the native loadings require PV-bus reactive control that the PQ-bus benchmark removes by construction

## Quickstart

```bash
# 1. Install the package and its deps (Python 3.10+ recommended)
pip install -e .

# 2. Install solvers
#    * Gurobi: https://www.gurobi.com/academia/
#    * MOSEK : https://www.mosek.com/academic/
#    * IPOPT : e.g. `pip install idaes-pse && idaes get-extensions`

# 3. Run
python -m scp_opf_bench.cli ieee33 ieee14       # subset
python -m scp_opf_bench.cli                     # all 6 networks
scp-opf-bench -m DOPF,SOCP,SCP1 ieee123         # specific methods

# 4. Tests
pytest
```

Output columns per method per network:

- **Obj.Est**  — method's own loss estimate
- **Obj.Eval** — NR-validated loss after plugging the method's dispatch into pandapower's Newton–Raphson power flow
- **ΔV %**     — average **per-bus magnitude** deviation between method and NR
- **Time**     — sum of solver `model.optimize` calls (s); excludes Python model build
- **Status**   — `ok`, `nr_fail`, `infeasible`, `not_conv`, `skip:*`

## Repository layout

```
.
├── LICENSE
├── README.md
├── pyproject.toml
├── src/scp_opf_bench/
│   ├── __init__.py
│   ├── networks.py         # single source of truth for per-network policy
│   ├── nr.py               # NR validation helpers
│   ├── result.py           # MethodResult dataclass
│   ├── cli.py              # python -m scp_opf_bench.cli
│   ├── methods/            # one file per method
│   │   ├── dopf.py  socp.py  sdp.py  lm_opf.py  nlp.py  scp.py
│   │   └── _runner.py      # wraps legacy code with CWD context manager
│   └── _legacy/            # unmodified research code (imported by methods/)
└── tests/
    └── test_smoke.py
```

## Citation

```bibtex
@article{Jo2026SCPOPF,
  author  = {Jo, Seungchan and Oh, Jae-Young and Lee, Gyu-Sub},
  title   = {A Successive Linear/Convex Programming Framework under a
             Voltage--Current Representation for Loss-Minimizing Optimal
             Power Flow},
  journal = {International Journal of Electrical Power \& Energy Systems},
  year    = {2026},
  note    = {Manuscript IJEPES-D-26-00259}
}
```

## License

MIT.
