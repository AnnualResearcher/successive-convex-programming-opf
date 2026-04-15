"""Common result type returned by every measure_* function."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MethodResult:
    method: str           # "DCOPF", "SOCP", ...
    network: str          # "ieee33", ...
    est_mw: float         # method's own loss estimate [MW]
    eval_mw: float        # NR-validated loss [MW]
    delta_v_pct: float    # avg |V| deviation from NR [%]
    runtime_s: float      # wall-clock seconds
    status: str           # "ok" | "nr_fail" | "infeasible" | "not_conv" | "skip:..."

    def fmt_row(self) -> str:
        def f(x, prec=4):
            return f"{x:.{prec}f}" if x == x else "—"   # NaN -> em-dash
        return (
            f"{self.network:<10} {self.method:<8} "
            f"{f(self.est_mw):<12} {f(self.eval_mw):<12} "
            f"{(f'{self.delta_v_pct:.2e}' if self.delta_v_pct == self.delta_v_pct else '—'):<10} "
            f"{self.runtime_s:>7.2f}  {self.status}"
        )


HEADER = (
    f"{'Net':<10} {'Method':<8} {'Obj.Est':<12} {'Obj.Eval':<12} "
    f"{'ΔV %':<10} {'Time':>7}  Status"
)
