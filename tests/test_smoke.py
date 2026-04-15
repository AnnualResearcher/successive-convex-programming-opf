"""Smoke tests: verify each method at least runs on IEEE-33 without error."""

import pytest

from scp_opf_bench import build_net
from scp_opf_bench.methods import (
    measure_dopf, measure_socp, measure_sdp, measure_lm_opf,
    measure_nlp, measure_scp,
)


@pytest.mark.parametrize("case", ["ieee33", "ieee14"])
def test_network_builds(case):
    net = build_net(case)
    assert len(net.bus) > 0
    assert len(net.load) > 0
    assert net.ext_grid.vm_pu.iloc[0] > 0.9


def test_dopf_ieee33():
    r = measure_dopf("ieee33")
    assert r.status == "ok"
    assert 0.05 < r.eval_mw < 0.20  # ~85 kW


def test_socp_ieee33():
    r = measure_socp("ieee33")
    assert r.status == "ok"
    # SOCP should be tight on radial
    assert 0.08 < r.eval_mw < 0.10
    assert abs(r.est_mw - r.eval_mw) < 0.001


def test_scp_eta1_ieee33():
    r = measure_scp("ieee33", 1.0)
    assert r.status == "ok"
    assert 0.08 < r.eval_mw < 0.10
    # SCP iteration should drive ΔV to machine precision
    assert r.delta_v_pct < 0.01
