"""Command-line entry point: ``python -m scp_opf_bench.cli ieee33 ieee14 ...``."""

from __future__ import annotations

import argparse
import sys

from .networks import ALL_CASES, RADIAL, is_radial
from .result import HEADER
from .methods import (
    measure_dopf, measure_socp, measure_sdp, measure_lm_opf,
    measure_nlp, measure_scp,
)


def run_one(name: str, methods: tuple[str, ...]) -> None:
    results = []
    for m in methods:
        if m == "DOPF":
            results.append(measure_dopf(name))
        elif m == "SOCP":
            if is_radial(name):
                results.append(measure_socp(name))
        elif m == "SDP":
            results.append(measure_sdp(name))
        elif m == "LM-OPF":
            results.append(measure_lm_opf(name))
        elif m == "NLP":
            results.append(measure_nlp(name))
        elif m == "SCP0":
            results.append(measure_scp(name, 0.0))
        elif m == "SCP1":
            results.append(measure_scp(name, 1.0))
        else:
            raise ValueError(f"Unknown method: {m}")

    for r in results:
        print(r.fmt_row())


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="SCP-OPF benchmark runner")
    p.add_argument("networks", nargs="*", default=list(ALL_CASES),
                   help=f"networks to test (default: all). Choices: {ALL_CASES}")
    p.add_argument("--methods", "-m", default="DOPF,SOCP,SDP,LM-OPF,NLP,SCP0,SCP1",
                   help="comma-separated methods to run")
    args = p.parse_args(argv)

    unknown = [n for n in args.networks if n not in ALL_CASES]
    if unknown:
        print(f"Unknown networks: {unknown}")
        print(f"Available: {ALL_CASES}")
        return 1

    methods = tuple(m.strip() for m in args.methods.split(",") if m.strip())
    print(HEADER)
    print("-" * len(HEADER))
    for name in args.networks:
        run_one(name, methods)
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
