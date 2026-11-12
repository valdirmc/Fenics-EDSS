#!/usr/bin/env python3
"""Compara CG+GAMG vs LU no B1 n=8 e n=16. Nao inventa tolerancia de aceite.

Uso (container): python3 /app/tests/run_fc_cg_vs_lu.py --work /work --levels 8,16
"""
from __future__ import annotations

import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("FC_SRC", "/app/src"))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

import numpy as np
from fc_b1 import write_level_inp
from run_from_inp import run_from_inp


def run_one(n, work, solver):
    pack = write_level_inp(work, n)
    out = os.path.join(work, "out_n{0}_{1}".format(n, solver.replace("-", "_")))
    os.makedirs(out, exist_ok=True)
    t0 = time.perf_counter()
    results = run_from_inp(
        pack["inp"],
        output_dir=out,
        debug=1,
        profile="fc",
        units="mm-N-MPa-t-s",
        solver=solver,
        rtol=1e-10,
    )
    elapsed = time.perf_counter() - t0
    uh, extra = results[0]
    return {
        "n": n,
        "solver": solver,
        "elapsed_s": elapsed,
        "u": np.asarray(uh.x.array, dtype=float).copy(),
        "sigma": np.asarray(extra["stress_cells"], dtype=float).copy(),
        "solver_info": extra.get("solver_info"),
        "reactions": {k: v for k, v in (extra.get("reactions") or {}).items() if k != "entries"},
        "out": out,
    }


def compare(lu, cg):
    du = np.abs(cg["u"] - lu["u"])
    ds = np.abs(cg["sigma"] - lu["sigma"])
    u_scale = float(np.max(np.abs(lu["u"]))) or 1.0
    s_scale = float(np.max(np.abs(lu["sigma"]))) or 1.0
    return {
        "n": lu["n"],
        "n_dofs": int(lu["u"].size),
        "max_abs_dU": float(du.max()),
        "max_abs_dU_over_max_abs_U_lu": float(du.max() / u_scale),
        "rms_dU": float(np.sqrt(np.mean(du ** 2))),
        "max_abs_dS": float(ds.max()),
        "max_abs_dS_over_max_abs_S_lu": float(ds.max() / s_scale),
        "rms_dS": float(np.sqrt(np.mean(ds ** 2))),
        "lu_elapsed_s": lu["elapsed_s"],
        "cg_elapsed_s": cg["elapsed_s"],
        "cg_solver_info": cg["solver_info"],
        "lu_solver_info": lu["solver_info"],
        "note": "Diferencas brutas CG+GAMG vs LU; nenhuma tolerancia de aceite inventada.",
    }


def main():
    argv = sys.argv[1:]
    work = "/tmp/fc_cg_vs_lu"
    levels = [8, 16]
    if "--work" in argv:
        work = argv[argv.index("--work") + 1]
    if "--levels" in argv:
        levels = [int(x) for x in argv[argv.index("--levels") + 1].split(",") if x]
    os.makedirs(work, exist_ok=True)
    rows = []
    for n in levels:
        print("==== B1 n={0} LU ====".format(n), flush=True)
        lu = run_one(n, work, "lu")
        print("==== B1 n={0} CG+GAMG ====".format(n), flush=True)
        cg = run_one(n, work, "cg-gamg")
        cmp = compare(lu, cg)
        rows.append(cmp)
        print(json.dumps(cmp, indent=2, default=str), flush=True)
        # libera arrays grandes antes do proximo nivel
        lu["u"] = lu["sigma"] = None
        cg["u"] = cg["sigma"] = None
    outp = os.path.join(work, "cg_vs_lu_summary.json")
    payload = {"benchmark": "B1 CG+GAMG vs LU", "levels": rows}
    with open(outp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, default=str)
        fh.write("\n")
    print("SUMMARY", outp, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
