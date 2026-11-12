#!/usr/bin/env python3
"""B1 n=4/8/16 CG+GAMG: esfera geometrica vs baricentro; gera copia CCX com *EL PRINT."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src"))

from run_fc_lote3 import run_level  # noqa: E402


def _ccx_fmt_number(tok: str) -> str:
    """Compacta float para o leitor Fortran do CCX (linhas curtas, E maiusculo)."""
    s = tok.strip()
    if s == "":
        return tok
    try:
        v = float(s)
    except ValueError:
        return tok
    if v == 0.0:
        return "0"
    if 1e-4 <= abs(v) < 1e7:
        t = "{0:.12g}".format(v)
        return t
    return "{0:.8E}".format(v)


def _ccx_compact_data_line(line: str) -> str:
    # so linhas de dados (comecam com numero): compacta campos numericos
    raw = line.strip()
    if not raw or raw.startswith("*") or raw.startswith("**"):
        return line
    if not raw[0].isdigit() and raw[0] not in "+-":
        return line
    parts = [p.strip() for p in raw.split(",")]
    return ", ".join(_ccx_fmt_number(p) for p in parts)


def ccx_deck_from_fc_inp(text: str) -> str:
    extra = "*NODE PRINT, NSET=Nall\nU\n*EL PRINT, ELSET=Eall, GLOBAL=YES\nS\n"
    if "*END STEP" not in text:
        raise RuntimeError("INP sem *END STEP")
    body = text.replace("*END STEP", extra + "*END STEP", 1)
    # *HEADING evita alguns leitores sensiveis; numeros compactos (<72 cols)
    if not body.lstrip().upper().startswith("*HEADING"):
        body = "*HEADING\nB1 flexao pura C3D4 (copia CCX, nao e a roda)\n" + body
    return "\n".join(_ccx_compact_data_line(l) for l in body.splitlines()) + "\n"


def main():
    out = Path(os.environ.get("FC_B1_OUT", "/tmp/fc-b1-sphere"))
    out.mkdir(parents=True, exist_ok=True)
    levels = [int(x) for x in os.environ.get("FC_B1_LEVELS", "4,8,16").split(",") if x]
    solver = os.environ.get("FC_B1_SOLVER", "cg-gamg")
    rtol = float(os.environ.get("FC_B1_RTOL", "1e-10"))
    rows = []
    for n in levels:
        work = out / "n{0}".format(n)
        work.mkdir(parents=True, exist_ok=True)
        ev = run_level(n, str(work), debug=1, solver=solver, rtol=rtol)
        inp = Path(ev["inp"])
        ccx_path = work / "b1_n{0}_ccx.inp".format(n)
        ccx_path.write_text(ccx_deck_from_fc_inp(inp.read_text(encoding="utf-8")))
        ev["ccx_inp"] = str(ccx_path)
        slim = {
            "n_thick": ev["meta"]["n_thick"],
            "h": ev["meta"]["h"],
            "n_nodes": ev["meta"]["n_nodes"],
            "n_elements": ev["meta"]["n_elements"],
            "elapsed_s": ev["elapsed_s"],
            "solver": ev.get("solver"),
            "sphere_geometric": ev.get("sphere_geometric") or ev.get("sphere"),
            "sphere_barycenter": ev.get("sphere_barycenter"),
            "uy_mid_rel_err": ev.get("uy_mid_rel_err"),
            "phi_rel_rel_err": ev.get("phi_rel_rel_err"),
            "l2": ev.get("l2"),
            "inp": ev.get("inp"),
            "ccx_inp": ev.get("ccx_inp"),
            "output_dir": ev.get("output_dir"),
        }
        rows.append(slim)
        (work / "sphere_eval.json").write_text(json.dumps(slim, indent=2, default=str) + "\n")
        print("N={0} geo={1} bary={2} its={3} true_rel={4}".format(
            n,
            (slim["sphere_geometric"] or {}).get("rel_err"),
            (slim["sphere_barycenter"] or {}).get("rel_err"),
            (slim["solver"] or {}).get("iterations"),
            (slim["solver"] or {}).get("true_residual_rel"),
        ), flush=True)
    (out / "sphere_summary.json").write_text(json.dumps({
        "solver": solver, "rtol": rtol, "levels": rows,
        "note": "Esfera geometrica = volume tet ∩ esfera; barycenter = lote 3. Propostas, nao gate.",
    }, indent=2, default=str) + "\n")
    print("WROTE", out / "sphere_summary.json")


if __name__ == "__main__":
    main()
