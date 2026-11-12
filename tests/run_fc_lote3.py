#!/usr/bin/env python3
"""Lote 3: B1 (2/4/8/16) + serie de recursos LU, isolado.

Nao toca o INP da roda. Criterios C7/C8 sao PROPOSTAS, nao gate aprovado.
"""

from __future__ import annotations

import json
import os
import sys
import time
import traceback

sys.path.insert(0, os.environ.get("FC_SRC", "/app/src"))
if os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")) not in sys.path:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

import numpy as np
from fc_b1 import (
    B,
    E,
    H,
    L,
    M,
    PHI_REL_EXACT,
    SIGMA_SPHERE_EXACT,
    SPHERE_CENTER,
    SPHERE_RADIUS,
    UY_MID_EXACT,
    U_EXACT_ENERGY,
    apply_loads,
    face_rotation_z,
    l2_errors,
    mean_sigma_xx_sphere,
    observed_order,
    peak_rss_bytes,
    structured_c3d4,
    support_321,
    tet_volumes,
    write_level_inp,
)


# Escalas declaradas para C2/C3 neste benchmark (PROPOSTA de aplicacao, nao congelada).
# Carga liquida e nula; F_ref_B1 = M / (H/2) = 1e4 N; L_ref_B1 = L.
F_REF_B1 = M / (H / 2.0)
L_REF_B1 = L


def _geom_cells(domain):
    geom = domain.geometry
    gdmap = geom.dofmaps[0] if hasattr(geom, "dofmaps") else geom.dofmap
    ncells = int(domain.topology.index_map(domain.topology.dim).size_local)

    def links(c):
        try:
            return [int(i) for i in gdmap[c]]
        except Exception:
            return [int(i) for i in gdmap.links(c)]

    return np.asarray([links(c) for c in range(ncells)], dtype=np.int32)


def _loads_on_geom(coords_g, coords_inp, loads_inp, atol=1e-9):
    mapped = {}
    for nid, f in loads_inp.items():
        dist = np.linalg.norm(coords_g - coords_inp[nid], axis=1)
        hits = np.flatnonzero(dist <= atol)
        if hits.size != 1:
            raise RuntimeError("Carga no INP {0} nao mapeou 1:1 na geometria (hits={1})".format(nid, hits.size))
        mapped[int(hits[0])] = f
    return mapped


def evaluate(domain, uh, sigma_cells, coords_inp, loads_inp, extra):
    coords = np.asarray(domain.geometry.x, dtype=float)
    cells = _geom_cells(domain)
    uarr = np.asarray(uh.x.array, dtype=float).reshape((-1, 3))
    if uarr.shape[0] != coords.shape[0]:
        raise RuntimeError("Layout u/geometry inesperado: u={0} geom={1}".format(uarr.shape, coords.shape))
    loads_g = _loads_on_geom(coords, coords_inp, loads_inp)
    mid = np.array([L / 2.0, 0.0, 0.0])
    im = int(np.argmin(np.linalg.norm(coords - mid, axis=1)))
    if np.linalg.norm(coords[im] - mid) > 1e-9:
        raise RuntimeError("No (L/2,0,0) ausente: {0}".format(coords[im]))
    uy = float(uarr[im, 1])
    phi_rel = face_rotation_z(coords, uarr, L) - face_rotation_z(coords, uarr, 0.0)
    l2 = l2_errors(coords, cells, uarr, sigma_cells)
    sph_bary = mean_sigma_xx_sphere(coords, cells, sigma_cells, method="barycenter")
    sph_geo = mean_sigma_xx_sphere(coords, cells, sigma_cells, method="geometric")
    sph = sph_geo
    energy = 0.0
    for gid, f in loads_g.items():
        energy += 0.5 * float(np.dot(f, uarr[gid]))
    rx = extra.get("reactions") or {}
    sum_rf = rx.get("sum_R_plus_F") or {}
    mom = rx.get("moment_R_plus_F_origin") or {}
    def c2(key, val):
        return abs(float(val)) <= 1e-6 * F_REF_B1
    def c3(key, val):
        return abs(float(val)) <= 1e-6 * F_REF_B1 * L_REF_B1
    return {
        "uy_mid": uy,
        "uy_mid_exact": UY_MID_EXACT,
        "uy_mid_rel_err": (uy - UY_MID_EXACT) / UY_MID_EXACT,
        "stiffer_than_exact_uy": abs(uy) < abs(UY_MID_EXACT),
        "phi_rel": phi_rel,
        "phi_rel_exact": PHI_REL_EXACT,
        "phi_rel_rel_err": (phi_rel - PHI_REL_EXACT) / PHI_REL_EXACT,
        "energy_half_f_u": energy,
        "energy_exact": U_EXACT_ENERGY,
        "energy_ratio": energy / U_EXACT_ENERGY,
        "energy_le_exact": energy <= U_EXACT_ENERGY * (1.0 + 1e-6),
        "l2": l2,
        "sphere": sph,
        "sphere_geometric": sph_geo,
        "sphere_barycenter": sph_bary,
        "reactions": {
            "sum_R_plus_F": sum_rf,
            "moment_R_plus_F_origin": mom,
            "residual_free_dofs_max_abs": rx.get("residual_free_dofs_max_abs"),
            "totals": rx.get("totals"),
            "applied_totals": rx.get("applied_totals"),
        },
        "proposed_C2": {k: {"value": sum_rf.get(k), "pass_proposed": c2(k, sum_rf.get(k, 0.0))} for k in ("Fx", "Fy", "Fz")},
        "proposed_C3_origin": {k: {"value": mom.get(k), "pass_proposed": c3(k, mom.get(k, 0.0))} for k in ("Mx", "My", "Mz")},
        "F_ref_B1": F_REF_B1,
        "L_ref_B1": L_REF_B1,
        "note": "C2/C3/C7/C8 sao criterios PROPOSTOS, nao congelados nem Gate V aprovado.",
    }


def run_level(n_thick: int, work: str, debug: int = 0, solver: str = "lu", rtol: float = 1e-10):
    from run_from_inp import run_from_inp

    pack = write_level_inp(work, n_thick)
    out = os.path.join(work, "out_n{0}".format(n_thick))
    os.makedirs(out, exist_ok=True)
    rss0 = peak_rss_bytes()
    t0 = time.perf_counter()
    kw = dict(output_dir=out, debug=debug, profile="fc", units="mm-N-MPa-t-s")
    if solver and solver != "lu":
        kw["solver"] = solver
        kw["rtol"] = rtol
    results = run_from_inp(pack["inp"], **kw)
    elapsed = time.perf_counter() - t0
    rss1 = peak_rss_bytes()
    uh, extra = results[0]
    ev = evaluate(
        uh.function_space.mesh,
        uh,
        extra["stress_cells"],
        pack["coords"],
        pack["loads"],
        extra,
    )
    ev["meta"] = pack["meta"]
    ev["elapsed_s"] = elapsed
    ev["rss_before_bytes"] = rss0
    ev["rss_peak_bytes"] = rss1
    ev["rss_peak_MiB"] = rss1 / (1024 * 1024)
    ev["output_dir"] = out
    ev["inp"] = pack["inp"]
    ev["mesh_info"] = pack["info"]
    ev["solver"] = extra.get("solver_info") or extra.get("solver")
    return ev


def series_table(rows):
    rows = sorted(rows, key=lambda r: r["meta"]["h"], reverse=True)
    out = []
    for i, r in enumerate(rows):
        item = {
            "n_thick": r["meta"]["n_thick"],
            "h": r["meta"]["h"],
            "n_nodes": r["meta"]["n_nodes"],
            "n_elements": r["meta"]["n_elements"],
            "n_dofs": r["meta"]["n_dofs"],
            "uy_mid_rel_err": r["uy_mid_rel_err"],
            "phi_rel_rel_err": r["phi_rel_rel_err"],
            "L2_u_rel": r["l2"]["L2_u_rel"],
            "L2_sigma_rel": r["l2"]["L2_sigma_rel"],
            "L2_sxx_rel": r["l2"].get("L2_sxx_rel"),
            "max_abs_sxy_sxz_syz": r["l2"].get("max_abs_sxy_sxz_syz"),
            "max_abs_syy_szz": r["l2"].get("max_abs_syy_szz"),
            "sphere_sigma_rel_err": r["sphere"].get("rel_err"),
            "energy_ratio": r["energy_ratio"],
            "stiffer": r["stiffer_than_exact_uy"],
            "elapsed_s": r["elapsed_s"],
            "rss_peak_MiB": r["rss_peak_MiB"],
        }
        if i + 1 < len(rows):
            a, b = r, rows[i + 1]
            item["order_to_next"] = {
                "L2_u": observed_order(a["meta"]["h"], a["l2"]["L2_u_rel"], b["meta"]["h"], b["l2"]["L2_u_rel"]),
                "L2_sigma": observed_order(a["meta"]["h"], a["l2"]["L2_sigma_rel"], b["meta"]["h"], b["l2"]["L2_sigma_rel"]),
                "L2_sxx": observed_order(a["meta"]["h"], a["l2"].get("L2_sxx_rel") or 1e300, b["meta"]["h"], b["l2"].get("L2_sxx_rel") or 1e300),
                "uy_mid": observed_order(a["meta"]["h"], abs(a["uy_mid_rel_err"]), b["meta"]["h"], abs(b["uy_mid_rel_err"])),
                "phi_rel": observed_order(a["meta"]["h"], abs(a["phi_rel_rel_err"]), b["meta"]["h"], abs(b["phi_rel_rel_err"])),
                "sphere": observed_order(
                    a["meta"]["h"], abs(a["sphere"].get("rel_err") or 1e300),
                    b["meta"]["h"], abs(b["sphere"].get("rel_err") or 1e300),
                ),
            }
        out.append(item)
    return out


def proposed_c7_c8(rows):
    """Avalia propostas C7/C8 sem declara-las aprovadas como gate."""
    by_n = {r["meta"]["n_thick"]: r for r in rows}
    finest = by_n.get(max(by_n))
    c7_u = abs(finest["uy_mid_rel_err"]) <= 0.02 if finest else False
    c7_phi = abs(finest["phi_rel_rel_err"]) <= 0.02 if finest else False
    # Luiz: erro relativo do monitor <= 2% na malha mais fina; ordem L2(u) em [1.6, 2.4] nos 3 mais finos
    fine_sorted = sorted(rows, key=lambda r: r["meta"]["h"])[:3]
    fine_sorted = list(reversed(fine_sorted))  # coarse to fine among the finest 3
    orders_u = []
    orders_s = []
    tab = series_table(rows)
    for item in tab:
        if item["n_thick"] in {r["meta"]["n_thick"] for r in fine_sorted} and "order_to_next" in item:
            if item["n_thick"] != min(r["meta"]["n_thick"] for r in fine_sorted):
                # only intervals among the 3 finest: 4-8 and 8-16 if all present
                ns = sorted(r["meta"]["n_thick"] for r in fine_sorted)
                if item["n_thick"] in ns[:-1]:
                    orders_u.append(item["order_to_next"]["L2_u"])
                    orders_s.append(item["order_to_next"]["L2_sigma"])
    def in_band(vals, lo, hi):
        return bool(vals) and all(lo <= v <= hi for v in vals if v == v)
    sph_err = (finest.get("sphere") or {}).get("rel_err") if finest else None
    c8_sph = False
    if sph_err is not None and sph_err == sph_err:
        c8_sph = abs(sph_err) <= 0.05
    return {
        "finest_n": finest["meta"]["n_thick"] if finest else None,
        "proposed_C7_uy_le_2pct": c7_u,
        "proposed_C7_phi_le_2pct": c7_phi,
        "proposed_C7_L2u_order_band_1.6_2.4": in_band(orders_u, 1.6, 2.4),
        "orders_L2_u_finest": orders_u,
        "proposed_C8_sphere_le_5pct": c8_sph,
        "proposed_C8_L2s_order_band_0.7_1.3": in_band(orders_s, 0.7, 1.3),
        "orders_L2_sigma_finest": orders_s,
        "gate_v_approved": False,
        "note": "Propostas apenas. Gate V nao aprovado: criterios nao estao congelados.",
    }


def extrapolate_ram(rows, target_dofs=376185, safety=3.0):
    """Ajuste log-log rss vs dofs (N**p). Nao e medicao da roda."""
    pts = [(r["meta"]["n_dofs"], r["rss_peak_MiB"]) for r in rows if r.get("rss_peak_MiB", 0) > 1]
    if len(pts) < 2:
        return {"ok": False, "reason": "menos de 2 pontos medidos"}
    logN = np.log([p[0] for p in pts])
    logM = np.log([p[1] for p in pts])
    # reta nos dois (ou mais) ultimos pontos
    use = pts[-3:] if len(pts) >= 3 else pts
    x = np.log([p[0] for p in use])
    y = np.log([p[1] for p in use])
    A = np.vstack([x, np.ones_like(x)]).T
    p, b = np.linalg.lstsq(A, y, rcond=None)[0]
    pred = float(np.exp(p * np.log(target_dofs) + b))
    return {
        "ok": True,
        "points_MiB": pts,
        "power": float(p),
        "pred_MiB_at_target": pred,
        "target_dofs": target_dofs,
        "safety_factor": safety,
        "pred_with_safety_MiB": pred * safety,
        "model": "rss = exp(b) * N**p, ajuste nos pontos mais finos",
    }


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    work = "/tmp/fc_lote3"
    levels = [2, 4, 8, 16]
    if "--work" in argv:
        work = argv[argv.index("--work") + 1]
    if "--levels" in argv:
        levels = [int(x) for x in argv[argv.index("--levels") + 1].split(",") if x]
    os.makedirs(work, exist_ok=True)
    print("B1 work=", work, "levels=", levels, "VmHWM0_MiB=", peak_rss_bytes() / 1024 / 1024)
    rows = []
    errors = []
    for n in levels:
        print("==== B1 n_thick={0} ====".format(n), flush=True)
        try:
            ev = run_level(n, work, debug=0)
            rows.append(ev)
            print(json.dumps({
                "n": n,
                "dofs": ev["meta"]["n_dofs"],
                "elems": ev["meta"]["n_elements"],
                "uy_rel": ev["uy_mid_rel_err"],
                "phi_rel": ev["phi_rel_rel_err"],
                "L2u": ev["l2"]["L2_u_rel"],
                "L2s": ev["l2"]["L2_sigma_rel"],
                "L2sxx": ev["l2"].get("L2_sxx_rel"),
                "tau_max": ev["l2"].get("max_abs_sxy_sxz_syz"),
                "syy_szz_max": ev["l2"].get("max_abs_syy_szz"),
                "sphere": ev["sphere"].get("rel_err"),
                "energy_ratio": ev["energy_ratio"],
                "stiffer": ev["stiffer_than_exact_uy"],
                "s": ev["elapsed_s"],
                "MiB": ev["rss_peak_MiB"],
                "C2": ev["proposed_C2"],
                "C3": ev["proposed_C3_origin"],
            }, indent=2), flush=True)
            with open(os.path.join(work, "b1_n{0}_eval.json".format(n)), "w", encoding="utf-8") as handle:
                json.dump(ev, handle, indent=2, default=str)
                handle.write("\n")
        except Exception as exc:
            traceback.print_exc()
            errors.append({"n_thick": n, "error": "{0}: {1}".format(type(exc).__name__, exc)})
            break
    table = series_table(rows)
    proposal = proposed_c7_c8(rows) if rows else {}
    ram = extrapolate_ram(rows) if rows else {"ok": False}
    payload = {
        "benchmark": "B1-pure-bending-C3D4",
        "criteria": "PROPOSTAS Luiz 11/11/2026 17:00, nao congeladas, Gate V nao aprovado",
        "kappa_note": (
            "Luiz citou kappa=1.07e-6 /mm; as formulas kappa=M/(E I) com M=1e5 e "
            "sigma_max=75 MPa dao kappa=1.0714e-4 /mm. Usamos as formulas, nao o typo."
        ),
        "levels": table,
        "proposed_C7_C8": proposal,
        "ram_extrapolation": ram,
        "errors": errors,
        "constants": {
            "L": L, "B": B, "H": H, "E": E, "M": M,
            "uy_mid_exact": UY_MID_EXACT,
            "phi_rel_exact": PHI_REL_EXACT,
            "sigma_sphere_exact": SIGMA_SPHERE_EXACT,
            "energy_exact": U_EXACT_ENERGY,
            "F_ref_B1": F_REF_B1,
            "L_ref_B1": L_REF_B1,
        },
    }
    outp = os.path.join(work, "b1_summary.json")
    with open(outp, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, default=str)
        handle.write("\n")
    print("SUMMARY", outp, flush=True)
    print(json.dumps(payload, indent=2, default=str), flush=True)
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
