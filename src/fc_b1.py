"""Benchmark B1 (Luiz): flexao pura 3D em prisma C3D4, Saint-Venant.

Geometria L x b x h = 200 x 20 x 20 mm, eixo x, E=70000 MPa, nu=0.33.
Malha estruturada: cada hexaedro em 6 tetraedros (padrao de Freudenthal,
diagonal comum n000-n111). Fixtures geradas em codigo; nao e INP da roda.

Solucao exata (Saint-Venant, kappa=M/(E I), origem no centroide da secao x=0):
  u_x = kappa*x*y
  u_y = -(kappa/2)*(x**2 + nu*(y**2 - z**2))
  u_z = -kappa*nu*y*z
  sigma_xx = M*y/I, demais nulas.

Carga: tracao sigma*n nas faces x=0 e x=L, forças consistentes T3
(tracao linear em y; ∫ N_i t dA = (A/12)*(2 t_i + t_j + t_k)).
Apoios 3-2-1 na secao x=0, sem engaste da face.

Criterios C7/C8/C2/C3 citados aqui sao PROPOSTAS de Luiz, nao gate congelado.
"""

from __future__ import annotations

import json
import math
import os
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

L = 200.0
B = 20.0
H = 20.0
E = 70000.0
NU = 0.33
M = 1.0e5  # N·mm
Izz = B * H**3 / 12.0  # 40000/3 mm^4
KAPPA = M / (E * Izz)
SIGMA_MAX = M * (H / 2.0) / Izz  # 75 MPa
U_EXACT_ENERGY = M**2 * L / (2.0 * E * Izz)
UY_MID_EXACT = -0.5 * KAPPA * (L / 2.0) ** 2  # em (L/2, 0, 0)
PHI_REL_EXACT = KAPPA * L
# Esfera monitor (C8): raio 2 mm em (L/2, h/2-3, 0) = (100, 7, 0)
SPHERE_CENTER = (L / 2.0, H / 2.0 - 3.0, 0.0)
SPHERE_RADIUS = 2.0
# sigma_xx medio analitico na esfera inteira (dentro do prisma): M*y_c/I
SIGMA_SPHERE_EXACT = M * SPHERE_CENTER[1] / Izz

# Padrao Freudenthal: 6 tets por hex, todos com a diagonal n000-n111.
# Vertices do hex em ordem: n000, n100, n110, n010, n001, n101, n111, n011.
HEX_TETS = (
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
)


def u_exact(x, y, z):
    ux = KAPPA * x * y
    uy = -0.5 * KAPPA * (x * x + NU * (y * y - z * z))
    uz = -KAPPA * NU * y * z
    return np.stack((ux, uy, uz), axis=-1)


def sigma_exact(x, y, z):
    """Retorna array (..., 3, 3) com sigma_xx = M y / I."""
    y = np.asarray(y, dtype=float)
    out = np.zeros(y.shape + (3, 3), dtype=float)
    out[..., 0, 0] = M * y / Izz
    return out


def node_index(i: int, j: int, k: int, nx: int, ny: int) -> int:
    return i + j * (nx + 1) + k * (nx + 1) * (ny + 1)


def structured_c3d4(n_thick: int, nx: int = None, nz: int = None):
    """Malha homologa: n_thick elementos na espessura (y) e em z; nx=10*n_thick (hx=hy=hz).

    Retorna coords (N,3), cells (E,4) 0-based, e metadados do grid.
    """
    ny = int(n_thick)
    if nz is None:
        nz = ny
    if nx is None:
        nx = 10 * ny  # L/H = 10, elementos aproximadamente cubicos
    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(-H / 2.0, H / 2.0, ny + 1)
    zs = np.linspace(-B / 2.0, B / 2.0, nz + 1)
    coords = []
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                coords.append((xs[i], ys[j], zs[k]))
    coords = np.asarray(coords, dtype=float)
    cells = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                verts = [
                    node_index(i, j, k, nx, ny),
                    node_index(i + 1, j, k, nx, ny),
                    node_index(i + 1, j + 1, k, nx, ny),
                    node_index(i, j + 1, k, nx, ny),
                    node_index(i, j, k + 1, nx, ny),
                    node_index(i + 1, j, k + 1, nx, ny),
                    node_index(i + 1, j + 1, k + 1, nx, ny),
                    node_index(i, j + 1, k + 1, nx, ny),
                ]
                for tet in HEX_TETS:
                    cells.append([verts[a] for a in tet])
    cells = np.asarray(cells, dtype=np.int32)
    meta = {
        "n_thick": ny,
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "pattern": "Freudenthal-6-tet, diagonal n000-n111, hex verts n000,n100,n110,n010,n001,n101,n111,n011",
        "h": H / ny,
        "n_nodes": int(coords.shape[0]),
        "n_elements": int(cells.shape[0]),
        "n_dofs": int(coords.shape[0] * 3),
    }
    return coords, cells, meta


def tet_volumes(coords: np.ndarray, cells: np.ndarray) -> np.ndarray:
    a = coords[cells[:, 0]]
    b = coords[cells[:, 1]]
    c = coords[cells[:, 2]]
    d = coords[cells[:, 3]]
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def end_face_triangles(coords: np.ndarray, cells: np.ndarray, x_value: float, atol: float = 1e-9):
    """Triangulos de contorno com os 3 nos em x=x_value, orientados com normal +x ou -x."""
    on = np.isclose(coords[:, 0], x_value, atol=atol)
    faces = []
    seen = set()
    for cell in cells:
        pts = cell.tolist()
        for tri in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            ids = [pts[i] for i in tri]
            if not all(on[i] for i in ids):
                continue
            key = tuple(sorted(ids))
            if key in seen:
                continue
            seen.add(key)
            p0, p1, p2 = coords[ids[0]], coords[ids[1]], coords[ids[2]]
            n = np.cross(p1 - p0, p2 - p0)
            # Orientar: x=L normal +x; x=0 normal -x (externa).
            want = 1.0 if x_value > L / 2.0 else -1.0
            if n[0] * want < 0:
                ids = [ids[0], ids[2], ids[1]]
            faces.append(ids)
    return faces


def consistent_face_forces(coords: np.ndarray, triangles: Sequence[Sequence[int]], sign_xx: float) -> Dict[int, np.ndarray]:
    """f_i = ∫ N_i t dA, t = (sign_xx * M y / I, 0, 0). Linear: (A/12)*(2ti+tj+tk)."""
    acc: Dict[int, np.ndarray] = {}
    for ids in triangles:
        p = coords[list(ids)]
        area_vec = np.cross(p[1] - p[0], p[2] - p[0])
        area = 0.5 * np.linalg.norm(area_vec)
        if area <= 0:
            continue
        tx = [sign_xx * M * p[k, 1] / Izz for k in range(3)]
        for a in range(3):
            fx = (area / 12.0) * (2.0 * tx[a] + tx[(a + 1) % 3] + tx[(a + 2) % 3])
            nid = int(ids[a])
            acc[nid] = acc.get(nid, np.zeros(3)) + np.array([fx, 0.0, 0.0])
    return acc


def apply_loads(coords: np.ndarray, cells: np.ndarray) -> Dict[int, np.ndarray]:
    loads: Dict[int, np.ndarray] = {}
    for x_value, sign in ((0.0, -1.0), (L, 1.0)):
        tris = end_face_triangles(coords, cells, x_value)
        part = consistent_face_forces(coords, tris, sign)
        for nid, f in part.items():
            loads[nid] = loads.get(nid, np.zeros(3)) + f
    return loads


def load_resultants(coords: np.ndarray, loads: Dict[int, np.ndarray], center=(0.0, 0.0, 0.0)):
    F = np.zeros(3)
    Mvec = np.zeros(3)
    c = np.asarray(center, dtype=float)
    for nid, f in loads.items():
        F += f
        r = coords[nid] - c
        Mvec += np.cross(r, f)
    return F, Mvec


def support_321(coords: np.ndarray, atol: float = 1e-9):
    """Nos 3-2-1 na secao x=0: centro (0,0,0), (0,H/2,0), (0,0,B/2)."""
    def find(xyz):
        dist = np.linalg.norm(coords - np.asarray(xyz, dtype=float), axis=1)
        i = int(np.argmin(dist))
        if dist[i] > atol:
            raise RuntimeError("No 3-2-1 ausente em {0} (min dist={1})".format(xyz, dist[i]))
        return i

    a = find((0.0, 0.0, 0.0))
    b = find((0.0, H / 2.0, 0.0))
    c = find((0.0, 0.0, B / 2.0))
    # (no, dof 1-3, valor) — valores exatos sao 0 nestes pontos para ux/uz;
    # uy no centro e 0. 3+2+1 GDL.
    bcs = [
        (a, 1, 0.0),
        (a, 2, 0.0),
        (a, 3, 0.0),
        (b, 1, 0.0),
        (b, 3, 0.0),
        (c, 1, 0.0),
    ]
    return bcs, {"A_000": a, "B_0_h2_0": b, "C_0_0_b2": c}


def inp_text(coords: np.ndarray, cells: np.ndarray, loads: Dict[int, np.ndarray], bcs, meta: dict) -> str:
    lines = [
        "** B1 flexao pura 3D, gerado em codigo (fc_b1.py). Nao e a roda.",
        "** padrao: {0}".format(meta["pattern"]),
        "** n_thick={0} nx={1} ny={2} nz={3} h={4}".format(
            meta["n_thick"], meta["nx"], meta["ny"], meta["nz"], meta["h"]
        ),
        "*NODE, NSET=Nall",
    ]
    for i, (x, y, z) in enumerate(coords, start=1):
        lines.append("{0}, {1:.16e}, {2:.16e}, {3:.16e}".format(i, x, y, z))
    lines.append("*ELEMENT, TYPE=C3D4, ELSET=Eall")
    for e, c in enumerate(cells, start=1):
        lines.append("{0}, {1}, {2}, {3}, {4}".format(e, c[0] + 1, c[1] + 1, c[2] + 1, c[3] + 1))
    lines.extend([
        "*MATERIAL, NAME=Mat",
        "*ELASTIC",
        "{0}, {1}".format(E, NU),
        "*SOLID SECTION, ELSET=Eall, MATERIAL=Mat",
        "*STEP, INC=100",
        "*STATIC",
        "1.0, 1.0",
        "*BOUNDARY",
    ])
    for nid0, dof, val in bcs:
        lines.append("{0}, {1}, {1}, {2:.16e}".format(nid0 + 1, dof, val))
    lines.append("*CLOAD")
    for nid0, f in sorted(loads.items()):
        if abs(f[0]) > 0.0:
            lines.append("{0}, 1, {1:.16e}".format(nid0 + 1, f[0]))
        if abs(f[1]) > 0.0:
            lines.append("{0}, 2, {1:.16e}".format(nid0 + 1, f[1]))
        if abs(f[2]) > 0.0:
            lines.append("{0}, 3, {1:.16e}".format(nid0 + 1, f[2]))
    lines.extend([
        "*OUTPUT, FREQUENCY=1",
        "*NODE FILE",
        "U, RF",
        "*EL FILE",
        "S, E",
        "*END STEP",
        "",
    ])
    return "\n".join(lines)


def face_rotation_z(coords: np.ndarray, uh: np.ndarray, x_value: float, atol: float = 1e-9) -> float:
    """Rotacao media phi_z da face: minimos quadrados u_x ≈ u0 + phi * y."""
    mask = np.isclose(coords[:, 0], x_value, atol=atol)
    y = coords[mask, 1]
    ux = uh[mask, 0]
    y = y - y.mean()
    ux = ux - ux.mean()
    den = float(np.dot(y, y))
    if den <= 0:
        return 0.0
    return float(np.dot(y, ux) / den)


def tet_quadrature_4():
    """Regra de 4 pontos (grau 2) no tetraedro de referencia."""
    a = 0.5854101966249685
    b = 0.1381966011250105
    w = 0.25
    return (
        np.array([[a, b, b, b], [b, a, b, b], [b, b, a, b], [b, b, b, a]], dtype=float),
        np.array([w, w, w, w], dtype=float),
    )


def l2_errors(coords: np.ndarray, cells: np.ndarray, uh_nodal: np.ndarray, sigma_cells: np.ndarray):
    """Erros L2(u) e L2(sigma) com quadratura de grau 2, vetorizado. sigma_cells (E,3,3)."""
    vols = np.abs(tet_volumes(coords, cells))
    p0, p1, p2, p3 = (coords[cells[:, i]] for i in range(4))
    uh0, uh1, uh2, uh3 = (uh_nodal[cells[:, i]] for i in range(4))
    bary, weights = tet_quadrature_4()
    e_u2 = n_u2 = e_s2 = n_s2 = e_sxx = n_sxx = 0.0
    for lam, w in zip(bary, weights):
        x = lam[0] * p0 + lam[1] * p1 + lam[2] * p2 + lam[3] * p3
        u_ex = u_exact(x[:, 0], x[:, 1], x[:, 2])
        u_h = lam[0] * uh0 + lam[1] * uh1 + lam[2] * uh2 + lam[3] * uh3
        du = u_h - u_ex
        e_u2 += w * float(np.sum(vols * np.sum(du * du, axis=1)))
        n_u2 += w * float(np.sum(vols * np.sum(u_ex * u_ex, axis=1)))
        s_ex = sigma_exact(x[:, 0], x[:, 1], x[:, 2])
        ds = sigma_cells - s_ex
        e_s2 += w * float(np.sum(vols * np.sum(ds * ds, axis=(1, 2))))
        n_s2 += w * float(np.sum(vols * np.sum(s_ex * s_ex, axis=(1, 2))))
        dxx = sigma_cells[:, 0, 0] - s_ex[:, 0, 0]
        e_sxx += w * float(np.sum(vols * dxx * dxx))
        n_sxx += w * float(np.sum(vols * s_ex[:, 0, 0] ** 2))
    l2u = math.sqrt(e_u2)
    l2s = math.sqrt(e_s2)
    l2sxx = math.sqrt(e_sxx)
    off = np.abs(sigma_cells[:, [0, 0, 1, 1, 2, 2], [1, 2, 0, 2, 0, 1]])
    diags_parasite = np.abs(sigma_cells[:, [1, 2], [1, 2]])  # syy, szz should be 0
    return {
        "L2_u": l2u,
        "L2_u_rel": l2u / math.sqrt(n_u2) if n_u2 > 0 else float("nan"),
        "L2_sigma": l2s,
        "L2_sigma_rel": l2s / math.sqrt(n_s2) if n_s2 > 0 else float("nan"),
        "L2_sxx": l2sxx,
        "L2_sxx_rel": l2sxx / math.sqrt(n_sxx) if n_sxx > 0 else float("nan"),
        "max_abs_sxy_sxz_syz": float(off.max()) if off.size else 0.0,
        "max_abs_syy_szz": float(diags_parasite.max()) if diags_parasite.size else 0.0,
        "max_abs_sxx": float(np.max(np.abs(sigma_cells[:, 0, 0]))),
        "volume": float(np.sum(vols)),
        "min_volume": float(np.min(tet_volumes(coords, cells))),
        "n_neg_volume": int(np.sum(tet_volumes(coords, cells) <= 0)),
    }



def mean_sigma_xx_sphere(coords: np.ndarray, cells: np.ndarray, sigma_cells: np.ndarray,
                         center=SPHERE_CENTER, radius=SPHERE_RADIUS):
    """Media volumetrica de sigma_xx nos tets cujo baricentro cai na esfera."""
    c = np.asarray(center, dtype=float)
    bary = coords[cells].mean(axis=1)
    dist = np.linalg.norm(bary - c, axis=1)
    vols = np.abs(tet_volumes(coords, cells))
    mask = dist <= radius
    if not np.any(mask):
        return {"n_tets": 0, "mean_sigma_xx": float("nan"), "volume": 0.0, "exact": SIGMA_SPHERE_EXACT, "rel_err": float("nan")}
    v = vols[mask]
    sxx = sigma_cells[mask, 0, 0]
    return {
        "n_tets": int(np.count_nonzero(mask)),
        "mean_sigma_xx": float(np.sum(sxx * v) / np.sum(v)),
        "volume": float(np.sum(v)),
        "exact": SIGMA_SPHERE_EXACT,
        "rel_err": float(np.sum(sxx * v) / np.sum(v) / SIGMA_SPHERE_EXACT - 1.0),
    }


def monitors(coords, cells, uh_nodal, sigma_cells, loads, reactions=None):
    mid = np.array([L / 2.0, 0.0, 0.0])
    dist = np.linalg.norm(coords - mid, axis=1)
    im = int(np.argmin(dist))
    uy = float(uh_nodal[im, 1])
    phi0 = face_rotation_z(coords, uh_nodal, 0.0)
    phiL = face_rotation_z(coords, uh_nodal, L)
    phi_rel = phiL - phi0
    l2 = l2_errors(coords, cells, uh_nodal, sigma_cells)
    sph = mean_sigma_xx_sphere(coords, cells, sigma_cells)
    # energia 1/2 f·u nas cargas aplicadas
    energy = 0.0
    for nid, f in loads.items():
        energy += 0.5 * float(np.dot(f, uh_nodal[nid]))
    F, Mv = load_resultants(coords, loads, (0.0, 0.0, 0.0))
    uy_err = (uy - UY_MID_EXACT) / UY_MID_EXACT if UY_MID_EXACT != 0 else float("nan")
    phi_err = (phi_rel - PHI_REL_EXACT) / PHI_REL_EXACT if PHI_REL_EXACT != 0 else float("nan")
    return {
        "uy_mid": uy,
        "uy_mid_exact": UY_MID_EXACT,
        "uy_mid_rel_err": uy_err,
        "phi_rel": phi_rel,
        "phi_rel_exact": PHI_REL_EXACT,
        "phi_rel_rel_err": phi_err,
        "stiffer_than_exact_uy": abs(uy) < abs(UY_MID_EXACT),
        "energy_half_f_u": energy,
        "energy_exact": U_EXACT_ENERGY,
        "energy_ratio": energy / U_EXACT_ENERGY if U_EXACT_ENERGY else float("nan"),
        "applied_F": F.tolist(),
        "applied_M_origin": Mv.tolist(),
        "sphere": sph,
        "l2": l2,
        "reactions": reactions,
        "constants": {
            "Izz": Izz,
            "kappa": KAPPA,
            "sigma_max": SIGMA_MAX,
            "M": M,
            "E": E,
            "nu": NU,
        },
    }


def observed_order(h1, e1, h2, e2):
    if e1 <= 0 or e2 <= 0 or h1 <= 0 or h2 <= 0 or h1 == h2:
        return float("nan")
    return math.log(e1 / e2) / math.log(h1 / h2)


def peak_rss_bytes() -> int:
    try:
        with open("/proc/self/status", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("VmHWM:"):
                    parts = line.split()
                    return int(parts[1]) * 1024
    except Exception:
        pass
    try:
        import resource
        return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024
    except Exception:
        return 0


def write_level_inp(out_dir: str, n_thick: int) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    coords, cells, meta = structured_c3d4(n_thick)
    vols = tet_volumes(coords, cells)
    if np.any(vols <= 0):
        raise RuntimeError("Volumes nao positivos no nivel {0}: min={1}".format(n_thick, vols.min()))
    loads = apply_loads(coords, cells)
    bcs, support = support_321(coords)
    text = inp_text(coords, cells, loads, bcs, meta)
    path = os.path.join(out_dir, "b1_n{0}.inp".format(n_thick))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    F, Mv = load_resultants(coords, loads)
    info = {
        "inp": path,
        "meta": meta,
        "n_loads": len(loads),
        "applied_F": F.tolist(),
        "applied_M_origin": Mv.tolist(),
        "support": support,
        "min_volume": float(vols.min()),
        "sum_volume": float(vols.sum()),
        "volume_exact": L * B * H,
    }
    with open(os.path.join(out_dir, "b1_n{0}_mesh.json".format(n_thick)), "w", encoding="utf-8") as handle:
        json.dump(info, handle, indent=2)
        handle.write("\n")
    return {"coords": coords, "cells": cells, "loads": loads, "bcs": bcs, "meta": meta, "info": info, "inp": path}
