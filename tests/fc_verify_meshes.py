"""Malhas sinteticas C3D4 para T0 (hidrostatico) e T2 (patch test).

Blocos pequenos, alguns centos de tetraedros lineares, com nos interiores
perturbados (semente fixa). Nao toca o INP da roda e nao e caso definitivo.
O perfil FC rejeita *DLOAD: pressao e tracao entram como *CLOAD consistentes
(F_i = A t / 3 no triangulo linear).
"""

from __future__ import annotations

import itertools
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

E = 70000.0
NU = 0.33
# Cubo 10 mm, malha 4 hex por aresta -> 5^3 nos, 64 hex, 384 tets.
L = 10.0
N_HEX = 4
PERTURB_FRAC = 0.22
SEED = 20261113
PRESSURE = 10.0  # MPa (N/mm^2)

# Campo linear do patch: cisalhamentos nao nulos e distintos nas 3 componentes
# (requisito T2 / mapeamento C5). Tambem ha rotacao rigida (parte anti-simetrica).
A_VEC = np.array([0.11, -0.07, 0.05], dtype=float)
B_MAT = np.array(
    [
        [1.2e-3, 4.0e-3, 6.0e-3],
        [2.0e-3, -3.1e-3, 1.0e-2],
        [3.0e-3, 1.0e-3, 2.4e-3],
    ],
    dtype=float,
)

# Freudenthal, mesma convencao de fc_b1: diagonal n000-n111.
HEX_TETS = (
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
)


def node_index(i: int, j: int, k: int, n: int) -> int:
    return i + j * (n + 1) + k * (n + 1) * (n + 1)


def structured_cube(n_hex: int = N_HEX, length: float = L):
    xs = np.linspace(0.0, length, n_hex + 1)
    coords = []
    for k in range(n_hex + 1):
        for j in range(n_hex + 1):
            for i in range(n_hex + 1):
                coords.append((xs[i], xs[j], xs[k]))
    coords = np.asarray(coords, dtype=float)
    cells = []
    for k in range(n_hex):
        for j in range(n_hex):
            for i in range(n_hex):
                verts = [
                    node_index(i, j, k, n_hex),
                    node_index(i + 1, j, k, n_hex),
                    node_index(i + 1, j + 1, k, n_hex),
                    node_index(i, j + 1, k, n_hex),
                    node_index(i, j, k + 1, n_hex),
                    node_index(i + 1, j, k + 1, n_hex),
                    node_index(i + 1, j + 1, k + 1, n_hex),
                    node_index(i, j + 1, k + 1, n_hex),
                ]
                for tet in HEX_TETS:
                    cells.append([verts[a] for a in tet])
    cells = np.asarray(cells, dtype=np.int32)
    meta = {
        "n_hex": n_hex,
        "length": length,
        "h": length / n_hex,
        "n_nodes": int(coords.shape[0]),
        "n_elements": int(cells.shape[0]),
        "pattern": "Freudenthal-6-tet, cubo [0,L]^3, diagonal n000-n111",
    }
    return coords, cells, meta


def tet_volumes(coords: np.ndarray, cells: np.ndarray) -> np.ndarray:
    a = coords[cells[:, 0]]
    b = coords[cells[:, 1]]
    c = coords[cells[:, 2]]
    d = coords[cells[:, 3]]
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def classify_nodes(coords: np.ndarray, length: float = L, atol: float = 1e-12):
    """Retorna mascaras: interior, e faces no plano original do cubo."""
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    on = {
        "x0": np.isclose(x, 0.0, atol=atol),
        "xL": np.isclose(x, length, atol=atol),
        "y0": np.isclose(y, 0.0, atol=atol),
        "yL": np.isclose(y, length, atol=atol),
        "z0": np.isclose(z, 0.0, atol=atol),
        "zL": np.isclose(z, length, atol=atol),
    }
    boundary = on["x0"] | on["xL"] | on["y0"] | on["yL"] | on["z0"] | on["zL"]
    return on, boundary, ~boundary


def perturb_nodes(
    coords: np.ndarray,
    length: float = L,
    frac: float = PERTURB_FRAC,
    seed: int = SEED,
    keep_faces_planar: bool = False,
):
    """Perturba nos para distorcer os tets, semente fixa.

    keep_faces_planar=True (hidrostatico com simetria): so nos interiores,
    e nos de face se movem so no plano da face. Cantos ficam fixos.
    """
    rng = np.random.default_rng(seed)
    h = length / N_HEX
    amp = frac * h
    out = coords.copy()
    on, boundary, interior = classify_nodes(coords, length)
    delta = rng.uniform(-amp, amp, size=out.shape)
    if keep_faces_planar:
        out[interior] += delta[interior]
        # Face, mas nao aresta de canto: desliza no plano.
        for mask, freeze in (
            (on["x0"] | on["xL"], 0),
            (on["y0"] | on["yL"], 1),
            (on["z0"] | on["zL"], 2),
        ):
            only = mask & ~interior
            move = delta[only]
            move[:, freeze] = 0.0
            # Arestas/cantos (2+ planos): nao mover, para nao sair do cubo.
            nplanes = (
                on["x0"][only].astype(int)
                + on["xL"][only].astype(int)
                + on["y0"][only].astype(int)
                + on["yL"][only].astype(int)
                + on["z0"][only].astype(int)
                + on["zL"][only].astype(int)
            )
            move[nplanes >= 2] = 0.0
            out[only] += move
    else:
        # Patch: distorce interiores e tambem nos de face (nao cantos).
        nplanes = (
            on["x0"].astype(int)
            + on["xL"].astype(int)
            + on["y0"].astype(int)
            + on["yL"].astype(int)
            + on["z0"].astype(int)
            + on["zL"].astype(int)
        )
        movable = (nplanes < 3)
        out[movable] += delta[movable]
        # Mantem o bloco dentro de uma caixa folgada; cantos (nplanes==3) fixos.
    return out


def topological_boundary_mask(coords: np.ndarray, cells: np.ndarray) -> np.ndarray:
    mask = np.zeros(coords.shape[0], dtype=bool)
    for ids, _area, _nhat in boundary_triangles(coords, cells):
        for nid in ids:
            mask[int(nid)] = True
    return mask


def make_mesh(kind: str, n_hex: int = N_HEX, length: float = L):
    coords0, cells, meta = structured_cube(n_hex, length)
    on0, boundary_planes0, interior0 = classify_nodes(coords0, length)
    keep_planar = kind == "hydrostatic"
    coords = perturb_nodes(coords0, length=length, keep_faces_planar=keep_planar)
    vols = tet_volumes(coords, cells)
    if np.any(vols <= 0):
        raise RuntimeError(
            "Volume nao positivo apos perturbacao ({0}): min={1}".format(kind, float(vols.min()))
        )
    boundary_topo = topological_boundary_mask(coords, cells)
    meta = dict(meta)
    meta["kind"] = kind
    meta["seed"] = SEED
    meta["perturb_frac"] = PERTURB_FRAC
    meta["keep_faces_planar"] = keep_planar
    meta["min_volume"] = float(vols.min())
    meta["sum_volume"] = float(vols.sum())
    meta["volume_bbox"] = float(length ** 3)
    meta["n_boundary_topo"] = int(np.count_nonzero(boundary_topo))
    meta["n_interior_topo"] = int(np.count_nonzero(~boundary_topo))
    return {
        "coords": coords,
        "coords0": coords0,
        "cells": cells,
        "meta": meta,
        "vols": vols,
        "on0": on0,
        "boundary_planes0": boundary_planes0,
        "boundary_topo": boundary_topo,
    }


def even_permutations_4() -> List[Tuple[int, int, int, int]]:
    out = []
    for perm in itertools.permutations(range(4)):
        inversions = sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j])
        if inversions % 2 == 0:
            out.append(perm)
    if len(out) != 12:
        raise RuntimeError("Esperadas 12 permutacoes pares, obtidas {0}".format(len(out)))
    return out


def permute_cells(cells: np.ndarray, perm: Sequence[int]) -> np.ndarray:
    return np.ascontiguousarray(cells[:, list(perm)], dtype=np.int32)


def hooke_sigma(eps: np.ndarray, e: float = E, nu: float = NU) -> np.ndarray:
    lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = e / (2.0 * (1.0 + nu))
    eye = np.eye(3, dtype=float)
    if eps.ndim == 2:
        return lam * np.trace(eps) * eye + 2.0 * mu * eps
    traces = np.trace(eps, axis1=-2, axis2=-1)[..., None, None]
    return lam * traces * eye + 2.0 * mu * eps


def patch_strain() -> np.ndarray:
    return 0.5 * (B_MAT + B_MAT.T)


def patch_sigma() -> np.ndarray:
    return hooke_sigma(patch_strain())


def patch_u(xyz: np.ndarray) -> np.ndarray:
    x = np.asarray(xyz, dtype=float)
    if x.ndim == 1:
        return A_VEC + B_MAT @ x
    return A_VEC + x @ B_MAT.T


def hydrostatic_alpha(p: float = PRESSURE, e: float = E, nu: float = NU) -> float:
    """eps_ii = -p (1-2 nu)/E; u = alpha * x com vinculo de simetria na origem."""
    return -p * (1.0 - 2.0 * nu) / e


def hydrostatic_u(xyz: np.ndarray, p: float = PRESSURE) -> np.ndarray:
    return hydrostatic_alpha(p) * np.asarray(xyz, dtype=float)


def hydrostatic_eps_vol(p: float = PRESSURE, e: float = E, nu: float = NU) -> float:
    return -3.0 * p * (1.0 - 2.0 * nu) / e


def boundary_triangles(coords: np.ndarray, cells: np.ndarray):
    """Faces de contorno com normal unitaria externa e area."""
    from collections import defaultdict

    occ: Dict[Tuple[int, ...], List[Tuple[int, List[int]]]] = defaultdict(list)
    faces = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    for e, cell in enumerate(cells):
        for face in faces:
            ids = [int(cell[i]) for i in face]
            occ[tuple(sorted(ids))].append((e, ids))
    tris = []
    for _key, items in occ.items():
        if len(items) != 1:
            continue
        e, ids = items[0]
        pts = coords[ids]
        normal = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        cell = cells[e]
        centroid = coords[cell].mean(axis=0)
        face_c = pts.mean(axis=0)
        if np.dot(normal, centroid - face_c) > 0.0:
            ids = [ids[0], ids[2], ids[1]]
            normal = -normal
        nlen = float(np.linalg.norm(normal))
        if nlen <= 0.0:
            continue
        area = 0.5 * nlen
        nhat = normal / nlen
        tris.append((ids, area, nhat))
    return tris


def pressure_cloads(coords: np.ndarray, cells: np.ndarray, p: float) -> Dict[int, np.ndarray]:
    """Tracao t = -p n (n externa) com forcas consistentes F_i = A t / 3."""
    loads: Dict[int, np.ndarray] = {}
    for ids, area, nhat in boundary_triangles(coords, cells):
        force = (area / 3.0) * (-p * nhat)
        for nid in ids:
            loads[nid] = loads.get(nid, np.zeros(3, dtype=float)) + force
    return loads


def traction_cloads(
    coords: np.ndarray,
    cells: np.ndarray,
    sigma: np.ndarray,
    on_faces: dict,
    faces: Iterable[str] = ("xL", "yL", "zL"),
) -> Dict[int, np.ndarray]:
    """Tracao t = sigma n nas faces pedidas (n externa), F_i = A t / 3.

    A face e identificada pelos nos originais do cubo (on_faces), nao pelos
    planos apos a perturbacao.
    """
    want = set(faces)
    loads: Dict[int, np.ndarray] = {}
    for ids, area, nhat in boundary_triangles(coords, cells):
        labels = []
        for name in ("x0", "xL", "y0", "yL", "z0", "zL"):
            if all(bool(on_faces[name][int(nid)]) for nid in ids):
                labels.append(name)
        if not any(lab in want for lab in labels):
            continue
        t = sigma @ nhat
        force = (area / 3.0) * t
        for nid in ids:
            loads[nid] = loads.get(nid, np.zeros(3, dtype=float)) + force
    return loads


def symmetry_bcs(coords: np.ndarray, length: float = L, atol: float = 1e-9):
    """ux=0 em x=0, uy=0 em y=0, uz=0 em z=0. (no 0-based, dof 1-3, valor)."""
    bcs = []
    for i, (x, y, z) in enumerate(coords):
        if abs(x - 0.0) <= atol:
            bcs.append((i, 1, 0.0))
        if abs(y - 0.0) <= atol:
            bcs.append((i, 2, 0.0))
        if abs(z - 0.0) <= atol:
            bcs.append((i, 3, 0.0))
    return bcs


def dirichlet_all_boundary(coords: np.ndarray, boundary_mask: np.ndarray):
    bcs = []
    for i in np.nonzero(boundary_mask)[0]:
        u = patch_u(coords[i])
        bcs.append((int(i), 1, float(u[0])))
        bcs.append((int(i), 2, float(u[1])))
        bcs.append((int(i), 3, float(u[2])))
    return bcs, boundary_mask


def dirichlet_partial_faces(
    coords: np.ndarray,
    on_faces: dict,
    faces: Iterable[str] = ("x0", "y0", "z0"),
):
    """Usa as mascaras de face do cubo ANTES da perturbacao (ids de no fixos)."""
    mask = np.zeros(coords.shape[0], dtype=bool)
    for name in faces:
        mask |= on_faces[name]
    bcs = []
    for i in np.nonzero(mask)[0]:
        u = patch_u(coords[i])
        bcs.append((int(i), 1, float(u[0])))
        bcs.append((int(i), 2, float(u[1])))
        bcs.append((int(i), 3, float(u[2])))
    return bcs, mask


def inp_text(
    coords: np.ndarray,
    cells: np.ndarray,
    bcs,
    loads: Dict[int, np.ndarray],
    meta: dict,
    heading: str,
    e: float = E,
    nu: float = NU,
) -> str:
    lines = [
        "** {0}".format(heading),
        "** Fixture de verificacao C3D4 do Fenics-EDSS. Nao e a roda.",
        "** kind={0} n_hex={1} L={2} n_nodes={3} n_elements={4} seed={5}".format(
            meta.get("kind"), meta.get("n_hex"), meta.get("length"),
            meta.get("n_nodes"), meta.get("n_elements"), meta.get("seed"),
        ),
        "** {0}".format(meta.get("pattern")),
        # *HEADING e *EL PRINT nao entram no .inp do perfil FC (fora da whitelist).
        # As copias compactas .ccx.inp acrescentam esses cartoes para o CCX.
        "*NODE, NSET=Nall",
    ]
    for i, (x, y, z) in enumerate(coords, start=1):
        lines.append("{0}, {1:.16e}, {2:.16e}, {3:.16e}".format(i, x, y, z))
    lines.append("*ELEMENT, TYPE=C3D4, ELSET=Eall")
    for e_id, c in enumerate(cells, start=1):
        lines.append("{0}, {1}, {2}, {3}, {4}".format(e_id, c[0] + 1, c[1] + 1, c[2] + 1, c[3] + 1))
    lines.extend([
        "*MATERIAL, NAME=Mat",
        "*ELASTIC",
        "{0}, {1}".format(e, nu),
        "*SOLID SECTION, ELSET=Eall, MATERIAL=Mat",
        "*STEP, INC=100",
        "*STATIC",
        "1.0, 1.0",
        "*BOUNDARY",
    ])
    for nid0, dof, val in bcs:
        lines.append("{0}, {1}, {1}, {2:.16e}".format(nid0 + 1, dof, val))
    if loads:
        lines.append("*CLOAD")
        for nid0, f in sorted(loads.items()):
            for comp, val in enumerate(f, start=1):
                if val != 0.0:
                    lines.append("{0}, {1}, {2:.16e}".format(nid0 + 1, comp, val))
    lines.extend([
        "*OUTPUT, FREQUENCY=1",
        "*NODE FILE",
        "U, RF",
        "*EL FILE",
        "S, E",
        "*NODE PRINT, NSET=Nall",
        "U, RF",
        "*END STEP",
        "",
    ])
    return "\n".join(lines)


def ccx_decorate(text: str) -> str:
    """Acrescenta *HEADING e *EL PRINT (aceitaveis no CCX, rejeitados no perfil FC)."""
    lines = text.splitlines()
    out = []
    inserted_heading = False
    for line in lines:
        if not inserted_heading and line.startswith("*NODE"):
            out.append("*HEADING")
            out.append("Fenics-EDSS T0/T2 fixture compacta para CCX 2.23")
            inserted_heading = True
        if line.strip() == "*END STEP":
            out.append("*EL PRINT, ELSET=Eall, GLOBAL=YES")
            out.append("S, COORD")
        out.append(line)
    return "\n".join(out) + "\n"


def ccx_compact_text(text: str) -> str:
    """Encurta linhas de dados numericos para o leitor Fortran do CCX 2.23."""
    out = []
    for line in text.splitlines():
        raw = line.strip()
        if not raw or raw.startswith("*") or raw.startswith("**"):
            out.append(line)
            continue
        if not (raw[0].isdigit() or raw[0] in "+-"):
            out.append(line)
            continue
        parts = [p.strip() for p in raw.split(",")]
        compact = []
        for tok in parts:
            try:
                v = float(tok)
            except ValueError:
                compact.append(tok)
                continue
            if v == 0.0:
                compact.append("0")
            elif 1e-4 <= abs(v) < 1e7:
                compact.append("{0:.12g}".format(v))
            else:
                compact.append("{0:.8E}".format(v))
        out.append(", ".join(compact))
    return "\n".join(out) + "\n"


def write_inp(path: str, text: str) -> str:
    from pathlib import Path

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    return str(p)


def build_patch_dirichlet():
    mesh = make_mesh("patch")
    coords, cells, meta, vols = mesh["coords"], mesh["cells"], mesh["meta"], mesh["vols"]
    bcs, boundary = dirichlet_all_boundary(coords, mesh["boundary_topo"])
    text = inp_text(
        coords, cells, bcs, {}, meta,
        heading="T2 patch test C3D4, Dirichlet em todo o contorno, u=a+B x",
    )
    return {
        "coords": coords, "cells": cells, "meta": meta, "vols": vols,
        "bcs": bcs, "loads": {}, "boundary": boundary, "text": text,
    }


def build_patch_traction():
    mesh = make_mesh("patch")
    coords, cells, meta, vols = mesh["coords"], mesh["cells"], mesh["meta"], mesh["vols"]
    bcs, dmask = dirichlet_partial_faces(coords, mesh["on0"])
    loads = traction_cloads(coords, cells, patch_sigma(), mesh["on0"], faces=("xL", "yL", "zL"))
    text = inp_text(
        coords, cells, bcs, loads, meta,
        heading="T2 patch test C3D4, Dirichlet parcial (x=0,y=0,z=0) + tracao consistente",
    )
    return {
        "coords": coords, "cells": cells, "meta": meta, "vols": vols,
        "bcs": bcs, "loads": loads, "dirichlet": dmask, "text": text,
    }


def build_hydrostatic(p: float = PRESSURE):
    mesh = make_mesh("hydrostatic")
    coords, cells, meta, vols = mesh["coords"], mesh["cells"], mesh["meta"], mesh["vols"]
    # Simetria nos planos originais (faces permaneceram planas).
    bcs = symmetry_bcs(coords)
    loads = pressure_cloads(coords, cells, p)
    meta = dict(meta)
    meta["pressure"] = p
    text = inp_text(
        coords, cells, bcs, loads, meta,
        heading="T0 hidrostatico C3D4, pressao p={0} MPa em todas as faces, simetria 0".format(p),
    )
    return {
        "coords": coords, "cells": cells, "meta": meta, "vols": vols,
        "bcs": bcs, "loads": loads, "text": text, "pressure": p,
    }


if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    dest = os.path.join(here, "fixtures")
    os.makedirs(dest, exist_ok=True)
    for name, pack in (
        ("t2_patch_block.inp", build_patch_dirichlet()),
        ("t2_patch_tracao.inp", build_patch_traction()),
        ("t0_hydrostatic_block.inp", build_hydrostatic()),
    ):
        write_inp(os.path.join(dest, name), pack["text"])
        write_inp(os.path.join(dest, name.replace(".inp", ".ccx.inp")), ccx_compact_text(ccx_decorate(pack["text"])))
        print(name, pack["meta"])
