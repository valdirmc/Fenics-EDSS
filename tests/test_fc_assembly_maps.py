"""Mapas INP/DOF vetorizados e agrupamento de BC no perfil FC."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from fc_assembly import CoordLookup  # noqa: E402
from inp_parser import parse_inp, get_material_for_element  # noqa: E402


def test_coord_lookup_unique_and_loose_fallback():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    look = CoordLookup(coords, atol=1e-9, loose=1e-6)
    idx, status, dist = look.find((1.0, 0.0, 0.0))
    assert status == "ok" and idx == 1 and dist == pytest.approx(0.0)
    # dentro da tolerancia folgada, fora da estrita
    idx, status, dist = look.find((0.0, 0.0, 1.0 + 5e-7))
    assert status == "ok" and idx == 3
    idx, status, _ = look.find((8.0, 8.0, 8.0))
    assert status == "missing" and idx is None


def test_coord_lookup_rejects_duplicates():
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=float)
    look = CoordLookup(coords)
    idx, status, _ = look.find((0.0, 0.0, 0.0))
    assert status == "ambiguous" and idx is None


def test_get_material_for_element_is_not_quadratic(tmp_path):
    """List membership por elemento era O(N^2); o cache deve ficar linear."""
    n = 8000
    lines = ["*NODE"]
    # dois nos e um tet repetido so para o parser de material (elset grande)
    lines += [
        "1, 0, 0, 0",
        "2, 1, 0, 0",
        "3, 0, 1, 0",
        "4, 0, 0, 1",
        "*ELEMENT, TYPE=C3D4, ELSET=Eall",
    ]
    # um elemento real + elset com n ids (incluindo 1)
    lines.append("1, 1, 2, 3, 4")
    lines += ["*ELSET, ELSET=Ebig"]
    # chunk ids
    ids = [1] + list(range(2, n + 1))
    for i in range(0, len(ids), 20):
        lines.append(", ".join(str(x) for x in ids[i:i + 20]))
    lines += [
        "*MATERIAL, NAME=MAT",
        "*ELASTIC",
        "70000, 0.33",
        "*SOLID SECTION, ELSET=Eall, MATERIAL=MAT",
        "*BOUNDARY",
        "1, 1, 3, 0",
        "*STEP",
        "*STATIC",
        "*END STEP",
    ]
    inp = tmp_path / "mat.inp"
    inp.write_text("\n".join(lines) + "\n")
    model = parse_inp(str(inp), debug=0, profile="fc")
    import time
    t0 = time.perf_counter()
    for eid in range(1, n + 1):
        get_material_for_element(model, 1)
        get_material_for_element(model, eid)
    elapsed = time.perf_counter() - t0
    assert get_material_for_element(model, 1).E == pytest.approx(70000)
    # O(N^2) em 8000 levaria segundos; o cache deve ficar bem abaixo de 1 s
    assert elapsed < 1.0, "get_material_for_element ainda parece quadratico: {0:.3f}s".format(elapsed)


def test_locator_matches_geometrical_and_groups_bcs(tmp_path):
    pytest.importorskip("dolfinx")
    from dolfinx import fem
    from fc_assembly import collect_constrained_dofs, component_dof_locator
    from run_from_inp import build_mesh_from_inp

    text = """*NODE
1, 0, 0, 0
2, 1, 0, 0
3, 0, 1, 0
4, 0, 0, 1
5, 1, 1, 0
6, 1, 0, 1
*ELEMENT, TYPE=C3D4, ELSET=Eall
1, 1, 2, 3, 4
2, 2, 5, 3, 4
3, 2, 6, 5, 4
*MATERIAL, NAME=MAT
*ELASTIC
70000, 0.33
*SOLID SECTION, ELSET=Eall, MATERIAL=MAT
*BOUNDARY
1, 1, 3, 0
2, 3, 3, 0
2, 3, 3, 0.25
*STEP
*STATIC
*CLOAD
5, 3, -10
*END STEP
"""
    inp = tmp_path / "tiny.inp"
    inp.write_text(text)
    model = parse_inp(str(inp), debug=0, profile="fc")
    domain, _ = build_mesh_from_inp(model, debug=0)
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))
    locator = component_dof_locator(V)
    # compara com o criterio geometrico original (um no)
    from fc_assembly import locate_component_parent_dofs as loc
    # freeze old path? we replaced it. Reconstruct geometrical once:
    def geo(node, component):
        def marker(x, nx=node.x, ny=node.y, nz=node.z):
            return (
                np.isclose(x[0], nx, atol=1e-9)
                & np.isclose(x[1], ny, atol=1e-9)
                & np.isclose(x[2], nz, atol=1e-9)
            )
        Vcol, _map = V.sub(component).collapse()
        dofs = fem.locate_dofs_geometrical((V.sub(component), Vcol), marker)
        parent = np.asarray(dofs[0] if isinstance(dofs, (tuple, list)) and len(dofs) == 2 else dofs, dtype=np.int32)
        return parent

    for nid, node in model.nodes.items():
        for c in range(3):
            parent, pair, _v = locator.locate(node, c)
            g = geo(node, c)
            assert parent.size == 1 and g.size == 1
            assert int(parent[0]) == int(g[0])

    bcs, rows = collect_constrained_dofs(V, model, True)
    # 4 DOFs prescritos unicos: no 1 comps 0,1,2 e no 2 comp 2 (ultimo valor 0.25)
    assert len(rows) == 4
    # um bc por (componente, valor): (0,0), (1,0), (2,0), (2,0.25) => 4 objetos, nao 5 linhas
    assert len(bcs) == 4
    # valor que permanece no no 2 comp 2 e o ultimo
    assert any(r[0] == 2 and r[1] == 2 for r in rows)
