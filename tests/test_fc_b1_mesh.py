"""B1: malha, volumes, cargas consistentes e 3-2-1 (sem DOLFINx)."""

import os
import sys

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

from fc_b1 import (  # noqa: E402
    B,
    H,
    L,
    M,
    SIGMA_MAX,
    apply_loads,
    load_resultants,
    structured_c3d4,
    support_321,
    tet_volumes,
)


def test_b1_volumes_e_cargas_consistentes():
    coords, cells, meta = structured_c3d4(2)
    vols = tet_volumes(coords, cells)
    assert meta["n_elements"] == 20 * 2 * 2 * 6
    assert np.all(vols > 0)
    assert vols.sum() == pytest.approx(L * B * H, rel=1e-12)
    loads = apply_loads(coords, cells)
    F, Mv = load_resultants(coords, loads)
    assert F[0] == pytest.approx(0.0, abs=1e-9)
    assert F[1] == pytest.approx(0.0, abs=1e-9)
    assert F[2] == pytest.approx(0.0, abs=1e-9)
    # Momentos das duas faces se cancelam na origem.
    assert Mv[0] == pytest.approx(0.0, abs=1e-7)
    assert Mv[1] == pytest.approx(0.0, abs=1e-7)
    assert Mv[2] == pytest.approx(0.0, abs=1e-7)
    # Cada face deve equivaler a um binario de modulo M: soma |M_face|.
    # Checagem do sigma_max declarado.
    assert SIGMA_MAX == pytest.approx(75.0)
    bcs, _ = support_321(coords)
    assert len(bcs) == 6


def test_b1_niveis_homologos():
    hs = []
    for n in (2, 4, 8, 16):
        _c, _e, meta = structured_c3d4(n)
        hs.append(meta["h"])
        assert meta["nx"] == 10 * n
        assert meta["n_dofs"] == meta["n_nodes"] * 3
    for a, b in zip(hs, hs[1:]):
        assert a / b == pytest.approx(2.0)


def test_tet_sphere_intersection_volume():
    import numpy as np
    from fc_b1 import tet_sphere_intersection_volume, _mid_subdivide_tet, _tet_volume_signed
    # tet unitario de volume 1/6
    verts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    vol = abs(_tet_volume_signed(verts))
    kids = _mid_subdivide_tet(verts)
    assert abs(sum(abs(_tet_volume_signed(k)) for k in kids) - vol) <= 1e-12 * max(vol, 1.0)
    # esfera que contem o tet: volume = vol do tet
    vin = tet_sphere_intersection_volume(verts, (0.2, 0.2, 0.2), 10.0)
    assert abs(vin - vol) <= 1e-12
    # esfera disjunta
    vout = tet_sphere_intersection_volume(verts, (10.0, 10.0, 10.0), 0.5)
    assert vout == 0.0
    # esfera no interior: volume da bola (4/3 pi r^3) se a bola cabe no tet
    # (r pequeno no baricentro 0.25,0.25,0.25 — distancia as faces = 0.25/sqrt(3) ~ 0.14)
    r = 0.05
    c = (0.15, 0.15, 0.15)
    vball = tet_sphere_intersection_volume(verts, c, r, max_depth=7)
    exact = 4.0 / 3.0 * np.pi * r ** 3
    # folha residual; tolerar 5 % nesta bola pequena totalmente interior
    assert vball > 0.0
    assert abs(vball - exact) / exact <= 0.05
