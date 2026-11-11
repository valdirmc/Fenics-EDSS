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
