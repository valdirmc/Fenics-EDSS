"""Campos obrigatórios do perfil FC para C3D4 P1.

Localizacao:
- U (displacement): nodal, espaco Lagrange P1 / CG1.
- S (Cauchy), E (deformacao infinitesimal), principais e von Mises:
  por celula, espaco DG0, constantes no tetra linear C3D4 (P1).
  Media/suavizacao nodal NAO substitui o campo por celula.
- RF: nos presos e totais (arquivo JSON +, se possivel, campo nodal).

Ordenacao:
- Tensor 3x3 em row-major: XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ.
- Voigt auxiliar (documentado): XX, YY, ZZ, XY, YZ, XZ.
- Principais algebricos com sinal: sigma1 >= sigma2 >= sigma3
  (autovalores de Cauchy, ordem decrescente, nao por modulo).
- von Mises: sqrt(0.5*((s1-s2)^2+(s2-s3)^2+(s3-s1)^2)).

Elastica linear de pequenas deformacoes: Cauchy e identificada com
sigma = lambda*tr(eps)*I + 2*mu*eps. Nao e identidade com CalculiX
alegada; testes analiticos cobrem este recorte.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


def principals_algebraic(sigma: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """sigma: (n,3,3) simetrico. Retorna (n,3) com s1>=s2>=s3 e von Mises (n,)."""
    if sigma.ndim != 3 or sigma.shape[1:] != (3, 3):
        raise ValueError("sigma deve ter forma (n,3,3)")
    # C3D4 P1: tensao constante por celula. Simetriza residual numerico
    # (sigma = (S+S.T)/2) antes dos autovalores; eigvalsh exige simetria.
    sigma = 0.5 * (np.asarray(sigma, dtype=float) + np.swapaxes(sigma, -1, -2))
    eig = np.linalg.eigvalsh(sigma)
    principals = eig[:, ::-1].copy()
    s1, s2, s3 = principals[:, 0], principals[:, 1], principals[:, 2]
    vm = np.sqrt(0.5 * ((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s3 - s1) ** 2))
    return principals, vm


def voigt_from_tensor(sigma: np.ndarray) -> np.ndarray:
    """(n,3,3) -> (n,6) XX,YY,ZZ,XY,YZ,XZ."""
    return np.column_stack((
        sigma[:, 0, 0],
        sigma[:, 1, 1],
        sigma[:, 2, 2],
        sigma[:, 0, 1],
        sigma[:, 1, 2],
        sigma[:, 0, 2],
    ))


def interpolate_cell_tensor(domain, expr, name: str):
    """Projeta expressao tensorial 3x3 no DG0 (constante por celula)."""
    from dolfinx import fem

    W = fem.functionspace(domain, ("DG", 0, (3, 3)))
    func = fem.Function(W)
    func.name = name
    points = W.element.interpolation_points
    if callable(points):
        points = points()
    func.interpolate(fem.Expression(expr, points))
    func.x.scatter_forward()
    return func


def interpolate_cell_vector(domain, values: np.ndarray, name: str):
    from dolfinx import fem

    W = fem.functionspace(domain, ("DG", 0, (3,)))
    func = fem.Function(W)
    func.name = name
    arr = values.reshape(-1)
    if arr.size != func.x.array.size:
        raise RuntimeError(
            "Tamanho do campo {0} ({1}) != dofs DG0 ({2})".format(
                name, arr.size, func.x.array.size
            )
        )
    func.x.array[:] = arr
    func.x.scatter_forward()
    return func


def interpolate_cell_scalar(domain, values: np.ndarray, name: str):
    from dolfinx import fem

    W = fem.functionspace(domain, ("DG", 0))
    func = fem.Function(W)
    func.name = name
    if values.size != func.x.array.size:
        raise RuntimeError(
            "Tamanho do campo {0} ({1}) != dofs DG0 ({2})".format(
                name, values.size, func.x.array.size
            )
        )
    func.x.array[:] = values.reshape(-1)
    func.x.scatter_forward()
    return func


def tensor_array_from_function(func) -> np.ndarray:
    arr = np.asarray(func.x.array)
    n = arr.size // 9
    return arr.reshape((n, 3, 3))


def write_functions_vtk(comm, functions: Sequence, output_base: str) -> str:
    """Grava varias funcoes no mesmo .pvd; falha e RuntimeError (nunca sucesso)."""
    import os
    from dolfinx.io import VTKFile
    from job_io import ensure_relative_vtk_links

    out_dir = os.path.dirname(os.path.abspath(output_base))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    pvd_path = "{0}.pvd".format(output_base)
    try:
        with VTKFile(comm, pvd_path, "w") as vtk:
            for fn in functions:
                vtk.write_function(fn)
        ensure_relative_vtk_links(pvd_path)
    except Exception as exc:
        raise RuntimeError("Falha ao gravar VTK {0}: {1}".format(pvd_path, exc)) from exc
    if not os.path.isfile(pvd_path):
        raise RuntimeError("Falha ao gravar VTK: arquivo ausente {0}".format(pvd_path))
    return pvd_path
