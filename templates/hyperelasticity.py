#!/usr/bin/env python3
"""
Template: Hiperelasticidade parametrizada
"""

import json
import sys
import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io, default_scalar_type
from dolfinx.fem.petsc import NonlinearProblem
import ufl

def load_params(path=None):
    defaults = {
        "L": 1.0,
        "W": 0.1,
        "nx": 16,
        "ny": 4,
        "nz": 4,
        "E": 1.0e4,
        "nu": 0.3,
        "load_z": -0.5,
        "output": "/app/output/hyper_param"
    }
    if path and path.endswith(".json"):
        with open(path) as f:
            defaults.update(json.load(f))
    return defaults

def main():
    param_file = sys.argv[1] if len(sys.argv) > 1 else None
    p = load_params(param_file)

    domain = mesh.create_box(
        MPI.COMM_WORLD,
        [np.array([0.0, 0.0, 0.0]), np.array([p["L"], p["W"], p["W"]])],
        [p["nx"], p["ny"], p["nz"]],
        cell_type=mesh.CellType.hexahedron,
    )

    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

    def clamped(x):
        return np.isclose(x[0], 0.0)

    fdim = domain.topology.dim - 1
    facets = mesh.locate_entities_boundary(domain, fdim, clamped)
    u_D = np.array([0.0, 0.0, 0.0], dtype=default_scalar_type)
    bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, facets), V)

    E = default_scalar_type(p["E"])
    nu = default_scalar_type(p["nu"])
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    u = fem.Function(V, name="displacement")

    d = len(u)
    I = ufl.Identity(d)
    F = ufl.variable(I + ufl.grad(u))
    C = F.T * F
    Ic = ufl.tr(C)
    J = ufl.det(F)

    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2
    P = ufl.diff(psi, F)

    B = fem.Constant(domain, default_scalar_type((0.0, 0.0, p["load_z"])))
    v = ufl.TestFunction(V)
    residual = ufl.inner(ufl.grad(v), P) * ufl.dx - ufl.inner(v, B) * ufl.dx

    problem = NonlinearProblem(
        residual, u, bcs=[bc],
        petsc_options={
            "snes_type": "newtonls",
            "snes_linesearch_type": "none",
            "snes_atol": 1e-8,
            "snes_rtol": 1e-8,
            "ksp_type": "preonly",
            "pc_type": "lu",
        },
        petsc_options_prefix="hyper_param_",
    )
    problem.solve()

    import os
    os.makedirs("/app/output", exist_ok=True)
    out_base = p["output"]

    with io.XDMFFile(domain.comm, f"{out_base}.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(u)

    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, f"{out_base}.pvd", "w") as vtk:
            vtk.write_function(u)
        print(f"PVD gerado: {out_base}.pvd")
    except Exception as e:
        print(f"Aviso VTU: {e}")

    print(f"Hiperelasticidade parametrizada concluída.")
    print(f"Norma: {np.linalg.norm(u.x.array):.6e}")

if __name__ == "__main__":
    main()