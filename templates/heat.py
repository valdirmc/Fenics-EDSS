#!/usr/bin/env python3
"""
Template: Equação do Calor parametrizada
"""

import json
import sys
import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem.petsc import LinearProblem
import ufl

def load_params(path=None):
    defaults = {
        "nx": 32,
        "ny": 32,
        "source": 10.0,
        "output": "/app/output/heat_param"
    }
    if path and path.endswith(".json"):
        with open(path) as f:
            defaults.update(json.load(f))
    return defaults

def main():
    param_file = sys.argv[1] if len(sys.argv) > 1 else None
    p = load_params(param_file)

    domain = mesh.create_unit_square(MPI.COMM_WORLD, p["nx"], p["ny"], mesh.CellType.triangle)
    V = fem.functionspace(domain, ("Lagrange", 1))

    def boundary(x):
        return np.logical_or(
            np.isclose(x[0], 0.0),
            np.logical_or(
                np.isclose(x[0], 1.0),
                np.logical_or(np.isclose(x[1], 0.0), np.isclose(x[1], 1.0))
            )
        )

    fdim = domain.topology.dim - 1
    facets = mesh.locate_entities_boundary(domain, fdim, boundary)
    bc = fem.dirichletbc(
        fem.Constant(domain, 0.0),
        fem.locate_dofs_topological(V, fdim, facets),
        V,
    )

    f = fem.Constant(domain, p["source"])
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = f * v * ufl.dx

    problem = LinearProblem(
        a, L, bcs=[bc],
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        petsc_options_prefix="heat_param_",
    )
    uh = problem.solve()
    uh.name = "temperature"

    import os
    os.makedirs("/app/output", exist_ok=True)
    out_base = p["output"]

    with io.XDMFFile(domain.comm, f"{out_base}.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)

    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, f"{out_base}.pvd", "w") as vtk:
            vtk.write_function(uh)
        print(f"PVD gerado: {out_base}.pvd")
    except Exception as e:
        print(f"Aviso VTU: {e}")

    print(f"Calor parametrizado concluído. T_max = {uh.x.array.max():.4f}")

if __name__ == "__main__":
    main()