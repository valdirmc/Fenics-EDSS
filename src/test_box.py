#!/usr/bin/env python3
"""
Teste de validação: estática linear com malha create_box.
Confirma se o solver funciona independentemente do parser .inp.
"""

import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
import ufl


def main(debug=1):
    # Malha gerada pelo DOLFINx
    domain = mesh.create_box(
        MPI.COMM_WORLD,
        [np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.2, 0.2])],
        [4, 2, 2],
        cell_type=mesh.CellType.hexahedron,
    )
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

    # BC: fixa a face x=0
    def clamped(x):
        return np.isclose(x[0], 0.0)

    fdim = domain.topology.dim - 1
    facets = mesh.locate_entities_boundary(domain, fdim, clamped)
    u_D = np.array([0.0, 0.0, 0.0], dtype=default_scalar_type)
    bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, facets), V)

    # Material
    E = default_scalar_type(210000.0)
    nu = default_scalar_type(0.3)
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2.0 * mu * epsilon(u)

    # Carga de corpo (peso)
    f = fem.Constant(domain, default_scalar_type((0.0, 0.0, -1000.0)))

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    L = ufl.dot(f, v) * ufl.dx

    problem = LinearProblem(
        a, L, bcs=[bc],
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        petsc_options_prefix="box_test_",
    )
    uh = problem.solve()
    uh.name = "displacement"

    norma = np.linalg.norm(uh.x.array)
    if debug:
        print(f"Norma do deslocamento: {norma:.6e}")

    # Saída
    import os
    os.makedirs("/app/output", exist_ok=True)
    with io.XDMFFile(domain.comm, "/app/output/box_test.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)
    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, "/app/output/box_test.pvd", "w") as vtk:
            vtk.write_function(uh)
        if debug:
            print("PVD gerado: /app/output/box_test.pvd")
    except Exception as e:
        if debug:
            print(f"Aviso: {e}")

    return norma


if __name__ == "__main__":
    main()
