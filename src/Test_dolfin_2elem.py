#!/usr/bin/env python3
"""
Teste DOLFINx com malha de 2 elementos.
"""

import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
import ufl
import basix.ufl


def main():
    # Coordenadas (2 elementos ao longo de x)
    coords = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.2, 0.0],
        [0.5, 0.2, 0.0],
        [1.0, 0.2, 0.0],
        [0.0, 0.0, 0.2],
        [0.5, 0.0, 0.2],
        [1.0, 0.0, 0.2],
        [0.0, 0.2, 0.2],
        [0.5, 0.2, 0.2],
        [1.0, 0.2, 0.2],
    ], dtype=np.float64)

    # 2 hexaedros — ordenação DOLFINx/Basix
    # CalculiX: 0,1,2,3,4,5,6,7
    # DOLFINx:  0,1,3,2,4,5,7,6  (reordena faces)
    cells = np.array([
        [0, 1, 3, 4, 6, 7, 9, 10],   # reordenado
        [1, 2, 4, 5, 7, 8, 10, 11],  # reordenado
    ], dtype=np.int32)

    coord_element = basix.ufl.element("Lagrange", "hexahedron", 1, shape=(3,))
    ufl_mesh = ufl.Mesh(coord_element)
    domain = mesh.create_mesh(MPI.COMM_WORLD, cells, ufl_mesh, coords)

    domain.topology.create_connectivity(0, 3)
    domain.topology.create_connectivity(3, 0)
    domain.topology.create_connectivity(2, 3)

    V = fem.functionspace(domain, ("Lagrange", 1, (3,)))

    # BC: face x=0
    def clamped(x):
        return np.isclose(x[0], 0.0)

    fdim = domain.topology.dim - 1
    facets = mesh.locate_entities_boundary(domain, fdim, clamped)
    print(f"Facets fixas: {len(facets)}")

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
        return lmbda * ufl.nabla_div(u) * ufl.Identity(3) + 2.0 * mu * epsilon(u)

    # Carga de corpo
    f = fem.Constant(domain, default_scalar_type((0.0, 0.0, -1000.0)))

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    L = ufl.dot(f, v) * ufl.dx

    problem = LinearProblem(
        a, L, bcs=[bc],
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        petsc_options_prefix="beam2_",
    )
    uh = problem.solve()
    uh.name = "displacement"

    norma = np.linalg.norm(uh.x.array)
    print(f"Norma: {norma:.6e}")
    print(f"Min: {uh.x.array.min():.6e}, Max: {uh.x.array.max():.6e}")

    import os
    os.makedirs("/app/output", exist_ok=True)
    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, "/app/output/beam_2elem.pvd", "w") as vtk:
            vtk.write_function(uh)
        print("PVD: /app/output/beam_2elem.pvd")
    except Exception as e:
        print(f"Aviso: {e}")


if __name__ == "__main__":
    main()
