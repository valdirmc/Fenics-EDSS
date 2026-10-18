#!/usr/bin/env python3
"""
Exemplo mínimo: Equação de Poisson em FEniCSx
Baseado no dolfinx-tutorial (Jørgen Dokken)
"""

import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem.petsc import LinearProblem
import ufl

def main():
    # Malha
    domain = mesh.create_unit_square(MPI.COMM_WORLD, 16, 16, mesh.CellType.triangle)
    
    # Espaço de funções
    V = fem.functionspace(domain, ("Lagrange", 1))
    
    # Condições de contorno (Dirichlet homogênea)
    def boundary(x):
        return np.logical_or(
            np.isclose(x[0], 0.0),
            np.logical_or(
                np.isclose(x[0], 1.0),
                np.logical_or(np.isclose(x[1], 0.0), np.isclose(x[1], 1.0))
            )
        )
    
    boundary_facets = mesh.locate_entities_boundary(domain, domain.topology.dim - 1, boundary)
    bc = fem.dirichletbc(fem.Constant(domain, 0.0), 
                         fem.locate_dofs_topological(V, domain.topology.dim - 1, boundary_facets), V)
    
    # Formulação variacional
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    f = fem.Constant(domain, 1.0)
    
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = f * v * ufl.dx
    
    # Resolução
    problem = LinearProblem(
        a, L, bcs=[bc],
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        petsc_options_prefix="poisson_"
    )
    uh = problem.solve()
    uh.name = "solution"
    
    # Saída
    import os
    os.makedirs("/app/output", exist_ok=True)
    
    with io.XDMFFile(domain.comm, "/app/output/poisson.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)
    
    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, "/app/output/poisson.pvd", "w") as vtk:
            vtk.write_function(uh)
        print("Arquivo PVD gerado: /app/output/poisson.pvd")
    except Exception as e:
        print(f"Aviso VTU: {e}")
    
    print("Poisson resolvido com sucesso.")
    print("Saída em: /app/output/")

if __name__ == "__main__":
    main()
