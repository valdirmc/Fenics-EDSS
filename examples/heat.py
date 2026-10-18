#!/usr/bin/env python3
"""
Exemplo: Equação do Calor (estacionária simplificada)
Domínio unitário com fonte constante e Dirichlet homogêneo.
"""

import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io
from dolfinx.fem.petsc import LinearProblem
import ufl

def main():
    # Malha
    domain = mesh.create_unit_square(MPI.COMM_WORLD, 32, 32, mesh.CellType.triangle)
    
    # Espaço de funções
    V = fem.functionspace(domain, ("Lagrange", 1))
    
    # Condição de contorno (temperatura = 0 na borda)
    def boundary(x):
        return np.logical_or(
            np.isclose(x[0], 0.0),
            np.logical_or(
                np.isclose(x[0], 1.0),
                np.logical_or(np.isclose(x[1], 0.0), np.isclose(x[1], 1.0))
            )
        )
    
    fdim = domain.topology.dim - 1
    boundary_facets = mesh.locate_entities_boundary(domain, fdim, boundary)
    bc = fem.dirichletbc(
        fem.Constant(domain, 0.0),
        fem.locate_dofs_topological(V, fdim, boundary_facets),
        V,
    )
    
    # Fonte de calor
    f = fem.Constant(domain, 10.0)
    
    # Formulação
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = f * v * ufl.dx
    
    # Resolução
    problem = LinearProblem(
        a, L, bcs=[bc],
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        petsc_options_prefix="heat_",
    )
    uh = problem.solve()
    uh.name = "temperature"
    
    # Saída
    import os
    os.makedirs("/app/output", exist_ok=True)
    
    with io.XDMFFile(domain.comm, "/app/output/heat.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)
    
    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, "/app/output/heat.pvd", "w") as vtk:
            vtk.write_function(uh)
        print("Arquivo PVD gerado: /app/output/heat.pvd")
    except Exception as e:
        print(f"Aviso VTU: {e}")
    
    print("Equação do calor resolvida com sucesso.")
    print(f"Temperatura máxima: {uh.x.array.max():.4f}")

if __name__ == "__main__":
    main()
