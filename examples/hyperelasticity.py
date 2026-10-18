#!/usr/bin/env python3
"""
Exemplo: Hiperelasticidade (Neo-Hookean) - viga sob carga
Baseado no dolfinx-tutorial.
"""

import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io, default_scalar_type
from dolfinx.fem.petsc import NonlinearProblem
import ufl

def main():
    # Malha (viga 3D)
    L, W = 1.0, 0.1
    domain = mesh.create_box(
        MPI.COMM_WORLD,
        [np.array([0.0, 0.0, 0.0]), np.array([L, W, W])],
        [16, 4, 4],
        cell_type=mesh.CellType.hexahedron,
    )
    
    # Espaço de funções
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))
    
    # Condição de contorno: engastada em x=0
    def clamped(x):
        return np.isclose(x[0], 0.0)
    
    fdim = domain.topology.dim - 1
    facets = mesh.locate_entities_boundary(domain, fdim, clamped)
    u_D = np.array([0.0, 0.0, 0.0], dtype=default_scalar_type)
    bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, facets), V)
    
    # Parâmetros
    E = default_scalar_type(1.0e4)
    nu = default_scalar_type(0.3)
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    
    # Solução
    u = fem.Function(V, name="displacement")
    
    # Cinemática (F precisa ser variable para o diff)
    d = len(u)
    I = ufl.Identity(d)
    F = ufl.variable(I + ufl.grad(u))
    C = F.T * F
    Ic = ufl.tr(C)
    J = ufl.det(F)
    
    # Energia (Neo-Hookean)
    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2
    
    # Tensão de Piola-Kirchhoff
    P = ufl.diff(psi, F)
    
    # Força de corpo
    B = fem.Constant(domain, default_scalar_type((0.0, 0.0, -0.5)))
    
    # Residual
    v = ufl.TestFunction(V)
    residual = ufl.inner(ufl.grad(v), P) * ufl.dx - ufl.inner(v, B) * ufl.dx
    
    # Problema não-linear
    problem = NonlinearProblem(
        residual,
        u,
        bcs=[bc],
        petsc_options={
            "snes_type": "newtonls",
            "snes_linesearch_type": "none",
            "snes_atol": 1e-8,
            "snes_rtol": 1e-8,
            "ksp_type": "preonly",
            "pc_type": "lu",
        },
        petsc_options_prefix="hyper_",
    )
    
    problem.solve()
    
    # Saída
    import os
    os.makedirs("/app/output", exist_ok=True)
    
    with io.XDMFFile(domain.comm, "/app/output/hyperelasticity.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(u)
    
    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, "/app/output/hyperelasticity.pvd", "w") as vtk:
            vtk.write_function(u)
        print("Arquivo PVD gerado: /app/output/hyperelasticity.pvd")
    except Exception as e:
        print(f"Aviso VTU: {e}")
    
    print("Hiperelasticidade resolvida com sucesso.")
    print(f"Norma do deslocamento: {np.linalg.norm(u.x.array):.6e}")

if __name__ == "__main__":
    main()
