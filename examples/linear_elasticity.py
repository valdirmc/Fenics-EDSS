#!/usr/bin/env python3
"""
Exemplo: Elasticidade Linear 3D (viga engastada sob peso próprio)
Baseado no dolfinx-tutorial (Jørgen S. Dokken)
"""

import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
import ufl

def main():
    # Parâmetros da viga
    L, W = 1.0, 0.2
    
    # Malha 3D (caixa)
    domain = mesh.create_box(
        MPI.COMM_WORLD,
        [np.array([0.0, 0.0, 0.0]), np.array([L, W, W])],
        [20, 6, 6],
        cell_type=mesh.CellType.hexahedron,
    )
    
    # Espaço de funções vetorial (deslocamento)
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))
    
    # Condição de contorno: engastada em x=0
    def clamped_boundary(x):
        return np.isclose(x[0], 0.0)
    
    fdim = domain.topology.dim - 1
    boundary_facets = mesh.locate_entities_boundary(domain, fdim, clamped_boundary)
    u_D = np.array([0.0, 0.0, 0.0], dtype=default_scalar_type)
    bc = fem.dirichletbc(
        u_D,
        fem.locate_dofs_topological(V, fdim, boundary_facets),
        V,
    )
    
    # Parâmetros do material (Lamé)
    E = default_scalar_type(1.0e5)
    nu = default_scalar_type(0.3)
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    
    # Força de corpo (peso próprio)
    rho = default_scalar_type(1.0)
    g = default_scalar_type(0.4 * (W / L) ** 2)
    f = fem.Constant(domain, default_scalar_type((0.0, 0.0, -rho * g)))
    
    # Funções de deformação e tensão
    def epsilon(u):
        return ufl.sym(ufl.grad(u))
    
    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2.0 * mu * epsilon(u)
    
    # Formulação variacional
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    L_form = ufl.dot(f, v) * ufl.dx
    
    # Resolução
    problem = LinearProblem(
        a,
        L_form,
        bcs=[bc],
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        petsc_options_prefix="elasticity_",
    )
    uh = problem.solve()
    uh.name = "displacement"
    
    # Saída
    import os
    os.makedirs("/app/output", exist_ok=True)
    
    # XDMF (formato nativo do DOLFINx)
    with io.XDMFFile(domain.comm, "/app/output/elasticity.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)
    
    # VTU (mais estável no ParaView)
    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, "/app/output/elasticity.pvd", "w") as vtk:
            vtk.write_function(uh)
        print("Arquivo VTU/PVD gerado: /app/output/elasticity.pvd")
    except Exception as e:
        print(f"Aviso: não foi possível gerar VTU ({e})")
    
    print("Elasticidade linear resolvida com sucesso.")
    print("Saída principal: /app/output/elasticity.xdmf")
    print(f"Norma do deslocamento: {np.linalg.norm(uh.x.array):.6e}")

if __name__ == "__main__":
    main()
