#!/usr/bin/env python3
"""
Executa análise hiperelástica (Neo-Hooke) a partir de um arquivo .inp.
"""

import sys
import os
import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io, default_scalar_type, log
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from petsc4py import PETSc
import ufl
import basix.ufl

sys.path.insert(0, "/app/src")
from inp_parser import parse_inp, summary
from run_from_inp import build_mesh_from_inp


def run_hyperelastic(model, step, output_base: str, debug: int = 0):
    """Executa um step hiperelástico."""
    domain, node_map = build_mesh_from_inp(model, debug=debug)
    
    domain.topology.create_connectivity(0, domain.topology.dim)
    domain.topology.create_connectivity(domain.topology.dim, 0)
    
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

    # Parâmetros Neo-Hooke (do material Borracha ou padrão)
    mu = default_scalar_type(1.0)
    lmbda = default_scalar_type(0.1)
    for mat in model.materials.values():
        if "BORRACHA" in mat.name.upper() or mat.E == 0:
            # Usa valores do HYPERELASTIC se disponíveis
            mu = default_scalar_type(1.0)
            lmbda = default_scalar_type(0.1)
            break

    if debug:
        print(f"Neo-Hooke: mu={float(mu)}, lambda={float(lmbda)}")

    # BCs
    bcs = []
    coords = domain.geometry.x
    fixed_nodes = set(bc.node for bc in model.boundary if bc.node in model.nodes)
    
    for nid in fixed_nodes:
        n = model.nodes[nid]
        def marker(x, nx=n.x, ny=n.y, nz=n.z):
            return (np.isclose(x[0], nx, atol=1e-6) &
                    np.isclose(x[1], ny, atol=1e-6) &
                    np.isclose(x[2], nz, atol=1e-6))
        dofs = fem.locate_dofs_geometrical(V, marker)
        if len(dofs) > 0:
            u_bc = fem.Function(V)
            u_bc.x.array[:] = 0.0
            bcs.append(fem.dirichletbc(u_bc, dofs))

    # Formulação hiperelástica
    u = fem.Function(V)
    v = ufl.TestFunction(V)
    
    d = len(u)
    I = ufl.Identity(d)
    F = I + ufl.grad(u)
    C = F.T * F
    Ic = ufl.tr(C)
    J = ufl.det(F)
    
    # Energia Neo-Hooke
    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J))**2
    
    # Força de corpo (das cargas do step)
    f_body = fem.Constant(domain, default_scalar_type((0.0, 0.0, 0.0)))
    if step.loads:
        total = np.zeros(3)
        for cl in step.loads:
            total[cl.dof - 1] += cl.value
        f_body = fem.Constant(domain, default_scalar_type(tuple(total * 0.01)))
    
    # Residual
    F_form = ufl.derivative(psi * ufl.dx, u, v) - ufl.dot(f_body, v) * ufl.dx
    
    problem = NonlinearProblem(
        F_form, u, bcs=bcs,
        petsc_options={"snes_type": "newtonls", "snes_monitor": ""},
        petsc_options_prefix="hyper_"
    )
    
    problem.solve()
    u.x.scatter_forward()
    
    if debug:
        print("Newton concluído")

    # Saída
    os.makedirs(os.path.dirname(output_base) or ".", exist_ok=True)
    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, f"{output_base}.pvd", "w") as vtk:
            vtk.write_function(u)
        if debug:
            print(f"PVD: {output_base}.pvd")
    except Exception as e:
        if debug:
            print(f"Aviso: {e}")

    norma = np.linalg.norm(u.x.array)
    if debug:
        print(f"Norma: {norma:.6e}")
    return u
