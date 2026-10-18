#!/usr/bin/env python3
"""
Executa condução térmica a partir de um arquivo .inp.
"""

import sys
import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
import ufl
import basix.ufl

sys.path.insert(0, "/app/src")
from inp_parser import parse_inp, summary


def run_thermal(inp_path: str, output_base: str = "/app/output/thermal", debug: int = 0):
    model = parse_inp(inp_path, debug=debug)
    if debug:
        print(summary(model))
        print("---")

    if model.analysis_type != "heat":
        print("Aviso: análise não é térmica")

    # Monta malha (reusa a lógica do run_from_inp)
    from run_from_inp import build_mesh_from_inp
    domain, node_map = build_mesh_from_inp(model, debug=debug)
    V = fem.functionspace(domain, ("Lagrange", 1))

    # Material (pega o primeiro com k)
    mat = None
    for m in model.materials.values():
        if m.k > 0:
            mat = m
            break
    if mat is None:
        mat = next(iter(model.materials.values()))
    k = default_scalar_type(mat.k if mat.k > 0 else 1.0)
    if debug:
        print(f"Condutividade: {float(k)}")

    # BCs de temperatura
    bcs = []
    coords = domain.geometry.x
    for bc in model.boundary:
        if bc.dof != 11:  # DOF 11 = temperatura no CalculiX
            continue
        if bc.node not in model.nodes:
            continue
        n = model.nodes[bc.node]
        def marker(x, nx=n.x, ny=n.y, nz=n.z):
            return (np.isclose(x[0], nx, atol=1e-6) &
                    np.isclose(x[1], ny, atol=1e-6) &
                    np.isclose(x[2], nz, atol=1e-6))
        dofs = fem.locate_dofs_geometrical(V, marker)
        if len(dofs) > 0:
            bc_val = fem.Constant(domain, default_scalar_type(bc.value))
            bcs.append(fem.dirichletbc(bc_val, dofs, V))

    if debug:
        print(f"BCs de temperatura: {len(bcs)}")

    # Formulação
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = k * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = fem.Constant(domain, default_scalar_type(0.0)) * v * ufl.dx

    problem = LinearProblem(
        a, L, bcs=bcs,
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        petsc_options_prefix="thermal_",
    )
    uh = problem.solve()
    uh.name = "temperature"

    # Saída
    import os
    os.makedirs("/app/output", exist_ok=True)
    with io.XDMFFile(domain.comm, f"{output_base}.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)
    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, f"{output_base}.pvd", "w") as vtk:
            vtk.write_function(uh)
        if debug:
            print(f"PVD gerado: {output_base}.pvd")
    except Exception as e:
        if debug:
            print(f"Aviso: {e}")

    tmax = uh.x.array.max()
    tmin = uh.x.array.min()
    if debug:
        print(f"Temperatura: min={tmin:.2f}, max={tmax:.2f}")
    return uh


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python3 run_thermal.py arquivo.inp")
        sys.exit(1)
    run_thermal(sys.argv[1])
