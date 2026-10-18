#!/usr/bin/env python3
"""
Template: Elasticidade Linear parametrizada
Com validação de parâmetros e logging.
"""

import json
import sys
import logging
import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
import ufl

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("elasticity")

def load_params(path=None):
    defaults = {
        "L": 1.0, "W": 0.2,
        "nx": 20, "ny": 6, "nz": 6,
        "E": 1.0e5, "nu": 0.3,
        "rho": 1.0, "g_scale": 0.4,
        "output": "/app/output/elasticity_param"
    }
    if path and path.endswith(".json"):
        log.info(f"Carregando parâmetros de {path}")
        with open(path) as f:
            defaults.update(json.load(f))
    return defaults

def validate(p):
    if p["L"] <= 0 or p["W"] <= 0:
        raise ValueError("L e W devem ser positivos")
    if p["nx"] < 1 or p["ny"] < 1 or p["nz"] < 1:
        raise ValueError("Divisões da malha devem ser >= 1")
    if p["E"] <= 0:
        raise ValueError("Módulo de Young (E) deve ser positivo")
    if not (0 < p["nu"] < 0.5):
        raise ValueError("Poisson (nu) deve estar entre 0 e 0.5")
    if p["rho"] <= 0:
        raise ValueError("Densidade (rho) deve ser positiva")
    log.info("Parâmetros validados")

def main():
    param_file = sys.argv[1] if len(sys.argv) > 1 else None
    p = load_params(param_file)
    validate(p)
    log.info(f"Parâmetros: L={p['L']}, W={p['W']}, E={p['E']}, nu={p['nu']}")

    domain = mesh.create_box(
        MPI.COMM_WORLD,
        [np.array([0.0, 0.0, 0.0]), np.array([p["L"], p["W"], p["W"]])],
        [p["nx"], p["ny"], p["nz"]],
        cell_type=mesh.CellType.hexahedron,
    )
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

    def clamped_boundary(x):
        return np.isclose(x[0], 0.0)

    fdim = domain.topology.dim - 1
    facets = mesh.locate_entities_boundary(domain, fdim, clamped_boundary)
    u_D = np.array([0.0, 0.0, 0.0], dtype=default_scalar_type)
    bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, facets), V)

    E = default_scalar_type(p["E"])
    nu = default_scalar_type(p["nu"])
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    rho = default_scalar_type(p["rho"])
    g = default_scalar_type(p["g_scale"] * (p["W"] / p["L"]) ** 2)
    f = fem.Constant(domain, default_scalar_type((0.0, 0.0, -rho * g)))

    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2.0 * mu * epsilon(u)

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    L_form = ufl.dot(f, v) * ufl.dx

    log.info("Resolvendo sistema linear...")
    problem = LinearProblem(
        a, L_form, bcs=[bc],
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        petsc_options_prefix="elast_param_",
    )
    uh = problem.solve()
    uh.name = "displacement"

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
        log.info(f"PVD gerado: {out_base}.pvd")
    except Exception as e:
        log.warning(f"Não foi possível gerar VTU: {e}")

    norma = np.linalg.norm(uh.x.array)
    log.info(f"Concluído. Norma do deslocamento: {norma:.6e}")

if __name__ == "__main__":
    main()
