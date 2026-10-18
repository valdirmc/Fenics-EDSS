#!/usr/bin/env python3
"""
Teste DOLFINx direto com os dados do simple_beam.inp.
Gerado a partir da saída do run_from_inp_Text.py.
"""

import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io, default_scalar_type
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, set_bc
from petsc4py import PETSc
import ufl
import basix.ufl


def main():
    # Coordenadas (do .inp)
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.2, 0.0],
        [0.0, 0.2, 0.0],
        [0.0, 0.0, 0.2],
        [1.0, 0.0, 0.2],
        [1.0, 0.2, 0.2],
        [0.0, 0.2, 0.2],
    ], dtype=np.float64)

    # Conectividade (1 elemento C3D8)
    cells = np.array([[0, 1, 2, 3, 4, 5, 6, 7]], dtype=np.int32)

    # Cria malha
    coord_element = basix.ufl.element("Lagrange", "hexahedron", 1, shape=(3,))
    ufl_mesh = ufl.Mesh(coord_element)
    domain = mesh.create_mesh(MPI.COMM_WORLD, cells, ufl_mesh, coords)

    # Cria conectividades
    domain.topology.create_connectivity(0, 3)
    domain.topology.create_connectivity(3, 0)
    domain.topology.create_connectivity(2, 3)

    V = fem.functionspace(domain, ("Lagrange", 1, (3,)))

    # Material
    E = default_scalar_type(210000.0)
    nu = default_scalar_type(0.3)
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    # BCs: fixa a face x=0 usando locate_dofs_topological
    def clamped(x):
        return np.isclose(x[0], 0.0)

    fdim = domain.topology.dim - 1
    facets = mesh.locate_entities_boundary(domain, fdim, clamped)
    print(f"Facets fixas encontradas: {len(facets)}")

    if len(facets) > 0:
        u_D = np.array([0.0, 0.0, 0.0], dtype=default_scalar_type)
        dofs = fem.locate_dofs_topological(V, fdim, facets)
        bc = fem.dirichletbc(u_D, dofs, V)
        bcs = [bc]
    else:
        bcs = []
        print("AVISO: nenhuma face fixa encontrada")

    print(f"BCs aplicadas: {len(bcs)}")

    # Formulação
    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(3) + 2.0 * mu * epsilon(u)

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    L_form = ufl.dot(fem.Constant(domain, default_scalar_type((0.0, 0.0, 0.0))), v) * ufl.dx

    a_form = fem.form(a)
    L_compiled = fem.form(L_form)
    A = assemble_matrix(a_form, bcs=bcs)
    A.assemble()
    b = assemble_vector(L_compiled)

    # Cargas: nós 7 e 8 do .inp → índices 6 e 7, DOF 3 (z)
    load_nodes = [6, 7]
    load_value = -100.0
    for idx in load_nodes:
        def marker(x, i=idx):
            return (np.isclose(x[0], coords[i, 0], atol=1e-6) &
                    np.isclose(x[1], coords[i, 1], atol=1e-6) &
                    np.isclose(x[2], coords[i, 2], atol=1e-6))
        dofs = fem.locate_dofs_geometrical(V, marker)
        if len(dofs) > 0:
            # Componente z é o terceiro dof do nó
            dof_z = dofs[2] if len(dofs) > 2 else dofs[-1]
            with b.localForm() as local:
                local[dof_z] += load_value
            print(f"Carga {load_value} aplicada no dof {dof_z} (nó índice {idx})")

    apply_lifting(b, [a_form], [bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b, bcs)

    # Resolve
    solver = PETSc.KSP().create(domain.comm)
    solver.setOperators(A)
    solver.setType("preonly")
    solver.getPC().setType("lu")
    uh = fem.Function(V)
    solver.solve(b, uh.x.petsc_vec)
    uh.x.scatter_forward()
    uh.name = "displacement"

    norma = np.linalg.norm(uh.x.array)
    print(f"Norma do deslocamento: {norma:.6e}")
    print(f"Min: {uh.x.array.min():.6e}, Max: {uh.x.array.max():.6e}")

    # Saída
    import os
    os.makedirs("/app/output", exist_ok=True)
    with io.XDMFFile(domain.comm, "/app/output/test_dolfin.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)
    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, "/app/output/test_dolfin.pvd", "w") as vtk:
            vtk.write_function(uh)
        print("PVD gerado: /app/output/test_dolfin.pvd")
    except Exception as e:
        print(f"Aviso: {e}")


if __name__ == "__main__":
    main()
