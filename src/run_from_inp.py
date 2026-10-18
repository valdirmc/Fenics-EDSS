#!/usr/bin/env python3
"""
Executa estática linear a partir de um arquivo .inp.
Usa o parser mínimo e monta o problema no DOLFINx.
"""

import sys
import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io, default_scalar_type
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, set_bc
from petsc4py import PETSc
import ufl
import basix.ufl

sys.path.insert(0, "/app/src")
from inp_parser import parse_inp, summary


def build_mesh_from_inp(model, debug=0):
    """Constrói uma malha DOLFINx a partir do modelo INP."""
    if not model.elements:
        raise ValueError("Nenhum elemento encontrado no .inp")

    first_el = next(iter(model.elements.values()))
    etype = first_el.type.upper()

    node_ids = sorted(model.nodes.keys())
    node_map = {nid: i for i, nid in enumerate(node_ids)}

    coords = np.array([[model.nodes[nid].x, model.nodes[nid].y, model.nodes[nid].z]
                       for nid in node_ids], dtype=np.float64)

    cells = []
    for el in model.elements.values():
        if el.type.upper() != etype:
            raise ValueError(f"Elementos mistos não suportados ({el.type} vs {etype})")
        nodes = [node_map[n] for n in el.nodes]
        # Reordena hexaedro de CalculiX para DOLFINx/Basix
        # CalculiX: 0,1,2,3,4,5,6,7
        # DOLFINx:  0,1,3,2,4,5,7,6
        if etype == "C3D8" and len(nodes) == 8:
            nodes = [nodes[0], nodes[1], nodes[3], nodes[2],
                     nodes[4], nodes[5], nodes[7], nodes[6]]
        cells.append(nodes)
    cells = np.array(cells, dtype=np.int32)

    if etype == "C3D8":
        cell_name = "hexahedron"
    elif etype == "C3D4":
        cell_name = "tetrahedron"
    else:
        raise ValueError(f"Tipo de elemento não suportado: {etype}")

    if debug:
        print(f"Malha: {len(node_ids)} nós, {len(cells)} elementos ({etype})")

    coord_element = basix.ufl.element("Lagrange", cell_name, 1, shape=(coords.shape[1],))
    ufl_mesh = ufl.Mesh(coord_element)
    domain = mesh.create_mesh(MPI.COMM_WORLD, cells, ufl_mesh, coords)
    return domain, node_map


def run_linear_static(inp_path: str, output_base: str = "/app/output/from_inp", debug: int = 1):
    model = parse_inp(inp_path)
    # Material
    mat = next(iter(model.materials.values()))
    E = default_scalar_type(mat.E)
    nu = default_scalar_type(mat.nu)
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    if debug:
        print(summary(model))
        print("---")
        print(f"Material: E={mat.E}, nu={mat.nu}")
        print(f"mu={float(mu):.2f}, lambda={float(lmbda):.2f}")

    domain, node_map = build_mesh_from_inp(model, debug=debug)
    
    # Cria conectividades
    domain.topology.create_connectivity(0, domain.topology.dim)
    domain.topology.create_connectivity(domain.topology.dim, 0)
    domain.topology.create_connectivity(domain.topology.dim - 1, domain.topology.dim)
    
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

    # Material
    mat = next(iter(model.materials.values()))
    E = default_scalar_type(mat.E)
    nu = default_scalar_type(mat.nu)
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    # Condições de contorno por nó/set
    bcs = []
    coords = domain.geometry.x
    
    # Coleta nós únicos com BC
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
    
    if debug:
        print(f"BCs aplicadas: {len(bcs)} (nós fixos: {len(fixed_nodes)})")

    # Formulação
    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2.0 * mu * epsilon(u)

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    L_form = ufl.dot(fem.Constant(domain, default_scalar_type((0.0, 0.0, 0.0))), v) * ufl.dx

    a_form = fem.form(a)
    L_compiled = fem.form(L_form)
    A = assemble_matrix(a_form, bcs=bcs)
    A.assemble()
    b = assemble_vector(L_compiled)

    if debug:
        print(f"Vetor b (antes): norma = {b.norm():.6e}")

    # Aplica CLOAD
    if model.cloads:
        coords = domain.geometry.x
        for cl in model.cloads:
            if cl.node not in model.nodes:
                continue
            n = model.nodes[cl.node]
            # Encontra o índice do nó na geometria
            dist = np.linalg.norm(coords - np.array([n.x, n.y, n.z]), axis=1)
            node_idx = np.argmin(dist)
            if dist[node_idx] > 1e-6:
                continue
            # Para espaço vetorial Lagrange, os dofs são: node_idx * 3 + componente
            comp = cl.dof - 1
            dof_comp = node_idx * 3 + comp
            with b.localForm() as local:
                if dof_comp < local.size:
                    local[dof_comp] += cl.value
                    if debug:
                        print(f"  CLOAD nó {cl.node}, comp {comp}: dof {dof_comp}")
        if debug:
            print(f"Cargas aplicadas: {len(model.cloads)}")

    if debug:
        print(f"Vetor b (depois): norma = {b.norm():.6e}")

    apply_lifting(b, [a_form], [bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b, bcs)

    solver = PETSc.KSP().create(domain.comm)
    solver.setOperators(A)
    solver.setType("preonly")
    solver.getPC().setType("lu")
    uh = fem.Function(V)
    solver.solve(b, uh.x.petsc_vec)
    uh.x.scatter_forward()
    uh.name = "displacement"

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
            print(f"Aviso VTU: {e}")

    norma = np.linalg.norm(uh.x.array)
    if debug:
        print(f"Norma do deslocamento: {norma:.6e}")
    return uh


if __name__ == "__main__":
    debug = 1
    if len(sys.argv) > 2 and sys.argv[2] == "--quiet":
        debug = 0
    if len(sys.argv) < 2:
        print("Uso: python3 run_from_inp.py arquivo.inp [--quiet]")
        sys.exit(1)
    run_linear_static(sys.argv[1], debug=debug)
