#!/usr/bin/env python3
"""
Executa análise a partir de um arquivo .inp.
Detecta o tipo (estática ou térmica) e chama a função apropriada.
"""

import sys
import os
import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, io, default_scalar_type
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, set_bc
from petsc4py import PETSc
import ufl
import basix.ufl

sys.path.insert(0, "/app/src")
from inp_parser import parse_inp, summary, Step


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


def get_material_for_element(model, elem_id):
    """Retorna o material associado a um elemento via elset/solid section."""
    for elset_name, elem_ids in model.elsets.items():
        if elem_id in elem_ids:
            mat_name = model.solid_sections.get(elset_name.upper())
            if mat_name and mat_name.upper() in [k.upper() for k in model.materials.keys()]:
                for name, mat in model.materials.items():
                    if name.upper() == mat_name.upper():
                        return mat
    # Fallback: primeiro material
def run_linear_static_step(model, step, output_base: str, debug: int = 0):
    """Executa um step de estática linear com múltiplos materiais."""
    domain, node_map = build_mesh_from_inp(model, debug=debug)
    
    domain.topology.create_connectivity(0, domain.topology.dim)
    domain.topology.create_connectivity(domain.topology.dim, 0)
    
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

    # Materiais por célula (DG-0)
    Q = fem.functionspace(domain, ("DG", 0))
    E_func = fem.Function(Q)
    nu_func = fem.Function(Q)
    
    # Mapeia elementos para células
    elem_to_cell = {el.id: i for i, el in enumerate(model.elements.values())}
    
    for elset_name, elem_ids in model.elsets.items():
        mat = get_material_for_element(model, elem_ids[0] if elem_ids else 0)
        for eid in elem_ids:
            if eid in elem_to_cell:
                cell_idx = elem_to_cell[eid]
                if cell_idx < len(E_func.x.array):
                    E_func.x.array[cell_idx] = mat.E if mat.E > 0 else 210000
                    nu_func.x.array[cell_idx] = mat.nu if mat.nu > 0 else 0.3
    
    E_func.x.scatter_forward()
    nu_func.x.scatter_forward()
    
    if debug:
        print(f"Materiais aplicados: E min={E_func.x.array.min():.0f}, max={E_func.x.array.max():.0f}")

    # Propriedades como funções
    E = E_func
    nu = nu_func
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

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

    # Cargas do step
    for cl in step.loads:
        if cl.node not in model.nodes:
            continue
        n = model.nodes[cl.node]
        dist = np.linalg.norm(coords - np.array([n.x, n.y, n.z]), axis=1)
        node_idx = np.argmin(dist)
        if dist[node_idx] > 1e-6:
            continue
        comp = cl.dof - 1
        dof_comp = node_idx * 3 + comp
        with b.localForm() as local:
            if dof_comp < local.size:
                local[dof_comp] += cl.value

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

    os.makedirs(os.path.dirname(output_base) or ".", exist_ok=True)
    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, f"{output_base}.pvd", "w") as vtk:
            vtk.write_function(uh)
        if debug:
            print(f"PVD: {output_base}.pvd")
    except Exception as e:
        if debug:
            print(f"Aviso: {e}")

    norma = np.linalg.norm(uh.x.array)
    if debug:
        print(f"Norma: {norma:.6e}")
    return uh
    """Executa estática linear."""
    model = parse_inp(inp_path)
    domain, node_map = build_mesh_from_inp(model, debug=debug)
    
    domain.topology.create_connectivity(0, domain.topology.dim)
    domain.topology.create_connectivity(domain.topology.dim, 0)
    
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

    mat = next(iter(model.materials.values()))
    E = default_scalar_type(mat.E)
    nu = default_scalar_type(mat.nu)
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))

    if debug:
        print(f"Material: E={mat.E}, nu={mat.nu}")

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
    
    if debug:
        print(f"BCs aplicadas: {len(bcs)}")

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

    # CLOAD
    if model.cloads:
        for cl in model.cloads:
            if cl.node not in model.nodes:
                continue
            n = model.nodes[cl.node]
            dist = np.linalg.norm(coords - np.array([n.x, n.y, n.z]), axis=1)
            node_idx = np.argmin(dist)
            if dist[node_idx] > 1e-6:
                continue
            comp = cl.dof - 1
            dof_comp = node_idx * 3 + comp
            with b.localForm() as local:
                if dof_comp < local.size:
                    local[dof_comp] += cl.value
        if debug:
            print(f"CLOADs: {len(model.cloads)}")

    # DLOAD
    if model.dloads:
        for dl in model.dloads:
            for elset_name, elem_ids in model.elsets.items():
                if elset_name.upper() == dl.elset.upper():
                    for eid in elem_ids:
                        if eid in model.elements:
                            for nid in model.elements[eid].nodes:
                                if nid in model.nodes:
                                    n = model.nodes[nid]
                                    dist = np.linalg.norm(coords - np.array([n.x, n.y, n.z]), axis=1)
                                    node_idx = np.argmin(dist)
                                    if dist[node_idx] < 1e-6:
                                        dof_z = node_idx * 3 + 2
                                        with b.localForm() as local:
                                            if dof_z < local.size:
                                                local[dof_z] += dl.value * 0.05
                    break
        if debug:
            print(f"DLOADs: {len(model.dloads)}")

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
    os.makedirs(os.path.dirname(output_base) or ".", exist_ok=True)
    try:
        from dolfinx.io import VTKFile
        with VTKFile(domain.comm, f"{output_base}.pvd", "w") as vtk:
            vtk.write_function(uh)
        if debug:
            print(f"PVD: {output_base}.pvd")
    except Exception as e:
        if debug:
            print(f"Aviso: {e}")

    norma = np.linalg.norm(uh.x.array)
    if debug:
        print(f"Norma: {norma:.6e}")
    return uh


def run_from_inp(inp_path: str, output_dir: str = "/app/output", debug: int = 0):
    """Ponto de entrada principal. Executa cada step separadamente."""
    if not os.path.exists(inp_path):
        print(f"ERRO: arquivo não encontrado: {inp_path}")
        sys.exit(1)
    
    model = parse_inp(inp_path, debug=debug)
    if debug:
        print(summary(model))
        print("---")

    base_name = os.path.splitext(os.path.basename(inp_path))[0]
    
    # Se não houver steps, cria um padrão
    if not model.steps:
        model.steps.append(Step())
    
    results = []
    for i, step in enumerate(model.steps):
        step_num = i + 1
        output_base = os.path.join(output_dir, f"{base_name}_{step_num}")
        
        if debug:
            print(f"\n=== Step {step_num} ===")
        
        if model.analysis_type == "heat":
            from run_thermal import run_thermal
            result = run_thermal(inp_path, output_base, debug)
        elif step.nlgeom or model.nlgeom:
            from run_hyper import run_hyperelastic
            result = run_hyperelastic(model, step, output_base, debug)
        else:
            result = run_linear_static_step(model, step, output_base, debug)
        results.append(result)
    
    return results


if __name__ == "__main__":
    debug = 0
    args = sys.argv[1:]
    inp_file = None
    for arg in args:
        if arg.startswith("-debug="):
            debug = int(arg.split("=")[1])
        elif not arg.startswith("-"):
            inp_file = arg
    if not inp_file:
        print("Uso: python3 run_from_inp.py arquivo.inp [-debug=1]")
        sys.exit(1)
    run_from_inp(inp_file, debug=debug)
