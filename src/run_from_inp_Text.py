#!/usr/bin/env python3
"""
Gera e imprime a entrada que seria enviada ao DOLFINx a partir de um .inp.
Não executa o solver — apenas mostra os dados.
"""

import sys
import numpy as np

sys.path.insert(0, "/app/src")
from inp_parser import parse_inp, summary


def generate_dolfinx_input(inp_path: str):
    model = parse_inp(inp_path)
    print("=" * 60)
    print("ENTRADA PARA DOLFINx")
    print("=" * 60)
    print()
    print(summary(model))
    print()

    # Nós
    print("--- NÓS ---")
    node_ids = sorted(model.nodes.keys())
    for nid in node_ids:
        n = model.nodes[nid]
        print(f"  {nid}: ({n.x}, {n.y}, {n.z})")
    print()

    # Elementos
    print("--- ELEMENTOS ---")
    for eid, el in model.elements.items():
        print(f"  {eid}: {el.type} -> {el.nodes}")
    print()

    # Material
    print("--- MATERIAL ---")
    for name, mat in model.materials.items():
        print(f"  {name}: E={mat.E}, nu={mat.nu}")
    print()

    # Nós fixos (após expansão de nsets)
    print("--- NÓS FIXOS (após expansão) ---")
    fixed = set(bc.node for bc in model.boundary)
    for nid in sorted(fixed):
        n = model.nodes.get(nid)
        if n:
            print(f"  Nó {nid}: ({n.x}, {n.y}, {n.z})")
    print(f"  Total: {len(fixed)} nós")
    print()

    # Node sets
    if model.nsets:
        print("--- NODE SETS ---")
        for name, nodes in model.nsets.items():
            print(f"  {name}: {nodes}")
        print()

    # Cargas
    print("--- CARGAS CONCENTRADAS ---")
    for cl in model.cloads:
        n = model.nodes.get(cl.node)
        coords = f"({n.x}, {n.y}, {n.z})" if n else "?"
        print(f"  Nó {cl.node} {coords}, DOF {cl.dof}, valor={cl.value}")
    print()

    # Coordenadas para malha
    print("--- COORDENADAS (array) ---")
    coords = np.array([[model.nodes[nid].x, model.nodes[nid].y, model.nodes[nid].z]
                       for nid in node_ids])
    print(coords)
    print()

    # Conectividade (com reordenação)
    print("--- CONECTIVIDADE REORDENADA (DOLFINx) ---")
    node_map = {nid: i for i, nid in enumerate(node_ids)}
    cells = []
    for el in model.elements.values():
        nodes = [node_map[n] for n in el.nodes]
        if el.type.upper() == "C3D8" and len(nodes) == 8:
            nodes = [nodes[0], nodes[1], nodes[3], nodes[2],
                     nodes[4], nodes[5], nodes[7], nodes[6]]
        cells.append(nodes)
    print(np.array(cells))
    print()

    # Verifica coordenadas mínimas/máximas
    coords = np.array([[model.nodes[nid].x, model.nodes[nid].y, model.nodes[nid].z]
                       for nid in node_ids])
    print("--- LIMITES DA MALHA ---")
    print(f"  X: {coords[:,0].min():.3f} a {coords[:,0].max():.3f}")
    print(f"  Y: {coords[:,1].min():.3f} a {coords[:,1].max():.3f}")
    print(f"  Z: {coords[:,2].min():.3f} a {coords[:,2].max():.3f}")
    print()

    # Nós na face x=0
    print("--- NÓS NA FACE X=0 ---")
    for i, nid in enumerate(node_ids):
        n = model.nodes[nid]
        if abs(n.x) < 1e-6:
            print(f"  Índice {i}, Nó {nid}: ({n.x}, {n.y}, {n.z})")
    print()

    print("=" * 60)
    print("FIM DA ENTRADA")
    print("=" * 60)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python3 run_from_inp_Text.py arquivo.inp")
        sys.exit(1)
    generate_dolfinx_input(sys.argv[1])
