#!/usr/bin/env python3
"""
Parser mínimo de arquivo .inp (CalculiX/Abaqus) para FeniCS-EDSS.
Suporta: *NODE, *ELEMENT (C3D8, C3D4), *MATERIAL, *ELASTIC,
         *SOLID SECTION, *BOUNDARY, *CLOAD, *STATIC
"""

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class Node:
    id: int
    x: float
    y: float
    z: float


@dataclass
class Element:
    id: int
    type: str
    nodes: List[int]


@dataclass
class Material:
    name: str
    E: float = 0.0
    nu: float = 0.0
    k: float = 0.0  # condutividade térmica


@dataclass
class BoundaryCondition:
    node: int
    dof: int
    value: float


@dataclass
class ConcentratedLoad:
    node: int
    dof: int
    value: float


@dataclass
class SurfaceLoad:
    elset: str
    face: str
    value: float


@dataclass
class Step:
    loads: List[ConcentratedLoad] = field(default_factory=list)
    dloads: List[SurfaceLoad] = field(default_factory=list)
    nlgeom: bool = False
    analysis_type: str = "static"


@dataclass
class INPModel:
    nodes: Dict[int, Node] = field(default_factory=dict)
    elements: Dict[int, Element] = field(default_factory=dict)
    materials: Dict[str, Material] = field(default_factory=dict)
    solid_sections: Dict[str, str] = field(default_factory=dict)
    elsets: Dict[str, List[int]] = field(default_factory=dict)
    nsets: Dict[str, List[int]] = field(default_factory=dict)
    boundary: List[BoundaryCondition] = field(default_factory=list)
    steps: List[Step] = field(default_factory=list)
    analysis_type: str = "static"
    nlgeom: bool = False


def parse_inp(path: str, debug: int = 0) -> INPModel:
    """Lê um arquivo .inp e retorna um INPModel."""
    model = INPModel()
    current_material = None
    current_step = None

    with open(path, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("**"):
            i += 1
            continue

        upper = line.upper()

        if upper.startswith("*NODE"):
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("*"):
                parts = [p.strip() for p in lines[i].split(",")]
                if len(parts) >= 4:
                    nid = int(parts[0])
                    model.nodes[nid] = Node(nid, float(parts[1]), float(parts[2]), float(parts[3]))
                i += 1
            continue

        if upper.startswith("*ELEMENT"):
            match = re.search(r"TYPE\s*=\s*(\w+)", upper)
            etype = match.group(1) if match else "C3D8"
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("*"):
                parts = [p.strip() for p in lines[i].split(",") if p.strip()]
                if len(parts) >= 2:
                    eid = int(parts[0])
                    nodes = [int(n) for n in parts[1:]]
                    model.elements[eid] = Element(eid, etype, nodes)
                i += 1
            continue

        if upper.startswith("*MATERIAL"):
            match = re.search(r"NAME\s*=\s*(\w+)", upper)
            if match:
                current_material = match.group(1)
                model.materials[current_material] = Material(name=current_material)
            i += 1
            continue

        if upper.startswith("*ELASTIC") and current_material:
            i += 1
            if i < len(lines):
                parts = [p.strip() for p in lines[i].split(",")]
                if len(parts) >= 2:
                    model.materials[current_material].E = float(parts[0])
                    model.materials[current_material].nu = float(parts[1])
            i += 1
            continue

        if upper.startswith("*CONDUCTIVITY") and current_material:
            i += 1
            if i < len(lines):
                parts = [p.strip() for p in lines[i].split(",")]
                if parts:
                    model.materials[current_material].k = float(parts[0])
            i += 1
            continue

        if upper.startswith("*SOLID SECTION"):
            match_elset = re.search(r"ELSET\s*=\s*(\w+)", upper)
            match_mat = re.search(r"MATERIAL\s*=\s*(\w+)", upper)
            if match_elset and match_mat:
                model.solid_sections[match_elset.group(1)] = match_mat.group(1)
            i += 1
            continue

        if upper.startswith("*NSET"):
            match = re.search(r"NSET\s*=\s*(\w+)", upper)
            if match:
                nset_name = match.group(1)
                model.nsets[nset_name] = []
                i += 1
                while i < len(lines) and not lines[i].strip().startswith("*"):
                    parts = [p.strip() for p in lines[i].split(",") if p.strip()]
                    for p in parts:
                        try:
                            model.nsets[nset_name].append(int(p))
                        except ValueError:
                            pass
                    i += 1
            continue

        if upper.startswith("*ELSET"):
            match = re.search(r"ELSET\s*=\s*(\w+)", upper)
            if match:
                elset_name = match.group(1)
                model.elsets[elset_name] = []
                i += 1
                while i < len(lines) and not lines[i].strip().startswith("*"):
                    parts = [p.strip() for p in lines[i].split(",") if p.strip()]
                    for p in parts:
                        try:
                            model.elsets[elset_name].append(int(p))
                        except ValueError:
                            pass
                    i += 1
            continue

        if upper.startswith("*DLOAD"):
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("*"):
                parts = [p.strip() for p in lines[i].split(",")]
                if len(parts) >= 3:
                    elset = parts[0]
                    face = parts[1]
                    value = float(parts[2])
                    dl = SurfaceLoad(elset, face, value)
                    if current_step:
                        current_step.dloads.append(dl)
                    else:
                        if not model.steps:
                            model.steps.append(Step())
                        model.steps[0].dloads.append(dl)
                i += 1
            continue

        if upper.startswith("*NLGEOM"):
            model.nlgeom = True
            i += 1
            continue

        if upper.startswith("*BOUNDARY"):
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("*"):
                parts = [p.strip() for p in lines[i].split(",")]
                if len(parts) >= 3:
                    target = parts[0]
                    dof_start = int(parts[1])
                    dof_end = int(parts[2]) if len(parts) > 2 and parts[2] else dof_start
                    value = float(parts[3]) if len(parts) > 3 and parts[3] else 0.0
                    # Expande o range de DOFs
                    nodes_to_apply = []
                    if target.isdigit():
                        nodes_to_apply = [int(target)]
                    else:
                        for nset_name, nodes in model.nsets.items():
                            if nset_name.upper() == target.upper():
                                nodes_to_apply = nodes
                                break
                    for nid in nodes_to_apply:
                        for dof in range(dof_start, dof_end + 1):
                            model.boundary.append(BoundaryCondition(nid, dof, value))
                i += 1
            continue

        if upper.startswith("*STEP"):
            current_step = Step()
            model.steps.append(current_step)
            i += 1
            continue

        if upper.startswith("*END STEP"):
            current_step = None
            i += 1
            continue

        if upper.startswith("*STATIC"):
            if "NLGEOM" in upper:
                model.nlgeom = True
                if current_step:
                    current_step.nlgeom = True
            model.analysis_type = "static"
            if current_step:
                current_step.analysis_type = "static"
            i += 1
            continue

        if upper.startswith("*CLOAD"):
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("*"):
                parts = [p.strip() for p in lines[i].split(",")]
                if len(parts) >= 3:
                    node = int(parts[0])
                    dof = int(parts[1])
                    value = float(parts[2])
                    cl = ConcentratedLoad(node, dof, value)
                    if current_step:
                        current_step.loads.append(cl)
                    else:
                        # Fora de step, adiciona ao primeiro step ou cria um
                        if not model.steps:
                            model.steps.append(Step())
                        model.steps[0].loads.append(cl)
                i += 1
            continue

        i += 1

    return model


def summary(model: INPModel) -> str:
    """Retorna um resumo do modelo lido."""
    lines = [
        f"Nós: {len(model.nodes)}",
        f"Elementos: {len(model.elements)}",
        f"Materiais: {list(model.materials.keys())}",
        f"Seções sólidas: {model.solid_sections}",
        f"Element sets: {list(model.elsets.keys())}",
        f"Node sets: {list(model.nsets.keys())}",
        f"Condições de contorno: {len(model.boundary)}",
        f"Steps: {len(model.steps)}",
        f"Tipo de análise: {model.analysis_type}",
        f"NLGEOM: {model.nlgeom}",
    ]
    for i, step in enumerate(model.steps):
        lines.append(f"  Step {i+1}: {len(step.loads)} cargas, {len(step.dloads)} dloads, nlgeom={step.nlgeom}")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Uso: python3 inp_parser.py arquivo.inp")
        sys.exit(1)
    model = parse_inp(sys.argv[1])
    print(summary(model))
