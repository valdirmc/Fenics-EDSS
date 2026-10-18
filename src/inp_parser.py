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
class INPModel:
    nodes: Dict[int, Node] = field(default_factory=dict)
    elements: Dict[int, Element] = field(default_factory=dict)
    materials: Dict[str, Material] = field(default_factory=dict)
    solid_sections: Dict[str, str] = field(default_factory=dict)  # elset -> material
    boundary: List[BoundaryCondition] = field(default_factory=list)
    cloads: List[ConcentratedLoad] = field(default_factory=list)
    analysis_type: str = "static"


def parse_inp(path: str) -> INPModel:
    """Lê um arquivo .inp e retorna um INPModel."""
    model = INPModel()
    current_material = None
    current_section_elset = None

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

        if upper.startswith("*SOLID SECTION"):
            match_elset = re.search(r"ELSET\s*=\s*(\w+)", upper)
            match_mat = re.search(r"MATERIAL\s*=\s*(\w+)", upper)
            if match_elset and match_mat:
                model.solid_sections[match_elset.group(1)] = match_mat.group(1)
            i += 1
            continue

        if upper.startswith("*BOUNDARY"):
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("*"):
                parts = [p.strip() for p in lines[i].split(",")]
                if len(parts) >= 4:
                    node = int(parts[0])
                    dof = int(parts[1])
                    value = float(parts[3]) if parts[3] else 0.0
                    model.boundary.append(BoundaryCondition(node, dof, value))
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
                    model.cloads.append(ConcentratedLoad(node, dof, value))
                i += 1
            continue

        if upper.startswith("*STATIC"):
            model.analysis_type = "static"
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
        f"Condições de contorno: {len(model.boundary)}",
        f"Cargas concentradas: {len(model.cloads)}",
        f"Tipo de análise: {model.analysis_type}",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Uso: python3 inp_parser.py arquivo.inp")
        sys.exit(1)
    model = parse_inp(sys.argv[1])
    print(summary(model))
