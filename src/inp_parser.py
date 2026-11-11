#!/usr/bin/env python3
"""
Parser minimo de arquivo .inp (CalculiX/Abaqus) para FeniCS-EDSS.

Fora do perfil FC o recorte anterior e preservado (C3D8/C3D4, hiperelastico,
termico, etc.). No perfil FC (--profile fc) a whitelist e fechada: ver fc_profile.py.
"""

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from fc_profile import (
    EXPLICITLY_IGNORED_KEYWORDS,
    EXPLICITLY_IGNORED_PARAMS,
    FCParseRecord,
    FCProfileError,
    IgnoredItem,
    OutputRequest,
    keyword_params,
    normalize_keyword,
    validate_fc_cards,
    parse_report_payload,
    PARSE_NOT_SOLVED_BANNER,
)


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
    # Neo-Hooke lido de *HYPERELASTIC (mu, lambda). None = não informado.
    mu: Optional[float] = None
    lmbda: Optional[float] = None
    hyperelastic: bool = False


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
    static_increments: Optional[List[float]] = None
    step_params: List[str] = field(default_factory=list)
    start_time: float = 0.0


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
    hyperelastic: bool = False
    fc_record: Optional[FCParseRecord] = None
    element_nset_from_node_card: Optional[str] = None



def _is_keyword_line(raw: str) -> bool:
    s = raw.strip()
    return s.startswith("*") and not s.startswith("**")



def _resolve_set_aliases(sets, aliases):
    """Expande referencias a outros sets (ex.: ELSET Eall -> Evolumes)."""
    pending = {k: list(v) for k, v in aliases.items()}
    for _ in range(16):
        progressed = False
        for name, refs in list(pending.items()):
            still = []
            for ref in refs:
                hit = None
                for key in sets:
                    if key.upper() == ref.upper():
                        hit = key
                        break
                if hit is not None and (hit.upper() != name.upper()):
                    sets.setdefault(name, []).extend(sets.get(hit, []))
                    progressed = True
                else:
                    still.append(ref)
            pending[name] = still
        if not progressed:
            break


def parse_inp(path: str, debug: int = 0, profile: Optional[str] = None) -> INPModel:
    """Le um arquivo .inp e retorna um INPModel."""
    model = INPModel()
    current_material = None
    current_step = None
    fc_record = FCParseRecord()
    model.fc_record = fc_record
    elset_aliases = {}
    nset_aliases = {}
    profile_fc = (profile or "").strip().lower() == "fc"


    with open(path, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("**"):
            i += 1
            continue

        if not line.startswith("*"):
            i += 1
            continue

        key = normalize_keyword(line)
        lineno = i + 1
        fc_record.keywords.append((lineno, line))
        upper = line.upper()

        if key == "NODE":
            nset_m = re.search(r"NSET\s*=\s*(\w+)", line, flags=re.IGNORECASE)
            nset_name = nset_m.group(1) if nset_m else None
            i += 1
            while i < len(lines) and not _is_keyword_line(lines[i]):
                parts = [p.strip() for p in lines[i].split(",")]
                if len(parts) >= 4 and not lines[i].strip().startswith("**"):
                    nid = int(parts[0])
                    model.nodes[nid] = Node(nid, float(parts[1]), float(parts[2]), float(parts[3]))
                    if nset_name:
                        model.nsets.setdefault(nset_name, []).append(nid)
                i += 1
            continue

        if key == "ELEMENT":
            match = re.search(r"TYPE\s*=\s*(\w+)", line, flags=re.IGNORECASE)
            if match:
                etype = match.group(1)
            elif profile_fc:
                raise FCProfileError(
                    "L{0}: *ELEMENT sem TYPE e rejeitado no perfil FC (sem default C3D8).".format(lineno)
                )
            else:
                etype = "C3D8"
            elset_m = re.search(r"ELSET\s*=\s*(\w+)", line, flags=re.IGNORECASE)
            elset_name = elset_m.group(1) if elset_m else None
            i += 1
            while i < len(lines) and not _is_keyword_line(lines[i]):
                parts = [p.strip() for p in lines[i].split(",") if p.strip()]
                if len(parts) >= 2 and not lines[i].strip().startswith("**"):
                    eid = int(parts[0])
                    nodes = [int(n) for n in parts[1:]]
                    model.elements[eid] = Element(eid, etype, nodes)
                    if elset_name:
                        model.elsets.setdefault(elset_name, []).append(eid)
                i += 1
            continue

        if key == "MATERIAL":
            match = re.search(r"NAME\s*=\s*(\w+)", line, flags=re.IGNORECASE)
            if match:
                current_material = match.group(1)
                model.materials[current_material] = Material(name=current_material)
            i += 1
            continue

        if key == "ELASTIC" and current_material:
            i += 1
            while i < len(lines) and not _is_keyword_line(lines[i]):
                s = lines[i].strip()
                if s and not s.startswith("**"):
                    parts = [p.strip() for p in s.split(",")]
                    if len(parts) >= 2:
                        model.materials[current_material].E = float(parts[0])
                        model.materials[current_material].nu = float(parts[1])
                    i += 1
                    break
                i += 1
            continue

        if key == "HYPERELASTIC" and current_material:
            mat = model.materials[current_material]
            mat.hyperelastic = True
            model.hyperelastic = True
            i += 1
            while i < len(lines) and not _is_keyword_line(lines[i]):
                data = lines[i].strip()
                if data and not data.startswith("**"):
                    parts = [p.strip() for p in data.split(",") if p.strip()]
                    if parts:
                        mat.mu = float(parts[0])
                    if len(parts) >= 2:
                        mat.lmbda = float(parts[1])
                    i += 1
                    break
                i += 1
            continue

        if key == "CONDUCTIVITY" and current_material:
            i += 1
            while i < len(lines) and not _is_keyword_line(lines[i]):
                s = lines[i].strip()
                if s and not s.startswith("**"):
                    parts = [p.strip() for p in s.split(",")]
                    if parts:
                        model.materials[current_material].k = float(parts[0])
                    i += 1
                    break
                i += 1
            continue

        if key == "SOLID SECTION":
            match_elset = re.search(r"ELSET\s*=\s*(\w+)", line, flags=re.IGNORECASE)
            match_mat = re.search(r"MATERIAL\s*=\s*(\w+)", line, flags=re.IGNORECASE)
            if match_elset and match_mat:
                model.solid_sections[match_elset.group(1)] = match_mat.group(1)
            i += 1
            continue

        if key == "NSET":
            match = re.search(r"NSET\s*=\s*(\w+)", line, flags=re.IGNORECASE)
            if match:
                nset_name = match.group(1)
                model.nsets.setdefault(nset_name, [])
                i += 1
                while i < len(lines) and not _is_keyword_line(lines[i]):
                    parts = [p.strip() for p in lines[i].split(",") if p.strip()]
                    for tok in parts:
                        if tok.startswith("**"):
                            continue
                        try:
                            model.nsets[nset_name].append(int(tok))
                        except ValueError:
                            nset_aliases.setdefault(nset_name, []).append(tok)
                    i += 1
                continue
            i += 1
            continue

        if key == "ELSET":
            match = re.search(r"ELSET\s*=\s*(\w+)", line, flags=re.IGNORECASE)
            if match:
                elset_name = match.group(1)
                model.elsets.setdefault(elset_name, [])
                i += 1
                while i < len(lines) and not _is_keyword_line(lines[i]):
                    parts = [p.strip() for p in lines[i].split(",") if p.strip()]
                    for tok in parts:
                        if tok.startswith("**"):
                            continue
                        try:
                            model.elsets[elset_name].append(int(tok))
                        except ValueError:
                            elset_aliases.setdefault(elset_name, []).append(tok)
                    i += 1
                continue
            i += 1
            continue

        if key == "DLOAD":
            i += 1
            while i < len(lines) and not _is_keyword_line(lines[i]):
                parts = [p.strip() for p in lines[i].split(",")]
                if len(parts) >= 3 and not lines[i].strip().startswith("**"):
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

        if key == "NLGEOM":
            model.nlgeom = True
            i += 1
            continue

        if key == "BOUNDARY":
            i += 1
            while i < len(lines) and not _is_keyword_line(lines[i]):
                raw_parts = [p.strip() for p in lines[i].split(",")]
                while raw_parts and raw_parts[-1] == "":
                    raw_parts.pop()
                parts = raw_parts
                if len(parts) >= 2 and not lines[i].strip().startswith("**"):
                    target = parts[0]
                    if profile_fc and len(parts) == 3 and parts[2] != "":
                        third = parts[2]
                        looks_value = ("." in third) or ("e" in third.lower()) or (
                            third.lstrip("+-").isdigit() and int(third) not in range(1, 7)
                        )
                        if looks_value:
                            raise FCProfileError(
                                "L{0}: *BOUNDARY com 3 colunas e ambiguo no perfil FC. "
                                "CalculiX usa no, dof_ini, dof_fim[, valor]. "
                                "Use 2 colunas (valor 0) ou 4 colunas (no, ini, fim, valor). "
                                "Linha: {1}".format(i + 1, lines[i].strip())
                            )
                    try:
                        dof_start = int(parts[1])
                        if len(parts) >= 3 and parts[2] != "":
                            dof_end = int(parts[2])
                        else:
                            dof_end = dof_start
                        value = float(parts[3]) if len(parts) > 3 and parts[3] else 0.0
                    except ValueError as exc:
                        raise FCProfileError(
                            "L{0}: *BOUNDARY invalido no perfil FC ({1}). "
                            "Formato: no, dof_ini[, dof_fim[, valor]]. Linha: {2}".format(
                                i + 1, exc, lines[i].strip()
                            )
                        ) from exc
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

        if key == "STEP":
            current_step = Step()
            current_step.step_params = keyword_params(line)
            current_step.start_time = 0.0
            for param in current_step.step_params:
                pname = param.split("=", 1)[0].strip().upper()
                if (key, pname) in EXPLICITLY_IGNORED_PARAMS:
                    fc_record.ignored.append(
                        IgnoredItem(
                            "STEP",
                            "parametro INC ignorado (perfil FC: um incremento, tempo inicial 0, sem heranca)",
                            lineno,
                            param,
                        )
                    )
            model.steps.append(current_step)
            fc_record.step_params = list(current_step.step_params)
            i += 1
            continue

        if key == "END STEP":
            current_step = None
            i += 1
            continue

        if key == "STATIC":
            if "NLGEOM" in upper:
                model.nlgeom = True
                if current_step:
                    current_step.nlgeom = True
            model.analysis_type = "static"
            if current_step:
                current_step.analysis_type = "static"
            i += 1
            while i < len(lines) and not _is_keyword_line(lines[i]):
                data = lines[i].strip()
                if data and not data.startswith("**"):
                    try:
                        inc = [float(x.strip()) for x in data.split(",") if x.strip()]
                    except ValueError:
                        inc = None
                    if inc:
                        if current_step is not None:
                            current_step.static_increments = inc
                        fc_record.static_increments = inc
                    i += 1
                    break
                i += 1
            continue

        if key == "HEAT TRANSFER":
            model.analysis_type = "heat"
            if current_step:
                current_step.analysis_type = "heat"
            i += 1
            continue

        if key == "CLOAD":
            i += 1
            while i < len(lines) and not _is_keyword_line(lines[i]):
                parts = [p.strip() for p in lines[i].split(",")]
                if len(parts) >= 3 and not lines[i].strip().startswith("**"):
                    node = int(parts[0])
                    dof = int(parts[1])
                    value = float(parts[2])
                    cl = ConcentratedLoad(node, dof, value)
                    if current_step:
                        current_step.loads.append(cl)
                    else:
                        if not model.steps:
                            model.steps.append(Step())
                        model.steps[0].loads.append(cl)
                i += 1
            continue

        if key == "NODE FILE":
            req_line = lineno
            i += 1
            fields = []
            while i < len(lines) and not _is_keyword_line(lines[i]):
                s = lines[i].strip()
                if s and not s.startswith("**"):
                    fields.extend([x.strip().upper() for x in s.split(",") if x.strip()])
                i += 1
            fc_record.requests.append(OutputRequest("node_file", fields, line=req_line))
            continue

        if key == "EL FILE":
            req_line = lineno
            i += 1
            fields = []
            while i < len(lines) and not _is_keyword_line(lines[i]):
                s = lines[i].strip()
                if s and not s.startswith("**"):
                    fields.extend([x.strip().upper() for x in s.split(",") if x.strip()])
                i += 1
            fc_record.requests.append(OutputRequest("el_file", fields, line=req_line))
            continue

        if key == "NODE PRINT":
            req_line = lineno
            nset_m = re.search(r"NSET\s*=\s*(\w+)", line, flags=re.IGNORECASE)
            nset_name = nset_m.group(1) if nset_m else None
            totals = "TOTALS" in upper and "ONLY" in upper
            i += 1
            fields = []
            while i < len(lines) and not _is_keyword_line(lines[i]):
                s = lines[i].strip()
                if s and not s.startswith("**"):
                    fields.extend([x.strip().upper() for x in s.split(",") if x.strip()])
                i += 1
            fc_record.requests.append(
                OutputRequest("node_print", fields, nset=nset_name, totals_only=totals, line=req_line)
            )
            continue

        if key == "OUTPUT":
            # Parametros sao validados em validate_fc_cards: so FREQUENCY=1
            # e honrado. Reconhecer o cartao nao implica saida honrada.
            i += 1
            while i < len(lines) and not _is_keyword_line(lines[i]):
                s = lines[i].strip()
                if s and not s.startswith("**") and profile_fc:
                    raise FCProfileError(
                        "L{0}: dados sob *OUTPUT nao sao honrados no perfil FC: {1}".format(
                            i + 1, s
                        )
                    )
                i += 1
            continue

        if key in EXPLICITLY_IGNORED_KEYWORDS:
            if key == "EXPANSION":
                reason = (
                    "EXPANSION ignorado somente no perfil estritamente mecanico, "
                    "sem temperatura: o coeficiente nao entra na formulacao. "
                    "Cartao termico concomitante e rejeitado."
                )
            else:
                reason = (
                    "metadado ignorado explicitamente (nao altera a estatica "
                    "mecanica do perfil FC)"
                )
            fc_record.ignored.append(
                IgnoredItem(
                    key,
                    reason,
                    lineno,
                    line,
                )
            )
            i += 1
            while i < len(lines) and not _is_keyword_line(lines[i]):
                i += 1
            continue

        # Keyword nao despachada: avanca o cartao e os dados (sem interpretar).
        i += 1
        while i < len(lines) and not _is_keyword_line(lines[i]):
            i += 1
        continue

    _resolve_set_aliases(model.elsets, elset_aliases)
    _resolve_set_aliases(model.nsets, nset_aliases)

    if profile_fc:
        types = {el.type for el in model.elements.values()}
        validate_fc_cards(fc_record, types, len(model.steps))

    return model


def select_analysis(model: INPModel, step: Step) -> str:
    """Escolhe o solver do step: ``heat``, ``hyper`` ou ``linear``."""
    if step.analysis_type == "heat" or model.analysis_type == "heat":
        return "heat"
    if step.nlgeom or model.nlgeom or model.hyperelastic:
        return "hyper"
    return "linear"


def get_material_for_element(model: INPModel, elem_id: int) -> Optional[Material]:
    """Retorna o material do elemento via elset/solid section.

    Se o elemento não estiver em nenhum elset com material, devolve o
    primeiro material do modelo (ou ``None`` se não houver materiais).
    """
    for elset_name, elem_ids in model.elsets.items():
        if elem_id in elem_ids:
            mat_name = None
            for sname, mname in model.solid_sections.items():
                if sname.upper() == elset_name.upper():
                    mat_name = mname
                    break
            if not mat_name:
                continue
            for name, mat in model.materials.items():
                if name.upper() == mat_name.upper():
                    return mat
    if model.materials:
        return next(iter(model.materials.values()))
    return None


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
        f"Hiperelástico: {model.hyperelastic}",
    ]
    for i, step in enumerate(model.steps):
        lines.append(f"  Step {i+1}: {len(step.loads)} cargas, {len(step.dloads)} dloads, nlgeom={step.nlgeom}")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys
    args = [a for a in sys.argv[1:] if a]
    profile = None
    if "--profile" in args:
        ix = args.index("--profile")
        if ix + 1 >= len(args):
            print("faltou o valor de --profile")
            sys.exit(2)
        profile = args[ix + 1]
        del args[ix:ix + 2]
    if args and args[0].startswith("--profile="):
        profile = args[0].split("=", 1)[1]
        args = args[1:]
    if len(args) < 1:
        print("Uso: python3 inp_parser.py arquivo.inp [--profile fc]")
        sys.exit(1)
    try:
        model = parse_inp(args[0], profile=profile)
    except FCProfileError as exc:
        print("ERRO FC: {0}".format(exc))
        sys.exit(2)
    print(summary(model))
    if (profile or "").strip().lower() == "fc":
        print(PARSE_NOT_SOLVED_BANNER)
        rec = model.fc_record
        if rec:
            print("Keywords:", len(rec.keywords))
            print("Ignorados explicitos:", [(i.keyword, i.reason) for i in rec.ignored])
            print("Pedidos de saida:", [(r.kind, r.fields, r.nset) for r in rec.requests])
    elif model.fc_record:
        rec = model.fc_record
        print("Keywords:", len(rec.keywords))
        print("Ignorados explicitos:", [(i.keyword, i.reason) for i in rec.ignored])
        print("Pedidos de saida:", [(r.kind, r.fields, r.nset) for r in rec.requests])
