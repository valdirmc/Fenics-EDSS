"""Perfil FC (First Customer): whitelist fechada e contrato de saidas.

Opt-in via --profile fc. Fora deste perfil o comportamento anterior
e preservado (exemplos termicos/hiperelasticos e cards nao listadas).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Sequence, Set, Tuple


PROFILE_NAME = "fc"

FC_UNITS_ALIASES = {
    "mm-n-mpa-t-s": "mm-N-MPa-t-s",
    "mm-n-s-t": "mm-N-MPa-t-s",
    "mm-n-mpa-t": "mm-N-MPa-t-s",
}

ALLOWED_KEYWORDS: Set[str] = {
    "NODE",
    "ELEMENT",
    "NSET",
    "ELSET",
    "MATERIAL",
    "ELASTIC",
    "SOLID SECTION",
    "BOUNDARY",
    "CLOAD",
    "STEP",
    "STATIC",
    "END STEP",
    "NODE FILE",
    "EL FILE",
    "NODE PRINT",
    "EL PRINT",
    "OUTPUT",
}

EXPLICITLY_IGNORED_KEYWORDS: Set[str] = {
    "PHYSICAL CONSTANTS",
    "EXPANSION",
    "DENSITY",
}

EXPLICITLY_IGNORED_PARAMS = {
    ("STEP", "INC"),
}

BLOCKED_KEYWORDS: Set[str] = {
    "DLOAD",
    "NLGEOM",
    "HYPERELASTIC",
    "HEAT TRANSFER",
    "CONDUCTIVITY",
    "PLASTIC",
    "FREQUENCY",
    "DYNAMIC",
    "MODAL DYNAMIC",
    "STEADY STATE DYNAMICS",
    "CONTACT PAIR",
    "SURFACE INTERACTION",
    "FRICTION",
    "TIE",
    "AMPLITUDE",
    "INCLUDE",
    "TEMPERATURE",
    "INITIAL CONDITIONS",
    "CFLUX",
    "DFLUX",
    "FILM",
    "RADIATE",
    "SPECIFIC HEAT",
}

# Cartões térmicos: *EXPANSION só pode ser ignorado se nenhum destes aparecer.
THERMAL_PHYSICS_KEYWORDS: Set[str] = {
    "TEMPERATURE",
    "INITIAL CONDITIONS",
    "HEAT TRANSFER",
    "CONDUCTIVITY",
    "CFLUX",
    "DFLUX",
    "FILM",
    "RADIATE",
    "SPECIFIC HEAT",
}

PARSE_NOT_SOLVED_BANNER = (
    "STATUS: parse_ok_not_solved\n"
    "AVISO: relatorio de parse apenas. Nao implica solver valido, "
    "equilibrio, honra de saidas na execucao, nem artefatos VTK."
)

REQUIRED_FIELDS_ALWAYS = ("U", "S", "PRINCIPAL", "VON_MISES")

HONORABLE_OUTPUT_FIELDS = {
    "U": "displacement (nodal, CG1 / P1)",
    "S": "Cauchy stress (celula, DG0 / constante no C3D4 P1)",
    "E": "infinitesimal strain (celula, DG0 / constante no C3D4 P1)",
    "RF": "reaction force (nos do NSET pedido e totais)",
}


class FCProfileError(ValueError):
    """Rejeicao explicita do perfil FC (exit code 2 no CLI)."""


def normalize_units(raw: Optional[str]) -> Optional[str]:
    if raw is None or not str(raw).strip():
        return None
    key = str(raw).strip().lower()
    if key in FC_UNITS_ALIASES:
        return FC_UNITS_ALIASES[key]
    raise FCProfileError(
        "Sistema de unidades nao reconhecido: {0!r}. "
        "Use mm-N-MPa-t-s (alias: mm-N-s-t).".format(raw)
    )


def normalize_keyword(header: str) -> str:
    text = header.strip()
    if text.startswith("*"):
        text = text[1:]
    return text.split(",")[0].strip().upper()


def keyword_params(header: str) -> List[str]:
    text = header.strip()
    if text.startswith("*"):
        text = text[1:]
    return [p.strip() for p in text.split(",")[1:] if p.strip()]


@dataclass
class OutputRequest:
    kind: str
    fields: List[str]
    nset: Optional[str] = None
    elset: Optional[str] = None
    totals_only: bool = False
    line: int = 0


@dataclass
class IgnoredItem:
    keyword: str
    reason: str
    line: int = 0
    detail: str = ""


@dataclass
class FCParseRecord:
    keywords: List[Tuple[int, str]] = field(default_factory=list)
    ignored: List[IgnoredItem] = field(default_factory=list)
    requests: List[OutputRequest] = field(default_factory=list)
    static_increments: Optional[List[float]] = None
    step_params: List[str] = field(default_factory=list)


def validate_fc_cards(record: FCParseRecord, element_types: Iterable[str], n_steps: int) -> None:
    errors: List[str] = []

    for lineno, header in record.keywords:
        key = normalize_keyword(header)
        if key in EXPLICITLY_IGNORED_KEYWORDS or key in ALLOWED_KEYWORDS:
            continue
        if key in BLOCKED_KEYWORDS:
            errors.append(
                "L{0}: *{1} e bloqueado no perfil FC (sem fallback). Cartao: {2}".format(
                    lineno, key, header.strip()
                )
            )
        elif key:
            errors.append(
                "L{0}: keyword nao esta na whitelist do perfil FC: {1}".format(
                    lineno, header.strip()
                )
            )

    if n_steps != 1:
        errors.append(
            "Perfil FC exige exatamente um *STEP (demo um step; "
            "multistep e rejeitado). Encontrados: {0}.".format(n_steps)
        )

    types = {t.upper() for t in element_types if t}
    if not types:
        errors.append("Perfil FC exige malha com *ELEMENT, TYPE=C3D4.")
    elif types != {"C3D4"}:
        errors.append(
            "Perfil FC aceita somente C3D4 (tetra linear de 4 nos), "
            "sem fallback C3D8/C3D10. Tipos encontrados: {0}.".format(sorted(types))
        )

    seen = [normalize_keyword(h) for _, h in record.keywords]
    if "STATIC" not in seen:
        errors.append("Perfil FC exige um procedimento *STATIC (sem NLGEOM).")
    if any(normalize_keyword(h) == "STATIC" and "NLGEOM" in h.upper() for _, h in record.keywords):
        errors.append("*STATIC, NLGEOM e bloqueado no perfil FC.")

    if record.static_increments:
        inc = record.static_increments
        initial = inc[0]
        period = inc[1] if len(inc) > 1 else 1.0
        if abs(period - 1.0) > 1e-12 or initial > period + 1e-12:
            errors.append(
                "Perfil FC nao implementa incrementacao de *STATIC alem de "
                "um incremento com periodo 1 (lido: {0}).".format(inc)
            )

    seen = [normalize_keyword(h) for _, h in record.keywords]
    if "EXPANSION" in seen and any(k in THERMAL_PHYSICS_KEYWORDS for k in seen):
        errors.append(
            "*EXPANSION so pode ser ignorado num perfil estritamente mecanico, "
            "sem temperatura/fluxo termico. Ha cartao termico no .inp; "
            "ignorar a dilatacao perderia fisica."
        )

    for lineno, header in record.keywords:
        if normalize_keyword(header) != "OUTPUT":
            continue
        params = keyword_params(header)
        for param in params:
            compact = param.upper().replace(" ", "")
            if compact in {"FREQUENCY=1", "FREQ=1"}:
                continue
            errors.append(
                "L{0}: *OUTPUT nao honrado (parametro {1}). "
                "Reconhecer o cartao nao basta; o perfil FC so aceita "
                "FREQUENCY=1 com um unico incremento.".format(lineno, param)
            )

    if errors:
        raise FCProfileError("Rejeicao do perfil FC:\n- " + "\n- ".join(errors))


def required_output_fields(record: FCParseRecord) -> List[str]:
    fields: List[str] = list(REQUIRED_FIELDS_ALWAYS)
    for req in record.requests:
        for name in req.fields:
            if name not in fields:
                fields.append(name)
    return fields


def check_honored_outputs(requested: Sequence[str], honored: Sequence[str]) -> None:
    missing = [f for f in requested if f not in set(honored)]
    if missing:
        raise FCProfileError(
            "Perfil FC recusa a run: pedidos de saida nao honrados: "
            + ", ".join(missing)
            + ". Nao basta o parser aceitar o .inp."
        )


def parse_report_payload(model, inp_path: str = "") -> dict:
    """Relatorio de parse. Nunca afirma que o solver e valido."""
    rec = getattr(model, "fc_record", None)
    return {
        "status": "parse_ok_not_solved",
        "solver_run": False,
        "profile": "fc",
        "inp": inp_path,
        "n_nodes": len(model.nodes),
        "n_elements": len(model.elements),
        "n_steps": len(model.steps),
        "element_types": sorted({el.type.upper() for el in model.elements.values()}),
        "ignored_cards": [
            {"keyword": it.keyword, "reason": it.reason, "line": it.line, "detail": it.detail}
            for it in (rec.ignored if rec else [])
        ],
        "output_requests": [
            {"kind": r.kind, "fields": r.fields, "nset": r.nset, "totals_only": r.totals_only, "line": r.line}
            for r in (rec.requests if rec else [])
        ],
        "warning": (
            "Este relatorio e somente de parse. Nao significa que o solver, "
            "o equilibrio ou os artefatos VTK sejam validos."
        ),
    }
