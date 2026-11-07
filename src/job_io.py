#!/usr/bin/env python3
"""
Contrato de entrada/saída por job do FeniCS-EDSS.

O FemEDSS-app isola cada job em
``uploads/<CNPJ>/<id do usuário>/<id do job>/{input,output}``.
Este módulo interpreta a linha de comando e garante que os arquivos
``.pvd`` / ``.pvtu`` / ``.vtu`` fiquem somente no diretório pedido,
com referências relativas válidas dentro dessa pasta.
"""

import os
import re
from typing import List, Optional, Sequence, Tuple

# Compatível com chamadas antigas que não passam diretório de saída.
DEFAULT_OUTPUT_DIR = "/app/output"

USAGE = (
    "Uso: python3 run_from_inp.py arquivo.inp "
    "[--output-dir DIR | -output DIR] [-debug=1]\n"
    "  arquivo.inp          caminho do .inp deste job\n"
    "  --output-dir, -output  diretório de saída deste job "
    f"(padrão: {DEFAULT_OUTPUT_DIR})\n"
    "  -debug=1             imprime diagnóstico"
)

_FILE_ATTR = re.compile(
    r'(\b(?:file|Source)\s*=\s*")([^"]*)(")',
    re.IGNORECASE,
)


def parse_run_args(argv: Sequence[str]) -> Tuple[Optional[str], Optional[str], int]:
    """Interpreta a CLI.

    Retorna ``(inp, output_dir, debug)``. ``output_dir`` é ``None`` quando
    a opção não foi informada (o chamador aplica o padrão).
    """
    debug = 0
    output_dir: Optional[str] = None
    inp_file: Optional[str] = None
    args: List[str] = list(argv)
    i = 0
    while i < len(args):
        arg = args[i]
        if arg.startswith("-debug="):
            debug = _parse_debug(arg.split("=", 1)[1])
        elif arg == "-debug":
            i += 1
            if i >= len(args):
                raise ValueError("faltou o valor de -debug")
            debug = _parse_debug(args[i])
        elif arg in ("--output-dir", "-output", "--output", "-o"):
            i += 1
            if i >= len(args):
                raise ValueError(f"faltou o diretório depois de {arg}")
            output_dir = args[i]
        elif arg.startswith("--output-dir="):
            output_dir = arg.split("=", 1)[1]
        elif arg.startswith("-output="):
            output_dir = arg.split("=", 1)[1]
        elif arg.startswith("--output="):
            output_dir = arg.split("=", 1)[1]
        elif arg.startswith("-"):
            raise ValueError(f"opção desconhecida: {arg}")
        else:
            if inp_file is not None:
                raise ValueError("informe apenas um arquivo .inp")
            inp_file = arg
        i += 1
    if output_dir is not None and not str(output_dir).strip():
        raise ValueError("diretório de saída vazio")
    return inp_file, output_dir, debug


def resolve_output_dir(output_dir: Optional[str]) -> str:
    """Aplica o diretório padrão quando a opção não foi passada."""
    if output_dir is None or not str(output_dir).strip():
        return DEFAULT_OUTPUT_DIR
    return str(output_dir)


def prepare_output_dir(output_dir: str) -> str:
    """Cria somente o diretório de saída do job. Não apaga nem lista outras pastas."""
    path = os.path.abspath(resolve_output_dir(output_dir))
    os.makedirs(path, exist_ok=True)
    return path


def step_output_base(output_dir: str, inp_path: str, step_num: int) -> str:
    """Caminho sem extensão de ``{nome}_{step}`` dentro do diretório do job."""
    directory = prepare_output_dir(output_dir)
    base_name = os.path.splitext(os.path.basename(inp_path))[0] or "resultado"
    return os.path.join(directory, f"{base_name}_{step_num}")


def _parse_debug(raw: str) -> int:
    try:
        return int(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"valor de debug inválido: {raw}") from exc


def _is_within(root: str, path: str) -> bool:
    root_real = os.path.realpath(root)
    path_real = os.path.realpath(path)
    return path_real == root_real or path_real.startswith(root_real + os.sep)


def _relative_ref(ref: str, out_dir: str, xml_dir: str) -> str:
    """Converte uma referência VTK para um caminho relativo ao XML, dentro de out_dir."""
    ref = ref.strip().replace("\\", "/")
    if not ref or "://" in ref:
        return ref
    basename = os.path.basename(ref)
    sibling = os.path.join(out_dir, basename)
    if basename and os.path.isfile(sibling) and _is_within(out_dir, sibling):
        return os.path.relpath(sibling, xml_dir).replace("\\", "/")
    if not os.path.isabs(ref):
        joined = os.path.normpath(os.path.join(xml_dir, ref))
        if _is_within(out_dir, joined):
            return ref
    return basename or ref


def _rewrite_xml_refs(path: str, out_dir: str) -> None:
    with open(path, encoding="utf-8", errors="replace") as handle:
        original = handle.read()

    xml_dir = os.path.dirname(os.path.abspath(path))

    def repl(match: re.Match) -> str:
        return f"{match.group(1)}{_relative_ref(match.group(2), out_dir, xml_dir)}{match.group(3)}"

    updated = _FILE_ATTR.sub(repl, original)
    if updated != original:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(updated)


def ensure_relative_vtk_links(pvd_path: str) -> None:
    """Ajusta ``file``/``Source`` do .pvd e dos .pvtu irmãos para caminhos relativos.

    Só abre o .pvd indicado e os .pvtu que ele cita, e somente se esses
    arquivos estiverem dentro do mesmo diretório. Não há varredura nem
    limpeza de outras pastas.
    """
    if not pvd_path or not os.path.isfile(pvd_path):
        return
    pvd_path = os.path.abspath(pvd_path)
    out_dir = os.path.dirname(pvd_path)
    _rewrite_xml_refs(pvd_path, out_dir)

    with open(pvd_path, encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    for ref in re.findall(r'\bfile\s*=\s*"([^"]+)"', text, flags=re.IGNORECASE):
        candidate = os.path.normpath(os.path.join(out_dir, ref))
        if (
            candidate.lower().endswith(".pvtu")
            and os.path.isfile(candidate)
            and _is_within(out_dir, candidate)
        ):
            _rewrite_xml_refs(candidate, out_dir)


def write_vtk_function(comm, function, output_base: str, debug: int = 0) -> str:
    """Grava ``{output_base}.pvd`` e os .vtu/.pvtu ao lado, com links relativos."""
    out_dir = os.path.dirname(os.path.abspath(output_base))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    pvd_path = f"{output_base}.pvd"
    try:
        from dolfinx.io import VTKFile

        with VTKFile(comm, pvd_path, "w") as vtk:
            vtk.write_function(function)
        ensure_relative_vtk_links(pvd_path)
        if debug:
            print(f"PVD: {pvd_path}")
    except Exception as exc:
        if debug:
            print(f"Aviso: {exc}")
    return pvd_path
