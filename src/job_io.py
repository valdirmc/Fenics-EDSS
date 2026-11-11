#!/usr/bin/env python3
"""
Contrato de entrada/saída por job do FeniCS-EDSS.

O FemEDSS-app isola cada job em
``uploads/<CNPJ>/<id do usuário>/<id do job>/{input,output}``.
Este módulo interpreta a linha de comando e garante que os arquivos
``.pvd`` / ``.pvtu`` / ``.vtu`` fiquem somente no diretório pedido,
com referências relativas válidas dentro dessa pasta.
"""

import hashlib
import json
import os
import re
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

# Compatível com chamadas antigas que não passam diretório de saída.
DEFAULT_OUTPUT_DIR = "/app/output"

USAGE = (
    "Uso: python3 run_from_inp.py arquivo.inp "
    "[--output-dir DIR | -output DIR] [--profile fc] "
    "[--units mm-N-MPa-t-s] [-debug=1]\n"
    "  arquivo.inp            caminho do .inp deste job\n"
    "  --output-dir, -output  diretorio de saida deste job "
    f"(padrao: {DEFAULT_OUTPUT_DIR})\n"
    "  --profile fc           whitelist fechada e contrato de saidas FC\n"
    "  --units                obrigatorio com --profile fc (mm-N-MPa-t-s)\n"
    "  -debug=1               imprime diagnostico"
)

_FILE_ATTR = re.compile(
    r'(\b(?:file|Source)\s*=\s*")([^"]*)(")',
    re.IGNORECASE,
)


class RunArgs:
    """Argumentos da CLI. Desempacota como (inp, output_dir, debug) para testes antigos."""

    def __init__(
        self,
        inp_file: Optional[str],
        output_dir: Optional[str],
        debug: int,
        profile: Optional[str] = None,
        units: Optional[str] = None,
    ):
        self.inp_file = inp_file
        self.output_dir = output_dir
        self.debug = debug
        self.profile = profile
        self.units = units

    def __iter__(self):
        yield self.inp_file
        yield self.output_dir
        yield self.debug

    def __getitem__(self, index):
        return (self.inp_file, self.output_dir, self.debug)[index]


def parse_run_args(argv: Sequence[str]) -> RunArgs:
    """Interpreta a CLI.

    Retorna RunArgs desempacotavel como ``(inp, output_dir, debug)``.
    ``output_dir`` e ``None`` quando a opcao nao foi informada.
    """
    debug = 0
    output_dir: Optional[str] = None
    inp_file: Optional[str] = None
    profile: Optional[str] = None
    units: Optional[str] = None
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
                raise ValueError(f"faltou o diretorio depois de {arg}")
            output_dir = args[i]
        elif arg.startswith("--output-dir="):
            output_dir = arg.split("=", 1)[1]
        elif arg.startswith("-output="):
            output_dir = arg.split("=", 1)[1]
        elif arg.startswith("--output="):
            output_dir = arg.split("=", 1)[1]
        elif arg in ("--profile",):
            i += 1
            if i >= len(args):
                raise ValueError("faltou o valor de --profile")
            profile = args[i]
        elif arg.startswith("--profile="):
            profile = arg.split("=", 1)[1]
        elif arg in ("--units",):
            i += 1
            if i >= len(args):
                raise ValueError("faltou o valor de --units")
            units = args[i]
        elif arg.startswith("--units="):
            units = arg.split("=", 1)[1]
        elif arg.startswith("-"):
            raise ValueError(f"opcao desconhecida: {arg}")
        else:
            if inp_file is not None:
                raise ValueError("informe apenas um arquivo .inp")
            inp_file = arg
        i += 1
    if output_dir is not None and not str(output_dir).strip():
        raise ValueError("diretorio de saida vazio")
    if profile is not None:
        profile = str(profile).strip().lower()
        if profile not in ("fc",):
            raise ValueError("perfil desconhecido: {0} (use fc)".format(profile))
    return RunArgs(inp_file, output_dir, debug, profile=profile, units=units)


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
        raise RuntimeError("Falha ao gravar VTK {0}: {1}".format(pvd_path, exc)) from exc
    if not os.path.isfile(pvd_path):
        raise RuntimeError("Falha ao gravar VTK: arquivo ausente {0}".format(pvd_path))
    return pvd_path


def file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_run_manifest(output_dir: str, payload: Dict[str, Any]) -> str:
    """Grava run_manifest.json. Sucesso so conta artefatos listados neste run."""
    directory = os.path.abspath(output_dir)
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, "run_manifest.json")
    artifacts = []
    for raw in payload.get("artifacts", []):
        item = dict(raw)
        fpath = item.get("path")
        if fpath and os.path.isfile(fpath):
            item["sha256"] = file_sha256(fpath)
            item["bytes"] = os.path.getsize(fpath)
            item["mtime"] = os.path.getmtime(fpath)
        artifacts.append(item)
    payload = dict(payload)
    payload["artifacts"] = artifacts
    payload["written_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return path


def verify_run_artifacts(manifest_path: str, required_names: Sequence[str]) -> None:
    """Nao aceita lixo velho: so arquivos deste manifesto, com hash conferido."""
    with open(manifest_path, encoding="utf-8") as handle:
        payload = json.load(handle)
    by_name = {os.path.basename(a.get("path", "")): a for a in payload.get("artifacts", [])}
    missing = []
    stale = []
    for name in required_names:
        item = by_name.get(name)
        if not item:
            missing.append(name)
            continue
        path = item.get("path")
        expected = item.get("sha256")
        if not path or not os.path.isfile(path) or not expected:
            missing.append(name)
            continue
        if file_sha256(path) != expected:
            stale.append(name)
    if missing or stale:
        raise RuntimeError(
            "Artefatos desta run ausentes ou nao conferem com o manifesto "
            "(lixo velho nao conta como sucesso). ausentes={0} divergentes={1}".format(
                missing, stale
            )
        )
