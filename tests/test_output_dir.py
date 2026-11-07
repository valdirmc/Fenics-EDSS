"""Testes do diretório de saída por job, sem DOLFINx.

Rodar: python3 -m pytest tests/test_output_dir.py tests/test_parser_cards.py -v
"""

import os
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

from job_io import (  # noqa: E402
    DEFAULT_OUTPUT_DIR,
    ensure_relative_vtk_links,
    parse_run_args,
    prepare_output_dir,
    resolve_output_dir,
    step_output_base,
)


def test_cli_padrao_nao_cria_saida_global():
    antes = os.path.exists(DEFAULT_OUTPUT_DIR)
    inp, output_dir, debug = parse_run_args(["/job/input/caso.inp", "-debug=1"])
    assert inp == "/job/input/caso.inp"
    assert output_dir is None
    assert debug == 1
    assert resolve_output_dir(output_dir) == DEFAULT_OUTPUT_DIR
    assert os.path.exists(DEFAULT_OUTPUT_DIR) is antes


def test_cli_output_dir_formas_equivalentes():
    esperado = "/dados/uploads/cnpj/user/job/output"
    variantes = [
        ["caso.inp", "--output-dir", esperado],
        ["caso.inp", "--output-dir=" + esperado],
        ["-output", esperado, "caso.inp"],
        ["caso.inp", "-output=" + esperado],
        ["--output", esperado, "caso.inp", "-debug=0"],
    ]
    for argv in variantes:
        inp, output_dir, debug = parse_run_args(argv)
        assert inp == "caso.inp"
        assert output_dir == esperado
        assert debug == 0


def test_cli_rejeita_opcao_desconhecida_e_diretorio_vazio():
    with pytest.raises(ValueError):
        parse_run_args(["caso.inp", "--saida", "/tmp/job"])
    with pytest.raises(ValueError):
        parse_run_args(["caso.inp", "--output-dir"])
    with pytest.raises(ValueError):
        parse_run_args(["a.inp", "b.inp"])


def test_prepare_grava_so_no_diretorio_do_job(tmp_path):
    job = tmp_path / "uploads" / "12345678000199" / "user-7" / "job-9" / "output"
    compartilhada = tmp_path / "output"
    compartilhada.mkdir()
    marcador = compartilhada / "antigo.pvd"
    marcador.write_text("manter", encoding="utf-8")
    dentro = job / "ja_existia.txt"
    job.mkdir(parents=True)
    dentro.write_text("preservar", encoding="utf-8")

    base = step_output_base(str(job), "pasta/caso.inp", 1)
    assert base == os.path.join(os.path.abspath(job), "caso_1")
    assert dentro.read_text(encoding="utf-8") == "preservar"
    assert marcador.read_text(encoding="utf-8") == "manter"
    assert list(compartilhada.iterdir()) == [marcador]


def test_dois_jobs_nao_compartilham_diretorio(tmp_path):
    job_a = tmp_path / "job-a" / "output"
    job_b = tmp_path / "job-b" / "output"
    base_a = step_output_base(str(job_a), "viga.inp", 1)
    base_b = step_output_base(str(job_b), "viga.inp", 2)
    assert os.path.dirname(base_a) == os.path.abspath(job_a)
    assert os.path.dirname(base_b) == os.path.abspath(job_b)
    assert base_a != base_b
    aberto = prepare_output_dir(str(job_a))
    assert aberto == os.path.abspath(job_a)


def test_pvd_fica_com_caminhos_relativos_dentro_do_diretorio(tmp_path):
    out = tmp_path / "output"
    out.mkdir()
    pvtu = out / "caso_1000000.pvtu"
    vtu = out / "caso_1_p0_000000.vtu"
    pvd = out / "caso_1.pvd"
    vtu.write_text("<VTKFile></VTKFile>", encoding="utf-8")
    pvtu.write_text(
        '<VTKFile><PUnstructuredGrid>'
        '<Piece Source="/tmp/absolut/caso_1_p0_000000.vtu"/>'
        "</PUnstructuredGrid></VTKFile>",
        encoding="utf-8",
    )
    pvd.write_text(
        '<?xml version="1.0"?>\n'
        '<VTKFile type="Collection" version="1.0">\n'
        "  <Collection>\n"
        '    <DataSet timestep="0" part="0" '
        'file="/tmp/absolut/caso_1000000.pvtu"/>\n'
        "  </Collection>\n"
        "</VTKFile>\n",
        encoding="utf-8",
    )
    fora = tmp_path / "fora.vtu"
    fora.write_text("FORA", encoding="utf-8")
    outro = out / "outro_job.pvd"
    outro.write_text('file="/tmp/nao/mexer.pvtu"', encoding="utf-8")

    ensure_relative_vtk_links(str(pvd))

    texto = pvd.read_text(encoding="utf-8")
    assert 'file="caso_1000000.pvtu"' in texto
    assert "/tmp/absolut" not in texto
    pvtu_texto = pvtu.read_text(encoding="utf-8")
    assert 'Source="caso_1_p0_000000.vtu"' in pvtu_texto
    assert "/tmp/" not in pvtu_texto
    assert fora.read_text(encoding="utf-8") == "FORA"
    assert outro.read_text(encoding="utf-8") == 'file="/tmp/nao/mexer.pvtu"'


def test_pvd_nao_segue_referencia_fora_do_diretorio(tmp_path):
    out = tmp_path / "output"
    out.mkdir()
    secreto = tmp_path / "secreto"
    secreto.mkdir()
    alvo = secreto / "caso.pvtu"
    alvo.write_text('Source="/etc/passwd"', encoding="utf-8")
    pvd = out / "caso_1.pvd"
    pvd.write_text(
        '<VTKFile><Collection>'
        '<DataSet timestep="0" part="0" file="../secreto/caso.pvtu"/>'
        "</Collection></VTKFile>",
        encoding="utf-8",
    )

    ensure_relative_vtk_links(str(pvd))

    assert 'file="caso.pvtu"' in pvd.read_text(encoding="utf-8")
    assert alvo.read_text(encoding="utf-8") == 'Source="/etc/passwd"'
