# SPRINT 3 — FenicsEDSS × CalculiX 2.23 na roda A02 (candidato)

**Rotulo: CANDIDATO / nao congelado.** Nao e resultado definitivo. **Gate V nao aprovado. Gate M nao aprovado.** Nenhuma tolerancia de aceite foi inventada. As formulas C4/C5/C6 aparecem abaixo so como **PROPOSTA**.

Data: 13/11/2026, America/Sao_Paulo (UTC-3). Autor: Mike. Branch local Nitro `fc/demo-alwheel-static-c3d4-p0`. Sem push neste lote, sem Avell. Containers e programas do usuario (`femedss-app-app-1`, `mngsys-app-1`, `ubuntu-calculix`, `femedss-app-db-1`) nao foram parados nem alterados. Swap da maquina nao foi criado nem modificado.

## 1. INP candidato (somente leitura)

| | |
|---|---|
| Caminho | `/home/valdir/GROK_BOTS/Luiz (CAE specialist)/RODAS/R02-GENERICA-EDSS/analises/A02-estatica-pisadas30-refino/DEMO-RODAGEN-STATIC-02.inp` |
| sha256 (antes e depois) | `dd9466db39471c76f8baee741cc1bc798e80bb8aef0f7281604189992ea6b06c` |
| Tamanho | 161 527 377 B (161,53 MB) |
| Malha | 724 613 nos, 3 019 484 C3D4, E = 70 000 MPa, ν = 0,33 |
| Cargas | so `*CLOAD` no GDL 2; ΣFy = −5000 N exactos; FIX nos 5 furos |
| Resultados CCX de Luiz | mesma pasta `run/` (`.frd`, `.dat`, `post_results.json`, `VMe.npy`, `fixn.npy`, `ELN.npy`) — tambem somente leitura |

O INP **nao foi copiado, movido nem reescrito**. Foi montado somente leitura no container.

## 2. Perfil e solver

Pipeline FC habitual: `python3 src/run_from_inp.py --profile fc --units mm-N-MPa-t-s --solver cg-gamg --rtol 1e-10`, imagem `edss/fenicsx:dev`, `docker run --rm` ephemeral, 4 CPUs, teto `--memory` = `--memory-swap` = 13,520 GiB (MemAvailable 17,520 GiB menos 4 GiB reservados para o desktop). Sem fallback para LU, C3D8 ou C3D10.

O perfil FC da branch `a6e13c9` recusou o INP na primeira tentativa util por dois cartoes fora da whitelist: `*DENSITY` e `*EL PRINT, ELSET=Eall`. **O INP nao foi editado.** O codigo do perfil passou a:

- ignorar `*DENSITY` de forma **explicita** no manifesto (metadado; a massa nao entra na estatica linear, sem gravitacao nem inercia);
- aceitar `*EL PRINT` como pedido de saida de tensao, no mesmo espirito de `*EL FILE` / `*NODE PRINT`.

Teste `tests/test_fc_profile.py`: **29 passed** no container efemero da mesma imagem. Docs: `docs/FC_PROFILE.md`.

(Uma tentativa anterior morreu em 1 s porque a imagem nao tem `/usr/bin/time`; o pico passou a ser medido pelo cgroup e por `resource.ru_maxrss`. Logs em `_run/container.log.attempt1-notime` e `_run/attempt2-whitelist/`.)

## 3. Recursos e tempos (medidos)

Estimativa previa (escala linear do pico 1,485 GiB da roda auxiliar, 476 960 C3D4): 8,58 GiB (por DOF) a 9,40 GiB (por elementos), 11,75 GiB com +25 %. A corrida coube no teto de 13,52 GiB.

| | |
|---|---|
| Inicio / fim (America/Sao_Paulo) | 13/11/2026 12:02:31 – 12:06:29 |
| Parede do container | 3 min 58 s (238 s) |
| Solve + dump (processo) | 215,6 s (`time.txt`) |
| CG+GAMG (estagio `solve` − `apply_cloads`) | 99,2 s |
| Pos-processamento | 22 s (12:06:07–12:06:29) |
| Pico cgroup (`memory.peak`) | **9 188 057 088 B = 8,555 GiB** |
| HWM do processo (`FC_MEM` / `ru_maxrss`) | 7990,4 MiB = 7,803 GiB, no estagio `fc_fields_rf` |
| Iteracoes CG | **122** (1 reposicao de residual, apos 100 its com residual verdadeiro 1,399×10⁻¹⁰) |
| Residual verdadeiro \|\|Au−b\|\|/\|\|b\|\| | **7,609×10⁻¹¹** ≤ 1×10⁻¹⁰ (`true_residual_meets_rtol=true`) |
| Razao PETSc | `CONVERGED_RTOL` (codigo 2) |
| Near-nullspace | 6 modos de corpo rigido |

Estagios `FC_MEM` (RSS / HWM em MiB, tempo acumulado):

| Estagio | t (s) | RSS | HWM |
|---|---:|---:|---:|
| mesh | 17,7 | 1972 | 2733 |
| inp_node_to_geom | 31,0 | 2225 | 2733 |
| inp_element_to_cell | 42,9 | 2900 | 3251 |
| assemble_A | 67,6 | 5477 | 5477 |
| apply_cloads | 67,7 | 5477 | 5477 |
| solve (CG+GAMG) | 167,0 | 5692 | 6968 |
| fc_fields_rf | 205,2 | 6009 | 7990 |

## 4. Equilibrio (ΣRF)

Carga aplicada no Fenics (CLOAD do INP): Fy = −5000,000000000 N (Fx = Fz = 0).

| | FenicsEDSS | CCX 2.23 `.dat` (set FIX) | Δ (Fenics − CCX) |
|---|---:|---:|---:|
| RFx (N) | −2,11×10⁻⁹ | −3,76×10⁻⁹ | — |
| **RFy (N)** | **4999,999999993** | **5000,000** | **−6,78×10⁻⁹** |
| RFz (N) | +1,25×10⁻¹⁰ | −8,42×10⁻⁹ | — |

Σ(R+F) no Fenics: (−2,11×10⁻⁹, −6,77×10⁻⁹, +1,25×10⁻¹⁰) N. Momentos Σ(R+F) na origem: da ordem de 10⁻⁶ N·mm. Residual maximo nos DOFs livres: 2,78×10⁻¹⁰.

## 5. Deslocamento nodal (Fenics P1 × CCX `.frd` DISP)

Casamento por ID de no INP, 724 613 nos comuns (0 ausentes). O mapa dof↔geometria coincide (max \|x_dof − x_geom\| = 6,4×10⁻¹⁴ mm).

| | FenicsEDSS | CCX `.frd` |
|---|---|---|
| \|U\|max | **0,353107476 mm** | **0,353107603 mm** |
| No | **373552** (o mesmo) | **373552** |
| xyz (mm) | (0,92980, 260,09632, −128,34945) | identico |
| Componentes (mm) | (−5,28427×10⁻⁵, −0,26752276, −0,23047009) | (−5,28427×10⁻⁵, −0,267523, −0,23047) |

Δ(\|U\|max) = −1,27×10⁻⁷ mm (Fenics − CCX). O valor do CCX bate com o `post_results.json` do Luiz (0,353107603 mm, no 373552).

Diferenca nodal \|Δu\| (mm):

- max = **7,10×10⁻⁷ mm** (2,01×10⁻⁶ · \|U\|max_CCX), no 65803, (107,515, 191,879, −71,771) mm
- RMS = **1,87×10⁻⁷ mm** (5,28×10⁻⁷ · \|U\|max_CCX)
- p50 = 4,76×10⁻⁸ mm; p99 = 5,60×10⁻⁷ mm

A magnitude do maximo e a do arredondamento de impressao do `.frd` (formato ~7 algarismos).

**C4 PROPOSTA** (nao aceite, nao congelada): \|Δu_k\| ≤ 1×10⁻⁶·u_ref + 1×10⁻⁵·\|u_CCX,k\| por componente, com u_ref = \|U\|max do CCX = 0,353107603 mm. **0 violacoes** em 724 613 nos (max erro/tol = 0,371).

## 6. Tensao (base comparavel)

- **FenicsEDSS:** Cauchy **DG0**, constante no C3D4 P1 (1 valor por elemento), nao suavizado nos nos; tensor simetrizado antes dos principais.
- **CCX (base desta comparacao):** `*EL PRINT` no `.dat`, **1 ponto de integracao por C3D4**, ordem Voigt xx, yy, zz, xy, xz, yz. E a mesma base que o Luiz usou em `VMe.npy` (diferenca Luiz vs `.dat` relida = 0).
- **CCX `.frd`:** extrapola/media tensao para os nos; **nao** e a base principal desta comparacao.

3 019 484 elementos comparados. Destes, 2138 tocam o FIX; 3 017 346 ficam fora do FIX.

| | Fenics DG0 | CCX `.dat` IP | mesmo eid? |
|---|---:|---:|---|
| VM global max (MPa) | **33,689592** | **33,689591** | sim, **2821722** (toca o FIX) |
| VM fora do FIX (MPa) | **22,716863** | **22,716865** | sim, **2817879** |
| S1 global (min / max, MPa) | −9,191621 / 31,597008 | −9,191623 / 31,597012 | — |
| S3 global (min / max, MPa) | −39,768375 / 9,325591 | −39,768372 / 9,325593 | — |

Os maximos do CCX reproduzem o `post_results.json` do Luiz (VM 33,69 MPa no FIX e 22,72 MPa fora do FIX).

Diferencas brutas (Fenics − CCX), todos os elementos:

- max \|ΔS_ij\| = **5,00×10⁻⁶ MPa** (eid 1331115, fora do FIX; o desvio cai no xx, no ultimo algarismo impresso do `.dat`)
- RMS de ‖ΔS‖_F = 3,46×10⁻⁷ MPa
- max \|ΔVM\| = **8,06×10⁻⁶ MPa** (eid 505495, fora do FIX)
- RMS \|ΔVM\| = 2,80×10⁻⁷ MPa (2,79×10⁻⁷ fora do FIX; 7,26×10⁻⁷ nos que tocam o FIX)

**C5 PROPOSTA** (nao aceite): \|Δσ_ij\| ≤ 1×10⁻⁵·σ_ref + 1×10⁻⁴·\|σ_ij,CCX\|, com σ_ref = max \|S_ij\| do CCX IP = 39,43733 MPa. **0 violacoes** (max erro/tol = 0,0036).

**C6 PROPOSTA** (nao aceite; o CCX nao e julgado por C6): a mesma folga com max(\|S1\|, \|S3\|)_CCX, para principais e VM. **0 violacoes**.

## 7. Artefatos

Pasta de resultados (Nitro):

`/home/valdir/GROK_BOTS/Mike (Dev FemEDSS)/FC/A02-fenics/`

| Ficheiro | Conteudo |
|---|---|
| `_run/out/run_status.json`, `run_manifest.json`, `DEMO-RODAGEN-STATIC-02_1_reactions.json` | solver, residual verdadeiro, RF |
| `_run/out/compare.json`, `compare_brief.json` | comparacao completa e resumo |
| `_run/out/*.pvd` / `*.vtu` | VTK do pipeline FC (U nodal, S/E/principais/VM em celula DG0, RF) |
| `_run/out/xdmf/DEMO-RODAGEN-STATIC-02.xdmf` + `.h5` | XDMF do DOLFINx (malha + U) |
| `_run/out/xdmf/DEMO-RODAGEN-STATIC-02-fields.xdmf` + `*.bin` | XDMF binario para ParaView (U, VM, principais) |
| `_run/out/npy/` | campos e mapas INP↔celula (nao vao para o git) |
| `_run/out/time.txt`, `_run/stats.csv`, `_run/peak.snapshot`, `_run/container.log` | tempos e pico de memoria |
| `mem_estimate.json` | estimativa previa |

Os VTK/XDMF/npy (~3,4 GiB) **nao entram no git**.

## 8. Conclusao (sem gate)

O FenicsEDSS (CG+GAMG, residual verdadeiro 7,61×10⁻¹¹, 8,56 GiB de pico) reproduz o CCX 2.23 de Luiz neste INP candidato ate a precisao com que o CalculiX imprime o `.frd` e o `.dat`: equilibrio em 10⁻⁹ N, \|U\|max no mesmo no, max \|Δu\| = 7,1×10⁻⁷ mm, max \|ΔS\| = 5×10⁻⁶ MPa, maximos de VM no mesmo elemento (33,69 MPa no FIX e 22,72 MPa fora). As formulas C4/C5/C6 propostas fecham com folga, **mas isso nao e aceite oficial**. O INP continua **candidato, nao congelado**. Gate V e Gate M permanecem fechados ate o Valdir congelar os criterios e a referencia (u_ref, σ_ref, L_ref).
