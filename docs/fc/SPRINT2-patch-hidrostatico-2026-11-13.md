# Sprint 2 — patch test (T2) e teste hidrostático (T0)

**Data:** 13/11/2026, America/Sao_Paulo (UTC−3).  
**Branch local:** `fc/demo-alwheel-static-c3d4-p0` (a partir de `a5344d2`).  
**Status:** **PROPOSTA.** Este documento **não congela** critérios, **não** declara os Gates V e M aprovados e **não** é o caso definitivo da roda.

Autorização do Valdir (13/11 10:09): «podemos seguir» com os testes da sprint 2, na branch local da Nitro, sem push e sem Avell.

---

## 1. O que foi feito

Foram implementados e rodados, pelo leitor `.inp` do perfil FC (C3D4 linear elástico, E = 70 000 MPa, ν = 0,33, unidades mm-N-MPa-t-s):

- **T2 — patch test** em um bloco distorcido de 125 nós / 384 C3D4 (27 nós interiores), com campo linear `u = a + B x` cujo tensor de deformação tem cisalhamentos **não nulos e distintos** nas três componentes (εxy = 0,003, εxz = 0,0045, εyz = 0,0055). Três variantes: Dirichlet em todo o contorno; as 12 permutações pares da conectividade do tetraedro; Dirichlet parcial (faces x = 0, y = 0, z = 0) com tração consistente nas faces opostas.
- **T0 — hidrostático** no mesmo recorte (faces mantidas planas, volume = 1 000 mm³), pressão uniforme p = 10 MPa em todas as faces (forças nodais consistentes `F_i = A t / 3`; o perfil FC **rejeita** `*DLOAD`) e vínculos de simetria ux = 0 em x = 0, uy = 0 em y = 0, uz = 0 em z = 0. Também foi refeito o teste sintético do pós-processador (já existente em `test_t0_hidrostatico_von_mises`).

O motor **não** gera `.inp` de produção. As fixtures de verificação estão em `tests/fixtures/`. Não houve fallback automático de elemento (somente C3D4). Não foi tocado o INP original de Luiz, nem `docker-compose.yaml`, nem os containers persistentes da Nitro (`femedss-app-app-1`, `mngsys-app-1`, `ubuntu-calculix`, `femedss-app-db-1`). Imagem efêmera: `edss/fenicsx:dev` (image ID local `33180f8f7f1d`, **não** é RepoDigest).

Arquivos desta entrega:

- `tests/fc_verify_meshes.py`
- `tests/test_fc_patch_hydrostatic.py`
- `tests/fixtures/t2_patch_block.inp` (SHA-256 `4e0cadb4…d01076`)
- `tests/fixtures/t2_patch_tracao.inp` (SHA-256 `1391c7b2…461e11`)
- `tests/fixtures/t0_hydrostatic_block.inp` (SHA-256 `a198735c…22f617`)
- cópias compactas `.ccx.inp` (com `*HEADING` e `*EL PRINT`, aceitos no CCX e **rejeitados** pela whitelist do perfil FC)
- este documento

Os arquivos já sujos no tree (`input/beam_nset.inp`, `input/teste.inp`, `docs/Information_TBD.odt`) **não** foram incluídos.

---

## 2. Resultados MEASURED (FenicsEDSS)

Tolerâncias usadas **somente** nos asserts do pytest: abs 1e-8 mm e 1e-6 MPa, no mesmo espírito do teste axial já existente. São consistência numérica (LU / campo no espaço P1), **não** aceite de Luiz e **não** critério congelado.

### T2 — patch test (LU)

| Variante | max\|du\| (mm) | max\|du\| interiores / livres (mm) | max\|dS\| (MPa) | Σ(R+F) máx (N) | momento máx (N·mm) | resíduo livre máx (N) |
|---|---:|---:|---:|---:|---:|---:|
| Dirichlet no contorno | 1,388e-16 | 1,388e-16 | 8,768e-12 | 1,381e-10 | 1,745e-9 | 8,819e-11 |
| 12 permutações pares (pior) | 2,220e-16 | — | 1,111e-11 | — | — | — |
| Dirichlet parcial + tração | 2,776e-16 | 2,776e-16 | 1,252e-11 | 1,892e-10 | 1,208e-9 | 1,039e-10 |

O campo linear é recuperado nos nós interiores (e nos livres da variante de tração) no nível de arredondamento. A tensão é constante e igual à analítica, componente a componente (σxx, σyy, σzz, σxy, σxz, σyz), o que detectaria uma troca de cisalhamentos. Volume mínimo 0,913 mm³, todos os 384 tets com volume > 0, inclusive nas 12 permutações pares.

### T0 — hidrostático (p = 10 MPa)

Esperado: σ = −p I, S1 = S2 = S3 = −10 MPa, σvm = 0, εvol = −3p(1−2ν)/E = −1,45714e-4, u = α x com α = −4,85714e-5.

| Solver | max\|du\| (mm) | max\|dS\| (MPa) | max σvm (MPa) | max\|Δ(S1,S2,S3)\| (MPa) | max\|Δεvol\| | Σ(R+F) máx (N) | momento máx (N·mm) | residual verdadeiro \|Au−b\|/\|b\| |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| LU | 2,168e-18 | 4,441e-14 | 4,868e-14 | 4,796e-14 | 4,337e-19 | 8,455e-13 | 8,672e-12 | 2,884e-15 |
| CG+GAMG (16 its) | 1,209e-16 | 5,704e-12 | 3,803e-12 | 6,056e-12 | 6,535e-17 | 1,597e-10 | 1,020e-9 | **2,050e-13** ≤ 1e-10 |

Critério **proposto** de Luiz para T0 no pós-processador: σvm ≤ 1e-12·\|p\| = 1e-11 MPa neste caso. O pós-processador sintético atende. O solver FenicsEDSS (LU e CG+GAMG) **também atendeu** neste bloco.

Volume 1 000 mm³, mínimo 0,960 mm³, 100 % dos tets com volume > 0. Campos finitos, S1 ≥ S2 ≥ S3 com sinal.

### Suíte completa

`python3 -m pytest /app/tests -v` na imagem `edss/fenicsx:dev` (container efêmero `fenics-fc-sprint2-fullsuite`): **69 passed** em 6,20 s (America/Sao_Paulo, 13/11/2026 ~10:17). Inclui os 8 testes novos e o T0 sintético já existente.

---

## 3. Checagem cruzada opcional no CalculiX 2.23

Mesmos decks, em cópia compacta `.ccx.inp`, no container **efêmero** `fenics-fc-sprint2-ccx` (`ubuntu-calculix:latest`, `/opt/CalculiX/bin/ccx`, SPOOLES). O container vivo `ubuntu-calculix` permaneceu Up e **não** foi tocado. Três jobs com `ccx_exit=0` (~0,01–0,03 s).

| Caso | max\|du\| vs analítico (mm) | max\|dS\| (MPa) | max σvm ou \|ΔVM\| (MPa) |
|---|---:|---:|---:|
| T2 Dirichlet | 4,967e-8 | 3,684e-5 | 6,195e-6 (\|ΔVM\|) |
| T2 tração | 4,967e-8 | 3,684e-5 | 6,195e-6 (\|ΔVM\|) |
| T0 hidrostático | 4,994e-11 | 1,040e-11 | **2,184e-11** (σvm) |

As diferenças de U/S no T2 estão na casa da impressão do `.dat` do CCX (~7 algarismos). No T0, o σvm impresso pelo CCX (2,184e-11 MPa) fica **ligeiramente acima** de 1e-12·\|p\| = 1e-11 MPa. Isso **não** é uma tolerância nova: é o valor medido, e entra na decisão abaixo. Não se afirma identidade dos solvers além do que foi impresso.

---

## 4. Mapeamento C0–C11 (PROPOSTAS da revisão Luiz 11/11 17:00)

Nenhuma tolerância nova foi inventada. Onde o critério não se aplica a estes testes, isso é dito explicitamente.

| ID | Papel proposto | O que estes testes mostram | Lacuna / decisão |
|---|---|---|---|
| **C0** | PASS (carga Pisada) | **Não se aplica.** Não há conjuntos Pisada1/Pisada2. | Continua pendente no caso da roda, depois de congelar as pisadas. |
| **C1** | PASS (bench) | T0 com CG+GAMG: residual verdadeiro 2,050e-13 ≤ **1e-10** (critério já proposto / D16). LU: 2,884e-15. | Declarar o rtol no manifesto da run oficial da roda. |
| **C2** | PASS | T0 LU: Σ(R+F) = 8,46e-13 N; CG: 1,60e-10 N. Ambos ≪ 5e-3 N (1e-6·F_ref com F_ref = 5 000 N da roda). T2 também fecha (~1e-10 N). | F_ref é da roda, não deste cubo. Repetir no caso definitivo. |
| **C3** | PASS | T0 LU: momento máx 8,67e-12 N·mm; CG: 1,02e-9 N·mm. Ambos ≪ 1,23 N·mm (1e-6·F_ref·L_ref com L_ref = 245,486 mm da roda). | L_ref é da roda (o cubo tem L = 10 mm). Congelar **c**; repetir no definitivo. |
| **C3b** | só reportar | Não feito (é comparação de momento com o CCX da roda). | Fora desta sprint. |
| **C4** | PASS em B1 e **T2** | Fenics vs campo analítico no T2: max\|du\| ~ 1e-16 mm (arredondamento). CCX vs analítico: 4,97e-8 mm (impressão do `.dat`). A fórmula proposta é \|Δu_k\| ≤ 1e-6·u_ref + 1e-5·\|u_CCX,k\| e foi escrita **contra o CCX**, não contra o campo exato. | **Decisão do Valdir:** aplicar essa forma no T2 contra o analítico? Qual u_ref? (Se u_ref = max\|u_ex\| = 0,222 mm, **não congelado**, o CCX teria 0 violações — só ilustração, não é escolha deste documento.) |
| **C5** | PASS em B1 e **T2** | Fenics vs analítico no T2: max\|dS\| = 8,77e-12 a 1,25e-11 MPa, por **nome** (σxx, σyy, σzz, σxy, σxz, σyz). CCX vs analítico: 3,68e-5 MPa. A fórmula proposta é \|Δσ_ij\| ≤ 1e-5·σ_ref + 1e-4·\|σ_ij,CCX\|. σ_ref do B1 (75 MPa) **não** é deste estado e **não está congelado**. | **Decisão do Valdir:** não foi escolhido σ_ref para o patch. Só se reportam os valores medidos. |
| **C6** | PASS em **T0** e B1 | Pós-processador: S1=S2=S3, σvm ≤ 1e-12·\|p\|. Solver Fenics (T0): max σvm LU 4,87e-14 e CG 3,80e-12, ambos ≤ 1e-11. Ordenação 100 %. CCX impresso: σvm = 2,18e-11, **acima** de 1e-11. | **Decisão do Valdir:** o 1e-12·\|p\| vale só para o pós-processador, ou também para o solver e para o valor impresso pelo CCX? |
| **C7** | PASS (B1) | **Não se aplica** (é convergência em flexão do B1, já evidenciada). | Congelar os números do B1. |
| **C8** | PASS (B1) | **Não se aplica.** | Congelar raio/σ_ref do B1. |
| **C9** | só reportar (roda) | **Não se aplica.** | Família homóloga da geometria congelada. |
| **C10** | PASS (volume > 0) | 100 % dos 384 tets com volume > 0 no T0 e no T2 (e nas 12 permutações). | O q da roda continua só reportado; remalha do definitivo (D9) inalterada. |
| **C11** | PASS | U nodal, Cauchy / S1–S3 / VM por célula DG0, finitos, ordenados com sinal, 384 células, nos dois testes. | Repetir no `.inp` definitivo da roda. |

---

## 5. Decisões pedidas ao Valdir (para congelar)

Enquanto estas decisões não forem escritas, os números acima permanecem **PROPOSTAS** e os Gates V e M **não** estão aprovados.

1. **Congelar C0–C11** (item D8), com as fórmulas da revisão 17:00 (κ já corrigido) e os papéis PASS / só reportar confirmados pelo Luiz em 12/11.
2. **T0 / C6:** o limite proposto σvm ≤ 1e-12·\|p\| é atendido pelo pós-processador e pelo solver FenicsEDSS neste bloco. O mesmo `.inp` no CCX 2.23 imprime σvm = 2,184e-11 MPa, um pouco acima de 1e-11 MPa (compatível com a precisão do `.dat`). Vale para o pós-processador isolado, para o solver, para os dois, ou também para a impressão do CCX?
3. **T2 / C4–C5:** o C0–C11 **não traz** uma tolerância de arredondamento contra o campo linear exato. Os erros do FenicsEDSS estão em ~1e-16 mm e ~1e-11 MPa. Não foi escolhido u_ref nem σ_ref para este estado. Se o Valdir quiser uma forma mista igual à do B1 vs CCX, precisa declarar as referências.
4. Nada disto substitui o `.inp` congelado da roda genérica nem a comparação Fenics × CCX nesse caso (checkpoint 20/11).

---

## 6. Como reproduzir

Na Nitro, sem alterar containers persistentes:

```bash
docker run --rm --name fenics-fc-sprint2-fullsuite \
  --cpus=2 --memory=4g --memory-swap=4g \
  -v /mnt/MAIN/000_VALDIR/GIT_PJR/FemEDSS/Fenics-EDSS:/app \
  -w /app edss/fenicsx:dev \
  python3 -m pytest /app/tests -v --tb=short
```

---

## 7. O que isto não é

Não é Gate V aprovado, não é Gate M aprovado, não é a roda da demo, não é aceite de tensão e não é identidade demonstrada com o CalculiX além das diferenças medidas acima.
