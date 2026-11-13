# Sprint 2 — patch test (T2) e teste hidrostático (T0)

**Data da entrega dos testes:** 13/11/2026, America/Sao_Paulo (UTC−3).  
**Branch local:** `fc/demo-alwheel-static-c3d4-p0`.  
**Status:** **PROPOSTA.** Este documento **não congela** critérios, **não** declara os Gates V e M aprovados e **não** é o caso definitivo da roda.

Autorização do Valdir para rodar os testes (13/11 10:09): «podemos seguir».

**Decisão do Valdir (13/11/2026 às 10:22, America/Sao_Paulo), citação exata:** «vamos usar o CalculiX como referencia».

Consequência, até o congelamento (D8): nos testes T2 (patch) e T0 (hidrostático), o FenicsEDSS é comparado ao **CalculiX 2.23** pelas fórmulas **C4 e C5** da revisão Luiz 11/11 17:00, **não** contra a solução exata. O CCX é a referência e **não é julgado por C6** (nem pelo limite proposto σvm ≤ 1e-12·|p|). `u_ref`, `σ_ref` e `L_ref` **não foram escolhidos**; permanecem pendentes, sem número inventado neste documento, e as violações oficiais de C4/C5/C6/COORD **não são contabilizadas**.

---

## 1. O que foi feito

Foram implementados e rodados, pelo leitor `.inp` do perfil FC (C3D4 linear elástico, E = 70 000 MPa, ν = 0,33, unidades mm-N-MPa-t-s):

- **T2 — patch test** em um bloco distorcido de 125 nós / 384 C3D4 (27 nós interiores), com campo linear `u = a + B x` cujo tensor de deformação tem cisalhamentos **não nulos e distintos** nas três componentes (εxy = 0,003, εxz = 0,0045, εyz = 0,0055). Três variantes: Dirichlet em todo o contorno; as 12 permutações pares da conectividade do tetraedro; Dirichlet parcial (faces x = 0, y = 0, z = 0) com tração consistente nas faces opostas.
- **T0 — hidrostático** no mesmo recorte (faces mantidas planas, volume = 1 000 mm³), pressão uniforme p = 10 MPa em todas as faces (forças nodais consistentes `F_i = A t / 3`; o perfil FC **rejeita** `*DLOAD`) e vínculos de simetria ux = 0 em x = 0, uy = 0 em y = 0, uz = 0 em z = 0. Também permanece o teste sintético do pós-processador (`test_t0_hidrostatico_von_mises`).

O motor **não** gera `.inp` de produção. As fixtures de verificação estão em `tests/fixtures/`. Não houve fallback automático de elemento (somente C3D4). Não foi tocado o INP original de Luiz, nem `docker-compose.yaml`, nem os containers persistentes da Nitro (`femedss-app-app-1`, `mngsys-app-1`, `ubuntu-calculix`, `femedss-app-db-1`). Imagem efêmera do Fenics: `edss/fenicsx:dev` (image ID local `33180f8f7f1d`, **não** é RepoDigest). CCX 2.23 em container efêmero `fenics-fc-sprint2-ccx` (`ubuntu-calculix:latest`, `/opt/CalculiX/bin/ccx`, SPOOLES); o container vivo `ubuntu-calculix` permaneceu Up e **não** foi tocado.

Arquivos desta entrega (testes no commit `2514891`; este texto atualizado no commit seguinte da mesma branch):

- `tests/fc_verify_meshes.py`
- `tests/test_fc_patch_hydrostatic.py`
- `tests/fixtures/t2_patch_block.inp` (SHA-256 `4e0cadb4…d01076`) e cópia compacta `.ccx.inp`
- `tests/fixtures/t2_patch_tracao.inp` (SHA-256 `1391c7b2…461e11`) e cópia compacta `.ccx.inp`
- `tests/fixtures/t0_hydrostatic_block.inp` (SHA-256 `a198735c…22f617`) e cópia compacta `.ccx.inp`
- este documento

As cópias `.ccx.inp` acrescentam `*HEADING` e `*EL PRINT`, aceitos no CCX e **rejeitados** pela whitelist do perfil FC. Os arquivos já sujos no tree (`input/beam_nset.inp`, `input/teste.inp`, `docs/Information_TBD.odt`) **não** foram incluídos.

---

## 2. Comparação FenicsEDSS × CCX 2.23 (referência)

Fórmulas propostas (revisão Luiz 11/11 17:00), **não congeladas**:

- **C4:** |Δu_k| ≤ 1e-6·u_ref + 1e-5·|u_CCX,k| (por componente, casamento por ID de nó)
- **C5:** |Δσ_ij| ≤ 1e-5·σ_ref + 1e-4·|σ_ij,CCX| (casamento **por nome**: σxx, σyy, σzz, σxy, σxz, σyz; cisalhamento tensorial)
- **C6:** |Δσ_k| ≤ 1e-5·σ_ref + 1e-4·máx(|σ1|,|σ3|)_CCX, e a mesma regra para VM; ordenação S1 ≥ S2 ≥ S3 com sinal. O CCX, sendo a referência, **não** é julgado por C6.
- **COORD:** |Δx| ≤ 1e-9·L_ref (PI do CCX versus baricentro da célula no Fenics)

Solver Fenics desta comparação: **LU**, like-for-like com o SPOOLES do CCX. 125 nós e 384 elementos casados em todos os casos; nenhum ID ausente.

### Medido — T2 patch (Dirichlet no contorno)

| Grandeza | Valor |
|---|---|
| max\|Δu\| | **4,967e-8 mm** (ux 4,967e-8; uy 4,959e-9; uz 4,719e-8) |
| max\|u_CCX\| | 0,222 mm |
| max\|Δσ\| | **3,684e-5 MPa** (σxy) |
| max\|Δσ\| por nome (MPa) | σxx 4,025e-7; σyy 9,288e-7; σzz 1,486e-5; σxy 3,684e-5; σxz 5,263e-6; σyz 1,579e-5 |
| max\|Δ(S1,S2,S3)\| | 2,732e-5 MPa (em S2) |
| max\|ΔVM\| | 6,195e-6 MPa |
| max\|ΔCOORD\| | 4,996e-7 mm |

A variante de **Dirichlet parcial + tração** deu os mesmos máximos até a precisão reportada (max\|Δu\| = 4,967e-8 mm; max\|Δσ\| = 3,684e-5 MPa; max\|ΔVM\| = 6,195e-6 MPa; max\|ΔCOORD\| = 4,996e-7 mm).

Estas diferenças estão na casa da impressão do `.dat` do CCX (~7 algarismos). Não se afirma identidade dos solvers além do que foi impresso.

### Medido — T0 hidrostático (p = 10 MPa)

| Grandeza | Valor |
|---|---|
| max\|Δu\| | **4,994e-11 mm** (ux 4,994e-11; uy 4,981e-11; uz 4,929e-11) |
| max\|u_CCX\| | 4,857e-4 mm |
| max\|Δσ\| | **1,041e-11 MPa** (σxy; diagonais ~4e-14) |
| max\|Δ(S1,S2,S3)\| | 1,440e-11 MPa |
| max\|ΔVM\| | 2,183e-11 MPa |
| max\|ΔCOORD\| | 4,994e-7 mm |

O CCX **não** é avaliado contra σvm ≤ 1e-12·|p|. Esse limite proposto permanece, quando couber, como checagem do **pós-processador** do FenicsEDSS, não como juízo sobre a referência.

`u_ref`, `σ_ref` e `L_ref` estão **pendentes**. Por isso este documento **não conta violações** de C4, C5, C6 ou COORD. Só para ordem de grandeza (não é aceite): o termo rtol já escrito nas fórmulas, 1e-5·|u_CCX|, vale no máximo 2,220e-6 mm no T2 e 4,857e-9 mm no T0; o termo 1e-4·|σ_ij,CCX| no T2 chega a 2,895e-2 MPa em σyz.

---

## 3. Consistência do campo no espaço P1 (FenicsEDSS, não é o aceite)

Os asserts do pytest usam abs 1e-8 mm e 1e-6 MPa, no mesmo espírito do teste axial já existente. São consistência numérica (LU / campo no espaço P1), **não** aceite de Luiz, **não** a comparação de referência e **não** critério congelado. A comparação de aceite desta sprint, por decisão de 13/11 10:22, é a da seção 2 (Fenics × CCX).

### T2 — patch test (LU, versus o campo linear exato)

| Variante | max\|du\| (mm) | max\|du\| interiores / livres (mm) | max\|dS\| (MPa) | Σ(R+F) máx (N) | momento máx (N·mm) | resíduo livre máx (N) |
|---|---:|---:|---:|---:|---:|---:|
| Dirichlet no contorno | 1,388e-16 | 1,388e-16 | 8,768e-12 | 1,381e-10 | 1,745e-9 | 8,819e-11 |
| 12 permutações pares (pior) | 2,220e-16 | — | 1,111e-11 | — | — | — |
| Dirichlet parcial + tração | 2,776e-16 | 2,776e-16 | 1,252e-11 | 1,892e-10 | 1,208e-9 | 1,039e-10 |

O campo linear é recuperado nos nós interiores no nível de arredondamento. A tensão é constante, componente a componente. Volume mínimo 0,913 mm³, todos os 384 tets com volume > 0, inclusive nas 12 permutações pares.

### T0 — hidrostático (versus σ = −p I)

Esperado analiticamente: σ = −p I, S1 = S2 = S3 = −10 MPa, σvm = 0, εvol = −3p(1−2ν)/E = −1,45714e-4, u = α x com α = −4,85714e-5.

| Solver | max\|du\| (mm) | max\|dS\| (MPa) | max σvm (MPa) | max\|Δ(S1,S2,S3)\| (MPa) | max\|Δεvol\| | Σ(R+F) máx (N) | momento máx (N·mm) | residual verdadeiro \|Au−b\|/\|b\| |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| LU | 2,168e-18 | 4,441e-14 | 4,868e-14 | 4,796e-14 | 4,337e-19 | 8,455e-13 | 8,672e-12 | 2,884e-15 |
| CG+GAMG (16 its) | 1,209e-16 | 5,704e-12 | 3,803e-12 | 6,056e-12 | 6,535e-17 | 1,597e-10 | 1,020e-9 | **2,050e-13** ≤ 1e-10 |

O pós-processador sintético atende ao limite proposto σvm ≤ 1e-12·|p| = 1e-11 MPa. O solver FenicsEDSS também produziu σvm abaixo desse valor neste bloco; isso **não** se usa para julgar o CCX.

Volume 1 000 mm³, mínimo 0,960 mm³, 100 % dos tets com volume > 0. Campos finitos, S1 ≥ S2 ≥ S3 com sinal.

### Suíte completa

`python3 -m pytest /app/tests -v` na imagem `edss/fenicsx:dev` (container efêmero `fenics-fc-sprint2-fullsuite`): **69 passed** em 6,20 s (America/Sao_Paulo, 13/11/2026 ~10:17). Inclui os 8 testes novos e o T0 sintético já existente.

---

## 4. Mapeamento C0–C11 (PROPOSTAS da revisão Luiz 11/11 17:00)

Nenhuma tolerância nova foi inventada. Onde o critério não se aplica, isso é dito explicitamente. A referência de T2 e T0, por decisão de 13/11 10:22, é o CCX 2.23.

| ID | Papel proposto | O que estes testes mostram | Lacuna / decisão |
|---|---|---|---|
| **C0** | PASS (carga Pisada) | **Não se aplica.** Não há conjuntos Pisada1/Pisada2. | Continua pendente no caso da roda, depois de congelar as pisadas. |
| **C1** | PASS (bench) | T0 com CG+GAMG: residual verdadeiro 2,050e-13 ≤ **1e-10** (critério já proposto / D16). LU: 2,884e-15. | Declarar o rtol no manifesto da run oficial da roda. |
| **C2** | PASS | T0 LU: Σ(R+F) = 8,46e-13 N; CG: 1,60e-10 N. Ambos ≪ 5e-3 N (1e-6·F_ref com F_ref = 5 000 N da roda). T2 também fecha (~1e-10 N). | F_ref é da roda, não deste cubo. Repetir no caso definitivo. |
| **C3** | PASS | T0 LU: momento máx 8,67e-12 N·mm; CG: 1,02e-9 N·mm. Ambos ≪ 1,23 N·mm (1e-6·F_ref·L_ref com L_ref = 245,486 mm da roda). | L_ref é da roda (o cubo tem L = 10 mm). Congelar **c**; repetir no definitivo. |
| **C3b** | só reportar | Não feito (é comparação de momento com o CCX da roda). | Fora desta sprint. |
| **C4** | PASS em B1 e **T2** | Fenics × CCX no T2 (125 nós): max\|Δu\| = **4,967e-8 mm**. A fórmula proposta precisa de **u_ref**, que **não foi escolhido**. Violações **não contabilizadas**. | **Pendente:** o Valdir congelar u_ref (e os demais números da revisão 17:00). |
| **C5** | PASS em B1 e **T2** | Fenics × CCX no T2, casamento por nome e por ID de elemento: max\|Δσ\| = **3,684e-5 MPa** (σxy). σ_ref **não foi escolhido** (o 75 MPa do B1 não é deste estado). Violações **não contabilizadas**. | **Pendente:** congelar σ_ref. |
| **C6** | PASS em **T0** e B1 | Fenics × CCX no T0: max\|Δ(S1,S2,S3)\| = 1,440e-11 MPa; max\|ΔVM\| = 2,183e-11 MPa; ordenação 100 % no Fenics. O **CCX não é julgado** por C6 nem por σvm ≤ 1e-12·\|p\|. O pós-processador sintético do Fenics atende a 1e-12·\|p\|. σ_ref pendente; violações da forma mista **não contabilizadas**. | **Pendente:** congelar σ_ref e confirmar que C6 compara o Fenics à referência CCX, sem julgar o CCX. |
| **C7** | PASS (B1) | **Não se aplica** (é convergência em flexão do B1, já evidenciada). | Congelar os números do B1. |
| **C8** | PASS (B1) | **Não se aplica.** | Congelar raio/σ_ref do B1. |
| **C9** | só reportar (roda) | **Não se aplica.** | Família homóloga da geometria congelada. |
| **C10** | PASS (volume > 0) | 100 % dos 384 tets com volume > 0 no T0 e no T2 (e nas 12 permutações). | O q da roda continua só reportado; remalha do definitivo (D9) inalterada. |
| **C11** | PASS | U nodal, Cauchy / S1–S3 / VM por célula DG0, finitos, ordenados com sinal, 384 células, nos dois testes. COORD do CCX casou com o baricentro (max\|Δx\| = 4,996e-7 mm); L_ref pendente, violações de COORD não contabilizadas. | Repetir no `.inp` definitivo da roda. |

---

## 5. Decisões pedidas ao Valdir (para congelar)

A referência **CalculiX** para T0 e T2 já foi decidida (13/11 10:22). Enquanto o restante não for escrito, os números acima permanecem **PROPOSTAS** e os Gates V e M **não** estão aprovados.

1. **Congelar C0–C11** (item D8), com as fórmulas da revisão 17:00 (κ já corrigido) e os papéis PASS / só reportar confirmados pelo Luiz em 12/11.
2. **Declarar u_ref, σ_ref e L_ref** aplicáveis a estes testes (ou declarar que se usam só as referências já propostas da roda / do B1). Este documento **não escolhe** esses valores e **não conta** violações.
3. Nada disto substitui o `.inp` congelado da roda genérica nem a comparação Fenics × CCX nesse caso (checkpoint 20/11).

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

A comparação Fenics × CCX desta nota usou as fixtures de `tests/fixtures/`, o CCX em container efêmero (`ubuntu-calculix:latest`, sem tocar o container vivo) e o Fenics em `edss/fenicsx:dev` com `--solver lu`.

---

## 7. O que isto não é

Não é Gate V aprovado, não é Gate M aprovado, não é a roda da demo, não é aceite de tensão e não é identidade demonstrada com o CalculiX além das diferenças medidas na seção 2. Os critérios continuam **PROPOSTAS** até o Valdir congelar.
