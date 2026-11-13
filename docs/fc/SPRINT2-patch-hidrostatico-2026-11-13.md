# Sprint 2 — patch test (T2) e teste hidrostático (T0)

**Data da entrega dos testes:** 13/11/2026, America/Sao_Paulo (UTC−3).  
**Branch local:** `fc/demo-alwheel-static-c3d4-p0`.  
**Status:** **PROPOSTA.** Este documento **não congela** critérios, **não** declara os Gates V e M aprovados e **não** é o caso definitivo da roda.

Autorização do Valdir para rodar os testes (13/11 10:09): «podemos seguir».

**Decisão do Valdir (13/11/2026 às 10:22, America/Sao_Paulo), citação exata:** «vamos usar o CalculiX como referencia».

Consequência, até o congelamento (D8): nos testes T2 (patch) e T0 (hidrostático), o FenicsEDSS é comparado ao **CalculiX 2.23** pelas fórmulas **C4 e C5** da revisão Luiz 11/11 17:00, **não** contra a solução exata. O CCX é a referência e **não é julgado por C6** (nem pelo limite proposto σvm ≤ 1e-12·|p|).

Luiz (revisão 13/11 da sprint 2) propôs referências para estes testes, ainda **pendentes de decisão do Valdir**:

| Teste | u_ref (PROPOSTA) | σ_ref (PROPOSTA) |
|---|---|---|
| T2 patch | 0,222 mm | 289,473684 MPa |
| T0 hidrostático | 4,857e-4 mm | 10 MPa |

Estas propostas são as usadas nos asserts do pytest e nas colunas `C4_*` / `C5_*` das tabelas ponto a ponto. **Não estão congeladas.** Para COORD, Luiz propôs \|Δx\| ≤ 1e-9·L_ref + 2·ulp (ver §8 abaixo); isso permanece **proposta pendente** e **não** é imposta como critério congelado.

---

## 1. O que foi feito

Foram implementados e rodados, pelo leitor `.inp` do perfil FC (C3D4 linear elástico, E = 70 000 MPa, ν = 0,33, unidades mm-N-MPa-t-s):

- **T2 — patch test** em um bloco distorcido de 125 nós / 384 C3D4 (27 nós interiores), com campo linear `u = a + B x` cujo tensor de deformação tem cisalhamentos **não nulos e distintos** nas três componentes (εxy = 0,003, εxz = 0,0045, εyz = 0,0055). Três variantes: Dirichlet em todo o contorno; as 12 permutações pares da conectividade do tetraedro; Dirichlet parcial (faces x = 0, y = 0, z = 0) com tração consistente nas faces opostas. A variante de tração **é resolvida a partir do arquivo** `tests/fixtures/t2_patch_tracao.inp`.
- **T0 — hidrostático** no mesmo recorte (faces mantidas planas, volume = 1 000 mm³), pressão uniforme p = 10 MPa em todas as faces (forças nodais consistentes `F_i = A t / 3`; o perfil FC **rejeita** `*DLOAD`) e vínculos de simetria ux = 0 em x = 0, uy = 0 em y = 0, uz = 0 em z = 0. Também permanece o teste sintético do pós-processador (`test_t0_hidrostatico_von_mises`).

O motor **não** gera `.inp` de produção. As fixtures de verificação estão em `tests/fixtures/`. Não houve fallback automático de elemento (somente C3D4). Não foi tocado o INP original de Luiz, nem `docker-compose.yaml`, nem os containers persistentes da Nitro. Imagem efêmera do Fenics: `edss/fenicsx:dev` (image ID local `33180f8f7f1d`, **não** é RepoDigest). CCX 2.23 em container efêmero (`ubuntu-calculix:latest`, `/opt/CalculiX/bin/ccx`, SPOOLES); o container vivo `ubuntu-calculix` permaneceu Up e **não** foi tocado.

Arquivos desta entrega:

- `tests/fc_verify_meshes.py`
- `tests/test_fc_patch_hydrostatic.py`
- `tests/fixtures/t2_patch_block.inp` (SHA-256 `4e0cadb41273f7015e2ea5db621f014f1b55410943fa3e617c448853dbd01076`) e cópia compacta `.ccx.inp` (`cadc0d7e23d24717670f7891039cc49a533c691408b490de43b8057ea4d3ae98`)
- `tests/fixtures/t2_patch_tracao.inp` (SHA-256 `ee8f6f63cd7a3da9ae9cd693486853fbaf5a59923f27a4c90c36dc4307be7b32`) e cópia compacta `.ccx.inp` (`a44c7dfb2acf04270d20376f3dddd26941c0e34de7cb38b8d1102854edb2c888`)
- `tests/fixtures/t0_hydrostatic_block.inp` (SHA-256 `6d81d1be36a8de1c9f39f30b24d2e974e77439b26e4a91b5c1537b0a1fe4d763`) e cópia compacta `.ccx.inp` (`1c59044a97e0e475dd7b5f753fae529c987d0facb23f5539819e3b525514aff2`)
- tabelas ponto a ponto em `docs/fc/sprint2-fenics-vs-ccx/`
- este documento

As cópias `.ccx.inp` acrescentam `*HEADING` e `*EL PRINT`, aceitos no CCX e **rejeitados** pela whitelist do perfil FC. Os `*CLOAD` das fixtures FC foram realinhados ao gerador (soma determinística + formato `.15e`); as cópias compactas `.ccx.inp` permaneceram com o mesmo hash. Os arquivos já sujos no tree (`input/beam_nset.inp`, `input/teste.inp`, `docs/Information_TBD.odt`) **não** foram incluídos.

---

## 2. Comparação FenicsEDSS × CCX 2.23 (referência)

Fórmulas propostas (revisão Luiz 11/11 17:00), **não congeladas**:

- **C4:** |Δu_k| ≤ 1e-6·u_ref + 1e-5·|u_CCX,k| (por componente, casamento por ID de nó)
- **C5:** |Δσ_ij| ≤ 1e-5·σ_ref + 1e-4·|σ_ij,CCX| (casamento **por nome**: σxx, σyy, σzz, σxy, σxz, σyz; cisalhamento tensorial)
- **C6:** |Δσ_k| ≤ 1e-5·σ_ref + 1e-4·máx(|σ1|,|σ3|)_CCX, e a mesma regra para VM; ordenação S1 ≥ S2 ≥ S3 com sinal. O CCX, sendo a referência, **não** é julgado por C6.
- **COORD:** a fórmula escrita na revisão 17:00 é |Δx| ≤ 1e-9·L_ref. Luiz (13/11) propôs acrescer um piso de formato: |Δx| ≤ 1e-9·L_ref + 2·ulp₁₃.₆(\|x\|), ou de forma equivalente ≤ 1e-9·L_ref + 1e-6·(1 mm + \|x\|), com L_ref = 10 mm (aresta do cubo) ou o L_ref geral que o Valdir congelar. **Proposta pendente: não é aplicada como critério congelado.**

Solver Fenics desta comparação: **LU**, like-for-like com o SPOOLES do CCX. 125 nós e 384 elementos casados em todos os casos; nenhum ID ausente.

**Tabelas ponto a ponto** (um registro por nó e por elemento), referenciadas aqui e gravadas em:

- `docs/fc/sprint2-fenics-vs-ccx/t2_patch_dirichlet_nos.csv` e `t2_patch_dirichlet_elementos.csv`
- `docs/fc/sprint2-fenics-vs-ccx/t2_patch_tracao_nos.csv` e `t2_patch_tracao_elementos.csv`
- `docs/fc/sprint2-fenics-vs-ccx/t0_hydrostatic_nos.csv` e `t0_hydrostatic_elementos.csv`
- resumo: `docs/fc/sprint2-fenics-vs-ccx/resumo.json` (e `README.md` das colunas)

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
| C4/C5 com refs **PROPOSTAS** de Luiz | 0 violações; máx (erro/tol) C4 = 0,038; C5 = 0,002 (**não** é aceite oficial) |

A variante de **Dirichlet parcial + tração** deu os mesmos máximos até a precisão reportada (max\|Δu\| = 4,967e-8 mm; max\|Δσ\| = 3,684e-5 MPa; max\|ΔVM\| = 6,195e-6 MPa; max\|ΔCOORD\| = 4,996e-7 mm; C4/C5 propostos: 0 violações).

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
| C4/C5 com refs **PROPOSTAS** de Luiz | 0 violações; máx (erro/tol) C4 = 0,030; C5 ≪ 1 (**não** é aceite oficial) |

O CCX **não** é avaliado contra σvm ≤ 1e-12·|p|. Esse limite proposto permanece, quando couber, como checagem do **pós-processador** e da solução **em memória** do FenicsEDSS, não como juízo sobre a referência.

---

## 3. Consistência do campo no espaço P1 (FenicsEDSS)

Os asserts do pytest passaram a usar as fórmulas C4/C5/C6 com as referências **PROPOSTAS** por Luiz (13/11), listadas no cabeçalho. Não são mais folgados que essas fórmulas. Continuam **não** sendo aceite oficial nem critério congelado. No T0, σvm em memória é assertado contra 1e-12·|p|.

### T2 — patch test (LU, versus o campo linear exato)

| Variante | max\|du\| (mm) | max\|du\| interiores / livres (mm) | max\|dS\| (MPa) | Σ(R+F) máx (N) | momento máx (N·mm) | resíduo livre máx (N) |
|---|---:|---:|---:|---:|---:|---:|
| Dirichlet no contorno | 1,388e-16 | 1,388e-16 | 8,768e-12 | 1,381e-10 | 1,745e-9 | 8,819e-11 |
| 12 permutações pares (pior) | 2,220e-16 | — | 1,111e-11 | — | — | — |
| Dirichlet parcial + tração (arquivo) | 2,776e-16 | 2,776e-16 | 1,252e-11 | 1,892e-10 | 1,208e-9 | 1,039e-10 |

O campo linear é recuperado nos nós interiores no nível de arredondamento. Volume mínimo 0,913 mm³, todos os 384 tets com volume > 0, inclusive nas 12 permutações pares.

### T0 — hidrostático (versus σ = −p I, em memória)

Esperado analiticamente: σ = −p I, S1 = S2 = S3 = −10 MPa, σvm = 0, α = −p(1−2ν)/E = −4,85714e-5, u = α x, εvol = 3α = −1,45714e-4.

| Solver | max\|du\| (mm) | max\|dS\| (MPa) | max σvm (MPa) | max\|Δ(S1,S2,S3)\| (MPa) | max\|Δεvol\| | Σ(R+F) máx (N) | momento máx (N·mm) | residual verdadeiro \|Au−b\|/\|b\| |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| LU | 2,168e-18 | 4,441e-14 | 4,868e-14 | 4,796e-14 | 4,337e-19 | 8,455e-13 | 8,672e-12 | 2,884e-15 |
| CG+GAMG (16 its) | 1,209e-16 | 5,704e-12 | 3,803e-12 | 6,056e-12 | 6,535e-17 | 1,597e-10 | 1,020e-9 | **2,050e-13** ≤ 1e-10 |

O pós-processador sintético e o solver FenicsEDSS (LU e CG+GAMG), **em memória**, atendem σvm ≤ 1e-12·|p| = 1e-11 MPa neste bloco. Isso **não** se usa para julgar o CCX.

Volume 1 000 mm³, mínimo 0,960 mm³, 100 % dos tets com volume > 0. Campos finitos, S1 ≥ S2 ≥ S3 com sinal.

---

## 4. Mapeamento C0–C11 (PROPOSTAS da revisão Luiz 11/11 17:00)

Nenhuma tolerância nova foi inventada além de registrar as referências **propostas** por Luiz em 13/11, ainda pendentes do Valdir.

| ID | Papel proposto | O que estes testes mostram | Lacuna / decisão |
|---|---|---|---|
| **C0** | PASS (carga Pisada) | **Não se aplica.** | Pendente no caso da roda. |
| **C1** | PASS (bench) | T0 CG+GAMG: residual verdadeiro 2,050e-13 ≤ **1e-10**. LU: 2,884e-15. | Declarar o rtol no manifesto da run oficial. |
| **C2** | PASS | T0 LU: Σ(R+F) = 8,46e-13 N; CG: 1,60e-10 N. Ambos ≪ 5e-3 N (F_ref da roda). | F_ref é da roda. Repetir no definitivo. |
| **C3** | PASS | T0 LU: momento máx 8,67e-12 N·mm; CG: 1,02e-9 N·mm. Ambos ≪ 1,23 N·mm (L_ref da roda). | Congelar **c**; repetir no definitivo. |
| **C3b** | só reportar | Não feito. | Fora desta sprint. |
| **C4** | PASS em B1 e **T2** | Fenics × CCX no T2 (125 nós, tabela ponto a ponto): max\|Δu\| = **4,967e-8 mm**. Com u_ref **PROPOSTO** 0,222 mm: 0 violações (máx erro/tol = 0,038). **Não congelado.** | **Pendente:** o Valdir aceitar ou corrigir u_ref. |
| **C5** | PASS em B1 e **T2** | Fenics × CCX no T2, por nome e por ID: max\|Δσ\| = **3,684e-5 MPa** (σxy). Com σ_ref **PROPOSTO** 289,473684 MPa: 0 violações (máx erro/tol = 0,002). **Não congelado.** | **Pendente:** congelar σ_ref. |
| **C6** | PASS em **T0** e B1 | Fenics × CCX no T0: max\|Δ(S1,S2,S3)\| = 1,440e-11 MPa; max\|ΔVM\| = 2,183e-11 MPa. O **CCX não é julgado** por C6 nem por σvm ≤ 1e-12·\|p\|. O Fenics, em memória e no pós-processador, atende 1e-12·\|p\|. | **Pendente:** congelar σ_ref e confirmar o papel do C6. |
| **C7** | PASS (B1) | **Não se aplica.** | Congelar os números do B1. |
| **C8** | PASS (B1) | **Não se aplica.** | Congelar raio/σ_ref do B1. |
| **C9** | só reportar (roda) | **Não se aplica.** | Família homóloga da geometria congelada. |
| **C10** | PASS (volume > 0) | 100 % dos 384 tets com volume > 0 no T0 e no T2 (e nas 12 permutações). | Remalha do definitivo (D9) inalterada. |
| **C11** | PASS | U nodal, Cauchy / S1–S3 / VM por célula DG0, finitos, 384 células. COORD medido (max\|Δx\| = 4,996e-7 mm); a fórmula 1e-9·L_ref + 2·ulp permanece proposta pendente, não é imposta. | Repetir no `.inp` definitivo. |

---

## 5. Revisão Luiz 13/11 — achados e resolução

Fonte (somente leitura): `/home/valdir/GROK_BOTS/Luiz (CAE specialist)/FC-revisoes/SPRINT2-revisao-Luiz-2026-11-13.md`, sobre o commit `2514891`. Veredito dele: aprovado com ressalvas. Esta seção registra o que foi feito depois, na mesma branch local.

| # | Achado | Resolução |
|---|---|---|
| 1 (média) | Comparação Fenics × CCX like-for-like, por ID, não estava tabelada. | Incluídas as tabelas ponto a ponto em `docs/fc/sprint2-fenics-vs-ccx/` (nós e elementos dos três casos), referenciadas na §2. Números coincidem com o resumo já medido no commit `7c0bd3b`. |
| 2 (média) | O teste de COORD (\|Δx\| ≤ 1e-9·L_ref) fica abaixo do piso de impressão do `.dat`. | Documentada a proposta de Luiz: \|Δx\| ≤ 1e-9·L_ref + 2·ulp₁₃.₆(\|x\|), ou ≤ 1e-9·L_ref + 1e-6·(1 mm + \|x\|). **Não** é imposta como critério congelado; o max\|ΔCOORD\| medido (4,996e-7 mm) fica só reportado. |
| 3 (baixa) | Suíte Fenics não reproduzida no computador do Luiz. | Sem ação de código. A suíte é reexecutada nesta correção (ver §7). |
| 4 (baixa) | Asserts do pytest mais folgados que C4/C5/C6; C6 do solver só impresso. | Asserts passaram a usar C4/C5/C6 com as refs **PROPOSTAS** de Luiz (T2: 0,222 mm e 289,473684 MPa; T0: 4,857e-4 mm e 10 MPa). No T0, σvm em memória é assertado contra 1e-12·\|p\|. Tudo rotulado PROPOSTA, pendente do Valdir. |
| 5 (baixa) | `test_t2_e_t0_via_fixture_arquivo` não carregava a fixture de tração; a variante de tração rodava só do gerador em memória. | A variante de tração e o teste via arquivo **leem** `tests/fixtures/t2_patch_tracao.inp`. Há também assert de que as três fixtures batem com o gerador atual. |
| 6 (baixa) | Sete linhas de `*CLOAD` em `t2_patch_tracao.inp` diferiam do gerador por 1 ULP. | Gerador passou a somar com `math.fsum`, em ordem determinística, e a gravar `*CLOAD` com `.15e`. Fixtures de tração e hidrostático realinhadas (hashes atualizados na §1). As cópias `.ccx.inp` compactas não mudaram de hash. |
| 7 (baixa) | Aviso de incremento mínimo nulo no `*STATIC` dos decks CCX. | Inócuo neste passo único; não alterado. |

---

## 6. Decisões pedidas ao Valdir (para congelar)

A referência **CalculiX** para T0 e T2 já foi decidida (13/11 10:22). Enquanto o restante não for escrito, os números acima permanecem **PROPOSTAS** e os Gates V e M **não** estão aprovados.

1. **Congelar C0–C11** (item D8), com as fórmulas da revisão 17:00 (κ já corrigido) e os papéis PASS / só reportar confirmados pelo Luiz em 12/11.
2. **Declarar u_ref, σ_ref e L_ref** para estes testes. Luiz propôs os valores da tabela do cabeçalho, e para COORD o piso 1e-9·L_ref + 2·ulp. Este documento usa essas refs só como **PROPOSTA** nos testes e nas tabelas; **não** as dá como congeladas.
3. Nada disto substitui o `.inp` congelado da roda genérica nem a comparação Fenics × CCX nesse caso (checkpoint 20/11).

---

## 7. Como reproduzir

Na Nitro, sem alterar containers persistentes:

```bash
docker run --rm --name fenics-fc-sprint2-fullsuite \
  --cpus=2 --memory=4g --memory-swap=4g \
  -v /mnt/MAIN/000_VALDIR/GIT_PJR/FemEDSS/Fenics-EDSS:/app \
  -w /app edss/fenicsx:dev \
  python3 -m pytest /app/tests -v --tb=short
```

A comparação Fenics × CCX desta nota usou as fixtures de `tests/fixtures/`, o CCX em container efêmero (`ubuntu-calculix:latest`, sem tocar o container vivo; decks compactos de hash inalterado) e o Fenics em `edss/fenicsx:dev` com `--solver lu`.

---

## 8. O que isto não é

Não é Gate V aprovado, não é Gate M aprovado, não é a roda da demo, não é aceite de tensão e não é identidade demonstrada com o CalculiX além das diferenças medidas na seção 2. Os critérios e as referências u_ref / σ_ref / piso de COORD continuam **PROPOSTAS** até o Valdir congelar.
