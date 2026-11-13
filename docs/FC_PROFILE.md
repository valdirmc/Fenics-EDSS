# Perfil FC (Fenics-EDSS) — o que o codigo faz agora

Documento do repositorio, alinhado ao codigo na branch `fc/demo-alwheel-static-c3d4-p0`.
Nao substitui os canonicos do Alex. Nao e validacao de roda.

## Opt-in

```
python3 /app/src/run_from_inp.py arquivo.inp --output-dir DIR --profile fc --units mm-N-MPa-t-s
```

Sem `--profile fc` o comportamento anterior permanece (C3D8, termico, hiperelastico, cards desconhecidas ignoradas).

`--units` e obrigatorio no perfil FC. Alias aceito: `mm-N-s-t` -> `mm-N-MPa-t-s` (N/mm^2 = MPa, massa em t).

## Whitelist e rejeicoes

Aceita: *NODE, *ELEMENT (TYPE=C3D4 obrigatorio), *NSET, *ELSET (incluindo alias), *MATERIAL, *ELASTIC, *SOLID SECTION, *BOUNDARY (forma `nset,dof` ou `nset,dof1,dof2,valor`), *CLOAD, um *STEP, *STATIC (um incremento, periodo 1, tempo inicial 0, sem heranca), *END STEP, *NODE FILE, *EL FILE, *NODE PRINT, *EL PRINT, *OUTPUT.

Ignorados **explicitamente** (manifesto, nunca silencio): *PHYSICAL CONSTANTS, *EXPANSION, *DENSITY, parametro `INC` de *STEP.

Rejeita (exit 2): keyword fora da lista, *ELEMENT sem TYPE, C3D8/C3D10, multistep, *DLOAD, *NLGEOM, *HYPERELASTIC, *HEAT TRANSFER, *CONDUCTIVITY, *PLASTIC, *FREQUENCY e demais bloqueados, incrementacao *STATIC alem de um passo de periodo 1.

## Campos

Sempre obrigatorios no perfil FC: U, S (Cauchy), principais algebricos sigma1>=sigma2>=sigma3 **com sinal**, von Mises.
Se o .inp pede E ou RF, tambem tem de ser honrados — parser sozinho nao e sucesso (exit 2 se faltar; exit 3 se falhar gravacao VTK).

Localizacao: U nodal (Lagrange P1); S, E, principais e von Mises **por celula DG0**, constantes no C3D4 P1. Suavizacao nodal nao substitui. RF nos nos presos e totais (`*_reactions.json`).

Ordenacao do tensor no manifesto: row-major XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ.

## CLI / artefatos

- `--output-dir` e o contrato de saida.
- Falha VTK nunca e sucesso (nao ha mais `except` que engole e devolve 0).
- `run_manifest.json` lista artefatos desta run com SHA-256. Lixo velho no diretorio nao conta como sucesso.
- Exit: 1 arquivo/uso; 2 CLI ou rejeicao do perfil FC; 3 falha de solver/VTK/artefato.

## Imagem

A run registra **versoes importadas** (dolfinx etc.). Image ID local (`docker image inspect .Id`) identifica o artefato nesta maquina e **nao** e RepoDigest de registry, nem referencia de pull. Nao usar `repo:tag@sha256:<imageId>` como pin. Tag `edss/fenicsx:dev` e `dolfinx/dolfinx:stable` continuam flutuantes.

## Material 71700 vs 70000

O plano falava E=70000 MPa, nu=0.33. O INP de referencia de Luiz (`FEMMeshGmsh.inp`, 11/11/2026) exporta E=71700, nu=0.33, comentario Aluminum-7075-T6. O codigo **nao corrige** nem presume equivalencia. A decisao e de Alex/Valdir.

## O que este lote nao faz

- Nao gera/exporta INP de roda; a entrada e .inp fornecido.
- Nao implementa envelope de pisada ±30° nem eixos/areas (dados de geometria ainda insuficientes; so CLOAD/BOUNDARY genericos).
- Nao roda o modelo de 125395 nos / 476960 C3D4 (sem estimativa de recurso nesta passagem).
- Nao afirma validacao de roda, benchmark independente, nem identidade com CalculiX.
- Go/no-go alvo do usuario: 13/11. Nao promete roda aceita.

## Lote 2 (validacao/robustez, 11/11/2026)

Opt-in inalterado. Correções e contratos adicionais:

- Indexacao: `create_mesh` pode reordenar vertices. CLOAD/RF/material nao usam `node_idx * 3` nem o `coords` de entrada (no lote 1 o CLOAD referencia `coords` fora de escopo — regressao corrigida). Mapas: no INP → vertice `geometry.x` (coordenadas unicas) e elemento INP → celula (conjunto de vertices).
- CLOAD acumula no DOF localizado pelas coordenadas do no (mesmo criterio das BC, agora via indice vetorizado; um `collapse` por componente). Duas linhas no mesmo no/DOF somam.
- `*BOUNDARY`: CalculiX e `no, dof_ini, dof_fim[, valor]`. Tres colunas com valor flutuante sao ambiguas e rejeitadas no FC; valor nao nulo exige 4 colunas.
- Matriz com BC so no solve. RF = `K u - f_ext` (`assemble_vector` esparso da forma, apenas nos DOFs prescritos; nao ha K denso nem `A_full` unused). JSON traz `sum_R_plus_F` e momentos na origem. Desequilibrio acima da tolerancia numerica de consistencia recusa sucesso.
- MPI: perfil FC so em serial (`size=1`). `size>1` e rejeitado (ghosts nao validados).
- Tensor de Cauchy simetrizado `(S+S.T)/2` antes dos principais; ordem algebrica `σ1≥σ2≥σ3` com sinal; von Mises a partir desses principais.
- `*EXPANSION` so pode ser ignorado no perfil estritamente mecanico, sem temperatura; `*TEMPERATURE` / `*INITIAL CONDITIONS` / fluxos termicos sao bloqueados.
- `INC` em `*STEP` permanece ignorado de forma explicita no estatico linear de um incremento.
- `*OUTPUT` so e honrado com `FREQUENCY=1` (ou sem parametro). Outro parametro e rejeitado: reconhecer o cartao nao basta.
- CLI: `--parse-only` grava `parse_report.json` com `STATUS: parse_ok_not_solved` e **nao** executa o solver nem gera manifesto de sucesso. O parser CLI imprime o mesmo aviso.
- `run_status.json` (`started` / `ok` / `failed` / `parse_ok_not_solved`) com `run_id`. Manifesto antigo no diretorio nao conta como sucesso da run atual. Falha de parse/solver/VTK nao reescreve manifesto `ok`.
- Sistema sem `*BOUNDARY` e recusado (LU que nao acusa singularidade nao e sucesso). Solucao com NaN/Inf ou KSP nao convergido e recusada.

### Fixture analitica (nao e a roda)

Prisma triangular de 6 nos / 3 C3D4, L=1 mm, A=0,5 mm², E=70000 MPa, ν=0,33, Fz=−50 N (tracao constante no topo, cargas consistentes F_i=A t/3; Zienkiewicz & Taylor, vetor de carga do triangulo linear). Estado uniaxial no espaco P1:

- σ_zz=F/A=−100 MPa, demais nulas; principais (0, 0, −100); von Mises=100 MPa
- ε_zz=σ_zz/E, ε_lat=−ν σ_zz/E
- ΣR+ΣF=0 e momentos na origem ~0, com DOFs livres no topo

Tolerancia dos testes: abs 1e-6 MPa/N ou 1e-8 mm (LU, campo no espaco FE). **Nao** e aceite de Luiz nem GateV aprovado. Ha tambem caso de deslocamento prescrito nao nulo (sem CLOAD) e casos de erro (sem BC, parse-only, lixo velho).

Flexao/convergencia formal ficou para o proximo lote (o axial+compressao fechou neste).

### Caso auxiliar E=71700 (decisao Valdir, 11/11/2026)

Autorizado usar o INP original de Luiz (E=71700 MPa, cargas simplificadas) **apenas como teste auxiliar**, preservando o arquivo e o hash, com saidas rotuladas AUXILIAR/E71700/cargas simplificadas. O caso definitivo permanece E=70000 e o perfil combinado. Este lote **nao rodou** o modelo de 125395 nos / 476960 C3D4: o solver FC ainda usa LU direto, e nao ha medicao que mostre caber com margem no teto isolado de 2 GB sem risco a outros processos da Nitro. Comparacoes CCX so com campos/localizacao compativeis e tolerancias fornecidas; nenhuma tolerancia de aceite foi inventada.

## Lote 3 (11/11/2026) — B1 e recursos

Benchmark B1 (proposta de Luiz): prisma 200×20×20 mm, flexão pura Saint-Venant, E=70000, ν=0,33, M=1,0e5 N·mm, σ_máx=75 MPa. Malha estruturada, 6 tets/hex (Freudenthal, diagonal n000–n111), 2/4/8/16 elementos na espessura com nx=10 n (elementos aproximadamente cúbicos). Apoios 3-2-1 em x=0, sem engaste. Cargas consistentes T3 nas faces x=0 e x=L.

Nota: Luiz citou κ=1,07e-6 /mm; as fórmulas κ=M/(E I) com σ_máx=75 MPa dão **κ=1,0714e-4 /mm**. O código usa as fórmulas.

C3D4 é rígido em flexão (esperado). Flecha e rotação menores que a exata em todos os níveis; energia ½f·u_h sobe monotonicamente em direção à exata (0,49 → 0,78 → 0,93 → 0,98). Tensões parasitárias (σyy, σzz, cisalhamentos) presentes e decaem com o refinamento.

Critérios C7/C8/C2/C3 são **PROPOSTAS**, não congelados. **Gate V não está aprovado.** Contra as propostas, na malha de 16 na espessura: |erro uy| e |erro rotação| ≤ 2 %; σxx médio na esfera (baricentros) ≤ 5 %; ordens L2(u) em 4→8 e 8→16 dentro de [1,6; 2,4]; ordens L2(σ) dentro de [0,7; 1,3]. A esfera usa tets com baricentro no raio de 2 mm (o volume do conjunto muda com h; não é recorte geométrico exato).

Recursos (LU direto, container isolado, processo fresco): ver nota de execução do lote 3. A extrapolação incremental (N^1,35) para ~376k DOFs cai na casa de 10 GiB sem margem, e 20–30 GiB com fator 2–3. Com ~11 GiB disponíveis na Nitro e serviços no ar, o caso auxiliar da roda **não foi executado**.

### Caso auxiliar E=71700 — aborto 11/11 18:56 America/Sao_Paulo (opcao B)

Valdir autorizou (18:56 America/Sao_Paulo, retransmitido por Alex) parar temporariamente só `mngsys-app-1` e `ubuntu-calculix` para liberar RAM, religando depois, e rodar o INP original (E=71700, cargas simplificadas Face119/121) como teste **AUXILIAR**, com critério MemAvailable ≥ ~20 GiB após a parada (margem 2× sobre a previsão LU de ~10,4 GiB).

Inventário e checagem de atividade (18:57–18:59 America/Sao_Paulo): sem requests mngsys nos últimos 15 min; Guto-MngSys.log parado em 10/11 18:46; nenhum `ccx`; Postgres do mngsys é o servidor 18 do host em :5432 (não o `femedss-app-db-1` em :5434) e não foi parado. Os dois containers autorizados somavam **108 MiB** (94 + 14). MemAvailable = **10,40 GiB**; após uma parada a estimativa seria ≈ **10,50 GiB**, ainda ~9,5 GiB abaixo do critério.

**Decisão:** abortar **sem** parar os containers e **sem** executar o solver (evitar interrupção inútil). Hash do INP conferido somente leitura: `9380e87b6073d3180b85591bcd677cd5be7d15f52a870fcc2fba6a2c3fce7e92`. Sem RF, equilíbrio, U, tensões nem comparação CCX. Gate V não aprovado. Nota: `nota-execucao-FC-auxiliar-2026-11-11.md`.

## Correção de memoria (11/11/2026, apos OOM auxiliar)

A rodada auxiliar E=71700 (125395 nos, 476960 C3D4, 8367 BC, 8582 CLOAD) foi OOMKilled a 18,7 GiB ~55 min depois de "Materiais aplicados", com crescimento linear de ~1,22 GiB / 6 s. B1 n=16 (140 k DOFs, poucos BC) tinha pico de só 2,8 GiB: o estouro nao era o LU.

Causa: `collect_constrained_dofs` em `src/fc_assembly.py` chamava `V.sub(c).collapse()`, `locate_dofs_geometrical` e `fem.Function(Vcol)` **por cada linha *BOUNDARY**. No modelo da roda isso retém ~2 MiB por BC (Function + mapa do collapse) e chega a dezenas de GiB. `apply_cloads` tinha o mesmo locate no laço (tempo O(N·nCLOAD), memoria liberada). Alem disso, `inp_node_to_geom` era O(N²) em alocações e `get_material_for_element` fazia `elem_id in list` por elemento (O(N²) no tempo, ~48 min antes de "Materiais aplicados").

Correção (resultados numericos inalterados; conferido no B1 2/4/8/16 e nos testes axiais):

- Um `ComponentDofLocator` (um collapse por componente + indice de coordenadas O(N)); BC/CLOAD em O(n) lookups.
- Um `dirichletbc` + `Function` por (componente, valor), nunca por no.
- `inp_node_to_geom` com hash de coordenadas (mesmas tolerâncias 1e-9 / 1e-6).
- Cache linear de material por elemento.
- Removido o `A_full` não usado. RF continua `assemble_vector` esparso só nos DOFs prescritos.
- `FC_MEM_PROFILE=1` imprime RSS/HWM por estágio.

A roda **nao foi reexecutada** neste lote. O vazamento de BC nao e mais o limitante; o pico esperado continua sendo o LU (~10 GiB na extrapolação incremental N^1,35 a partir do B1 n=16, ou ~20 GiB com fator 2). Ver nota de execucao `nota-execucao-FC-fix-memoria-2026-11-11.md`.

## Solver iterativo opt-in (11/11/2026) — CG+GAMG

CLI: `--solver lu|cg-gamg` (default `lu`) e `--solver cg-gamg --rtol 1e-10`. O `rtol` maior que 1e-10 e rejeitado. Nao ha fallback para LU: se o KSP nao convergir, a run termina com erro (exit ≠ 0) e `run_status.json` em `failed`.

Implementacao (`src/fc_solver.py`): PETSc KSP CG + PC GAMG (aggregation), residual verdadeiro (UNPRECONDITIONED), `max_it=10000`, near-nullspace elastico com 6 modos de corpo rigido (3 translações + 3 rotações em torno da origem, ortonormalizados) via `A.setNearNullSpace` (nao `setNullSpace`, que seria incorreto apos Dirichlet). Iteracoes, `||r||`, `||b||`, `||r||/||b||` e `converged_reason` vao para `run_status.json` / `run_manifest.json`.

Validacao rapida B1 (mesmas malhas do lote 3, n=8 e n=16; diferencas brutas, sem tolerancia de aceite inventada), evidencia `FC/cg-vs-lu-b1/cg_vs_lu_summary.json`:

- n=8 (19683 DOFs): max|dU|=1,019e-9 (rel. 5,09e-10), max|dS|=4,97e-7 (rel. 6,40e-9); CG 27 its, `||r||/||b||`=8,53e-11, 1,12 s vs LU 3,19 s.
- n=16 (139587 DOFs): max|dU|=6,44e-9 (rel. 3,06e-9), max|dS|=6,78e-6 (rel. 8,74e-8); CG 30 its, `||r||/||b||`=5,32e-11, 9,33 s vs LU 58,78 s.

Testes do perfil FC no container `edss/fenicsx:dev`: 60 passaram. Testes de host (parser + output_dir): 12 passaram.

### Execucao auxiliar E=71700 com CG+GAMG (11/11 21:18 America/Sao_Paulo)

Autorizacao de Valdir 11/11 21:10 America/Sao_Paulo, retransmitida por Alex: tentar o solver iterativo, maquina liberada exceto o grok-bot. Job detached `fenics-fc-auxiliar-e71700-cggamg` (cid `3653ba2082d2…`), imagem local `edss/fenicsx:dev`, INP copia (original intacto, SHA-256 `9380e87b6073d3180b85591bcd677cd5be7d15f52a870fcc2fba6a2c3fce7e92` antes e depois). Teto `--memory`=`--memory-swap`=19,592 GiB (MemAvailable 23,592 GiB − 4 GiB), 4 CPUs. Sem push, sem Avell, sem contato com outros containers.

O solve CG+GAMG convergiu (200 its, `CONVERGED_RTOL`, `||r||/||b||` verdadeiro = 2,326e-10; `solve_exit=0`, pico cgroup 1,479 GiB, ~41 s ate o manifesto). O pos-processamento *dentro* do container falhou (exit 1): a imagem nao tem `meshio`/`vtk`. O pos-processamento foi refeito no host sobre os VTU/JSON ja gravados (leitor ASCII VTU 2.2; glob do stem do PVD foi corrigido para nao pegar `_RF`/`_cell` como U). **Rotulo AUXILIAR / nao definitivo / nao validado. Gate V nao aprovado.** Tolerancias de Luiz (C2/C3/C4/energia) sao propostas, nao aceite. Detalhe e comparacao CCX (so U nodal e totais do `.dat`): `nota-execucao-FC-auxiliar-iterativo-2026-11-11.md`.

## Noite 11/11/2026 — residual verdadeiro, esfera geometrica, calibracao CCX B1

Criterio oficial do CG+GAMG: **TRUE** `||Au−b||/||b|| ≤ rtol` (autorizado 1e-10), verificado com `A.mult` apos o solve. O rnorm UNPRECONDITIONED recursivo do CG pode ficar abaixo do rtol enquanto o residual verdadeiro ainda nao (roda auxiliar, primeira run: KSP ~8,76e-11 vs verdadeiro 2,326e-10). Nao ha `setConvergenceTest` (uma tentativa anterior nao parava). O KSP usa rtol interno `min(rtol*0,1, 1e-12)`; se o residual verdadeiro ainda exceder o rtol, o solver reinicia a partir do `u` atual (residual replacement, r0 verdadeiro, ciclos de ate 100 its, no maximo 40). Se ao final `rel > rtol` ou o KSP divergir, a run termina com `RuntimeError` / exit ≠ 0, **sem fallback para LU**.

A imagem compartilhada `edss/fenicsx:dev` **nao foi alterada** (meshio/VTU nao entraram nela). O pos-processamento da roda usa o script no host (leitor ASCII VTU 2.2), bind-mount em `/scripts/postprocess.py`.

### Roda auxiliar CG+GAMG com criterio verdadeiro (21:42:47–21:44:37 America/Sao_Paulo)

Evidencia da run anterior (200 its, residual verdadeiro 2,326e-10 > 1e-10) preservada em `_run_iterativo/`. Run nova em `_run_iterativo_true_rtol/` (tentativas 1–2 arquivadas no mesmo diretorio). Job `fenics-fc-auxiliar-e71700-cggamg-tr`, teto 19,543 GiB, 4 CPUs, INP copia com SHA-256 `9380e87b6073d3180b85591bcd677cd5be7d15f52a870fcc2fba6a2c3fce7e92` antes e depois. **AUXILIAR, nao definitivo, nao validado. Gate V nao aprovado.**

MEASURED: `solve_exit=0`, 668 its, 11 residual replacements, `true ||r||/||b|| = 9,578e-11 ≤ 1e-10`, `true_residual_meets_rtol=true`, pico cgroup 1,485 GiB (1594101760 B) sob o teto de 19,543 GiB. max|U|=0,39003889 mm no no 3233. RFy=4999,927425905 N; Σ(R+F) ~ 1e-9 N. C4 proposta vs .dat Luiz: max‖du‖=7,02e-7 mm, 0 violacoes; energia rel 4,32e-8. Extremos S1/S2/S3/VM globais **nao** caem nos 183 elementos com q<0,1 (os valores fora dos pobres coincidem com o global). Detalhe na nota `nota-execucao-FC-noite-2026-11-11.md`.

### Esfera B1 (recorte geometrico)

O monitor da esfera (r=2 mm em (100, 7, 0)) passou a reportar volume tet ∩ esfera (subdivisao mid-edge, profundidade 5) alem do criterio por baricentro. Volumes geometricos aproximam 4/3 π R³ ≈ 33,5103 mm³. Erros de σ_xx medio (CG+GAMG, residual verdadeiro ≤ 1e-10), MEASURED:

- n=4: geo rel_err=−0,18913 (vol=33,4924, 12 tets); bary=−0,14701 (vol=41,6667)
- n=8: geo rel_err=−0,05907 (vol=33,5000, 38 tets); bary=−0,08874 (vol=41,6667)
- n=16: geo rel_err=−0,01620 (vol=33,5075, 248 tets); bary=−0,01628 (vol=33,8542)

Em n=4/8 o conjunto por baricentro tinha volume 41,67 mm³ (nao e uma bola geometrica). Propostas C7/C8, nao gate.

### Calibracao CCX B1 (C4–C6)

INP B1 n=8 e n=16 rodados no CCX 2.23 em container **temporario** `ubuntu-calculix:latest` (o container vivo `ubuntu-calculix` nao foi parado nem alterado). Decks com notacao cientifica longa (16 digitos) falharam (`*ERROR reading *NODE`, exit 201); a copia CCX foi compactada (E maiusculo / decimais curtos, linhas ≪ 72). `*EL PRINT` no ponto de integracao (C3D4: 1 IP) vs Cauchy DG0.

MEASURED, propostas, nao aceite. Tolerancia numerica C5/C6 de Luiz **nao encontrada** nos docs locais/Drive disponiveis; nao foi inventada.

- n=8: C4 proposta 0 violacoes, max‖du‖=5,049e-7 mm (u_ref=2,011 mm); max|dS|=5,234e-6 MPa (s_ref=77,636 MPa), 6561 nos / 30720 elems
- n=16: C4 proposta 0 violacoes, max‖du‖=5,364e-7 mm (u_ref=2,116 mm); max|dS|=4,379e-5 MPa (s_ref=77,542 MPa), 46529 nos / 245760 elems
