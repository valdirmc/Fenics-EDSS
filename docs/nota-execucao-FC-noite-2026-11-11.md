# Nota de execucao FC — noite 11/11/2026

America/Sao_Paulo (UTC−3). Branch local `fc/demo-alwheel-static-c3d4-p0`. Sem push, sem Avell. Containers do usuario (`femedss-app-app-1`, `mngsys-app-1`, `ubuntu-calculix`, `femedss-app-db-1`) nao foram parados nem alterados. INP original de Luiz intacto (SHA-256 `9380e87b6073d3180b85591bcd677cd5be7d15f52a870fcc2fba6a2c3fce7e92` antes e depois). Caso AUXILIAR/E71700/cargas-simplificadas: **nao definitivo, nao validado. Gate V nao aprovado.** Criterios de Luiz (C2–C8, C4–C6) sao propostas, nao aceite.

Autorizacao desta noite: Valdir dormindo; Alex concordou em continuar nas regras ja dadas, sem nova decisao de Valdir (pendencias listadas no fim).

---

## MEASURED

### 0) Criterio de parada com residual verdadeiro + rerun da roda

O rnorm UNPRECONDITIONED recursivo do CG pode cumprir o rtol do KSP enquanto `||Au−b||/||b||` (matvec independente) ainda nao. Implementacao em `src/fc_solver.py`:

- rtol interno do KSP = `min(rtol*0,1, 1e-12)`
- apos o solve, verifica TRUE `||Au−b||/||b||` com `A.mult`
- se ainda `rel > rtol`, residual replacement: reinicia o CG a partir do `u` atual (r0 verdadeiro), ciclos de no maximo 100 its, ate 40 ciclos
- se no fim `rel > rtol` ou o KSP divergir: `RuntimeError`, exit ≠ 0, **sem fallback para LU**
- um `setConvergenceTest` com residual mal avaliado no callback foi tentado e abandonado (nao parava)

A imagem compartilhada `edss/fenicsx:dev` **nao foi modificada** (incluir meshio/VTU exigiria rebuild da tag compartilhada). O pos-processamento permanece o script do host (leitor ASCII VTU 2.2; meshio 5.3.5 nao le VTKFile 2.2), bind-mount em `/scripts/postprocess.py` no job. Neste rerun, `postprocess_exit=0`.

Evidencia da run anterior (21:18:19–21:19:14 America/Sao_Paulo, 200 its, residual verdadeiro 2,326e-10 > 1e-10, pico 1,479 GiB) **preservada** em `FC/AUXILIAR/E71700/cargas-simplificadas/_run_iterativo/`.

Run nova: `_run_iterativo_true_rtol/` (tentativas intermediarias arquivadas no mesmo diretorio).

| Tentativa | Resultado |
|---|---|
| 1 (21:36, `setConvergenceTest`) | Travou (~it 1023, rnorm recursivo ~1e-68). Container do job parado (`docker stop` so deste nome); exit 137, OOMKilled=false. Arquivo: `attempt1_convtest_stuck/` |
| 2 (21:39:09–21:39:47) | KSP CONVERGED_RTOL em 229 its; residual verdadeiro 2,219e-10 > 1e-10; o post-check recusou corretamente (`solve_exit=3`). Pico 1,174 GiB. Arquivo: `attempt2_true_residual_rejected/` |
| 3 (21:42:47–21:44:37) | Sucesso. cid `c7961b7b3d6e4f61…`. MemAvailable 23,543 GiB; teto `--memory`=`--memory-swap`=19,543 GiB; 4 CPUs |

**Tentativa 3 (MEASURED):**

- `solve_exit=0`, `postprocess_exit=0`, `true_residual_meets_rtol=true`
- TRUE `||Au−b||/||b||` = 9,577896e-11 ≤ 1e-10
- 668 iteracoes acumuladas, 11 residual replacements (apos 310 its o residual verdadeiro oscilava em ~1,00–1,10e-10; caiu abaixo do rtol no 12º ciclo)
- `CONVERGED_RTOL` (reason=2), rtol interno 1e-12, rtol pedido 1e-10
- pico cgroup 1594101760 B = 1,485 GiB (HWM processo 1460,4 MiB), bem abaixo do teto
- parede 21:42:47–21:44:37 America/Sao_Paulo (~110 s); solve 72,5 s ate o manifesto de campos/RF
- max\|U\| = 0,390038893665 mm no no INP 3233, xyz = (4,81980, 245,43863, −199,70251) mm
- RF: RFx=−1,0248e-6 N, RFy=4999,927425905 N, RFz=−1,5389e-9 N
- aplicado Fy=−4999,927425905 N; Σ(R+F) = (6,97e-10, 2,08e-10, −1,54e-9) N
- momentos Σ(R+F) na origem ~ 1e-7 N·mm
- vs totais do .dat CCX de Luiz: ΔFy=4,26e-4 N (CCX imprime 4999,927)
- C4 proposta (nao aceite) vs U do .dat: max‖du‖=7,0198e-7 mm, u_ref=0,390039 mm, 0 violacoes em 125390 nos
- energia ½f·u (so CLOAD): E_fc=320,01213862, E_ccx=320,01212479, rel=4,324e-8

**Extremos S1/S2/S3/VM (MPa, DG0), com localizacao (centroide mm) e eid INP.** Nenhum extremo global cai nos 183 elementos com q<0,1; portanto GLOBAL = fora dos pobres.

GLOBAL (e fora de q<0,1):

- S1 max = 47,74775, eid 648176, (−154,751, −155,954, 11,509)
- S1 min = −18,22391, eid 662534, (33,589, 218,946, 10,595)
- S2 max = 18,50995, eid 645914, (−154,769, −157,545, 11,048)
- S2 min = −28,11341, eid 662534, (33,589, 218,946, 10,595)
- S3 max = 14,04284, eid 651104, (218,594, 35,634, 10,613)
- S3 min = −69,93191, eid 662534, (33,589, 218,946, 10,595)
- VM max = 62,75435, eid 645534, (33,362, 216,988, 11,572)
- VM min = 0, eid 305065, (31,874, −47,855, 13,642)

So dentro dos 183 q<0,1:

- S1 max = 28,28257, eid 648887, (219,083, 33,184, 10,401)
- S1 min = −8,49287, eid 481414, (26,856, 232,407, 44,647)
- S2 max = 12,88929, eid 674181, (158,519, −154,734, 10,434)
- S2 min = −18,71272, eid 673604, (−36,557, 218,437, 10,490)
- S3 max = 7,76503, eid 674183, (158,370, −154,934, 10,516)
- S3 min = −33,03094, eid 649176, (−196,554, 102,230, 10,358)
- VM max = 43,95603, eid 649148, (−37,328, 218,654, 10,187)
- VM min = 0,88785, eid 673557, (53,499, 56,797, 63,455)

### 1) Calibracao CCX C4–C6 no B1

Gerados INP B1 n=8 e n=16 (o n=16 e leve). CCX 2.23 em `/opt/CalculiX/bin/ccx`, imagem `ubuntu-calculix:latest`, container temporario `fenics-fc-b1-ccx-noite` (docker run, pasta temporaria, removido depois). O container vivo `ubuntu-calculix` permaneceu Up (2 dias) e nao recebeu exec/stop.

Decks com `%.16e` (minusculo, linhas longas) falharam: `*ERROR reading *NODE`, ccx_exit=201, `.dat` vazio. Copia CCX compactada (decimais curtos / E maiusculo, linhas ≤ 34 caracteres), `*HEADING`, `*NODE PRINT U` e `*EL PRINT, ELSET=Eall, GLOBAL=YES S`. Depois: ccx_exit=0, n=8 em ~1 s (21:42:47–21:42:48), n=16 em ~12 s (21:42:48–21:43:00). Tensao no `.dat` em uma linha por elemento: `eid ipt sxx syy szz sxy sxz syz` (1 IP no C3D4).

C4 conhecida (proposta): `|du_k| ≤ 1e-6*u_ref + 1e-5*|u_CCX_k|` por componente. C5 qualitativo = like-for-like DG0 vs *EL PRINT no IP. Tolerancia numerica C5/C6 **nao encontrada**; so diferencas medidas.

| nivel | nos / elems | u_ref (mm) | max‖du‖ (mm) | max‖du‖/u_ref | violacoes C4 proposta | s_ref \|S\|_∞ CCX (MPa) | max\|dS\| (MPa) | max\|dS\|/s_ref |
|---|---|---|---|---|---|---|---|---|
| n=8 | 6561 / 30720 | 2,010855 | 5,049e-7 | 2,511e-7 | 0 (passa a proposta) | 77,63637 | 5,234e-6 | 6,741e-8 |
| n=16 | 46529 / 245760 | 2,116073 | 5,364e-7 | 2,535e-7 | 0 (passa a proposta) | 77,54243 | 4,379e-5 | 5,647e-7 |

Ordem Voigt usada: xx,yy,zz,sxy,sxz,syz. vtk unmatched = 0. As diferencas de U/S estao na casa da precisao de impressao do `.dat` do CCX (~7 algarismos). JSON: `FC/noite-2026-11-11/b1/ccx_fc_compare.json`.

### 2) Esfera B1 — recorte geometrico

Selecao por baricentro substituida (o monitor agora reporta os dois metodos). Volume tet ∩ esfera por classificacao in/out + subdivisao mid-edge (profundidade 5; folha se hmax ≤ 0,05 R). CG+GAMG, `true_residual_meets_rtol=true` nos tres niveis. σ_xx exato no centro da esfera = 52,5 MPa; volume exato 4/3 π R³ ≈ 33,5103 mm³.

| n | its | true rel | geo rel_err | geo vol (mm³) | geo tets | bary rel_err | bary vol (mm³) |
|---|---|---|---|---|---|---|---|
| 4 | 29 | 4,981e-11 | −0,18913 | 33,4924 | 12 | −0,14701 | 41,6667 |
| 8 | 27 | 8,533e-11 | −0,05907 | 33,5000 | 38 | −0,08874 | 41,6667 |
| 16 | 30 | 5,321e-11 | −0,01620 | 33,5075 | 248 | −0,01628 | 33,8542 |

Em n=4 o erro geometrico e pior que o de baricentro porque o conjunto bary (vol 41,67 mm³) nao e a bola. A partir de n=8 o geometrico e melhor ou equivalente, e os volumes geometricos aproximam a bola exata. Propostas, nao gate. JSON: `FC/noite-2026-11-11/b1/sphere_summary.json`.

### 3) Diagnostico de malha da roda (so relatorio)

Sobre a malha da run auxiliar (476960 C3D4). q = 6√2 V / l_rms³, limiar 0,1. JSON: `FC/noite-2026-11-11/diag/mesh_quality_diag-run_iterativo.json`.

- 183 elementos com q<0,1 (**bate com a contagem de Luiz**)
- q min=2,619e-4, media=0,7062, p50=0,7431, p99=0,9588, max=0,99997
- histograma: 76 com q<0,05; 107 em [0,05, 0,1); 690 em [0,1, 0,2)
- os 183 estao todos no ELSET `Evolumes` (nenhum outro elset os contem)
- espaciais dos pobres: r media 212,6 mm (faixa 42,8–234,6; aros mais externos; r medio global=141,5 mm); z media 27,2 mm (faixa −4,46 a 64,09); theta quase volta completa
- nenhum extremo global de S1/S2/S3/VM cai em elemento pobre
- distancia (centroides) do extremo global ao pobre mais proximo: S1 max 3,515 mm (eid pobre 673180, q=0,080); S3 min / S1 min 3,806 mm (eid 649198, q=0,0018); S3 max 2,507 mm (eid 648887, q=0,095, 11 pobres a menos de 5 mm); VM max 4,547 mm (eid 649198). Mediana das distancias a todos os pobres: ~264–303 mm

### 4) Commits e documentos

Commit local apenas (sem push). `docs/FC_PROFILE.md` atualizado. Esta nota em `FC/nota-execucao-FC-noite-2026-11-11.md` e copia em `docs/`, tambem no Drive (pasta FC `1EWCNZ5WEzVoJjhOkP7IswSS2ZiwoiwtF`). Arquivos sujos do usuario no tree (`D input/beam_nset.inp`, `D input/teste.inp`, `?? docs/Information_TBD.odt`) nao foram incluídos.

---

## ESTIMATED / interpretacao (nao medido como fato novo)

- O piso ~1,00e-10 do residual verdadeiro durante varias reposicoes, neste mesh com 183 tets muito pobres, e compativel com deriva do residual recursivo do CG + condicionamento. Nao e medicao de cond(A).
- As diferencas U/S FC vs CCX no B1, na casa de 1e-7 mm e 1e-6–1e-5 MPa, sao compativeis com o arredondamento do `.dat` do CCX. Nao e prova de identidade dos solvers alem do que foi impresso.
- Extrapolacoes de memoria LU do lote 3 permanecem estimativas; esta run iterativa mediu 1,485 GiB.

---

## Pendente (sem decisao tomada nesta noite)

- Tolerancias numericas C5/C6 da revisao de Luiz (`revisao-Luiz-CAE-FC-criterios-gate-e-INP-teste-2026-11-11-1700.md` nao foi encontrada localmente nem por busca no Drive). Nao inventadas.
- Congelar ou nao qualquer criterio (C2–C8, C4, energia, esfera).
- O rtol autorizado 1e-10 foi cumprido, mas ficou justos nesta malha (muitas reposicoes em torno de 1,00e-10). Qualquer folga ou limpeza dos 183 pobres e decisao de Valdir/Luiz.
- Incluir meshio/VTU na imagem compartilhada `edss/fenicsx:dev` (nao feito).
- Gate V / caso definitivo E=70000: nao aprovado, nao executado.
