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

Aceita: *NODE, *ELEMENT (TYPE=C3D4 obrigatorio), *NSET, *ELSET (incluindo alias), *MATERIAL, *ELASTIC, *SOLID SECTION, *BOUNDARY (forma `nset,dof` ou `nset,dof1,dof2,valor`), *CLOAD, um *STEP, *STATIC (um incremento, periodo 1, tempo inicial 0, sem heranca), *END STEP, *NODE FILE, *EL FILE, *NODE PRINT, *OUTPUT.

Ignorados **explicitamente** (manifesto, nunca silencio): *PHYSICAL CONSTANTS, *EXPANSION, parametro `INC` de *STEP.

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
- CLOAD acumula no DOF localizado por `locate_dofs_geometrical` (mesmo criterio das BC). Duas linhas no mesmo no/DOF somam.
- `*BOUNDARY`: CalculiX e `no, dof_ini, dof_fim[, valor]`. Tres colunas com valor flutuante sao ambiguas e rejeitadas no FC; valor nao nulo exige 4 colunas.
- Matriz completa (sem BC) para residual; matriz com BC so no solve. RF = `K_full u - f_ext` apenas nos DOFs prescritos. JSON traz `sum_R_plus_F` e momentos na origem. Desequilibrio acima da tolerancia numerica de consistencia recusa sucesso.
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
