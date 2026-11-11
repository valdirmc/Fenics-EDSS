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
