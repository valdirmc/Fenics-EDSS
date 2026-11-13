# Tabelas ponto a ponto FenicsEDSS × CCX 2.23 (Sprint 2)

Comparação like-for-like, por ID de nó (`*_nos.csv`) e de elemento (`*_elementos.csv`), do FenicsEDSS (LU) contra o CalculiX 2.23 (SPOOLES, `*NODE PRINT` / `*EL PRINT` no `.dat`). Casamento por nome nas tensões (σxx, σyy, σzz, σxy, σxz, σyz).

As colunas `C4_*_proposta` e `C5_*_proposta` usam as referências **PROPOSTAS** por Luiz em 13/11 (T2: u_ref = 0,222 mm, σ_ref = 289,473684 MPa; T0: u_ref = 4,857e-4 mm, σ_ref = 10 MPa). **Não estão congeladas** e não aprovam o Gate V.

COORD aparece nas colunas `dx_coord`, `dy_coord`, `dz_coord` (baricentro Fenics versus `*EL PRINT COORD` do CCX). A fórmula proposta por Luiz, 1e-9·L_ref + 2·ulp, é documentada no relatório da sprint e **não** é aplicada como critério congelado nestes CSV.
