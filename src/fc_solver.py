"""Solvers do perfil FC: LU (padrao) e CG+GAMG opt-in. Sem fallback silencioso.

CG+GAMG exige near-nullspace elastico (6 modos de corpo rigido) na matriz.
Nao convergiu => RuntimeError, exit != 0; nunca cai para LU.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from job_io import normalize_rtol, normalize_solver

SOLVER_LU = "lu"
SOLVER_CG_GAMG = "cg-gamg"
ALLOWED_SOLVERS = (SOLVER_LU, SOLVER_CG_GAMG)
DEFAULT_RTOL = 1e-10
MAX_RTOL = 1e-10
DEFAULT_MAX_IT = 10000

_REASON_NAMES = {
    1: "CONVERGED_RTOL_NORMAL",
    2: "CONVERGED_ATOL",
    3: "CONVERGED_ITOL",
    4: "CONVERGED_ITS",
    9: "CONVERGED_ITERATING",  # placeholder; real codes below
    -2: "DIVERGED_NULL",
    -3: "DIVERGED_ITS",
    -4: "DIVERGED_DTOL",
    -5: "DIVERGED_BREAKDOWN",
    -6: "DIVERGED_BREAKDOWN_BICG",
    -7: "DIVERGED_NONSYMMETRIC",
    -8: "DIVERGED_INDEFINITE_PC",
    -9: "DIVERGED_NANORINF",
    -10: "DIVERGED_INDEFINITE_MAT",
    -11: "DIVERGED_PC_FAILED",
}


def reason_name(reason: int) -> str:
    try:
        from petsc4py import PETSc
        mapping = {int(getattr(PETSc.KSP.ConvergedReason, n)): n
                   for n in dir(PETSc.KSP.ConvergedReason) if n.isupper()}
        if int(reason) in mapping:
            return mapping[int(reason)]
    except Exception:
        pass
    return _REASON_NAMES.get(int(reason), "REASON_{0}".format(reason))


def _petsc_array(vec):
    if hasattr(vec, "array"):
        return vec.array
    return vec.getArray()


def build_elastic_near_nullspace(V, A):
    """6 modos de corpo rigido (3 translacoes + 3 rotacoes) como near-nullspace.

    Com *BOUNDARY a nulidade exata some; GAMG precisa de setNearNullSpace, nao
    setNullSpace. Serial (perfil FC).
    """
    from petsc4py import PETSc

    bs = int(V.dofmap.index_map_bs)
    nblock = int(V.dofmap.index_map.size_local)
    ndof = nblock * bs
    coords = np.asarray(V.tabulate_dof_coordinates(), dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] < 3:
        raise RuntimeError("tabulate_dof_coordinates forma inesperada: {0}".format(coords.shape))
    if coords.shape[0] == nblock:
        xyz = coords[:, :3]
    elif coords.shape[0] == ndof and bs > 1:
        xyz = coords[::bs, :3]
        if xyz.shape[0] != nblock:
            raise RuntimeError("coordenadas vs blocos: {0} vs {1}".format(xyz.shape[0], nblock))
    else:
        raise RuntimeError(
            "coordenadas ({0}) nao batem com blocos ({1}) * bs ({2})".format(
                coords.shape[0], nblock, bs
            )
        )
    if bs != 3:
        raise RuntimeError("near-nullspace elastico 3D exige block size 3 (bs={0}).".format(bs))

    mat_size = A.getLocalSize()[0] if hasattr(A, "getLocalSize") else A.getSize()[0]
    if mat_size != ndof:
        # alguns Mat reportam so linhas globais em serial
        gsize = A.getSize()[0]
        if gsize != ndof:
            raise RuntimeError(
                "tamanho da matriz ({0}/{1}) != ndof {2}".format(mat_size, gsize, ndof)
            )

    modes = []
    # translacoes
    for axis in range(3):
        v = A.createVecLeft()
        v.set(0.0)
        arr = np.asarray(_petsc_array(v))
        if arr.size < ndof:
            raise RuntimeError("vec size {0} < ndof {1}".format(arr.size, ndof))
        arr[:ndof] = 0.0
        arr[axis:ndof:3] = 1.0
        v.assemble()
        modes.append(v)
    # rotacoes omega x r no origem
    rot = (
        lambda x, y, z: (0.0, -z, y),
        lambda x, y, z: (z, 0.0, -x),
        lambda x, y, z: (-y, x, 0.0),
    )
    for fn in rot:
        v = A.createVecLeft()
        v.set(0.0)
        arr = np.asarray(_petsc_array(v))
        arr[:ndof] = 0.0
        for i, (x, y, z) in enumerate(xyz):
            vx, vy, vz = fn(float(x), float(y), float(z))
            arr[3 * i] = vx
            arr[3 * i + 1] = vy
            arr[3 * i + 2] = vz
        v.assemble()
        modes.append(v)

    # ortonormaliza (PETSc nao faz isso sozinho em NullSpaceCreate)
    basis = []
    for v in modes:
        w = v.copy()
        for u in basis:
            w.axpy(-float(w.dot(u)), u)
        nrm = float(w.norm())
        if nrm < 1e-14:
            w.destroy()
            continue
        w.scale(1.0 / nrm)
        basis.append(w)
    if len(basis) < 6:
        # malha degenerada ou coordenadas nulas: ainda tenta com o que houver
        pass
    ns = PETSc.NullSpace().create(constant=False, vectors=basis, comm=A.getComm())
    A.setNearNullSpace(ns)
    # PETSc.Mat nao aceita atributos Python; devolve refs para o caller reter.
    return {"nullspace": ns, "basis": basis, "modes": modes}


def true_residual(A, uvec, bvec) -> Tuple[float, float, float]:
    """Retorna (||r||, ||b||, ||r||/||b||) com r = A u - b."""
    r = bvec.duplicate()
    A.mult(uvec, r)
    r.axpy(-1.0, bvec)
    rn = float(r.norm())
    bn = float(bvec.norm())
    r.destroy()
    rel = rn / bn if bn > 0.0 else rn
    return rn, bn, rel


def _configure_cg_gamg(ksp, A, V, rtol: float, max_it: int) -> None:
    from petsc4py import PETSc

    refs = build_elastic_near_nullspace(V, A)
    ksp.setOperators(A)
    ksp.setType(PETSc.KSP.Type.CG)
    ksp.setTolerances(rtol=rtol, atol=0.0, divtol=1e5, max_it=int(max_it))
    # residual verdadeira ||b-Ax||, nao a precondicionada
    try:
        ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
    except Exception:
        ksp.setNormType(PETSc.KSP.NormType.NORM_UNPRECONDITIONED)
    pc = ksp.getPC()
    pc.setType(PETSc.PC.Type.GAMG)
    try:
        pc.setGAMGType("agg")
    except Exception:
        pass
    ksp.setOptionsPrefix("fc_")
    opts = PETSc.Options()
    opts["fc_pc_gamg_agg_nsmooths"] = 1
    opts["fc_pc_gamg_coarse_eq_limit"] = 100
    opts["fc_pc_gamg_threshold"] = 0.001
    ksp.setFromOptions()
    try:
        ksp.setErrorIfNotConverged(False)  # nos mesmos checamos e erramos
    except Exception:
        pass
    return refs


def solve_linear(A, b, uh, V, *, solver: str = SOLVER_LU, rtol: float = DEFAULT_RTOL,
                 max_it: int = DEFAULT_MAX_IT, debug: int = 0) -> Dict[str, Any]:
    """Resolve A u = b. LU padrao; cg-gamg opt-in. Sem fallback."""
    from petsc4py import PETSc

    kind = normalize_solver(solver)
    rtol = normalize_rtol(rtol, kind)
    ksp = PETSc.KSP().create(V.mesh.comm)
    history: List[Dict[str, Any]] = []
    ns_refs = None  # retido ate o fim do solve (near-nullspace Vecs)

    def _monitor(ks, it, rnorm):
        rec = {"iteration": int(it), "rnorm": float(rnorm)}
        history.append(rec)
        if debug or it == 0 or it % 25 == 0:
            print("KSP it={0} rnorm={1:.6e}".format(it, rnorm), flush=True)

    if kind == SOLVER_LU:
        ksp.setOperators(A)
        ksp.setType("preonly")
        ksp.getPC().setType("lu")
        try:
            ksp.setErrorIfNotConverged(True)
        except Exception:
            pass
    else:
        ns_refs = _configure_cg_gamg(ksp, A, V, rtol, max_it)
        ksp.setMonitor(_monitor)

    uvec = uh.x.petsc_vec
    try:
        ksp.solve(b, uvec)
    except Exception as exc:
        info = {
            "type": kind,
            "rtol": rtol,
            "max_it": int(max_it) if kind == SOLVER_CG_GAMG else None,
            "error": str(exc),
            "history": history,
        }
        err = RuntimeError(
            "Solver PETSc falhou ({0}): {1}. Sem fallback para LU.".format(kind, exc)
        )
        err.solver_info = info
        raise err from exc

    uh.x.scatter_forward()
    reason = int(ksp.getConvergedReason())
    its = int(ksp.getIterationNumber())
    rn, bn, rel = true_residual(A, uvec, b)
    rname = reason_name(reason)
    info = {
        "type": kind,
        "pc": "lu" if kind == SOLVER_LU else "gamg",
        "ksp": "preonly" if kind == SOLVER_LU else "cg",
        "rtol": rtol if kind == SOLVER_CG_GAMG else None,
        "max_it": int(max_it) if kind == SOLVER_CG_GAMG else None,
        "iterations": its,
        "converged_reason": reason,
        "converged_reason_name": rname,
        "residual_norm": rn,
        "rhs_norm": bn,
        "true_residual_rel": rel,
        "near_nullspace": "elastic-6-rbm" if kind == SOLVER_CG_GAMG else None,
        "history_len": len(history),
        "history_tail": history[-5:],
    }
    if debug or kind == SOLVER_CG_GAMG:
        print(
            "Solver {0}: reason={1} ({2}) its={3} ||r||={4:.6e} ||b||={5:.6e} ||r||/||b||={6:.6e}".format(
                kind, reason, rname, its, rn, bn, rel
            ),
            flush=True,
        )
    if kind == SOLVER_CG_GAMG and reason <= 0:
        err = RuntimeError(
            "CG+GAMG nao convergiu (reason={0} {1}, its={2}/{3}, "
            "||r||/||b||={4:.6e}, rtol={5:g}). Sem fallback para LU.".format(
                reason, rname, its, max_it, rel, rtol
            )
        )
        err.solver_info = info
        raise err
    from fc_assembly import check_solution_finite
    check_solution_finite(uh, ksp)
    info["ok"] = True
    if ns_refs:
        info["near_nullspace_basis"] = len(ns_refs.get("basis") or [])
    return info
