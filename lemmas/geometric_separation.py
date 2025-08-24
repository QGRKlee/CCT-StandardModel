# === [LEMMA-7 GEOMETRIC-SEPARATION STUB v1] ===
# Minimal utilities to test "geometric separation ⇒ emergent commutativity".
from __future__ import annotations
import numpy as np

# --- expected hooks from your codebase (adjust names if different) ---
try:
    from shells import build_shell_projectors
except Exception:
    build_shell_projectors = None

try:
    from isoclinic import get_sigma, get_S
except Exception:
    def get_sigma(): raise NotImplementedError("Wire isoclinic.get_sigma()")
    def get_S(): raise NotImplementedError("Wire isoclinic.get_S()")

def _comm_norm(A: np.ndarray, B: np.ndarray, ord: str = "fro") -> float:
    return np.linalg.norm(A @ B - B @ A, ord=ord)

def _is_orth_proj(P: np.ndarray, tol: float = 1e-10) -> bool:
    return (
        np.allclose(P.T.conj(), P, atol=tol) and
        np.allclose(P @ P, P, atol=tol)
    )

def check_geometric_separation(tol: float = 1e-12) -> dict:
    if build_shell_projectors is None:
        raise RuntimeError("Provide shells.build_shell_projectors() that returns list of projectors")

    P_list = build_shell_projectors()
    assert len(P_list) > 0, "No shell projectors returned"

    for idx, P in enumerate(P_list):
        assert _is_orth_proj(P), f"Projector {idx} is not orthogonal idempotent"

    S = get_S()
    sigma = get_sigma()

    def _pow_eq(M, k):
        Mk = np.eye(M.shape[0], dtype=M.dtype)
        for _ in range(k): Mk = Mk @ M
        return Mk

    orders_ok = np.allclose(_pow_eq(S, 5), np.eye(S.shape[0]), atol=tol) and \
                np.allclose(_pow_eq(sigma, 4), np.eye(sigma.shape[0]), atol=tol)

    shell_stats = []
    worst = 0.0
    for i, P in enumerate(P_list):
        SPS = P @ S @ P
        PsigP = P @ sigma @ P
        c = _comm_norm(SPS, PsigP, ord="fro")
        worst = max(worst, float(c))
        shell_stats.append({"shell": i, "comm_norm": float(c)})

    passed = (worst <= tol)

    return {
        "orders_ok": bool(orders_ok),
        "per_shell_comm_norms": shell_stats,
        "max_comm_norm": float(worst),
        "passed_geometric_separation": bool(passed),
        "tol": float(tol),
    }

if __name__ == "__main__":
    out = check_geometric_separation()
    print(out)
# === [LEMMA-7 GEOMETRIC-SEPARATION STUB v1 END] ===
