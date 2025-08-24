"""
Lemma (Geometric separation ⇒ effective commutation).
If P is the projector onto any E8 24‑cell shell and S (order 5) and σ (order 4) are the
Cl(8) rotors used in the paper, then on the projected subspace we have
  || [P S P, P σ P] ||_F ≤ tol
for the same numerical tolerance used in isoclinic.py.

This formalizes the intuition that disjoint shells yield commuting actions “in practice”,
matching the twist/discrete‑curvature equivalence used in the construction.  See also
Fang–Clawson–Irwin for the twist ↔ discrete curvature angle matching. 
"""


# geometric_separation.py
# Machine check for “geometric separation”: S (order 5) and σ (order 4) commute
# and, when restricted to each 24‑cell shell, their effective commutator is ~0.

# --- path shim (so we can import from ../src) ---
import os, sys, numpy as np
_here = os.path.dirname(__file__)
src_path = os.path.abspath(os.path.join(_here, "..", "..", "src"))
if src_path not in sys.path:
    sys.path.append(src_path)

from shells import partition_shells, validate_shells
from isoclinic import (
    construct_standard_isoclinic_S,
    construct_standard_isoclinic_sigma,
)
# optional data dir
DATA_DIR = os.path.abspath(os.path.join(src_path, "data"))

def _orth_proj_from_vectors(vs: np.ndarray) -> np.ndarray:
    """
    vs: array of shape (m, 8) with m vectors in R^8 (e.g., a 24‑cell's 24 roots).
    returns the 8x8 orthogonal projector onto span(vs).
    """
    # Orthonormal columns spanning span(vs)
    # Work with the transpose to QR the column space of vs^T
    # vs.T has shape (8, m). QR gives Q (8 x k) with orthonormal columns.
    Q, _ = np.linalg.qr(vs.T)  # economic QR
    return Q @ Q.T  # 8x8 projector

def _mat_pow(M: np.ndarray, k: int) -> np.ndarray:
    return np.linalg.matrix_power(M, k)

def check_geometric_separation(tol: float = 1e-9):
    # 1) Build shells (10 shells, each expected to have 24 root vectors in R^8)
    shells = partition_shells()
    validate_shells(shells)  # raises if anything is off

    # 2) Get the canonical rotors/matrices S (order 5) and sigma (order 4)
    S = construct_standard_isoclinic_S()
    sigma = construct_standard_isoclinic_sigma()

    I = np.eye(8)

    # 3) Global algebraic checks
    orders_ok = (
        np.allclose(_mat_pow(S, 5), I, atol=tol) and
        np.allclose(_mat_pow(sigma, 4), I, atol=tol)
    )
    comm_global = np.linalg.norm(S @ sigma - sigma @ S, ord="fro")

    # 4) Per‑shell projectors and effective commutator
    shell_stats = []
    worst = 0.0

    for idx, sh in enumerate(shells):
        # sh expected shape: (24, 8)
        P = _orth_proj_from_vectors(sh)  # 8x8 projector onto that 24‑cell's span
        # Restricted (effective) commutator on the shell’s subspace:
        C_eff = P @ S @ P @ sigma @ P - P @ sigma @ P @ S @ P
        c = np.linalg.norm(C_eff, ord="fro")
        worst = max(worst, float(c))
        shell_stats.append({"shell": idx, "comm_norm": float(c)})

    passed = (orders_ok and comm_global <= 1e-12 and worst <= 1e-9)

    # 5) Report
    out = {
        "orders_ok": bool(orders_ok),
        "global_comm_norm": float(comm_global),
        "per_shell_comm_norms": shell_stats,
        "max_comm_norm_over_shells": float(worst),
        "passed_geometric_separation": bool(passed),
        "tol": float(tol),
    }
    return out
def run_geometric_separation(tol: float = 1e-12):
    out = check_geometric_separation(tol=tol)
    # normalize a stable, machine‑parsable summary
    return {
        "orders_ok": bool(out.get("orders_ok", False)),
        "passed_geometric_separation": bool(out.get("passed_geometric_separation", False)),
        "max_comm_norm": float(out.get("max_comm_norm", float("inf"))),
        "tol": float(out.get("tol", tol)),
        "per_shell_comm_norms": out.get("per_shell_comm_norms", []),
    }


if __name__ == "__main__":
    result = check_geometric_separation()
    print(result)
